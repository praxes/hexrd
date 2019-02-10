from __future__ import print_function

import cPickle
import logging
import multiprocessing as mp
import os
import time
from collections import namedtuple

import numpy as np
#np.seterr(over='ignore', invalid='ignore')

import scipy.cluster as cluster
from scipy import ndimage

from hexrd import matrixutil as mutil
from hexrd.constants import sqrt_epsf
from hexrd.xrd import experiment as expt
from hexrd.xrd import indexer as idx
from hexrd.xrd import rotations as rot
from hexrd.xrd import symmetry as sym
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.coreutil import initialize_experiment
from hexrd.xrd import xrdutil
from hexrd.xrd import distortion as dFuncs
from hexrd.fitgrains import get_instrument_parameters

logger = logging.getLogger(__name__)

save_as_ascii = False # FIX LATER...

# just require scikit-learn?
have_sklearn = False
try:
    import sklearn
    vstring = sklearn.__version__.split('.')
    if vstring[0] == '0' and int(vstring[1]) >= 14:
        import sklearn.cluster
        from sklearn.metrics.pairwise import pairwise_distances
        have_sklearn = True
except ImportError:
    pass

have_parallel_dbscan = False
try:
    import parallel_dbscan
    have_parallel_dbscan = True
except ImportError:
    pass

def generate_orientation_fibers(eta_ome, chi, threshold, seed_hkl_ids,
                                fiber_ndiv, filt_stdev=1.0, ncpus=1):
    """
    From ome-eta maps and hklid spec, generate list of
    quaternions from fibers
    """
    # seed_hkl_ids must be consistent with this...
    pd_hkl_ids = eta_ome.iHKLList[seed_hkl_ids]

    # grab angular grid infor from maps
    del_ome = eta_ome.omegas[1] - eta_ome.omegas[0]
    del_eta = eta_ome.etas[1] - eta_ome.etas[0]

    # labeling mask
    structureNDI_label = ndimage.generate_binary_structure(2, 1)

    # crystallography data from the pd object
    pd = eta_ome.planeData
    hkls = pd.hkls
    tTh  = pd.getTTh()
    bMat = pd.latVecOps['B']
    csym = pd.getLaueGroup()

    params = {
        'bMat': bMat,
        'chi': chi,
        'csym': csym,
        'fiber_ndiv': fiber_ndiv
    }

    ############################################
    ##    Labeling of spots from seed hkls    ##
    ############################################

    qfib     = []
    input_p  = []
    numSpots = []
    coms     = []
    for i in seed_hkl_ids:
        # First apply filter
        this_map_f = -ndimage.filters.gaussian_laplace(eta_ome.dataStore[i],
                                                       filt_stdev)

        labels_t, numSpots_t = ndimage.label(this_map_f > threshold,
                                             structureNDI_label )
        coms_t = np.atleast_2d(
            ndimage.center_of_mass(
                this_map_f,
                labels=labels_t,
                index=np.arange(1, np.amax(labels_t)+1)
                )
            )
        numSpots.append(numSpots_t)
        coms.append(coms_t)
        pass

    for i in range(len(pd_hkl_ids)):
        for ispot in range(numSpots[i]):
            if not np.isnan(coms[i][ispot][0]):
                ome_c = eta_ome.omeEdges[0] + (0.5 + coms[i][ispot][0])*del_ome
                eta_c = eta_ome.etaEdges[0] + (0.5 + coms[i][ispot][1])*del_eta
                input_p.append(np.hstack([hkls[:, pd_hkl_ids[i]],
                                          tTh[pd_hkl_ids[i]],
                                          eta_c, ome_c]))
                pass
            pass
        pass

    # do the mapping
    start = time.time()
    qfib = None
    if ncpus > 1:
        # multiple process version
        pool = mp.Pool(ncpus, discrete_fiber_init, (params, ))
        qfib = pool.map(discrete_fiber_reduced, input_p) # chunksize=chunksize)
        pool.close()
    else:
        # single process version
        global paramMP
        discrete_fiber_init(params)
        qfib = map(discrete_fiber_reduced, input_p)
        paramMP = None # clear paramMP

    elapsed = (time.time()-start)
    logger.info('fiber generation took %.3f seconds', elapsed)

    return np.hstack(qfib)


def discrete_fiber_init(params):
    global paramMP
    paramMP = params

def discrete_fiber_reduced(params_in):
    """
    input parameters are [hkl_id, com_ome, com_eta]
    """
    bMat       = paramMP['bMat']
    chi        = paramMP['chi']
    csym       = paramMP['csym']
    fiber_ndiv = paramMP['fiber_ndiv']

    hkl = params_in[:3].reshape(3, 1)

    gVec_s = xfcapi.anglesToGVec(
        np.atleast_2d(params_in[3:]),
        chi=chi,
        ).T

    tmp = mutil.uniqueVectors(
        rot.discreteFiber(
            hkl,
            gVec_s,
            B=bMat,
            ndiv=fiber_ndiv,
            invert=False,
            csym=csym)[0]
    )

    return tmp


_clustering_option = namedtuple('_clustering_option', 'fn')
_clustering_algorithm_dict = {}


class ClusteringError(Exception):
    """A ClusteringError that contains a reason (message) of failure as well as an
    optional alternative method of clustering to use."""
    def __init__(self, message, alternative=None):
        super(ClusteringError, self).__init__(message)
        self.alternative = alternative


def get_supported_clustering_algorithms():
    """get a list of the supported clustering algorithms"""
    # note: this is used by the yaml parser to know available
    #       options
    return _clustering_algorithm_dict.keys()


def clustering_algorithm(key):
    """A decorator that registers clustering algorithms automagically.

    A valid cluster algorithm must return an array of [int] indices
    that map each fiber to its corresponding cluster id.

    Valid clusters id are positive integers, with 0 reserved for noise
    fibers where applicable.

    """
    def wrapper(fn):
        assert key not in _clustering_algorithm_dict
        val = _clustering_option(fn)
        _clustering_algorithm_dict.update({key: val })
        return fn

    return wrapper


def _check_dbscan():
    if not have_sklearn:
        msg = 'sklearn >= 0.14 required to use dbscan; using fclusterdata'
        raise ClusteringError(msg, alternative='fclusterdata')


def _normalize_labels_from_dbscan(labels):
    """returns labels normalized as a numpy array with 1-based indices and -1 for
    noise. Input are the labels as returned by sklearn.cluster.dbscan

    """
    cl = np.array(labels, dtype=int)
    noise_points = cl == -1
    cl += 1
    cl[noise_points] = -1 # remark it as -1
    logger.info("dbscan found %d noise points", sum(noise_points))
    return cl


def _compute_centroids_split(cl, qfib_r, qsym):
    """generic compute centroids (clusters may be split)"""
    if np.any(cl == -1):
        nblobs = len(np.unique(cl)) - 1
    else:
        nblobs = len(np.unique(cl))

    qbar = np.zeros((4, nblobs))
    for i in range(nblobs):
        cluster_indices = (cl == i + 1)
        this_cluster = qfib_r[:, cluster_indices]
        npts = sum(cluster_indices)
        this_cluster = this_cluster.reshape(4, npts)
        qbar[:, i] = rot.quatAverageCluster(this_cluster, qsym).flatten()

    return np.atleast_2d(qbar)


def _compute_centroids_dense(cl, qfib_r, qsym):
    """compute centroids when clusters are compact"""
    if np.any(cl == -1):
        nblobs = len(np.unique(cl)) - 1
    else:
        nblobs = len(np.unique(cl))

    qbar = np.zeros((4, nblobs))
    for i in range(nblobs):
        cluster_indices = (cl == i + 1)
        this_cluster = qfib_r[:, cluster_indices]
        qbar[:, i] = np.average(np.atleast_2d(this_cluster), axis=1)
    qbar = sym.toFundamentalRegion(mutil.unitVector(qbar), crysSym=qsym)
    return np.atleast_2d(qbar)


def _handle_duplicate_orientations(qbar, qsym, cl_radius):
    """removes duplicate orienations within a tolerance"""

    if qbar.size > 4:
        # need nblobs coming in to detect duplicates
        # WARNING: qbat assumed to be 2-d with shape (4, n)
        nblobs = qbar.shape[1]
        logger.info('\tchecking for duplicate orientations...')
        def quat_distance(x, y):
            return xfcapi.quat_distance(np.array(x, order='C'), np.array(y, order='C'), qsym)
        cl = cluster.hierarchy.fclusterdata(qbar.T, np.radians(cl_radius),
                                            criterion='distance',
                                            metric=quat_distance)
        nblobs_new = len(np.unique(cl))
        if nblobs_new < nblobs:
            logger.info("\tfound %d duplicates within %f degrees"
                        % (nblobs - nblobs_new, cl_radius))

            # if duplicates found, average the duplicates
            tmp = np.zeros((4, nblobs_new))
            for i in range(nblobs_new):
                # TODO: this could be simplified and made faster.
                npts = sum(cl == i + 1)
                duplicates = qbar[:, cl == i+1].reshape(4, npts)
                tmp[:,i] = rot.quatAverageCluster(duplicates, qsym).flatten()
                pass
            qbar = tmp
            pass
        pass
    return qbar


@clustering_algorithm('fclusterdata')
def cluster_fclusterdata(qfib_r, qsym, cl_radius, min_samples):
    num_ors = qfib_r.shape[1]
    if num_ors > 25000:
        if have_sklearn:
            msg = 'Size too big for fclusterdata. Defaulting to ort-dbscan.'
            raise ClusteringError(msg, alternative='ort-dbscan')
        else:
            msg = 'Size too big for fclusterdata and no dbscan present.'
            raise ClusteringError(msg)

    pts = qfib_r.T
    qsym = np.array(qsym.T, order='C').T
    def quat_distance(x, y):
        return xfcapi.quat_distance(np.array(x, order='C'), np.array(y, order='C'), qsym)
    labels = cluster.hierarchy.fclusterdata(pts, np.radians(cl_radius),
                                            criterion='distance', metric = quat_distance)
    return _compute_centroids_split(labels, qfib_r, qsym), labels


@clustering_algorithm('dbscan')
def cluster_dbscan(qfib_r, qsym, cl_radius, min_samples):
    # CAVEAT: the euclidean misorientation of two quaternions  
    # is ~2x smaller than the true quaternion misorientation
    # for magnitudes < ~5deg
    _check_dbscan()
    dbscan = sklearn.cluster.dbscan
    pts = qfib_r.T
    _, labels = dbscan(pts, eps=0.5*np.radians(cl_radius),
                       min_samples=min_samples, metric='minkowski', p=2)
    labels = _normalize_labels_from_dbscan(labels)
    qbar = _compute_centroids_dense(labels, qfib_r, qsym)
    qbar = _handle_duplicate_orientations(qbar, qsym, cl_radius)
    return qbar, labels


@clustering_algorithm('ort-dbscan')
def cluster_ort_dbscan(qfib_r, qsym, cl_radius, min_samples):
    # CAVEAT: the euclidean misorientation of the vector parts of two
    # quaternion is ~(2+eps)x smaller than the true quaternion misorientation 
    # for magnitudes < ~5deg.  The distribution is generally larger than the
    # full quaternion 2-norm, however!!!
    _check_dbscan()
    dbscan = sklearn.cluster.dbscan
    pts = qfib_r[1:, :].T
    _, labels = dbscan(pts, eps=0.5*np.radians(cl_radius),
                       min_samples=min_samples, metric='minkowski', p=2)
    labels = _normalize_labels_from_dbscan(labels)
    qbar = _compute_centroids_dense(labels, qfib_r, qsym)
    qbar = _handle_duplicate_orientations(qbar, qsym, cl_radius)
    return qbar, labels


@clustering_algorithm('sph-dbscan')
def cluster_sph_dbscan(qfib_r, qsym, cl_radius, min_samples):
    _check_dbscan()
    num_ors = qfib_r.shape[1]
    if num_ors > 25000:
        msg = 'Size too big for sph-dbscan. Defaulting to dbscan.'
        raise ClusteringError(msg, alternative='dbscan')

    dbscan = sklearn.cluster.dbscan
    pts = qfib_r.T
    qsym = np.array(qsym.T, order='C').T
    def quat_distance(x, y):
        return xfcapi.quat_distance(np.array(x, order='C'), np.array(y, order='C'), qsym)

    pdist = pairwise_distances(pts, metric=quat_distance, n_jobs=1)
    _, labels = dbscan(pdist, eps=np.radians(cl_radius),
                       min_samples=min_samples, metric='precomputed')
    labels = _normalize_labels_from_dbscan(labels)
    return _compute_centroids_split(labels, qfib_r, qsym), labels


def run_cluster(compl, qfib, qsym, cfg, min_samples=None, compl_thresh=None, radius=None):
    """
    """
    algorithm = cfg.find_orientations.clustering.algorithm

    cl_radius = cfg.find_orientations.clustering.radius
    min_compl = cfg.find_orientations.clustering.completeness

    # check for override on completeness threshold
    if compl_thresh is not None:
        min_compl = compl_thresh

    # check for override on radius
    if radius is not None:
        cl_radius = radius

    t0 = time.clock() # time this

    num_above = sum(np.array(compl) > min_compl)
    if num_above == 0:
        # nothing to cluster
        qbar = cl = np.array([])
    elif num_above == 1:
        # short circuit
        qbar = qfib[:, np.array(compl) > min_compl]
        cl = [1]
    else:
        qfib_r = np.ascontiguousarray(qfib[:, np.array(compl) > min_compl])

        logger.info(
            "Feeding %d orientations above %.1f%% to clustering",
            qfib_r.shape[1], 100*min_compl
            )

        cl_dict = _clustering_algorithm_dict
        if min_samples is None or cfg.find_orientations.use_quaternion_grid is not None:
            min_samples = 1
        cluster_args = [qfib_r, qsym, cl_radius, min_samples]
        while algorithm is not None:
            if algorithm not in cl_dict:
                raise RuntimeError(
                    "Clustering '{0}' not recognized".format(algorithm)
                )
            try:
                logger.info("Trying '%s' over %d orientations",
                            algorithm, qfib_r.shape[1])
                qbar, cl = cl_dict[algorithm].fn(*cluster_args)
                algorithm_used = algorithm
                algorithm = None
            except ClusteringError as error:
                fb = error.alternative
                if fb is None:
                    msg = "Clustering '{0}' failed: {1}\n no fallback."
                    raise RuntimeError(msg.format(algorithm, error))
                msg = "Clustering '{0}' failed: {2}\ntrying '{1}'"
                logger.info(msg.format(algorithm, fb, error))
                algorithm = fb

    logger.info("clustering took %f seconds", time.clock() - t0)
    logger.info(
        "Found %d orientation clusters with >=%.1f%% completeness"
        " and %2f misorientation",
        qbar.size/4,
        100.*min_compl,
        cl_radius
        )

    return qbar, cl


def load_eta_ome_maps(cfg, pd, reader, detector, hkls=None, clean=False):
    fn = os.path.join(
        cfg.working_dir,
        cfg.find_orientations.orientation_maps.file
        )

    if not clean:
        try:
            res = cPickle.load(open(fn, 'r'))
            pd = res.planeData
            available_hkls = pd.hkls.T
            logger.info('loaded eta/ome orientation maps from %s', fn)
            hkls = [str(i) for i in available_hkls[res.iHKLList]]
            logger.info(
                'hkls used to generate orientation maps: %s', hkls)
            return res
        except (AttributeError, IOError):
            return generate_eta_ome_maps(cfg, pd, reader, detector, hkls)
    else:
        logger.info('clean option specified; recomputing eta/ome orientation maps')
        return generate_eta_ome_maps(cfg, pd, reader, detector, hkls)


def generate_eta_ome_maps(cfg, pd, reader, detector, hkls=None):

    available_hkls = pd.hkls.T
    # default to all hkls defined for material
    active_hkls = range(available_hkls.shape[0])
    # override with hkls from config, if specified
    temp = cfg.find_orientations.orientation_maps.active_hkls
    active_hkls = active_hkls if temp == 'all' else temp
    # override with hkls from command line, if specified
    active_hkls = hkls if hkls is not None else active_hkls

    logger.info(
        "using hkls to generate orientation maps: %s",
        ', '.join([str(i) for i in available_hkls[active_hkls]])
        )

    # not ready # eta_ome = xrdutil.EtaOmeMaps(cfg, reader=reader, eta_step=None)
    bin_frames = cfg.find_orientations.orientation_maps.bin_frames
    eta_bins = np.int(2*np.pi / abs(reader.getDeltaOmega())) / bin_frames
    eta_ome = xrdutil.CollapseOmeEta(
        reader,
        pd,
        pd.hkls[:, active_hkls],
        detector,
        nframesLump=bin_frames,
        nEtaBins=eta_bins,
        debug=False,
        threshold=cfg.find_orientations.orientation_maps.threshold
        ).getEtaOmeMaps()

    fn = os.path.join(
        cfg.working_dir,
        cfg.find_orientations.orientation_maps.file
        )
    fd = os.path.split(fn)[0]
    if not os.path.isdir(fd):
        os.makedirs(fd)
    with open(fn, 'w') as f:
        cPickle.dump(eta_ome, f)
    logger.info("saved eta/ome orientation maps to %s", fn)
    return eta_ome


def find_orientations(cfg, hkls=None, clean=False, profile=False):
    """
    Takes a config dict as input, generally a yml document

    NOTE: single cfg instance, not iterator!
    """

    # ...make this an attribute in cfg?
    analysis_id = '%s_%s' % (
        cfg.analysis_name.strip().replace(' ', '-'),
        cfg.material.active.strip().replace(' ', '-')
    )

    # a goofy call, could be replaced with two more targeted calls
    pd, reader, detector = initialize_experiment(cfg)

    # need instrument cfg later on down...
    instr_cfg = get_instrument_parameters(cfg)
    detector_params = np.hstack([
        instr_cfg['detector']['transform']['tilt_angles'],
        instr_cfg['detector']['transform']['t_vec_d'],
        instr_cfg['oscillation_stage']['chi'],
        instr_cfg['oscillation_stage']['t_vec_s'],
        ])
    rdim = cfg.instrument.detector.pixels.size[0]*cfg.instrument.detector.pixels.rows
    cdim = cfg.instrument.detector.pixels.size[1]*cfg.instrument.detector.pixels.columns
    panel_dims = ((-0.5*cdim, -0.5*rdim),
                  ( 0.5*cdim,  0.5*rdim),
                  )
    # UGH! hard-coded distortion...
    if instr_cfg['detector']['distortion']['function_name'] == 'GE_41RT':
        distortion = (dFuncs.GE_41RT,
                      instr_cfg['detector']['distortion']['parameters'],
                      )
    else:
        distortion = None

    # start logger
    logger.info("beginning analysis '%s'", cfg.analysis_name)

    # load the eta_ome orientation maps
    eta_ome = load_eta_ome_maps(cfg, pd, reader, detector, hkls=hkls, clean=clean)

    # KLUDGE: need to enforce + delta omega in maps for optimized paintGrid
    # <JVB 2017/08/16>
    del_ome = np.degrees(eta_ome.omeEdges[1] - eta_ome.omeEdges[0])
    if np.sign(del_ome) == -1:
        assert abs(cfg.image_series.omega.step - del_ome) < sqrt_epsf, \
          "inconsistency with omega spec"
        eta_ome.dataStore = eta_ome.dataStore[:, ::-1, :]    # flip omega axis
        eta_ome.omegas = eta_ome.omegas[::-1]
        eta_ome.omeEdges = eta_ome.omeEdges[::-1]
        pass

    ome_range = (np.min(eta_ome.omeEdges),
                 np.max(eta_ome.omeEdges)
                 )
    try:
        # are we searching the full grid of orientation space?
        qgrid_f = cfg.find_orientations.use_quaternion_grid
        quats = np.load(qgrid_f).T
        logger.info("Using %s for full quaternion search", qgrid_f)
        hkl_ids = None
    except (IOError, ValueError, AttributeError):
        # or doing a seeded search?
        logger.info("Defaulting to seeded search")
        hkl_seeds = cfg.find_orientations.seed_search.hkl_seeds
        hkl_ids = [eta_ome.planeData.hklDataList[i]['hklID'] for i in hkl_seeds]
        hklseedstr = ', '.join(
            [str(i) for i in eta_ome.planeData.hkls.T[hkl_seeds]]
            )
        logger.info(
            "Seeding search using hkls from %s: %s",
            cfg.find_orientations.orientation_maps.file,
            hklseedstr
            )
        quats = generate_orientation_fibers(
            eta_ome,
            detector_params[6],
            cfg.find_orientations.threshold,
            cfg.find_orientations.seed_search.hkl_seeds,
            cfg.find_orientations.seed_search.fiber_ndiv,
            ncpus=cfg.multiprocessing
            )
        if save_as_ascii:
            np.savetxt(
                os.path.join(cfg.working_dir, 'trial_orientations.dat'),
                quats.T,
                fmt="%.18e",
                delimiter="\t"
                )
            pass
        pass # close conditional on grid search

    # generate the completion maps
    logger.info("Running paintgrid on %d trial orientations", quats.shape[1])
    if profile:
        logger.info("Profiling mode active, forcing ncpus to 1")
        ncpus = 1
    else:
        ncpus = cfg.multiprocessing
        logger.info(
            "%d of %d available processors requested", ncpus, mp.cpu_count()
            )
    compl = idx.paintGrid(
        quats,
        eta_ome,
        etaRange=np.radians(cfg.find_orientations.eta.range),
        omeTol=np.radians(cfg.find_orientations.omega.tolerance),
        etaTol=np.radians(cfg.find_orientations.eta.tolerance),
        omePeriod=np.radians(cfg.find_orientations.omega.period),
        threshold=cfg.find_orientations.threshold,
        doMultiProc=ncpus > 1,
        nCPUs=ncpus
        )

    if save_as_ascii:
        np.savetxt(os.path.join(cfg.working_dir, 'completeness.dat'), compl)
    else:
        np.save(
            os.path.join(cfg.working_dir,
                         'scored_orientations_%s.npy' % analysis_id),
            np.vstack([quats, compl])
        )

    ##########################################################
    ##   Simulate N random grains to get neighborhood size  ##
    ##########################################################
    if hkl_ids is not None:
        ngrains = 100
        rand_q = mutil.unitVector(np.random.randn(4, ngrains))
        rand_e = np.tile(2.*np.arccos(rand_q[0, :]), (3, 1)) \
          * mutil.unitVector(rand_q[1:, :])
        refl_per_grain = np.zeros(ngrains)
        num_seed_refls = np.zeros(ngrains)
        for i in range(ngrains):
            grain_params = np.hstack([rand_e[:, i],
                                      xf.zeroVec.flatten(),
                                      xf.vInv_ref.flatten()
                                      ])

            eta_range = np.radians(cfg.find_orientations.eta.range)
            pixel_pitch = cfg.instrument.detector.pixels.size
            sim_results = xrdutil.simulateGVecs(
                pd,
                detector_params,
                grain_params,
                ome_range=(ome_range,),
                ome_period=(ome_range[0], ome_range[0]+2*np.pi),
                eta_range=eta_range,
                panel_dims=panel_dims,
                pixel_pitch=pixel_pitch,
                distortion=distortion,
            )
            refl_per_grain[i] = len(sim_results[0])
            num_seed_refls[i] = np.sum([sum(sim_results[0] == hkl_id) for hkl_id in hkl_ids])
            pass

        cfg_completeness = cfg.find_orientations.clustering.completeness
        min_samples = max(np.floor(cfg_completeness*np.average(num_seed_refls)), 2)
        mean_rpg = int(np.round(np.average(refl_per_grain)))
    else:
        min_samples = 1
        mean_rpg = 1

    logger.info("mean number of reflections per grain is %d", mean_rpg)
    logger.info("neighborhood size estimate is %d points", min_samples)

    # cluster analysis to identify orientation blobs, the final output:
    qbar, cl = run_cluster(compl, quats, pd.getQSym(), cfg, min_samples=min_samples)

    np.savetxt(
        os.path.join(cfg.working_dir, 'accepted_orientations_%s.dat' % analysis_id),
        qbar.T,
        fmt="%.18e",
        delimiter="\t"
        )
    return
