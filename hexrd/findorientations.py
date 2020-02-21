from __future__ import print_function

import cPickle
import logging
import multiprocessing as mp
import os
import time

import yaml

import numpy as np
# np.seterr(over='ignore', invalid='ignore')

import scipy.cluster as cluster
from scipy import ndimage

from hexrd import matrixutil as mutil
from hexrd.xrd import indexer as idx
from hexrd.xrd import rotations as rot
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi

from hexrd.xrd.xrdutil import GenerateEtaOmeMaps, EtaOmeMaps, simulateGVecs

from hexrd.xrd import distortion as dFuncs

from hexrd.fitgrains import get_instrument_parameters

from skimage.feature import blob_dog, blob_log

# just require scikit-learn?
have_sklearn = False
try:
    import sklearn
    vstring = sklearn.__version__.split('.')
    if vstring[0] == '0' and int(vstring[1]) >= 14:
        from sklearn.cluster import dbscan
        from sklearn.metrics.pairwise import pairwise_distances
        have_sklearn = True
except ImportError:
    pass


method = "blob_dog"  # !!! have to get this from the config
save_as_ascii = False  # FIX LATER...

logger = logging.getLogger(__name__)


# =============================================================================
# FUNCTIONS
# =============================================================================


def generate_orientation_fibers(
        eta_ome, chi, threshold, seed_hkl_ids, fiber_ndiv,
        method='blob_dog', filt_stdev=0.8, ncpus=1):
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
    tTh = pd.getTTh()
    bMat = pd.latVecOps['B']
    csym = pd.getLaueGroup()

    params = dict(
        bMat=bMat,
        chi=chi,
        csym=csym,
        fiber_ndiv=fiber_ndiv)

    # =========================================================================
    # Labeling of spots from seed hkls
    # =========================================================================

    qfib = []
    input_p = []
    numSpots = []
    coms = []
    for i in seed_hkl_ids:
        if method == 'label':
            # First apply filter
            this_map_f = -ndimage.filters.gaussian_laplace(
                eta_ome.dataStore[i], filt_stdev)
            
            labels_t, numSpots_t = ndimage.label(
                this_map_f > threshold,
                structureNDI_label
                )
            coms_t = np.atleast_2d(
                ndimage.center_of_mass(
                    this_map_f,
                    labels=labels_t,
                    index=np.arange(1, np.amax(labels_t) + 1)
                    )
                )
        elif method in ['blob_log', 'blob_dog']:
            # must scale map
            this_map = eta_ome.dataStore[i]
            this_map[np.isnan(this_map)] = 0.
            this_map -= np.min(this_map)
            scl_map = 2*this_map/np.max(this_map) - 1.

            # FIXME: need to expose the parameters to config options.
            if method == 'blob_log':
                blobs = np.atleast_2d(
                    blob_log(scl_map, min_sigma=0.5, max_sigma=5,
                             num_sigma=10, threshold=0.01, overlap=0.1)
                )
            else:
                blobs = np.atleast_2d(
                    blob_dog(scl_map, min_sigma=0.5, max_sigma=5,
                             sigma_ratio=1.6, threshold=0.01, overlap=0.1)
                    )
            numSpots_t = len(blobs)
            coms_t = blobs[:, :2]
        numSpots.append(numSpots_t)
        coms.append(coms_t)
        pass

    for i in range(len(pd_hkl_ids)):
        for ispot in range(numSpots[i]):
            if not np.isnan(coms[i][ispot][0]):
                ome_c = eta_ome.omeEdges[0] + (0.5 + coms[i][ispot][0])*del_ome
                eta_c = eta_ome.etaEdges[0] + (0.5 + coms[i][ispot][1])*del_eta
                input_p.append(
                    np.hstack(
                        [hkls[:, pd_hkl_ids[i]],
                         tTh[pd_hkl_ids[i]], eta_c, ome_c]
                    )
                )
                pass
            pass
        pass

    # do the mapping
    start = time.time()
    qfib = None
    if ncpus > 1:
        # multiple process version
        # QUESTION: Need a chunksize?
        pool = mp.Pool(ncpus, discretefiber_init, (params, ))
        qfib = pool.map(discretefiber_reduced, input_p)  # chunksize=chunksize)
        pool.close()
    else:
        # single process version.
        global paramMP
        discretefiber_init(params)  # sets paramMP
        qfib = map(discretefiber_reduced, input_p)
        paramMP = None  # clear paramMP
    elapsed = (time.time() - start)
    logger.info("fiber generation took %.3f seconds", elapsed)
    return np.hstack(qfib)


def discretefiber_init(params):
    global paramMP
    paramMP = params


def discretefiber_reduced(params_in):
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
            csym=csym
            )[0]
        )
    return tmp


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

    start = time.clock() # time this

    num_above = sum(np.array(compl) > min_compl)
    if num_above == 0:
        # nothing to cluster
        qbar = cl = np.array([])
    elif num_above == 1:
        # short circuit
        qbar = qfib[:, np.array(compl) > min_compl]
        cl = [1]
    else:
        # use compiled module for distance
        # just to be safe, must order qsym as C-contiguous
        qsym  = np.array(qsym.T, order='C').T
        def quat_distance(x, y):
            return xfcapi.quat_distance(np.array(x, order='C'), np.array(y, order='C'), qsym)

        qfib_r = qfib[:, np.array(compl) > min_compl]

        num_ors = qfib_r.shape[1]

        if num_ors > 25000:
            if algorithm == 'sph-dbscan' or algorithm == 'fclusterdata':
                logger.info("falling back to euclidean DBSCAN")
                algorithm = 'ort-dbscan'
            #raise RuntimeError, \
            #    "Requested clustering of %d orientations, which would be too slow!" %qfib_r.shape[1]

        logger.info(
            "Feeding %d orientations above %.1f%% to clustering",
            num_ors, 100*min_compl
            )

        if algorithm == 'dbscan' and not have_sklearn:
            algorithm = 'fclusterdata'
            logger.warning(
                "sklearn >= 0.14 required for dbscan; using fclusterdata"
                )

        if algorithm == 'dbscan' or algorithm == 'ort-dbscan' or algorithm == 'sph-dbscan':
            # munge min_samples according to options
            if min_samples is None or cfg.find_orientations.use_quaternion_grid is not None:
                min_samples = 1

            if algorithm == 'sph-dbscan':
                logger.info("using spherical DBSCAN")
                # compute distance matrix
                pdist = pairwise_distances(
                    qfib_r.T, metric=quat_distance, n_jobs=1
                    )

                # run dbscan
                core_samples, labels = dbscan(
                    pdist,
                    eps=np.radians(cl_radius),
                    min_samples=min_samples,
                    metric='precomputed'
                    )
            else:
                if algorithm == 'ort-dbscan':
                    logger.info("using euclidean orthographic DBSCAN")
                    pts = qfib_r[1:, :].T
                    eps = 0.25*np.radians(cl_radius)
                else:
                    logger.info("using euclidean DBSCAN")
                    pts = qfib_r.T
                    eps = 0.5*np.radians(cl_radius)

                # run dbscan
                core_samples, labels = dbscan(
                    pts,
                    eps=eps,
                    min_samples=min_samples,
                    metric='minkowski', p=2,
                    )

            # extract cluster labels
            cl = np.array(labels, dtype=int) # convert to array
            noise_points = cl == -1 # index for marking noise
            cl += 1 # move index to 1-based instead of 0
            cl[noise_points] = -1 # re-mark noise as -1
            logger.info("dbscan found %d noise points", sum(noise_points))
        elif algorithm == 'fclusterdata':
            logger.info("using spherical fclusetrdata")
            cl = cluster.hierarchy.fclusterdata(
                qfib_r.T,
                np.radians(cl_radius),
                criterion='distance',
                metric=quat_distance
                )
        else:
            raise RuntimeError(
                "Clustering algorithm %s not recognized" % algorithm
                )

        # extract number of clusters
        if np.any(cl == -1):
            nblobs = len(np.unique(cl)) - 1
        else:
            nblobs = len(np.unique(cl))

        """ PERFORM AVERAGING TO GET CLUSTER CENTROIDS """
        qbar = np.zeros((4, nblobs))
        for i in range(nblobs):
            npts = sum(cl == i + 1)
            qbar[:, i] = rot.quatAverageCluster(
                qfib_r[:, cl == i + 1].reshape(4, npts), qsym
            ).flatten()
            pass
        pass

    if (algorithm == 'dbscan' or algorithm == 'ort-dbscan') \
      and qbar.size/4 > 1:
        logger.info("\tchecking for duplicate orientations...")
        cl = cluster.hierarchy.fclusterdata(
            qbar.T,
            np.radians(cl_radius),
            criterion='distance',
            metric=quat_distance)
        nblobs_new = len(np.unique(cl))
        if nblobs_new < nblobs:
            logger.info("\tfound %d duplicates within %f degrees" \
                        %(nblobs-nblobs_new, cl_radius))
            tmp = np.zeros((4, nblobs_new))
            for i in range(nblobs_new):
                npts = sum(cl == i + 1)
                tmp[:, i] = rot.quatAverageCluster(
                    qbar[:, cl == i + 1].reshape(4, npts), qsym
                ).flatten()
                pass
            qbar = tmp
            pass
        pass

    logger.info("clustering took %f seconds", time.clock() - start)
    logger.info(
        "Found %d orientation clusters with >=%.1f%% completeness"
        " and %2f misorientation",
        qbar.size/4,
        100.*min_compl,
        cl_radius
        )

    return np.atleast_2d(qbar), cl


def load_eta_ome_maps(cfg, pd, image_series, hkls=None, clean=False):
    fn = os.path.join(
        cfg.working_dir,
        cfg.find_orientations.orientation_maps.file
        )

    # ...necessary?
    if fn.split('.')[-1] != 'npz':
        fn = fn + '.npz'

    if not clean:
        try:
            res = EtaOmeMaps(fn)
            pd = res.planeData
            available_hkls = pd.hkls.T
            logger.info('loaded eta/ome orientation maps from %s', fn)
            hkls = [str(i) for i in available_hkls[res.iHKLList]]
            logger.info(
                'hkls used to generate orientation maps: %s', hkls)
            return res
        except (AttributeError, IOError):
            return generate_eta_ome_maps(cfg, pd, image_series, hkls)
    else:
        logger.info('clean option specified; recomputing eta/ome orientation maps')
        return generate_eta_ome_maps(cfg, pd, image_series, hkls)

def generate_eta_ome_maps(cfg, pd, image_series, hkls=None):

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

    bin_frames = cfg.find_orientations.orientation_maps.bin_frames
    ome_step = cfg.image_series.omega.step*bin_frames
    instrument_params = yaml.load(open(cfg.instrument.parameters, 'r'))

    # generate maps
    eta_ome = GenerateEtaOmeMaps(
        image_series, instrument_params, pd, active_hkls,
        ome_step=ome_step,
        threshold=cfg.find_orientations.orientation_maps.threshold
        )

    fn = os.path.join(
        cfg.working_dir,
        cfg.find_orientations.orientation_maps.file
        )
    fd = os.path.split(fn)[0]
    if not os.path.isdir(fd):
        os.makedirs(fd)
    eta_ome.save(fn)
    logger.info("saved eta/ome orientation maps to %s", fn)
    return eta_ome


def find_orientations(cfg, hkls=None, clean=False, profile=False):
    """
    Takes a config dict as input, generally a yml document

    NOTE: single cfg instance, not iterator!
    """

    # ...make this an attribute in cfg?
    analysis_id = '%s_%s' %(
        cfg.analysis_name.strip().replace(' ', '-'),
        cfg.material.active.strip().replace(' ', '-'),
        )

    # grab planeData object
    matl = cPickle.load(open('materials.cpl', 'r'))
    md = dict(zip([matl[i].name for i in range(len(matl))], matl))
    pd = md[cfg.material.active].planeData

    # make image_series
    image_series = cfg.image_series.omegaseries

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

    min_compl = cfg.find_orientations.clustering.completeness

    # start logger
    logger.info("beginning analysis '%s'", cfg.analysis_name)

    # load the eta_ome orientation maps
    eta_ome = load_eta_ome_maps(cfg, pd, image_series, hkls=hkls, clean=clean)

    ome_range = (
        np.min(eta_ome.omeEdges),
        np.max(eta_ome.omeEdges)
        )
    try:
        # are we searching the full grid of orientation space?
        qgrid_f = cfg.find_orientations.use_quaternion_grid
        quats = np.load(qgrid_f)
        logger.info("Using %s for full quaternion search", qgrid_f)
        hkl_ids = None
    except (IOError, ValueError, AttributeError):
        # or doing a seeded search?
        logger.info("Defaulting to seeded search")
        hkl_seeds = cfg.find_orientations.seed_search.hkl_seeds
        hkl_ids = [
            eta_ome.planeData.hklDataList[i]['hklID'] for i in hkl_seeds
        ]
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
            ncpus=cfg.multiprocessing,
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
            os.path.join(
                cfg.working_dir,
                'scored_orientations_%s.npy' %analysis_id
                ),
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
        print('fo: hklids = ', hkl_ids)
        for i in range(ngrains):
            grain_params = np.hstack([rand_e[:, i],
                                      xf.zeroVec.flatten(),
                                      xf.vInv_ref.flatten()
                                      ])
            sim_results = simulateGVecs(pd,
                                        detector_params,
                                        grain_params,
                                        ome_range=(ome_range,),
                                        ome_period=(ome_range[0], ome_range[0]+2*np.pi),
                                        eta_range=np.radians(cfg.find_orientations.eta.range),
                                        panel_dims=panel_dims,
                                        pixel_pitch=cfg.instrument.detector.pixels.size,
                                        distortion=distortion,
                                        )
            refl_per_grain[i] = len(sim_results[0])
            # lines below fix bug when sim_results[0] is empty
            if refl_per_grain[i] > 0:
                num_seed_refls[i] = np.sum([sum(sim_results[0] == hkl_id) for hkl_id in hkl_ids])
            else:
                num_seed_refls[i] = 0
        #min_samples = 2
        min_samples = max(
            int(np.floor(0.5*min_compl*min(num_seed_refls))),
            2
            )
        mean_rpg = int(np.round(np.average(refl_per_grain)))
    else:
        min_samples = 1
        mean_rpg = 1

    logger.info("mean number of reflections per grain is %d", mean_rpg)
    logger.info("neighborhood size estimate is %d points", min_samples)

    # cluster analysis to identify orientation blobs, the final output:
    qbar, cl = run_cluster(compl, quats, pd.getQSym(), cfg, min_samples=min_samples)

    analysis_id = '%s_%s' %(
        cfg.analysis_name.strip().replace(' ', '-'),
        cfg.material.active.strip().replace(' ', '-'),
        )

    np.savetxt(
        os.path.join(
            cfg.working_dir,
            'accepted_orientations_%s.dat' %analysis_id
            ),
        qbar.T,
        fmt="%.18e",
        delimiter="\t")

    return
