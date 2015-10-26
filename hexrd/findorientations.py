
from __future__ import print_function

import cPickle
import logging
import multiprocessing as mp
import os
import shelve
import sys
import time

import numpy as np
#np.seterr(over='ignore', invalid='ignore')

import scipy.cluster as cluster
import scipy.optimize as opt
from scipy import ndimage

import yaml

from hexrd import matrixutil as mutil
from hexrd.xrd import experiment as expt
from hexrd.xrd import indexer as idx
from hexrd.xrd import rotations as rot
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.coreutil import initialize_experiment

from hexrd.xrd import xrdutil
from hexrd.xrd.detector import ReadGE
from hexrd.xrd.xrdutil import simulateGVecs

from hexrd.xrd import distortion as dFuncs

from hexrd.fitgrains import get_instrument_parameters

logger = logging.getLogger(__name__)


# TODO: just require scikit-learn?
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

have_parallel_dbscan = False
try:
    from parallel_dbscan import omp_dbscan
    have_parallel_dbscan = True
except ImportError:
    pass

def generate_orientation_fibers(eta_ome, threshold, seed_hkl_ids, fiber_ndiv):
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
    structureNDI_label = ndimage.generate_binary_structure(2, 2)

    # crystallography data from the pd object
    pd = eta_ome.planeData
    tTh  = pd.getTTh()
    bMat = pd.latVecOps['B']
    csym = pd.getLaueGroup()
    qsym = pd.getQSym()
    
    ############################################
    ##    Labeling of spots from seed hkls    ##
    ############################################

    qfib     = []
    labels   = []
    numSpots = []
    coms     = []
    for i in seed_hkl_ids:
        labels_t, numSpots_t = ndimage.label(
            eta_ome.dataStore[i] > threshold,
            structureNDI_label
            )
        coms_t = np.atleast_2d(
            ndimage.center_of_mass(
                eta_ome.dataStore[i],
                labels=labels_t,
                index=np.arange(1, np.amax(labels_t)+1)
                )
            )
        labels.append(labels_t)
        numSpots.append(numSpots_t)
        coms.append(coms_t)
        pass

    ############################################
    ##  Generate discrete fibers from labels  ##
    ############################################

    for i in range(len(pd_hkl_ids)):
        ii = 0
        qfib_tmp = np.empty((4, fiber_ndiv*numSpots[i]))
        for ispot in range(numSpots[i]):
            if not np.isnan(coms[i][ispot][0]):
                ome_c = eta_ome.omeEdges[0] \
                        + (0.5 + coms[i][ispot][0])*del_ome
                eta_c = eta_ome.etaEdges[0] \
                        + (0.5 + coms[i][ispot][1])*del_eta

                gVec_s = xrdutil.makeMeasuredScatteringVectors(
                    tTh[pd_hkl_ids[i]], eta_c, ome_c
                    )
                tmp = mutil.uniqueVectors(
                    rot.discreteFiber(
                        pd.hkls[:, pd_hkl_ids[i]].reshape(3, 1),
                        gVec_s,
                        B=bMat,
                        ndiv=fiber_ndiv,
                        invert=False,
                        csym=csym
                        )[0]
                    )
                jj = ii + tmp.shape[1]
                qfib_tmp[:, ii:jj] = tmp
                ii += tmp.shape[1]
                pass
            pass
        qfib.append(qfib_tmp[:, :ii])
        pass
    return np.hstack(qfib)


def cluster_parallel_dbscan(qfib_r, cl_radius, min_samples):
    homochoric_coords = xfcapi.homochoricOfQuat(qfib_r.T)
    labels = parallel_dbscan(
        homochoric_coords,
        eps=cl_radius,
        min_samples=min_samples)
    return labels + 1


def cluster_homochoric_dbscan(qfib_r, cl_radius, min_samples):
    homochoric_coords = xfcapi.homochoricOfQuat(qfib_r.T)
    _, labels = dbscan(
        homochoric_coords,
        eps=cl_radius,
        min_samples=min_samples)
    return labels + 1


def cluster_fcluster(qfib_r, cl_radius, min_samples):
    cl = cluster.hierarchy.fclusterdata(
        qfib_r.T,
        np.radians(cl_radius),
        criterion='distance',
        metric=xfcapi.quat_distance
    )
    return cl


def cluster_dbscan(qfib_r, cl_radius, min_samples):
    def quat_distance(x, y):
        qsym  = np.array(qsym.T, order='C').T
        return xfcapi.quat_distance(np.array(x, order='C'), np.array(y, order='C'), qsym)

    pdist = pairwise_distances(
        qfib_r.T, metric=xfcapi.quat_distance, n_jobs=-1
    )
    _, labels = dbscan(
        pdist,
        eps=np.radians(cl_radius),
        min_samples=min_samples,
        metric='precomputed'
    )
    cl = np.array(labels, dtype=int) + 1

def run_cluster(compl, qfib, qsym, cfg, min_samples=None):
    """
    """
    cl_radius = cfg.find_orientations.clustering.radius
    min_compl = cfg.find_orientations.clustering.completeness
    algorithm = cfg.find_orientations.clustering.algorithm

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

        qfib_r = qfib[:, np.array(compl) > min_compl].ascontiguousarray()
        logger.info(
            "Feeding %d orientations above %.1f%% to clustering",
            qfib_r.shape[1], 100*min_compl
            )

        if algorithm == 'parallel-dbscan' and not have_parallel_dbscan:
            algorithm = 'homochoric-dbscan'
            logger.warning(
                "parallel_dbscan not found, trying sklearn dbscan"
                )

        if algorithm in ('dbscan', 'homochoric-dbscan') and not have_sklearn:
            algorithm = 'fclusterdata'
            logger.warning(
                "sklearn >= 0.14 required for dbscan, using fclusterdata"
                )

        if algorithm == 'parallel-dbscan':
            cl = cluster_parallel_dbscan(qfib_r, cl_radius, min_samples)
        elif algorithm == 'homochoric-dbscan':
            cl = cluster_homochoric_dbscan(qfib_r, cl_radius, min_samples)
        elif algorithm == 'dbscan':
            cl = cluster_dbscan(qfib_r, cl_radius, min_samples)
        elif algorithm == 'fclusterdata':
            cl = cluster_fclusterdata(qfib_r, cl_radius, min_samples)
        else:
            raise RuntimeError(
                "Clustering algorithm %s not recognized" % algorithm
                )

        nblobs = len(np.unique(cl))

        qbar = np.zeros((4, nblobs))
        for i in range(nblobs):
            npts = sum(cl == i + 1)
            qbar[:, i] = rot.quatAverage(
                qfib_r[:, cl == i + 1].reshape(4, npts), qsym
                ).flatten()

    logger.info("clustering took %f seconds", time.clock() - start)
    logger.info(
        "Found %d orientation clusters with >=%.1f%% completeness"
        " and %2f misorientation",
        qbar.size/4,
        100.*min_compl,
        cl_radius
        )

    return np.atleast_2d(qbar), cl


def load_eta_ome_maps(cfg, pd, reader, detector, hkls=None):
    fn = os.path.join(
        cfg.working_dir,
        cfg.find_orientations.orientation_maps.file
        )
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


def find_orientations(cfg, hkls=None, profile=False):
    """
    Takes a config dict as input, generally a yml document
    
    NOTE: single cfg instance, not iterator!
    """

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
    eta_ome = load_eta_ome_maps(cfg, pd, reader, detector, hkls)

    ome_range = (np.min(eta_ome.omeEdges),
                 np.max(eta_ome.omeEdges)
                 )
    try:
        # are we searching the full grid of orientation space?
        qgrid_f = cfg.find_orientations.use_quaternion_grid
        quats = np.loadtxt(qgrid_f).T
        logger.info("Using %s for full quaternion search", qgrid_f)
        hkl_ids = None
    except (IOError, ValueError):
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
            cfg.find_orientations.threshold,
            cfg.find_orientations.seed_search.hkl_seeds,
            cfg.find_orientations.seed_search.fiber_ndiv
            )
        np.savetxt(
            os.path.join(cfg.working_dir, 'trial_orientations.dat'),
            quats.T,
            fmt="%.18e",
            delimiter="\t"
            )

    # generate the completion maps
    logger.info("Running paintgrid on %d trial orientations", (quats.shape[1]))
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
    np.savetxt(os.path.join(cfg.working_dir, 'completeness.dat'), compl)

    ##########################################################
    ##   Simulate N random grains to get neighborhood size  ##
    ##########################################################
    if hkl_ids is not None:
        ngrains = 10
        rand_q = mutil.unitVector(np.random.randn(4, ngrains))
        rand_e = np.tile(2.*np.arccos(rand_q[0, :]), (3, 1)) \
          * mutil.unitVector(rand_q[1:, :])
        num_seed_refls = np.zeros(ngrains)
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
            num_seed_refls[i] = np.sum([sum(sim_results[0] == hkl_id) for hkl_id in hkl_ids])
            pass
        min_samples = max(
            cfg.find_orientations.clustering.completeness*np.floor(np.average(num_seed_refls)),
            4
            )
    else:
        min_samples = 1

    logger.info("neighborhood size estimate is %d points", (min_samples))
    
    # cluster analysis to identify orientation blobs, the final output:
    qbar, cl = run_cluster(compl, quats, pd.getQSym(), cfg, min_samples=min_samples)
    np.savetxt(
        os.path.join(cfg.working_dir, 'accepted_orientations.dat'),
        qbar.T,
        fmt="%.18e",
        delimiter="\t"
        )
