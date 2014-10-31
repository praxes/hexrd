from __future__ import print_function

import cPickle
import logging
import multiprocessing as mp
import os
import shelve
import sys
import time

import numpy as np
np.seterr(over='ignore', invalid='ignore')

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

from hexrd.xrd import distortion as dFuncs

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except:
    pass


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


def generate_orientation_fibers(eta_ome, threshold, seed_hkl_ids, fiber_ndiv):
    """ From ome-eta maps and hklid spec, generate list of
    quaternions from fibers

    ** passing of BOTH pd and eta_ome object is redundant; flag for fix!
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

    ii       = 0
    jj       = fiber_ndiv
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

    qfib_tmp = np.empty((4, fiber_ndiv*sum(numSpots)))

    for i in range(len(pd_hkl_ids)):
        for ispot in range(numSpots[i]):
            if not np.isnan(coms[i][ispot][0]):
                ome_c = eta_ome.omeEdges[0] \
                        + (0.5 + coms[i][ispot][0])*del_ome
                eta_c = eta_ome.etaEdges[0] \
                        + (0.5 + coms[i][ispot][1])*del_eta

                gVec_s = xrdutil.makeMeasuredScatteringVectors(
                    tTh[pd_hkl_ids[i]], eta_c, ome_c
                    )

                qfib_tmp[:, ii:jj] = rot.discreteFiber(
                    pd.hkls[:, pd_hkl_ids[i]].reshape(3, 1),
                    gVec_s,
                    B=bMat,
                    ndiv=fiber_ndiv,
                    invert=False,
                    csym=csym
                    )[0]
                ii  = jj
                jj += fiber_ndiv
                pass
            pass
        qfib.append(mutil.uniqueVectors(qfib_tmp))
        pass
    return np.hstack(qfib)


def run_cluster(compl, qfib, qsym, cfg):
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
        qsym  = np.array(qsym.T, order='C').T
        quat_distance = lambda x, y: xfcapi.quat_distance(
            np.array(x, order='C'),
            np.array(y, order='C'),
            qsym
            )

        qfib_r = qfib[:, np.array(compl) > min_compl]

        logger.info(
            "Feeding %d orientations above %.1f%% to clustering",
            qfib_r.shape[1], 100*min_compl
            )

        if algorithm == 'dbscan' and not have_sklearn:
            algorithm = 'fclusterdata'
            logger.warning(
                "sklearn >= 0.14 required for dbscan, using fclusterdata"
                )
        if algorithm == 'dbscan':
            pdist = pairwise_distances(
                qfib_r.T, metric=quat_distance, n_jobs=-1
                )
            core_samples, labels = dbscan(
                pdist,
                eps=np.radians(cl_radius),
                min_samples=1,
                metric='precomputed'
                )
            cl = np.array(labels, dtype=int) + 1
        elif algorithm == 'fclusterdata':
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

    bin_frames = cfg.find_orientations.orientation_maps.bin_frames
    eta_bins = np.int(2*np.pi / abs(reader.getDeltaOmega())) / bin_frames
    eta_ome = xrdutil.CollapseOmeEta(
        reader,
        pd,
        active_hkls,
        detector,
        nframesLump=bin_frames,
        nEtaBins=eta_bins,
        debug=False,
        threshold=cfg.find_orientations.orientation_maps.threshold
        )

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


def find_orientations(
    cfg, hkls=None, force=False
    ):
    """Takes a config dict as input, generally a yml document"""

    # a goofy call, could be replaced with two more targeted calls
    pd, reader, detector = initialize_experiment(cfg)

    if os.path.exists(cfg.analysis_dir) and not force:
        logger.error(
            'analysis "%s" already exists, change yml file or specify "force"',
            cfg.analysis_name
            )
        sys.exit()
    if not os.path.exists(cfg.analysis_dir):
        os.makedirs(cfg.analysis_dir)
    logger.info("beginning analysis '%s'", cfg.analysis_name)

    # load the eta_ome orientation maps
    eta_ome = load_eta_ome_maps(cfg, pd, reader, detector, hkls)

    try:
        # are we searching the full grid of orientation space?
        qgrid_f = cfg.find_orientations.use_quaternion_grid
        quats = np.loadtxt(qgrid_f)
        logger.info("Using %s for full quaternian search", qgrid_f)
    except (IOError, ValueError):
        # or doing a seeded search?
        logger.info("Defaulting to seeded search")
        hkl_seeds = cfg.find_orientations.seed_search.hkl_seeds
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
            os.path.join(cfg.analysis_dir, 'testq.out'),
            quats.T,
            fmt="%.18e",
            delimiter="\t"
            )

    # generate the completion maps
    logger.info("Running paintgrid on %d trial orientations", (quats.shape[1]))
    ncpus = cfg.multiprocessing
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
    np.savetxt(os.path.join(cfg.analysis_dir, 'compl.out'), compl)

    # cluster analysis to identify orientation blobs, the final output:
    qbar, cl = run_cluster(compl, quats, pd.getQSym(), cfg)
    np.savetxt(
        os.path.join(cfg.analysis_dir, 'quats.out'),
        qbar.T,
        fmt="%.18e",
        delimiter="\t"
        )

    # do the peak extraction now?
    if cfg.find_orientations.extract_measured_g_vectors:
        raise ImplementationError('TODO: implement extract gvecs')
        #extract_measured_g_vectors(cfg)
