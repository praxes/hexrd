import cPickle
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
from hexrd.coreutil import initialize_experiment, make_eta_ranges, merge_dicts

from hexrd.xrd import xrdutil
from hexrd.xrd.xrdbase import multiprocessing
from hexrd.xrd.detector import ReadGE

from hexrd.xrd import distortion as dFuncs

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except:
    pass

# TODO: just require scikit-learn?
have_sklearn = False
try:
    import sklearn
    vstring = sklearn.__version__.split('.')
    if vstring[0] == '0' and int(vstring[1]) >= 14:
        from sklearn.cluster import dbscan
        from sklearn.metrics.pairwise import pairwise_distances
        have_sklearn = True
    else:
        print "Installed scikit-learn is too old (<0.14), using scipy fallback"
except:
    print "System does not have SciKit installed, using scipy fallback"


def generate_orientation_fibers(
    eta_ome, threshold, seed_hkl_ids, fiber_ndiv, verbose=False
    ):
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

    if verbose:
        print "labeling maps..."
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

    if verbose:
        print "generating quaternions..."

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


def paintgrid(pd, eta_ome, quats, threshold,
                  omeTol=None, etaTol=None,
                  omeRange=None, etaRange=None,
                  omePeriod=(-np.pi, np.pi),
                  qTol=1e-7,
                  ncpus=multiprocessing.cpu_count(),
                  verbose=False
                  ):
    """
    wrapper for indexer.paintGrid
    """

    # tolerances in degrees...  I know, pathological
    if omeTol is None:
        omeTol = 360. / float(fiber_ndiv)
    if etaTol is None:
        etaTol = 360. / float(fiber_ndiv)

    if verbose:
        print "Running paintgrid on %d trial orientations" % (quats.shape[1])
    complPG = idx.paintGrid(
        quats,
        eta_ome,
        omegaRange=omeRange, etaRange=etaRange,
        omeTol=np.radians(omeTol), etaTol=np.radians(etaTol),
        omePeriod=omePeriod, threshold=threshold,
        doMultiProc=ncpus>1,
        nCPUs=ncpus)
    return complPG


def run_cluster(compl, qfib, qsym, cfg, verbose=False):
    """
    """
    clcfg = cfg['find_orientations']['clustering']
    cl_radius = clcfg['radius']
    min_compl = clcfg['completeness']
    algorithm = clcfg.get('algorithm', 'dbscan')

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

        if verbose:
            print (
                "Feeding %d orientations above %.1f%% to clustering"
                % (qfib_r.shape[1], 100*min_compl)
                )

        if algorithm == 'dbscan' and not have_sklearn:
            algorithm = 'fclusterdata'
            if verbose:
                print "sklearn >= 0.14 required for dbscan"
                print "falling back to fclusterdata"
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
            qbar[:, i] = rot.quatAverage(qfib_r[:, cl == i + 1].reshape(4, npts),
                                         qsym).flatten()
        elapsed = (time.clock() - start)

    if verbose:
        print "clustering took %f seconds" % (elapsed)

    if verbose:
        print (
            "Found %d orientation clusters"
            " with >=%.1f%% completeness"
            " and %2f misorientation"
            % (qbar.size/4, 100.*min_compl, cl_radius)
            )

    return np.atleast_2d(qbar), cl


def load_eta_ome_maps(cfg, pd, reader, detector, verbose=False, hkls=None):
    cwd = cfg.get('working_dir', os.getcwd())
    fn = cfg['find_orientations']['orientation_maps'].get('file', None)
    try:
        fn = os.path.join(cwd, fn)
        res = cPickle.load(open(fn, 'r'))
        pd = res.planeData
        available_hkls = pd.hkls.T
        active_hkls = range(available_hkls.shape[0])
        if verbose:
            print "loaded eta/ome orientation maps from %s" % fn
            print "hkls used to generate orientation maps:"
            for i in available_hkls[active_hkls]:
                print i
        return res
    except (AttributeError, IOError):
        return generate_eta_ome_maps(cfg, pd, reader, detector, verbose, hkls)


def generate_eta_ome_maps(cfg, pd, reader, detector, verbose=False, hkls=None):
    omcfg = cfg['find_orientations']['orientation_maps']
    threshold = omcfg['threshold']
    bin_frames = omcfg.get('bin_frames', 1)

    available_hkls = pd.hkls.T
    # default to all hkls defined for material
    active_hkls = range(available_hkls.shape[0])
    # override with hkls from config, if specified
    active_hkls = omcfg.get('active_hkls', active_hkls)
    # override with hkls from command line, if specified
    active_hkls = hkls if hkls is not None else active_hkls

    if verbose:
        print "using hkls to generate orientation maps:"
        for i in available_hkls[active_hkls]:
            print i

    eta_bins = np.int(2*np.pi / abs(reader.getDeltaOmega())) / bin_frames
    if verbose:
        print "Using %d eta bins" % (eta_bins)
        print "loading data...",
    #import pdb; pdb.set_trace()
    eta_ome = xrdutil.CollapseOmeEta(
        reader,
        pd,
        active_hkls,
        detector,
        nframesLump=bin_frames,
        nEtaBins=eta_bins,
        debug=False,
        threshold=threshold
        )
    if verbose:
        print "done"

    cwd = cfg.get('working_dir', os.getcwd())
    outfile = omcfg.get('file', None)
    if outfile is not None:
        fn = os.path.join(cwd, outfile)
        fd = os.path.split(fn)[0]
        if not os.path.isdir(fd):
            os.makedirs(fd)
        with open(fn, 'w') as f:
            cPickle.dump(eta_ome, f)
        if verbose:
            print "saved eta/ome orientation maps to %s" % fn
    return eta_ome


def find_orientations(
    cfg, verbose=False, hkls=None, force=False
    ):
    """Takes a config dict as input, generally a yml document"""

    # a goofy call, could be replaced with two more targeted calls
    pd, reader, detector = initialize_experiment(cfg, verbose)

    cwd = cfg.get('working_dir', os.getcwd())
    analysis_name = cfg['analysis_name']
    analysis_root = os.path.join(cwd, analysis_name)
    if os.path.exists(analysis_root) and not force:
        print (
            'analysis "%s" already exists, change yml file or specify "force"'
            % analysis_name
            )
        sys.exit()
    if not os.path.exists(analysis_root):
        os.makedirs(analysis_root)
    if verbose:
        print "beginning analysis '%s'" % analysis_name

    # load the eta_ome orientation maps
    eta_ome = load_eta_ome_maps(cfg, pd, reader, detector, verbose, hkls)

    ############################################
    ## load parameters required by paint grid ##
    ############################################

    pgcfg = cfg['find_orientations']

    try:
        ome_tol = pgcfg['ome'].get('tolerance', None)
        ome_period = pgcfg['ome'].get('period', None)
    except (KeyError, AttributeError):
        ome_tol = None
        ome_period = None
    if ome_tol is None:
        ome_tol = abs(cfg['image_series']['ome']['step'])
    if ome_period is None:
        temp = cfg['image_series']['ome']['start']
        if cfg['image_series']['ome']['step'] > 0:
            ome_period = [temp, temp + 360]
        else:
            ome_period = [temp, temp - 360]
        if verbose:
            print "Omega tolerance: %g" % ome_tol
            print "Omega period: %s" % ome_period
    ome_period = np.radians(ome_period)
    try:
        eta_tol = pgcfg['eta'].get('tolerance', None)
        eta_mask = abs(pgcfg.get('mask', 5))
    except (KeyError, AttributeError):
        eta_tol = None
        eta_mask = 5
    if eta_tol is None:
        eta_tol = 2*ome_tol     # ome tol is half, eta is full
    if verbose:
        print "Eta tolerance: %g" % eta_tol
    eta_range = None
    if eta_mask:
        eta_range = make_eta_ranges(eta_mask)
        if verbose:
            print (
                "Masking eta angles within %g of ome rotation axis"
                % eta_mask
                )
    else:
        if verbose:
            print "Using full eta range"

    threshold = pgcfg.get('threshold', 1)

    # determine number of processes to run in parallel
    multiproc = pgcfg.get('multiprocessing', -1)
    ncpus = multiprocessing.cpu_count()
    if multiproc == 'all':
        pass
    elif multiproc == -1:
        ncpus -= 1
    elif int(ncpus) == 'half':
        ncpus /= 2
    elif isinstance(multiproc, int):
        if multiproc < ncpus:
            ncpus = multiproc
    else:
        ncpus -= 1
        if verbose:
            print (
                "Invalid value %s for find_orientations:multiprocessing"
                % multiproc
                )
    ncpus = ncpus if ncpus else 1

    # are we searching the full grid of orientation space, or a seeded search:
    seeded = False
    try:
        qgrid_f = pgcfg.get('use_quaternian_grid', None)
        quats = np.loadtxt(qgrid_f)
        if verbose:
            print "Using %s for full quaternian search" % qgrid_f
    except (IOError, ValueError):
        seeded = True
        if verbose:
            print "Could not open quaternian grid file %s" % qgrid_f
            print "Defaulting to seeded search"
        try:
            seed_hkl_ids = pgcfg['seed_search'].get('hkl_seeds', [0])
            fiber_step = pgcfg.get('fiber_step', None)
        except KeyError:
            raise RuntimeError(
                "if use_quaternian_grid is undefined, you must specify"
                " seed_search:hkl_seeds"
                )
        if fiber_step is None:
            fiber_step = ome_tol
        fiber_ndiv = int(360.0 / fiber_step)
        if verbose:
            print (
                "Seeding search using hkls from %s:"
                % pgcfg['orientation_maps']['file']
                )
            print eta_ome.planeData.hkls.T[seed_hkl_ids]
        quats = generate_orientation_fibers(
            eta_ome, threshold, seed_hkl_ids, fiber_ndiv, verbose
            )

    # run paintgrid
    compl = paintgrid(
        pd, eta_ome, quats, threshold,
        omeTol=ome_tol, etaTol=eta_tol, etaRange=eta_range,
        omePeriod=ome_period,
        qTol=1e-7,
        ncpus=ncpus,
        verbose=verbose
        )

    # cluster analysis to identify orientation blobs
    qbar, cl = run_cluster(compl, quats, pd.getQSym(), cfg, verbose)

    #################
    ## save output ##
    #################

    if seeded:
        # all of the orientations tested
        testq_f = os.path.join(analysis_root, 'testq.out')
        np.savetxt(testq_f, quats.T, fmt="%.18e", delimiter="\t")
    # raw completeness
    np.savetxt(os.path.join(analysis_root, 'compl.out'), compl)
    # main output, the list of quaternian orientation clusters
    # the result of cluster analysis on the thresholded completion map
    quats_f = os.path.join(analysis_root, 'quats.out')
    np.savetxt(quats_f, qbar.T, fmt="%.18e", delimiter="\t")

    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # from hexrd.xrd import rotations as rot
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # phis, ns = rot.angleAxisOfRotMat(rot.rotMatOfQuat(quats[:, np.r_[compl] > min_compl]))
    # rod = np.tile(np.tan(0.5*phis), (3, 1)) * ns
    # ax.scatter(rod[0, :], rod[1, :], rod[2, :], c='r', marker='o')
    #
    # phis, ns = rot.angleAxisOfRotMat(rot.rotMatOfQuat(qbar))
    # rod = np.tile(np.tan(0.5*phis), (3, 1)) * ns
    # ax.scatter(rod[0, :], rod[1, :], rod[2, :], c='b', marker='*')
    #
    # plt.show()

    ########################
    ## do extraction now? ##
    ########################

    if pgcfg.get('extract_measured_g_vectors', False):
        raise ImplementationError('TODO: implement extract gvecs')
        #extract_measured_g_vectors(cfg, verbose)
