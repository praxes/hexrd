import sys, os, time

import shelve, cPickle
import numpy as np
np.seterr(over='ignore', invalid='ignore')

import scipy.cluster  as cluster
import scipy.optimize as opt
from scipy import ndimage

from ConfigParser import SafeConfigParser

# check PYTHONPATH
try:
    import hexrd
except:
    cfgfile = 'default.cfg'
    parser = SafeConfigParser()
    parser.read(cfgfile)

    HEXRD_ROOT = parser.get('base', 'hexrd_root')
    sys.path.append(HEXRD_ROOT)

# now can import hexrd modules
from hexrd     import matrixutil      as mutil
from hexrd.xrd import experiment      as expt
from hexrd.xrd import indexer         as idx
from hexrd.xrd import rotations       as rot
from hexrd.xrd import transforms      as xf
from hexrd.xrd import transforms_CAPI as xfcapi

from hexrd.xrd          import xrdutil
from hexrd.xrd.xrdbase  import multiprocessing
from hexrd.xrd.detector import ReadGE

from hexrd.xrd import distortion as dFuncs

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except:
    pass

haveScikit = False
try:
    import sklearn
    vstring = sklearn.__version__.split('.') 
    if vstring[0] == '0' and int(vstring[1]) >= 15:
        from sklearn.cluster          import dbscan
        from sklearn.metrics.pairwise import pairwise_distances
        haveScikit = True
    else:
        print "Installed SciKit is too old (<0.14), using scipy fallback"        
except:
    print "System does not have SciKit installed, using scipy fallback"

"""
##################### BEGIN COMPUTATION HERE #####################
"""

# constants
d2r = np.pi/180.
r2d = 180./np.pi

bVec_ref = xf.bVec_ref

min_compl = 0.5
cl_radius = 1.0

# #################################################
# LEGACY FUCTIONS FOR WORKING WITH OLD HEXRD PAR
def tVec_d_from_old_parfile(old_par, detOrigin):
    beamXYD = old_par[:3, 0]
    rMat_d  = xf.makeDetectorRotMat(old_par[3:6, 0])
    args=(rMat_d, beamXYD, detOrigin, bVec_ref)
    tvd_xy = opt.leastsq(objFun_tVec_d, -beamXYD[:2], args=args)[0]
    return np.hstack([tvd_xy, -beamXYD[2]]).reshape(3, 1)

def objFun_tVec_d(tvd_xy, rMat_d, beamXYD, detOrigin, bHat_l):
    """
    """
    xformed_xy = beamXYD[:2] - detOrigin
    tVec_d = np.hstack([tvd_xy, -beamXYD[2]]).T
    n_d    = rMat_d[:, 2]

    bVec_l = (np.dot(n_d, tVec_d) / np.dot(n_d, bHat_l)) * bHat_l
    bVec_d = np.hstack([xformed_xy, 0.]).T

    return np.dot(rMat_d, bVec_d).flatten() + tVec_d.flatten() - bVec_l.flatten()

def beamXYD_from_tVec_d(rMat_d, tVec_d, bVec_ref, detOrigin):
    # calculate beam position
    n_d = np.dot(rMat_d, np.c_[0., 0., 1.].T)

    u = np.dot(n_d.flatten(), tVec_d) / np.dot(n_d.flatten(), bVec_ref)

    detOrigin - tVec_d[:2].flatten()
    return np.hstack([detOrigin - tVec_d[:2].flatten(), u.flatten()])

def write_old_parfile(filename, results):
    if isinstance(filename, file):
        fid = filename
    elif isinstance(filename, str) or isinstance(filename, unicode):
        fid = open(filename, 'w')
        pass
    rMat_d = xf.makeDetectorRotMat(results['tiltAngles'])
    tVec_d = results['tVec_d'] - results['tVec_s']
    beamXYD = beamXYD_from_tVec_d(rMat_d, tVec_d, bVec_ref, detOrigin)
    det_plist = np.zeros(12)
    det_plist[:3]  = beamXYD.flatten()
    det_plist[3:6] = results['tiltAngles']
    det_plist[6:]  = results['dParams']
    print >> fid, "# DETECTOR PARAMETERS (from new geometry model fit)"
    print >> fid, "# \n# <class 'hexrd.xrd.detector.DetectorGeomGE'>\n#"
    for i in range(len(det_plist)):
        print >> fid, "%1.8e\t%d" % (det_plist[i], 0)
    fid.close()
    return
# #################################################

def initialize_experiment(cfg_file):
    """
    """
    parser = SafeConfigParser()
    parser.read(cfg_file)

    hexrd_root = parser.get('base', 'hexrd_root')

    # make experiment
    ws = expt.Experiment(cfgFile=os.path.join(hexrd_root, "hexrd/data/materials.cfg"),
                         matFile=os.path.join(hexrd_root, "hexrd/data/all_materials.cfg"))

    working_dir = parser.get('base', 'working_dir')

    materials_fname = parser.get('material', 'materials_fname')
    material_name   = parser.get('material', 'material_name')
    detector_fname  = parser.get('detector', 'parfile_name')

    # MATERIALS
    ws.loadMaterialList(os.path.join(working_dir, materials_fname))
    ws.activeMaterial = material_name
    print "setting active material to '%s'" % (material_name)

    pd = ws.activeMaterial.planeData

    image_dir = parser.get('reader', 'image_dir')

    # number of files ASSUMING SEQUENTIAL SCAN NUMBERS
    file_start  = parser.getint('reader', 'file_start')
    file_stop   = parser.getint('reader', 'file_stop')
    file_suffix = parser.get('reader', 'file_suffix')
    nfiles      = file_stop - file_start + 1
    zpad_str    = '%0' + parser.get('reader', 'file_zpad') + 'd' # int
    fileInfo = []
    for i in [file_start + i for i in range(nfiles)]:
        if file_suffix == '':
            image_filename = parser.get('reader', 'file_root') + '_' + zpad_str % (i)
        else:
            image_filename = parser.get('reader', 'file_root') + '_' + zpad_str % (i) + '.' + file_suffix
        fileInfo.append( ( os.path.join(image_dir, image_filename), parser.getint('reader', 'nempty') ) )
    ome_start = parser.getfloat('reader', 'ome_start') * d2r
    ome_delta = parser.getfloat('reader', 'ome_delta') * d2r
    darkName  = parser.get('reader', 'dark')
    if darkName.strip() == '':
        dark         = None
        subtractDark = False
    else:
        # dark = os.path.join(image_dir, darkName)
        dark = darkName
        subtractDark = True
    doFlip  = parser.getboolean('reader', 'doFlip')
    flipArg = parser.get('reader', 'flipArg')

    # make frame reader
    reader   = ReadGE(fileInfo, ome_start, ome_delta,
                      subtractDark=subtractDark, dark=dark,
                      doFlip=doFlip, flipArg=flipArg)

    # DETECTOR
    ws.loadDetector(os.path.join(working_dir, detector_fname))

    return pd, reader, ws.detector

def make_maps(pd, reader, detector, hkl_ids, threshold, nframesLump, output=None):
    """
    OME-ETA MAPS
    """
    nEtaBins = np.int(2*np.pi / abs(reader.getDeltaOmega())) / nframesLump
    print "Using %d eta bins" % (nEtaBins)
    omeEta   = xrdutil.CollapseOmeEta(reader,
                                      pd, hkl_ids, detector,
                                      nframesLump=nframesLump, nEtaBins=nEtaBins,
                                      debug=False, threshold=threshold)
    if output is not None:
        if isinstance(output, str):
            fid = open(output, 'w')
            cPickle.dump(omeEta, fid)
            fid.close()
        elif isinstance(output, file):
            cPickle.dump(omeEta, output)
        else:
            raise RuntimeError, "output specifier must be a string or file"
    return omeEta

def run_paintGrid(pd, omeEta, seed_hkl_ids, threshold, fiber_ndiv,
                  omeTol=None, etaTol=None,
                  omeRange=None, etaRange=None,
                  omePeriod=(-np.pi, np.pi),
                  qTol=1e-7,
                  doMultiProc=True, nCPUs=multiprocessing.cpu_count(),
                  useGrid=None):
    """
    wrapper for indexer.paintGrid
    """
    del_ome = omeEta.omegas[1] - omeEta.omegas[0]
    del_eta = omeEta.etas[1] - omeEta.etas[0]

    # tolerances in degrees...  I know, pathological
    if omeTol is None:
        omeTol = 360. / float(fiber_ndiv)
    if etaTol is None:
        etaTol = 360. / float(fiber_ndiv)

    # must be consistent
    pd_hkl_ids = omeEta.iHKLList[seed_hkl_ids]

    tTh  = pd.getTTh()
    bMat = pd.latVecOps['B']
    csym = pd.getLaueGroup()
    qsym = pd.getQSym()

    if useGrid is not None:
        try:
            print "loading quaternion grid file: %s" % (useGrid)
            qfib = np.loadtxt(useGrid).T
        except:
            raise RuntimeError, "unable to load quaternion grid file"
    else:
        structureNDI_label = np.array([
                [0,1,0],
                [1,1,1],
                [0,1,0]
                ])
        qfib = []
        ii = 0
        jj = fiber_ndiv
        print "labeling maps..."
        labels   = []
        numSpots = []
        coms     = []
        for i in seed_hkl_ids:
            labels_t, numSpots_t = ndimage.label(omeEta.dataStore[i] > threshold, structureNDI_label)
            coms_t = np.atleast_2d(ndimage.center_of_mass(omeEta.dataStore[i],
                                                          labels=labels_t,
                                                          index=np.arange(1, np.amax(labels_t)+1)))
            labels.append(labels_t)
            numSpots.append(numSpots_t)
            coms.append(coms_t)
            pass

        # second pass for generation
        print "generating quaternions..."
        qfib_tmp = np.empty((4, fiber_ndiv*sum(numSpots)))
        for i in range(len(pd_hkl_ids)):
            for ispot in range(numSpots[i]):
                if not np.isnan(coms[i][ispot][0]):
                    ome_c = omeEta.omeEdges[0] + (0.5 + coms[i][ispot][0])*del_ome
                    eta_c = omeEta.etaEdges[0] + (0.5 + coms[i][ispot][1])*del_eta

                    gVec_s = xrdutil.makeMeasuredScatteringVectors(tTh[pd_hkl_ids[i]], eta_c, ome_c)

                    qfib_tmp[:, ii:jj] = rot.discreteFiber(pd.hkls[:, pd_hkl_ids[i]].reshape(3, 1),
                                                           gVec_s, B=bMat, ndiv=fiber_ndiv,
                                                           invert=False, csym=csym)[0]
                    ii  = jj
                    jj += fiber_ndiv
                    pass
                pass
            qfib.append(mutil.uniqueVectors(qfib_tmp))
            pass
        qfib = np.hstack(qfib)
    print "Running paintGrid on %d orientations" % (qfib.shape[1])
    complPG = idx.paintGrid(qfib,
                            omeEta,
                            omegaRange=omeRange, etaRange=etaRange,
                            omeTol=d2r*omeTol, etaTol=d2r*etaTol,
                            omePeriod=omePeriod, threshold=threshold,
                            doMultiProc=doMultiProc,
                            nCPUs=nCPUs)
    return complPG, qfib

def run_cluster(complPG, qfib, qsym,
                cl_radius=cl_radius, min_compl=min_compl):
    """
    """
    start = time.clock()                      # time this

    # # use transforms module for distance
    # quatDistance = lambda x, y: xf.quat_distance(x, y, qsym)

    # use compiled module for distance
    # just to be safe, must order qsym as C-contiguous
    qsym  = np.array(qsym.T, order='C').T
    quatDistance = lambda x, y: xfcapi.quat_distance(np.array(x, order='C'), \
                                                     np.array(y, order='C'), \
                                                     qsym)

    qfib_r = qfib[:, np.r_[complPG] > min_compl]

    print "Feeding %d orientations above %.1f%% to clustering" % (qfib_r.shape[1], 100*min_compl)

    if haveScikit:
        print "Using scikit..."
        pdist = pairwise_distances(qfib_r.T, metric=quatDistance, n_jobs=-1)
        core_samples, labels = dbscan(pdist, eps=d2r*cl_radius, min_samples=1, metric='precomputed')
        cl = np.array(labels, dtype=int) + 1
    else:
        print "Using fclusterdata with a tolerance of %f degrees..." % (cl_radius)
        cl = cluster.hierarchy.fclusterdata(qfib_r.T, d2r*cl_radius, criterion='distance', metric=quatDistance)

    nblobs = len(np.unique(cl))

    qbar = np.zeros((4, nblobs))
    for i in range(nblobs):
        npts = sum(cl == i + 1)
        # qbar[:, i] = mutil.unitVector(
        #     np.sum(qfib_r[:, cl == i + 1].reshape(4, npts), axis=1).reshape(4, 1)).flatten()
        qbar[:, i] = rot.quatAverage(qfib_r[:, cl == i + 1].reshape(4, npts), qsym)
    elapsed = (time.clock() - start)

    print "clustering took %f seconds" % (elapsed)
    return qbar, cl

if __name__ == "__main__":
    cfg_filename = sys.argv[1]
    hkl_ids = np.array(sys.argv[2:], dtype=int)

    print "Using cfg file '%s'" % (cfg_filename)

    parser = SafeConfigParser()
    parser.read(cfg_filename)

    # output for eta-ome maps as pickles
    working_dir   = parser.get('base', 'working_dir')
    analysis_name = parser.get('base', 'analysis_name')

    eta_ome_filename = os.path.join(working_dir, analysis_name + '-eta_ome.cpl')
    quats_filename   = os.path.join(working_dir, analysis_name + '-quats.out')
    compl_filename   = os.path.join(working_dir, analysis_name + '-compl.out')
    pull_filename    = os.path.join(working_dir, analysis_name + '-spots_%05d.out')

    print "analysis name is '%s'" %analysis_name

    pd, reader, detector = initialize_experiment(cfg_filename)

    if len(hkl_ids) == 0:
        hkl_ids = range(pd.hkls.shape[1])
        print "hkl ids not specified; grabbing from materials file..."
        print hkl_ids

    # some ome-eta parameters
    threshold   = parser.getfloat('ome_eta', 'threshold')
    nframesLump = parser.getint('ome_eta', 'nframesLump')

    # load stored maps ("if possible")
    load_maps_str = parser.get('ome_eta', 'load_maps')
    save_maps_str = parser.get('ome_eta', 'save_maps')
    if load_maps_str.strip() == '1' or load_maps_str.strip().lower() == 'true':
        load_maps = True
    elif load_maps_str.strip() == '' or load_maps_str.strip() == '0' or load_maps_str.strip().lower() == 'false':
        load_maps = False
    else:
        eta_ome_filename = os.path.join(working_dir, load_maps_str.strip())
        load_maps = True
    if load_maps:
        print "attempting to load stored maps from '%s'" %(eta_ome_filename)
        try:
            eta_ome_file = open(eta_ome_filename ,'r')
        except:
            load_maps = False
            raise RuntimeWarning, "load of eta ome maps failed...  making them instead"

    if load_maps:
        ome_eta = cPickle.load(eta_ome_file)
        pd = ome_eta.planeData
        hkl_ids = range(pd.hkls.shape[1])
        print "loaded stored maps, forcing overwrite of planeData and hkl_ids; overriding save preferences"
        save_maps_str = '0'
    else:
        if save_maps_str.strip() == '1' or save_maps_str.strip().lower() == 'true':
            eta_ome_outputname = eta_ome_filename
            save_maps = True
        elif save_maps_str.strip() == '' or save_maps_str.strip() == '0' or save_maps_str.strip().lower() == 'false':
            save_maps = False
        else:
            eta_ome_outputname = os.path.join(working_dir, save_maps_str.strip())
            save_maps = True
        if save_maps:
            try:
                eta_ome_out = open(eta_ome_outputname ,'w')
            except:
                eta_ome_out = None
                raise RuntimeWarning, "can't open output file; skipping save"
        # if you got here, make the maps
        ome_eta = make_maps(pd, reader, detector, hkl_ids, threshold, nframesLump, output=eta_ome_out)

    # more parameters
    seed_hkl_ids = [0,]
    threshold_pg = parser.getfloat('paint_grid', 'threshold')
    fiber_ndiv   = parser.getint('paint_grid', 'fiber_ndiv')
    multiproc    = parser.getboolean('paint_grid', 'multiproc')
    ncpus        = parser.get('paint_grid', 'ncpus')
    use_qgrid    = parser.getboolean('paint_grid', 'use_qgrid')
    omepd_str    = parser.get('paint_grid', 'ome_period')
    ome_period   = tuple(d2r*np.array(omepd_str.split(','), dtype=float))
    if use_qgrid:
        qgrid_file = parser.get('paint_grid', 'qgrid_file')
        print "using the file '%s' for quaternion search" % (qgrid_file)
    else:
        qgrid_file = None
        seed_hkl_str = parser.get('paint_grid', 'hkl_seeds')
        seed_hkl_ids = np.array(seed_hkl_str.split(','), dtype=int)
        print "using the following for seed hkls:"
        # print np.r_[hkl_ids][seed_hkl_ids]
        print pd.hkls[:, np.r_[hkl_ids][seed_hkl_ids]]
    # tolerances go in IN DEGREES
    ome_tol      = parser.getfloat('paint_grid', 'ome_tol')
    eta_tol      = parser.getfloat('paint_grid', 'eta_tol')
    restrict_eta = parser.getfloat('paint_grid', 'restrict_eta')

    etaRange = None
    if restrict_eta > 0:
        eta_del = d2r*abs(restrict_eta)
        etaRange = [[-0.5*np.pi + eta_del, 0.5*np.pi - eta_del],
                    [ 0.5*np.pi + eta_del, 1.5*np.pi - eta_del]]
        print "eta ranges restricted to:"
        print r2d*np.array(etaRange)
    else:
        print "using full eta range"
    print "omega period is:"
    print r2d*np.r_[ome_period]
    if ncpus.strip() == '':
        ncpus = multiprocessing.cpu_count()
    else:
        ncpus = int(ncpus)
    compl, qfib = run_paintGrid(pd, ome_eta, seed_hkl_ids, threshold_pg, fiber_ndiv,
                                omeTol=ome_tol, etaTol=eta_tol, etaRange=etaRange,
                                omePeriod=ome_period,
                                qTol=1e-7,
                                doMultiProc=multiproc,
                                nCPUs=ncpus,
                                useGrid=qgrid_file)

    np.savetxt(os.path.join(working_dir, compl_filename), compl)

    cl_radius = parser.getfloat('clustering', 'cl_radius')
    min_compl = parser.getfloat('clustering', 'min_compl')
    num_above = sum(np.r_[compl] > min_compl)
    if num_above == 0:
        raise RuntimeError, "No orientations above specified threshold of %.1f%%" % (100.*min_compl)
    elif num_above == 1:
        qbar = qfib[:, np.r_[compl] > min_compl]
    else:
        qbar, cl = run_cluster(compl, qfib, pd.getQSym(), cl_radius=cl_radius, min_compl=min_compl)
    qbar = np.atleast_2d(qbar)

    print "Found %d orientation clusters with >=%.1f%% completeness and %2f misorientation" %(qbar.size/4, 100.*min_compl, cl_radius)

    # SAVE OUTPUT
    np.savetxt(quats_filename, qbar.T, fmt="%1.12e", delimiter="\t")

    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # from hexrd.xrd import rotations as rot
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # phis, ns = rot.angleAxisOfRotMat(rot.rotMatOfQuat(qfib[:, np.r_[compl] > min_compl]))
    # rod = np.tile(np.tan(0.5*phis), (3, 1)) * ns
    # ax.scatter(rod[0, :], rod[1, :], rod[2, :], c='r', marker='o')
    #
    # phis, ns = rot.angleAxisOfRotMat(rot.rotMatOfQuat(qbar))
    # rod = np.tile(np.tan(0.5*phis), (3, 1)) * ns
    # ax.scatter(rod[0, :], rod[1, :], rod[2, :], c='b', marker='*')
    #
    # plt.show()

    # PULL SPOTS
    do_pull_str = parser.get('pull_spots', 'do_pull')
    if do_pull_str.strip() == '1' or do_pull_str.strip().lower() == 'true':

        pthresh = parser.getfloat('pull_spots', 'threshold')
        tth_tol = parser.getfloat('pull_spots', 'tth_tol')

        det_origin_str = parser.get('pull_spots', 'det_origin')
        det_origin = np.array(det_origin_str.split(','), dtype=float)

        geomParams = np.vstack([detector.getParams(allParams=True)[:6],
                                np.zeros(6)]).T

        distortion = (dFuncs.GE_41RT, detector.getParams(allParams=True)[6:])

        tVec_d = tVec_d_from_old_parfile(geomParams, det_origin)
        detector_params = np.hstack([geomParams[3:6, 0], tVec_d.flatten(), 0., np.zeros(3)])

        use_tth_max = parser.get('pull_spots', 'use_tth_max')
        if use_tth_max.strip() == '1' or use_tth_max.strip().lower() == 'true':
            excl = np.zeros_like(pd.exclusions, dtype=bool)
            pd.exclusions = excl
            excl = pd.getTTh() > detector.getTThMax()
            pd.exclusions = excl
            pass
        phi, n = rot.angleAxisOfRotMat(rot.rotMatOfQuat(qbar))
        if have_progBar:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=len(qbar.T)).start()
            pass
        print "pulling spots for %d orientations..." %len(qbar.T)
        for iq, quat in enumerate(qbar.T):
            if have_progBar:
                pbar.update(iq)
            exp_map = phi[iq]*n[:, iq]
            grain_params = np.hstack([exp_map.flatten(), 0., 0., 0., 1., 1., 1., 0., 0., 0.])
            sd = xrdutil.pullSpots(pd, detector_params, grain_params, reader,
                                   filename=pull_filename %iq,
                                   eta_range=etaRange, ome_period=ome_period,
                                   eta_tol=eta_tol, ome_tol=ome_tol,
                                   threshold=pthresh, tth_tol=tth_tol,
                                   distortion=distortion)
            pass
        if have_progBar:
            pbar.finish()
            pass
        pass # condidional on pull spots
    pass
