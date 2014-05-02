import sys, os, time

import shelve, cPickle
import numpy as np

import scipy.cluster as cluster
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
from hexrd     import matrixutil as mutil
from hexrd.xrd import experiment as expt
from hexrd.xrd import indexer    as idx
from hexrd.xrd import rotations  as rot
from hexrd.xrd import transforms as xf

from hexrd.xrd          import xrdutil
from hexrd.xrd.xrdbase  import multiprocessing
from hexrd.xrd.detector import ReadGE

haveScikit = False
# try:
#     from sklearn.cluster          import dbscan
#     from sklearn.metrics.pairwise import pairwise_distances
#     haveScikit = True
# except:
#     print "System does not have SciKit installed, using scipy fallback"

"""
##################### BEGIN COMPUTATION HERE #####################
"""

# constants
d2r = np.pi/180.
r2d = 180./np.pi

min_compl = 0.5
cl_radius = 1.0

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
    dark      = os.path.join(image_dir, parser.get('reader', 'dark')) 
    doFlip    = parser.getboolean('reader', 'doFlip')
    flipArg   = parser.get('reader', 'flipArg')
    
    # make frame reader
    reader   = ReadGE(fileInfo, ome_start, ome_delta,
                      subtractDark=True, dark=dark,
                      doFlip=doFlip, flipArg=flipArg)
    
    # DETECTOR
    ws.loadDetector(os.path.join(working_dir, detector_fname))
    
    return pd, reader, ws.detector

def make_maps(pd, reader, detector, hkl_ids, threshold, nframesLump, output=None):
    """
    OME-ETA MAPS
    """
    nEtaBins = np.int(2*np.pi / reader.getDeltaOmega()) / nframesLump
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
                  omeTol=None, etaTol=None, qTol=1e-7, 
                  doMultiProc=True, nCPUs=multiprocessing.cpu_count(), 
                  useGrid=None):
    """
    """
    del_ome = omeEta.omegas[1] - omeEta.omegas[0]
    del_eta = omeEta.etas[1] - omeEta.etas[0]
    
    omegaRange = [ ( np.amin(omeEta.omeEdges), np.amax(omeEta.omeEdges) ), ]
    
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
            qfib = np.loadtxt(useGrid).T
            print "loading quaternion grid file: %s" % (useGrid)
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
            coms_t               = ndimage.measurements.center_of_mass(omeEta.dataStore[i], labels_t, range(np.amax(labels_t) + 1))
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
                    ome_c = omeEta.omeEdges[0] + coms[i][ispot][0]*del_ome
                    eta_c = omeEta.etaEdges[0] + coms[i][ispot][1]*del_eta
                    
                    gVec_s = xrdutil.makeMeasuredScatteringVectors(tTh[pd_hkl_ids[i]], eta_c, ome_c)
                    
                    qfib_tmp[:, ii:jj] = rot.discreteFiber(pd.hkls[:, pd_hkl_ids[i]].reshape(3, 1), 
                                                           gVec_s, B=bMat, ndiv=fiber_ndiv, 
                                                           invert=False, csym=csym)[0]
                    ii  = jj
                    jj += fiber_ndiv
                    pass
                pass
            qfib.append(qfib_tmp)
            pass
        print "condensing quaternions..."
        try:
            qfib = mutil.uniqueVectors(np.hstack(qfib), tol=qTol)
        except:
            raise RuntimeError, "No seed points found (numSpots = 0)"
    
    print "Running paintGrid on %d orientations" % (qfib.shape[1])
    complPG = idx.paintGrid(qfib,
                            omeEta,
                            omegaRange=omegaRange,
                            omeTol=d2r*omeTol, etaTol=d2r*etaTol,
                            threshold=threshold,
                            doMultiProc=doMultiProc,
                            nCPUs=nCPUs)
    return complPG, qfib

def run_cluster(complPG, qfib, qsym, 
                cl_radius=cl_radius, min_compl=min_compl):
    """
    """
    start = time.clock()                      # time this
    
    quatDistance = lambda x, y: xf.quat_distance(x, y, qsym)
    
    qfib_r = qfib[:, np.r_[complPG] > min_compl]
    
    print "Feeding %d orientations above %.1f%% to clustering" % (qfib_r.shape[1], 100*min_compl)
    
    if haveScikit:
        print "Using scikit..."
        pdist = pairwise_distances(qfib_r.T, Y=None, metric=quatDistance, n_jobs=-2)
        core_samples, labels = dbscan(pdist, eps=d2r*cl_radius, min_samples=2, metric='precomputed')
        cl = np.array(labels, dtype=int) + 1
    else:
        print "Using fclusterdata with a tolerance of %f degrees..." % (cl_radius)
        cl = cluster.hierarchy.fclusterdata(qfib_r.T, d2r*cl_radius, criterion='distance', metric=quatDistance)
        
    nblobs = len(np.unique(cl))
    
    qbar = np.zeros((4, nblobs))
    for i in range(nblobs):
        npts = sum(cl == i + 1)
        qbar[:, i] = mutil.unitVector(
            np.sum(qfib[:, np.r_[complPG] > min_compl][:, cl == i + 1].reshape(4, npts), axis=1).reshape(4, 1)).flatten()
    elapsed = (time.clock() - start)
    
    print "clustering took %f seconds" % (elapsed)
    return qbar, cl

if __name__ == "__main__":
    cfg_filename = sys.argv[1]
    out_filename = sys.argv[2]
    hkl_ids = np.array(sys.argv[3:], dtype=int)
    
    print "Using cfg file '%s'" % (cfg_filename)

    parser = SafeConfigParser()
    parser.read(cfg_filename)

    # output for eta-ome maps as pickles
    working_dir   = parser.get('base', 'working_dir')
    analysis_name = parser.get('base', 'analysis_name')
     
    eta_ome_file = open(
        os.path.join(working_dir, analysis_name + '-eta_ome.cpl'),
        'w')
    
    pd, reader, detector = initialize_experiment(cfg_filename)

    if len(hkl_ids) == 0:
        hkl_ids = range(pd.hkls.shape[1])
        print "hkl ids not specified; grabbing from materials file..."
        print hkl_ids
    
    threshold   = parser.getfloat('ome_eta', 'threshold')
    nframesLump = parser.getint('ome_eta', 'nframesLump')
    ome_eta = make_maps(pd, reader, detector, hkl_ids, threshold, nframesLump, output=eta_ome_file)
    
    seed_hkl_ids = [0,]
    threshold_pg = parser.getfloat('paint_grid', 'threshold')
    fiber_ndiv   = parser.getint('paint_grid', 'fiber_ndiv')
    multiproc    = parser.getboolean('paint_grid', 'multiproc')
    ncpus        = parser.get('paint_grid', 'ncpus')
    use_qgrid    = parser.getboolean('paint_grid', 'use_qgrid')
    if use_qgrid:
        qgrid_file = parser.get('paint_grid', 'qgrid_file')
    else:
        qgrid_file = None
        seed_hkl_str = parser.get('paint_grid', 'hkl_seeds')
        seed_hkl_ids = np.array(seed_hkl_str.split(','), dtype=int)
        print "using the following for seed hkls:"
        print np.r_[hkl_ids][hkl_ids]
    
    # tolerances go in IN DEGREES 
    ome_tol      = parser.getfloat('paint_grid', 'ome_tol') 
    eta_tol      = parser.getfloat('paint_grid', 'eta_tol') 

    if ncpus.strip() == '':
        ncpus = multiprocessing.cpu_count()
    else:
        ncpus = int(ncpus)
    compl, qfib = run_paintGrid(pd, ome_eta, seed_hkl_ids, threshold_pg, fiber_ndiv, 
                                omeTol=ome_tol, etaTol=eta_tol, qTol=1e-7, 
                                doMultiProc=multiproc, 
                                nCPUs=ncpus, 
                                useGrid=qgrid_file)

    cl_radius = parser.getfloat('clustering', 'cl_radius')
    min_compl = parser.getfloat('clustering', 'min_compl')
    num_above = sum(np.r_[compl] > min_compl)
    if num_above == 0:
        import pdb;pdb.set_trace()
        # raise RuntimeError, "No orientations above specified threshold of %.1f%%" % (100.*min_compl)
    elif num_above == 1:
        qbar = qfib[:, np.r_[compl] > min_compl]
    else:
        qbar, cl = run_cluster(compl, qfib, pd.getQSym(), cl_radius=cl_radius, min_compl=min_compl)
    
    # SAVE OUTPUT
    np.savetxt(os.path.join(working_dir, out_filename), qbar.T, fmt="%1.12e", delimiter="\t")
    
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
    pass
