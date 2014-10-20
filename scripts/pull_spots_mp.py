# !/usr/bin/env python
#
# Multiprocessing pull_spots script.
#

import os, sys, time, datetime
import multiprocessing
import itertools

from ConfigParser import SafeConfigParser

import numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg.matfuncs import logm

from hexrd.xrd import fitting
from hexrd.xrd import material
from hexrd.xrd import xrdutil

from hexrd     import matrixutil as mutil
from hexrd     import coreutil
from hexrd.xrd import distortion as dFuncs
from hexrd.xrd import rotations  as rot
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi

from hexrd.xrd.detector import ReadGE

try:
    from progressbar import ProgressBar, Bar, ETA, Percentage
except:
    # Dummy no-op progress bar to simplify code using ProgressBar
    class ProgressBar(object):
        def __init__(*args, **kwargs):
            pass

        def start(self):
            return self

        def finish(self):
            pass

        def update(self, x):
            pass

    class Bar(object):
        pass

    class ETA(object):
        pass

    class Percentage(object):
        pass



class Config(object):
    def __init__(self, cfg_filename):
        self.pd, reader, self.detector = coreutil.initialize_experiment(cfg_filename)

        parser = SafeConfigParser()
        parser.read(cfg_filename)

        self.working_dir = parser.get('base', 'working_dir')
        self.analysis_name = parser.get('base', 'analysis_name')
        self.multiproc = parser.get('paint_grid', 'multiproc')
        ncpus_str = parser.get('paint_grid', 'ncpus')
        self.ncpus = get_ncpus(ncpus_str)

        self.threshold = parser.getfloat('pull_spots', 'threshold')

        self.threshold = parser.getfloat('pull_spots', 'threshold')
        det_origin_str = parser.get('pull_spots', 'det_origin')
        self.det_origin = np.array(det_origin_str.split(','), dtype=float)

        self.tth_tol = parser.getfloat('pull_spots', 'tth_tol')
        self.eta_tol = parser.getfloat('pull_spots', 'eta_tol')
        self.ome_tol = parser.getfloat('pull_spots', 'ome_tol')
        self.tth_tol_r = parser.getfloat('pull_spots', 'tth_tol_r')
        self.eta_tol_r = parser.getfloat('pull_spots', 'eta_tol_r')
        self.ome_tol_r = parser.getfloat('pull_spots', 'ome_tol_r')

        maxTTh_str = parser.get('pull_spots', 'use_tth_max')
        self.maxTTh = float(maxTTh_str) # in DEGREES

        material_name = parser.get('material', 'material_name')
        self.matl = material.loadMaterialList(os.path.join(self.working_dir, material_name+'.ini'))[0]
        self.old_par = np.loadtxt(parser.get('detector', 'parfile_name'))
        self.new_par = np.loadtxt(parser.get('pull_spots', 'parfile_name'))
        self.detector_params = self.new_par[:10]

        self.distortion = (dFuncs.GE_41RT, self.old_par[-6:,0])

        restrict_eta = parser.getfloat('paint_grid', 'restrict_eta')
        eta_del = d2r*abs(restrict_eta)
        self.etaRange = [[-0.5*np.pi + eta_del, 0.5*np.pi - eta_del],
                         [ 0.5*np.pi + eta_del, 1.5*np.pi - eta_del]]

        self.reader = read_frames(reader, parser)

        self.ome_start = parser.getfloat('reader', 'ome_start')
        self.ome_delta = parser.getfloat('reader', 'ome_delta')
        self.ome_stop = self.ome_start + len(self.reader[0])*self.ome_delta
        omepd_str = parser.get('paint_grid', 'ome_period')
        self.ome_period = tuple(d2r*np.array(omepd_str.split(','), dtype=float))


d2r = np.pi/180.
r2d = 180./np.pi

bVec_ref = xf.bVec_ref # reference beam vector (propagation) [0, 0, -1]
eta_ref  = xf.eta_ref  # eta=0 reference vector [1, 0, 0]
vInv_ref = xf.vInv_ref # reference inverse stretch [1, 1, 1, 0, 0, 0]

# grain parameter refinement flags
gFlag = np.array([1, 1, 1,
                  1, 1, 1,
                  1, 1, 1, 1, 1, 1], dtype=bool)
# grain parameter scalings
gScl  = np.array([1., 1., 1.,
                  1., 1., 1.,
                  1., 1., 1., 0.01, 0.01, 0.01])

def read_frames(reader, parser):
    start = time.time()                      # time this

    threshold = parser.getfloat('pull_spots', 'threshold')
    ome_start = parser.getfloat('reader', 'ome_start')     # in DEGREES
    ome_delta = parser.getfloat('reader', 'ome_delta')     # in DEGREES

    frame_list = []
    nframes = reader.getNFrames()
    print "Reading %d frames:" % nframes
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nframes).start()
    for i in range(nframes):
        frame = reader.read()
        frame[frame <= threshold] = 0
        frame_list.append(coo_matrix(frame))
        pbar.update(i+1)
    pbar.finish()
    # frame_list = np.array(frame_list)
    reader = [frame_list, [ome_start*d2r, ome_delta*d2r]]

    elapsed = (time.time() - start)
    print "Reading %d frames took %.2f seconds" % (nframes, elapsed)
    return reader

def get_ncpus(ncpus_str):
    cpucount = multiprocessing.cpu_count()
    if ncpus_str.strip() == '':
        ncpus = cpucount
    elif int(ncpus_str) == -1:
        ncpus = cpucount - 1
    elif int(ncpus_str) == -2:
        ncpus = cpucount / 2
    else:
        ncpus = int(ncpus_str)
    return ncpus

def process_grain(jobdata):
    # Unpack the job data
    jobnum = jobdata['job']
    quat = jobdata['quat']
    config = jobdata['config']

    if config.multiproc:
        # multiprocessing...
        # Redirect output to a process-specific logfile
        procnum = multiprocessing.current_process()._identity[0]
        logfile = open(config.analysis_name + '-proc%02d-log.out' % procnum, 'a')
        saved_stdout = sys.stdout
        saved_stderr = sys.stderr
        sys.stderr = sys.stdout = logfile

        print '\nTimestamp %s' % (datetime.datetime.utcnow().isoformat())
        print 'Process %d, job %d, quat %s' % (procnum, jobnum, quat)

    # put the job num on it
    fileroot = config.analysis_name + '_job_%05d' %jobnum
    filename = config.analysis_name + '-spots_%05d.out' %jobnum

    """
    ####### INITIALIZATION
    """
    ome_start, ome_delta = config.ome_start, config.ome_delta
    ome_period, ome_stop = config.ome_period, config.ome_stop
    pd = config.matl.planeData
    pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
    pd.exclusions = pd.getTTh() >= d2r*config.maxTTh

    bMat = np.ascontiguousarray(pd.latVecOps['B']) # hexrd convention; necessary to re-cast (?)
    wlen = pd.wavelength                           # Angstroms

    """
    ####### PULL SPOTS
    """
    results = np.zeros((1, 21))

    start = time.time()                      # time this
    print "fitting %d" %jobnum
    #
    phi   = 2*np.arccos(quat[0])
    n     = xf.unitVector(quat[1:].reshape(3, 1))
    #
    grain_params = np.hstack([phi*n.flatten(), np.zeros(3), np.ones(3), np.zeros(3)])
    #
    fid = open(filename, 'w')
    sd = xrdutil.pullSpots(pd, config.detector_params, grain_params, config.reader,
                           distortion=config.distortion,
                           eta_range=config.etaRange,
                           ome_period=ome_period,
                           tth_tol=config.tth_tol, eta_tol=config.eta_tol, ome_tol=config.ome_tol,
                           panel_buff=[10, 10],
                           npdiv=2, threshold=config.threshold, doClipping=False,
                           filename=fid)
    fid.close()

    # strain fitting
    for i in range(2):
        gtable  = np.loadtxt(filename) # load pull_spots output table
        idx0    = gtable[:, 0] >= 0             # select valid reflections
        #
        pred_ome = gtable[:, 6]
        if np.sign(ome_delta) < 0:
            idx_ome  = np.logical_and(pred_ome < d2r*(ome_start + 2*ome_delta),
                                      pred_ome > d2r*(ome_stop  - 2*ome_delta))
        else:
            idx_ome  = np.logical_and(pred_ome > d2r*(ome_start + 2*ome_delta),
                                      pred_ome < d2r*(ome_stop  - 2*ome_delta))
        #
        idx     = np.logical_and(idx0, idx_ome)
        hkls    = gtable[idx, 1:4].T            # must be column vectors
        xyo_det = gtable[idx, -3:]              # these are the cartesian centroids + ome
        xyo_det[:, 2] = xf.mapAngle(xyo_det[:, 2], ome_period)
        print "completeness: %f%%" %(100. * sum(idx)/float(len(idx)))
        if sum(idx) > 12:
            g_initial = grain_params
            g_refined = fitting.fitGrain(xyo_det, hkls, bMat, wlen,
                                         config.detector_params,
                                         g_initial[:3], g_initial[3:6], g_initial[6:],
                                         beamVec=bVec_ref, etaVec=eta_ref,
                                         distortion=config.distortion,
                                         gFlag=gFlag, gScl=gScl,
                                         omePeriod=ome_period)
            if i == 0:
                fid = open(filename, 'w')
                sd = xrdutil.pullSpots(pd, config.detector_params, g_refined, config.reader,
                                       distortion=config.distortion,
                                       eta_range=config.etaRange, ome_period=ome_period,
                                       tth_tol=config.tth_tol_r,
                                       eta_tol=config.eta_tol_r,
                                       ome_tol=config.ome_tol_r,
                                       panel_buff=[10, 10],
                                       npdiv=2, threshold=config.threshold,
                                       use_closest=True, doClipping=False,
                                       filename=fid)
                fid.close()
            pass
        else:
            g_refined = grain_params
            break
        pass
    eMat = logm(np.linalg.inv(mutil.vecMVToSymm(g_refined[6:])))

    resd_f2 = fitting.objFuncFitGrain(g_refined[gFlag], g_refined, gFlag,
                                      config.detector_params,
                                      xyo_det, hkls, bMat, wlen,
                                      bVec_ref, eta_ref,
                                      config.distortion[0], config.distortion[1],
                                      ome_period,
                                      simOnly=False)

    # Save the intermediate grain data as an npy file
    graindata = np.empty(21)
    graindata[:3] = (jobnum, sum(idx)/float(len(idx)), sum(resd_f2**2))
    graindata[3:15] = g_refined
    graindata[15:] = (eMat[0, 0], eMat[1, 1], eMat[2, 2], eMat[1, 2], eMat[0, 2], eMat[0, 1])
    np.save(os.path.join('pstmp', fileroot + '-grains.npy'), graindata)

    elapsed = (time.time() - start)
    print "grain %d took %.2f seconds" %(jobnum, elapsed)

    if config.multiproc:
        # Restore output
        sys.stdout = saved_stdout
        sys.stderr = saved_stderr
        logfile.close()

    return True

def process_all_grains(pool, jobdata, chunksize):
    print "Processing %d grains:" % len(jobdata)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=len(jobdata)).start()
    if pool is not None:
        for i, result in enumerate(pool.imap_unordered(process_grain, jobdata, chunksize)):
            pbar.update(i + 1)
    else:
        for i, jd in enumerate(jobdata):
            process_grain(jd)
            pbar.update(i + 1)
    pbar.finish()

"""
####### INPUT GOES HERE
"""
# def pull_spots_block(cfg_filename, blockID, pd, reader, detector):
if __name__ == "__main__":
    total_start = time.time()                      # time this
    cfg_filename = sys.argv[1]

    print "Using cfg file '%s'" % (cfg_filename)

    config = Config(cfg_filename)

    if len(sys.argv) < 3:
        quats_filename = config.analysis_name+'-quats.out'
    else:
        quats_filename = sys.argv[2]
    quats = np.loadtxt(os.path.join(config.working_dir, quats_filename))

    nquats = len(quats)
    jobdata = [{'job':i, 'quat':quat, 'config': config}
               for i, quat in enumerate(quats)]


    # material class
    material_name = parser.get('material', 'material_name')
    matl = material.loadMaterialList(os.path.join(working_dir, material_name+'.ini'))[0]
    
    # Temporary directory for intermediate files
    if not os.path.exists('pstmp'):
        os.mkdir('pstmp')

    # config.multiproc = False

    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nquats).start()
    if config.multiproc:
        print
        print 'Creating pool with %d processes' % config.ncpus
        pool = multiprocessing.Pool(config.ncpus)

        start = time.time()                      # time this

        # Controls how many jobs each process gets at once
        chunksize = 4

        print "Processing %d grains:" % nquats
        for i, result in enumerate(pool.imap_unordered(process_grain, jobdata, chunksize)):
            pbar.update(i + 1)

        elapsed = (time.time() - start)
        print "Processing %d grains took %.2f seconds" % (nquats, elapsed)
    else:
        print 'Running sequential'
        start = time.time()
        print "Processing %d grains:" % nquats
        for i, result in enumerate(itertools.imap(process_grain, jobdata)):
            pbar.update(i+1)
        elapsed = time.time() - start
        print "Processing %d grains took %.2f seconds" % (nquats, elapsed)

    pbar.finish()

    print "\nCopying grains to a single file"
    grains_file = open(config.analysis_name + '-grains.out', 'w')
    print >> grains_file, \
        "# grain ID\tcompleteness\tsum(resd**2)/n_refl\t" + \
        "xi[0]\txi[1]\txi[2]\t" + \
        "tVec_c[0]\ttVec_c[1]\ttVec_c[2]\t" + \
        "vInv_s[0]\tvInv_s[1]\tvInv_s[2]\tvInv_s[4]*sqrt(2)\tvInv_s[5]*sqrt(2)\tvInv_s[6]*sqrt(2)\t" + \
        "ln(V[0,0])\tln(V[1,1])\tln(V[2,2])\tln(V[1,2])\tln(V[0,2])\tln(V[0,1])"
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nquats).start()
    for jobnum in range(nquats):
        fileroot = config.analysis_name + '_job_%05d' %jobnum
        graindata = np.load(os.path.join('pstmp', fileroot + '-grains.npy'))
        print >> grains_file, \
            ("%d\t%1.7e\t%1.7e\t"
             "%1.7e\t%1.7e\t%1.7e\t"
             "%1.7e\t%1.7e\t%1.7e\t"
             "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t"
             "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e") % tuple(graindata)
        pbar.update(jobnum + 1)

    pbar.finish()
    grains_file.close()

    total_elapsed = (time.time() - total_start)
    print "\nTotal processing time %.2f seconds" % (total_elapsed)
