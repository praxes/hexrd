# !/usr/bin/env python
import os, sys, time

from ConfigParser import SafeConfigParser

import numpy as np
from scipy.sparse import dok_matrix
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

"""
####### INPUT GOES HERE
"""
# def pull_spots_block(cfg_filename, blockID, pd, reader, detector):
if __name__ == "__main__":
    cfg_filename = sys.argv[1]
    blockID      = int(sys.argv[2])
    gp_fileroot  = sys.argv[3]

    print "Using cfg file '%s'" % (cfg_filename)

    pd, reader, detector = coreutil.initialize_experiment(cfg_filename)

    parser = SafeConfigParser()
    parser.read(cfg_filename)

    # output for eta-ome maps as pickles
    working_dir   = parser.get('base', 'working_dir')
    analysis_name = parser.get('base', 'analysis_name')

    restrict_eta = parser.getfloat('paint_grid', 'restrict_eta')
    omepd_str    = parser.get('paint_grid', 'ome_period')
    ome_period   = tuple(d2r*np.array(omepd_str.split(','), dtype=float))

    threshold      = parser.getfloat('pull_spots', 'threshold')
    det_origin_str = parser.get('pull_spots', 'det_origin')
    det_origin     = np.array(det_origin_str.split(','), dtype=float)

    # for spot pulling; got this from GUI
    tth_tol    = parser.getfloat('pull_spots', 'tth_tol')
    eta_tol    = parser.getfloat('pull_spots', 'eta_tol')
    ome_tol    = parser.getfloat('pull_spots', 'ome_tol')
    tth_tol_r  = parser.getfloat('pull_spots', 'tth_tol_r')
    eta_tol_r  = parser.getfloat('pull_spots', 'eta_tol_r')
    ome_tol_r  = parser.getfloat('pull_spots', 'ome_tol_r')

    maxTTh_str = parser.get('pull_spots', 'use_tth_max')
    maxTTh = float(maxTTh_str)            # in DEGREES
    # if maxTTh_str.strip() == '1' or maxTTh_str.strip().lower() == 'true':
    #     maxTTh = detector.getTThMax()*r2d
    # elif maxTTh_str.strip() != '0' or maxTTh_str.strip().lower() == 'false':
    #     maxTTh = float(maxTTh_str)            # in DEGREES

    # put a block ID on it
    fileroot = analysis_name + '_block_%03d' %blockID
    filename = fileroot + '-spots_%05d.out'

    """
    ####### INITIALIZATION
    """
    # material class
    material_name = parser.get('material', 'material_name')
    matl = material.loadMaterialList(os.path.join(working_dir, material_name+'.ini'))[0]

    # planeData and reader
    pd = matl.planeData
    pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
    pd.exclusions = pd.getTTh() >= d2r*maxTTh

    bMat = np.ascontiguousarray(pd.latVecOps['B']) # hexrd convention; necessary to re-cast (?)
    wlen = pd.wavelength                           # Angstroms

    # parameters for detector
    old_par = np.loadtxt(parser.get('detector',   'parfile_name'))
    new_par = np.loadtxt(parser.get('pull_spots', 'parfile_name'))

    detector_params = new_par[:10]

    dFunc   = dFuncs.GE_41RT
    dParams = old_par[-6:, 0]                 # MUST CHANGE THIS

    # need this below for cases where full 360 isn't here
    ome_start = parser.getfloat('reader', 'ome_start')     # in DEGREES
    ome_delta = parser.getfloat('reader', 'ome_delta')     # in DEGREES

    frame_list = []
    for i in range(reader.getNFrames()):
        frame = reader.read()
        frame[frame <= threshold] = 0
        frame_list.append(dok_matrix(frame))
    frame_list = np.array(frame_list)
    reader = [frame_list, [ome_start*d2r, ome_delta*d2r]]

    ome_stop = ome_start + len(reader[0])*ome_delta

    # restrict eta range
    #  - important for killng edge cases near eta=+/-90
    eta_del = d2r*abs(restrict_eta)
    etaRange = [[-0.5*np.pi + eta_del, 0.5*np.pi - eta_del],
                [ 0.5*np.pi + eta_del, 1.5*np.pi - eta_del]]

    """
    ####### PULL SPOTS
    """
    gp_table = np.loadtxt(gp_fileroot + '_block_%03d-grains.out' %blockID)
    n_grains = len(gp_table)
    full_results = np.zeros((n_grains, 21))
    grains_file = open(fileroot + '-grains.out', 'w')
    print >> grains_file, \
      "# grain ID\tcompleteness\tsum(resd**2)/n_refl\t" + \
      "xi[0]\txi[1]\txi[2]\t" + \
      "tVec_c[0]\ttVec_c[1]\ttVec_c[2]\t" + \
      "vInv_s[0]\tvInv_s[1]\tvInv_s[2]\tvInv_s[4]*sqrt(2)\tvInv_s[5]*sqrt(2)\tvInv_s[6]*sqrt(2)\t" + \
      "ln(V[0,0])\tln(V[1,1])\tln(V[2,2])\tln(V[1,2])\tln(V[0,2])\tln(V[0,1])"

    # loop over grains from previous fit
    for grainID in range(n_grains):
        start = time.clock()                      # time this
        print "fitting %d" %grainID
        #
        grain_params = gp_table[grainID, 3:15]
        #
        fid = open(filename %grainID, 'w')
        sd = xrdutil.pullSpots(pd, detector_params, grain_params, reader,
                               distortion=(dFunc, dParams),
                               eta_range=etaRange, ome_period=ome_period,
                               tth_tol=tth_tol, eta_tol=eta_tol, ome_tol=ome_tol,
                               panel_buff=[10, 10],
                               npdiv=2, threshold=threshold, doClipping=False,
                               filename=fid)
        fid.close()

        # strain fitting
        for i in range(2):
            gtable  = np.loadtxt(filename %grainID) # load pull_spots output table
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
                                             detector_params,
                                             g_initial[:3], g_initial[3:6], g_initial[6:],
                                             beamVec=bVec_ref, etaVec=eta_ref,
                                             distortion=(dFunc, dParams),
                                             gFlag=gFlag, gScl=gScl,
                                             omePeriod=ome_period)
                if i == 0:
                    fid = open(filename %grainID, 'w')
                    sd = xrdutil.pullSpots(pd, detector_params, g_refined, reader,
                                           distortion=(dFunc, dParams),
                                           eta_range=etaRange, ome_period=ome_period,
                                           tth_tol=tth_tol_r, eta_tol=eta_tol_r, ome_tol=ome_tol_r,
                                           panel_buff=[10, 10],
                                           npdiv=2, threshold=threshold,
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
                                          detector_params,
                                          xyo_det, hkls, bMat, wlen,
                                          bVec_ref, eta_ref,
                                          dFunc, dParams,
                                          ome_period,
                                          simOnly=False)
        print >> grains_file, \
          "%d\t%1.7e\t%1.7e\t"                         % (grainID, sum(idx)/float(len(idx)), sum(resd_f2**2)) + \
          "%1.7e\t%1.7e\t%1.7e\t"                      % tuple(g_refined[:3]) + \
          "%1.7e\t%1.7e\t%1.7e\t"                      % tuple(g_refined[3:6]) + \
          "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t" % tuple(g_refined[6:]) + \
          "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e"   % (eMat[0, 0], eMat[1, 1], eMat[2, 2], eMat[1, 2], eMat[0, 2], eMat[0, 1])
        elapsed = (time.clock() - start)
        print "grain %d took %.2f seconds" %(grainID, elapsed)
        pass
    grains_file.close()
    # return
