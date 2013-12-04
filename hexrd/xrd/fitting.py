import numpy as np
from scipy import optimize as opt

from hexrd.xrd import transforms as xf
from hexrd.xrd import distortion as dFuncs

import pdb

# ######################################################################
# Module Data
d2r = np.pi/180.
r2d = 180./np.pi

dFunc_ref   = dFuncs.GE_41RT
dParams_ref = [0., 0., 0., 2., 2., 2]
dFlag_ref   = np.array([0, 0, 0, 0, 0, 0], dtype=bool)
dScl_ref    = np.array([1, 1, 1, 1, 1, 1], dtype=float)
pFlag_ref   = np.array([1, 1, 1, 
                        1, 1, 1, 
                        0, 
                        0, 0, 0, 
                        0, 0, 0, 
                        0, 0, 0], dtype=bool)
pScl_ref    = np.array([1, 1, 1, 
                        1, 1, 1, 
                        1, 
                        1, 1, 1, 
                        1, 1, 1, 
                        1, 1, 1])
bVec_ref    = xf.bVec_ref
eta_ref     = xf.eta_ref

# ######################################################################
# Funtions

def calibrateDetectorFromSX(xyo_det, hkls_idx, bMat, wavelength, 
                            tiltAngles, chi, expMap_c,
                            tVec_d, tVec_s, tVec_c, 
                            beamVec=bVec_ref, etaVec=eta_ref, 
                            distortion=(dFunc_ref, dParams_ref, dFlag_ref, dScl_ref), 
                            pFlag=pFlag_ref, pScl=pScl_ref):
    """
    """
    dFunc   = distortion[0]
    dParams = distortion[1]
    dFlag   = distortion[2]
    dScl    = distortion[3]
    
    # p = np.zeros(16)
    # 
    # p[0]  = tiltAngles[0]
    # p[1]  = tiltAngles[1]
    # p[2]  = tiltAngles[2]
    # p[3]  = tVec_d[0]
    # p[4]  = tVec_d[1]
    # p[5]  = tVec_d[2]
    # p[6]  = chi
    # p[7]  = tVec_s[0]
    # p[8]  = tVec_s[1]
    # p[9]  = tVec_s[2]
    # p[10] = expMap_c[0]
    # p[11] = expMap_c[1]
    # p[12] = expMap_c[2]
    # p[13] = tVec_c[0]
    # p[14] = tVec_c[1]
    # p[15] = tVec_c[2]
    # 
    # pFull = np.hstack([p, dParams])

    pFull = geomParamsToInput(tiltAngles, chi, expMap_c,
                              tVec_d, tVec_s, tVec_c, 
                              dParams)

    refineFlag = np.hstack([pFlag, dFlag])
    scl        = np.hstack([pScl, dScl])
    pFit       = pFull[refineFlag] 
    fitArgs    = (pFull, pFlag, dFunc, dFlag, xyo_det, hkls_idx, bMat, wavelength, beamVec, etaVec)
    
    results = opt.leastsq(objFuncSX, pFit, args=fitArgs, diag=scl[refineFlag].flatten(), factor=1, ftol=1e-10)
    
    pFit_opt = results[0]

    retval = pFull
    retval[refineFlag] = pFit_opt
    return retval

# pdb.set_trace()
def objFuncSX(pFit, pFull, pFlag, dFunc, dFlag, 
              xyo_det, hkls_idx, bMat, wavelength, 
              bVec, eVec):

    npts   = len(xyo_det)
    
    refineFlag = np.hstack([pFlag, dFlag])
    
    # pFull[refineFlag] = pFit/scl[refineFlag]
    pFull[refineFlag] = pFit
    
    dParams = pFull[-len(dFlag):]
    xy_unwarped = dFunc(xyo_det[:, :2], dParams)
    
    rMat_d = xf.makeDetectorRotMat(pFull[:3])
    tVec_d = pFull[3:6].reshape(3, 1)

    chi = pFull[6]

    tVec_s = pFull[7:10].reshape(3, 1)
    
    gVec_c = np.dot(bMat, hkls_idx)
    rMat_c = xf.makeRotMatOfExpMap(pFull[10:13])
    tVec_c = pFull[13:16].reshape(3, 1)
    
    # get omegas for rMat_s calculation
    meas_omes  = xyo_det[:, 2] 
    
    oangs0, oangs1 = xf.oscillAnglesOfHKLs(hkls_idx, chi, rMat_c, bMat, wavelength, 
                                           beamVec=bVec, etaVec=eVec)

    calc_omes  = np.vstack( [xf.mapAngle(oangs0[2, :]), 
                             xf.mapAngle(oangs1[2, :]) ] )
    match_omes = np.argsort(abs(np.tile(meas_omes, (2, 1)) - calc_omes), axis=0) == 0
    calc_omes  = calc_omes.T.flatten()[match_omes.T.flatten()]
    
    xy_det = np.zeros((npts, 2))
    for i in range(npts):
        rMat_s = xf.makeOscillRotMat([chi, calc_omes[i]])
        xy_det[i, :] = xf.gvecToDetectorXY(gVec_c[:, i].reshape(3, 1), 
                                           rMat_d, rMat_s, rMat_c, 
                                           tVec_d, tVec_s, tVec_c, 
                                           beamVec=bVec).flatten()
        pass
    if np.any(np.isnan(xy_det)):
        print "infeasible pFull"
    
    # return values
    # retval = np.sum((xy_det - xy_unwarped[:, :2])**2)
    retval = np.sum( (xy_det - xy_unwarped[:, :2])**2, axis=1)
    return retval

def geomParamsToInput(tiltAngles, chi, expMap_c,
                      tVec_d, tVec_s, tVec_c, 
                      dParams):
    """
    """
    p = np.zeros(16)

    p[0]  = tiltAngles[0]
    p[1]  = tiltAngles[1]
    p[2]  = tiltAngles[2]
    p[3]  = tVec_d[0]
    p[4]  = tVec_d[1]
    p[5]  = tVec_d[2]
    p[6]  = chi
    p[7]  = tVec_s[0]
    p[8]  = tVec_s[1]
    p[9]  = tVec_s[2]
    p[10] = expMap_c[0]
    p[11] = expMap_c[1]
    p[12] = expMap_c[2]
    p[13] = tVec_c[0]
    p[14] = tVec_c[1]
    p[15] = tVec_c[2]

    return np.hstack([p, dParams])

def inputToGeomParams(p):
    """
    """
    retval = {}

    retval['tiltAngles'] = (p[0], p[1], p[2])
    retval['tVec_d']     = np.c_[p[3], p[4], p[5]].T
    retval['chi']        = p[6]
    retval['tVec_s']     = np.c_[p[7], p[8], p[9]].T
    retval['expMap_c']   = np.c_[p[10], p[11], p[12]].T
    retval['tVec_c']     = np.c_[p[13], p[14], p[15]].T
    retval['dParams']    = p[16:]
              
    return retval

