import sys, os, time
from numbapro import vectorize, float64, jit, guvectorize, autojit
import numpy as np
from hexrd.xrd import nbdistortion as dFuncs
from hexrd.xrd import transforms as xf 
from hexrd.xrd.transforms import makeDetectorRotMat, makeRotMatOfExpMap, makeOscillRotMat, makeEtaFrameRotMat, rotate_vecs_about_axis, makeBinaryRotMat

from hexrd.xrd import pycfuncs_transforms as xfcapi

#epsf = np.finfo(float).eps      # ~2.2e-16
#
#sqrt_epsf = np.sqrt(epsf)            # ~1.5e-8
## basis vector
#I3 = np.eye(3) 
#Xl = np.array([[1., 0., 0.]], order='C').T     # X in the lab frame
#
#Zl = np.array([[0., 0., 1.]]).T     # Z in the lab frame
#
#bVec_ref = -Zl 
#eta_ref = Xl
#
#
## distortion for warping detector coords
#dFunc_ref   = dFuncs.GE_41RT
#dParams_ref = [0., 0., 0., 2., 2., 2]
#

 


if __name__ == '__main__':
    bVec_ref = xf.bVec_ref

    rMat_d = makeDetectorRotMat( ( 0.0011546340766314521,
                                     -0.0040527538387122993,
                                     -0.0026221336905160211 ) ) 

#
    tVec_d = np.array( [ [   -1.44904 ],
                         [   -3.235616],
                         [-1050.74026 ] ] )

    chi    = -0.0011591608938627839
    tVec_s = np.array([ [-0.15354144],
                        [ 0.        ],
                        [-0.23294777] ] )

    rMat_c = makeRotMatOfExpMap(np.array( [ [ 0.66931818],
                                               [-0.98578066],
                                               [ 0.73593251] ] ) )
    tVec_c = np.array( [ [ 0.07547626],
                         [ 0.08827523],
                         [-0.02131205] ] )

    rMat_s = makeOscillRotMat([chi, 0.])

    # ######################################################################
    # Calculate pixel coordinates
    #
    #pvec  = 204.8 * np.linspace(-1, 1, 2048)
    pvec  = 204.8 * np.linspace(-1, 1, 512)
    dcrds = np.meshgrid(pvec, pvec)
    XY    = np.vstack([dcrds[0].flatten(), dcrds[1].flatten()]).T

    # Check the timings
    start1 = time.clock()                      # time this
    dangs1 = xf.detectorXYToGvec(XY, rMat_d, rMat_s,
                                 tVec_d, tVec_s, tVec_c, 
                                 beamVec=bVec_ref)
    tTh_d1   = dangs1[0][0]
    eta_d1   = dangs1[0][1]
    gVec_l1  = dangs1[1]
    elapsed1 = (time.clock() - start1)
    print "Time for Python detectorXYToGvec: %f"%(elapsed1)
    start2 = time.clock()                      # time this
    #dangs2 = xfcapi.detectorXYToGvec(XY, rMat_d, rMat_s,
   #                                  tVec_d, tVec_s, tVec_c,
    #                                 beamVec=bVec_ref)
    #tTh_d2   = dangs2[0][0]
    #eta_d2   = dangs2[0][1]
    #gVec_l2  = dangs1[1]
   # elapsed2 = (time.clock() - start2)
    #print "Time for CAPI detectorXYToGvec: %f"%(elapsed2)

   # maxDiff_tTh = np.linalg.norm(tTh_d1-tTh_d2,np.inf)
   # print "Maximum disagreement in tTh:  %f"%maxDiff_tTh
   # maxDiff_eta = np.linalg.norm(eta_d1-eta_d2,np.inf)
   ## print "Maximum disagreement in eta:  %f"%maxDiff_eta
   # maxDiff_gVec = np.linalg.norm(np.sqrt(np.sum(np.asarray(gVec_l1.T-gVec_l2)**2,1)),np.inf)
   # print "Maximum disagreement in gVec: %f"%maxDiff_gVec



    gVec_c1 = np.dot(rMat_c.T,np.dot(rMat_s.T,gVec_l1))
    #gVec_c1 = nb_dot(rMat_c.T,np.dot(rMat_s.T,gVec_l1))
#    gVec_c2 = np.ascontiguousarray(np.dot(rMat_c.T,np.dot(rMat_s.T,gVec_l2.T)).T)
    #gVec_c2 = np.ascontiguousarray(nb_dot(rMat_c.T,np.dot(rMat_s.T,gVec_l2.T)).T)

    start3 = time.clock()                      # time this

    #xy1 = xfcapi.gvecToDetectorXY(gVec_c1,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,bVec_ref)
    bHat_l = np.zeros(3)
    nVec_l = np.zeros(3)
    P0_l = np.zeros(3)
    P2_l = np.zeros(3)
    P2_d = np.zeros(3)
    P3_l = np.zeros(3)
    gHat_c = np.zeros(3)
    gVec_l = np.zeros(3)
    dVec_l = np.zeros(3)
    rMat_sc = np.zeros(9)
    brMat = np.zeros(9) 
    result = np.empty((gVec_c1.T.shape[0], 3))
    bVec_ref= bVec_ref.flatten()

#
    xy1 = xfcapi.gvecToDetectorXY(gVec_c1.T,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,bVec_ref, bHat_l, nVec_l, P0_l, P2_l, P2_d, P3_l, gHat_c, gVec_l, dVec_l, rMat_sc, brMat, result)
    #xy1 = xfcapi.gvecToDetectorXY(gVec_c2,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,bVec_ref)
    elapsed3 = (time.clock() - start3)
    print "Time for Python gvecToDetectorXY: %f"%(elapsed3)
#    import pdb
#    pdb.set_trace()

#    maxDiff_xy = np.linalg.norm(np.sqrt(np.sum(np.asarray(XY-xy1)**2,1)),np.inf)
#    print "Maximum disagreement in gVec: %f"%maxDiff_xy

    start4 = time.clock()                      # time this
    #xy2 = xfcapi.gvecToDetectorXY(gVec_c2,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,beamVec=bVec_ref)
#    xy2 = xfcapi.gvecToDetectorXY(gVec_c2,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,beamVec=bVec_ref)
    elapsed4 = (time.clock() - start4)
    print "Time for CAPI gvecToDetectorXY: %f"%(elapsed4)

 #   maxDiff_xy = np.linalg.norm(np.sqrt(np.sum(np.asarray(XY-xy2)**2,1)),np.inf)
    #print "Maximum disagreement in gVec: %f"%maxDiff_xy




