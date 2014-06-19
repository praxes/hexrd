from timeit import default_timer as timer
import sys, os
#from numbapro import vectorize, float64, jit, guvectorize, autojit
import numpy as np
#from hexrd.xrd import nbdistortion as dFuncs
from hexrd.xrd import transforms as xf 
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import pycfuncs_transforms as pycfuncs 
import numba.cuda

#input parameters
bVec_ref = xf.bVec_ref

rMat_d = xf.makeDetectorRotMat( ( 0.0011546340766314521,
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

rMat_c = xf.makeRotMatOfExpMap(np.array( [ [ 0.66931818],
                                           [-0.98578066],
                                           [ 0.73593251] ] ) )
tVec_c = np.array( [ [ 0.07547626],
                     [ 0.08827523],
                     [-0.02131205] ] )

rMat_s = xf.makeOscillRotMat([chi, 0.])

# ######################################################################
# Calculate pixel coordinates
#
pvec  = 204.8 * np.linspace(-1, 1, 2048)
#pvec  = 204.8 * np.linspace(-1, 1, 512)
dcrds = np.meshgrid(pvec, pvec)
XY    = np.vstack([dcrds[0].flatten(), dcrds[1].flatten()]).T

# Check the timings
start1 = timer()                      # time this
dangs1 = xf.detectorXYToGvec(XY, rMat_d, rMat_s,
                             tVec_d, tVec_s, tVec_c, 
                             beamVec=bVec_ref)
tTh_d1   = dangs1[0][0]
eta_d1   = dangs1[0][1]
gVec_l1  = dangs1[1]
elapsed1 = (timer() - start1)
print "Time for Python detectorXYToGvec: %f"%(elapsed1)

start2 = timer()                      # time this
dangs2 = xfcapi.detectorXYToGvec(XY, rMat_d, rMat_s,
                                 tVec_d, tVec_s, tVec_c,
                                 beamVec=bVec_ref)
tTh_d2   = dangs2[0][0]
eta_d2   = dangs2[0][1]
gVec_l2  = dangs2[1]
elapsed2 = (timer() - start2)
print "Time for CAPI detectorXYToGvec: %f"%(elapsed2)

# maxDiff_tTh = np.linalg.norm(tTh_d1-tTh_d2,np.inf)
# print "Maximum disagreement in tTh:  %f"%maxDiff_tTh
# maxDiff_eta = np.linalg.norm(eta_d1-eta_d2,np.inf)
## print "Maximum disagreement in eta:  %f"%maxDiff_eta
# maxDiff_gVec = np.linalg.norm(np.sqrt(np.sum(np.asarray(gVec_l1.T-gVec_l2)**2,1)),np.inf)
# print "Maximum disagreement in gVec: %f"%maxDiff_gVec



gVec_c1 = np.dot(rMat_c.T,np.dot(rMat_s.T,gVec_l1))
gVec_c2 = np.ascontiguousarray(np.dot(rMat_c.T,np.dot(rMat_s.T,gVec_l2.T)).T)

start3 = timer()                      # time this
xy1 = xf.gvecToDetectorXY(gVec_c1,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,beamVec=bVec_ref)
elapsed3 = (timer() - start3)
print "Time for Python gvecToDetectorXY: %f"%(elapsed3)

start4 = timer()                      # time this
xy2 = xfcapi.gvecToDetectorXY(gVec_c2,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,beamVec=bVec_ref)
elapsed4 = (timer() - start4)
print "Time for CAPI gvecToDetectorXY: %f"%(elapsed4)


print 'cudadevice: ', numba.cuda.get_current_device().name


# setup or numba version
# should be able to run in nopython mode
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
result = np.empty((gVec_c2.shape[0], 3))
bVec_ref= bVec_ref.flatten()

start5 = timer()
pycfuncs.gvecToDetectorXY(gVec_c2,rMat_d,rMat_s,rMat_c,tVec_d,tVec_s,tVec_c,bVec_ref, bHat_l, nVec_l, P0_l, P2_l, P2_d, P3_l, gHat_c, gVec_l, dVec_l, rMat_sc, brMat, result)
elapsed5 = (timer() - start5)
print "Time for Numba PyCFuncs gvecToDetectorXY: %f"%(elapsed5)

xy3 = result[:, 0:2] # result should be a 2 col matrix but an issue with numba

# debug
#    import pdb
#    pdb.set_trace()

#    maxDiff_xy = np.linalg.norm(np.sqrt(np.sum(np.asarray(XY-xy1)**2,1)),np.inf)
#    print "Maximum disagreement in gVec: %f"%maxDiff_xy


print 'testing data ...'

print xy2
print '*' * 50
print xy3


assert np.allclose(xy1, xy2), 'xy1 not close to xy3'
assert np.allclose(xy2, xy3), 'C and PyC not close: xy2, xy3'
print 'all tests passed'

