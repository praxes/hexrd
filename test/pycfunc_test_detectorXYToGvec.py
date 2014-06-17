from timeit import default_timer as timer
import sys, os, time
import numpy as np

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import pycfuncs_transforms as pycfuncs 

# input parameters
bVec_ref = xf.bVec_ref

rMat_d = xf.makeDetectorRotMat( ( 0.0011546340766314521,
                                 -0.0040527538387122993,
                                 -0.0026221336905160211 ) ) 
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
#pvec  = 204.8 * np.linspace(-1, 1, 2048)
pvec  = 204.8 * np.linspace(-1, 1, 512)
dcrds = np.meshgrid(pvec, pvec)
XY    = np.ascontiguousarray(np.vstack([dcrds[0].flatten(), dcrds[1].flatten()]).T)

# Check the timings
start1 = timer()
dangs1 = xf.detectorXYToGvec(XY, rMat_d, rMat_s,
                             tVec_d, tVec_s, tVec_c, 
                             beamVec=bVec_ref)
tTh_d1 = dangs1[0][0]
eta_d1 = dangs1[0][1]
gVec1  = dangs1[1]
elapsed1 = (timer() - start1)
print "Time for Python detectorXYToGvec: %f"%(elapsed1)

# start2 = timer()                      # time this
# dangs2 = xfc.detectorXYToGvecC(XY, rMat_d, rMat_s,
#                                tVec_d, tVec_s, tVec_c, 
#                                beamVec=bVec_ref.flatten())
# tTh_d2 = dangs2[0][0]
# eta_d2 = dangs2[0][1]
# gVec2  = dangs2[1]
# elapsed2 = (timer() - start2)
# print "Time for Cython detectorXYToGvec: %f"%(elapsed2)

#maxDiff_tTh = np.linalg.norm(tTh_d1-tTh_d2,np.inf)
#print "Maximum disagreement in tTh:  %f"%maxDiff_tTh
#maxDiff_eta = np.linalg.norm(eta_d1-eta_d2,np.inf)
#print "Maximum disagreement in eta:  %f"%maxDiff_eta
#maxDiff_gVec = np.linalg.norm(np.sqrt(np.sum(np.asarray(gVec1.T-gVec2)**2,1)),np.inf)
#print "Maximum disagreement in gVec: %f"%maxDiff_gVec

start3 = timer()                      # time this
dangs3 = xfcapi.detectorXYToGvec(XY, rMat_d, rMat_s,
                                 tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(), 
                                 beamVec=bVec_ref.flatten(),etaVec=np.array([1.0,0.0,0.0]))
tTh_d3 = dangs3[0][0]
eta_d3 = dangs3[0][1]
gVec3  = dangs3[1]
elapsed3 = (timer() - start3)
print "Time for CAPI detectorXYToGvec: %f"%(elapsed3)

maxDiff_tTh = np.linalg.norm(tTh_d1-tTh_d3,np.inf)
print "Maximum disagreement in tTh:  %f"%maxDiff_tTh
maxDiff_eta = np.linalg.norm(eta_d1-eta_d3,np.inf)
print "Maximum disagreement in eta:  %f"%maxDiff_eta
maxDiff_gVec = np.linalg.norm(np.sqrt(np.sum(np.asarray(gVec1.T-gVec3)**2,1)),np.inf)
print "Maximum disagreement in gVec: %f"%maxDiff_gVec

#setup for detectorXYToGVec
rMat_e = np.zeros(9)
bVec = np.zeros(3)
tVec1 = np.zeros(3)
tVec2 = np.zeros(3)
dHat_l = np.zeros(3)
n_g = np.zeros(3)
npts = XY.shape[0]
#return values
tTh = np.zeros(npts)
eta = np.zeros(npts)
gVec_l = np.zeros((npts, 3))


start4 = timer()                      # time this



pycfuncs.detectorXYToGvec(XY, rMat_d, rMat_s,
                                 tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(), 
                                 bVec_ref.flatten(),np.array([1.0,0.0,0.0]),
                                 rMat_e, bVec, tVec1, tVec2, dHat_l, n_g, npts, tTh, eta, gVec_l)
                                 
tTh_d4 = tTh  #dangs4[0][0]
eta_d4 = eta  # dangs4[0][1]
gVec4  = gVec_l  # dangs4[1]
elapsed4 = (timer() - start4)
print "Time for Numba PyCFuncs detectorXYToGvec: %f"%(elapsed4)
maxDiff_tTh = np.linalg.norm(tTh_d3 - tTh_d4,np.inf)
print "Maximum disagreement in tTh:  %f"%maxDiff_tTh
maxDiff_eta = np.linalg.norm(eta_d3 - eta_d4,np.inf)
print "Maximum disagreement in eta:  %f"%maxDiff_eta
maxDiff_gVec = np.linalg.norm(np.sqrt(np.sum(np.asarray(gVec3 - gVec4)**2,1)),np.inf)
print "Maximum disagreement in gVec: %f"%maxDiff_gVec




print 'testing data ...'
#test tTh
assert np.allclose(tTh_d1, tTh_d3), 'tTh_d1 not close to tTh_d3'
#assert np.array_equal(tTh_d1, tTh_d3),  'tTh_d1 not equal to tTh_d3'
assert np.allclose(tTh_d3, tTh_d4), 'C and PyC not close: tTh_d3, tTh_d4'
#test eta
assert np.allclose(eta_d1, eta_d3), 'eta_d1 not close to eta_d3'
assert np.allclose(eta_d3, eta_d4), 'C and PyC not close: eta_d3, eta_d4'
#test gVec 
assert np.allclose(gVec1.T, gVec3), 'eta_d1 not close to eta_d3'
assert np.allclose(gVec3, gVec4), 'C and PyC not close: eta_d3, eta_d4'
print 'all tests passed'
