from timeit import default_timer as timer
import sys, os
import numpy as np

import numba.cuda

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import pycfuncs_transforms as pycfuncs

bVec_ref    = np.ascontiguousarray(xf.bVec_ref)
eta_ref     = np.ascontiguousarray(xf.eta_ref)

#idxFile = './ruby_4537-8_log.txt'

idxFile = './ruby_triple.txt'
gtable  = np.loadtxt(idxFile, delimiter='\t')
idx     = gtable[:, 0] >= 0
hklsT   = np.ascontiguousarray(gtable[idx, 2:5].T)
hkls    = np.ascontiguousarray(gtable[idx, 2:5])

# input parameters
wavelength = 0.153588                     # Angstroms (80.725keV)

chi    = -0.0011591608938627839

bMat = np.array( [ [  2.10048731e-01,   0.00000000e+00,   0.00000000e+00],
                   [  1.21271692e-01,   2.42543383e-01,   0.00000000e+00],
                   [  0.00000000e+00,   0.00000000e+00,   7.69486476e-02] ] )

rMat_c = xf.makeRotMatOfExpMap(np.array( [ [ 0.66931818],
                                           [-0.98578066],
                                           [ 0.73593251] ] ) )

# ######################################################################

# oscillation angle arrays
n=22
start1 = timer()                      # time this
for i in range(n):
    oangs01, oangs11 = xf.oscillAnglesOfHKLs(hklsT, chi, rMat_c, bMat, wavelength, 
                                             beamVec=bVec_ref, etaVec=eta_ref)
elapsed1 = (timer() - start1)
print "Time for Python oscillAnglesOfHKLs: %f"%elapsed1

start2 = timer()                      # time this
for i in range(n):
    oangs02, oangs12 = xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
    #xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
elapsed2 = (timer() - start2)
print "Time for CAPI oscillAnglesOfHKLs:   %f"%elapsed2
#print oangs01.shape, oangs11.shape
#print oangs02.shape, oangs12.shape
#print np.linalg.norm(oangs01[:,0]),np.linalg.norm(oangs01[:,1]),np.linalg.norm(oangs01[:,2])
#print np.linalg.norm(oangs11[:,0]),np.linalg.norm(oangs11[:,1]),np.linalg.norm(oangs11[:,2])
#print "Maximum Relative Differences: %f, %f"%(np.linalg.norm(oangs01-oangs02.T)/np.linalg.norm(oangs01),np.linalg.norm(oangs11-oangs12.T)/np.linalg.norm(oangs11))

#print "  Speedup: %f"%(elapsed1/elapsed2)


print 'cudadevice: ', numba.cuda.get_current_device().name

gVec_e = np.zeros(3)
gHat_c = np.zeros(3)
gHat_s = np.zeros(3)
bHat_l = np.zeros(3)
eHat_l = np.zeros(3) 
oVec = np.zeros(2)
tVec0 = np.zeros(3)
rMat_e = np.zeros(9)
rMat_s = np.zeros(9)
npts = hkls.shape[0]
#return arrays
oangs0 = np.zeros((npts, 3))
oangs1 = np.zeros((npts, 3))


start3 = timer()                      # time this
#for i in range(n):
pycfuncs.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, bVec_ref, eta_ref,
                            gVec_e, gHat_c, gHat_s, 
                            bHat_l, eHat_l, oVec, tVec0, 
                            rMat_e, rMat_s, npts,
                            oangs0, oangs1)
    #xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
elapsed3 = (timer() - start3)
print "Time for Numba PyCFuncs oscillAnglesOfHKLs:   %f" %elapsed3
oangs03 = oangs0
oangs13 = oangs1

#print oangs01.shape, oangs11.shape
#print oangs02.shape, oangs12.shape
#print np.linalg.norm(oangs01[:,0]),np.linalg.norm(oangs01[:,1]),np.linalg.norm(oangs01[:,2])
#print np.linalg.norm(oangs11[:,0]),np.linalg.norm(oangs11[:,1]),np.linalg.norm(oangs11[:,2])
#print "Maximum Relative Differences: %f, %f"%(np.linalg.norm(oangs01-oangs02.T)/np.linalg.norm(oangs01),np.linalg.norm(oangs11-oangs12.T)/np.linalg.norm(oangs11))

print 'testing data ...'
#test oangs
#print oangs02
#print '*' * 50
#print oangs03

assert np.allclose(oangs01.T, oangs02), 'oangs01 and oangs02 not close'
assert np.allclose(oangs02, oangs03),  'oangs02 and oangs03 not close'
assert np.allclose(oangs11.T, oangs12), 'oangs11 and oangs12 not close'
assert np.allclose(oangs12, oangs13),  'oangs12 and oangs13 not close'
print 'all tests passed'
