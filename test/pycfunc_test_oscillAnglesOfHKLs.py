import sys, os, time
import numpy as np

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import pycfuncs_transforms as pycfuncs

bVec_ref    = np.ascontiguousarray(xf.bVec_ref)
eta_ref     = np.ascontiguousarray(xf.eta_ref)

idxFile = './ruby_4537-8_log.txt'

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
start1 = time.clock()                      # time this
for i in range(n):
    oangs01, oangs11 = xf.oscillAnglesOfHKLs(hklsT, chi, rMat_c, bMat, wavelength, 
                                             beamVec=bVec_ref, etaVec=eta_ref)
elapsed1 = (time.clock() - start1)
print "Time for Python oscillAnglesOfHKLs: %f"%elapsed1

start2 = time.clock()                      # time this
for i in range(n):
    oangs02, oangs12 = xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
    #xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
elapsed2 = (time.clock() - start2)
print "Time for CAPI oscillAnglesOfHKLs:   %f"%elapsed2
#print oangs01.shape, oangs11.shape
#print oangs02.shape, oangs12.shape
#print np.linalg.norm(oangs01[:,0]),np.linalg.norm(oangs01[:,1]),np.linalg.norm(oangs01[:,2])
#print np.linalg.norm(oangs11[:,0]),np.linalg.norm(oangs11[:,1]),np.linalg.norm(oangs11[:,2])
#print "Maximum Relative Differences: %f, %f"%(np.linalg.norm(oangs01-oangs02.T)/np.linalg.norm(oangs01),np.linalg.norm(oangs11-oangs12.T)/np.linalg.norm(oangs11))

#print "  Speedup: %f"%(elapsed1/elapsed2)


start3 = time.clock()                      # time this
for i in range(n):
    oangs03, oangs13 = pycfuncs.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
    #xfcapi.oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength, beamVec=bVec_ref, etaVec=eta_ref)
elapsed3 = (time.clock() - start3)
print "Time for Numba PyCFuncs oscillAnglesOfHKLs:   %f" %elapsed3
#print oangs01.shape, oangs11.shape
#print oangs02.shape, oangs12.shape
#print np.linalg.norm(oangs01[:,0]),np.linalg.norm(oangs01[:,1]),np.linalg.norm(oangs01[:,2])
#print np.linalg.norm(oangs11[:,0]),np.linalg.norm(oangs11[:,1]),np.linalg.norm(oangs11[:,2])
#print "Maximum Relative Differences: %f, %f"%(np.linalg.norm(oangs01-oangs02.T)/np.linalg.norm(oangs01),np.linalg.norm(oangs11-oangs12.T)/np.linalg.norm(oangs11))

print 'testing data ...'
#test oangs


assert np.allclose(oangs01.T, oangs02), 'oangs01 and oangs02 not close'
assert np.allclose(oangs02, oangs03),  'oangs02 and oangs03 not close'
assert np.allclose(oangs11.T, oangs12), 'oangs11 and oangs12 not close'
assert np.allclose(oangs12, oangs13),  'oangs12 and oangs13 not close'
print 'all tests passed'
