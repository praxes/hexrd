import sys, os, time

import numpy as np

from hexrd import matrixutil as mutil

from hexrd.xrd import rotations as rot
from hexrd.xrd import transforms_CAPI as xfcapi

n_quats = 1000000
rq = mutil.unitVector(np.random.randn(4, n_quats))
quats = np.array(rq.T, dtype=float, order='C')

# phi = 2. * rot.arccosSafe(rq[0, :])
# n   = np.tile(1. / np.sin(0.5*phi), (3, 1)) * rq[1:, :]

start0 = time.clock()                      # time this
rMats0 = rot.rotMatOfQuat(rq)
elapsed0 = (time.clock() - start0)

# rMats1 = np.zeros((n_quats, 3, 3))
# for i in range(n_quats):
#     rMats1[i, :, :] = xfcapi.makeRotMatOfQuat(quats[i, :])

start1 = time.time()                      # time this
rMats1 = xfcapi.makeRotMatOfQuat(quats)
elapsed1 = (time.time() - start1)
print "Time for %d quats:\t%g v. %g (%f)"%(n_quats, elapsed0/float(n_quats), elapsed1/float(n_quats), elapsed0/elapsed1)
print "Maximum discrepancy:\t%f" % (np.amax(abs(rMats0 - rMats1)))
