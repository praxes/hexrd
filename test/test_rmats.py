import sys, os, time

import numpy as np

from hexrd import matrixutil as mutil

from hexrd.xrd import rotations as rot
from hexrd.xrd import transforms_CAPI as xfcapi

n_quats = 1e6
rq = mutil.unitVector(np.random.randn(4, n_quats))
quats = np.array(rq.T, dtype=float, order='C')

""" NUMPY """
start0 = time.clock()                      # time this
rMats0 = rot.rotMatOfQuat(rq)
elapsed0 = (time.clock() - start0)

""" CAPI """
start1 = time.time()                      # time this
rMats1 = xfcapi.makeRotMatOfQuat(quats)
# rMats1 = np.zeros((n_quats, 3, 3))
# for i in range(n_quats):
#     rMats1[i, :, :] = xfcapi.makeRotMatOfQuat(quats[i, :])
elapsed1 = (time.time() - start1)
print "Time for %d quats:\t%g v. %g (%f)"%(n_quats, elapsed0/float(n_quats), elapsed1/float(n_quats), elapsed0/elapsed1)
print "Maximum discrepancy:\t%f" % (np.amax(abs(rMats0 - rMats1)))

""" NOTES

If I call xfcapi.makeRotMatOfQuat(quats) I get:

[bernier2@photonjb2] test > python -i test_rmats.py 
Time for 10000000 quats: 3.90575e-07 v. 1.85224e-07 (2.108666)
Maximum discrepancy:    0.000000

Strange, but if I call xfcapi.makeRotMatOfQuat(rq.T) instead, I get:

[bernier2@photonjb2] test > python -i test_rmats.py 
Time for 10000000 quats: 3.88884e-07 v. 2.22061e-07 (1.751246)
Maximum discrepancy:    0.000000
"""
