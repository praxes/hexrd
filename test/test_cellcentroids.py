#
# This script runs alternative cellCentroid implementations against
# the one in the hexrd. package. The results are compared against the
# original using np.allclose() to ensure correctness. Timing
# information is also gathered used the numbapro.nvtx explicit
# instrumentation profiler (which happens to be integrated with the
# nvidia nvvp timeline profiler application).
#
# In order to run this script, a data file with the arguments for the
# cellCentroid function must be passed in (defaults to DEFAULT_ARG_FILE
# below).
#
# The arg file is just a series of pickled tuples containing the
# arguments to the calls. The best way to build that file is running a
# dataset in hexrd with the cellCentroids modified so that for each call it
# pickles out its parameters. An example file is NOT included due to its
# relative big size.
#


from __future__ import print_function

import os, sys, time

import cPickle as pickle
from scipy import optimize as opt
import numpy as np
import numba
from numbapro import nvtx
# import the reference implementation
from hexrd.gridutil import cellCentroids


DEFAULT_ARG_FILE = 'cell_centroids_args.pickled'

def cell_centroids_original(crd, con):
    """
    con.shape = (nele, 4)
    crd.shape = (ncrd, 2)

    con.shape = (nele, 8)
    crd.shape = (ncrd, 3)
    """
 
    nele = con.shape[0]
    dim  = crd.shape[1]
    centroid_xy = np.zeros((nele, dim))
    for i in range(len(con)):
        el_crds = crd[con[i, :], :] # (4, 2)
        centroid_xy[i, :] = (el_crds).mean(axis=0)
    return centroid_xy

@numba.njit
def _cell_centroids_opt1(crd, con, out):
    nele, conn_count = con.shape
    dim = crd.shape[1]
    inv_conn = 1.0/conn_count
    for i in range(nele):
        for j in range(dim):
            acc = 0.0
            for k in range(conn_count):
                acc += crd[con[i,k], j]
            out[i,j] = acc * inv_conn
    return out

def cell_centroids_opt1(crd, con):
    nele = con.shape[0]
    dim = crd.shape[1]
    result_xy = np.empty((nele, dim))
    return _cell_centroids_opt1(crd, con, result_xy)

def test_cellcentroids(in_file):
    count = 0
    with open(in_file, 'rb') as f:
        while True:
            try:
                args = pickle.load(f)
            except EOFError:
                break

            print(args[0].shape, args[0].dtype, args[1].shape, args[1].dtype)
            with nvtx.profile_range('original'):
                res_orig = cell_centroids_original(*args)
            with nvtx.profile_range('current in hexrd'):
                res_current = cellCentroids(*args)
            with nvtx.profile_range('numba1'):
                res_numba1 = cell_centroids_opt1(*args)

            assert np.allclose(res_current, res_orig)
            assert np.allclose(res_current, res_numba1)
            count += 1


if __name__=='__main__':
    in_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ARG_FILE
    test_cellcentroids(in_file)

    print(" STATS ".center(72, '='))
    fmt = "{2:>14}, {1:>8}, {0:<40}"
    print(fmt.format("FUNCTION", "CALLS", "TIME"))
    fmt = "{2:>14F}, {1:>8}, {0:<40}"
    sorted_by_time = sorted(nvtx.getstats().iteritems(), key=lambda tup: tup[1][1])
    for key, val in sorted_by_time:
        print(fmt.format(key, *val))
