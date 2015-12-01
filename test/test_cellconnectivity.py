#
# This script runs alternative cellConnectivity implementations against the one
# in the hexrd. package. The results are compared against the original using
# np.allclose() to ensure correctness. Timing information is also gathered used
# the nvtxpy explicit instrumentation profiler (which happens to be integrated
# with the nvidia nvvp timeline profiler application).
#
# In order to run this script, a data file with the arguments for the
# cellConnectivity function must be passed in (defaults to DEFAULT_ARG_FILE
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
import nvtxpy as nvtx
# import the reference implementation
from hexrd.gridutil import cellConnectivity

DEFAULT_ARG_FILE = 'cell_connectivity_args.pickled'


def cell_connectivity_original(m, n, p=1, origin='ul'):
    """
    p x m x n (layers x rows x cols)

    origin can be upper left -- 'ul' <default> or lower left -- 'll'

    choice will affect handedness (cw or ccw)
    """
    nele = p*m*n
    con  = np.zeros((nele, 4), dtype=int)
    i_con = 0
    for k in range(p):
        for j in range(m):
            for i in range(n):
                con[i_con, :] = np.array([ i + j*(n + 1) + 1, 
                                           i + j*(n + 1), 
                                           i + j + n*(j + 1) + 1,
                                           i + j + n*(j + 1) + 2 ]) + k*(n + 1)*(m + 1)
                i_con += 1
    if p > 1:
        nele = m*n*(p-1)
        tmp_con3 = con.reshape(p, m*n, 4)
        hex_con = []
        for layer in range(p - 1):
            hex_con.append(np.hstack([tmp_con3[layer], tmp_con3[layer + 1]]))
        con = vstack(hex_con)
        pass
    if origin.lower().strip() == 'll':
        con = con[:, ::-1]
    return con


@numba.njit
def _fill_connectivity(out, m, n, p):
    i_con = 0
    for k in range(p):
        for j in range(m):
            for i in range(n):
                extra = k*(n+1)*(m+1)
                out[i_con, 0] = i + j*(n + 1) + 1 + extra
                out[i_con, 1] = i + j*(n + 1) + extra
                out[i_con, 2] = i + j + n*(j+1) + 1 + extra
                out[i_con, 3] = i + j + n*(j+1) + 2 + extra
                i_con += 1


def cell_connectivity_numba(m, n, p=1, origin='ul'):
    """
    p x m x n (layers x rows x cols)

    origin can be upper left -- 'ul' <default> or lower left -- 'll'

    choice will affect handedness (cw or ccw)
    """
    nele = p*m*n
    con  = np.empty((nele, 4), dtype=int)
    _fill_connectivity(con, m, n, p)

    if p > 1:
        nele = m*n*(p-1)
        tmp_con3 = con.reshape(p, m*n, 4)
        hex_con = []
        for layer in range(p - 1):
            hex_con.append(np.hstack([tmp_con3[layer], tmp_con3[layer + 1]]))
        con = vstack(hex_con)
        pass
    if origin.lower().strip() == 'll':
        con = con[:, ::-1]
    return con


def run_test(in_file):
    count = 0
    with open(in_file, 'rb') as f:
        while True:
            try:
                args = pickle.load(f)
            except EOFError:
                break

            with nvtx.profile_range('original'):
                res_orig = cell_connectivity_original(*args)
            with nvtx.profile_range('current in hexrd'):
                res_current = cellConnectivity(*args)
            with nvtx.profile_range('numba'):
                res_numba = cell_connectivity_numba(*args)
            with nvtx.profile_range('numba_bis'):
                res_numba = cell_connectivity_numba(*args)

            assert np.allclose(res_current, res_orig)
            assert np.allclose(res_current, res_numba)
            count += 1


if __name__=='__main__':
    in_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ARG_FILE
    run_test(in_file)

    print(" STATS ".center(72, '='))
    fmt = "{2:>14}, {1:>8}, {0:<40}"
    print(fmt.format("FUNCTION", "CALLS", "TIME"))
    fmt = "{2:>14F}, {1:>8}, {0:<40}"
    sorted_by_time = sorted(nvtx.getstats().iteritems(), key=lambda tup: tup[1][1])
    for key, val in sorted_by_time:
        print(fmt.format(key, *val))
