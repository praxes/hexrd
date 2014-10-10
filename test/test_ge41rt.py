#
# This script runs alternative GE_41RT implementation against the one
# in the hexrd.xrd package. The results are compared against the
# original using np.allclose() to ensure correctness. Timing information
# is also gathered used the numbapro.nvtx explicit instrumentation
# profiler (which happens to be integrated with the nvidia nvvp timeline
# profiler application).
#
# In order to run this script, a data file with the arguments for the
# GE_41RT function must be passed in (defaults to DEFAULT_ARG_FILE
# below).
#
# The arg file is just a series of pickled tuples containing the
# arguments to the calls. The best way to build that file is running a
# dataset in hexrd with the GE_41RT modified so that for each call it
# pickles out its parameters. An example file is NOT include due to its
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
from hexrd.xrd.distortion import GE_41RT


DEFAULT_ARG_FILE = 'ge_41rt_args'

@numba.njit
def GE_41RT_rho_scl_func_inv_opt(ri, ni, ro, rx, p, res):
    ri = ri[0]
    ratio = ri/rx
    a = p[0]*ratio**p[3] * np.cos(2.0 * ni)
    b = p[1]*ratio**p[4] * np.cos(4.0 * ni)
    c = p[2]*ratio**p[5]

    return (a + b + c + 1.0)*ri - ro


@numba.njit
def GE_41RT_rho_scl_fi_prime_opt(ri, ni, ro, rx, p, res):
    ri = ri[0]
    ratio = ri/rx
    res[0] = p[0]*ratio**p[3] * np.cos(2.0 * ni) * (p[3] + 1) + \
             p[1]*ratio**p[4] * np.cos(4.0 * ni) * (p[4] + 1) + \
             p[2]*ratio**p[5] * (p[5] + 1) + 1
    return res



@numba.guvectorize(['void(float64[:], float64[:])'], '(n)->(n)')
def to_polar(_in, _out):
    _out[0] = np.sqrt(_in[0]*_in[0] + _in[1]*_in[1])
    _out[1] = np.arctan2(_in[1], _in[0])


def GE_41RT_opt(xy_in, params, invert=False):
    """
    Apply radial distortion to polar coordinates on GE detector

    xin, yin are 1D arrays or scalars, assumed to be relative to self.xc, self.yc
    Units are [mm, radians].  This is the power-law based function of Bernier.

    Available Keyword Arguments :

    invert = True or >False< :: apply inverse warping
    """

    if params[0] == 0 and params[1] == 0 and params[2] == 0:
        return xy_in
    else:
        # canonical max radius based on perfectly centered beam
        rhoMax = 204.8

        polar = to_polar(xy_in)
        npts = len(polar)

        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = polar[:,0].flatten()
        eta0 = polar[:,1].flatten()

        if invert:
            # in here must do nonlinear solve for distortion
            # must loop to call fsolve individually for each point
            rhoOut = np.zeros(npts, dtype=float)

            rhoSclFuncInv = GE_41RT_rho_scl_func_inv_opt
            rhoSclFIprime = GE_41RT_rho_scl_fi_prime_opt
            res = np.array((1,), dtype=np.double)
            for iRho in range(len(rho0)):
                rhoOut[iRho] = opt.fsolve(rhoSclFuncInv, rho0[iRho],
                                          fprime=rhoSclFIprime,
                                          args=(eta0[iRho], rho0[iRho], rhoMax, params, res) )
                pass
        else:
            # usual case: calculate scaling to take you from image to detector plane
            # 1 + p[0]*(ri/rx)**p[2] * np.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
            rhoSclFunc = lambda ri, rx=rhoMax, p=params, ni=eta0: \
                         p[0]*(ri/rx)**p[3] * np.cos(2.0 * ni) + \
                         p[1]*(ri/rx)**p[4] * np.cos(4.0 * ni) + \
                         p[2]*(ri/rx)**p[5] + 1

            rhoOut = np.squeeze( rho0 * rhoSclFunc(rho0) )
            pass

        xout = rhoOut * np.cos(eta0)
        yout = rhoOut * np.sin(eta0)
    return np.vstack([xout, yout]).T


def newton(x0, f, fp, extra, prec=1e-16, maxiter=100):
    for i in range(maxiter):
        x = x0 - f(x0, *extra) / fp(x0, *extra)
        if np.max(np.abs(x - x0)) < prec:
            return x
        x0 = x
    return x0


@numba.njit
def GE_41RT_rho_scl_func_inv_newton(ri, ni, ro, rx, p):
    ratio = ri/rx
    a = p[0]*ratio**p[3] * np.cos(2.0 * ni)
    b = p[1]*ratio**p[4] * np.cos(4.0 * ni)
    c = p[2]*ratio**p[5]

    return (a + b + c + 1.0)*ri - ro


@numba.njit
def GE_41RT_rho_scl_fi_prime_newton(ri, ni, ro, rx, p):
    ratio = ri/rx
    res = p[0]*ratio**p[3] * np.cos(2.0 * ni) * (p[3] + 1) + \
          p[1]*ratio**p[4] * np.cos(4.0 * ni) * (p[4] + 1) + \
          p[2]*ratio**p[5] * (p[5] + 1) + 1
    return res



def GE_41RT_newton(xy_in, params, invert=False):
    """
    Apply radial distortion to polar coordinates on GE detector

    xin, yin are 1D arrays or scalars, assumed to be relative to self.xc, self.yc
    Units are [mm, radians].  This is the power-law based function of Bernier.

    Available Keyword Arguments :

    invert = True or >False< :: apply inverse warping
    """

    if params[0] == 0 and params[1] == 0 and params[2] == 0:
        return xy_in
    else:
        # canonical max radius based on perfectly centered beam
        rhoMax = 204.8

        polar = to_polar(xy_in)
        npts = len(polar)

        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = polar[:,0].flatten()
        eta0 = polar[:,1].flatten()

        if invert:
            # in here must do nonlinear solve for distortion
            # must loop to call fsolve individually for each point
            rhoOut = np.zeros(npts, dtype=float)

            rhoSclFuncInv = GE_41RT_rho_scl_func_inv_newton
            rhoSclFIprime = GE_41RT_rho_scl_fi_prime_newton
            for iRho in range(len(rho0)):
                rhoOut[iRho] = newton(float(rho0[iRho]), rhoSclFuncInv, rhoSclFIprime,
                                      (float(eta0[iRho]), float(rho0[iRho]), rhoMax, params) )
                pass
        else:
            # usual case: calculate scaling to take you from image to detector plane
            # 1 + p[0]*(ri/rx)**p[2] * np.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
            rhoSclFunc = lambda ri, rx=rhoMax, p=params, ni=eta0: \
                         p[0]*(ri/rx)**p[3] * np.cos(2.0 * ni) + \
                         p[1]*(ri/rx)**p[4] * np.cos(4.0 * ni) + \
                         p[2]*(ri/rx)**p[5] + 1

            rhoOut = np.squeeze( rho0 * rhoSclFunc(rho0) )
            pass

        xout = rhoOut * np.cos(eta0)
        yout = rhoOut * np.sin(eta0)
    return np.vstack([xout, yout]).T


@numba.njit
def inverse_distortion(rhoOut, rho0, eta0, rhoMax, params):
    maxiter = 100
    prec = 1e-16

    p0, p1, p2, p3, p4, p5 = params[0:6]
    rx = rhoMax
    for el in range(len(rhoOut)):
        ri = rho0[el]
        ni = eta0[el]
        ro = ri
        cos2ni = np.cos(2.0*ni)
        cos4ni = np.cos(4.0*ni)
        for i in range(maxiter):
            ratio = ri/rx
            fx = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri - ro
            fxp = (p0*ratio**p3*cos2ni*(p3+1) +
                   p1*ratio**p4*cos4ni*(p4+1) +
                   p2*ratio**p5*(p5+1) + 1)

            delta = fx/fxp
            ri = ri - delta
            if np.abs(delta) < prec:
                break

        rhoOut[el] = ri

def GE_41RT_newton_numba(xy_in, params, invert=False):
    """
    Apply radial distortion to polar coordinates on GE detector

    xin, yin are 1D arrays or scalars, assumed to be relative to self.xc, self.yc
    Units are [mm, radians].  This is the power-law based function of Bernier.

    Available Keyword Arguments :

    invert = True or >False< :: apply inverse warping
    """

    if params[0] == 0 and params[1] == 0 and params[2] == 0:
        return xy_in
    else:
        # canonical max radius based on perfectly centered beam
        rhoMax = 204.8

        polar = to_polar(xy_in)
        npts = len(polar)

        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = polar[:,0].flatten()
        eta0 = polar[:,1].flatten()

        if invert:
            # in here must do nonlinear solve for distortion
            # must loop to call fsolve individually for each point
            rhoOut = np.empty(npts, dtype=float)
            inverse_distortion(rhoOut, rho0, eta0, rhoMax, params)
        else:
            # usual case: calculate scaling to take you from image to detector plane
            # 1 + p[0]*(ri/rx)**p[2] * np.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
            rhoSclFunc = lambda ri, rx=rhoMax, p=params, ni=eta0: \
                         p[0]*(ri/rx)**p[3] * np.cos(2.0 * ni) + \
                         p[1]*(ri/rx)**p[4] * np.cos(4.0 * ni) + \
                         p[2]*(ri/rx)**p[5] + 1

            rhoOut = np.squeeze( rho0 * rhoSclFunc(rho0) )
            pass

        xout = rhoOut * np.cos(eta0)
        yout = rhoOut * np.sin(eta0)
    return np.vstack([xout, yout]).T


def test_ge41rt(in_file):
    _colors = [nvtx.colors.cyan, nvtx.colors.magenta,
               nvtx.colors.yellow, nvtx.colors.black]
    count = 0
    with open(in_file, 'rb') as f:
        while True:
            try:
                args = pickle.load(f)
            except EOFError:
                break
            with nvtx.profile_range('ge_41rt',
                                    color=_colors[(count*2)%len(_colors)],
                                    payload=args[1].shape[0]):
                res = GE_41RT(*args)
            with nvtx.profile_range('ge_41rt_opt',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_opt = GE_41RT_opt(*args)

            with nvtx.profile_range('ge_41rt_newton',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton = GE_41RT_newton(*args)

            with nvtx.profile_range('ge_41rt_newton_numba',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba = GE_41RT_newton_numba(*args)

            assert np.allclose(res, res_opt)
            assert np.allclose(res, res_newton)
            assert np.allclose(res, res_newton_numba)
            count += 1


if __name__=='__main__':
    in_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ARG_FILE
    test_ge41rt(in_file)

    print(" STATS ".center(72, '='))
    fmt = "{2:>14}, {1:>8}, {0:<40}"
    print(fmt.format("FUNCTION", "CALLS", "TIME"))
    fmt = "{2:>14F}, {1:>8}, {0:<40}"
    sorted_by_time = sorted(nvtx.getstats().iteritems(), key=lambda tup: tup[1][1])
    for key, val in sorted_by_time:
        print(fmt.format(key, *val))
