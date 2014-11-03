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

@numba.njit
def to_polar(_in_out):
    "inplace conversion from rectangular coordinates to polar"
    for i in range(len(_in_out)):
        x, y = _in_out[i, 0:2]
        _in_out[i, 0] = np.sqrt(x*x + y*y)
        _in_out[i, 1] = np.arctan2(y,x)
    return _in_out

@numba.njit
def to_rectangular(_in_out):
    "inplace conversion form polar coordinates to rectangular"
    for i in range(len(_in_out)):
        rho, theta = _in_out[i, 0:2]
        _in_out[i, 0] = rho * np.cos(theta)
        _in_out[i, 1] = rho * np.sin(theta)
    return _in_out

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

        polar = to_polar(np.copy(xy_in, order='F'))
        npts = len(polar)
        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = polar[:,0]
        eta0 = polar[:,1]

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
        if np.max(np.abs(x - x0)) <= prec:
            return x
        x0 = x
    return x0


def GE_41RT_rho_scl_func_inv_newton(ri, ni, ro, rx, p):
    ratio = ri/rx
    a = p[0]*ratio**p[3] * np.cos(2.0 * ni)
    b = p[1]*ratio**p[4] * np.cos(4.0 * ni)
    c = p[2]*ratio**p[5]

    return (a + b + c + 1.0)*ri - ro


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

        polar = to_polar(np.copy(xy_in, order='F'))
        npts = len(polar)

        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = polar[:,0]
        eta0 = polar[:,1]

        if invert:
            # in here must do nonlinear solve for distortion
            # must loop to call fsolve individually for each point
            rhoOut = np.zeros(npts, dtype=float)

            rhoSclFuncInv = GE_41RT_rho_scl_func_inv_newton
            rhoSclFIprime = GE_41RT_rho_scl_fi_prime_newton
            rhoOut = newton(rho0, rhoSclFuncInv, rhoSclFIprime,
                            (eta0, rho0, rhoMax, params))
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
def inverse_distortion(polar, rhoMax, params):
    maxiter = 100
    prec = 1e-16

    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax
    for el in range(len(polar)):
        ri, ni = polar[el, 0:2]
        ro = ri
        cos2ni = np.cos(2.0*ni)
        cos4ni = np.cos(4.0*ni)
        for i in range(maxiter):
            ratio = ri*rxi
            fx = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri - ro
            fxp = (p0*ratio**p3*cos2ni*(p3+1) +
                   p1*ratio**p4*cos4ni*(p4+1) +
                   p2*ratio**p5*(p5+1) + 1)

            delta = fx/fxp
            ri = ri - delta
            if np.abs(delta) <= prec*np.abs(ri):
                break

        polar[el, 0] = ri


@numba.njit
def direct_distortion(polar, rhoMax, params):
    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax
    for el in range(len(polar)):
        ri, ni = polar[el, 0:2]
        cos2ni = np.cos(2.0*ni)
        cos4ni = np.cos(4.0*ni)
        ratio = ri*rxi

        polar[el,0] = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri


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
        polar = to_polar(np.copy(xy_in))

        if invert:
            # in here must do nonlinear solve for distortion
            # must loop to call fsolve individually for each point
            inverse_distortion(polar, rhoMax, params)
        else:
            # usual case: calculate scaling to take you from image to detector plane
            # 1 + p[0]*(ri/rx)**p[2] * np.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
            direct_distortion(polar, rhoMax, params)

        return to_rectangular(polar)




@numba.njit
def inverse_distortion_full(out, in_, rhoMax, params):
    maxiter = 100
    prec = 1e-16

    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax
    for el in range(len(in_)):
        xi, yi = in_[el, 0:2]
        ri = np.sqrt(xi*xi + yi*yi)
        ri_inv = 1.0/ri
        sinni = yi*ri_inv
        cosni = xi*ri_inv
        ro = ri
        cos2ni = cosni*cosni - sinni*sinni
        sin2ni = 2*sinni*cosni
        cos4ni = cos2ni*cos2ni - sin2ni*sin2ni
        for i in range(maxiter):
            ratio = ri*rxi
            fx = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri - ro
            fxp = (p0*ratio**p3*cos2ni*(p3+1) +
                   p1*ratio**p4*cos4ni*(p4+1) +
                   p2*ratio**p5*(p5+1) + 1)

            delta = fx/fxp
            ri = ri - delta
            if np.abs(delta) <= prec*np.abs(ri):
                break
            
        xi = ri*cosni
        yi = ri*sinni
        out[el, 0] = xi
        out[el, 1] = yi

    return out


@numba.njit
def direct_distortion_full(out, in_, rhoMax, params):
    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax

    for el in range(len(in_)):
        xi, yi = in_[el, 0:2]
        ri = np.sqrt(xi*xi + yi*yi)
        ri_inv = 1.0/ri
        sinni = yi*ri_inv
        cosni = xi*ri_inv
        cos2ni = cosni*cosni - sinni*sinni
        sin2ni = 2*sinni*cosni
        cos4ni = cos2ni*cos2ni - sin2ni*sin2ni
        ratio = ri*rxi

        ri = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri
        xi = ri*cosni
        yi = ri*sinni
        out[el, 0] = xi
        out[el, 1] = yi

    return out


def GE_41RT_newton_numba_full(xy_in, params, invert=False):
    if params[0] == 0 and params[1] == 0 and params[2] == 0:
        return xy_in
    else:
        rhoMax = 204.8
        if invert:
            xy_out = np.empty_like(xy_in)
            return inverse_distortion_full(xy_out, xy_in, rhoMax, params)
        else:
            xy_out = np.empty_like(xy_in)
            return direct_distortion_full(xy_out, xy_in, rhoMax, params)

@numba.njit('f8[:,:](f8[:,:], f8[:,:], f8, f8[6])')
def inverse_distortion_full_exp(out, in_, rhoMax, params):
    maxiter = 100
    prec = 1e-16

    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax
    for el in range(len(in_)):
        xi, yi = in_[el, 0:2]
        ri = np.sqrt(xi*xi + yi*yi)
        ri_inv = 1.0/ri
        sinni = yi*ri_inv
        cosni = xi*ri_inv
        ro = ri
        cos2ni = cosni*cosni - sinni*sinni
        sin2ni = 2*sinni*cosni
        cos4ni = cos2ni*cos2ni - sin2ni*sin2ni
        for i in range(maxiter):
            ratio = ri*rxi
            fx = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri - ro
            fxp = (p0*ratio**p3*cos2ni*(p3+1) +
                   p1*ratio**p4*cos4ni*(p4+1) +
                   p2*ratio**p5*(p5+1) + 1)

            delta = fx/fxp
            ri = ri - delta
            if np.abs(delta) <= prec*np.abs(ri):
                break
            
        xi = ri*cosni
        yi = ri*sinni
        out[el, 0] = xi
        out[el, 1] = yi

    return out


@numba.njit('f8[:,:](f8[:,:], f8[:,:], f8, f8[6])')
def direct_distortion_full_exp(out, in_, rhoMax, params):
    p0, p1, p2, p3, p4, p5 = params[0:6]
    rxi = 1.0/rhoMax

    for el in range(len(in_)):
        xi, yi = in_[el, 0:2]
        ri = np.sqrt(xi*xi + yi*yi)
        ri_inv = 1.0/ri
        sinni = yi*ri_inv
        cosni = xi*ri_inv
        cos2ni = cosni*cosni - sinni*sinni
        sin2ni = 2*sinni*cosni
        cos4ni = cos2ni*cos2ni - sin2ni*sin2ni
        ratio = ri*rxi

        ri = (p0*ratio**p3*cos2ni + p1*ratio**p4*cos4ni + p2*ratio**p5 + 1)*ri
        xi = ri*cosni
        yi = ri*sinni
        out[el, 0] = xi
        out[el, 1] = yi

    return out


def GE_41RT_newton_numba_full_exp(xy_in, params, invert=False):
    if params[0] == 0 and params[1] == 0 and params[2] == 0:
        return xy_in
    else:
        rhoMax = 204.8
        if invert:
            xy_out = np.empty_like(xy_in)
            return inverse_distortion_full(xy_out, xy_in, rhoMax, params)
        else:
            xy_out = np.empty_like(xy_in)
            return direct_distortion_full(xy_out, xy_in, rhoMax, params)


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
            with nvtx.profile_range('ge_41rt_bis',
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

            with nvtx.profile_range('ge_41rt_newton_numba_bis',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba = GE_41RT_newton_numba(*args)

            with nvtx.profile_range('ge_41rt_newton_numba_full',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba_full = GE_41RT_newton_numba_full(*args)

            with nvtx.profile_range('ge_41rt_newton_numba_full_bis',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba_full = GE_41RT_newton_numba_full(*args)

            with nvtx.profile_range('ge_41rt_newton_numba_full_exp',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba_full = GE_41RT_newton_numba_full_exp(*args)

            with nvtx.profile_range('ge_41rt_newton_numba_full_exp_bis',
                                    color=_colors[(count*2+1)%len(_colors)],
                                    payload=args[1].shape[0]):
                res_newton_numba_full = GE_41RT_newton_numba_full_exp(*args)

            assert np.allclose(res, res_opt)
            assert np.allclose(res, res_newton)
            assert np.allclose(res, res_newton_numba)
            assert np.allclose(res, res_newton_numba_full)
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
