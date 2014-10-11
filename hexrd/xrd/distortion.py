import numpy as np
from scipy import optimize as opt
import numba

def dummy(xy_in, params, invert=False):
    """
    """
    return xy_in

def newton(x0, f, fp, extra, prec=3e-16, maxiter=100):
    for i in range(maxiter):
        x = x0 - f(x0, *extra) / fp(x0, *extra)
        relerr = np.max(np.abs(x - x0)) / np.max(np.abs(x))
        if relerr < prec:
            # print 'stopping at %d iters' % i
            return x
        x0 = x
    return x0

@numba.njit('void(f8[:], f8[:], f8[:], f8, f8[6])')
def inverse_distortion_numba(rhoOut, rho0, eta0, rhoMax, params):
    # Apply Newton's method to invert the GE_41RT distortion,
    # inlining the function and its inverse to help Numba's JIT
    # produce good machine code.
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
            # Stop when relative error reaches threshold
            if np.abs(delta) <= prec * np.abs(ri):
                break

        rhoOut[el] = ri

def inverse_distortion_numpy(rho0, eta0, rhoMax, params):
    rhoSclFuncInv = lambda ri, ni, ro, rx, p: \
        (p[0]*(ri/rx)**p[3] * np.cos(2.0 * ni) + \
         p[1]*(ri/rx)**p[4] * np.cos(4.0 * ni) + \
         p[2]*(ri/rx)**p[5] + 1)*ri - ro

    rhoSclFIprime = lambda ri, ni, ro, rx, p: \
        p[0]*(ri/rx)**p[3] * np.cos(2.0 * ni) * (p[3] + 1) + \
        p[1]*(ri/rx)**p[4] * np.cos(4.0 * ni) * (p[4] + 1) + \
        p[2]*(ri/rx)**p[5] * (p[5] + 1) + 1

    return newton(rho0, rhoSclFuncInv, rhoSclFIprime,
                  (eta0, rho0, rhoMax, params))

def GE_41RT(xy_in, params, invert=False):
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

        x0 = xy_in[:, 0].flatten()
        y0 = xy_in[:, 1].flatten()

        npts = len(x0)

        # detector relative polar coordinates
        #   - this is the radius that gets rescaled
        rho0 = np.sqrt( x0*x0 + y0*y0 )
        eta0 = np.arctan2( y0, x0 )

        if invert:
            # in here must do nonlinear solve for distortion
            rhoOut = np.empty(npts, dtype=float)
            inverse_distortion(rhoOut, rho0, eta0, rhoMax, params)
            #rhoOut = inverse_distortion_numpy(rhoOut, rho0, eta0, rhoMax, params)
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
