#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Joel Bernier <bernier2@llnl.gov> and others.
# LLNL-CODE-529294.
# All rights reserved.
#
# This file is part of HEXRD. For details on dowloading the source,
# see the file COPYING.
#
# Please also see the file LICENSE.
#
# This program is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program (see file LICENSE); if not, write to
# the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA or visit <http://www.gnu.org/licenses/>.
# ============================================================

import numpy as np
import sys

from hexrd.xrd import _transforms_CAPI

from numpy import float_ as nFloat
from numpy import int_ as nInt

# ######################################################################
# Module Data
epsf      = np.finfo(float).eps      # ~2.2e-16
ten_epsf  = 10 * epsf                # ~2.2e-15
sqrt_epsf = np.sqrt(epsf)            # ~1.5e-8

periodDict   = {'degrees': 360.0, 'radians': 2*np.pi}
angularUnits = 'radians'        # module-level angle units

# basis vectors
I3 = np.eye(3)                                        # (3, 3) identity
Xl = np.ascontiguousarray(I3[:, 0].reshape(3, 1))     # X in the lab frame
Yl = np.ascontiguousarray(I3[:, 1].reshape(3, 1))     # Y in the lab frame
Zl = np.ascontiguousarray(I3[:, 2].reshape(3, 1))     # Z in the lab frame

# reference stretch
vInv_ref = np.array([[1., 1., 1., 0., 0., 0.]], order='C').T

# reference beam direction and eta=0 ref in LAB FRAME for standard geometry
bVec_ref = -Zl
eta_ref  =  Xl

# ######################################################################
# Funtions

def makeGVector(hkl, bMat):
    """
    take a CRYSTAL RELATIVE B matrix onto a list of hkls to output unit
    reciprocal lattice vectors (a.k.a. lattice plane normals)

    Required Arguments:
    hkls -- (3, n) ndarray of n hstacked reciprocal lattice vector component
            triplets
    bMat -- (3, 3) ndarray representing the matirix taking reciprocal lattice
            vectors to the crystal reference frame

    Output:
    gVecs -- (3, n) ndarray of n unit reciprocal lattice vectors
             (a.k.a. lattice plane normals)

    To Do:
    * might benefit from some assert statements to catch improperly shaped
      input.
    """
    assert hkl.shape[0] == 3, 'hkl input must be (3, n)'
    return unitVector(np.dot(bMat, hkl))

def gvecToDetectorXY(gVec_c,
                     rMat_d, rMat_s, rMat_c,
                     tVec_d, tVec_s, tVec_c,
                     beamVec=bVec_ref):
    """
    Takes a list of unit reciprocal lattice vectors in crystal frame to the
    specified detector-relative frame, subject to the conditions:

    1) the reciprocal lattice vector must be able to satisfy a bragg condition
    2) the associated diffracted beam must intersect the detector plane

    Required Arguments:
    gVec_c -- (n, 3) ndarray of n reciprocal lattice vectors in the CRYSTAL FRAME
    rMat_d -- (3, 3) ndarray, the COB taking DETECTOR FRAME components to LAB FRAME
    rMat_s -- (3, 3) ndarray, the COB taking SAMPLE FRAME components to LAB FRAME
    rMat_c -- (3, 3) ndarray, the COB taking CRYSTAL FRAME components to SAMPLE FRAME
    tVec_d -- (3, 1) ndarray, the translation vector connecting LAB to DETECTOR
    tVec_s -- (3, 1) ndarray, the translation vector connecting LAB to SAMPLE
    tVec_c -- (3, 1) ndarray, the translation vector connecting SAMPLE to CRYSTAL

    Outputs:
    (m, 2) ndarray containing the intersections of m <= n diffracted beams
    associated with gVecs
    """
    return _transforms_CAPI.gvecToDetectorXY(np.ascontiguousarray(gVec_c),
                                             rMat_d, rMat_s, rMat_c,
                                             tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(),
                                             beamVec.flatten())
    # ztol = epsf

    # nVec_l = np.dot(rMat_d, Zl)                # detector plane normal
    # bHat_l = unitRowVector(beamVec.reshape(1, 3)) # make sure beam vector is unit
    # P0_l   = tVec_s + np.dot(rMat_s, tVec_c)   # origin of CRYSTAL FRAME
    # P3_l   = tVec_d                            # origin of DETECTOR FRAME

    # # form unit reciprocal lattice vectors in lab frame (w/o translation)
    # rMat_sc = rMat_s.dot(rMat_c)
    # gVec_l = np.dot(unitRowVector(gVec_c), rMat_sc.T)

    # # dot with beam vector (upstream direction)
    # bDot   = np.dot(gVec_l, -bHat_l.T).squeeze()

    # # see who can diffract; initialize output array with NaNs
    # canDiffract = np.atleast_1d( np.logical_and( bDot >= ztol, bDot <= 1. - ztol ) )
    # npts        = sum(canDiffract)
    # retval      = np.nan * np.ones_like(gVec_l)
    # if np.any(canDiffract):  # subset of admissable reciprocal lattice vectors
    #     adm_gVec_l = gVec_l[canDiffract, :].reshape(npts, 3)
    #     dVec_l = np.empty((npts, 3)) # initialize diffracted beam vector array
    #     for ipt in range(npts):
    #         dVec_l[ipt, :] = np.dot(makeBinaryRotMat(adm_gVec_l[ipt, :]), -bHat_l.T).squeeze()
    #         # tmp_op = np.dot(adm_gVec_l[:, ipt].reshape(3, 1),
    #         #                 adm_gVec_l[:, ipt].reshape(1, 3))
    #         # dVec_l[:, ipt] = np.dot(2*tmp_op - I3, -bHat_l).squeeze()
    #         pass
    #     # ###############################################################
    #     # displacement vector calculation

    #     # first check for non-instersections
    #     denom = np.dot(dVec_l, nVec_l).flatten()
    #     dzero = abs(denom) < ztol
    #     denom[dzero] = 1.          # mitigate divide-by-zero
    #     cantIntersect = denom > 0. # index to dVec_l that can't hit det

    #     # displacement scaling (along dVec_l)
    #     u = np.dot(nVec_l.T, P3_l - P0_l).flatten() / denom
        
    #     # filter out non-intersections, fill with NaNs
    #     u[np.logical_or(dzero, cantIntersect)] = np.nan

    #     # diffracted beam points IN DETECTOR FRAME
    #     P2_l = np.empty((npts, 3))
    #     for ipt in range(npts):
    #         P2_l[ipt,:] = P0_l.T + u[ipt] * dVec_l[ipt,:]
    #     P2_d = np.dot(P2_l - tVec_d.T, rMat_d)

    #     # put feasible transformed gVecs into return array
    #     retval[canDiffract, :] = P2_d
    #     pass
    # return retval[:, :2]

def detectorXYToGvec(xy_det,
                     rMat_d, rMat_s,
                     tVec_d, tVec_s, tVec_c,
                     beamVec=bVec_ref, etaVec=eta_ref):
    """
    Takes a list cartesian (x, y) pairs in the detector coordinates and calculates
    the associated reciprocal lattice (G) vectors and (bragg angle, azimuth) pairs
    with respect to the specified beam and azimth (eta) reference directions

    Required Arguments:
    xy_det -- (n, 2) ndarray or list-like input of n detector (x, y) points
    rMat_d -- (3, 3) ndarray, the COB taking DETECTOR FRAME components to LAB FRAME
    rMat_s -- (3, 3) ndarray, the COB taking SAMPLE FRAME components to LAB FRAME
    tVec_d -- (3, 1) ndarray, the translation vector connecting LAB to DETECTOR
    tVec_s -- (3, 1) ndarray, the translation vector connecting LAB to SAMPLE
    tVec_c -- (3, 1) ndarray, the translation vector connecting SAMPLE to CRYSTAL

    Optional Keyword Arguments:
    beamVec -- (3, 1) mdarray containing the incident beam direction components in the LAB FRAME
    etaVec  -- (3, 1) mdarray containing the reference azimuth direction components in the LAB FRAME

    Outputs:
    (n, 2) ndarray containing the (tTh, eta) pairs associated with each (x, y)
    (n, 3) ndarray containing the associated G vector directions in the LAB FRAME
    associated with gVecs
    """
    return _transforms_CAPI.detectorXYToGvec(np.ascontiguousarray(xy_det),
                                             rMat_d, rMat_s,
                                             tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(),
                                             beamVec.flatten(),etaVec.flatten())

def oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength,
                       vInv=None, beamVec=bVec_ref, etaVec=eta_ref):
    """
    Takes a list of unit reciprocal lattice vectors in crystal frame to the
    specified detector-relative frame, subject to the conditions:

    1) the reciprocal lattice vector must be able to satisfy a bragg condition
    2) the associated diffracted beam must intersect the detector plane

    Required Arguments:
    hkls       -- (n, 3) ndarray of n reciprocal lattice vectors in the CRYSTAL FRAME
    chi        -- float representing the inclination angle of the oscillation axis (std coords)
    rMat_c     -- (3, 3) ndarray, the COB taking CRYSTAL FRAME components to SAMPLE FRAME
    bMat       -- (3, 3) ndarray, the COB taking RECIPROCAL LATTICE components to CRYSTAL FRAME
    wavelength -- float representing the x-ray wavelength in Angstroms

    Optional Keyword Arguments:
    beamVec -- (3, 1) mdarray containing the incident beam direction components in the LAB FRAME
    etaVec  -- (3, 1) mdarray containing the reference azimuth direction components in the LAB FRAME

    Outputs:
    ome0 -- (n, 3) ndarray containing the feasible (tTh, eta, ome) triplets for each input hkl (first solution)
    ome1 -- (n, 3) ndarray containing the feasible (tTh, eta, ome) triplets for each input hkl (second solution)

    Notes:
    ------------------------------------------------------------------------
    The reciprocal lattice vector, G, will satisfy the the Bragg condition
    when:

        b.T * G / ||G|| = -sin(theta)

    where b is the incident beam direction (k_i) and theta is the Bragg
    angle consistent with G and the specified wavelength. The components of
    G in the lab frame in this case are obtained using the crystal
    orientation, Rc, and the single-parameter oscillation matrix, Rs(ome):

        Rs(ome) * Rc * G / ||G||

    The equation above can be rearranged to yield an expression of the form:

        a*sin(ome) + b*cos(ome) = c

    which is solved using the relation:

        a*sin(x) + b*cos(x) = sqrt(a**2 + b**2) * sin(x + alpha)

        --> sin(x + alpha) = c / sqrt(a**2 + b**2)

    where:

        alpha = atan2(b, a)

     The solutions are:

                /
                |       arcsin(c / sqrt(a**2 + b**2)) - alpha
            x = <
                |  pi - arcsin(c / sqrt(a**2 + b**2)) - alpha
                \

    There is a double root in the case the reflection is tangent to the
    Debye-Scherrer cone (c**2 = a**2 + b**2), and no solution if the
    Laue condition cannot be satisfied (filled with NaNs in the results
    array here)
    """
    if vInv is None:
        vInv = vInv_ref
    else:
        vInv = np.ascontiguousarray(vInv.flatten())
    return _transforms_CAPI.oscillAnglesOfHKLs(np.ascontiguousarray(hkls),chi,rMat_c,bMat,wavelength,
                                               vInv,beamVec.flatten(),etaVec.flatten())

"""
#######################################################################
######                  Utility Functions                        ######
#######################################################################

"""

def arccosSafe(temp):
    """
    Protect against numbers slightly larger than 1 in magnitude due to round-off
    """

    # Oh, the tricks we must play to make this overloaded and robust...
    if type(temp) is list:
        temp = nd.asarray(temp)
    elif type(temp) is ndarray:
        if len(temp.shape) == 0:
            temp = temp.reshape(1)

    if (temp > 1.00001).any():
        print >> sys.stderr, "attempt to take arccos of %s" % temp
        raise RuntimeError, "unrecoverable error"
    elif (temp < -1.00001).any():
        print >> sys.stderr, "attempt to take arccos of %s" % temp
        raise RuntimeError, "unrecoverable error"

    gte1 = temp >=  1.
    lte1 = temp <= -1.

    temp[gte1] =  1
    temp[lte1] = -1

    ang = arccos(temp)

    return ang

def angularDifference(angList0, angList1, units=angularUnits):
    """
    Do the proper (acute) angular difference in the context of a branch cut.

    *) Default angular range is [-pi, pi]
    """
    period = periodDict[units]
    # take difference as arrays
    diffAngles = np.atleast_1d(angList0) - np.atleast_1d(angList1)

    return abs(np.remainder(diffAngles + 0.5*period, period) - 0.5*period)

def mapAngle(ang, *args, **kwargs):
    """
    Utility routine to map an angle into a specified period
    """
    units  = angularUnits
    period = periodDict[units]

    kwargKeys = kwargs.keys()
    for iArg in range(len(kwargKeys)):
        if kwargKeys[iArg] == 'units':
            units = kwargs[ kwargKeys[iArg] ]
        else:
            raise RuntimeError, "Unknown keyword argument: " + str(kwargKeys[iArg])

    try:
        period = periodDict[units.lower()]
    except:
        raise RuntimeError, "unknown angular units: " + str( kwargs[ kwargKeys[iArg] ] )

    ang = np.atleast_1d(nFloat( ang ) )

    # if we have a specified angular range, use that
    if len(args) > 0:
        angRange = np.atleast_1d(nFloat( args[0] ) )

        # divide of multiples of period
        ang = ang - nInt(ang / period) * period

        lb = angRange.min()
        ub = angRange.max()

        if abs(ub - lb) != period:
            raise RuntimeError, 'range is incomplete!'

        lbi = ang < lb
        while lbi.sum() > 0:
            ang[lbi] = ang[lbi] + period
            lbi = ang < lb
            pass
        ubi = ang > ub
        while ubi.sum() > 0:
            ang[ubi] = ang[ubi] - period
            ubi = ang > ub
            pass
        retval = ang
    else:
        retval = np.mod(ang + 0.5*period, period) - 0.5*period
    return retval

def columnNorm(a):
    """
    normalize array of column vectors (hstacked, axis = 0)
    """
    if len(a.shape) > 2:
        raise RuntimeError, "incorrect shape: arg must be 1-d or 2-d, yours is %d" %(len(a.shape))

    cnrma = np.sqrt(np.sum(np.asarray(a)**2, 0))

    return cnrma

def rowNorm(a):
    """
    normalize array of row vectors (vstacked, axis = 1)
    """
    if len(a.shape) > 2:
        raise RuntimeError, "incorrect shape: arg must be 1-d or 2-d, yours is %d" %(len(a.shape))

    cnrma = np.sqrt(np.sum(np.asarray(a)**2, 1))

    return cnrma

def unitRowVector(vecIn):
    if vecIn.ndim == 1:
        return _transforms_CAPI.unitRowVector(vecIn)
    elif vecIn.ndim == 2:
        return _transforms_CAPI.unitRowVectors(vecIn)
    else:
        assert vecIn.ndim in [1,2], "incorrect arg shape; must be 1-d or 2-d, yours is %d-d" % (a.ndim)

def makeDetectorRotMat(tiltAngles):
    """
    Form the (3, 3) tilt rotations from the tilt angle list:

    tiltAngles = [gamma_Xl, gamma_Yl, gamma_Zl] in radians
    """
    return _transforms_CAPI.makeDetectorRotMat(tiltAngles.flatten())

def makeOscillRotMat(oscillAngles):
    """
    oscillAngles = [chi, ome]
    """
    return _transforms_CAPI.makeOscillRotMat(oscillAngles.flatten())

def makeRotMatOfExpMap(expMap):
    """
    make a rotation matrix from an exponential map
    """
    return _transforms_CAPI.makeRotMatOfExpMap(expMap.flatten())

def makeRotMatOfQuat(quat):
    """
    make rotation matrix from a unit quaternion

    ...check to set if input is unit magnitude?
    """
    return _transforms_CAPI.makeRotMatOfQuat(quat)

def makeBinaryRotMat(axis):
    return _transforms_CAPI.makeBinaryRotMat(axis.flatten())

def makeEtaFrameRotMat(bHat_l, eHat_l):
    return _transforms_CAPI.makeEtaFrameRotMat(bHat_l.flatten(),eHat_l.flatten())

def validateAngleRanges(angList, angMin, angMax, ccw=True):
    return _transforms_CAPI.validateAngleRanges(angList,angMin,angMax,ccw)

def rotate_vecs_about_axis(angle, axis, vecs):
    return _transforms_CAPI.rotate_vecs_about_axis(angle, axis, vecs)

def quat_distance(q1, q2, qsym):
    return _transforms_CAPI.quat_distance(q1, q2, qsym)

#def rotateVecsAboutAxis(angle, axis, vecs):
#    return _transforms_CAPI.rotateVecsAboutAxis(angle, axis, vecs)
