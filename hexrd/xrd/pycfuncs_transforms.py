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
from numbapro import vectorize, jit, autojit, guvectorize, njit, cuda, float64, int64 
import math

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

Z1 = np.array([0., 0., 1.]) 

# reference beam direction and eta=0 ref in LAB FRAME for standard geometry
bVec_ref = -Zl
eta_ref  =  Xl


#numba doesn't handle NaN's correctly
not_a_num = float('nan')

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

#@jit('void(i8, f8[:,:], f8[:,:])')
def nb_unitRowVector_cfunc(n, cIn, cOut):
    #n = len(cIn)
    nrm = 0.0
    for j in range(n):
        nrm += cIn[j] * cIn[j]
    nrm = math.sqrt(nrm)
    if nrm > epsf:
        for j in range(n):
            cOut[j] = cIn[j] / nrm
    else:
        for j in range(n):
            cOut[j] = cIn[j]

#@autojit # 21.9 secs
#@jit('void(f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:] , f8[:], f8[:] , f8[:, :]  )', nopython=True)

#@jit('void(f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:] , f8[:], f8[:] , f8[:, :]  )')
#@guvectorize(['void(f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:,:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:] , f8[:], f8[:] , f8[:, :]  )'], '(m,n),(n,n),(n,n),(n,n),(n,p),(n,p),(n,p),(n,p),(n),(n),(n),(n),(n),(n),(n),(n),(n),(q),(q)->(m,n)')
def gvecToDetectorXY(gVec_c,
                     rMat_d, rMat_s, rMat_c,
                     tVec_d, tVec_s, tVec_c,
                     beamVec,
                     bHat_l, 
                     nVec_l, 
                     P0_l, 
                     P2_l, 
                     P2_d, 
                     P3_l, 
                     gHat_c,
                     gVec_l, 
                     dVec_l,
                     rMat_sc,
                     brMat,  
                     result
                     ):
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
#    return _transforms_CAPI.gvecToDetectorXY(np.ascontiguousarray(gVec_c),
#                                             rMat_d, rMat_s, rMat_c,
#                                             tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(),
#                                             beamVec.flatten())
    ztol = epsf
#    bHat_l = np.zeros(3)
#    nVec_l = np.zeros(3)
#    P0_l = np.zeros(3)
#    P2_l = np.zeros(3)
#    P2_d = np.zeros(3)
#    P3_l = np.zeros(3)
#    gHat_c = np.zeros(3)
#    gVec_l = np.zeros(3)
#    dVec_l = np.zeros(3)
#    rMat_sc = np.zeros(9)
#    brMat = np.zeros(9) 
#    result = np.empty((gVec_c.shape[0], 2)) 

    # Normalize the beam vector 
    nb_unitRowVector_cfunc(3, beamVec, bHat_l)
    
    # Initialize the detector normal and frame origins
    num = 0.0
    for j in range(3): 
        nVec_l[j] = 0.0
        P0_l[j] = tVec_s[j]
        for k in range(3):
            nVec_l[j] += rMat_d[j, k] * Z1[k]
            P0_l[j]   += rMat_s[j, k] * tVec_c[k]
        P3_l[j] = tVec_d[j]
        num += nVec_l[j] * (P3_l[j] - P0_l[j])

    # Compute the matrix product of rMat_s and rMat_c */
    for j in range(3): 
        for k in range(3): 
            rMat_sc[3 * j + k] = 0.0
            for l in range(3): 
	        rMat_sc[3 * j + k] += rMat_s[j, l] * rMat_c[l, k]


    # args = (gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat, dVec_l,
    #           nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result)

    # gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat, dVec_l,
    #           nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result)
    # print(map(typeof, args))
    # cpu_gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat, dVec_l,
    #               nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result)



    gpu_gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat, dVec_l,
                  nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result)


# @guvectorize(["float64[:], float64[:], float64[:], float64[:], float64[:], "
#               "float64, float64[:], float64[:], float64[:], float64, float64[:], "
#               "float64[:], float64[:], float64[:,:], float64[:,:], float64[:]"],
#              "(a,),(b,),(c,),(d,),(e,),"
#              "(), (f,), (g,), (h,), (), (i,),"
#              "(j,), (k,), (l, m), (n, o) -> (c,)")
# def gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat,
#                    dVec_l,
#               nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result):



'''
[array(float64, 2d, C), array(float64, 1d, C), array(float64, 1d, C), array(float64, 1d, C), array(float64, 1d, C), float64, array(float64, 1d, C), array(float64, 1d, C), array(float64, 1d, C), float64, array(float64, 1d, C), array(float64, 1d, C), array(float64, 1d, C), array(float64, 2d, C), array(float64, 2d, C), array(float64, 2d, C)]
gVec_c,
result
'''

def gpu_gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat,
                  dVec_l, nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d,
                  result):

    # gpu_core_loop_kernel.forall(result.shape[0])(gVec_c, rMat_sc, bHat_l, ztol, nVec_l, num, P0_l,
    #                      rMat_d, tVec_d, result)

    # test_loop_kernel(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat,
    #                  dVec_l, nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d,
    #                  result)

    dev_gVec_c = cuda.to_device(gVec_c)
    dev_rMat_sc = cuda.to_device(rMat_sc)
    dev_bHat_l = cuda.to_device(bHat_l)
    dev_nVec_l = cuda.to_device(nVec_l)
    dev_P0_l = cuda.to_device(P0_l)
    dev_rMat_d = cuda.to_device(rMat_d)
    dev_tVec_d = cuda.to_device(tVec_d)
    dev_result = cuda.to_device(result)

    gpu_gvec_core_loop_kernel.forall(result.shape[0])(dev_gVec_c, dev_rMat_sc, dev_bHat_l,
                                                 ztol, dev_nVec_l, num, dev_P0_l,
                                                 dev_rMat_d, dev_tVec_d, dev_result)

    dev_result.copy_to_host(ary=result)


def gvec_test_loop_kernel(gVec_c, rMat_sc, bHat_l, ztol, nVec_l, num, P0_l,
                     rMat_d, tVec_d, result):

    P2_d = np.zeros(3, dtype=np.float64)
    P2_l = np.zeros(3, dtype=np.float64)
    dVec_l = np.zeros(3, dtype=np.float64)
    gVec_l = np.zeros(3, dtype=np.float64)
    gHat_c = np.zeros(3, dtype=np.float64)
    brMat = np.zeros(3 * 3, dtype=np.float64)

    for i in range(gVec_c.shape[0]): # npts
        #nb_unitRowVector_cfunc(3, gVec_c[i], gHat_c)
        nrm = 0.0
        for j in range(3):
            nrm += gVec_c[i, j] * gVec_c[i, j]
        nrm = math.sqrt(nrm)
        if nrm > epsf:
            for j in range(3):
                gHat_c[j] = gVec_c[i, j] / nrm
        else:
            for j in range(3):
                gHat_c[j] = gVec_c[i, j]

        bDot = 0.0
        for j in range(3):
            gVec_l[j] = 0.0
            for k in range(3):
                gVec_l[j] += rMat_sc[3 * j + k] * gHat_c[k]
            bDot -= bHat_l[j] * gVec_l[j]

        if bDot >= ztol and bDot <= 1.0 - ztol:
            #If we are here diffraction is possible so increment the number of admissable vectors */
            #nb_makeBinaryRotMat_cfunc(gVec_l, brMat);
            for j in range(3):
                for k in range(3):
                    brMat[3 *j + k] = 2.0 * gVec_l[j] * gVec_l[k]
            j = 2
            brMat[3 * j + j] -= 1.0

            denom = 0.0

            for j in range(3):
                dVec_l[j] = 0.0
                for k in range(3):
                    dVec_l[j] -= brMat[3 * j + k] * bHat_l[k]
                denom += nVec_l[j] * dVec_l[j]

            if denom < -ztol:
                u = num / denom
                for j in range(3):
                    P2_l[j] = P0_l[j] + u * dVec_l[j]
                for j in range(2):
                    P2_d[j] = 0.0
                    for k in range(3):
                        P2_d[j] += rMat_d[k, j] * (P2_l[k] - tVec_d[k, 0])
                    result[i, j] = P2_d[j]
            else:
                result[i, 0] = not_a_num
                result[i, 1] = not_a_num
        else:
            result[i, 0] = not_a_num
            result[i, 1] = not_a_num


# ("float64[:,:], float64[:], float64[:], float64[:], float64[:], "
#           "float64, float64[:], float64[:], float64[:], float64, float64[:], "
#           "float64[:], float64[:], float64[:,:], float64[:,:], float64[:,:]",
#           debug=True)





@cuda.jit("float64[:,:], float64[:], float64[:], float64, float64[:], "
          "float64, float64[:], float64[:,:], float64[:,:], float64[:,:]")
def gpu_gvec_core_loop_kernel(gVec_c, rMat_sc, bHat_l, ztol, nVec_l, num, P0_l,
                         rMat_d, tVec_d, result):

    P2_d = cuda.local.array(3, dtype=float64)
    P2_l = cuda.local.array(3, dtype=float64)
    dVec_l = cuda.local.array(3, dtype=float64)
    gVec_l = cuda.local.array(3, dtype=float64)
    gHat_c = cuda.local.array(3, dtype=float64)
    brMat = cuda.local.array(3 * 3, dtype=float64)

    i = cuda.grid(1)
    if i >= result.shape[0]:
        return

    nrm = 0.0
    for j in range(3):
        nrm += gVec_c[i, j] * gVec_c[i, j]
    nrm = math.sqrt(nrm)
    if nrm > epsf:
        for j in range(3):
            gHat_c[j] = gVec_c[i, j] / nrm
    else:
        for j in range(3):
            gHat_c[j] = gVec_c[i, j]

    bDot = 0.0
    for j in range(3):
        gVec_l[j] = 0.0
        for k in range(3):
            gVec_l[j] += rMat_sc[3 * j + k] * gHat_c[k]
        bDot -= bHat_l[j] * gVec_l[j]

    if bDot >= ztol and bDot <= 1.0 - ztol:
        #If we are here diffraction is possible so increment the number of admissable vectors */
        #nb_makeBinaryRotMat_cfunc(gVec_l, brMat);
        for j in range(3):
            for k in range(3):
                brMat[3 *j + k] = 2.0 * gVec_l[j] * gVec_l[k]
        j = 2
        brMat[3 * j + j] -= 1.0

        denom = 0.0

        for j in range(3):
            dVec_l[j] = 0.0
            for k in range(3):
                dVec_l[j] -= brMat[3 * j + k] * bHat_l[k]
            denom += nVec_l[j] * dVec_l[j]

        if denom < -ztol:
            u = num / denom
            for j in range(3):
                P2_l[j] = P0_l[j] + u * dVec_l[j]
            for j in range(2):
                P2_d[j] = 0.0
                for k in range(3):
                    P2_d[j] += rMat_d[k, j] * (P2_l[k] - tVec_d[k, 0])
                result[i, j] = P2_d[j]
        else:
            result[i, 0] = not_a_num
            result[i, 1] = not_a_num
    else:
        result[i, 0] = not_a_num
        result[i, 1] = not_a_num                


@njit
def cpu_gvec_core_loop(gVec_c, gHat_c, gVec_l, rMat_sc, bHat_l, ztol, brMat, dVec_l,
              nVec_l, num, P2_l, P0_l, P2_d, rMat_d, tVec_d, result):

    for i in range(gVec_c.shape[0]): # npts
        #nb_unitRowVector_cfunc(3, gVec_c[i], gHat_c)
        nrm = 0.0
        for j in range(3):
            nrm += gVec_c[i, j] * gVec_c[i, j]
        nrm = math.sqrt(nrm)
        if nrm > epsf:
            for j in range(3):
                gHat_c[j] = gVec_c[i, j] / nrm
        else:
            for j in range(3):
                gHat_c[j] = gVec_c[i, j]

        bDot = 0.0
        for j in range(3):
            gVec_l[j] = 0.0
            for k in range(3):
                gVec_l[j] += rMat_sc[3 * j + k] * gHat_c[k]
            bDot -= bHat_l[j] * gVec_l[j]

        if bDot >= ztol and bDot <= 1.0 - ztol:
            #If we are here diffraction is possible so increment the number of admissable vectors */
            #nb_makeBinaryRotMat_cfunc(gVec_l, brMat);
            for j in range(3):
                for k in range(3):
                    brMat[3 *j + k] = 2.0 * gVec_l[j] * gVec_l[k]
            j = 2
            brMat[3 * j + j] -= 1.0

            denom = 0.0

            for j in range(3):
                dVec_l[j] = 0.0
                for k in range(3):
                    dVec_l[j] -= brMat[3 * j + k] * bHat_l[k]
                denom += nVec_l[j] * dVec_l[j]

            if denom < -ztol:
                u = num / denom;
                for j in range(3):
                    P2_l[j] = P0_l[j] + u * dVec_l[j]
                for j in range(2):
                    P2_d[j] = 0.0
                    for k in range(3):
                        P2_d[j] += rMat_d[k, j] * (P2_l[k] - tVec_d[k, 0])
                    result[i, j] = P2_d[j]
            else:
                result[i, 0] = not_a_num
                result[i, 1] = not_a_num
        else:
            result[i, 0] = not_a_num
            result[i, 1] = not_a_num

def nb_makeEtaFrameRotMat_cfunc(bPtr, ePtr, rPtr):
    # matrices dim
    # bPtr (3,)
    # ePtr (3,)
    # rPtr 9,
   
    nrmZ = 0.0;
    for i in range(3):
        nrmZ += bPtr[i] * bPtr[i]
    nrmZ = math.sqrt(nrmZ);

    # Assign Z column */

    for i in range(3):
        rPtr[3 * i + 2] = -bPtr[i] / nrmZ
 
    # Determine dot product of Z column and eHat */
    dotZe = 0.0;
    for i in range(3):
        dotZe += rPtr[3 * i + 2] * ePtr[i]

    #Assign X column */
    for i in range(3):
        rPtr[3 * i + 0] = ePtr[i] - dotZe * rPtr[3 * i + 2]

    #Normalize X column */
    nrmX = 0.0;
    for i in range(3):
        nrmX += rPtr[3 * i + 0] * rPtr[3 * i + 0]
       
    nrmX = math.sqrt(nrmX)

    #Assign Y column
    for i in range(3):
        rPtr[3*i+1] = rPtr[3 * ((i+1) % 3) + 2] * rPtr[3 * ((i+2) %3) + 0] - rPtr[3 * ((i+2) % 3) + 2] * rPtr[3 * ((i+1) % 3) + 0]
  
def nb_rotate_vecs_about_axis_cfunc(na, angles, nax, axes, nv, vecs, rVecs, row):

#  int i, j, sa, sax;
#  double c, s, nrm, proj, aCrossV[3];

    aCrossV = np.zeros(3)

    if na == 1:
        sa = 0
    else: 
        sa = 1
    if nax == 1:
        sax = 0
    else: 
        sax = 3

    for i in range(nv): 
        # Rotate using the Rodrigues' Rotation Formula 
        # nb the C code was written with angles being an array but
        # in reality it is just a scalar and i = 0 and nv = 1 
        # 
        c = math.cos(angles)
        s = math.sin(angles)

        # Compute projection of vec along axis 
        proj = 0.0
        for j in range(3): 
            proj += axes[sax * i + j] * vecs[3 * i + j]

        # Compute norm of axis 
        if nax > 1 or i == 0: 
            nrm = 0.0
            for j in range(3): 
	        nrm += axes[sax * i + j] * axes[sax * i + j]
            nrm = math.sqrt(nrm);

        # Compute projection of vec along axis */
        proj = 0.0
        for j in range(3):
            proj += axes[sax * i + j] * vecs[3 * i + j]

        # Compute the cross product of the axis with vec */
        for j in range(3): 
            aCrossV[j] = axes[sax * i + (j + 1) % 3] * vecs[3 * i + (j + 2) % 3] - axes[sax * i + (j + 2) % 3] * vecs[3 * i+ (j + 1) % 3]

        # Combine the three terms to compute the rotated vector */
        for j in range(3):
            rVecs[row, j] = c * vecs[3 * i + j] + (s / nrm) * aCrossV[j] + (1.0 - c) * proj * axes[sax * i + j] / (nrm * nrm)


def detectorXYToGvec(xy_det,
                     rMat_d, rMat_s,
                     tVec_d, tVec_s, tVec_c,
                     beamVec, etaVec,
                     rMat_e,
                     bVec,
                     tVec1,
                     tVec2,
                     dHat_l, 
                     n_g,
                     npts,
                     tTh,
                     eta,
                     gVec_l):
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
    
#    rMat_e = np.zeros(9)
#    bVec = np.zeros(3)
#    tVec1 = np.zeros(3)
#    tVec2 = np.zeros(3)
#    dHat_l = np.zeros(3)
#    n_g = np.zeros(3)
#    npts = xy_det.shape[0]
#    #return values
#    tTh = np.zeros(npts)
#    eta = np.zeros(npts)
#    gVec_l = np.zeros((npts, 3))

    # Fill rMat_e */
    nb_makeEtaFrameRotMat_cfunc(beamVec, etaVec, rMat_e);

    # Normalize the beam vector
    nrm = 0.0
    for j in range(3):
        nrm += beamVec[j] * beamVec[j]
    nrm = math.sqrt(nrm)
    if nrm > epsf: 
        for j in range(3):
            bVec[j] = beamVec[j] / nrm
    else:
        for j in range(3):
            bVec[j] = beamVec[j]

    # Compute shift vector 
    for j in range(3):
        tVec1[j] = tVec_d[j] - tVec_s[j]
        for k in range(3):
            tVec1[j] -= rMat_s[j, k] * tVec_c[k]

    gpu_detector_core_loop(xy_det, rMat_d, rMat_e, bVec, tVec1, tTh, eta, gVec_l)

    
    # cpu  
    #detector_core_loop(xy_det, rMat_d, rMat_e, bVec, tVec1, npts, tTh, eta, gVec_l)

def gpu_detector_core_loop(xy_det, rMat_d,  
        rMat_e, bVec, tVec1, tTh, eta, gVec_l):

    dev_xy_det = cuda.to_device(xy_det)
    dev_rMat_d = cuda.to_device(rMat_d)
    dev_rMat_e = cuda.to_device(rMat_e)
    dev_bVec = cuda.to_device(bVec)
    dev_tVec1 = cuda.to_device(tVec1)
    dev_tTh = cuda.to_device(tTh)
    dev_eta = cuda.to_device(eta)
    dev_gVec_l = cuda.to_device(gVec_l)
      
    gpu_detector_core_loop_kernel.forall(xy_det.shape[0])(dev_xy_det, dev_rMat_d, dev_rMat_e, dev_bVec, dev_tVec1, dev_tTh, dev_eta, dev_gVec_l)

    dev_tTh.copy_to_host(ary=tTh)
    dev_eta.copy_to_host(ary=eta)
    dev_gVec_l.copy_to_host(ary=gVec_l)

 
@cuda.jit("float64[:,:], float64[:,:], float64[:], float64[:], float64[:], "
        "float64[:], float64[:], float64[:,:]")
def gpu_detector_core_loop_kernel(xy_det, rMat_d,
                       rMat_e, bVec, tVec1, tTh, eta, gVec_l):


    dHat_l = cuda.local.array(3, dtype=float64)
    tVec2 = cuda.local.array(3, dtype=float64)
    n_g = cuda.local.array(3, dtype=float64)
    aCrossV = cuda.local.array(3, dtype=float64)

    i = cuda.grid(1)
    if i >= xy_det.shape[0]:
        return


    # Compute dHat_l vector
    nrm = 0.0;
    for j in range(3): 
        dHat_l[j] = tVec1[j]
        for k in range(2): 
            dHat_l[j] += rMat_d[j, k] * xy_det[i, k]
        nrm += dHat_l[j] * dHat_l[j]
    if nrm > epsf:
        for j in range(3):
            dHat_l[j] /= math.sqrt(nrm)

    # Compute tTh 
    nrm = 0.0;
    for j in range(3):
        nrm += bVec[j] * dHat_l[j]

    tTh[i] = math.acos(nrm)

    # Compute eta 
    for j in range(2):
        tVec2[j] = 0.0
        for k in range(3):
            tVec2[j] += rMat_e[3 * k + j] * dHat_l[k]

    eta[i] = math.atan2(tVec2[1], tVec2[0])

   # Compute n_g vector
    nrm = 0.0
    for j in range(3):
        n_g[j] = bVec[(j + 1) % 3] * dHat_l[(j + 2) % 3] - bVec[(j + 2) % 3] * dHat_l[(j + 1) % 3]
        nrm += n_g[j] * n_g[j]
    nrm = math.sqrt(nrm)
    for j in range(3):
        n_g[j] /= nrm

    # Rotate dHat_l vector
    phi = 0.5 * (math.pi - tTh[i])

    #nb_rotate_vecs_about_axis_cfunc(1, phi, 1, n_g, 1, dHat_l, gVec_l, i)
    # rewrote this as local function 
    # removed for loop which went from for j in range(1)
    # and removed beginning if since na == 1 and nax == 1 always in this case
    c = math.cos(phi)
    s = math.sin(phi)
    
    # Compute projection of vec along axis 
    proj = 0.0
    for j in range(3): 
        proj += n_g[j] * dHat_l[j]
  

    # Compute norm of axis 
    nrm = 0.0
    for j in range(3): 
        nrm += n_g[j] * n_g[j]
    nrm = math.sqrt(nrm);


    # Compute projection of vec along axis 
    # not sure why this is done again but this is a copy of the
    # original code
    proj = 0.0
    for j in range(3):
        proj += n_g[j] * dHat_l[j]
   
    # Compute the cross product of the axis with vec */
    for j in range(3): 
        aCrossV[j] = n_g[(j + 1) % 3] * dHat_l[(j + 2) % 3] - n_g[(j + 2) % 3] * dHat_l[(j + 1) % 3]

    # Combine the three terms to compute the rotated vector */
    for j in range(3):
        gVec_l[i, j] = c * dHat_l[j] + (s / nrm) * aCrossV[j] + (1.0 - c) * proj * n_g[j] / (nrm * nrm)


@jit
def detector_core_loop(xy_det, rMat_d,
                       rMat_e, bVec, tVec1, npts, tTh, eta, gVec_l):

    dHat_l = np.zeros(3)
    tVec2 = np.zeros(3)
    n_g = np.zeros(3)
    aCrossV = np.zeros(3)

    for i in range(npts): 
        # Compute dHat_l vector
        nrm = 0.0;
        for j in range(3): 
            dHat_l[j] = tVec1[j]
            for k in range(2): 
	        dHat_l[j] += rMat_d[j, k] * xy_det[i, k]
            nrm += dHat_l[j] * dHat_l[j]
        if nrm > epsf:
            for j in range(3):
	        dHat_l[j] /= math.sqrt(nrm)

        # Compute tTh 
        nrm = 0.0;
        for j in range(3):
            nrm += bVec[j] * dHat_l[j]
    
        tTh[i] = math.acos(nrm)

        # Compute eta 
        for j in range(2):
            tVec2[j] = 0.0
            for k in range(3):
	        tVec2[j] += rMat_e[3 * k + j] * dHat_l[k]
   
        eta[i] = math.atan2(tVec2[1], tVec2[0])

       # Compute n_g vector
        nrm = 0.0
        for j in range(3):
            n_g[j] = bVec[(j + 1) % 3] * dHat_l[(j + 2) % 3] - bVec[(j + 2) % 3] * dHat_l[(j + 1) % 3]
            nrm += n_g[j] * n_g[j]
        nrm = math.sqrt(nrm)
        for j in range(3):
            n_g[j] /= nrm

        # Rotate dHat_l vector
        phi = 0.5 * (math.pi - tTh[i])


        #nb_rotate_vecs_about_axis_cfunc(1, phi, 1, n_g, 1, dHat_l, gVec_l, i)
        # rewrote this as local function 
        # removed for loop which went from for j in range(1)
        # and removed beginning if since na == 1 and nax == 1 always in this case
        sa = 0
        sax = 0
        c = math.cos(phi)
        s = math.sin(phi)
        
        # Compute projection of vec along axis 
        proj = 0.0
        for j in range(3): 
            proj += n_g[j] * dHat_l[j]
      

        # Compute norm of axis 
        nrm = 0.0
        for j in range(3): 
	    nrm += n_g[j] * n_g[j]
        nrm = math.sqrt(nrm);

 
        # Compute projection of vec along axis 
        # not sure why this is done again but this is a copy of the
        # original code
        proj = 0.0
        for j in range(3):
            proj += n_g[j] * dHat_l[j]
       
        # Compute the cross product of the axis with vec */
        for j in range(3): 
            aCrossV[j] = n_g[(j + 1) % 3] * dHat_l[(j + 2) % 3] - n_g[(j + 2) % 3] * dHat_l[(j + 1) % 3]

        # Combine the three terms to compute the rotated vector */
        for j in range(3):
            gVec_l[i, j] = c * dHat_l[j] + (s / nrm) * aCrossV[j] + (1.0 - c) * proj * n_g[j] / (nrm * nrm)



        
#    return ((tTh, eta), gVec_l)

#    return _transforms_CAPI.detectorXYToGvec(np.ascontiguousarray(xy_det),
#                                             rMat_d, rMat_s,
#                                             tVec_d.flatten(), tVec_s.flatten(), tVec_c.flatten(),
#                                             beamVec.flatten(),etaVec.flatten())


def nb_makeOscillRotMat_cfunc(oPtr, rPtr):
#  int i;
#  double c[2],s[2];
#
#  for (i=0; i<2; i++) {
#    c[i] = cos(oPtr[i]);
#    s[i] = sin(oPtr[i]);
#  }

  rPtr[0] =  math.cos(oPtr[1])  # c[1];
  rPtr[1] =  0.0;
  rPtr[2] =  math.sin(oPtr[1])  # s[1];
  rPtr[3] =  math.sin(oPtr[0]) * math.sin(oPtr[1]) # s[0]*s[1];
  rPtr[4] =  math.cos(oPtr[0])  # c[0];
  rPtr[5] =  -math.sin(oPtr[0]) * math.cos(oPtr[1])  # -s[0]*c[1];
  rPtr[6] = -math.cos(oPtr[0]) * math.sin(oPtr[1])  # -c[0]*s[1];
  rPtr[7] = math.sin(oPtr[0])  # s[0];
  rPtr[8] =  math.cos(oPtr[0]) * math.cos(oPtr[1]) # c[0]*c[1];



def oscillAnglesOfHKLs(hkls, chi, rMat_c, bMat, wavelength,
                       beamVec, etaVec, 
                       gVec_e, gHat_c, gHat_s, 
                       bHat_l, eHat_l, oVec, tVec0, 
                       rMat_e, rMat_s, npts,
                       oangs0, oangs1):
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
#    gVec_e = np.zeros(3)
#    gHat_c = np.zeros(3)
#    gHat_s = np.zeros(3)
#    bHat_l = np.zeros(3)
#    eHat_l = np.zeros(3) 
#    oVec = np.zeros(2)
#    tVec0 = np.zeros(3)
#    rMat_e = np.zeros(9)
#    rMat_s = np.zeros(9)
#    npts = hkls.shape[0]
#    #return arrays
#    oangs0 = np.zeros((npts, 3))
#    oangs1 = np.zeros((npts, 3))

    crc = 0
    # Normalize the beam vector
    nrm0 = 0.0
    for j in range(3):
        nrm0 += beamVec[j] * beamVec[j]
    nrm0 = math.sqrt(nrm0)
    if nrm0 > epsf: 
        for j in range(3): 
            bHat_l[j] = beamVec[j] / nrm0
    else:
        for j in range(3):
            bHat_l[j] = beamVec[j]

    # Normalize the eta vector 
    nrm0 = 0.0
    for j in range(3): 
        nrm0 += etaVec[j] * etaVec[j]
    nrm0 = math.sqrt(nrm0)
    if nrm0 > epsf:
        for j in range(3):
            eHat_l[j] = etaVec[j] / nrm0
    else:
        for j in range(3): 
            eHat_l[j] = etaVec[j]

    # Check for consistent reference coordiantes
    nrm0 = 0.0
    for j in range(3): 
        nrm0 += bHat_l[j] * eHat_l[j]
  
    if math.fabs(nrm0) < 1.0 - sqrt_epsf:
        crc = 1

    # Compute the sine and cosine of the oscillation axis tilt
    cchi = math.cos(chi);
    schi = math.sin(chi);

    # move main loop to its own function so it can be numbaized
    #oscill_core_loop(hkls, chi, rMat_c, bMat, wavelength,
    #                   beamVec, etaVec, 
    #                   crc, cchi, schi,
    #                   bHat_l, eHat_l, 
    #                   oangs0, oangs1)

    gpu_oscill_core_loop(hkls, chi, rMat_c, bMat, wavelength,
                       beamVec, etaVec, 
                       crc, cchi, schi,
                       bHat_l, eHat_l, 
                       oangs0, oangs1)



def gpu_oscill_core_loop(hkls, chi, rMat_c, bMat, wavelength,
                       beamVec, etaVec, 
                       crc, cchi, schi,
                       bHat_l, eHat_l, 
                       oangs0, oangs1):

    dev_hkls = cuda.to_device(hkls)
    dev_rMat_c = cuda.to_device(rMat_c)
    dev_bMat = cuda.to_device(bMat)
    dev_beamVec = cuda.to_device(beamVec)
    dev_etaVec = cuda.to_device(etaVec)
    dev_bHat_l = cuda.to_device(bHat_l)
    dev_eHat_l = cuda.to_device(eHat_l)
    dev_oangs0 = cuda.to_device(oangs0)
    dev_oangs1 = cuda.to_device(oangs1)

    gpu_oscill_core_loop_kernel.forall(hkls.shape[0])(dev_hkls, chi, dev_rMat_c, dev_bMat, wavelength,
                       dev_beamVec, dev_etaVec, 
                       crc, cchi, schi,
                       dev_bHat_l, dev_eHat_l,  
                       dev_oangs0, dev_oangs1)

    dev_oangs0.copy_to_host(ary=oangs0)
    dev_oangs1.copy_to_host(ary=oangs1)


@cuda.jit("float64[:,:], float64, float64[:,:], float64[:,:], float64, float64[:], float64[:], int64, float64, float64, float64[:], float64[:], float64[:,:], float64[:,:]")
def gpu_oscill_core_loop_kernel(hkls, chi, rMat_c, bMat, wavelength,
                       beamVec, etaVec, 
                       crc, cchi, schi,
                       bHat_l, eHat_l,  
                       oangs0, oangs1):

    gHat_c = cuda.local.array(3, dtype=float64)
    gHat_s = cuda.local.array(3, dtype=float64)
    oVec = cuda.local.array(2, dtype=float64)
    tVec0 = cuda.local.array(3, dtype=float64)
    gVec_e = cuda.local.array(3, dtype=float64)
    #c = cuda.local.array(2, dtype=float64)
    #s = cuda.local.array(2, dtype=float64)
    rMat_s = cuda.local.array(9, dtype=float64)
    rMat_e = cuda.local.array(9, dtype=float64)

    i = cuda.grid(1)
    if i >= hkls.shape[0]:
        return

    # Compute gVec_c 
    nrm0 = 0.0
    for j in range(3):
        gHat_c[j] = 0.0
        for k in range(3):
            gHat_c[j] += bMat[j, k] * hkls[i, k]
        nrm0 += gHat_c[j] * gHat_c[j]
    nrm0 = math.sqrt(nrm0)

    # Compute gVec_s
    nrm1 = 0.0
    for j in range(3):
        gHat_s[j] = 0.0
        for k in range(3): 
            gHat_s[j] += rMat_c[j, k] * gHat_c[k]
        nrm1 += gHat_s[j] * gHat_s[j]
    nrm1 = math.sqrt(nrm1)

    # Normalize gVec_c to make gHat_c 
    if nrm0 > epsf:
        for j in range(3): 
            gHat_c[j] /= nrm0

    # Normalize gVec_s to make gHat_s */
    if nrm1 > epsf:
        for j in range(3): 
            gHat_s[j] /= nrm1

    # Compute the sine of the Bragg angle */
    sintht = 0.5 * wavelength * nrm0

    # Compute the coefficients of the harmonic equation
    a = gHat_s[2] * bHat_l[0] + schi * gHat_s[0] * bHat_l[1] - cchi * gHat_s[0] * bHat_l[2]
    b = gHat_s[0] * bHat_l[0] - schi * gHat_s[2] * bHat_l[1] + cchi * gHat_s[2] * bHat_l[2]
    c = -sintht - cchi * gHat_s[1] * bHat_l[1] - schi * gHat_s[1] * bHat_l[2]

    # Form solution
    abMag = math.sqrt(a * a + b * b)
    if abMag > 0.0:
        #print('abMag should be <= 0')
        #print('bailing')
        exit
    #assert abMag > 0.0, "abMag <= 0.0" 
    phaseAng = math.atan2(b, a)
    rhs = c / abMag

    if math.fabs(rhs) > 1.0: 
        for j in range(3): 
            oangs0[i, j] = not_a_num
            oangs1[i, j] = not_a_num
        #continue
        return 
    
    #try:
    if rhs >= -1.0 and rhs <= 1:  # check domain
        rhsAng = math.asin(rhs)
    #except ValueError:
    else:
        rhsAng = not_a_num

   # Write ome angles
    oangs0[i, 2] = rhsAng - phaseAng
    oangs1[i, 2] = math.pi - rhsAng - phaseAng

    if crc == 1:
        #nb_makeEtaFrameRotMat_cfunc(bHat_l, eHat_l, rMat_e)

        nrmZ = 0.0
        for j in range(3): 
            nrmZ += bHat_l[j] * bHat_l[j]
        nrmZ = math.sqrt(nrmZ)

        # Assign Z column */
        for j in range(3):
            rMat_e[3 * j + 2] = -bHat_l[j] / nrmZ

        # Determine dot product of Z column and eHat
        dotZe = 0.0
        for j in range(3):
            dotZe += rMat_e[3 * j + 2] * eHat_l[j]

        # Assign X column
        for j in range(3): 
           rMat_e[3 * j + 0] = eHat_l[j] - dotZe * rMat_e[3 * j + 2]

        # Normalize X column 
        nrmX = 0.0
        for j in range(3): 
            nrmX += rMat_e[3 * j] * rMat_e[3 * j]
        nrmX = math.sqrt(nrmX)

        # Assign Y column
        for j in range(3):
            rMat_e[3 * j + 1] = rMat_e[3 * ((j + 1) % 3) + 2] * rMat_e[3 * ((j + 2) % 3) + 0] - rMat_e[3 * ((j +2) % 3) + 2] * rMat_e[3 *((j + 1) % 3) + 0];


        oVec[0] = chi
        oVec[1] = oangs0[i, 2]

    #    nb_makeOscillRotMat_cfunc(oVec, rMat_s)
        #inlined this function
        rMat_s[0] =  math.cos(oVec[1])  # c[1];
        rMat_s[1] =  0.0
        rMat_s[2] =  math.sin(oVec[1])  # s[1];
        rMat_s[3] =  math.sin(oVec[0]) * math.sin(oVec[1]) # s[0]*s[1];
        rMat_s[4] =  math.cos(oVec[0])  # c[0];
        rMat_s[5] =  -math.sin(oVec[0]) * math.cos(oVec[1])  # -s[0]*c[1];
        rMat_s[6] = -math.cos(oVec[0]) * math.sin(oVec[1])  # -c[0]*s[1];
        rMat_s[7] = math.sin(oVec[0])  # s[0];
        rMat_s[8] =  math.cos(oVec[0]) * math.cos(oVec[1]) # c[0]*c[1];


        for j in range(3):
            tVec0[j] = 0.0
            for k in range(3): 
                tVec0[j] += rMat_s[3 * j + k] * gHat_s[k]
        for j in range(2):
            gVec_e[j] = 0.0
            for k in range(3): 
                gVec_e[j] += rMat_e[3 * k + j] * tVec0[k]
      
        oangs0[i, 1] = math.atan2(gVec_e[1], gVec_e[0])

        oVec[1] = oangs1[i, 2]
        #nb_makeOscillRotMat_cfunc(oVec, rMat_s)
        rMat_s[0] =  math.cos(oVec[1])  # c[1];
        rMat_s[1] =  0.0
        rMat_s[2] =  math.sin(oVec[1])  # s[1];
        rMat_s[3] =  math.sin(oVec[0]) * math.sin(oVec[1]) # s[0]*s[1];
        rMat_s[4] =  math.cos(oVec[0])  # c[0];
        rMat_s[5] =  -math.sin(oVec[0]) * math.cos(oVec[1])  # -s[0]*c[1];
        rMat_s[6] = -math.cos(oVec[0]) * math.sin(oVec[1])  # -c[0]*s[1];
        rMat_s[7] = math.sin(oVec[0])  # s[0];
        rMat_s[8] =  math.cos(oVec[0]) * math.cos(oVec[1]) # c[0]*c[1];

        for j in range(3): 
            tVec0[j] = 0.0
            for k in range(3): 
                tVec0[j] += rMat_s[3 * j + k] * gHat_s[k]
        for j in range(2):
            gVec_e[j] = 0.0
            for k in range(3):
                gVec_e[j] += rMat_e[3 * k + j] * tVec0[k]
        oangs1[i, 1] = math.atan2(gVec_e[1], gVec_e[0]);

        oangs0[i, 0] = 2.0 * math.asin(sintht)
        oangs1[i, 0] = oangs0[i, 0]


def oscill_core_loop(hkls, chi, rMat_c, bMat, wavelength,
                   beamVec, etaVec, 
                   crc, cchi, schi,
                   bHat_l, eHat_l, 
                   oangs0, oangs1):

    gHat_c = np.zeros(3)
    gHat_s = np.zeros(3) 
    oVec = np.zeros(2) 
    tVec0 = np.zeros(3)
    gVec_e = np.zeros(3) 
    rMat_s = np.zeros(9)
    rMat_e = np.zeros(9)

    for i in range(hkls.shape[0]): # npts

        # Compute gVec_c 
        nrm0 = 0.0
        for j in range(3):
            gHat_c[j] = 0.0
            for k in range(3):
                gHat_c[j] += bMat[j, k] * hkls[i, k]
            nrm0 += gHat_c[j] * gHat_c[j]
        nrm0 = math.sqrt(nrm0)

        # Compute gVec_s
        nrm1 = 0.0
        for j in range(3):
            gHat_s[j] = 0.0
            for k in range(3): 
                gHat_s[j] += rMat_c[j, k] * gHat_c[k]
            nrm1 += gHat_s[j] * gHat_s[j]
        nrm1 = math.sqrt(nrm1)

        # Normalize gVec_c to make gHat_c 
        if nrm0 > epsf:
            for j in range(3): 
                gHat_c[j] /= nrm0

        # Normalize gVec_s to make gHat_s */
        if nrm1 > epsf:
            for j in range(3): 
                gHat_s[j] /= nrm1

        # Compute the sine of the Bragg angle */
        sintht = 0.5 * wavelength * nrm0

        # Compute the coefficients of the harmonic equation
        a = gHat_s[2] * bHat_l[0] + schi * gHat_s[0] * bHat_l[1] - cchi * gHat_s[0] * bHat_l[2]
        b = gHat_s[0] * bHat_l[0] - schi * gHat_s[2] * bHat_l[1] + cchi * gHat_s[2] * bHat_l[2]
        c = -sintht - cchi * gHat_s[1] * bHat_l[1] - schi * gHat_s[1] * bHat_l[2]

        # Form solution
        abMag = math.sqrt(a * a + b * b) 
        assert abMag > 0.0, "abMag <= 0.0" 
        phaseAng = math.atan2(b, a)
        rhs = c / abMag

        if math.fabs(rhs) > 1.0: 
            for j in range(3): 
                oangs0[i, j] = not_a_num
                oangs1[i, j] = not_a_num
            continue

        if rhs >= -1.0 and rhs <= 1:  # check domain
            rhsAng = math.asin(rhs)
        else:
            rhsAng = not_a_num

       # try:
       #     rhsAng = math.asin(rhs)
       # except ValueError:
       #     rhsAng = not_a_num

       # Write ome angles
        oangs0[i, 2] = rhsAng - phaseAng
        oangs1[i, 2] = math.pi - rhsAng - phaseAng

        if crc:
            #nb_makeEtaFrameRotMat_cfunc(bHat_l, eHat_l, rMat_e)
            nrmZ = 0.0
            for j in range(3): 
                nrmZ += bHat_l[j] * bHat_l[j]
            nrmZ = math.sqrt(nrmZ)

            # Assign Z column */
            for j in range(3):
                rMat_e[3 * j + 2] = -bHat_l[j] / nrmZ

            # Determine dot product of Z column and eHat
            dotZe = 0.0
            for j in range(3):
                dotZe += rMat_e[3 * j + 2] * eHat_l[j]

            # Assign X column
            for j in range(3): 
               rMat_e[3 * j + 0] = eHat_l[j] - dotZe * rMat_e[3 * j + 2]

            # Normalize X column 
            nrmX = 0.0
            for j in range(3): 
                nrmX += rMat_e[3 * j] * rMat_e[3 * j]
            nrmX = math.sqrt(nrmX)

            # Assign Y column
            for j in range(3):
                rMat_e[3 * j + 1] = rMat_e[3 * ((j + 1) % 3) + 2] * rMat_e[3 * ((j + 2) % 3) + 0] - rMat_e[3 * ((j +2) % 3) + 2] * rMat_e[3 *((j + 1) % 3) + 0];

            oVec[0] = chi
            oVec[1] = oangs0[i, 2]
            #nb_makeOscillRotMat_cfunc(oVec, rMat_s)

            rMat_s[0] =  math.cos(oVec[1])  # c[1];
            rMat_s[1] =  0.0
            rMat_s[2] =  math.sin(oVec[1])  # s[1];
            rMat_s[3] =  math.sin(oVec[0]) * math.sin(oVec[1]) # s[0]*s[1];
            rMat_s[4] =  math.cos(oVec[0])  # c[0];
            rMat_s[5] =  -math.sin(oVec[0]) * math.cos(oVec[1])  # -s[0]*c[1];
            rMat_s[6] = -math.cos(oVec[0]) * math.sin(oVec[1])  # -c[0]*s[1];
            rMat_s[7] = math.sin(oVec[0])  # s[0];
            rMat_s[8] =  math.cos(oVec[0]) * math.cos(oVec[1]) # c[0]*c[1];




            for j in range(3):
                tVec0[j] = 0.0
                for k in range(3): 
                    tVec0[j] += rMat_s[3 * j + k] * gHat_s[k]
            for j in range(2):
                gVec_e[j] = 0.0
                for k in range(3): 
                    gVec_e[j] += rMat_e[3 * k + j] * tVec0[k]
          
            oangs0[i, 1] = math.atan2(gVec_e[1], gVec_e[0])

            oVec[1] = oangs1[i, 2]
            #nb_makeOscillRotMat_cfunc(oVec, rMat_s)

            rMat_s[0] =  math.cos(oVec[1])  # c[1];
            rMat_s[1] =  0.0
            rMat_s[2] =  math.sin(oVec[1])  # s[1];
            rMat_s[3] =  math.sin(oVec[0]) * math.sin(oVec[1]) # s[0]*s[1];
            rMat_s[4] =  math.cos(oVec[0])  # c[0];
            rMat_s[5] =  -math.sin(oVec[0]) * math.cos(oVec[1])  # -s[0]*c[1];
            rMat_s[6] = -math.cos(oVec[0]) * math.sin(oVec[1])  # -c[0]*s[1];
            rMat_s[7] = math.sin(oVec[0])  # s[0];
            rMat_s[8] =  math.cos(oVec[0]) * math.cos(oVec[1]) # c[0]*c[1];



            for j in range(3): 
                tVec0[j] = 0.0
                for k in range(3): 
                    tVec0[j] += rMat_s[3 * j + k] * gHat_s[k]
            for j in range(2):
                gVec_e[j] = 0.0
                for k in range(3):
                    gVec_e[j] += rMat_e[3 * k + j] * tVec0[k]
            oangs1[i, 1] = math.atan2(gVec_e[1], gVec_e[0]);

            oangs0[i, 0] = 2.0 * math.asin(sintht)
            oangs1[i, 0] = oangs0[i, 0]


#    return _transforms_CAPI.oscillAnglesOfHKLs(np.ascontiguousarray(hkls),chi,rMat_c,bMat,wavelength,
#                                               beamVec.flatten(),etaVec.flatten())

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
    d = abs(angList1 - angList0)
    return np.minimum(d, period - d)

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
