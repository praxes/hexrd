#! /usr/bin/env python
# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details, see https://github.com/joelvbernier/hexrd.
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
# orientation space stuff: orientations.py

import sys
import copy
import math
from math import sqrt

from numpy import *
import numpy as num
import scipy.linalg

import hexrd.matrixUtils as mU
from hexrd.matrixUtils import normvec3, normalized, cross, normvec

if __name__ != '__main__':
    debug = 0

# some parameters
# machAr = num.MachAr()
epsHedge     = 10.*num.finfo(float).eps
epsSqrt      = math.sqrt(num.finfo(float).eps)
tinyRotAng   = num.finfo(float).eps # 2e-16
nearPiRotMod = pi - epsSqrt # 1.0e-7
skewSymmMeth = 0.95 * pi
piby2 = pi * 0.5
piby3 = pi / 3.0
I3 = eye(3, dtype='float64')

def arccosSafe(temp):
    # protect against numbers slightly larger than 1 in magnitude 
    # due to round-off
    if temp > 1.00001:
        print >> sys.stderr, "attempt to take arccos of %s" % temp
        raise RuntimeError, "unrecoverable error"
    elif temp < -1.00001:
        print >> sys.stderr, "attempt to take arccos of %s" % temp
        raise RuntimeError, "unrecoverable error"
    elif temp >= 1.:
        ang = 0.
    elif temp <= -1.:
        ang = pi
    else:
        ang = math.acos(temp)
    return ang

def orthogonalize(rMatIn):
    UU = num.matrix(num.dot(rMatIn.T, rMatIn))
    U = scipy.linalg.sqrtm(UU).real
    rMat = num.array(scipy.linalg.solve(U.T, rMatIn.T), dtype=float).T
    return rMat

def traceToAng(tr):
    "given trace of a rotation matrix find the angle or rotation"
    temp = (tr-1.)*0.5
    return arccosSafe(temp)

def matToCanova(r):
    euler = zeros(3, dtype='float64')
                
    euler[1] = arccosSafe(r[2,2])

    if (abs(abs(r[2,2]) - 1.0) > epsHedge):
        sth = sin(euler[1])
        euler[2] = arctan2(r[0,2]/sth,  r[1,2]/sth)
        euler[0] = arctan2(r[2,0]/sth, -r[2,1]/sth)
    else:
        euler[2] = 0.
        euler[0] = arctan2(r[0,1],r[0,0])
    return euler

def invToRodr(inv):
    'do not check for divide-by-zero'
    vect = False
    if hasattr(inv, 'shape'):
        if len(inv.shape) == 2: vect = True
    if vect:
        r = zeros((3,inv.shape[1]), dtype='float64')
        r[:] = num.tan(inv[0,:] * 0.5) * inv[1:,:]
    else:
        r = zeros(3, dtype='float64')
        r[:] = tan(inv[0] * 0.5) * inv[1:]
    return r

def rodrToInv(rodr): 
    'do not check for divide-by-zero'
    vect = False
    sqr3i = 1./sqrt(3.)
    if hasattr(rodr, 'shape'):
        if len(rodr.shape) == 2: vect = True
    if vect:
        inv = num.zeros((4,rodr.shape[1]), dtype='float64')
        inv[1,:] = 1. # for cases really close to the identity
        #factor = num.apply_along_axis(normvec, 0, rodr)
        factor = num.apply_along_axis(num.linalg.norm, 0, rodr)
        inv[0,:] = 2.0 * num.arctan(factor[:])
        b = num.abs(factor) > tinyRotAng
        indices = b.nonzero()[0]
        inv[1:,indices] = rodr[:,indices] / factor[indices]
    else:
        inv = num.zeros(4, dtype='float64')
        factor = num.sqrt(num.sum(rodr[:] * rodr[:]))
        inv[0] = 2.0 * num.arctan(factor)
        if abs(inv[0]) <= tinyRotAng:
            inv[1:] = sqr3i
        else:
            inv[1:] = rodr[:] / factor
    return inv

def rodrToQuat(rodr): 
    return invToQuat(rodrToInv(rodr))

def invToQuat(inv):
    vect = False
    if hasattr(inv, 'shape'):
        if len(inv.shape) == 2: vect = True
    if vect:
        q = num.zeros((4,inv.shape[1]), dtype='float64')
        q[0,:]   = num.cos(inv[0,:] * 0.5)
        sth      = num.sin(inv[0,:] * 0.5)
        q[1:4,:] = sth * inv[1:,:]
    else:
        q = num.zeros(4, dtype='float64')
        q[0]   = cos(inv[0] * 0.5)
        sth    = sin(inv[0] * 0.5)
        q[1:4] = sth * inv[1:]
    return q
    
def bungeToMat(euler):
    
    mat = zeros([3,3], dtype='float64')

    c = num.cos(euler[:])
    s = num.sin(euler[:])
    
    mat[0, 0]  =  c[0]*c[2] - s[0]*c[1]*s[2]
    mat[1, 0]  =  s[0]*c[2] + c[0]*c[1]*s[2]
    mat[2, 0]  =  s[1]*s[2]
    mat[0, 1]  = -c[0]*s[2] - s[0]*c[1]*c[2]
    mat[1, 1]  = -s[0]*s[2] + c[0]*c[1]*c[2]
    mat[2, 1]  =  s[1]*c[2]
    mat[0, 2]  =  s[0]*s[1]
    mat[1, 2]  = -c[0]*s[1]
    mat[2, 2]  =  c[1]
    
    return mat

def matToQuat(r):
    """
    based on Spurrier's algorithm for quaternion extraction, as
    described in \cite{sim-vuq-85a} 

    returns a 4-vector, not a Quat instance

    BibTeX:
    @TechReport{sim-vuq-85a,
    author      = {J. C. Simo and L. {Vu Quoc}},
    title       = {Three dimensional finite strain rod model part
                   {II}: computational aspects, Memorandum
                   No. {UCB/ERL M85/31}},  
    institution = {Electronics Research Laboratory, College of
                   Engineering, University of California, Berkeley},
    year        = {1985}
    }
    """
    quat = zeros(4, dtype='float64')

    tr  = r[0,0]+r[1,1]+r[2,2]
    rDiag = array([r[0,0],r[1,1],r[2,2]])
    maxDiag = rDiag.max()
    
    if tr > maxDiag:
        quat[0] = 0.5 * sqrt(1 + tr)
        temp = 4.0 * quat[0]
        quat[1] = (r[2,1] - r[1,2]) / temp
        quat[2] = (r[0,2] - r[2,0]) / temp
        quat[3] = (r[1,0] - r[0,1]) / temp
    else:
        mIndex = where(rDiag == rDiag.max())
        mi = mIndex[0][0]
        mj = mi+1
        mk = mi+2
        if mj > 2: mj = mj - 3
        if mk > 2: mk = mk - 3
        miq = mi+1
        mjq = mj+1
        mkq = mk+1
        quat[miq] = sqrt(r[mi,mi]*0.5 + (1.0 - tr)*.25)
        temp = 4.0*quat[miq]
        quat[0]   = (r[mk,mj] - r[mj,mk])/temp
        quat[mjq] = (r[mj,mi] + r[mi,mj])/temp
        quat[mkq] = (r[mk,mi] + r[mi,mk])/temp
    
    return quat

def matToThetaN(r):
    """2nd order tensor => angle/axis

    to see that angle is right take R in a basis so that R11=1;
    get tr(R) = 1 + 2 cos(theta), solve for theta and use argument
    about invariance of tr(R)

    references:
    1) box 4 in Simo and VuQuoc, 1985, ERL Berkeley memorandum 
    no. UCB/ERL M85/31 
    2) Marin and Dawson 98 part 1, 
    equation for update d_rstar (exponential mapping)

    n is in vector notation of a skew tensor according to
    w_i = 1/2 epsilon_jik W_jk"""

    theta = 0.
    n     = zeros(3, dtype='float64')

    sqr3i = 1./sqrt(3.)

    ra = zeros(3, dtype='float64')

    #REAL(idp) :: factor, trR1, ra(DIMS), temp
    #LOGICAL :: sign_pos
    #INTEGER :: max_ax

    trR1  = r[0,0]+r[1,1]+r[2,2]
    theta = traceToAng(trR1)
    trR1  = trR1 - 1.

    if theta <= tinyRotAng:
        # R is or is nearly the identity
        n[:] = sqr3i

    elif theta <= skewSymmMeth:
        # determine axis of R from skew part of R

        factor = 1.0 / sin(theta)
        n[0] = 0.5 * (r[2,1] - r[1,2]) * factor
        n[1] = 0.5 * (r[0,2] - r[2,0]) * factor
        n[2] = 0.5 * (r[1,0] - r[0,1]) * factor

    else:
        # determine magnitude of axis components from symm part of R
        #
        # have already ruled out the identity and small rotations 
        # (close to identity) handled by theta .LE. skew_symm_meth 
        # case above, so this factor is okay
        factor = 1./(2.- trR1) # 1/(3-tr(R))

        # protect against taking sqrt of a -zero number
        temp = (2.*r[0,0] - trR1)*factor
        if temp < 1e-8:
            n[0] = 0.
        else:
            n[0] = sqrt(temp)

        temp = (2.*r[1,1] - trR1)*factor
        if temp < 1e-8:
            n[1] = 0.
        else:
            n[1] = sqrt(temp)

        temp = (2.*r[2,2] - trR1)*factor
        if temp < 1e-8:
            n[2] = 0.
        else:
            n[2] = sqrt(temp)

        # determine sign of components of axis, use skew(R) if large
        # enough, otherwise assume that rotation is by pi and thus
        # the direction of the axis is arbitrary so set the largest
        # component of the axis positive
        ra[0] = (r[2,1] - r[1,2])
        ra[1] = (r[0,2] - r[2,0])
        ra[2] = (r[1,0] - r[0,1])

        max_ax = argmax(abs(ra))
        if abs(ra[max_ax]) > tinyRotAng:
            if ra[max_ax] < 0.:
                n[max_ax] = -n[max_ax]
                sign_pos = 0
            else:
                sign_pos = 1
        else:
            max_ax = argmax(n)
            sign_pos = 1

        # will have a nonzero symmetric part of R for use below unless R=I
        # which has already been taken care of or R is a rotation by
        # pi about one of the coordinate directions in which case
        # only one of the axis components should be nonzero and the
        # sign of it does not matter anyway
        if (max_ax != 0) and (n[0] > 0.):
            # get sign for axis(1)
            if sign_pos != (r[max_ax,0]+r[0,max_ax] > 0.) :
                n[0] = -n[0]
        if (max_ax != 1) and (n[1] > 0.):
            # get sign for axis(2)
            if sign_pos != (r[max_ax,1]+r[1,max_ax] > 0.) :
                n[1] = -n[1]
        if (max_ax != 2) and (n[2] > 0.):
            # get sign for axis(3)
            if sign_pos != (r[max_ax,2]+r[2,max_ax] > 0.) :
                n[2] = -n[2]

    return theta, n

def quatToInv(q):
    invParams = zeros(4, dtype='float64')
    invParams[1:] = q[1:4]
    period = num.pi * 2.
    mag = sqrt(sum(invParams[1:] * invParams[1:]))
    if abs(q[0]) > 0.9:
        theta = 2.0 * arcsin(mag)
    else:
        theta = 2.0 * arccosSafe(q[0])
    theta = mod(theta + 0.5*period, period) - 0.5*period # put in -pi to pi, useful if have thetaScale stuff going on
    invParams[0] = theta
    if mag > 0:
        invParams[1:] = invParams[1:] / mag
    else:
        invParams[1:] = [1., 0., 0.]
    return invParams

def quatToMat(quat):
    """
    Take an array of n quats (numpy ndarray, 4 x n) and generate an
    array of rotation matrices (n x 3 x 3)

    Uses the truncated series expansion for the exponential map;
    divide-by-zero is checked using the global 'tinyRotAng'
    """
    if quat.shape[0] != 4:
        raise RuntimeError, "input is the wrong shape"
    
    n    = quat.shape[1]
    rmat = zeros((n, 3, 3), dtype='float64')
    
    for i in range(n):
        theta = 2. * arccosSafe(quat[0, i])
        
        # find axial vector
        if (theta > tinyRotAng):
            a = sin(theta) / theta
            b = (1. - cos(theta)) / (theta*theta)
            w = (theta / sin(0.5 * theta)) * quat[1:4, i]
            
            wskew = array([[   0., -w[2],  w[1]],
                           [ w[2],    0., -w[0]],
                           [-w[1],  w[0],    0.]])
            
            rmat[i, :, :] = I3 + a * wskew + b * dot(wskew, wskew)
        else:
            rmat[i, :, :] = I3
    
    return rmat

def quatToProdMat(quat, mult='right'):
    """
    Form 4 x 4 arrays to perform the quaternion product
    
    USAGE
        qmats = quatToProdMat(quats, mult='right')
        
    INPUTS
        1) quats is (4, n), a numpy ndarray array of n quaternions
           horizontally concatenated 
        2) mult is a keyword arg, either 'left' or 'right', denoting
           the sense of the multiplication:
            
                      / quatToProdMat(h, 'right') * q
           q * h  --> <
                      \ quatToProdMat(q, 'left') * h 

    OUTPUTS
        1) qmats is (n, 4, 4), the left or right quaternion product
           operator 

    NOTES
       *) This function is intended to replace a cross-product based
          routine for products of quaternions with large arrays of
          quaternions (e.g. applying symmetries to a large set of
          orientations).
    """

    if quats.shape[0] != 4:
        raise RuntimeError, "input is the wrong size along the 0-axis"
    
    nq = quats.shape[1]
    
    q0 = quats[0, :].copy()
    q1 = quats[1, :].copy()
    q2 = quats[2, :].copy()
    q3 = quats[3, :].copy()
    
    if mult == 'right':
        qmats = array([[ q0], [ q1], [ q2], [ q3], 
                       [-q1], [ q0], [-q3], [ q2],
                       [-q2], [ q3], [ q0], [-q1],
                       [-q3], [-q2], [ q1], [ q0]])
    elif mult == 'left':
        qmats = array([[ q0], [ q1], [ q2], [ q3], 
                       [-q1], [ q0], [ q3], [-q2],
                       [-q2], [-q3], [ q0], [ q1],
                       [-q3], [ q2], [-q1], [ q0]])
    
    # some fancy reshuffling:
    qmats = transpose(reshape(qmats.T, (nq, 4, 4)), (0, 2, 1))
    
    return qmats

def sampleToLatticeT2(A_sm, C):
    '''
                        T
    [A_sm]=[C][A_lat][C]
    '''
    A_lat = num.matrix(C).T * num.matrix(A_sm) * num.matrix(C)
    return A_lat

def latticeToSampleT2(A_lat, C):
    '''
                        T
    [A_sm]=[C][A_lat][C]
    '''
    A_sm = num.matrix(C) * num.matrix(A_lat) * num.matrix(C).T
    return A_sm

def latticeToSampleV(V_lat, C):
    V_sm = num.dot(C, V_lat)
    return V_sm

class RotationParameterization:
    "template for rotation parameterization class"
    def __init__(self, args):
        raise RuntimeError, "need implementation"
    def toMatrix(self):
        """use matrix as common representation
        all classes should have a constructor that works using this
        
            **  Construct [C] matrix (Kocks convention)
            **         {a}      = [C] {a}
            **            sm             cr
            """
        raise RuntimeError, "need implementation"

def rotMatrixFromCrystalVectors(cvs1=None, cvs2=None, cvs3=None):
    """\
    Make a rotation matrix in the RotationParameterization convention
    from components of crystal vectors that are along given sample directions
    """
    mat = zeros([3,3], dtype='float64')
    if cvs1 == None:
        if ((cvs2 == None) or (cvs3 == None)):
            print >> sys.stderr, "need more inputs"
            raise RuntimeError, "unrecoverable error"
        cvs1 = cross(cvs2,cvs3)
        cvs1 = normalized(cvs1)
    if cvs2 == None:
        if ((cvs1 == None) or (cvs3 == None)):
            print >> sys.stderr, "need more inputs"
            raise RuntimeError, "unrecoverable error"
        cvs2 = cross(cvs1,cvs3)
        cvs2 = normalized(cvs2)
    if cvs3 == None:
        if ((cvs1 == None) or (cvs2 == None)):
            print >> sys.stderr, "need more inputs"
            raise RuntimeError, "unrecoverable error"
        cvs3 = cross(cvs1,cvs2)
        cvs3 = normalized(cvs3)
    mat[0,:] = cvs1[:]
    mat[1,:] = cvs2[:]
    mat[2,:] = cvs3[:]
    det = determinant3(mat)
    if (abs(det-1.0) < 1e-12):
        return mat
    elif (abs(det+1.0) < 1e-12):
        mat = -mat
        return mat
    else:
        print >> sys.stderr, "vectors not close enough to orthonormal"
        raise RuntimeError, "unrecoverable error"

class RotInv(RotationParameterization):
    "rotation invariants"
    #__rvNorm = stats.norm().rsv # scipy
    __rvNorm = random.normal # numpy
    def __init__(self, *args):
        if isinstance(args[0], str):
            a = args[0]
            if a == 'rand':
                if len(args) > 1:
                    if hasattr(args[1], '__len__'):
                        thetaScale = args[1]
                        if len(thetaScale) == 3:
                            self.n    = zeros(3, dtype='float64')
                            self.n[0] = RotInv.__rvNorm() * thetaScale[0]
                            self.n[1] = RotInv.__rvNorm() * thetaScale[1]
                            self.n[2] = RotInv.__rvNorm() * thetaScale[2]
                            mag = normvec3(self.n)
                            if mag > 0.:
                                self.n[:] = self.n[:] / mag
                                self.theta = mag
                            else:
                                newInv = Quat('rand').toInv()
                                self.n[:] = newInv.n[:]
                                self.theta = mag
                        else:
                            raise RuntimeError, "unrecoverable error: bad thetaScale length"
                    else:
                        newInv = Quat('rand').toInv()
                        self.n     = newInv.n
                        thetaScale = float(args[1])
                        self.theta = RotInv.__rvNorm() * thetaScale
                else:
                    newInv = Quat('rand').toInv()
                    self.n     = newInv.n
                    self.theta = newInv.theta
            elif a == 'align':
                # align given crystal vector with given sample vector
                vl = num.array(args[1])
                vs = num.array(args[2])
                # normalize
                mag = normvec3(vl)
                vl /= mag
                mag = normvec3(vs)
                vs /= mag
                self.theta = math.acos(num.dot(vs,vl)) # arccosine
                w = cross(vl,vs)
                mag = normvec3(w)
                if mag > 1e-14:
                    self.n = w / mag
                else:
                    sqr3i = 1./sqrt(3.)
                    self.n = num.array([sqr3i, sqr3i, sqr3i])
            else:
                print >> sys.stderr, "bad string"
                raise RuntimeError, "unrecoverable error"
        elif len(args) == 1:
            a = args[0]
            if (hasattr(a,"toInv")): 
                "use direct conversion, assuming it is more efficient"
                newInv = a.toInv()
                self.theta = newInv.theta
                self.n     = newInv.n
            elif isinstance(a,RotationParameterization):
                "from 3x3 matrix of components"
                (self.theta, self.n) = matToThetaN(a.toMatrix())
            else:
                a = num.atleast_1d(args[0])
                if size(shape(a)) == 2:
                    "from 3x3 matrix of components"
                    assert a.shape[0] == 3 and a.shape[1] == 3, 'wrong shape for matrix'
                    (self.theta, self.n) = matToThetaN(a)
                elif size(shape(a)) == 1:
                    "exponential map parameters"
                    assert a.shape[0] == 3, 'wrong shape for exp map parameters'
                    self.n     = a
                    self.theta = sqrt(sum(self.n[:] * self.n[:]))
                    if self.theta > 0.:
                        self.n     = self.n / self.theta
                    else:
                        self.n     = num.array([1.,0.,0.]) # n not really well defined
                else:
                    print >> sys.stderr, "wrong shape of args[0]", shape(a)
                    raise RuntimeError, "unrecoverable error"
        elif len(args) == 2:
            "from theta and normal vector (not necessarily normalized)"
            self.theta = float(args[0])
            self.n     = zeros(3, dtype='float64')
            self.n[:]  = args[1]
            mag = sqrt(sum(self.n[:] * self.n[:]))
            self.n[:] = self.n[:] / mag
        elif len(args) == 4:
            "from theta and normal vector (not necessarily normalized)"
            self.theta = float(args[0])
            self.n     = zeros(3, dtype='float64')
            self.n[:]  = args[1:4]
            mag = sqrt(sum(self.n[:] * self.n[:]))
            self.n[:] = self.n[:] / mag
        else:
            print >> sys.stderr, "wrong number of args %g" % len(args)
            raise RuntimeError, "unrecoverable error"
    def __repr__(self):
        return "RotInv(%s,[%s,%s,%s])" % (self.theta,self.n[0],self.n[1],self.n[2])
    def __str__(self):
        return "(%g, [%g, %g, %g])" % (self.theta,self.n[0],self.n[1],self.n[2])
    def toQuat(self):
        q = Quat(invToQuat(num.hstack((self.theta, self.n)))) # [self.theta]+list(self.n)
        if debug > 1: print "from theta, n:", self.theta, self.n
        if debug > 1: print "q:", q
        return q
    def toMatrix(self):
        mat = zeros([3,3], dtype='float64')
        a = cos(self.theta)
        b = sin(self.theta)
        am1 = (1. - a)
        mat[0,0] = am1 * self.n[0] * self.n[0] + a
        mat[0,1] = am1 * self.n[0] * self.n[1] - self.n[2] * b
        mat[0,2] = am1 * self.n[0] * self.n[2] + self.n[1] * b
        mat[1,0] = am1 * self.n[0] * self.n[1] + self.n[2] * b
        mat[1,1] = am1 * self.n[1] * self.n[1] + a  
        mat[1,2] = am1 * self.n[1] * self.n[2] - self.n[0] * b
        mat[2,0] = am1 * self.n[0] * self.n[2] - self.n[1] * b
        mat[2,1] = am1 * self.n[1] * self.n[2] + self.n[0] * b
        mat[2,2] = am1 * self.n[2] * self.n[2] + a
        return mat

class CanovaEuler(RotationParameterization):
    def __init__(self, *args):
        if len(args) == 3:
            self.euler = zeros(3, dtype='float64')
            self.euler[:] = args[:]
        elif len(args) == 1:
            a = args[0]
            self.euler = zeros(3, dtype='float64')
            if isinstance(a,RotationParameterization):
                if hasattr(a, 'toCanova'):
                    self.euler = a.toCanova()
                else:
                    "from 3x3 matrix of components"
                    self.euler[:] = CanovaEuler(a.toMatrix())
            if size(shape(args[0])) == 2:
                "from 3x3 matrix of components"
                self.euler[:] = matToCanova(args[0])
            else:
                print >> sys.stderr, "wrong shape of args[0]", shape(args[0])
                raise RuntimeError, "unrecoverable error"
        else:
            print >> sys.stderr, "wrong number of args %g" % len(args)
            raise RuntimeError, "unrecoverable error"
    def __repr__(self):
        return "CanovaEuler(%s,%s,%s)" % (self.euler[0],self.euler[1],self.euler[2])
    def __str__(self):
        return "(%g %g %g)" % (self.euler[0],self.euler[1],self.euler[2])
    def toKocks(self):
        e_k1 = piby2 - self.euler[2]
        e_k2 = self.euler[1]
        e_k3 = self.euler[0] - piby2
        kocks = KocksEuler(e_k1, e_k2, e_k3)
        return kocks
    def toMatrix(self):
        return (self.toKocks()).toMatrix()
    
class KocksEuler(RotationParameterization):
    "Kocks Euler angles see equation 6 page 65 of koc-tom-wen-98a (Kocks, Tome, & Wenk; Texture and Anisotropy)"
    def __init__(self, *args):
        if len(args) == 3:
            self.euler = zeros(3, dtype='float64')
            self.euler[:] = args[:]
        elif len(args) == 1:
            a = args[0]
            if isinstance(a,RotationParameterization):
                self.euler = zeros(3, dtype='float64')
                if hasattr(a, 'toKocks'):
                    self.euler[:] = a.toKocks()
                else:
                    "from 3x3 matrix of components"
                    canova = CanovaEuler(a.toMatrix())
                    kocks = canova.toKocks()
                    self.euler[:] = kocks.euler
            elif size(shape(a)) == 2:
                "from 3x3 matrix of components"
                self.euler = zeros(3, dtype='float64')
                canova = CanovaEuler(a)
                kocks = canova.toKocks()
                self.euler = kocks.euler
            else:
                print >> sys.stderr, "wrong shape of args[0]", shape(a)
                raise RuntimeError, "unrecoverable error"
        else:
            print >> sys.stderr, "wrong number of args %g" % len(args)
            raise RuntimeError, "unrecoverable error"
    def __repr__(self):
        return "KocksEuler(%s,%s,%s)" % (self.euler[0],self.euler[1],self.euler[2])
    def __str__(self):
        return "(%g %g %g)" % (self.euler[0],self.euler[1],self.euler[2])
    def toBunge(self):
        'trusting Table 1a of koc-tom-wen-98a (Kocks, Tome, & Wenk; Texture and Anisotropy)'
        e_b1 =  self.euler[0] + piby2 # phi1
        e_b2 =  self.euler[1] # Phi
        e_b3 = -self.euler[2] + piby2 # phi2
        bunge = BungeEuler(e_b1, e_b2, e_b3)
        return bunge
    def toMatrix(self):
        mat = zeros([3,3], dtype='float64')
        
        sps = sin(self.euler[0])
        cps = cos(self.euler[0])
        sth = sin(self.euler[1])
        cth = cos(self.euler[1])
        sph = sin(self.euler[2])
        cph = cos(self.euler[2])

        mat[0,0] = -sps * sph - cps * cph * cth
        mat[1,0] =  cps * sph - sps * cph * cth
        mat[2,0] =  cph * sth
        mat[0,1] =  cph * sps - sph * cps * cth
        mat[1,1] = -cps * cph - sps * sph * cth
        mat[2,1] =  sph * sth
        mat[0,2] =  cps * sth
        mat[1,2] =  sps * sth
        mat[2,2] =  cth

        return mat
        
class BungeEuler(RotationParameterization):
    '''Bunge Euler angles
    see koc-tom-wen-98a (Kocks, Tome, & Wenk; Texture and Anisotropy)
    '''
    def __init__(self, *args):
        if len(args) == 3:
            self.euler = zeros(3, dtype='float64')
            self.euler[:] = args[:]
        elif len(args) == 1:
            a = args[0]
            if isinstance(a,RotationParameterization):
                self.euler = zeros(3, dtype='float64')
                if hasattr(a, 'toBunge'):
                    self.euler[:] = a.toBunge()
                else:
                    "from 3x3 matrix of components"
                    bunge = CanovaEuler(a.toMatrix()).toKocks().toBunge()
                    self.euler[:] = bunge.euler
            elif size(shape(a)) == 2:
                "from 3x3 matrix of components"
                self.euler = zeros(3, dtype='float64')
                bunge = CanovaEuler(a).toKocks().toBunge()
                self.euler = bunge.euler
            else:
                print >> sys.stderr, "wrong shape of args[0]", shape(a)
                raise RuntimeError, "unrecoverable error"
        else:
            print >> sys.stderr, "wrong number of args %g" % len(args)
            raise RuntimeError, "unrecoverable error"
    def __repr__(self):
        return "BungeEuler(%s,%s,%s)" % (self.euler[0],self.euler[1],self.euler[2])
    def __str__(self):
        return "%24.16e %24.16e %24.16e" % (self.euler[0],self.euler[1],self.euler[2])
    def toKocks(self):
        'trusting Table 1a of koc-tom-wen-98a (Kocks, Tome, & Wenk; Texture and Anisotropy)'
        e_k1 = self.euler[0] - piby2
        e_k2 = self.euler[1]
        e_k3 = piby2         - self.euler[2] 
        kocks = KocksEuler(e_k1, e_k2, e_k3)
        return kocks
    def toMatrix(self):
        #return self.toKocks().toMatrix()
        mat = bungeToMatrix(self.euler)
        return mat
        
# // Find the quaternion for this q vector and its canonical representation
# Quaternion rot = quatFinder.find(qhat, L);
# rot = findCanonical(rot);
# 
# // Quaternion components are (r, i, j, k).  A Vector is a standard
# // vector in 3-space.
# Vector v(rot.i(), rot.j(), rot.k());
# v /= rot.r();
# double alpha = 2.0* std::atan(v.length());
# v.normalize();
# v *= alpha *0.5;
# 
# v += Vector(0.5, 0.5, 0.5);
# 
# PIX color;
# 
# color.r = short( 255*fabs(v[0]*v[0]) );
# color.g = short( 255*fabs(v[1]*v[1]) );
# color.b = short( 255*fabs(v[2]*v[2]) );
class Quat(RotationParameterization):
    "quaternions, normalized; parameterization of SO(3); do not bother making unique (is 2-to-1)"
    #__rvNorm = stats.norm().rsv # scipy
    __rvNorm = random.normal # numpy
    def __init__(self, *args):
        self.q = zeros(4, dtype='float64')
        # self.symmPath = []
        self.symmVariant = []
        self.symmVariantPath = []
        if len(args) == 4:
            self.q[0] = float(args[0])
            self.q[1] = float(args[1])
            self.q[2] = float(args[2])
            self.q[3] = float(args[3])
            self.normalize()
        elif len(args) == 1:
            a = args[0]
            if isinstance(a,RotationParameterization):
                if (hasattr(a,"toQuat")):
                    "use direct conversion, assuming it is more efficient"
                    newQ = a.toQuat()
                    self.q = newQ.q
                else:
                    "from 3x3 matrix of components"
                    #inv = RotInv(a.toMatrix())
                    #newQ = inv.toQuat()
                    #self.q = newQ.q
                    self.q = matToQuat(a.toMatrix())
            elif isinstance(a, str):
                if a == 'rand':
                    "generate a random quaternion -- a sample from a uniform distribution in SO(3)"
                    self.q = Quat.getRandQuat()
                    #self.q = Quat.__rvNorm(size=4)
                    # #self.q[0] = float(Quat.__rvNorm.rvs())
                    # #self.q[1] = float(Quat.__rvNorm.rvs())
                    # #self.q[2] = float(Quat.__rvNorm.rvs())
                    # #self.q[3] = float(Quat.__rvNorm.rvs())
                else:
                    print >> sys.stderr, "bad string"
                    raise RuntimeError, "unrecoverable error"
            else:
                if size(a) == 4:
                    self.q[0] = float(a[0])
                    self.q[1] = float(a[1])
                    self.q[2] = float(a[2])
                    self.q[3] = float(a[3])
                elif size(shape(a)) == 2:
                    # assume 3x3 matrix of components
                    inv = RotInv(a)
                    if debug > 1: print "got: ",inv
                    newQ = inv.toQuat()
                    self.q = newQ.q
                else:
                    print >> sys.stderr, "bad size"
                    raise RuntimeError, "unrecoverable error"
            self.normalize()
        elif len(args) == 0:
            self.q[0] = 1.
            self.q[1] = 0.
            self.q[2] = 0.
            self.q[3] = 0.
        else:
            print >> sys.stderr, "wrong number of args"
            raise RuntimeError, "unrecoverable error"
    @staticmethod
    def normalizeQuat(q):
        'normalize in place'
        temp = q[:] * q[:]
        mag = sqrt(sum(temp))
        if (mag < 1e-6):
            print >> sys.stderr, "quat magnitude too small"
            raise RuntimeError, "unrecoverable error"
        q = q / mag
        return q
    @staticmethod
    def getRandQuat(n=1):
        'sample uniform orientation distribution to get quaternion parameters'
        if n == 1:
            q = Quat.__rvNorm(size=4)
            retval = Quat.normalizeQuat(q)
        else:
            retval = []
            for iN in range(n):
                q = Quat.__rvNorm(size=4)
                q = Quat.normalizeQuat(q)
                retval.append(q)
            retval = num.array(retval)
        return retval
    def __repr__(self):
        return "Quat(%s,%s,%s,%s)" % (self.q[0],self.q[1],self.q[2],self.q[3])
    def __str__(self):
        #qString = "(%g, %g, %g, %g)" % (self.q[0],self.q[1],self.q[2],self.q[3])
        #qString = "%s %s %s %s" % (self.q[0],self.q[1],self.q[2],self.q[3])
        qString = "%24.16e %24.16e %24.16e %24.16e" % (self.q[0],self.q[1],self.q[2],self.q[3])
        if len(self.symmVariant) > 0:
            qString += " symmVariant: " + self.symmVariant.__str__()
        if len(self.symmVariantPath) > 0:
            qString += " symmVariantPath: " + self.symmVariantPath.__str__()
        return qString
    def __mul__(self,other):
        "multiplication; if multiply a list of Quats then list keeps metadata"
        if isinstance(other, Quat):
            new = Quat()
            new.q[0] = self.q[0]*other.q[0] - self.q[1]*other.q[1] - self.q[2]*other.q[2] - self.q[3]*other.q[3]
            new.q[1] = self.q[0]*other.q[1] + self.q[1]*other.q[0] + self.q[2]*other.q[3] - self.q[3]*other.q[2]
            new.q[2] = self.q[0]*other.q[2] - self.q[1]*other.q[3] + self.q[2]*other.q[0] + self.q[3]*other.q[1]
            new.q[3] = self.q[0]*other.q[3] + self.q[1]*other.q[2] - self.q[2]*other.q[1] + self.q[3]*other.q[0]
            return new
        elif isinstance(other, SymmGroup):
            return other.__rmul__(self)
        elif hasattr(other, '__len__'):
            if ((len(self.symmVariant) > 0) or (len(self.symmVariantPath) > 0)):
                print >> sys.stderr, "quat multiplying list and quat has metadata"
                raise RuntimeError, "unrecoverable error"
            result = copy.deepcopy(other)
            for iR in range(len(other)):
                tempQuat = self * result[iR]
                result[iR].q[:] = tempQuat.q[:]
            return result
        else:
            raise RuntimeError, "mul with unsupported operand"
    def __coerce__(self, other):
        if isinstance(other, Quat):
            return self, other
        if hasattr(other, '__len__'):
            if len(other) == 4:
                return self,Quat(other[0],other[1],other[2],other[3])
            else:
                # do not know what else to do
                return self, other
        else:
            # do not know what else to do
            return self, other
    def strip(self):
        self.symmVariant = []
        self.symmVariantPath = []
    def T(self):
        "return transposed quaternion"
        return Quat(self.q[0],-self.q[1],-self.q[2],-self.q[3])
    def transposed(self):
        "return transposed quaternion"
        return self.T()
    def transpose(self):
        "transpose the quaternion in place; retains metadata"
        self.q[1:4] = -self.q[1:4]
        return 
    def normalized(self):
        "return normalize quaternion"
        temp = self.q[:] * self.q[:]
        mag = sqrt(sum(temp))
        if (mag < 1e-6):
            print >> sys.stderr, "quat magnitude too small"
            raise RuntimeError, "unrecoverable error"
        new = Quat()
        new.q = self.q[:] / mag
        return new
    def normalize(self):
        "normalize the quaternion"
        #temp = self.normalized()
        #self.q[:] = temp.q[:]
        self.q = Quat.normalizeQuat(self.q)
        return
    def misorAng(self, other):
        "compute the misorientation angle, in radians, in [0,2*pi]"
        oT = other.T()
        q0 = self.q[0]*oT.q[0] - self.q[1]*oT.q[1] - self.q[2]*oT.q[2] - self.q[3]*oT.q[3]
        ang = 2.0 * arccosSafe(q0)
        return ang
    def toInv(self):
        invParams = quatToInv(self.q)
        return RotInv(*invParams)
    def toMatrix(self):
        mat = zeros([3,3], dtype='float64')
        
        x0sq = self.q[0]*self.q[0]
        x1sq = self.q[1]*self.q[1]
        x2sq = self.q[2]*self.q[2]
        x3sq = self.q[3]*self.q[3]
        
        x0x1 = self.q[0]*self.q[1]
        x0x2 = self.q[0]*self.q[2]
        x0x3 = self.q[0]*self.q[3]
        
        x1x2 = self.q[1]*self.q[2]
        x1x3 = self.q[1]*self.q[3]
        
        x2x3 = self.q[2]*self.q[3]
        
        mat[0,0] = x0sq+x1sq-x2sq-x3sq
        mat[0,1] = 2.0*(x1x2-x0x3)
        mat[0,2] = 2.0*(x1x3+x0x2)
        mat[1,0] = 2.0*(x1x2+x0x3)
        mat[1,1] = x0sq-x1sq+x2sq-x3sq
        mat[1,2] = 2.0*(x2x3-x0x1)
        mat[2,0] = 2.0*(x1x3-x0x2)
        mat[2,1] = 2.0*(x2x3+x0x1)
        mat[2,2] = x0sq-x1sq-x2sq+x3sq
        
        return mat
        
class SymmGroup:
    "symmetry group"
    __misorZeroTol = 1e-6
    def __init__(self, *args):
        if len(args) == 1:
            a = args[0]
            if isinstance(a, str):
                if a == 'cubic':
                    self.nSymm = int(24)
                    self.qSymm = [ Quat() for i in range(self.nSymm)]

                    self.setQuatI(0,  (  1.,  0.,  0.,  0.  )  ) #  11 
                    self.setQuatI(1,  (  1.,  1.,  1.,  1.  )  ) #  12 
                    self.setQuatI(2,  ( -1.,  1.,  1.,  1.  )  ) #  13 
                    self.setQuatI(3,  (  1.,  0.,  0.,  1.  )  ) #  21 
                    self.setQuatI(4,  (  0.,  0.,  1.,  1.  )  ) #  22 
                    self.setQuatI(5,  ( -1.,  0.,  1.,  0.  )  ) #  23 
                    self.setQuatI(6,  (  0.,  0.,  0.,  1.  )  ) #  31 
                    self.setQuatI(7,  ( -1., -1.,  1.,  1.  )  ) #  32 
                    self.setQuatI(8,  ( -1., -1.,  1., -1.  )  ) #  33 
                    self.setQuatI(9,  ( -1.,  0.,  0.,  1.  )  ) #  41 
                    self.setQuatI(10, ( -1., -1.,  0.,  0.  )  ) #  42 
                    self.setQuatI(11, (  0., -1.,  0., -1.  )  ) #  43 
                    self.setQuatI(12, (  0.,  1.,  1.,  0.  )  ) #  51 
                    self.setQuatI(13, ( -1.,  1.,  0.,  0.  )  ) #  52 
                    self.setQuatI(14, (  1.,  0.,  1.,  0.  )  ) #  53 
                    self.setQuatI(15, (  0.,  0., -1.,  0.  )  ) #  61 
                    self.setQuatI(16, (  1., -1., -1.,  1.  )  ) #  62 
                    self.setQuatI(17, (  1., -1.,  1.,  1.  )  ) #  63 
                    self.setQuatI(18, (  0.,  1., -1.,  0.  )  ) #  71 
                    self.setQuatI(19, (  0.,  0., -1.,  1.  )  ) #  72 
                    self.setQuatI(20, (  0., -1.,  0.,  1.  )  ) #  73 
                    self.setQuatI(21, (  0.,  1.,  0.,  0.  )  ) #  81 
                    self.setQuatI(22, ( -1.,  1., -1.,  1.  )  ) #  82 
                    self.setQuatI(23, ( -1., -1., -1.,  1.  )  ) #  83 
    
                    self.symmName = "cubic"

                elif a == 'hexag':
                    self.nSymm = int(12)
                    self.qSymm = [ Quat() for i in range(self.nSymm)]

                    halfsqr3 = sqrt(3.)*0.5
    
                    self.setQuatI(0,  ( 1.,     0.,      0.,      0.      )  ) #  00
                    self.setQuatI(1,  ( 0.,     1.,      0.,      0.      )  ) #  20
                    self.setQuatI(2,  ( 0.,     halfsqr3,  0.5,   0.      )  ) #  21
                    self.setQuatI(3,  ( 0.,     0.5,   halfsqr3,  0.      )  ) #  22
                    self.setQuatI(4,  ( 0.,     0.,      1.,      0.      )  ) #  23
                    self.setQuatI(5,  ( 0.,    -0.5,   halfsqr3,  0.      )  ) #  24
                    self.setQuatI(6,  ( 0.,    -halfsqr3,  0.5,   0.      )  ) #  25
                    self.setQuatI(7,  ( halfsqr3, 0.,      0.,      0.5   )  ) #  60
                    self.setQuatI(8,  ( 0.5,  0.,      0.,      halfsqr3  )  ) #  61
                    self.setQuatI(9,  ( 0.,    0.,      0.,      1.       )  ) #  62
                    self.setQuatI(10, ( halfsqr3, 0.,      0.,     -0.5   )  ) #  63
                    self.setQuatI(11, ( 0.5,  0.,      0.,     -halfsqr3  )  ) #  64
    
                    self.symmName = "hexag"

                else:
                    print >> sys.stderr, "unknown symmetry type"
                    raise RuntimeError, "unrecoverable error"
            else:
                self.nSymm = int(a)
                self.qSymm = [ Quat() for i in range(a)]
                self.symmName = "unInit"
        else:
            print >> sys.stderr, "wrong number of args"
            raise RuntimeError, "unrecoverable error"
    def setQuatI(self,i,a):
        self.qSymm[i].q[:] = a[:]
        self.qSymm[i].normalize()
    def __rmul__(self, q):
        "post-multiply a quaternion or list of quaternions by all symmetry elements"
        if hasattr(q, '__len__'):
            result = [ Quat() for i in range(self.nSymm*len(q))]
            i = 0
            for iSymm in range(self.nSymm):
                qSymmi = self.qSymm[iSymm]
                for qi in q:
                    result[i] = qi * qSymmi
                    # result[i].symmPath.append([(self, iSymm)])
                    if (len(qi.symmVariant)>0):
                        result[i].symmVariantPath.append((qi.symmVariant[-1:][0], iSymm))
                    i += 1
            return result
        else:
            result = [ Quat() for i in range(self.nSymm)]
            for iSymm in range(self.nSymm):
                result[iSymm] = q * self.qSymm[iSymm]
                # result[iSymm].symmPath.append([(self, iSymm)])
                if (len(q.symmVariant)>0):
                    result[iSymm].symmVariantPath.append((q.symmVariant[-1:][0], iSymm))
            return result
    def findDistinct(self, qList):
        "given a list of quaternions, find those which are distinct under the symmetry group"
        distinct =  ones(len(qList), dtype='i')
        slavedTo = zeros(len(qList), dtype='i')
        for i in range(len(qList)):
            if not distinct[i]: continue
            for j in range(i+1,len(qList)):
                if not distinct[j]: continue
                for k in range(self.nSymm):
                    # if not distinct[j]: break # not needed
                    """
                    tempa = qList[j] * self.qSymm[k]
                    tempb = qList[i] # * self.qSymm[l]
                    misAng = tempb.misorAng(tempa)
                    if debug :
                        #print "misAng between %g and %g with symms %g,%g: %g" % (i, j, k, l, misAng)
                        print "misAng between %g and %g with symm %g: %g" % (i, j, k, misAng)
                    if misAng < SymmGroup.__misorZeroTol:
                        distinct[j] = 0
                        slavedTo[j] = i
                        break
                    """
                    for l in range(self.nSymm):
                        tempa = qList[j] * self.qSymm[k]
                        tempb = qList[i] * self.qSymm[l]
                        misAng = tempb.misorAng(tempa)
                        if debug :
                            print "misAng between %g and %g with symms %g,%g: %g" % (i, j, k, l, misAng)
                        if misAng < SymmGroup.__misorZeroTol:
                            distinct[j] = 0
                            slavedTo[j] = i
                            break
        result = []
        if debug >= 0:
            print "distinct: ", distinct
        for i in range(len(qList)):
            if distinct[i]:
                distinct[i] = len(result)
                result.append(copy.deepcopy(qList[i]))
                result[distinct[i]].symmVariant.append(distinct[i])
            else:
                # should already have changed distinct[slavedTo[i]] to have mapping
                if len(qList[i].symmVariantPath) == 1:
                    result[distinct[slavedTo[i]]].symmVariantPath.append(qList[i].symmVariantPath[0])
                elif len(qList[i].symmVariantPath) == 0:
                    pass # nothing to do
                else:
                    # do not know what to do with this case
                    print >> sys.stderr, "internal error"
                    raise RuntimeError, "unrecoverable error"
        return result
    def getQArray(self):
        retval = num.array(map(lambda x:x.q, self.qSymm))
        return retval 
    def checkForBinary(self, qList):
        binaryRelatedTo = zeros((len(qList),len(qList)), dtype=bool)
        for iQ, q_i in enumerate(qList):
            for jQ in range(iQ+1,len(qList)):
                q_j = qList[jQ]
                for k in range(self.nSymm):
                    for l in range(self.nSymm):
                        tempa = q_j * self.qSymm[k]
                        tempb = q_i * self.qSymm[l]
                        misAng = tempb.misorAng(tempa)
                        'misAng is in range [0,2*pi]'
                        if abs(misAng - math.pi) < epsSqrt:
                            binaryRelatedTo[iQ,jQ] = binaryRelatedTo[jQ,iQ] = True
        return binaryRelatedTo

def transposeQuats(qList):
    for qTemp in qList:
        qTemp.transpose()

def stripQuatList(qList):
    for qTemp in qList:
        qTemp.strip()

def writeQuats(qList, f):
    for q in qList:
        print >> f, "%s %s %s %s" % (q.q[0],q.q[1],q.q[2],q.q[3])
#
# from mmod_set_texture_evppg:
#    CALL mat_x_mat_3(state%c, sc_rel%B_matx(:,:,i_cnst), c_cnst) 
#    CALL tensor_to_quat(texture(:,i_cnst), c_cnst)
# texture(:,i_cnst) then has texture quat for the constituent
#
# from setup_pt_wsymm:
#        CALL mat_x_mat_3(&
#             & symm%rsymm(:,:,i_symm), &
#             & lat_OR, &
#             & sc_rel%B_matx(:,:,i_cnst))
# lat_OR one of PT_lat_rel
# symm is for reference crystal, i_symm is in a subset of total symm

def makeQuatsBall(qRef, thetaScale, n):
    """make n quaternions in a ball around qRef, with ball size scaled by thetaScale
    """
    result = []
    for i in range(n):
        qRand = Quat(RotInv('rand',thetaScale))
        qCur = qRand * qRef
        result.append(qCur)
    return result

didSeedSet = False
def makeQuatsComponents(nGrain, scale=None):
    global didSeedSet
    if not didSeedSet:
        num.random.seed(0) # for reproducable results
        didSeedSet = True
    quatsBall = num.ones([nGrain,4])
    if scale is not None:
        quatsBall[:,1:] = num.random.normal(loc=0.0, scale=scale, size=[nGrain,3])
    else:
        'ignore set of measure zero over which mag will be zero'
        quatsBall[:,:] = num.random.normal(loc=0.0, scale=1.0, size=[nGrain,4])
    for iQ in range(quatsBall.shape[0]):
        mag = math.sqrt(num.sum(quatsBall[iQ,:]*quatsBall[iQ,:]))
        quatsBall[iQ,:] = quatsBall[iQ,:] / mag
    return quatsBall

def millerBravais2Normal(invec, *args):
    """
    Generate the normal(s) for a plane(s) given in 
       the Miller-Bravais convention for the
       hexagonal basis {a1, a2, a3, c}.  The 
       basis for the output {o1, o2, o3}
       is chosen such that:

       o1 || a1
       o3 || c
       o2 = o3 ^ o1

    returns a (3, n) array of horizontally concatenated
    unit vectors
    """
    aspect = 1.
    if len(args) > 0:
        aspect = args[0]

    r1 = invec[0]
    r2 = (2*invec[1] + invec[0]) / sqrt(3.)
    r3 = invec[3]/aspect
    
    outvec = num.array([r1,r2,r3])
    outvec = outvec / math.sqrt(num.dot(outvec,outvec))

    return outvec

class Fiber:
    """
    Like John Edmiston's MakeFiber class, but with the implementation more tightly coupled to the rest of the code base
    """
    def __init__(self, latVec, recipVec):
        
        self.latVec   = mU.normalized(latVec)
        self.recipVec = mU.normalized(recipVec)

        self.qBase  = q1  = Quat(RotInv('align', 
                                                self.latVec,
                                                self.recipVec)
                                     )
        e1 = q1.q
        rFiber = RotInv(0.5, self.recipVec) # 0.5 is a very arbitary distance along the fiber
        q2 = Quat(rFiber)
        q12 = (q2*q1).q
        e2 = mU.normalized(q12 - num.dot(q12,e1) * e1)
        #
        self.e1 = e1
        self.e2 = e2
        
        return
    def __sub__(self,other):
        return self.distBetweenFibers(other)        
    def distBetweenFibers(self,other):
        """
Compute the distance between two fibers using the polar decomposition of the projection operator taking one geodesic plane to the other. 
input: instance of MakeFiber class
output: (max_eigenvalue, Rotation at max_eigenvalue), intersecting fibers would have max_eigenvalue ~ 1. Rotation at max_eigenvalue would be the 'closest' Rotation which would relate the two fibers.
        """
        
        e1,e2 = self.e1,self.e2
        e1_,e2_ = other.e1,other.e2
        F11 = num.dot(e1,e1_);F12 = num.dot(e1,e2_);
        F21 = num.dot(e2,e1_);F22 = num.dot(e2,e2_);
        #not guaranteed to be non singular... most of the time seems to work though
        F = num.array([[F11,F12],[F21,F22]])
        #want to use F = R.U for F non singular
        #eigenvalues of U give measure of fiber distance
        #R*eigenvector(U) gives quaternion at min/max fiber distance on e1,e2 (self's) basis
        C = num.dot(F.T,F)
        eigval,eigvec = num.linalg.eig(C)
        l1,l2 = map(float,eigval);
        u1,u2 = eigvec[:,0],eigvec[:,1]        
        if(l1>l2):
            max_eigval = l1
            max_eigvec = u1
        else:
            max_eigval = l2
            max_eigvec = u2
        #had non orthogonal u1,u2 once with C++ implementation using lapack::geev which is probably what scipy uses, so keep this "if ..." safeguard unless there is a better safeguard against finding non orthogonal eigenvectors
        if num.dot(u1,u2)>1.0e-4:
            tmp = num.zeros(3)
            tmp[0] = u1[0]
            tmp[1] = u1[1]
            tmp[2] = 0.
            tmp2 = num.cross(num.array([0.,0.,1.]),tmp)
            u2[0] = tmp2[0]
            u2[1] = tmp2[1]
        #polar composition F = R.U, F^TF = U^2 = C, R = F.U^-1
        #U = sqrt(abs(l1))*num.outer(u1,u1)+sqrt(abs(l2))*num.outer(u2,u2)
        U_inv = (1./sqrt(abs(l1)))*num.outer(u1,u1) + (1./sqrt(abs(l2)))*num.outer(u2,u2)
        R = num.dot(F,U_inv)
        Ru1 = num.dot(R,max_eigvec)
        q_min_dist = Ru1[0]*e1 + Ru1[1]*e2
        #check for consistency
        q1_ = num.dot(q_min_dist,e1_)
        q2_ = num.dot(q_min_dist,e2_)
        #max_eigval = max(l1,l2)
        #just for initial testing of the fiber method, trying to understand the properties
        if(abs(sqrt(q1_*q1_+q2_*q2_) - mU.normvec(q_min_dist))>1.0e-3 and max_eigval>=.9999):
            raise Exception, 'min dist quaternion is not in both planes'
           
        R_min_dist = Quat(q_min_dist) # da.quaternion_map(q_min_dist)
        return max_eigval,R_min_dist

    def constructOrientation(self, angle):
        qAxis = Quat(ors.RotInv(angle, self.latVec))
        qCur  = self.qBase * qAxis
        return qCur.q



    
