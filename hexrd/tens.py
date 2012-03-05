#! /usr/bin/env python
# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC. 
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
# stuff for representations of tensor components: tens.py
import sys
from math import sqrt

from numpy import *
import numpy as num

from hexrd.constants import *

def vecdvToVecds(vecdv):
    """convert from [t1,...,t5,v] to vecds[:] representation,
    where v is the relative volume"""
    vecds = zeros(6,dtype='float64')
    vecds[:-1] = vecdv[:-1]
    vecds[-1] = log(vecdv[-1])
    return vecds

def vecdsToSymm(vecds):
    """convert from vecds representation to symmetry matrix"""
    A = zeros([3,3],dtype='float64')
    Akk_by_3 = sqr3i * vecds[5] # -p
    t1 = sqr2i*vecds[0]
    t2 = sqr6i*vecds[1]

    A[0, 0] =    t1 - t2 + Akk_by_3     
    A[1, 1] =   -t1 - t2 + Akk_by_3     
    A[2, 2] = sqr2b3*vecds[1] + Akk_by_3
    A[1, 0] = vecds[2] * sqr2i          
    A[2, 0] = vecds[3] * sqr2i          
    A[2, 1] = vecds[4] * sqr2i          

    A[0, 1] = A[1, 0]
    A[0, 2] = A[2, 0]
    A[1, 2] = A[2, 1]
    return A

def traceToVecdsS(Akk):
    return sqr3i * Akk

def vecdsSToTrace(vecdsS):
    return vecdsS * sqr3

def trace3(A):
    return A[0,0]+A[1,1]+A[2,2]

def symmToVecds(A):
    """convert from symmetry matrix to vecds representation"""
    vecds = zeros(6,dtype='float64')
    vecds[0] = sqr2i * (A[0,0] - A[1,1])
    vecds[1] = sqr6i * (2. * A[2,2] - A[0,0] - A[1,1]) 
    vecds[2] = sqr2 * A[1,0]
    vecds[3] = sqr2 * A[2,0]
    vecds[4] = sqr2 * A[2,1]
    vecds[5] = traceToVecdsS(trace3(A))
    return vecds

def matxToSkew(A):
    return skewOfMatx(A)
def skewOfMatx(A):
    W = zeros([3,3],dtype='float64')
    W[0,1] =  0.5*(A[0,1] - A[1,0])
    W[0,2] =  0.5*(A[0,2] - A[2,0])
    W[1,2] =  0.5*(A[1,2] - A[2,1])
    W[1,0] = -W[0,1]
    W[2,0] = -W[0,2]
    W[2,1] = -W[1,2]
    return W

def matxToSymm(A):
    return symmOfMatx(A)
def symmOfMatx(A):
    D = zeros([3,3],dtype='float64')
    D = 0.5*(A + A.T)
    return D

def symmToMVvec(A):
    """
    convert from symmetric matrix to Mandel-Voigt vector
    representation (JVB)
    """ 
    mvvec = zeros(6, dtype='float64')
    mvvec[0] = A[0,0]
    mvvec[1] = A[1,1]
    mvvec[2] = A[2,2]
    mvvec[3] = sqr2 * A[1,2]
    mvvec[4] = sqr2 * A[0,2]
    mvvec[5] = sqr2 * A[0,1] 
    return mvvec

def MVvecToSymm(A):
    """
    convert from Mandel-Voigt vector to symmetric matrix 
    representation (JVB) 
    """ 
    symm = zeros((3, 3), dtype='float64')
    symm[0, 0] = A[0]
    symm[1, 1] = A[1]
    symm[2, 2] = A[2]
    symm[1, 2] = A[3]/sqr2
    symm[0, 2] = A[4]/sqr2
    symm[0, 1] = A[5]/sqr2 
    symm[2, 1] = A[3]/sqr2
    symm[2, 0] = A[4]/sqr2
    symm[1, 0] = A[5]/sqr2 
    return symm

def MVCOBMatrix(R):
    """
    GenerateS array of 6 x 6 basis transformation matrices for the
    Mandel-Voigt tensor representation in 3-D given by: 
    
    [A] = [[A_11, A_12, A_13],
           [A_12, A_22, A_23],
           [A_13, A_23, A_33]]
        |
        |
        V
    {A} = [A_11, A_22, A_33, sqrt(2)*A_23, sqrt(2)*A_13, sqrt(2)*A_12]
              
    where the operation R * A *R.T (in tensor notation) is obtained by
    the matrix-vector product [T]*{A}.
    
    USAGE
    
        T = MVCOBMatrix(R)
    
    INPUTS
    
        1) R is (3, 3) an ndarray representing a change of basis matrix
    
    OUTPUTS
    
        1) T is (6, 6), an ndarray of transformation matrices as
           described above
    
    NOTES
    
        1) Compoments of symmetric 4th-rank tensors transform in a
           manner analogous to symmetric 2nd-rank tensors in full
           matrix notation. 

    SEE ALSO

    symmToMVvec, quatToMat
    """
    T = zeros((6, 6), dtype='float64')
    
    T[0, 0] = R[0, 0]**2;    
    T[0, 1] = R[0, 1]**2;
    T[0, 2] = R[0, 2]**2;
    T[0, 3] = sqr2 * R[0, 1] * R[0, 2];
    T[0, 4] = sqr2 * R[0, 0] * R[0, 2];
    T[0, 5] = sqr2 * R[0, 0] * R[0, 1];
    T[1, 0] = R[1, 0]**2;
    T[1, 1] = R[1, 1]**2;
    T[1, 2] = R[1, 2]**2;
    T[1, 3] = sqr2 * R[1, 1] * R[1, 2];
    T[1, 4] = sqr2 * R[1, 0] * R[1, 2];
    T[1, 5] = sqr2 * R[1, 0] * R[1, 1];
    T[2, 0] = R[2, 0]**2;
    T[2, 1] = R[2, 1]**2;
    T[2, 2] = R[2, 2]**2;
    T[2, 3] = sqr2 * R[2, 1] * R[2, 2];
    T[2, 4] = sqr2 * R[2, 0] * R[2, 2];
    T[2, 5] = sqr2 * R[2, 0] * R[2, 1];
    T[3, 0] = sqr2 * R[1, 0] * R[2, 0];
    T[3, 1] = sqr2 * R[1, 1] * R[2, 1];
    T[3, 2] = sqr2 * R[1, 2] * R[2, 2];
    T[3, 3] = R[1, 2] * R[2, 1] + R[1, 1] * R[2, 2];
    T[3, 4] = R[1, 2] * R[2, 0] + R[1, 0] * R[2, 2];
    T[3, 5] = R[1, 1] * R[2, 0] + R[1, 0] * R[2, 1];
    T[4, 0] = sqr2 * R[0, 0] * R[2, 0];
    T[4, 1] = sqr2 * R[0, 1] * R[2, 1];
    T[4, 2] = sqr2 * R[0, 2] * R[2, 2];
    T[4, 3] = R[0, 2] * R[2, 1] + R[0, 1] * R[2, 2];
    T[4, 4] = R[0, 2] * R[2, 0] + R[0, 0] * R[2, 2];
    T[4, 5] = R[0, 1] * R[2, 0] + R[0, 0] * R[2, 1];
    T[5, 0] = sqr2 * R[0, 0] * R[1, 0];
    T[5, 1] = sqr2 * R[0, 1] * R[1, 1];
    T[5, 2] = sqr2 * R[0, 2] * R[1, 2];
    T[5, 3] = R[0, 2] * R[1, 1] + R[0, 1] * R[1, 2];
    T[5, 4] = R[0, 0] * R[1, 2] + R[0, 2] * R[1, 0];
    T[5, 5] = R[0, 1] * R[1, 0] + R[0, 0] * R[1, 1];
    return T

def NormalProjectionOfMV(vec):
    # 
    # To perform n' * A * n as [N]*{A}
    #

    # normalize in place... col vectors!
    v2 = vec**2
    n  = vec / sqrt(tile(v2.sum(0), (vec.shape[0], 1)))
    
    nmat = array([
        n[0, :]**2, 
        n[1, :]**2, 
        n[2, :]**2, 
        sqr2*n[1, :]*n[2, :], 
        sqr2*n[0, :]*n[2, :], 
        sqr2*n[0, :]*n[1, :]])
    
    nmat = nmat.T
    return nmat

def svecToVecds(svec):
    """convert from svec to vecds representation"""
    vecds = zeros(6,dtype='float64')
    vecds[0] = sqr2i * (svec[0] - svec[1])
    vecds[1] = sqr6i * (2. * svec[2] - svec[0] - svec[1]) 
    vecds[2] = sqr2 * svec[5]
    vecds[3] = sqr2 * svec[4]
    vecds[4] = sqr2 * svec[3]
    vecds[5] = traceToVecdsS(svec[0]+svec[1]+svec[2])
    return vecds

def symmPlusI(Ain):
    """add the identity to a symmetric matrix"""
    A = zeros([3,3],dtype='float64')
    A[:,:] = Ain[:,:]
    A[0,0] = A[0,0] + 1.
    A[1,1] = A[1,1] + 1.
    A[2,2] = A[2,2] + 1.
    return A

def svecpToSvec(svecp):
    svec = zeros(6,dtype='float64')
    svec[0:6] = svecp[0:6]
    p         = svecp[6]
    svec[0] = svec[0] - p
    svec[1] = svec[1] - p
    svec[2] = svec[2] - p
    return svec

def symmToSvec(symm):
    svec = zeros(6,dtype='float64')
    svec[0] = symm[0,0]
    svec[1] = symm[1,1]
    svec[2] = symm[2,2]
    svec[3] = symm[1,2]
    svec[4] = symm[2,0]
    svec[5] = symm[0,1]
    return svec

def matxToSvec(matx):
    svec = zeros(6,dtype='float64')
    svec[0] = matx[0,0]
    svec[1] = matx[1,1]
    svec[2] = matx[2,2]
    svec[3] = (matx[1,2]+matx[2,1])*0.5
    svec[4] = (matx[2,0]+matx[0,2])*0.5
    svec[5] = (matx[0,1]+matx[1,0])*0.5
    return svec

svecGam = array([ 1.,  1.,  1.,  2.,  2.,  2.])

def svecToMatx(svec):
    return svecToSymm(svec)
def svecToSymm(svec):
    matx = zeros([3,3],dtype='float64')
    matx[0,0] = svec[0]
    matx[1,1] = svec[1]
    matx[2,2] = svec[2]
    matx[1,2] = svec[3]
    matx[2,0] = svec[4]
    matx[0,1] = svec[5]
    matx[2,1] = svec[3]
    matx[0,2] = svec[4]
    matx[1,0] = svec[5]
    return matx

def dAiAoH_svecList(aInv):
    """
    derivative of inverse of symmetric matrix wrt svec components of that matrix; aInv is the inverse of the matrix
    """
    retval = []
    for iSvec in range(6):
        svec = num.zeros(6)
        svec[iSvec] = 1.
        hMatx = svecToSymm(svec)
        retval.append(
            - num.dot( num.dot(aInv, hMatx), aInv )
              )
    return retval

def svecToSvecP(svec):
    svecp = zeros(7,dtype='float64')
    p          = -(svec[0]+svec[1]+svec[2])/3.0
    svecp[6]   = p
    svecp[0:3] = array(svec[0:3]) + p
    svecp[3:6] = svec[3:6]
    return svecp


class T2Symm:
    'template for symmetric second order tensor components'
    def __init__(self, args):
        raise RuntimeError, "need implementation"
    def toSymm(self):
        raise RuntimeError, "need implementation"

class T2Vecds(T2Symm):
    def __init__(self, vecds):
        self.vecds = num.array(vecds)
        return
    def toSymm(self):
        return vecdsToSymm(self.vecds)
    def toVecds(self):
        return self # self.vecds
    def toUnit(self):
        mag = sqrt(num.sum(self.vecds[:] * self.vecds[:]))
        return T2Vecds(self.vecds / mag)

class T2Svec(T2Symm):
    def __init__(self, val):
        # self.svec = num.array(svec)
        if hasattr(val,'toSymm'):
            self.svec = symmToSvec(val.toSymm())
        elif hasattr(val, '__len__'):
            if len(val) == 7:
                self.svec = svecpToSvec(val)
            elif len(val) == 6:
                self.svec = num.array(val)
            elif len(val.shape) == 2:
                self.svec = matxToSvec(val)
            else:
                print >> sys.stderr, "bad size "+str(val)
                raise RuntimeError, "unrecoverable error"
        else:
            raise RuntimeError, 'do not know what to do with '+str(val)
        return
    def toSymm(self):
        return svecToMatx(self.svec)
    def toVecds(self):
        return T2Vecds(svecToVecds(self.svec))
    def toSvec(self):
        return self.svec

# class TCompSvec:
#     def __init__(self, val):
#         if len(val) == 7:
#             self.c = svecpToSvec(val)
#         elif len(val) == 6:
#             self.c = array(val)
#         elif len(val.shape) == 2:
#             self.c = matxToSvec(val)
#         else:
#             print >> sys.stderr, "bad size"
#             raise RuntimeError, "unrecoverable error"
#         return
#     def __call__(self):
#         return self.c
class T2SvecP(T2Symm):
    def __init__(self, svecp):
        self.svecp = num.array(svecp)
        return
    def toSymm(self):
        return svecToMatx(svecpToSvec(self.svecp))
    def toVecds(self):
        return T2Vecds(svecToVecds(svecpToSvec(self.svecp)))

sProj = {}
sProj['sxxDev'] = T2Svec(num.array([ 2.0,-1.0,-1.0,0.,0.,0.])/3.0)
sProj['syyDev'] = T2Svec(num.array([-1.0, 2.0,-1.0,0.,0.,0.])/3.0)
sProj['szzDev'] = T2Svec(num.array([-1.0,-1.0, 2.0,0.,0.,0.])/3.0)
sProj['sxx']    = T2Svec(num.array([ 1.0, 0.0, 0.0,0.,0.,0.]))
sProj['syy']    = T2Svec(num.array([ 0.0, 1.0, 0.0,0.,0.,0.]))
sProj['szz']    = T2Svec(num.array([ 0.0, 0.0, 1.0,0.,0.,0.]))
sProj['p']      = T2Svec(num.array([-1.0,-1.0,-1.0,0.,0.,0.])/3.0)
#
sProj['sxxDevUnit'] = sProj['sxxDev'].toVecds().toUnit()
sProj['syyDevUnit'] = sProj['syyDev'].toVecds().toUnit()
sProj['szzDevUnit'] = sProj['szzDev'].toVecds().toUnit()

