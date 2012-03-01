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
import numpy as num

ndim = 3
from qloc3dData import *

def qloc1():
    nqpt = 1 
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])

    xi[0,:]= 0.5e0; w[0] = 1.0e0
    return xi, w

def qloc8():
    nqpt = 8
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])

    a = (1.0e0 - xia)*0.5e0
    b = (1.0e0 + xia)*0.5e0
    c = 0.5e0*0.5e0*0.5e0

    xi[0,:] = [ a, a, a]
    xi[1,:] = [ a, a, b]
    xi[2,:] = [ a, b, b]
    xi[3,:] = [ a, b, a]
    xi[4,:] = [ b, a, a]
    xi[5,:] = [ b, a, b]
    xi[6,:] = [ b, b, b]
    xi[7,:] = [ b, b, a]

    w[:] = c

    return xi, w
    
def qloc27():
    '3x3x3 quadrature, product of qloc1d03 in three directions'
    """
    import math
    nqpt = 27
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])
    
    xl = num.array(
        [0.5e0 * (1.0e0 - xis), 
         0.5e0,
         0.5e0 * (1.0e0 + xis)
         ])
    wl = num.array([
            math.sqrt(25.0e0 / 324.0e0),
            math.sqrt(64.0e0 / 324.0e0),
            math.sqrt(25.0e0 / 324.0e0)
            ])

    i_qpt = 0
    for i_k in range(3):
        for i_j in range(3):
            for i_i in range(3):
                xi[i_qpt,:] = [ xl[i_i],  xl[i_j],  xl[i_k] ]
                w[i_qpt]    =   wl[i_i] * wl[i_j] * wl[i_k]
                i_qpt += 1
    """
    xi, w = qLocFrom1D(3)
    return xi, w

def qLocFrom1D(quadr1d):
    """
    product of 1d quadrature rules;
    given accuracy may be available with fewer quadrature points
    using a native 3D rule
    """
    from quadrature import q1db
    
    if hasattr(quadr1d,'__len__'):
        assert len(quadr1d) == ndim, 'wrong length'
    else:
        quadr1d = num.tile(quadr1d,(ndim))
    
    xi1_i, w1_i = q1db.qLoc(quadr1d[0], promote=True)
    xi1_j, w1_j = q1db.qLoc(quadr1d[1], promote=True)
    xi1_k, w1_k = q1db.qLoc(quadr1d[2], promote=True)
    nqpt = len(w1_i)*len(w1_j)*len(w1_k)
    
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])
    i_qpt = 0
    for xi_i, w_i in zip(xi1_i, w1_i):
        for xi_j, w_j in zip(xi1_j, w1_j):
            for xi_k, w_k in zip(xi1_k, w1_k):
                xi[i_qpt] = [ xi_i, xi_j, xi_k ]
                w [i_qpt] =    w_i * w_j * w_k
                i_qpt += 1
    return xi, w
    
def qLoc(quadr):
    if isinstance(quadr,int):
        if quadr == 3:
            xi, w = qloc27()
        elif quadr == 2:
            xi, w = qloc8()
        elif quadr == 1:
            xi, w = qloc1()
        else:
            raise NotImplementedError, 'quadr rule %d not implemented' % (quadr)
    else:
        qsplit_x = quadr.split('x')
        qsplit_b = quadr.split('b')
        if len(qsplit_x) == 2:
            assert int(qsplit_x[0]) == ndim, \
                'bad quadr syntax : '+str(quadr)
            quadr1d = int(qsplit_x[1])
            xi, w = qLocFrom1D(quadr1d)
        elif len(qsplit_b) == ndim:
            xi, w = qLocFrom1D(map(int, qsplit_b))
        else:
            raise NotImplementedError, 'quadr rule %s not implemented' % (str(quadr))
    return xi, w
