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

from qloc2dData import *

def qloc1():
    nqpt = 1
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])

    xi[0,0]=xik; xi[0,1]=xik;  w[0] = wf;
    return xi, w

def qloc4():
    nqpt = 4
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])

    xi[0,0]=xil;  xi[0,1]=xil;  w[1] = wj;
    xi[1,0]=xil;  xi[1,1]=xim;  w[2] = wj;
    xi[2,0]=xim;  xi[2,1]=xil;  w[3] = wj;
    xi[3,0]=xim;  xi[3,1]=xim;  w[4] = wj;

    return xi, w
    
def qloc9():
    nqpt = 9
    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])
    xi[0,0]=xii;  xi[0,1]=xii;  w[0] = wg;
    xi[1,0]=xii;  xi[1,1]=xik;  w[1] = wh;
    xi[2,0]=xii;  xi[2,1]=xij;  w[2] = wg;
    xi[3,0]=xik;  xi[3,1]=xii;  w[3] = wh;
    xi[4,0]=xik;  xi[4,1]=xik;  w[4] = wi;
    xi[5,0]=xik;  xi[5,1]=xij;  w[5] = wh;
    xi[6,0]=xij;  xi[6,1]=xii;  w[6] = wg;
    xi[7,0]=xij;  xi[7,1]=xik;  w[7] = wh;
    xi[8,0]=xij;  xi[8,1]=xij;  w[8] = wg;
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
    nqpt = len(w1_i)*len(w1_j)

    xi = num.empty([nqpt,ndim])
    w  = num.empty([nqpt])
    i_qpt = 0
    for xi_i, w_i in zip(xi1_i, w1_i):
        for xi_j, w_j in zip(xi1_j, w1_j):
            xi[i_qpt] = [ xi_i, xi_j ]
            w [i_qpt] =    w_i * w_j 
            i_qpt += 1
    return xi, w
    
def qLoc(quadr):
    if isinstance(quadr,int):
        if quadr == 3:
            xi, w = qloc9()
        elif quadr == 2:
            xi, w = qloc4()
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
