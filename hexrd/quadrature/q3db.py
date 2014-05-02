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
import numpy as num

from . import q1db

ndim = 3
# formerly:  qloc3ddata

xia = 0.577350269189626e0
xib = 0.25e0
xic = 0.714285714285714285e-1
xid = 0.785714285714285714e0
xie = 0.399403576166799219e0
xif = 0.115964233833200785e0
xig = 0.5e0
xih = 0.0e0
xii = 0.100526765225204467e0
xij = 0.698419704324386603e0
xik = 0.314372873493192195e0
xil = 0.568813795204234229e-1
xim = 0.16666666666666666666e0
xin = 0.333333333333333333e0
xio = 0.909090909090909091e-1
xip = 0.727272727272727273e0
xiq = 0.665501535736642813e-1
xir = 0.433449846426335728e0
xis = 0.7745966692414834e0 # sqrt(0.6e0)
#
wa = -0.131555555555555550e-1
wb = 0.7622222222222222222e-2
wc = 0.2488888888888888880e-1
wd = 0.317460317460317450e-2
we = 0.147649707904967828e-1
wf = 0.221397911142651221e-1
wg = -0.8e0
wh = 0.45e0
wi = 0.602678571428571597e-2
wj = 0.302836780970891856e-1
wk = 0.116452490860289742e-1
wl = 0.109491415613864534e-1


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
