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
from math import sqrt

from .qloc1ddata import *

def qloc1():
    nqpt = 1
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])

    xi[0]=0.5e0; w[0] = 1.0e0;
    return xi, w

def qloc2():
    nqpt = 2
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])

    xi[0]=(1.0e0 - xi_a2)*0.5e0;  w[0] = 0.5e0;
    xi[1]=(1.0e0 + xi_a2)*0.5e0;  w[1] = 0.5e0;

    return xi, w

def qloc3():
    nqpt = 3
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])
    xi[0]=0.5e0 * (1.0e0 - xi_a3);  w[0] = sqrt(25.0e0 / 324.0e0);
    xi[1]=0.5e0                  ;  w[1] = sqrt(64.0e0 / 324.0e0);
    xi[2]=0.5e0 * (1.0e0 + xi_a3);  w[2] = w[0];
    return xi, w

def qloc4():
    nqpt = 4
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])
    xi_a = sqrt((3.+2.*sqrt(6./5.))/7.)
    xi_b = sqrt((3.-2.*sqrt(6./5.))/7.)
    xi[0]=0.5e0 * (1.0e0 - xi_a);  w[0] = 0.5 * (18.-sqrt(30.))/36.
    xi[1]=0.5e0 * (1.0e0 - xi_b);  w[1] = 0.5 * (18.+sqrt(30.))/36.
    xi[2]=0.5e0 * (1.0e0 + xi_b);  w[2] = w[1]
    xi[3]=0.5e0 * (1.0e0 + xi_a);  w[3] = w[0]
    return xi, w

def qloc5():
    nqpt = 5
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])
    xi_a = sqrt(5.+2.*sqrt(10./7.))/3.
    xi_b = sqrt(5.-2.*sqrt(10./7.))/3.
    xi[0]=0.5e0 * (1.0e0 - xi_a);  w[0] = 0.5 * (322. - 13.*sqrt(70.))/900.
    xi[1]=0.5e0 * (1.0e0 - xi_b);  w[1] = 0.5 * (322. + 13.*sqrt(70.))/900.
    xi[2]=0.5e0 * (1.0e0 + 0.  );  w[2] = 0.5 * 128. / 225.
    xi[3]=0.5e0 * (1.0e0 + xi_b);  w[3] = w[1]
    xi[4]=0.5e0 * (1.0e0 + xi_a);  w[4] = w[0]
    return xi, w

def qloc8():
    nqpt = 8
    xi = num.empty([nqpt]) # ,ndim
    w  = num.empty([nqpt])
    xi[0]=(1.0e0-xi_8a)*0.5e0;  w[0] = w_8a*0.5e0;
    xi[1]=(1.0e0-xi_8b)*0.5e0;  w[1] = w_8b*0.5e0;
    xi[2]=(1.0e0-xi_8c)*0.5e0;  w[2] = w_8c*0.5e0;
    xi[3]=(1.0e0-xi_8d)*0.5e0;  w[3] = w_8d*0.5e0;
    xi[4]=(1.0e0+xi_8d)*0.5e0;  w[4] = w_8d*0.5e0;
    xi[5]=(1.0e0+xi_8c)*0.5e0;  w[5] = w_8c*0.5e0;
    xi[6]=(1.0e0+xi_8b)*0.5e0;  w[6] = w_8b*0.5e0;
    xi[7]=(1.0e0+xi_8a)*0.5e0;  w[7] = w_8a*0.5e0;
    return xi, w

def qLoc(quadr, promote=False):
    if promote:
        nqp = num.array([1,2,3,4,5,8])
        'the following may raise an exception if quadr is too big:'
        quadr = nqp[num.where(nqp >= quadr)][0]
    if quadr == 8:
        xi, w = qloc8()
    elif quadr == 5:
        xi, w = qloc5()
    elif quadr == 4:
        xi, w = qloc4()
    elif quadr == 3:
        xi, w = qloc3()
    elif quadr == 2:
        xi, w = qloc2()
    elif quadr == 1:
        xi, w = qloc1()
    else:
        raise NotImplementedError, 'quadr rule %d not implemented' % (quadr)
    return xi, w
