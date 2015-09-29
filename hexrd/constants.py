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
from scipy import constants as sc

# pi related
pi    = np.pi
piby2 = 0.5 * pi
piby3 = pi / 3.
piby4 = 0.25 * pi
piby6 = pi / 6.

# misc radicals
sq2    = np.sqrt(2.)
sq3    = np.sqrt(3.)
sq3by2 = 0.5 * sq3

# tolerancing
epsf      = np.finfo(float).eps      # ~2.2e-16
ten_epsf  = 10 * epsf                # ~2.2e-15
sqrt_epsf = np.sqrt(epsf)            # ~1.5e-8

# for angles
periodDict = {'degrees': 360.0, 'radians': 2*pi}
angularUnits = 'radians' # module-level angle units
d2r = deg2rad = pi / 180.
r2d = rad2deg = 180. / pi

# basis vectors
I3    = np.eye(3) # (3, 3) identity
X = np.array(I3[:, 0].reshape(3, 1), order='C') # X in the lab frame
Y = np.array(I3[:, 1].reshape(3, 1), order='C') # Y in the lab frame
Z = np.array(I3[:, 2].reshape(3, 1), order='C') # Z in the lab frame

zeros_3x1 = np.zeros((3, 1))
zeros_6x1 = np.zeros((6, 1))

# reference beam direction and eta=0 ref in LAB FRAME for standard geometry
bVec = -Z
eVec =  X

# for strain
vInv = np.array(np.r_[1., 1., 1., 0., 0., 0.].reshape(6, 1), order='C')

# for energy/wavelength conversions
keVToAngstrom = lambda x: (1e7*sc.c*sc.h/sc.e) / x
