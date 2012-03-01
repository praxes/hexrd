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
from hexrd.orientations import *

kocks = KocksEuler(0., 0., 0.)
quat = Quat(kocks) # quat = kocks.toQuat()

qRand = Quat('rand')
invRand = RotInv(qRand) # qRand.toInv()
kocksRand = KocksEuler(qRand)
qRandB = Quat(kocksRand)
qI = qRand * qRandB.transposed()
qI.toMatrix()
# yes, is the identity to machine precision

# generate orientations perturbed randomly about a given orientation
vl = array([ 0., -1.,  0.])
vs = array([ 1.,  0.,  0.])
rotRef = Quat(RotInv('align',vl,vs))
quatsBall = makeQuatsBall(rotRef, 0.05, 10000)
f = open('inBall.quat','w')
for quat in quatsBall:
    print >> f, quat.__str__() + " 0.5 0.5 0.5"
f.close()











