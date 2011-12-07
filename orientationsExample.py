# DO-NOT-DELETE revisionify.begin() 
#
#   Copyright (c) 2007-2009 Lawrence Livermore National Security,
#   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
#   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
#   
#   Please also read the file NOTICES.
#   
#   This file is part of the mdef package (version 0.2) and is
#   free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   
#   A copy of the GNU Lesser General Public License may be found in the
#   file NOTICES. If this file is missing, see
#   <http://www.gnu.org/licenses/>.
#
# DO-NOT-DELETE revisionify.end() 
from orientations import *

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











