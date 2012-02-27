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
from grainIndex_MakeFiber import MakeFiber
import numpy
from numpy import dot



a1 = numpy.array([4,5,3])
R =  numpy.array([[ 0.95533649, -0.29552021 , 0.        ], [ 0.29552021 , 0.95533649,  0.        ], [ 0.   ,       0.   ,       1.        ]])
a2 = dot(R,a1)
f1 = MakeFiber(a1,a2)

b1 = numpy.array([-6,4,5])
b2 = dot(R,b1)
f2 = MakeFiber(b1,b2)
dist_,Closest_Rotation = f1.distBetweenFibers(f2)
print "original rotation", '\n',R
print "closest_rotation", '\n',Closest_Rotation


R2 = numpy.array([[ 0.69670671, -0.50724736,  0.50724736],       [ 0.50724736,  0.84835335,  0.15164665],       [-0.50724736,  0.15164665,  0.84835335]])

b2_prime = dot(R2,b1)
f2_prime = MakeFiber(b1,b2_prime)
dist_prime, Closest_Rotation_prime = f1.distBetweenFibers(f2_prime)
print "original rotation, different grain", '\n',R2
print "closest_rotation", '\n',Closest_Rotation_prime

#alternate method using rodrigues space
print "alternate method, rodrigues space distance, Closest Rotation"
dist_rodrigues, Closest_Rotation = f1.distBetweenFibers_Rodrigues(f2)
print dist_rodrigues, '\n', Closest_Rotation
