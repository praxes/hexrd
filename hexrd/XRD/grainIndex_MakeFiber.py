# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
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
from numpy import dot,sqrt
from scipy.linalg import eig

import hexrd.Vector_funcs as vf
import hexrd.Diffraction_Analysis_10 as da
from hexrd.Vector_funcs import Rodrigues3,e1func,e2func,e3func,Unit,Mag,\
    From_Rodrigues_Vector,Proj_func
import hexrd.orientations as ors

def dist_to_fiber_alternate(a1,a2,b1,b2):    
    """distance to fiber in Rodrigues space/straight line closest distance"""
    a1_ = Unit(a1); a2_ = Unit(a2)
    r0_1 = num.cross(a1_,a2_)/(1+dot(a1_,a2_))
    g1 = (a1_ + a2_ )/(1 + dot(a1_,a2_))
    g1_ = Unit(g1)
    b1_ = Unit(b1); b2_ = Unit(b2)
    r0_2 = num.cross(b1_,b2_)/(1+dot(b1_,b2_))
    g2 = (b1_ + b2_ )/(1 + dot(b1_,b2_))
    g2_ = Unit(g2)
    Proj = Proj_func(g1_)
    basis_1 = Unit(dot(Proj, g2))
    basis_2 = Unit(num.cross(g1_,basis_1))
    t_intersect = -dot(dot(Proj, r0_2-r0_1), basis_1)
    closest_dist = dot(dot(Proj, r0_2-r0_1), basis_2)
    intersect = r0_2 + t_intersect*g2
    return closest_dist, From_Rodrigues_Vector(intersect)

                                

class MakeFiber:
    """make fiber from two rotation related vectors including fancy basis stuff
    input: vector1, vector2
    attributes: .e1,.e2, the basis vectors in R4 which span the plane containing the fiber traversal geodesic on S3
    
example usage:
from grainIndex import MakeFiber
import num
from num import dot
a1 = num.array([4,5,3])

from Vector_funcs import *
R = Rodrigues3(e3func(),.3) #angle axis

a2 = dot(R,a1)
f1 = MakeFiber(a1,a2)
b1 = num.array([-6,4,5])
b2 = dot(R,b1)
f2 = MakeFiber(b1,b2)
dist_,Closest_Rotation = f1.distBetweenFibers(f2)
#or, equivalently
dist_,Closest_Rotation = f1-f2 #== f2-f1
#dist_ ~= 1 
#Closest_Rotation ~= R

    """
    def __init__(self, latVec, recipVec):

        '''
        this ors stuff in here because it gets used in constructOrientation,
        but ultimately it will probably go away
        '''
        self.latVec   = latVec
        self.recipVec = recipVec
        self.qBase  = q1  = ors.Quat(ors.RotInv('align', Unit(latVec), Unit(recipVec)))
        e1 = q1.q
        Rfiber = Rodrigues3(Unit(recipVec),.5) #.5 is a very arbitary distance along the fiber
        q2 = ors.Quat(Rfiber)
        q12 = (q2*q1).q
        e2 = Unit(q12 - dot(q12,e1)*e1)
        self.e1 = e1
        self.e2 = e2
    def distBetweenFibers_Rodrigues(self,other):
        return dist_to_fiber_alternate(self.latVec,self.recipVec, other.latVec, other.recipVec)
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
        F11 = dot(e1,e1_);F12 = dot(e1,e2_);
        F21 = dot(e2,e1_);F22 = dot(e2,e2_);
        #not guaranteed to be non singular... most of the time seems to work though
        F = num.array([[F11,F12],[F21,F22]])
        #want to use F = R.U for F non singular
        #eigenvalues of U give measure of fiber distance
        #R*eigenvector(U) gives quaternion at min/max fiber distance on e1,e2 (self's) basis
        C = dot(F.T,F)
        eval,evec = eig(C)
        l1,l2 = map(float,eval);
        u1,u2 = evec[:,0],evec[:,1]        
        if(l1>l2):
            max_eval = l1
            max_evec = u1
        else:
            max_eval = l2
            max_evec = u2
        #had non orthogonal u1,u2 once with C++ implementation using lapack::geev which is probably what scipy uses, so keep this "if ..." safeguard unless there is a better safeguard against finding non orthogonal eigenvectors
        if dot(u1,u2)>10**-4:
            tmp = num.zeros(3)
            tmp[0] = u1[0]
            tmp[1] = u1[1]
            tmp[2] = 0.
            tmp2 = num.cross(e3func(),tmp)
            u2[0] = tmp2[0]
            u2[1] = tmp2[1]
        #polar composition F = R.U, F^TF = U^2 = C, R = F.U^-1
        #U = sqrt(abs(l1))*num.outer(u1,u1)+sqrt(abs(l2))*num.outer(u2,u2)
        U_inv = (1./sqrt(abs(l1)))*num.outer(u1,u1) + (1./sqrt(abs(l2)))*num.outer(u2,u2)
        R = dot(F,U_inv)
        Ru1 = dot(R,max_evec)
        q_min_dist = Ru1[0]*e1 + Ru1[1]*e2
        #check for consistency
        q1_ = dot(q_min_dist,e1_)
        q2_ = dot(q_min_dist,e2_)
        #max_eval = max(l1,l2)
        #just for initial testing of the fiber method, trying to understand the properties
        if(abs(sqrt(q1_**2+q2_**2)-Mag(q_min_dist))>10**-3 and max_eval>=.9999):
            raise Exception, 'min dist quaternion is not in both planes'
           
        R_min_dist = da.quaternion_map(q_min_dist)
        return max_eval,R_min_dist

    def constructOrientation(self, angle):
        qAxis = ors.Quat(ors.RotInv(angle, self.latVec))
        qCur  = self.qBase * qAxis
        return qCur.q



