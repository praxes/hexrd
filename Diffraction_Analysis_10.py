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

import numpy
from Vector_funcs import *
from numpy import dot
from scipy.linalg import det,inv,eig
from scipy.optimize import fmin
import scipy.optimize
from math import cos,sin,atan,sqrt,tan,asin,acos
import math
from data_class import inv_dict
import random
from Misc_funcs import *
def deg(x):
    return x*180./math.pi
def rad(x):
    return x*math.pi/180.
def Construct_s_a_b(a,b):
    return sin(b)*e3func() + cos(b)*erfunc(a)
def Get_a_b(s):
    pass
def Find_Rotation_Angle(p,R,w0 = 0.1,tol = 10**-1,stop = 0):
    
        
    print 'rot_angle'
    def diffa(w,p,R):
        #print 'w,p,R',w,p,R
        d = Rodrigues3(p,w)-R
        return two_norm(d,d)
    #out = scipy.optimize.fmin(diffa,x0 = (w0,),args = (p,R), full_output = 1)
    
    #if out[1]>1:
    #    if stop==1:
    #        raise Exception
    #    out = Find_Rotation_Angle(p,R,w0 = w0-.57,stop = 1)
    #    return out
        
    out = scipy.optimize.fsolve(diffa,x0 = (w0,),args = (p,R), full_output = 1)
    return out

def Find_Rotation_Comps(R,call = 0,tol=10**-2):
    """finds the rotation components w,(r,theta,phi) making up a matrix, for use with Rodrigues3(Unit(erhofunc(theta,phi)),w)""" 
    #if call==10:
        #raise Exception,"unable to find correct oslution"
    eval,evec = eig(R)
    index = list(eval).index(max(list(eval)))
    axis = numpy.asarray(evec[:,index],dtype = 'float64')
    axis = Unit(axis)
    arb_vec = numpy.array([random.random(),random.random(),random.random()],dtype = 'float64')
    arb_vec = Unit(arb_vec)
    while Mag(arb_vec-axis)<tol:
        arb_vec = numpy.array([random.random(),random.random(),random.random()],dtype = 'float64')
    proj_arb_vec = Unit(dot(Proj_func(axis),arb_vec))
    arb_vec_prime = Unit(dot(R,proj_arb_vec))
    
    w = Compute_Angle(proj_arb_vec,arb_vec_prime)
    sign = numpy.sign(dot(axis,Cross(proj_arb_vec,arb_vec_prime)))
    if sign!=0:
        w = w*sign
    r,theta,phi = Get_Spherical_Coords(axis)
    
        
    
    if two_norm_diff(Rodrigues3(axis,w),R)>tol:
        print "two_norm",R,w,theta,phi,proj_arb_vec,arb_vec_prime
        #call +=1
        #w,(r,theta,phi) = Find_Rotation_Comps(R,call)
        
        raise AttributeError, "Rotations do not match, error"
    return w,(r,theta,phi)

    
def Find_Rotation(v1,v2,v3,v4,v5,v6,w0=0,tol = 10**-1):
    #Rv1 = v2
    print "find_rotation"
    def Rdiff(rparams,v1,v2,v3,v4,v5,v6):
        w,theta,phi = rparams
        R1 = Rodrigues3(e1func(),w)
        R2 = Rodrigues3(e2func(),theta)
        R3 = Rodrigues3(e3func(),phi)
        R = dot(R1,R2)
        R = dot(R,R3)
        #p =erhofunc(theta,phi)
        #R = Rodrigues3(p,w)
        out1 = (dot(R,v1)-v2)
        out2 = (dot(R,v3)-v4)
        out3 = (dot(R,v5)-v6)
        o1 = list(out1)
        o2 = list(out2)
        o3 = list(out3)
        o1.extend(o2)
        o1.extend(o3)
        
        o1 = numpy.array(o1)
        return Mag(o1)
    res = scipy.optimize.fmin(Rdiff,x0 = (w0,0,0),args = (v1,v2,v3,v4,v5,v6),full_output = 1,ftol = tol,xtol = 10**-5,maxfun = 100000)
    #res = scipy.optimize.fsolve(Rdiff,x0 = (w0,0,0),args = (v1,v2,v3,v4,v5,v6),full_output = 1,xtol = 10**-5)
    w,theta,phi = res[0]
    R1 = Rodrigues3(e1func(),w)
    R2 = Rodrigues3(e2func(),theta)
    R3 = Rodrigues3(e3func(),phi)
    R = dot(R1,R2)
    R = dot(R,R3)
    
    return res,R

    

def Compute_Recip_Vector(w,eta,theta,lam = 1.):
    R = Rodrigues(w,e3func(),e1func(),e2func())
    s_s0 = Construct_s(eta,theta)-e1func()
    return (1./lam)*dot(inv(R),s_s0)
def Quadratic(A,B,C):
    return (-B + sqrt(B**2 - 4*A*C))/(2*A), (-B - sqrt(B**2 - 4*A*C))/(2*A)

def Get_Recip_Vector_From_Pixel(pixel,w,center_px,L0,mm_per_px = .2,lam = .15338 * 10**-10):
    """0,0 in upper left corner, -x to the right, +y up"""

    loc = pixel-center_px
    print 'loc',loc
    e1 = numpy.array([1,0])
    e2 = numpy.array([0,1])
    Y = dot(loc,e2)*mm_per_px
    X = dot(loc,-e1)*mm_per_px
    print 'X,Y',X,Y
    r = L0*e1func()+X*e2func()+Y*e3func()
    s = Unit(r)
    print 's',s
    r,b,a = Get_Spherical_Coords(s)
    b = math.pi/2-b
    print 'r,b,a',r,b,a
    print 'eta,theta',alphabeta_to_etatheta(a,b)
    s_s0 = s-e1func()
    R = Rodrigues3(e3func(),w)
    gi = dot(R.T,s_s0)/lam
    return gi,(a,b)

def Get_Recip_Vector(w,a,b,wave_length = 1.):
    """compute recip vector from R(w)^-1 (s-s_0/lambda)(a,b)"""
    s = sin(b)*e3func()+cos(b)*erfunc(a)
    s0 = e1func()
    R = Rodrigues(w,e3func(),e1func(),e2func())
    R_inv = inv(R)
    g_I = (1./wave_length)*dot(R_inv,s-s0)
    return g_I

def find_w_a_b(v0,lam,w0 = 0,tol = 10**-1):
    #w  = 0 at e1
    def get_diff(arg1,arg2,arg3):
        #print 'arg1,arg2,arg3,arg4',arg1,arg2,arg3
        lam = arg3
        v = arg2
        w,a,b = arg1
        #print w,a,b,v
        g_I = Get_Recip_Vector(w,a,b,lam)
        return (dot(g_I-v,g_I-v))**.5
    def get_diff2(arg1,arg2,arg3):
        
        #tol = 10**-3
        print 'arg1,arg2,arg3',arg1,arg2,arg3
        lam = arg3
        v = arg2
        w,a,b = arg1
        #print w,a,b,v
        g_I = Get_Recip_Vector(w,a,b,lam)
        print w,'gI-v',g_I - v
        return (g_I-v)*lam
    #print 'fmin'
    res = scipy.optimize.fmin(get_diff,x0 = (w0,0,0),args = (v0,lam),full_output = 1,ftol = tol,xtol = 10**-5,maxfun = 100000)
    #res2 = scipy.optimize.fsolve(get_diff2,x0 = (0,0,0),args = (v0,lam),full_output = 1,xtol = 1**-12)
    return res[0]
def alphabeta_to_etatheta(a,b):
    theta = asin(sqrt((1 - cos(a)*cos(b))/2.))
    eta1 = asin(dot(Construct_s_a_b(a,b) - e1func(),e3func())/sin(2*theta))
    eta2 = numpy.sign(eta1)*(abs(eta1)+2*(math.pi/2. - abs(eta1)))
    #s1 = Construct_s(eta1,theta)
    #s2 = Construct_s(eta2,theta)
    s = Construct_s_a_b(a,b)
    if dot(s,e2func())>0:
        eta = eta2
    else:
        eta = eta1
    return eta,theta
            
def CCW(x_in_rad):
    if x_in_rad > 2*math.pi:
        revs = Int_Map(x_in_rad/(2*math.pi),rounding = 'ROUND_DOWN')
        return CCW(x_in_rad - revs*2*math.pi)
    if x_in_rad<0:
        return 2*math.pi - abs(x_in_rad)
    else:
        return x_in_rad
def find_pixel_orientation_polar(v,lam=1.):
    g1,g2,g3 = v
    d = Mag(v)**-1
    theta = asin(lam/(2*d))
    beta = acos(lam*g3)
    alpha1 = acos(cos(2*theta)/sin(beta))
    alpha2 = -alpha1
    alphas = [alpha1,alpha2]
    A = -g2
    B = g1
    if A>=0:
        phi = asin(B/sqrt(A**2+B**2))
    else:
        phi = math.pi - asin(B/sqrt(A**2+B**2))
    
    C = (cos(2*theta) - 1)/lam
    w1 = asin(C/(sqrt(A**2+B**2)))-phi
    #what quadrant is R.v
    w_max = math.pi/2. - phi
    w2 = w1 + 2*abs((w_max - w1))
    Rv = dot(Rodrigues3(e3func(),w1),v)
    quadrant = Quadrant2D(Rv)
    w1 = CCW(w1)
    w2 = CCW(w2)
    if quadrant == 2:
        pair1 = [w1,alpha1,beta]
        pair2 = [w2,alpha2,beta]
    elif quadrant == 3:
        #alpha should be negative here
        pair1 = [w1,alpha2,beta]
        pair2 = [w2,alpha1,beta]
    return pair1,pair2
    
        
    
def Construct_s_s0_2_(a,b):
    return e3func()*cos(b)+sin(b)*(cos(a)*e1func()+sin(a)*e2func())-e1func()

def Reverse_vec(w,a,b,lam):
    R = Rodrigues3(e3func(),w)
    s_s0 = Construct_s_s0_2_(a,b)
    return dot(R.T,s_s0)/lam

def Compute_Coords_From_Recip_Vector(v,lam=1.):
    
    d = sqrt(dot(v,v))**-1
    theta = asin(lam/(2*d)) #Should always work
    v1 = dot(v,e1func())
    v2 = dot(v,e2func())
    v3 = dot(v,e3func())
    a_ = -v2
    b_ = v1
    c_ = (1./lam)*(cos(2*theta)-1)
    c = 2*sin(theta)**2
    a = v2
    b = -v1    
    if a_<0:
        phi = (math.pi - asin((b_)/sqrt(a_**2+b_**2)))
        w = asin((c_)/sqrt(a_**2+b_**2)) - phi       
    else:
        phi = asin((b_)/sqrt(a_**2+b_**2))
        w = asin((c_)/sqrt(a_**2+b_**2)) - phi
    w1 = w
    wmax = math.pi/2. - phi
    w2 = w1 + 2*abs(w1-wmax)
    
    beta = asin(lam*v3) #should be ok with single value
    alpha1 = acos((2*sin(theta)**2 - 1 )/-cos(beta))
    alpha2 = -acos((2*sin(theta)**2 - 1 )/-cos(beta))
    s_s0_1 = Construct_s_a_b(alpha1,beta) - e1func()
    s_s0_2 = Construct_s_a_b(alpha2,beta) - e1func()
    R1 = Rodrigues(w1)
    R2 = Rodrigues(w2)
    ws = [w1,w2]
    alphas = [alpha1,alpha2]
    betas = [beta]
    sols = inv_dict()
    tol = 1
    ct = 0
    for i in range(len(ws)):
        for j in range(len(alphas)):
            for k in range(len(betas)):
                w,alpha,beta = ws[i],alphas[j],betas[k]
                diff = Mag(Get_Recip_Vector(w,alpha,beta,lam)-v)
                if diff < tol:
                    sols.Add(ct,(w,alpha,beta))
                    ct+=1
           
    try:

        if len(sols.keys())!=2:
            raise Exception, "only one solution" +w1.__str__()+'_'+w2.__str__()
        return sols
    except Exception:
        print "only 1 solution"
        return sols

def Construct_s(eta,theta):
    s = cos(2*theta)*e1func()+sin(2*theta)*(erfunc(eta,-e2func(),e3func()))
    return s
def Compute_Pixel_From_Recip_Vector(p0,v,L_0,center_px,mm_per_px = .2):
    raise Exception, 'not fixed'
    w,eta,theta = Compute_Coords_From_Recip_Vector(v)
    return w,eta,theta,Compute_Pixel_From_Coords(p0,w,eta,theta,L_0,center_px,mm_per_px = .2)
    
def keV2Angstrom(energy):
    h = 4.13566743
    c = 2.99792458
    cfac = h*c
    return cfac/energy
def Compute_Pixel_From_Coords_a_b(p0,w,a,b,L_0,center_px,mm_per_px = .2):
    R = Rodrigues(w,e3func(),e1func(),e2func())
    p = dot(R,p0)
    s = Construct_s_a_b(a,b)
    L = L_0 - dot(p,e1func())
    mag_r = L/(dot(s,e1func()))
    r = mag_r*s
    R_0 = r+p
    W_0 = dot(R_0,e2func())
    H_0 = dot(R_0,e3func())
    pix_coord = numpy.array([-H_0,-W_0])/mm_per_px
    #pix_coord = numpy.array([-W_0,-H_0])/mm_per_px
    return center_px + pix_coord
    
    
def Compute_Pixel_From_Coords(p0,w,eta,theta,L_0,center_px,mm_per_px = .2):
    R = Rodrigues(w,e3func(),e1func(),e2func())
    p = dot(R,p0)
    s = Construct_s(eta,theta)
    L = L_0 - dot(p,e1func()) #L = (R_0 - p).e1
    rho = L*tan(2*theta)
    r = sqrt(L**2 + rho**2)*s #r = |r|*s
    R_0 = r + p
    W_0 = dot(R_0,e2func())
    H_0 = dot(R_0,e3func())
    #pix_coord = numpy.array([-H_0,-W_0])/mm_per_px
    pix_coord = numpy.array([-W_0,-H_0])/mm_per_px
    return center_px + pix_coord
"""    
L=1000.
rho = 100.
lam = 1.
eta = rad(67)
w = rad(40)
theta = rad(10)
v = Compute_Recip_Vector(w,eta,theta,lam = 1.)

res = Compute_Coords_From_Recip_Vector(v)

"""
#center_px = numpy.array([1000,1000])
#px = numpy.array([1000,1100])
#p0=numpy.array([1,0,0])

#pix = Compute_Pixel_From_Recip_Vector(p0,numpy.array([.1,.2,0]),L,center_px)

def quaternion_map((q0,q1,q2,q3),tol = 10**-2):
    w = acos(q0)*2.
    q_ = numpy.array([q1,q2,q3])
    axis = Unit(q_)
    #w_ = asin(Mag(q_))*2.
    #if abs(w -w_)>tol :
    #    raise Exception, 'w does not match'
    return Rodrigues3(axis,w)
def rodrigues_from_quat(quat):
    return rodrigues_from_rot(quaternion_map(quat))
def get_spherical_coords_quat(qin):
    try:
        a = acos(dot(qin,e4func(4)))
    except ValueError:
        if abs(dot(qin,e4func(4)))<tol:
            a = acos(1*numpy.sign(dot(qin,e4func(4))))
    proj = Ifunc4() - numpy.outer(e4func(4),e4func(4))
    v_proj = dot(proj,qin)/sin(a)
    x,b,c = Get_Spherical_Coords(v_proj[0:3])
    return a,b,c
def construct_quat_angle(a,b,c):
    return cos(a)*e4func(4)+sin(a)*(cos(b)*e3func(4)+sin(b)*(e1func(4)*cos(c)+e2func(4)*sin(c)))


    
def quat_from_rot(rot):
    r = r1,r2,r3 = rodrigues_from_rot(rot)
    tan_phi_over_2 = Mag(r)
    n = Unit(r)
    phi = 2*atan(tan_phi_over_2)
    cphi = cos(phi/2.)
    return numpy.array([1,r1,r2,r3])*cphi


def permutation(i,j,k):
    if i==j or j==k or i==k:
        return 0
    elif i<j<k or j<k<i or k<i<j:
        return 1
    else:
        return -1
def rodrigues_from_rot(R):
    trR = numpy.trace(R)
    r0 = -(permutation(0,1,2)*R[1,2]+permutation(0,2,1)*R[2,1])/(1+trR)
    r1 = -(permutation(1,2,0)*R[2,0]+permutation(1,0,2)*R[0,2])/(1+trR)
    r2 = -(permutation(2,0,1)*R[0,1]+permutation(2,1,0)*R[1,0])/(1+trR)
    return numpy.array([r0,r1,r2])


def Unit_Quaternion(inq):
    Mag_q = sqrt(inq[0]**2+inq[1]**2+inq[2]**2+inq[3]**2)
    return inq/Mag_q
#def map_to_fundamental_v(sym_rots,R):
#routs = numpy.dot(R,sym_rots)
    
def map_to_fundamental(sym_rots,R):
    tmp = inv_dict()
    for i in range(len(sym_rots)):
        rot = sym_rots[i]
        Rmap = dot(R,rot.T)
        w,(r,theta,phi) = Find_Rotation_Comps(Rmap)
        tmp.Add(abs(w),(theta,phi,numpy.sign(w)))
    keys = tmp.keys()
    min_w = min(keys)
    if tmp.triggerdict[min_w]==0:
        theta,phi,sgn = tmp[min_w]
    else:
        theta,phi,sgn = tmp[min_w][0]
    axis = erhofunc(theta,phi)
    mag = tan(min_w*sgn/2.)
    r_pos = mag*axis
    x,y,z = r_pos
    return r_pos
def r_to_2theta(r,lam):
    theta = asin(lam*r/2.)
    return 2*theta
def Get_Rodrigues_Curve_From_Vectors(v1,v2):
    a = Unit(v1)
    b = Unit(v2)
    r0 = Cross(a,b)/(1+dot(a,b))
    v = (a+b)/(1+dot(a,b))
    return r0,v
