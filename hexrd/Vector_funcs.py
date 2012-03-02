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
import warnings
from math import cos,sin,acos,asin,atan
import math
import decimal

import numpy
from numpy import dot,sqrt
from scipy.linalg import det,inv,eig

from hexrd.data_class import inv_dict
from hexrd.Misc_funcs import *
from hexrd.Diffraction_Analysis_10 import permutation as perm

def Rod_param(a,b,t=0):
    a_ = Unit(a)
    b_ = Unit(b)
    r0 = Cross(a_,b_)/(1+dot(a_,b_))
    R = From_Rodrigues_Vector(r0)
    r1 = t*(a_+b_)/(1+dot(a_,b_))
    R = From_Rodrigues_Vector(r0+r1)
    return R
    
def Ifunc4():
    return numpy.array([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]],dtype = 'float')
def Ifunc():
    return numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],dtype = 'float')
def SymTensor((a11,a22,a33,a12,a23,a31)):
    return numpy.array([[a11,a12,a31],[a12,a22,a23],[a31,a23,a33]])
def R2Tensor((a11,a12,a13,a21,a22,a23,a31,a32,a33)):
    return numpy.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
def Sym_to_vec(F):
    return numpy.array([F[0,0],F[1,1],F[2,2],F[0,1],F[1,2],F[2,0]])
def R2_to_vec(F):
    return numpy.array([F[0,0],F[0,1],F[0,2],F[1,0],F[1,1],F[1,2],F[2,0],F[2,1],F[2,2]])

def e1func(dim = 3):
    tmp = numpy.zeros(dim)
    tmp[0] = 1
    return tmp

def e2func(dim = 3):
    tmp = numpy.zeros(dim)
    tmp[1]=1
    return tmp


    #return numpy.array([0,1,0])
def e3func(dim = 3):
    tmp = numpy.zeros(dim)
    tmp[2] = 1
    return tmp
    #return numpy.array([0,0,1])

def e4func(dim = 4):
    return numpy.array([0,0,0,1])

def egen_func(index,dim=3):
    tmp = numpy.zeros(dim)
    tmp[index]=1
    return tmp
    

def erfunc(phi,p=e1func(),q=e2func()):
    return p*cos(phi)+q*sin(phi)
def ephifunc(phi):
    return -1*e1func()*sin(phi)+e2func()*cos(phi)
def erhofunc(theta,phi,p = e3func(),q = e1func(),r = e2func()):
    return cos(theta)*p+sin(theta)*erfunc(phi,q,r)
def erhofunc4(a,b,c,s = e4func(4),p = e3func(4),q = e1func(4),r = e2func(4)):
    return cos(a)*s + sin(a)*erhofunc(b,c,p,q,r)
def Get_Metric(g1,g2,g3):
    return numpy.array([[dot(g1,g1),dot(g1,g2),dot(g1,g3)],[dot(g2,g1),dot(g2,g2),dot(g2,g3)],[dot(g3,g1),dot(g3,g2),dot(g3,g3)]])
def Rodrigues(theta,p=e3func(),q=e1func(),r=e2func()):
    p = Unit(p)
    q = Unit(q)
    r = Unit(r)
    rot = Bond(p,p) + cos(theta)*(Bond(q,q)+Bond(r,r))-sin(theta)*(Bond(q,r)-Bond(r,q))
    return rot
def Mag(v1):
    return numpy.sqrt(dot(v1,v1))


def Compute_Angle(v1,v2,tol = 10**-6):
    """compute the angle (in radians) between vector v1 and v2. (numpy 3 element arrays). Tolerance is for arccos protection """
    dotprod = dot(Unit(v1),Unit(v2))
    if abs(abs(dotprod)-1)<tol:
        dotprod = 1*numpy.sign(dotprod)
    if dotprod == 0:
        crossprod = Mag(Cross(Unit(v1),Unit(v2)))
        if abs(abs(crossprod)-1)<tol:
            crossprod = 1*numpy.sign(crossprod)
        angle = asin(crossprod)
    else:
        angle = acos(dotprod)
    return angle
def From_Rodrigues_Vector(v):
    tan_phi_over_2 = Mag(v)
    n = Unit(v)
    phi = 2*atan(tan_phi_over_2)
    return Rodrigues3(n,phi)

def Unit(v1):
    return v1/Mag(v1)
def Bond(v1,v2,e1 = e1func(),e2 = e2func(),e3 = e3func()):
    temp = numpy.zeros([3,3],dtype='float')
    #e1 = numpy.array([1,0,0])
    #e2 = numpy.array([0,1,0])
    #e3 = numpy.array([0,0,1])
    temp[0,0]=dot(v1,e1)*dot(v2,e1)
    temp[0,1]=dot(v1,e1)*dot(v2,e2)
    temp[0,2]=dot(v1,e1)*dot(v2,e3)
    temp[1,0]=dot(v1,e2)*dot(v2,e1)
    temp[1,1]=dot(v1,e2)*dot(v2,e2)
    temp[1,2]=dot(v1,e2)*dot(v2,e3)
    temp[2,0]=dot(v1,e3)*dot(v2,e1)
    temp[2,1]=dot(v1,e3)*dot(v2,e2)
    temp[2,2]=dot(v1,e3)*dot(v2,e3)
    return temp.copy()

class sph_vector:
    def __init__(self,(x,y,z)):
        v = numpy.array([x,y,z])
        self.v = v
        r,theta,phi = self.sph_coords = Get_Spherical_Coords(v)
        self.r = r
        self.phi = phi
        self.theta = theta
        #self.coords = r,phi,theta
        self.tol = min(10**-3,self.Mag()*10**-3)
    def __str__(self):
        return self.v.__str__()
    def __repr__(self):
        return 'sph_vector '+'\t'+ self.v.__str__()
    def __getitem__(self,key):
        return self.v[key]
    def set_eq_tol(self,tol):
        """set tolerance for __eq__ checking"""
        self.tol = tol
    def Mag(self):
        #print 'Mag call'
        return Mag(self.v)
    def Cross(self,other):
        return sph_vector(Cross(self.v[:],other.v[:]))
    
        
    def __eq__(self,other):
        tol = max([self.tol,other.tol])
        if (self-other).Mag()<tol:
            return True
        else:
            return False
    def __hash__(self):
        r = decimal.Decimal(self.r.__str__())
        phi = decimal.Decimal(self.phi.__str__())
        theta = decimal.Decimal(self.theta.__str__())
        rtol = decimal.Decimal(self.tol.__str__())
        intr = int((r/rtol).to_integral(rounding = 'ROUND_HALF_DOWN'))
        inttheta = int((theta/rtol).to_integral(rounding = 'ROUND_HALF_DOWN'))
        intphi = int((phi/rtol).to_integral(rounding = 'ROUND_HALF_DOWN'))
        return hash(tuple([intr,inttheta,intphi]))
    def __cmp__(self,other):        
        if (self-other).Mag()<self.tol:
            return 0
        else:
            if self.Mag()-other.Mag()<0:
                return -1
            else:
                return 1
                       
            
    def __len__(self):
        return len(self.v)
    def __sub__(self,other):
        return self+(-1*other)
    def __add__(self,other):
        return sph_vector(self.v + other.v)
    def __mul__(self,other):
        if hasattr(other,'__len__'):
            tmp=[]
            for i in range(len(other)):
                tmp.append(other[i]*self[i])
            return sph_vector(tmp)
        else:
            return sph_vector(self.v*other)
    def __rmul__(self,other):
        return self.__mul__(other)
def Find_A_Rotation(v1,v2):
    v1 = Unit(v1)
    v2 = Unit(v2)
    
    axis = Unit(Cross(v1,v2))
    costheta = dot(v1,v2)
    theta = acos(costheta)
    return Rodrigues3(axis,theta)
    
    
def From_Spherical_Coords((r,theta,phi),p_ = e3func(),q_ = e1func(),r_ = e2func()):
    
    return r*erhofunc(theta,phi,p=p_,q=q_,r=r_)
def Proj_func(axis):
    return Ifunc() - Bond(Unit(axis),Unit(axis))
def Get_Spherical_Coords(v,e3 = e3func(),e1 = e1func(),e2 = e2func(),tol = 10**-6):
    #print 'v',v
    r = Mag(v)   
    v = Unit(v)
    #tol = 10**-6
    try:
        
        theta = acos(dot(v,Unit(e3)))
    except ValueError:
        if abs(dot(v,Unit(e3)))-1<tol:
            theta = acos(1*numpy.sign(dot(v,Unit(e3))))
        #print v,e3,dot(v,e3),dot(v,Unit(e3)).__repr__()
        #theta = acos(dot(v,Unit(e3)))
    #    print 'trieddot',v,e3
    proj = Ifunc() - Bond(e3,e3)
    v_proj = dot(proj,v)
    if Mag(v_proj)<tol:
        phi = 0
        return r,theta,phi
        
    v_proj = Unit(v_proj)
    if dot(v_proj,e1)<0:
        if dot(v_proj,e2)<0:
            #quadrant 3
            phi = acos(dot(v_proj,-e1))+math.pi
        else:
            #quadrant 2
            phi = acos(dot(v_proj,e2))+math.pi/2.
    else:
        if dot(v_proj,e2)<0:
            #quadrant 4
            phi = 2*math.pi - acos(dot(v_proj,e1))
        else:
            #quadrant 1
            phi = acos(dot(v_proj,e1))
            
    return r,theta,phi
def Get_Cyl_Polar_Coords(arg,e1 = e1func(),e2 = e2func()):
    r = Mag(arg)
    tmp = Unit(arg)
    
    if dot(tmp,e1)<0:
        if dot(tmp,e2)<0:
            #quadrant 3
            phi = acos(dot(tmp,-e1))+math.pi
        else:
            #quadrant 2
            phi = acos(dot(tmp,e2))+math.pi/2.
    else:
        if dot(tmp,e2)<0:
            #quadrant 4
            phi = 2*math.pi - acos(dot(tmp,e1))
        else:
            #quadrant 1
            phi = acos(dot(tmp,e1))
    return r,phi
def Fstar_Rparams(Rparams,U):
    R_ = Rodrigues3_angaxis(Rparams)
    F = dot(R_,U)
    return Star(F)
def Star(A):
    return det(A)*(inv(A.T))
def Inner_Prod(A,a,b):
    return dot(dot(A,b),a)
def Cross(a,b):
    a1,a2,a3 = dot(a,e1func()),dot(a,e2func()),dot(a,e3func())
    b1,b2,b3 = dot(b,e1func()),dot(b,e2func()),dot(b,e3func())
    return (a2*b3 - a3*b2)*e1func() + (a3*b1 - a1*b3)*e2func() + (a1*b2 - a2*b1)*e3func()

def Rodrigues3_angaxis((thetapx,thetapy,thetapz)):
    npi = numpy.array([thetapx,thetapy,thetapz])
    px,py,pz = Unit(npi)
    theta = Mag(npi)
    w = numpy.array([[0,-pz,py],[pz,0,-px],[-py,px,0]])
    return Ifunc()+w*sin(theta)+dot(w,w)*(1-cos(theta))

def Rodrigues3_((px,py,pz,theta)):
    w = numpy.array([[0,-pz,py],[pz,0,-px],[-py,px,0]])
    return Ifunc()+w*sin(theta)+dot(w,w)*(1-cos(theta))

def Rodrigues3(p,theta):
    px,py,pz = p
    w = numpy.array([[0,-pz,py],[pz,0,-px],[-py,px,0]])
    return Ifunc()+w*sin(theta)+dot(w,w)*(1-cos(theta))
def dRdtheta(p,theta):
    px,py,pz = p
    w = numpy.array([[0,-pz,py],[pz,0,-px],[-py,px,0]])
    return w*cos(theta)+dot(w,w)*sin(theta)

def dRdn(p,theta):
    px,py,pz = p
    dRij = num.zeros([3,3,3])
    def perm_loop1(i1,i2,i3,n):
        ret = 0
        for i in range(3):
            for j in range(3):
                ret+=perm(i1,i,i3)*perm(i,i2,j)*n[j]
        return ret
    
    def perm_loop2(i1,i2,i3,n):
        ret = 0
        for i in range(3):
            for j in range(3):
                ret+=perm(i1,i,j)*n[j]*perm(i,i2,i3)
        return ret
    
    for i in range(3):
        for m in range(3):
            for l in range(3):
                dRij[i,m,l]=-perm(i,m,l)*sin(theta) + (perm_loop1(i,m,l,p) + perm_loop2(i,m,l,p))*(1-cos(theta))
    return dRij

def two_norm(A,B):
    return sqrt(numpy.trace(dot(A,B.T)))

def two_norm_diff(A,B):
    return two_norm(A-B,A-B)

def Rodrigues2(theta,phi,w):
    p = erhofunc(theta,phi)
    a = e1func()+e2func()-e3func()
    if dot(Unit(a),p)<.9999:
        q = Unit(Cross(p,a))
        r = Cross(q,p)
    else:
        a = e1func()
        q = Unit(Cross(p,a))
        r = Cross(q,p)
    return Rodrigues(w,p,q,r)
def EigenDict(C):
    eval,evec = eig(C)
    b = inv_dict()
    c = {}
    #convert to real, allocate eigenvectors
    for i in range(len(eval)):
        b[i] = float(eval[i])
        c[i] = evec[:,i]
    
    b.Flip()
    
    keys = b.keys()[:]
    keys.sort()
    
    order = []
    if b.Contains_Extended() == False:
        for j in range(len(eval)):
            order.append(b[keys[j]])
    else:
        raise AttributeError,"multiples"
        for j in range(len(b.keys())):
            order.append(b[keys[j]])
    final_lambda_dict = {}
    final_evec_dict = {}
    b.Flip()
    for k in range(len(order)):
        elem = order[k]
        final_lambda_dict[k] = b[elem]
        final_evec_dict[k] = c[elem]
    ret = {}
    ret['eval']=final_lambda_dict
    ret['evec']=final_evec_dict
    tol = 10**-8
    for i in range(len(ret['evec'].keys())):
        key = ret['evec'].keys()[i]
        for j in range(len(ret['evec'].keys())):
            key2 = ret['evec'].keys()[j]
            if i!=j:
                if dot(ret['evec'][key],ret['evec'][key2]) > tol:
                    warnings.warn('not orthogonal')
        
    return ret
  

def Euler_2(a,b,c,e1=e1func(),e2 = e2func(),e3=e3func()):
    basis = numpy.array([e1,e2,e3])
    R1 = Rodrigues3(e3,a)
    e1pr,e2pr,e3pr = newbasis = dot(R1,basis)
    R2 = Rodrigues3(e1pr,b)
    e1pr_,e2pr_,e3pr_ = newbasis_ = dot(R2,newbasis)
    R3 = Rodrigues3(e3pr_,c)
    return dot(dot(R3,R2),R1)
class Hash_Vector:
    def __init__(self,v,tol = 1e-6):
        self.array = v[:]
        self.tol = tol
    def __len__(self):
        return len(self.array)
    def __sub__(self,other):
        return self+(-1*other)
    def __add__(self,other):
        return Hash_3Vector(self.array + other.array)
    def __mul__(self,other):
        if hasattr(other,'__len__'):
            tmp=[]
            for i in range(len(other)):
                tmp.append(other[i]*self[i])
            return Hash_3Vector(tmp)
        else:
            return Hash_3Vector(self.array*other)
    def __rmul__(self,other):
        return self.__mul__(other)
        
    
    def Set_Tol(self,tol):
        self.tol = tol
    def Mag(self):
        return Mag(self.array)
    def __getitem__(self,key):
        return self.array[key]
    def __str__(self):
        return self.array.__str__()
    def __repr__(self):
        return 'Hash_3Vector'+'\t'+self.array.__repr__()
    def __eq__(self,other):
        tol = max([self.tol,other.tol])
        if (self-other).Mag()<tol:
            return True
        else:
            return False
        
    def __cmp__(self,other):
        if (self-other).Mag()<self.tol:
            return 0
        else:
            if self.Mag()-other.Mag()<0:
                return -1
            else:
                return 1
    
    def __hash__(self):
        
        #a,b,c = self.array
        tol = self.tol
        ints = []
        for element in self.array:
            ints.append(Int_Map(element/tol))            
        #v1 = Int_Map(a/tol)
        #v2 = Int_Map(b/tol)
        #v3 = Int_Map(c/tol)
        return hash(tuple(ints))
    

    
class Hash_3Vector:
    def __init__(self,v,tol = 10**-6):
        """for use with dictionaries, tolerance by default is 10**-6)"""
        
        self.array = v[:]
        self.tol = tol
    def __len__(self):
        return len(self.array)
    def __sub__(self,other):
        return self+(-1*other)
    def __add__(self,other):
        return Hash_3Vector(self.array + other.array)
    def __mul__(self,other):
        if hasattr(other,'__len__'):
            tmp=[]
            for i in range(len(other)):
                tmp.append(other[i]*self[i])
            return Hash_3Vector(tmp)
        else:
            return Hash_3Vector(self.array*other)
    def __rmul__(self,other):
        return self.__mul__(other)
        
    
    def Set_Tol(self,tol):
        self.tol = tol
    def Mag(self):
        return Mag(self.array)
    def __getitem__(self,key):
        return self.array[key]
    def __str__(self):
        return self.array.__str__()
    def __repr__(self):
        return 'Hash_3Vector'+'\t'+self.array.__repr__()
    def __eq__(self,other):
        tol = max([self.tol,other.tol])
        if (self-other).Mag()<tol:
            return True
        else:
            return False
        
    def __cmp__(self,other):
        if (self-other).Mag()<self.tol:
            return 0
        else:
            if self.Mag()-other.Mag()<0:
                return -1
            else:
                return 1
    
    def __hash__(self):
        a,b,c = self.array
        tol = self.tol
        v1 = Int_Map(a/tol)
        v2 = Int_Map(b/tol)
        v3 = Int_Map(c/tol)
        return hash(tuple([v1,v2,v3]))
    
    
def Quadrant2D(v,e1=e1func(),e2=e2func()):
    if dot(v,e1)>=0:
        if dot(v,e2)>=0:
            return 1
        else:
            return 4
    else:
        if dot(v,e2)>=0:
            return 2
        else:
            return 3



def polarDecomposition(F,C=False):
    if C==True:
        C = F
    else:
        C = dot(F.T,F)
    eval,evec = eig(C)
    l1_sqr,l2_sqr,l3_sqr = numpy.asarray(eval,dtype = 'float')
    u1,u2,u3 = evec[:,0],evec[:,1],evec[:,2]
    U = sqrt(l1_sqr)*numpy.outer(u1,u1) + sqrt(l2_sqr)*numpy.outer(u2,u2) + sqrt(l3_sqr)*numpy.outer(u3,u3)
    invU = 1./sqrt(l1_sqr)*numpy.outer(u1,u1) + 1./sqrt(l2_sqr)*numpy.outer(u2,u2) + 1./sqrt(l3_sqr)*numpy.outer(u3,u3)
    R = dot(F,invU)

    return R,U
    
def unit_vector(length_, index):
    a = numpy.zeros(length_)
    a[index]=1
    return a
def spot_to_vec(spot, tmpdg,wavelength):
    # circular import below
    import hexrd.XRD.grain_wrap_precession as gw
    #angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
    angCOM= spot.angCOM(useFit=True)
    xyoCOM = num.array(spot.detectorGeom.angToXYO(*angCOM)).flatten()
    new_angCOM = tmpdg.xyoToAng(*xyoCOM)
    tth,eta,ome = new_angCOM
    rI = gw.makeARecipVector((tth,eta,ome), wavelength)
    return rI
