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
#import code
import copy

gauss_width_fact=2.*np.sqrt(2.*np.log(2.))
lorentz_width_fact=2.

#### 1-D Gaussian Functions
def _unit_gaussian(p,x):#Split the unit gaussian so this can be called for 2d and 3d functions
    """
    Required Arguments:
    p -- (m) [x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 



    x0=p[0]
    FWHM=p[1]
    sigma=FWHM/gauss_width_fact 
    
    f=np.exp(-(x-x0)**2/(2.*sigma**2.))
    return f

def _gaussian1d_no_bg(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    A=p[0]
    f=A*_unit_gaussian(p[[1,2]],x)    
    return f


def gaussian1d(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    bg0=p[3]
    bg1=p[4]  
    
    f=_gaussian1d_no_bg(p[:3],x)+bg0+bg1*x    
   
    return f
    
    
def _gaussian1d_no_bg_deriv(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    d_mat -- (3 x n) ndarray of derivative values at positions x
    """ 


    x0=p[1]
    FWHM=p[2]
    
    sigma=FWHM/gauss_width_fact 

    dydx0=_gaussian1d_no_bg(p,x)*((x-x0)/(sigma**2.))
    dydA=_unit_gaussian(p[[1,2]],x) 
    dydFWHM=_gaussian1d_no_bg(p,x)*((x-x0)**2./(sigma**3.))/gauss_width_fact
    
    d_mat=np.zeros((len(p),len(x)))
    
    d_mat[0,:]=dydA
    d_mat[1,:]=dydx0
    d_mat[2,:]=dydFWHM
    
    return d_mat
    
def gaussian1d_deriv(p,x): 
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    d_mat -- (5 x n) ndarray of derivative values at positions x
    """ 


    d_mat=np.zeros((len(p),len(x)))
    d_mat[0:3,:]=_gaussian1d_no_bg_deriv(p[0:3],x)
    d_mat[3,:]=1.    
    d_mat[4,:]=x
    
    return d_mat


#### 1-D Lorentzian Functions

def _unit_lorentzian(p,x):#Split the unit function so this can be called for 2d and 3d functions
    """
    Required Arguments:
    p -- (m) [x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    x0=p[0]
    FWHM=p[1]
    gamma=FWHM/lorentz_width_fact   
    
    f= gamma**2 / ((x-x0)**2 + gamma**2)
    return f

def _lorentzian1d_no_bg(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    A=p[0]
    f= A*_unit_lorentzian(p[[1,2]],x) 
   
    return f
    
def lorentzian1d(p,x):
    """
    Required Arguments:
    p -- (m) [x0,FWHM,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """     
    
    bg0=p[3]
    bg1=p[4]

    f=_lorentzian1d_no_bg(p[:3],x)+bg0+bg1*x 
   
    return f


def _lorentzian1d_no_bg_deriv(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    d_mat -- (3 x n) ndarray of derivative values at positions x
    """     
    
    x0=p[1]
    FWHM=p[2]
    
    gamma=FWHM/lorentz_width_fact 

    dydx0=_lorentzian1d_no_bg(p,x)*((2.*(x-x0))/((x-x0)**2 + gamma**2))
    dydA=_unit_lorentzian(p[[1,2]],x) 
    dydFWHM=_lorentzian1d_no_bg(p,x)*((2.*(x-x0)**2.)/(gamma*((x-x0)**2 + gamma**2)))/lorentz_width_fact
    
    d_mat=np.zeros((len(p),len(x)))
    d_mat[0,:]=dydA    
    d_mat[1,:]=dydx0
    d_mat[2,:]=dydFWHM
    
    return d_mat
    
def lorentzian1d_deriv(p,x): 
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    d_mat -- (5 x n) ndarray of derivative values at positions x
    """ 

    d_mat=np.zeros((len(p),len(x)))
    d_mat[0:3,:]=_lorentzian1d_no_bg_deriv(p[0:3],x)
    d_mat[3,:]=1.    
    d_mat[4,:]=x
    
    return d_mat


#### 1-D Psuedo Voigt Functions
def _unit_pvoigt1d(p,x):#Split the unit function so this can be called for 2d and 3d functions
    """
    Required Arguments:
    p -- (m) [x0,FWHM,n]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 


    n=p[2] 
    
    f=(n*_unit_gaussian(p[:2],x)+(1.-n)*_unit_lorentzian(p[:2],x))
    return f

def _pvoigt1d_no_bg(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM,n]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    A=p[0]    
    f=A*_unit_pvoigt1d(p[[1,2,3]],x)
    return f

def pvoigt1d(p,x):   
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM,n,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """ 

    bg0=p[4]
    bg1=p[5]

    f=_pvoigt1d_no_bg(p[:4],x)+bg0+bg1*x 
   
    return f

#### 1-D Split Psuedo Voigt Functions
def _split_pvoigt1d_no_bg(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM-,FWHM+,n-,n+]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """

    A=p[0]
    x0=p[1]

    f=np.zeros(x.shape[0])    
    
    #Define halves, using gthanorequal and lthan, choice is arbitrary
    xr=x>=x0
    xl=x<x0      
            
    #+  
    r=np.where(xr)  
    f[xr]=A*_unit_pvoigt1d(p[[1,3,5]],x[xr])   
    
    #-
    l=np.where(xl)  
    f[xl]=A*_unit_pvoigt1d(p[[1,2,4]],x[xl])     

    return f

def split_pvoigt1d(p,x):   
    """
    Required Arguments:
    p -- (m) [A,x0,FWHM-,FWHM+,n-,n+,c0,c1]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """


    bg0=p[6]
    bg1=p[7]

    f=_split_pvoigt1d_no_bg(p[:6],x)+bg0+bg1*x 
   
    return f

#### Tanh Step Down
def tanh_stepdown_nobg(p,x):
    """
    Required Arguments:
    p -- (m) [A,x0,w]
    x -- (n) ndarray of coordinate positions 

    Outputs:
    f -- (n) ndarray of function values at positions x
    """

    A=p[0]
    x0=p[1]
    w=p[2]

    f=A*(0.5*(1.-np.tanh((x-x0)/w)))    
        

    return f


#### 2-D Rotation Coordinate Transform

def _2d_coord_transform(theta,x0,y0,x,y):
    xprime=np.cos(theta)*x+np.sin(theta)*y
    yprime=-np.sin(theta)*x+np.cos(theta)*y

    x0prime=np.cos(theta)*x0+np.sin(theta)*y0
    y0prime=-np.sin(theta)*x0+np.cos(theta)*y0
    
    return x0prime, y0prime, xprime, yprime



#### 2-D Gaussian Function
def _gaussian2d_no_bg(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx,FWHMy]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 1

    Outputs:
    f -- (n x 0) ndarray of function values at positions (x,y)
    """ 

    A=p[0]    
    f=A*_unit_gaussian(p[[1,3]],x)*_unit_gaussian(p[[2,4]],y)  
    return f
    
    
def _gaussian2d_rot_no_bg(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx,FWHMy,theta]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 

    theta=p[5]    
  
    x0prime, y0prime, xprime, yprime=_2d_coord_transform(theta,p[1],p[2],x,y)
    
    newp=copy.copy(p)#this copy was needed so original parameters set isn't changed     
    
    newp[1]=x0prime
    newp[2]=y0prime
    
    f=_gaussian2d_no_bg(newp[:5],xprime,yprime) 

    return f 
    
    
def gaussian2d_rot(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx,FWHMy,theta,c0,c1x,c1y]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 

    bg0=p[6]
    bg1x=p[7]
    bg1y=p[8]
 
    f=_gaussian2d_rot_no_bg(p[:6],x,y)+(bg0+bg1x*x+bg1y*y)  
    return f 

def gaussian2d(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx,FWHMy,c0,c1x,c1y]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 

    bg0=p[5]
    bg1x=p[6]
    bg1y=p[7]
 
    f=_gaussian2d_no_bg(p[:5],x,y)+(bg0+bg1x*x+bg1y*y)  
    return f 


#### 2-D Split Psuedo Voigt Function

def _split_pvoigt2d_no_bg(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx-,FWHMx+,FWHMy-,FWHMy+,nx-,nx+,ny-,ny+]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 

    A=p[0]
    x0=p[1]
    y0=p[2]


    f=np.zeros([x.shape[0],x.shape[1]])    
    
    #Define quadrants, using gthanorequal and lthan, choice is arbitrary
    xr=x>=x0
    xl=x<x0      
    yr=y>=y0
    yl=y<y0         
    
    #++  
    q1=np.where(xr & yr)  
    f[q1]=A*_unit_pvoigt1d(p[[1,4,8]],x[q1])*_unit_pvoigt1d(p[[2,6,10]],y[q1])    
    
    #+-
    q2=np.where(xr & yl)  
    f[q2]=A*_unit_pvoigt1d(p[[1,4,8]],x[q2])*_unit_pvoigt1d(p[[2,5,9]],y[q2])  
    
    #-+
    q3=np.where(xl & yr)  
    f[q3]=A*_unit_pvoigt1d(p[[1,3,7]],x[q3])*_unit_pvoigt1d(p[[2,6,10]],y[q3]) 
    
    #--    
    q4=np.where(xl & yl)  
    f[q4]=A*_unit_pvoigt1d(p[[1,3,7]],x[q4])*_unit_pvoigt1d(p[[2,5,9]],y[q4])     
    

    return f


def _split_pvoigt2d_rot_no_bg(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx-,FWHMx+,FWHMy-,FWHMy+,nx-,nx+,ny-,ny+,theta]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 


    theta=p[11]    
  
    x0prime, y0prime, xprime, yprime=_2d_coord_transform(theta,p[1],p[2],x,y)
    
    newp=copy.copy(p)#this copy was needed so original parameters set isn't changed   
    
    newp[1]=x0prime
    newp[2]=y0prime
    
    f=_split_pvoigt2d_no_bg(newp[:11],xprime,yprime) 

    return f 

def split_pvoigt2d_rot(p,x,y):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,FWHMx-,FWHMx+,FWHMy-,FWHMy+,nx-,nx+,ny-,ny+,theta,c0,c1x,c1y]
    x -- (n x o) ndarray of coordinate positions for dimension 1
    y -- (n x o) ndarray of coordinate positions for dimension 2

    Outputs:
    f -- (n x o) ndarray of function values at positions (x,y)
    """ 


    bg0=p[12]
    bg1x=p[13]
    bg1y=p[14]
    
    f=_split_pvoigt2d_rot_no_bg(p[:12],x,y)+(bg0+bg1x*x+bg1y*y) 

    return f 

#### 3-D Gaussian Function
def _gaussian3d_no_bg(p,x,y,z):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,z0,FWHMx,FWHMy,FWHMz]
    x -- (n x o x q) ndarray of coordinate positions for dimension 1
    y -- (n x o x q) ndarray of coordinate positions for dimension 2
    y -- (z x o x q) ndarray of coordinate positions for dimension 3

    Outputs:
    f -- (n x o x q) ndarray of function values at positions (x,y)
    """ 

    A=p[0]    
    f=A*_unit_gaussian(p[[1,4]],x)*_unit_gaussian(p[[2,5]],y)*_unit_gaussian(p[[3,6]],z)    
    return f
    
    


def gaussian3d(p,x,y,z):
    """
    Required Arguments:
    p -- (m) [A,x0,y0,z0,FWHMx,FWHMy,FWHMz,c0,c1x,c1y,c1z]
    x -- (n x o x q) ndarray of coordinate positions for dimension 1
    y -- (n x o x q) ndarray of coordinate positions for dimension 2
    y -- (z x o x q) ndarray of coordinate positions for dimension 3

    Outputs:
    f -- (n x o x q) ndarray of function values at positions (x,y,z)
    """ 

    bg0=p[7]
    bg1x=p[8]
    bg1y=p[9]
    bg1z=p[10]
 
    f=_gaussian3d_no_bg(p[:5],x,y)+(bg0+bg1x*x+bg1y*y+bg1z*z)  
    return f 
