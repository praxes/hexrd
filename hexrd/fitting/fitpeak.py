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
import scipy.optimize as optimize
import hexrd.fitting.peakfunctions as pkfuncs
import scipy.ndimage as imgproc
import copy



#### 1-D Peak Fitting
def estimate_pk_parms_1d(x,f,pktype):
    """
    Gives initial guess of parameters for analytic fit of one dimensional peak
    data.

    Required Arguments:
    x -- (n) ndarray of coordinate positions
    f -- (n) ndarray of intensity measurements at coordinate positions x
    pktype -- string, type of analytic function that will be used to fit the data,
    current options are "gaussian","lorentzian","pvoigt" (psuedo voigt), and
    "split_pvoigt" (split psuedo voigt)
    

    Outputs:
    p -- (m) ndarray containing initial guesses for parameters for the input peaktype 
    (see peak function help for what each parameters corresponds to)
    """

    
    data_max=np.max(f)
#    lbg=np.mean(f[:2])
#    rbg=np.mean(f[:2]) 
    if((f[0]> (0.25*data_max)) and (f[-1]> (0.25*data_max))):#heuristic for wide peaks  
        bg0=0.
    elif (f[0]> (0.25*data_max)): #peak cut off on the left 
        bg0=f[-1]
    elif (f[-1]> (0.25*data_max)): #peak cut off on the right 
        bg0=f[0]
    else:
        bg0=(f[0]+f[-1])/2.    
    #bg1=(rbg-lbg)/(x[-1]-x[0])
    
    cen_index=np.argmax(f) 
    x0=x[cen_index]    
    A=data_max-bg0#-(bg0+bg1*x0)  
    
    num_pts=len(f)
   
    #checks for peaks that are cut off
    if cen_index == (num_pts-1):
        FWHM=x[cen_index]-x[np.argmin(np.abs(f[:cen_index]-A/2.))]#peak cut off on the left 
    elif cen_index == 0:
        FWHM=x[cen_index+np.argmin(np.abs(f[cen_index+1:]-A/2.))]-x[0] #peak cut off on the right 
    else:    
        FWHM=x[cen_index+np.argmin(np.abs(f[cen_index+1:]-A/2.))]-x[np.argmin(np.abs(f[:cen_index]-A/2.))]
        
    if FWHM <=0:##uh,oh something went bad
        FWHM=(x[-1]-x[0])/4. #completely arbitrary, set peak width to 1/4 window size

    
    
    if pktype=='gaussian' or pktype=='lorentzian':
        p=[A,x0,FWHM,bg0,0.]
    elif pktype=='pvoigt':
        p=[A,x0,FWHM,0.5,bg0,0.]
    elif pktype=='split_pvoigt':
        p=[A,x0,FWHM,FWHM,0.5,0.5,bg0,0.]
        
    p=np.array(p)
    return p
      
    
def fit_pk_parms_1d(p0,x,f,pktype):
    """
    Performs least squares fit to find parameters for 1d analytic functions fit
    to diffraction data

    Required Arguments:
    p0 -- (m) ndarray containing initial guesses for parameters for the input peaktype
    x -- (n) ndarray of coordinate positions
    f -- (n) ndarray of intensity measurements at coordinate positions x
    pktype -- string, type of analytic function that will be used to fit the data,
    current options are "gaussian","lorentzian","pvoigt" (psuedo voigt), and
    "split_pvoigt" (split psuedo voigt)
    

    Outputs:
    p -- (m) ndarray containing fit parameters for the input peaktype (see peak function
    help for what each parameters corresponds to)
    
    
    Notes:
    1. Currently no checks are in place to make sure that the guess of parameters
    has a consistent number of parameters with the requested peak type
    """    
    
    
    fitArgs=(x,f,pktype)
    
    ftol=1e-6
    xtol=1e-6    
    
    weight=np.max(f)*10.#hard coded should be changed
    
    if pktype == 'gaussian':
        p, outflag = optimize.leastsq(fit_pk_obj_1d, p0, args=fitArgs,Dfun=eval_pk_deriv_1d,ftol=ftol,xtol=xtol) 
    elif pktype == 'lorentzian':
        p, outflag = optimize.leastsq(fit_pk_obj_1d, p0, args=fitArgs,Dfun=eval_pk_deriv_1d,ftol=ftol,xtol=xtol) 
    elif pktype == 'pvoigt':
        lb=[p0[0]*0.5,np.min(x),0.,      0., 0.,None]
        ub=[p0[0]*2.0,np.max(x),4.*p0[2],1., 2.*p0[4],None]         
                
        fitArgs=(x,f,pktype,weight,lb,ub)
        p, outflag = optimize.leastsq(fit_pk_obj_1d_bnded, p0, args=fitArgs,ftol=ftol,xtol=xtol) 
    elif pktype == 'split_pvoigt':
        lb=[p0[0]*0.5,np.min(x),0.,      0.,      0., 0., 0.,None]
        ub=[p0[0]*2.0,np.max(x),4.*p0[2],4.*p0[2],1., 1., 2.*p0[4],None]         
        
        p, outflag = optimize.leastsq(fit_pk_obj_1d_bnded, p0, args=fitArgs,ftol=ftol,xtol=xtol) 
        
    elif pktype == 'tanh_stepdown':
        p, outflag = optimize.leastsq(fit_pk_obj_1d, p0, args=fitArgs,ftol=ftol,xtol=xtol)
    else:
        p=p0        
        print('non-valid option, returning guess')
           
    
    if np.any(np.isnan(p)):
        p=p0
        print('failed fitting, returning guess')
    
    return p
 
def eval_pk_deriv_1d(p,x,y0,pktype):  

    if pktype == 'gaussian':
        d_mat=pkfuncs.gaussian1d_deriv(p,x)
    elif pktype == 'lorentzian':
        d_mat=pkfuncs.lorentzian1d_deriv(p,x)
    
    return d_mat.T

    
def fit_pk_obj_1d(p,x,f0,pktype):  
    if pktype == 'gaussian':
        f=pkfuncs.gaussian1d(p,x)
    elif pktype == 'lorentzian':
        f=pkfuncs.lorentzian1d(p,x)
    elif pktype == 'pvoigt':
        f=pkfuncs.pvoigt1d(p,x)
    elif pktype == 'split_pvoigt':
        f=pkfuncs.split_pvoigt1d(p,x)
    elif pktype == 'tanh_stepdown':
        f=pkfuncs.tanh_stepdown_nobg(p,x)
    
    resd = f-f0
    return resd


def fit_pk_obj_1d_bnded(p,x,f0,pktype,weight,lb,ub):  
    if pktype == 'gaussian':
        f=pkfuncs.gaussian1d(p,x)
    elif pktype == 'lorentzian':
        f=pkfuncs.lorentzian1d(p,x)
    elif pktype == 'pvoigt':
        f=pkfuncs.pvoigt1d(p,x)
    elif pktype == 'split_pvoigt':
        f=pkfuncs.split_pvoigt1d(p,x)
    
    num_data=len(f)
    num_parm=len(p)
    resd=np.zeros(num_data+num_parm)
    #tub bnds implementation    
    
    resd[:num_data] = f-f0
    for ii in range(num_parm):
        if lb[ii] is not None:        
            resd[num_data+ii]=weight*np.max([-(p[ii]-lb[ii]),0.,(p[ii]-ub[ii])])
    
    
    return resd

#### 2-D Peak Fitting

def estimate_pk_parms_2d(x,y,f,pktype):
    """
    Gives initial guess of parameters for analytic fit of two dimensional peak
    data.

    Required Arguments:
    x -- (n x 0) ndarray of coordinate positions for dimension 1 (numpy.meshgrid formatting)
    y -- (n x 0) ndarray of coordinate positions for dimension 2 (numpy.meshgrid formatting)
    f -- (n x 0) ndarray of intensity measurements at coordinate positions x and y
    pktype -- string, type of analytic function that will be used to fit the data,
    current options are "gaussian", "gaussian_rot" (gaussian with arbitrary axes) and 
    "split_pvoigt_rot" (split psuedo voigt with arbitrary axes)
    

    Outputs:
    p -- (m) ndarray containing initial guesses for parameters for the input peaktype
    (see peakfunction help for more information)
    """


    
    bg0=np.mean([f[0,0],f[-1,0],f[-1,-1],f[0,-1]])
    bg1x=(np.mean([f[-1,-1],f[0,-1]])-np.mean([f[0,0],f[-1,0]]))/(x[0,-1]-x[0,0])
    bg1y=(np.mean([f[-1,-1],f[-1,0]])-np.mean([f[0,0],f[0,-1]]))/(y[-1,0]-y[0,0])
    
    fnobg=f-(bg0+bg1x*x+bg1y*y)    
    
    labels,numlabels=imgproc.label(fnobg>np.max(fnobg)/2.)
    
    #looks for the largest peak
    areas=np.zeros(numlabels)
    for ii in np.arange(1,numlabels+1,1):
        areas[ii-1]= np.sum(labels==ii)
    
    peakIndex=np.argmax(areas)+1  
    
    
#    #currently looks for peak closest to center
#    dist=np.zeros(numlabels)
#    for ii in np.arange(1,numlabels+1,1):
#        dist[ii-1]= ######
#    
#    peakIndex=np.argmin(dist)+1
    
    FWHMx=np.max(x[labels==peakIndex])-np.min(x[labels==peakIndex])
    FWHMy=np.max(y[labels==peakIndex])-np.min(y[labels==peakIndex])
    
    coords=imgproc.maximum_position(fnobg, labels=labels, index=peakIndex)
    A=imgproc.maximum(fnobg, labels=labels, index=peakIndex)
    x0=x[coords]
    y0=y[coords]
    
    if pktype=='gaussian':
        p=[A,x0,y0,FWHMx,FWHMy,bg0,bg1x,bg1y]
    elif pktype=='gaussian_rot':
        p=[A,x0,y0,FWHMx,FWHMy,0.,bg0,bg1x,bg1y]
    elif pktype=='split_pvoigt_rot':
        p=[A,x0,y0,FWHMx,FWHMx,FWHMy,FWHMy,0.5,0.5,0.5,0.5,0.,bg0,bg1x,bg1y]
        
    p=np.array(p)
    return p


def fit_pk_parms_2d(p0,x,y,f,pktype):
    """
    Performs least squares fit to find parameters for 2d analytic functions fit
    to diffraction data

    Required Arguments:
    p0 -- (m) ndarray containing initial guesses for parameters for the input peaktype
    x -- (n x 0) ndarray of coordinate positions for dimension 1 (numpy.meshgrid formatting)
    y -- (n x 0) ndarray of coordinate positions for dimension 2 (numpy.meshgrid formatting)
    f -- (n x 0) ndarray of intensity measurements at coordinate positions x and y
    pktype -- string, type of analytic function that will be used to fit the data,
    current options are "gaussian", "gaussian_rot" (gaussian with arbitrary axes) and 
    "split_pvoigt_rot" (split psuedo voigt with arbitrary axes)
    

    Outputs:
    p -- (m) ndarray containing fit parameters for the input peaktype (see peak function
    help for what each parameters corresponds to)
    
    
    Notes:
    1. Currently no checks are in place to make sure that the guess of parameters
    has a consisten number of parameters with the requested peak type
    """    


    fitArgs=(x,y,f,pktype)
    ftol=1e-9
    xtol=1e-9
    
    if pktype == 'gaussian':
        p, outflag = optimize.leastsq(fit_pk_obj_2d, p0, args=fitArgs,ftol=ftol, xtol=xtol) 
    elif pktype == 'gaussian_rot':
        p, outflag = optimize.leastsq(fit_pk_obj_2d, p0, args=fitArgs,ftol=ftol, xtol=xtol)
    elif pktype == 'split_pvoigt_rot':
        p, outflag = optimize.leastsq(fit_pk_obj_2d, p0, args=fitArgs,ftol=ftol, xtol=xtol) 
        
    
    if np.any(np.isnan(p)):
        p=p0
    
    return p

def fit_pk_obj_2d(p,x,y,f0,pktype):    
    if pktype == 'gaussian':
        f=pkfuncs.gaussian2d(p,x,y)
    elif pktype == 'gaussian_rot':       
        f=pkfuncs.gaussian2d_rot(p,x,y)
    elif pktype == 'split_pvoigt_rot':       
        f=pkfuncs.split_pvoigt2d_rot(p,x,y)
         
    resd = f-f0
    return resd.flatten()
    


#### Extra Utilities

def goodness_of_fit(f,f0):
    """
    Calculates two scalar measures of goodness of fit

    Required Arguments:
    f0 -- (n x 0) ndarray of intensity measurements at coordinate positions
    f -- (n x 0) ndarray of fit intensity values at coordinate positions

    Outputs:
    R -- (1) goodness of fit measure which is sum(error^2)/sum(meas^2)
    Rw -- (1) goodness of fit measure weighted by intensity sum(meas*error^2)/sum(meas^3)
    
    
    """ 


    R=np.sum((f-f0)**2)/np.sum(f0**2)
    Rw=np.sum(np.abs(f0*(f-f0)**2))/np.sum(np.abs(f0**3))
    
    return R, Rw    
    
        
    