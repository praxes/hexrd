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
import sys
import os
import copy
import time
import math
import shelve
import exceptions
import StringIO
import ctypes
import multiprocessing

import numpy as num
from matplotlib import cm
from scipy import ndimage
from scipy import optimize
from scipy import sparse
try:
    import scikits.image
    haveSciKits = True
except:
    haveSciKits = False

from hexrd import arrayUtil
from hexrd.quadrature import q3db, q2db
from hexrd import tens
from hexrd import valUnits
from hexrd.valUnits import toFloat
from hexrd import matrixUtils as mUtil
from hexrd.matrixUtils import rowNorm
from hexrd import plotWrap
        
from hexrd import XRD
from hexrd.XRD import crystallography
from hexrd.XRD import detector
from hexrd.XRD import xrdBase
from hexrd.XRD.xrdBase import getGaussNDParams
from hexrd.XRD import Rotations
from hexrd.XRD.Rotations import mapAngle

# from hexrd import uncertainty_analysis
        

debugDflt = False
debugMulti = True
debugFit = False
debugMasked = False
debugFrameRefs = False
debugFinalize = False

structureNDI_label = num.array([
    [0,1,0],
    [1,1,1],
    [0,1,0]
    ])
structureNDI_dilate = num.array([
    [1,1,1],
    [1,1,1],
    [1,1,1]
    ])

spotFitFuncDflt = 'gauss'

def getBin(thisframe, threshold, padSpot):
    bin = num.zeros(thisframe.shape, dtype=bool)
    bin[thisframe > threshold] = True
    if hasattr(thisframe, 'mask'):
        'for some reason, this seems to be necessary'
        bin[thisframe.mask] = False
        if debugMasked:
            nBin = num.sum(bin)
            print 'bin has %d pixels' % (nBin)
            if nBin > 0.1 * thisframe.size:
                print 'going very badly wrong'
    if padSpot:
        ''' do padding out here instead of passing padSpot=True to
        getSpot and dilate after mask application so that dilation
        will fill across rows of dead pixels; '''
        structure=ndimage.generate_binary_structure(thisframe.ndim, thisframe.ndim) # structureNDI_dilate
        bin = ndimage.morphology.binary_dilation(bin, structure=structure, iterations=int(padSpot)) 
    return bin

def dilateObj(obj, shape, nDilate=1):
    newObj = (
        slice(
            max(0, obj[0].start-nDilate), 
            min(obj[0].stop+nDilate, shape[0]-1),
            ),
        slice(
            max(0, obj[1].start-nDilate), 
            min(obj[1].stop+nDilate, shape[1]-1),
            ),
        )
    return newObj

def emptyBox(inpt, obj, dtype=None):
    dtype = dtype or inpt.dtype
    vbox = num.empty(
        (obj[0].stop-obj[0].start, obj[1].stop-obj[1].start),
        dtype=dtype
        )
    return vbox
def copyBox(inpt, obj):
    'deepcopy does not seem to do what we want with masked arrays'
    # vbox   = copy.deepcopy(inpt[useObj])
    vbox = emptyBox(inpt, obj, dtype=inpt.dtype)
    vbox[:,:] = inpt[obj]
    return vbox

def getSpot(inpt, labels, objs, index, keepWithinBBox, padSpot, darkframe=None):
    """
    labels is from ndimage.label;
    objs is from ndimage.find_objects
    """
    
    spot = {}
    
    iSpot = index - 1
    obj = objs[iSpot]
    objDilated = dilateObj(obj, labels.shape)

    'use objDilated instead of obj for atBoundary in case of masked pixels or odd boundary effects'
    atBoundary = (objDilated[0].start == 0) or (objDilated[1].start == 0) or (objDilated[0].stop == labels.shape[0]) or (objDilated[1].stop == labels.shape[1])
    
    if padSpot:
        useObj = objDilated
    else:
        useObj = obj
    'do not want to modify original labels'
    useLabels = copy.deepcopy(labels[useObj]) # do not want to modify original labels!
    
    if padSpot:
        dilatedBool = ndimage.morphology.binary_dilation(useLabels, structure=structureNDI_dilate)
        useLabels[dilatedBool] = index
        
    if hasattr(inpt, 'mask'):
        'make sure mask is respected'
        if debugMasked:
            print 'number of masked pixels : %d' % (num.sum(inpt.mask[useObj]))
        useLabels[inpt.mask[useObj]] = -1
        
    'spot-local x and y values:'
    if keepWithinBBox:
        xl, yl = num.where(useLabels >= 0) # should give all unmasked pixels in bbox
    else:
        xl, yl = num.where(useLabels == index) # labels[useObj] > 0
    vbox = copyBox(inpt, useObj)
    if darkframe is None:
        raise RuntimeError, 'darkframe is None'
    else:
        darkbox = copyBox(darkframe, useObj)
    # v = vbox[xl,yl] # = inpt[x,y]
    x = xl + useObj[0].start
    y = yl + useObj[1].start
    
    iSpot = index - 1
    obj = objs[iSpot]
    spot['index'] = index
    spot['obj']   = useObj
    spot['xl']    = xl
    spot['yl']    = yl
    spot['vbox']  = vbox
    spot['x']     = x
    spot['y']     = y
    spot['atBoundary'] = atBoundary
    spot['darkbox']=darkbox
    #spot['v']    = v # careful if put this back in -- look for where 'vbox' is used everywhere
    return spot

def cullSpotUsingBin(spot, bin):
    '''
    remove data for pixels where bin is true;
    vbox does not change, just change the locations that are marked as belonging with the spot
    '''
    
    useObj = spot['obj']
    xl     = spot['xl']
    yl     = spot['yl']
    
    binSpot = emptyBox(bin, useObj, dtype=bool)
    binSpot[:,:] = False
    binSpot[xl,yl] = True
    binSpot[bin[useObj]] = False
    xl, yl = num.where(binSpot)
    x = xl + useObj[0].start
    y = yl + useObj[1].start
    spot['xl']    = xl
    spot['yl']    = yl
    spot['x']     = x
    spot['y']     = y
    return 

def getValuesOnly(inpt, labels, objs, index):
    """
    labels is from ndimage.label;
    objs is from ndimage.find_objects
    """
    iSpot = index - 1
    obj = objs[iSpot]
    
    xl, yl = num.where(labels[obj] == index)
    # want to avoid keeping a reference to inpt array, hope deepcopy does the trick
    v = copy.deepcopy(inpt[obj][xl,yl])
    
    return v

def getIndices(inpt, labels, obj, index, mode='default', minlabel=0):

    if obj is None:
        obj = tuple([slice(0,labels.shape[iDim]) for iDim in range(len(labels.shape))])
    
    indices = None
    if inpt is not None:
        if hasattr(inpt, 'mask'):
            'make sure mask is respected'
            useLabels = copy.deepcopy(labels[obj]) # do not want to modify original labels!
            useLabels[inpt.mask[obj]] = -1
            if mode == 'default':
                indices = num.where(useLabels == index)
            elif mode == 'exclude':
                indices = num.where(num.logical_and(useLabels >= minlabel, useLabels != index))
            else:
                raise RuntimeError, 'unknown mode : '+str(mode)
            if debugMasked:
                print 'number of masked pixels : %d' % (num.sum(inpt.mask[obj]))
                print 'number of non-masked pixels : %d' % (len(indices[0]))
    if indices is None:
        if mode == 'default':
            indices = num.where(labels[obj] == index) # xl, yl =
        elif mode == 'exclude':    
            indices = num.where(num.logical_and(labels[obj] >= minlabel, labels[obj] != index))
        else:
            raise RuntimeError, 'unknown mode : '+str(mode)
    return indices

def getImCOM(inpt, labels, objs, index, floor=None, getVSum=False):
    """
    labels is from ndimage.label;
    objs is from ndimage.find_objects
    
    set floor for a minimum intensity to use in intensity weighted COM
    
    return sum of intensity as well if getVSum is True
    """
    
    iSpot = index - 1
    obj = objs[iSpot]

    indices = getIndices(inpt, labels, objs[index-1], index)
    
    if inpt is None:
        v = num.ones(len(indices[0]))
        assert not getVSum, \
            'really want to get vSum when inpt is None?'
    else:
        v = inpt[obj][indices]
        if floor is not None:
            v[v < floor] = floor
    vSum = float(num.sum(v))
    if vSum <= 0:
        print 'vSum is %g, from v array %s' % (vSum, str(v))
        raise RuntimeError, 'vSum <= 0'
    
    com = []
    for xL, slc in zip(indices, obj):
        x = xL + slc.start
        # coordsList.append(x)
        com.append( float(num.dot(x, v)) / vSum )
    com = num.array(com)
    
    if getVSum:
        retval = com, vSum
    else:
        retval = com
    return retval

def getObjSize(labels, objs, index):
    """
    labels is from ndimage.label;
    objs is from ndimage.find_objects
    """
    iSpot = index - 1
    obj = objs[iSpot]
    
    objSize = num.sum(labels[obj] == index)
    
    return objSize

def spotFinderSingle(
    thisframe, threshold, 
    minPx,
    keepWithinBBox, padSpot,
    weightedCOM=True,
    pw=None, debug=False, darkframe=None
    ):
    """
    find spots in thisframe;
    threshold can be a scalar or defined over of same dimension as thisframe;
    minPx is the minimum number of pixels to be kept as a spot;
    weightedCOM is true to use thisframe data for weighting center-of-mass position;
    if pw is present, it is used for plotting (assumed to be plotWrap.PlotWrap instance)
    """
    
    needCOM = debug or pw is not None

    bin = getBin(thisframe, threshold, padSpot)
    labels, numSpots = ndimage.label(bin, structureNDI_label)
    
    # find objects gives slices which are like bounding boxes
    objs = ndimage.find_objects(labels)
    
    if needCOM:
        com = num.empty([numSpots,2])
    nPx = num.empty([numSpots],dtype=int)
    for iSpot in range(numSpots):
        index = iSpot+1
        if debug:
            print 'ojbect %d : %s' % (iSpot, str(thisframe[objs[iSpot]]))
        
        """sum on bin gives number of pixels in the spot
        take advantage of objs instead of doing ndimage.sum"""
        # nPx[iSpot] = ndimage.sum(bin, labels=labels, index=index)
        # nPx[iSpot] = len(getValuesOnly(bin, labels, objs, index))
        nPx[iSpot] = getObjSize(labels, objs, index)
        
        if needCOM:
            if weightedCOM:
                # com[iSpot,:] = ndimage.center_of_mass(thisframe, labels=labels, index=index) # slow
                com[iSpot,:] = getImCOM(thisframe, labels, objs, index)
            else:
                com[iSpot,:] = getImCOM(None, labels, objs, index)
        if debug:
            print '    center of mass : '+str(com[iSpot,:])
    
    keepers = nPx >= minPx
    
    # pw.a.axis('on')
    if pw is not None:
        notKeepers = num.ones([numSpots],dtype=bool) - keepers
        lines = pw.a.get_lines()
        for line in lines:
            line.remove()
        pw(   com[keepers,0],    com[keepers,1], marker='o', ls='None', mec='r', mfc='None')
        pw(com[notKeepers,0], com[notKeepers,1], marker='o', ls='None', mec='b', mfc='None')
        pw.show()
    
    # make output data for kept spots
    spotDataList = []
    for iSpot in range(numSpots):
        if not keepers[iSpot]: continue
        index = iSpot+1
        
        spot = getSpot(thisframe, labels, objs, index, keepWithinBBox, False, darkframe=darkframe) # padSpot
        spotDataList.append(spot)
    return labels, objs, spotDataList, bin

class IntensityFunc3D(object):
    """
    This is just a template for intensity distribution functions in 3D
    """
    def __init__(self, *args, **kwargs):
        self.debug = False
        raise NotImplementedError
    def guessXVec(self, x, y, z, w=None, v=None, noBkg=False): 
        raise NotImplementedError
    def getCenter(self, xVec):
        raise NotImplementedError
    def getFWHM(self, xVec):
        raise NotImplementedError
    def getNParams(self, noBkg=False):
        raise NotImplementedError
    def eval(self, xVec, x, y, z, w=None, vSub=None, vScale=None, diff=False, noBkg=False):
        """
        if w is None: x, y, and z are 1D;
        if w has 1D weights: x, y, and z are 2D
        
        if vSub is present, it is subtracted from results -- useful
        for forming least-squares residuals; and vScale is used to scale the results
        if it is present
        """
        raise NotImplementedError
    def deval(self, xVec, x, y, z, w=None, vSub=None, vScale=None, noBkg=False):
        return self.eval(xVec, x, y, z, w=w, vSub=vSub, vScale=vScale, diff=True,  noBkg=noBkg)
    def __call__(self, xVec, x, y, z, w=None, vSub=None, vScale=None, noBkg=False, **kwargs):
        return self.eval(xVec, x, y, z, w=w, vSub=vSub, vScale=vScale, diff=False, noBkg=noBkg, **kwargs)
#
class IntensityFuncGauss3D(IntensityFunc3D):
    """
    8 parameters:
    	centers (3)
        FWHMs   (3)
        scaling  (1)
        background (1)
    
    """
    __nParams = 8
    def __init__(self):
        self.debug = False
        self.expFact  = -4.0 * math.log(2.0)
        return
    def guessXVec(self, x, y, z, w=None, v=None, noBkg=False): 
        nX = self.getNParams(noBkg=noBkg)
        'leave the following check for elsewhere -- useful to be able to guess parameters regardless'
        # assert num.size(x) >= nX,\
        #     'not enough data to reasonably guess at parameters'
        xList = [x, y, z]
        xVec = getGaussNDParams(xList, w=w, v=v)
        if noBkg:
            'last entry from getGaussNDParams is background'
            xVec = xVec[:-1]
        return xVec
    def getCenter(self, xVec):
        return xVec[0:3]
    def getFWHM(self, xVec):
        return xVec[3:6]
    def getNParams(self, noBkg=False):
        retval = self.__nParams
        if noBkg:
            retval -= 1
        return retval
    def constructParams(self, center, fwhm, A, B=0.0):
        xVec = num.zeros(self.getNParams())
        xVec[0:3] = center[:]
        xVec[3:6] = fwhm[:]
        xVec[6]   = A
        xVec[7]   = B
        return xVec
    @classmethod
    def get2DFunc(cls):
        instance = IntensityFuncGauss2D()
        return instance
    def eval(self, xVec, x, y, z, w=None, vSub=None, vScale=None, diff=False, noBkg=False, minWidth=None):
        
        unscaledNorm = None
        
        if w is not None:
             assert x.shape == w.shape,\
                 'w not the same shape as other arguments; %s != %s' % (str(x.shape), str(w.shape))
        
        nX = self.getNParams(noBkg=noBkg)
            
        x0, y0, z0 = xVec[0:3]
        if minWidth is None:
            xW, yW, zW = num.array(xVec[3:6], dtype=float)
            usedMax = num.zeros(3, dtype=bool)
        else:
            xW, yW, zW = num.array(num.maximum(minWidth,xVec[3:6]), dtype=float)
            usedMax = minWidth > xVec[3:6]
        A = xVec[6]
        if noBkg:
            B = 0.0
        else:
            B = xVec[7]
        
        xDist = num.atleast_1d( (x - x0) / xW )
        yDist = num.atleast_1d( (y - y0) / yW )
        zDist = num.atleast_1d( (z - z0) / zW )

        if diff:
            
            inputShape = num.shape(xDist)
            xDist = xDist.flatten()
            yDist = yDist.flatten()
            zDist = zDist.flatten()
            dg_dx = num.zeros( (xDist.size, nX) )

            expEval = num.exp( self.expFact * (
                xDist * xDist + yDist * yDist + zDist * zDist
                ))

            factor = A * expEval * self.expFact * 2.0 
            dg_dx[:,0] = factor * xDist * (-1.0/xW)
            dg_dx[:,1] = factor * yDist * (-1.0/yW)
            dg_dx[:,2] = factor * zDist * (-1.0/zW)
            if not usedMax[0]:
                dg_dx[:,3] = factor * xDist * (-xDist/xW)
            if not usedMax[1]:
                dg_dx[:,4] = factor * yDist * (-yDist/yW)
            if not usedMax[2]:
                dg_dx[:,5] = factor * zDist * (-zDist/zW)
            dg_dx[:,6] = expEval
            if not noBkg:
                dg_dx[:,7] = 1.0 # B
            #
            dg_dx = dg_dx.reshape(num.hstack((inputShape,nX)))

            if w is None:
                retval = dg_dx
            else:
                if len(w.shape) == 2:
                    'retval = num.sum(gauss3d * w, axis=1)'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,:,xInd] * w, axis=1)
                elif len(w.shape) == 1:
                    'retval = gauss3d * w'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,xInd] * w)
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
                
            'nothing needs to be done with vSub'
        
        else:
            
            gauss3d = A * num.exp( self.expFact * (
                xDist * xDist + yDist * yDist + zDist * zDist
                )) + B
            
            if w is None:
                retval = gauss3d
            else:
                # retval = num.dot(gauss3d, w) # no longer the right way to do it
                if len(w.shape) == 2:
                    retval = num.sum(gauss3d * w, axis=1)
                elif len(w.shape) == 1:
                    retval = gauss3d * w
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
            if vSub is not None:
                retval = (retval - vSub)
                unscaledNorm = num.linalg.norm(retval)

        if vScale is not None:
            if len(retval.shape)==1:
                realw = vScale
            else:
                'only happens when diff is True?'
                realw = num.tile(vScale.reshape(len(vScale),1),retval.shape[1])
            retval = retval * realw
    
        if self.debug and vSub is not None: 
            if diff:
                print '%s : did diff eval' % (self.__class__.__name__)
            else:
                print '%s : eval with norm %24.14e unscaled: %24.14e' % (self.__class__.__name__, num.linalg.norm(retval), unscaledNorm)

        return retval

class IntensityFuncGauss3DGenEll(IntensityFunc3D):
    """
    generalization of IntensityFuncGauss3D to have principal axes generally aligned 
    11 parameters:
    	centers (3)
        diagonal "fwhm" (3)
        scaling  (1)
        off-diagonal (3)
        background (1)
    
    """
    __nParams = 11
    def __init__(self):
        self.debug = False
        self.expFact  = -4.0 * math.log(2.0)
        return
    def guessXVec(self, x, y, z, w=None, v=None, noBkg=False): 
        """
        guess is for aligned ellipsoid
        """
        nX = self.getNParams(noBkg=noBkg)
        'leave the following check for elsewhere -- useful to be able to guess parameters regardless'
        # assert num.size(x) >= nX,\
        #     'not enough data to reasonably guess at parameters'
        xList = [x, y, z]
        xVecAligned = getGaussNDParams(xList, w=w, v=v)
        if noBkg:
            'last entry from getGaussNDParams is background'
            xVec = num.hstack((xVecAligned[:-1], num.zeros(3)))
        else:
            xVec = num.hstack((xVecAligned[:-1], num.zeros(3), xVecAligned[-1]))
        'leave the following check for elsewhere -- useful to be able to guess parameters regardless'
        # if num.any(xVec[3:6] <= 0.):
        #    raise UnfitableError(10, 'guess xVec has a zero on the diagonal of the distance measure')
        return xVec
    def getCenter(self, xVec):
        return xVec[0:3]
    def getFWHM(self, xVec):
        return xVec[3:6]
    def getNParams(self, noBkg=False):
        retval = self.__nParams
        if noBkg:
            retval -= 1
        return retval
    @classmethod
    def get2DFunc(cls):
        """
        if drop to 2D, also drop back to aligned ellipsoid
        """
        instance = IntensityFuncGauss2D()
        return instance
    def eval(self, xVec, x, y, z, w=None, vSub=None, vScale=None, diff=False, noBkg=False):
        
        if w is not None:
             assert x.shape == w.shape,\
                 'w not the same shape as other arguments; %s != %s' % (str(x.shape), str(w.shape))
            
        nX = self.getNParams(noBkg=noBkg)
            
        x0, y0, z0 = xVec[0:3]
        A = xVec[6]
        wSvec = num.hstack( (num.array(xVec[3:6],dtype=float), num.array(xVec[7:10],dtype=float)) )
        if noBkg:
            B = 0.0
        else:
            B = xVec[10]
        
        wMatx = tens.svecToSymm( wSvec )
        try:
            dMatx = num.linalg.inv(wMatx)
        except:
            # (etype, eobj, etb) = sys.exc_info()
            raise RuntimeError, 'upon exception from inv, xVec is '+str(xVec)
        
        xDist = num.atleast_1d( (x - x0) ) # / xW
        yDist = num.atleast_1d( (y - y0) ) # / yW
        zDist = num.atleast_1d( (z - z0) ) # / zW
        del x0, y0, z0
        
        if diff:

            dDMatx_dxW_svecList = tens.dAiAoH_svecList(dMatx)
            
            inputShape = num.shape(xDist)
            xDist = xDist.flatten()
            yDist = yDist.flatten()
            zDist = zDist.flatten()
            xyzDist = num.vstack((xDist,yDist,zDist))
            dAll = num.dot( dMatx, xyzDist )
            # dSqr = num.apply_along_axis( lambda x:num.dot(x,x), 0, dAll)
            dSqr = num.sum(dAll * dAll, axis=0)

            dg_dx = num.zeros( (xDist.size, nX) )

            expEval = num.exp( self.expFact * dSqr )
            factor = A * expEval * self.expFact * 2.0 
            
            ddSqr_dx0 = -num.dot(dMatx.T, dAll)
            dg_dx[:,0:3] = num.vstack( (factor * ddSqr_dx0[0,:], factor * ddSqr_dx0[1,:], factor * ddSqr_dx0[2,:]) ).T
            
            # these need to agree with wSvec assembly above!
            dg_dx[:, 3] = factor * num.sum( num.dot(dDMatx_dxW_svecList[0], xyzDist) * dAll, axis=0)
            dg_dx[:, 4] = factor * num.sum( num.dot(dDMatx_dxW_svecList[1], xyzDist) * dAll, axis=0)
            dg_dx[:, 5] = factor * num.sum( num.dot(dDMatx_dxW_svecList[2], xyzDist) * dAll, axis=0)
            dg_dx[:, 6] = expEval
            dg_dx[:, 7] = factor * num.sum( num.dot(dDMatx_dxW_svecList[3], xyzDist) * dAll, axis=0)
            dg_dx[:, 8] = factor * num.sum( num.dot(dDMatx_dxW_svecList[4], xyzDist) * dAll, axis=0)
            dg_dx[:, 9] = factor * num.sum( num.dot(dDMatx_dxW_svecList[5], xyzDist) * dAll, axis=0)
            if not noBkg:
                dg_dx[:,10] = 1.0 # B
            #
            dg_dx = dg_dx.reshape(num.hstack((inputShape,nX)))

            if w is None:
                retval = dg_dx
            else:
                if len(w.shape) == 2:
                    'retval = num.sum(gauss3d * w, axis=1)'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,:,xInd] * w, axis=1)
                elif len(w.shape) == 1:
                    'retval = gauss3d * w'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,xInd] * w)
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
                
            'nothing needs to be done with vSub'
            if vScale is not None:
                for xInd in range(dg_dx.shape[-1]):
                    retval[:,xInd] = retval[:,xInd] * vScale[:]
        
        else:
            
            inputShape = num.shape(xDist)
            xDist = xDist.flatten()
            yDist = yDist.flatten()
            zDist = zDist.flatten()
            
            dAll = num.dot( dMatx, num.vstack((xDist,yDist,zDist)) )
            # dSqr = num.apply_along_axis( lambda x:num.dot(x,x), 0, dAll)
            dSqr = num.sum(dAll * dAll, axis=0)
            
            gauss3d = A * num.exp( self.expFact * dSqr ) + B
            
            gauss3d = gauss3d.reshape(inputShape)
            
            if w is None:
                retval = gauss3d
            else:
                # retval = num.dot(gauss3d, w) # no longer the right way to do it
                if len(w.shape) == 2:
                    retval = num.sum(gauss3d * w, axis=1)
                elif len(w.shape) == 1:
                    retval = gauss3d * w
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
            if vSub is not None:
                retval = retval - vSub
                if debugMulti > 1:
                    print 'norm of retval : %g' % (num.linalg.norm(retval))
            if vScale is not None:
                retval = retval * vScale
            
        return retval

def getWtrShd(inp, threshold, gaussFilterSigma=None, footPrint=None, fpi=5, numPeaks=None):
    "fpi only used if footPrint is None"
    if haveSciKits:
        morphology = scikits.image.morphology
    else:
        msg = 'scikits package has not been loaded and is needed for this function'
        raise NameError(msg)
    
    if gaussFilterSigma is not None:
        inpUse = ndimage.filters.gaussian_filter(inp, gaussFilterSigma)
        # inpUse = ndimage.gaussian_gradient_magnitude(inp, gaussFilterSigma)
    else:
        inpUse = inp
    
    if footPrint is None:
        footprint = num.ones([fpi] * inp.ndim) # num.ones((5,5,5))
    
    localMax = morphology.is_local_maximum(inpUse, footprint=footprint)
    labelLM, numLM = ndimage.label(localMax)
    # objs = ndimage.find_objects(labelLM)
    lmMaxList = []
    for index in range(1,numLM+1):
        lmMax = num.max(inpUse[num.where(labelLM == index)])
        lmMaxList.append(lmMax)
    sortedLM = num.argsort(lmMaxList).tolist()
    if debugMulti:
        print 'local maxima are : %s' % (str(lmMaxList))
    
    if numPeaks is None:
        numPeaks = num.sum(num.array(lmMaxList) >= threshold)
        assert numPeaks > 0,\
            'did not find any peaks over threshold'
    else:
        numLMPeaks = num.sum(num.array(lmMaxList) >= threshold)
        assert numLMPeaks >= numPeaks,\
            'did not find enough isolated peaks'
        if numLMPeaks > numPeaks:
            print 'WARNING : found %d local maxima, using the %d highest' % (numLMPeaks, numPeaks)
    markers = num.zeros(inpUse.shape, dtype=int)
    markersIJK = []
    for iMarker in range(1,numPeaks+1):
        index = sortedLM.pop() + 1
        # whereIndex = num.where(labelLM == index)
        whereIndex = zip(*num.where(labelLM == index))[0]
        markers[whereIndex] = iMarker
        markersIJK.append(whereIndex)
    wtrshd = morphology.watershed(-inpUse, markers)
    
    return wtrshd, markersIJK

    
class IntensityFuncMulti3D(IntensityFunc3D):
    """
    combination of multiple overlapped functions
    """
    def __init__(self, subFuncClass, nSubFuncs, minWidth=None ):
        self.debug = False
        self.funcs = [ subFuncClass() for iSubFunc in range(nSubFuncs) ]
        self.minWidth = copy.deepcopy(minWidth)
        return
    def __repr__(self):
        retval = str(self.__class__.__name__)+"("
        retval += str(self.funcs[0]) + ', ' + str(len(self.funcs))
        retval += ")"
        return retval
    def getNX(self, noBkg=False):
        nX = 0        
        for iFunc, func in enumerate(self.funcs):
            nXThis = func.getNParams(noBkg=True)
            nX += nXThis
        if not noBkg:
            nX += 1
        return nX
    def getBounds(self, noBkg=False):
        
        if self.minWidth is None:
            bounds = None
        else:

            nX = self.getNX(noBkg=noBkg)
            lb = -1e120 * num.ones(nX) # None would be better for some implementations
            ub =  1e120 * num.ones(nX) # None would be better for some implementations
            
            iXLB = 0        
            for iFunc, func in enumerate(self.funcs):
                nXThis = func.getNParams(noBkg=True)
                xVecThis = lb[iXLB:iXLB+nXThis]
                fwhm = func.getFWHM(xVecThis)
                if self.minWidth is not None:
                    'fwhm references xVecThis, so should not need to set explicitly if assign by index'
                    fwhm[:] = self.minWidth[:]
                iXLB += nXThis
            if not noBkg:
                iXLB += 1
                'no bound on background'
                pass
            bounds = zip(lb,ub)
            
        return bounds
    def guessXVecPureNDImage(self, x, y, z, pxlCenterList, w=None, v=None, pxl=None, gaussFilterSigma=None, noBkg=False): 
        """
        for now, rely on specified centers;
        assumes that v can be made into integers without too much loss of precision
        """
        assert pxl is not None,\
            'have only coded for case in which pxl are given'
        assert len(pxlCenterList) == len(self.funcs), \
            'must have len(pxlCenterList) == len(self.funcs)'
        
        if not len(x.shape) == 1:
            raise RuntimeError, 'have not coded case of x.shape != 1 -- too messy'

        shape   = [num.max(iPxL)+1 for iPxL in pxl]
        markers = num.zeros( shape, dtype='int8')
        inp     = num.zeros( shape, dtype='uint16')
        
        if v is None:
            inp[pxl] = 1
        else:
            vMax = v.max()
            inp[pxl] = v.max() - num.fmax(v, 0)
            assert inp.max()-inp.min() > 16,\
                'suspiciously low dynamic range in data'
        #
        if gaussFilterSigma is not None:
            inpUse = ndimage.filters.gaussian_filter(inp, gaussFilterSigma)
            # inpUse = ndimage.gaussian_gradient_magnitude(inp, gaussFilterSigma)
        else:
            inpUse = inp
        #
        'make markers'
        for iObj, pxlCenter in enumerate(pxlCenterList):
            label = iObj+1
            markers[pxlCenter] = label
        # markers.reshape(markers.size)[num.argmin(inp)] = -1 # marker for background
        #
        wtrshd = ndimage.watershed_ift(inpUse, markers)
        #
        objs = ndimage.find_objects(wtrshd)
        #
        'for convenience of indexing; does not work with quadrature'
        cList = [x,y,z,w,v]
        cPXLList = []
        for c in cList:
            if c is None:
                cPXLList.append(c)
            else:
                cPXL = num.zeros( shape, dtype=c.dtype)
                cPXL[pxl] = c
                cPXLList.append(cPXL)
        #
        xVecLists = []
        for iObj, obj, func in zip( range(len(objs)), objs, self.funcs):
            label = iObj+1
            objInd = num.where(wtrshd[obj] == label)
            vals = inp[obj][objInd]
            x_offset, y_offset, z_offset = obj[0].start, obj[1].start, obj[2].start
            'this does not work with quadrature'
            ind = (objInd[0]+x_offset, objInd[1]+y_offset, objInd[2]+z_offset)
            xList = [ num.atleast_1d(cPXL[ind]) for cPXL in cPXLList[0:3] ]
            vThis = None; wThis = None;
            if w is not None:
                wThis = cPXLList[3][ind]
            if v is not None:
                vThis = cPXLList[4][ind]
            nXThis = func.getNParams(noBkg=True)
            'do not check nPx < nXThis -- that ends up being too restrictive given vagaries of the watershed algorithm'
            nPx = num.size(xList[0])
            # if nPx < nXThis:
            #     raise UnfitableError(11, 
            #                          "nPx=%d < nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
            #                          )
            xVecThis = func.guessXVec(*xList, w=wThis, v=vThis, noBkg=True)
            fwhm = func.getFWHM(xVecThis)
            if self.minWidth is not None:
                'fwhm references xVecThis, so should not need to set explicitly if assign by index'
                fwhm0, = num.where(fwhm < self.minWidth)
                fwhm[fwhm0] = 1.05 * self.minWidth[fwhm0]
            fwhm0, = num.where(fwhm <= 0)
            if len(fwhm0) > 0:
                raise UnfitableError(12, 
                                     "got a zero width : nPx=%d, nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % 
                                     (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
                                     )
            xVecLists.append(xVecThis)
        if not noBkg:
            if v is None:
                xVecLists.append(0)
            else:
                xVecLists.append(v.min())
        
        retval = num.hstack( xVecLists )
        return retval
    def guessXVec(self, x, y, z, pxlCenterList, w=None, v=None, pxl=None, gaussFilterSigma=None, noBkg=False, wtrshd=None): 
        if not haveSciKits:
            msg = 'scikits package has not been loaded and is needed for this function'
            raise NameError(msg)

        assert pxl is not None,\
            'have only coded for case in which pxl are given'
        assert len(pxlCenterList) == len(self.funcs), \
            'must have len(pxlCenterList) == len(self.funcs)'
        
        if not len(x.shape) == 1:
            raise RuntimeError, 'have not coded case of x.shape != 1 -- too messy'

        numPeaks = len(pxlCenterList)
        
        if wtrshd is None:
            
            shape   = [num.max(iPxL)+1 for iPxL in pxl]
            inp     = num.zeros( shape, dtype='int16')
            if v is None:
                inp[pxl] = 1
            else:
                inp[pxl] = num.maximum(v,0)
                assert inp.max()-inp.min() > 16,\
                    'suspiciously low dynamic range in data'
            
            wtrshd, markersIJK = getWtrShd(inp, 0, gaussFilterSigma=gaussFilterSigma, numPeaks=numPeaks)
    
        # vWtrShd = wtrshd[(x-x.min(),y-y.min(),z-z.min())]
        vWtrShd = wtrshd[pxl]
        # pwWtrshd = spot.display(vAll=vWtrShd)

        xVecLists = []
        for iMarker, func in zip(range(1,numPeaks+1), self.funcs):
            thesePx, = num.where(vWtrShd == iMarker)
            xList = [ x[thesePx], y[thesePx], z[thesePx] ]
            vThis = None; wThis = None;
            if w is not None:
                wThis = w[thesePx]
            if v is not None:
                vThis = v[thesePx]
            nXThis = func.getNParams(noBkg=True)
            'do not check nPx < nXThis -- that ends up being too restrictive given vagaries of the watershed algorithm'
            nPx = num.size(xList[0])
            # if nPx < nXThis:
            #     raise UnfitableError(11, 
            #                          "nPx=%d < nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
            #                          )
            xVecThis = func.guessXVec(*xList, w=wThis, v=vThis, noBkg=True)            
            fwhm = func.getFWHM(xVecThis)
            if self.minWidth is not None:
                'fwhm references xVecThis, so should not need to set explicitly if assign by index'
                fwhm0, = num.where(fwhm < self.minWidth)
                fwhm[fwhm0] = 1.05 * self.minWidth[fwhm0]
            fwhm0, = num.where(fwhm <= 0)
            if len(fwhm0) > 0:
                raise UnfitableError(12, 
                                     "got a zero width : nPx=%d, nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % 
                                     (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
                                     )
            xVecLists.append(xVecThis)
        if not noBkg:
            if v is None:
                xVecLists.append(0)
            else:
                xVecLists.append(v.min())
        
        retval = num.hstack( xVecLists )
        return retval
    def getCenter(self, xVecAll, iSubSpot):
        iXLB = 0
        retval = None
        for iFunc, func in enumerate(self.funcs):
            nXThis = func.getNParams(noBkg=True)
            if iFunc == iSubSpot:
                xVec = xVecAll[iXLB:iXLB+nXThis]
                retval = func.getCenter(xVec)
                break
            iXLB += nXThis
        assert retval is not None, \
            'failed to set a return value, iSubSpot out of range?'
        return retval
    def getFWHM(self, xVecAll):
        iXLB = 0
        retval = None
        for func in self.funcs:
            nXThis = func.getNParams(noBkg=True)
            xVec = xVecAll[iXLB:iXLB+nXThis]
            fwhmThis = func.getFWHM(xVec)
            if retval is None:
                retval = fwhmThis
            else:
                retval = num.minimum(retval, fwhmThis)
            iXLB += nXThis
        # if not noBkg:
        #     B = xVecAll[iXLB] # = xVecAll[-1]
        #     iXLB += 1
        return retval
    def getNParams(self, noBkg=False):
        retval = 0
        for func in self.funcs:
            retval += func.getNParams(noBkg=True)
        if not noBkg:        
            retval += 1 # B # background
        return retval
    def eval(self, xVecAll, x, y, z, w=None, vSub=None, vScale=None, diff=False, noBkg=False, iFunc=None):
        
        nX = self.getNParams(noBkg=noBkg)
        
        # minWidth = None
        # if self.minWidth is not None:
        #     minWidth = 0.99 * self.minWidth

        if diff:
            assert iFunc is None,\
                'not coded: diff and iFunc'
            
            iXLB = 0
            retval = None
            for func in self.funcs:
                nXThis = func.getNParams(noBkg=True)
                xVec = xVecAll[iXLB:iXLB+nXThis]
                retvalThis = func.eval(xVec, x, y, z, w=w, vSub=None, vScale=vScale, diff=diff, noBkg=True)
                # minWidth=minWidth
                if retval is None:
                    retval = num.empty( (num.shape(retvalThis)[0], nX) )
                retval[:,iXLB:iXLB+nXThis] = retvalThis
                iXLB += nXThis
            if not noBkg:
                B = xVecAll[iXLB] # = xVecAll[-1]
                if vScale is not None:
                    retval[:,iXLB] = vScale # B
                else:
                    retval[:,iXLB] = 1.0 # B
                iXLB += 1
            'nothing needs to be done with vSub'
        else:
            iXLB = 0
            retval = None
            for iF, func in enumerate(self.funcs):
                nXThis = func.getNParams(noBkg=True)
                if iFunc is None or iF == iFunc:
                    xVec = xVecAll[iXLB:iXLB+nXThis]
                    retvalThis = func.eval(xVec, x, y, z, w=w, vSub=None, vScale=None, diff=diff, noBkg=True)
                    # minWidth=minWidth)
                    if retval is None:
                        retval = retvalThis
                    else:
                        retval = retval + retvalThis
                iXLB += nXThis
            assert retval is not None,\
                'retval remained as None'
            if not noBkg:
                B = xVecAll[iXLB] # = xVecAll[-1]
                retval = retval + B
                iXLB += 1
            if vSub is not None:
                retval = retval - vSub
                if debugMulti > 1:
                    print 'norm in eval is '+str(num.linalg.norm(retval))
            if vScale is not None:
                retval = retval * vScale
            
        return retval
    
class IntensityFunc2D(object):
    """
    This is just a template for intensity distribution functions in 2D;
    It could be unified with the 3D version, but that would potentially make the 
    interface harder to understand 
    """
    def __init__(self, *args, **kwargs):
        self.debug = False
        raise NotImplementedError
    def guessXVec(self, x, y, w=None, v=None, noBkg=False):
        raise NotImplementedError
    def getCenter(self, xVec):
        raise NotImplementedError
    def getFWHM(self, xVec):
        raise NotImplementedError
    def getNParams(self, noBkg=False):
        raise NotImplementedError
    def eval(self, xVec, x, y, w=None, vSub=None, vScale=None, diff=False, noBkg=False):
        """
        if w is None: x and y are 1D;
        if w has 1D weights: x and y are 2D
        
        if vSub is present, it is subtracted from results -- useful
        for forming least-squares residuals; and vScale is used to scale the results
        if it is present
        """
        raise NotImplementedError
    def deval(self, xVec, x, y, w=None, vSub=None, vScale=None, noBkg=False):
        return self.eval(xVec, x, y, w=w, vSub=vSub, vScale=vScale, diff=True,  noBkg=noBkg)
    def __call__(self, xVec, x, y, w=None, vSub=None, vScale=None, noBkg=False):
        return self.eval(xVec, x, y, w=w, vSub=vSub, vScale=vScale, diff=False, noBkg=noBkg)
#
class IntensityFuncGauss2D(IntensityFunc2D):
    """
    6 parameters:
    	centers (2)
        FWHMs   (2)
        scaling  (1)
        background (1)
    """
    __nParams = 6
    def __init__(self):
        self.debug = False
        self.expFact  = -4.0 * math.log(2.0)
        return
    def guessXVec(self, x, y, w=None, v=None, noBkg=False):
        nX = self.getNParams(noBkg=noBkg)
        'leave the following check for elsewhere -- useful to be able to guess parameters regardless'
        # assert num.size(x) >= nX,\
        #     'not enough data to reasonably guess at parameters'
        xList = [x, y]
        xVec = getGaussNDParams(xList, w=w, v=v)
        if noBkg:
            'last entry from getGaussNDParams is background'
            xVec = xVec[:-1]
        return xVec
    def getCenter(self, xVec):
        return xVec[0:2]
    def getFWHM(self, xVec):
        return xVec[2:4]
    def getNParams(self, noBkg=False):
        retval = self.__nParams
        if noBkg:
            retval -= 1
        return retval
    def eval(self, xVec, x, y, w=None, vSub=None, vScale=None, diff=False, noBkg=False):
        
        if w is not None:
             assert x.shape == w.shape,\
                 'w not the same shape as other arguments; %s != %s' % (str(x.shape), str(w.shape))
            
        nX = self.getNParams(noBkg=noBkg)

        x0, y0 = xVec[0:2]
        xW, yW = num.array(xVec[2:4],dtype=float)
        A = xVec[4]
        if noBkg:
            B = 0.0
        else:
            B = xVec[5]
        
        xDist = num.atleast_1d( (x - x0) / xW )
        yDist = num.atleast_1d( (y - y0) / yW )

        if diff:
            
            inputShape = num.shape(xDist)
            xDist = xDist.flatten()
            yDist = yDist.flatten()
            dg_dx = num.zeros( (xDist.size, nX) )

            expEval = num.exp( self.expFact * (
                xDist * xDist + yDist * yDist
                ))

            factor = A * expEval * self.expFact * 2.0 
            dg_dx[:,0] = factor * xDist * (-1.0/xW)
            dg_dx[:,1] = factor * yDist * (-1.0/yW)
            dg_dx[:,2] = factor * xDist * (-xDist/xW)
            dg_dx[:,3] = factor * yDist * (-yDist/yW)
            dg_dx[:,4] = expEval
            if not noBkg:
                dg_dx[:,5] = 1.0 # B
            #
            dg_dx = dg_dx.reshape(num.hstack((inputShape,nX)))

            if w is None:
                retval = dg_dx
            else:
                if len(w.shape) == 2:
                    'retval = num.sum(gauss2d * w, axis=1)'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,:,xInd] * w, axis=1)
                elif len(w.shape) == 1:
                    'retval = gauss2d * w'
                    retval = num.empty( (inputShape[0], nX) )
                    for xInd in range(dg_dx.shape[-1]):
                        retval[:,xInd] = num.sum(dg_dx[:,xInd] * w)
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
                
            'nothing needs to be done with vSub'
        
        else:
            
            gauss2d = A * num.exp( self.expFact * (
                xDist * xDist + yDist * yDist
                )) + B
            
            if w is None:
                retval = gauss2d
            else:
                # retval = num.dot(gauss2d, w) # no longer the right way to do it
                if len(w.shape) == 2:
                    retval = num.sum(gauss2d * w, axis=1)
                elif len(w.shape) == 1:
                    retval = gauss2d * w
                else:
                    raise NotImplementedError,\
                        'shape of w has length %d' % (len(shape(w)))
            if vSub is not None:
                retval = retval - vSub
            
        if vScale is not None:
            if len(retval.shape)==1:
                realw = vScale
            else:
                'only happens when diff is True?'
                realw = num.tile(vSscale.reshape(len(vScale),1),retval.shape[1])
            retval = retval * realw
    
        return retval

class IntensityFuncMulti2D(IntensityFunc2D):
    """
    combination of multiple overlapped functions
    """
    def __init__(self, subFuncClass, nSubFuncs):
        self.debug = False
        self.funcs = [ subFuncClass() for iSubFunc in range(nSubFuncs) ]
        return
    def getBounds(self):
        return None
    def guessXVec(self, x, y, pxlCenterList, w=None, v=None, pxl=None, gaussFilterSigma=None, noBkg=False, minWidth=None, wtrshd=None): # no z 2d/3d
        """
        for now, rely on specified centers;
        assumes that v can be made into integers without too much loss of precision
        """
        assert pxl is not None,\
            'have only coded for case in which pxl are given'
        assert len(pxlCenterList) == len(self.funcs), \
            'must have len(pxlCenterList) == len(self.funcs)'
        
        if not len(x.shape) == 1:
            raise RuntimeError, 'have not coded case of x.shape != 1 -- too messy'

        if wtrshd is None:
            shape   = [num.max(iPxL)+1 for iPxL in pxl]
            inp     = num.zeros( shape, dtype='uint16')
            if v is None:
                inp[pxl] = 1
            else:
                vMax = v.max()
                inp[pxl] = v.max() - num.fmax(v, 0)
                assert inp.max()-inp.min() > 16,\
                    'suspiciously low dynamic range in data'
            
            wtrshd, markersIJK = getWtrShd(inp, 0, gaussFilterSigma=gaussFilterSigma, numPeaks=len(pxlCenterList))

        objs = ndimage.find_objects(wtrshd)
        #
        'for convenience of indexing; does not work with quadrature'
        cList = [x,y,w,v] # no z 2d/3d
        cPXLList = []
        for c in cList:
            if c is None:
                cPXLList.append(c)
            else:
                cPXL = num.zeros( shape, dtype=c.dtype)
                cPXL[pxl] = c
                cPXLList.append(cPXL)
        #
        xVecLists = []
        for iObj, obj, func in zip( range(len(objs)), objs, self.funcs):
            label = iObj+1
            objInd = num.where(wtrshd[obj] == label)
            vals = inp[obj][objInd]
            x_offset, y_offset = obj[0].start, obj[1].start # 2d/3d 0:3, no z
            'this does not work with quadrature'
            ind = (objInd[0]+x_offset, objInd[1]+y_offset) # 2d/3d 0:3, no z
            xList = [ cPXL[ind] for cPXL in cPXLList[0:2] ] # 2d/3d 0:3, no z
            vThis = None; wThis = None;
            if w is not None:
                wThis = cPXLList[2][ind] # 2d/3d
            if v is not None:
                vThis = cPXLList[3][ind] # 2d/3d
            nXThis = func.getNParams(noBkg=True)
            'do not check nPx < nXThis -- that ends up being too restrictive given vagaries of the watershed algorithm'
            nPx = num.size(xList[0])
            # if nPx < nXThis:
            #     raise UnfitableError(11, 
            #                          "nPx=%d < nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
            #                          )
            xVecThis = func.guessXVec(*xList, w=wThis, v=vThis, noBkg=True)
            fwhm = func.getFWHM(xVecThis)
            if minWidth is not None:
                'fwhm references xVecThis, so should not need to set explicitly if assign by index'
                fwhm0, = num.where(fwhm < minWidth)
                fwhm[fwhm0] = minWidth[fwhm0]
            fwhm0, = num.where(fwhm <= 0)
            if len(fwhm0) > 0:
                raise UnfitableError(12, 
                                     "got a zero width : nPx=%d, nParams=%d for sub-function; ind : %s; vThis : %s; pxlCenterList : %s" % 
                                     (nPx, nXThis, str(ind), str(vThis), str(pxlCenterList))
                                     )
            xVecLists.append(xVecThis)
        if not noBkg:
            if v is None:
                xVecLists.append(0)
            else:
                xVecLists.append(v.min())
        
        retval = num.hstack( xVecLists )
        return retval
    def getCenter(self, xVecAll, iSubSpot):
        iXLB = 0
        retval = None
        for iFunc, func in enumerate(self.funcs):
            nXThis = func.getNParams(noBkg=True)
            if iFunc == iSubSpot:
                xVec = xVecAll[iXLB:iXLB+nXThis]
                retval = func.getCenter(xVec)
                break
            iXLB += nXThis
        assert retval is not None, \
            'failed to set a return value, iSubSpot out of range?'
        return retval
    def getFWHM(self, xVecAll):
        iXLB = 0
        retval = None
        for func in self.funcs:
            nXThis = func.getNParams(noBkg=True)
            xVec = xVecAll[iXLB:iXLB+nXThis]
            fwhmThis = func.getFWHM(xVec)
            if retval is None:
                retval = fwhmThis
            else:
                retval = num.minimum(retval, fwhmThis)
            iXLB += nXThis
        # if not noBkg:
        #     B = xVecAll[iXLB] # = xVecAll[-1]
        #     iXLB += 1
        return retval
    def getNParams(self, noBkg=False):
        retval = 0
        for func in self.funcs:
            retval += func.getNParams(noBkg=True)
        if not noBkg:        
            retval += 1 # B # background
        return retval
    def eval(self, xVecAll, x, y, w=None, vSub=None, vScale=None, diff=False, noBkg=False, iFunc=None): # no z 2d/3d
        
        nX = self.getNParams(noBkg=noBkg)
        
        if diff:
            assert iFunc is None,\
                'not coded: diff and iFunc'
            
            iXLB = 0
            retval = None
            for func in self.funcs:
                nXThis = func.getNParams(noBkg=True)
                xVec = xVecAll[iXLB:iXLB+nXThis]
                retvalThis = func.eval(xVec, x, y, w=w, vSub=None, vScale=vScale, diff=diff, noBkg=True) # no z 2d/3d
                if retval is None:
                    retval = num.empty( (num.shape(retvalThis)[0], nX) )
                retval[:,iXLB:iXLB+nXThis] = retvalThis
                iXLB += nXThis
            if not noBkg:
                B = xVecAll[iXLB] # = xVecAll[-1]
                if vScale is not None:
                    retval[:,iXLB] = vScale # B
                else:
                    retval[:,iXLB] = 1.0 # B
                iXLB += 1
            'nothing needs to be done with vSub'
        else:
            iXLB = 0
            retval = None
            for iF, func in enumerate(self.funcs):
                nXThis = func.getNParams(noBkg=True)
                if iFunc is None or iF == iFunc:
                    xVec = xVecAll[iXLB:iXLB+nXThis]
                    retvalThis = func.eval(xVec, x, y, w=w, vSub=None, vScale=None, diff=diff, noBkg=True) # no z 2d/3d
                    if retval is None:
                        retval = retvalThis
                    else:
                        retval = retval + retvalThis
                iXLB += nXThis
            assert retval is not None,\
                'retval remained as None'
            if not noBkg:
                B = xVecAll[iXLB] # = xVecAll[-1]
                retval = retval + B
                iXLB += 1
            if vSub is not None:
                retval = retval - vSub
            if vScale is not None:
                retval = retval * vScale
            
        return retval
    
class UnfitableError(exceptions.Exception):
    def __init__(self,err,msg):
        self.err = err
        self.msg = msg
    
class FitFailedError(exceptions.Exception):
    def __init__(self,err,msg):
        self.err = err
        self.msg = msg
    
class Spot(object):
    __dtype = [
        ('xl', num.ndarray), 
        ('yl', num.ndarray),
        ('vbox', num.ndarray), # all values in box # num.ma.core.MaskedArray
        ('v', num.ndarray), # only values in object
        ('x', num.ndarray), # len(x) == len(v)
        ('y', num.ndarray), # len(y) == len(v)
        ('omega', float),
        ('iFrame', int),
        ('darkval',num.ndarray)
        ]
    markAtOmegaLo = -1
    markAtOmegaHi = -2
    markAtBoundary = -3
    atBoundMarks = (markAtOmegaLo, markAtOmegaHi, markAtBoundary)
    validMarks = (markAtOmegaLo, markAtOmegaHi, markAtBoundary)
    nQP_FWHM = 5 # number of quadrature points desired for integrating across FWHM of a spot
    uniformQP = False
    minQ1DInv = 0.125 # = 1./8. # minimum quad point spacing along any given pixel direction
    # num.hstack(self.data['x'][:,0])
    def __init__(self, key, delta_omega, data=None, 
                 omega=None, iFrame=None, 
                 detectorGeom=None):
        """
        key should be unique among all spots in the collection;
        can call with initial data or not
        
        if have tuple from getDataMinimal call, should be able to init directly with that tuple
        """
        
        self.key = key
        self.delta_omega_abs = abs(delta_omega)

        self.data = None
        self.xAll = self.yAll = self.oAll = self.vAll = None
        self.darkAll = None
        self.finalized = False
        self.marks = []
        self.debug = debugDflt
        self.__detectorGeom = None
        
        self.setDetectorGeom(detectorGeom)
        
        if data is not None:
            if isinstance(data, dict):
                assert omega  is not None, 'if have data dictionary, also need omega'
                assert iFrame is not None, 'if have data dictionary, also need iFrame'
                self.append(data, omega, iFrame)
            elif hasattr(data, '__len__'):
                assert omega  is None, 'should not specify omega'
                assert iFrame is None, 'should not spedify iFrame'
                if isinstance(data, num.ndarray): # hasattr(data, 'dtype')
                    'data from another spot'
                    # spotFinder.Spot._Spot__dtype
                    assert data.dtype == self.__dtype, \
                        'data has bad dtype'
                    self.data = copy.deepcopy(data)
                elif len(data) == 4:
                    'data as from getDataMinimal()[2]'
                    self.xAll, self.yAll, self.oAll, self.vAll = data
                    # detectorGeom = self.detectorGeom # not needed here!
                    assert detectorGeom is not None, \
                        'without darkAll being given, need to provide a detectorGeom'
                    self.darkAll = num.ones_like(self.vAll) * detectorGeom.getVDark()
                    self.finalize()
                elif len(data) == 5:
                    'data as from getDataMinimal()[2]'
                    self.xAll, self.yAll, self.oAll, self.vAll, self.darkAll = data
                    self.finalize()
            else:
                raise RuntimeError, 'do not know what to do with data : '+str(data)
        
        self.__split = False
        self.__xyoCOMList = None
        self.__wtrshd = None
        
        self.cleanFit()
        
        return
    def unsetMulti(self):
        self.__split = False
        self.__xyoCOMList = None
        self.__wtrshd = None
        self.cleanFit()
        return
    def setupMulti(self, xyoList, wtrshd=None):
        """
        for now, rely on candidate centers having been provided; probably want to do a fit
        """
        self.__split = True
        self.__xyoCOMList = copy.deepcopy(xyoList)
        self.__wtrshd = wtrshd
        self.cleanFit()
        return
    def nPx(self):
        assert self.finalized, 'called on spot that was not finalized'
        retval = len(self.vAll)
        return retval
    def getDataMinimal(self):
        assert self.finalized, 'called on spot that was not finalized'
        retval = self.key, self.delta_omega_abs, (self.xAll, self.yAll, self.oAll, self.vAll, self.darkAll)
        return retval
    def getShape(self):
        assert self.finalized, 'called on spot that was not finalized'
        return self.shape
    def __len__(self):
        return self.nPx()
    def mark(self, val):
        self.marks.append(val)
        return
    def isMarked(self, marks):
        """return true of spot has been marked as indicated; marks can be an integer or list;
        if marks is a list, returns true if any of marks are true"""
        retval = False
        for mark in num.atleast_1d(marks):
            Spot.validMarks.index(mark) # make sure mark is valid
            try:
                self.marks.index(mark)
                retval = True # return true of any of marks 
                break 
            except:
                'not marks with mark'
        return retval
    def append(self, dataDict, omega, iFrame):
        'dataDict should be like one of the entries in spotDataList coming back from spotFinderSingle'
        if self.finalized:
            raise RuntimeError, 'cannot append to finalized spot'
        if dataDict['atBoundary']:
            self.mark(self.markAtBoundary)
        xl   = dataDict['xl']
        yl   = dataDict['yl']
        vbox = dataDict['vbox'] 
        v    = vbox[xl,yl] 
        if dataDict.has_key('darkbox'):
            darkbox = dataDict['darkbox']
            darkval = darkbox[xl,yl]
            pass
        else:
            raise RuntimeError, 'darkbox now required'
            pass
        toInsert = num.array(
            [( xl, yl, vbox, v, dataDict['x'], dataDict['y'], omega, iFrame, darkval )],
            dtype=self.__dtype)
        if self.data is None:
            self.data = toInsert
        else:
            self.data = num.vstack((self.data, toInsert))
        return
    def merge(self, other):
        """merge in other's data"""
        assert self.data  is not None, 'self has no data'
        assert other.data is not None, 'other has no data'
        assert not other.finalized, 'other is already finalized'
        self.data = num.vstack((self.data, other.data))
        other.data = None
        return
    def isFinalized(self):
        return self.finalized
    def finalize(self, flatten=False, modifyDeltaOmega=False, cullDupl=True):
        """
        Could potentially get rid of self.data once finalize is called, 
        but leave it around just in case it is useful for subsequent operations.
        It might, for example, be needed if we go back and read more intensity data
        in the vicinity of the spot.
        """
        if self.data is None:
            assert \
                self.xAll is not None and \
                self.yAll is not None and \
                self.oAll is not None and \
                self.vAll is not None and \
                self.darkAll is not None, \
                'do not have necessary data'
            if flatten:
                """
                designed for cases in which spot is effectively only
                one omega frame wide and want to just go ahead and
                collapse it to that frame; if modifyDeltaOmega then
                adjust delta_omega for uncertainty calculations and so
                forth
                """
                if modifyDeltaOmega:
                    self.delta_omega_abs = max(self.delta_omega_abs, self.oAll.max() - self.oAll.min())
                vSum = float(num.sum(self.vAll))
                omegaMean = float(num.dot(self.oAll, self.vAll)) / vSum
                pixelDict = {}
                darkDict={}
                for x_p, y_p, v_p, d_p in zip(self.xAll, self.yAll, self.vAll, self.darkAll):
                    key = (x_p, y_p) # (int(x_p), int(y_p))
                    if pixelDict.has_key(key):
                        pixelDict[key] = pixelDict[key] + v_p
                    else:
                        pixelDict[key] = v_p
                    if darkDict.has_key(key):
                        darkDict[key].append(d_p)# = darkDict[key] + d_p
                    else:
                        darkDict[key] = [d_p]
                avgDark = {}
                for key in darkDict.keys():
                    avgDark[key]=num.mean(darkDict[key])
                    pass
                darkDict = avgDark
                
                nPxFlat = len(pixelDict)
                self.xAll    = num.empty(nPxFlat, dtype=self.xAll.dtype)
                self.yAll    = num.empty(nPxFlat, dtype=self.yAll.dtype)
                self.oAll    = num.empty(nPxFlat, dtype=self.oAll.dtype)
                self.vAll    = num.empty(nPxFlat, dtype=self.vAll.dtype)
                self.darkAll = num.empty(nPxFlat, dtype=self.darkAll.dtype)
                iPx = 0
                for key, value in pixelDict.iteritems():
                    self.xAll[iPx] = key[0]
                    self.yAll[iPx] = key[1]
                    self.vAll[iPx] = value
                    self.darkAll[iPx]=darkDict[key]
                    iPx += 1
                self.oAll[:] = omegaMean
        else:
            if flatten:
                raise NotImplementedError, 'flatten with self.data not None'
            if len(self.data) > 1:
                self.vAll = num.hstack(self.data['v'][:,0])
                self.xAll = num.hstack(self.data['x'][:,0])
                self.yAll = num.hstack(self.data['y'][:,0])
                self.darkAll = num.hstack(self.data['darkval'][:,0])
                lengths   = map(len, self.data['v'][:,0])
                oAll = num.empty(0)
                for iOm, om in enumerate(self.data['omega']):
                    oAll = num.hstack((oAll, num.tile(om,lengths[iOm])))
                self.oAll = oAll
            else:
                'do not want vAll to stay a special class instance if it is'
                self.vAll = num.array(self.data['v'][:][0])
                self.darkAll = num.array(self.data['darkval'][:][0])
                #
                self.xAll = self.data['x'][:][0]
                self.yAll = self.data['y'][:][0]
                length    = len(self.vAll)
                self.oAll = num.tile(self.data['omega'][:][0], length)
        if cullDupl:
            'cull duplicates that could have happened due to merging of spots'
            xyo = num.empty( len(self.vAll), dtype=[('x',int),('y',int),('o',float)] )
            xyo['x'], xyo['y'], xyo['o'] = self.xAll[:], self.yAll[:], self.oAll[:]
            try:
                xyoUnique, indices = num.unique...(xyo, return_index=True)
            except TypeError:
                xyoUnique, indices = num.unique1d...(xyo, return_index=True)
            if debugFinalize:
                nDiscard = len(self.vAll) - len(indices)
                if nDiscard > 0:
                    print '   finalize : %d of %d pixels are duplicate, discarding' % (nDiscard, len(self.vAll))
            self.xAll, self.yAll, self.oAll = xyoUnique['x'], xyoUnique['y'], xyoUnique['o']
            self.vAll = self.vAll[indices]
            self.darkAll = self.darkAll[indices]
        self.shape = num.array([
                self.xAll.max()-self.xAll.min()+1, 
                self.yAll.max()-self.yAll.min()+1, 
                len(num.unique(self.oAll)) # len(self.data)
                ])
        assert num.all(self.shape > 0), 'internal error in shape calculation'
        self.xMM = (self.xAll.min(), self.xAll.max())
        self.yMM = (self.yAll.min(), self.yAll.max())
        self.oMM = (self.oAll.min(), self.oAll.max())
        self.finalized = True
        return
    def setDetectorGeom(self, detectorGeom, clobber=False):
        if self.__detectorGeom is not None:
            if detectorGeom != self.detectorGeom and not clobber:
                raise RuntimeError,\
                    'detectorGeom already set'
        self.__detectorGeom = detectorGeom
        if detectorGeom is not None:
            self.__checkDG()
        return
    def getDetectorGeom(self):
        return self.__detectorGeom
    detectorGeom = property(getDetectorGeom, setDetectorGeom, None)
    def __checkDG(self):
        assert self.__detectorGeom is not None, \
            'need self.__detectorGeom to have been set'
        assert self.__detectorGeom.pVec is None, \
            'Spot should not use a detectorGeom with a pVec'
        return
    def getUncertainties(self):
        raise RuntimeError, \
            'getUncertainties deprecated, call angCOM with getUncertainties=True'
        return
    def cleanFit(self):
        """
        clean up data associated with having done fit()
        """
        self.__fitFailed = False
        self.__fitX = None
        self.__fitFunc = None
        self.__fitNDim = None
        self.__fitUncertainties = None
        self.__fitCov_X = None
        self.__fitConf_level = None
        # self.__xyoCOMList = None
        return
    def fitHasFailed(self):
        retval = bool(self.__fitFailed)
        return retval
    def exceptionOnFit(self):
        """
        return the exception if stored on fit;
        note that not all exceptions get stored
        """
        # isEx = isinstance(self.__fitFailed, exceptions.Exception) # storing exception instances seems to cause trouble
        isEx = isinstance(self.__fitFailed, tuple)
        if isEx:
            # retval = self.__fitFailed
            retval = FitFailedError(*self.__fitFailed)
        else:
            retval = False
        return retval
    def fitWrap(self, *args, **kwargs):
        """
        useful if just want to make sure an attempt was made to fit the spot;
        only lets exceptions through if the spot is suspect
        """
        uncertainties = False
        if kwargs.has_key('uncertainties'):
            uncertainties = kwargs['uncertainties']
        if self.__split:
            assert not uncertainties, \
                'uncertainties set True in fitWrap -- case not coded'
        
        try: 
            results = self.fit(*args, **kwargs)
        except UnfitableError, e:
            """
            okay, angCOM will still return a meaningful result;
            fitX and vCalc not necessarily usable
            """
            if debugFit:
                print 'spot unfitable : %d, %s' % (e.err, str(e.msg))
            self.cleanFit()
            self.__fitFailed = True
            if uncertainties:
                'base uncertainies on how far appart the corners of a pixel fall in angular space'
                self.__checkDG()
                self.__fitUncertainties = self.angPixelSize
        except FitFailedError, e:
            self.cleanFit()
            # self.__fitFailed = e
            self.__fitFailed = (e.err, e.msg)
            if uncertainties:
                'not necessarily okay to have such a nasty spot'
                raise
        '''
        # used to catch other stuff here, but then that results in errors getting masked and 
        # looking like failures to fit when they are really errors that need to be fixed
        except:
            (etype, eobj, etb) = sys.exc_info()
            self.cleanFit()
            self.__fitFailed = (1001, str(eobj))
            if uncertainties:
                'not necessarily okay to have such a nasty spot'
                raise
        '''
        
        return
    def fit(self, 
            funcType = None,
            quadr = 'auto',
            full_output = False,
            uncertainties = False,
            confidence_level = 0.95,
            debug=None,
            detectorGeom = None,
            fout=sys.stdout,
            ):
        """
        fit a distribution to a spot;

        may throw a UnfitableError Exception if it is not possible to fit 
        the spot with the given function
        
        if want to re-fit a spot, call the cleanFit method first

        if just want to make sure an attempt was made to fit the spot, call fitWrap
        """
        if debug is None:
            debug = self.debug
        #
        if not self.finalized:
            raise RuntimeError, 'can only fit finalized spots'
        #
        if funcType is None:
            if spotFitFuncDflt == 'auto':
                nSplit = max(1, self.getNSplit())
                'if spot is bigger than about 6**3 pixels per peak, then use gaussGenEll'
                if self.nPx() > 200*nSplit:
                    funcType = 'gaussGenEll'
                else:
                    funcType = 'gauss'
            else:
                funcType = spotFitFuncDflt
        #
        if detectorGeom is not None:
            self.setDetectorGeom(detectorGeom)
        if self.__fitFailed or \
                (self.__fitFunc is not None and \
                     (not uncertainties or self.__fitUncertainties is not None)):
            e = self.exceptionOnFit()
            if e:
                raise e
            if self.debug:
                print 'have alread tried fit, not retrying'
            return
        if self.debug:
            print 'doing spot fit'
        
        #
        if self.__split:
            assert not uncertainties, 'not yet coded with uncertainties'
            fit_results = self.__fitSpotMulti(
                funcType,
                quadr,
                full_output, 
                fout=fout)
        else:
            fit_results = self.__fitSpot(
                funcType,
                quadr,
                full_output or uncertainties, 
                fout=fout)
        if full_output:
            retval = fit_results
        else:
            retval = fit_results[0:2]
        
        fitX = fit_results[0]
        fitV = fit_results[1]
        if num.any(num.isnan(fitV)):
            raise FitFailedError(5, 'fit has nan values')
        'store results'
        self.__fitX    = fitX
        self.__fitV    = fitV
        self.__fitNDim = fit_results[2]
        self.__fitFunc = fit_results[3]
        if uncertainties:
            raise NotImplementedError('uncertainty not implemented')
            params = fit_results[0]
            cov_x = fit_results[4]
            #u_is = uncertainty_analysis.computeAllparameterUncertainties(params, cov_x, confidence_level)
            self.__fitUncertainties = u_is
            self.__fitCov_X         = cov_x
            self.__fitConf_level    = confidence_level
            
        return retval
    def pixelDistXYO(self, x_a, y_a, o_a, x_b, y_b, o_b, asInt=False):
        if asInt:
            retval = \
                num.abs( num.round( x_a - x_b) ) + \
                num.abs( num.round( y_a - y_b) ) + \
                num.abs( num.round((o_a - o_b)/self.delta_omega_abs) )
            retval = num.array(retval, dtype=int)
        else:
            retval = \
                num.abs(  x_a - x_b) + \
                num.abs(  y_a - y_b) + \
                num.abs( (o_a - o_b)/self.delta_omega_abs )
        return retval
    def getPixelDists(self, x, y, o, asInt=False):
        retval = self.pixelDistXYO(
            self.xAll, self.yAll, self.oAll,
            x, y, o
            )
        return retval
    def xyoToIJKLocal(self,
                      x, y, o):
        xMin, yMin, oMin = self.xAll.min(), self.yAll.min(), self.oAll.min()
        k = num.array(num.round((o-oMin) / self.delta_omega_abs), dtype=int)
        retval = (x - xMin, y - yMin, k)
        return retval
    def ijkLocalToXYO(self,
                      i, j, k):
        xMin, yMin, oMin = self.xAll.min(), self.yAll.min(), self.oAll.min()
        retval = (i+xMin, j+yMin, k*self.delta_omega_abs+oMin)
        return retval
    def __fitSpotMulti(self, 
                       funcType,
                       quadr, 
                       full_output,
                       minXYODist=3.,
                       checkForValley=False,
                       fout=sys.stdout,
                       ):
        debug1 = self.debug > 1 or debugMulti
        
        quadrMesgForm = 'quadrature rule %s not available, peak is perhaps too tight for fitting; pixel widths : %s; requested quad points : %s'
        
        self.__checkDG()
        xyoToAng  = self.__detectorGeom.xyoToAngMap
        angToXYO  = self.__detectorGeom.angToXYO
        angPixelSize = self.angPixelSize
        
        xAng, yAng, zAng, w = xyoToAng(self.xAll, self.yAll, self.oAll, outputDV=True)
        
        nSubFuncs = len(self.__xyoCOMList)
        
        withOmega = True
        if self.shape[0] == 1 or self.shape[1] == 1:
            raise UnfitableError(1, 'fitSpotMulti does not support this spot shape '+str(self.shape))
        if self.shape[2] == 1:
            'this can happen, even with spot padding, if flattenOmega was called on the spot'
            withOmega = False
        'check for spots that should be flattened in omega'
        xVecNDP = getGaussNDParams( (xAng, yAng, zAng), w=w, v=self.vAll)
        fwhmOme = xVecNDP[3+2]
        if fwhmOme < angPixelSize[2]/8.:
            if debugMulti:
                print >> fout,  'omega width is estimated at %g compared to pixel size of %g, flattening' % (fwhmOme, angPixelSize[2])
            self.flattenOmega()
            xAng, yAng, zAng, w = xyoToAng(self.xAll, self.yAll, self.oAll, outputDV=True)
            withOmega = False
        
        # minWidth = angPixelSize
        minWidth = (1.1 * self.minQ1DInv * self.nQP_FWHM) * angPixelSize
        
        if withOmega:
            ndim = 3
            if funcType is 'gauss' or funcType is 'normal':
                intensityFuncClass = IntensityFuncGauss3D
            elif funcType is 'gaussGenEll' or funcType is 'normalGenEll':
                intensityFuncClass = IntensityFuncGauss3DGenEll
            else:
                raise NotImplementedError, 'non-gauss intensity function'
            intensityFunc = IntensityFuncMulti3D( intensityFuncClass, nSubFuncs, minWidth=minWidth )
        else:
            ndim = 2
            if funcType is 'gauss' or funcType is 'normal' or funcType is 'gaussGenEll' or funcType is 'normalGenEll':
                intensityFuncClass = IntensityFuncGauss2D
            else:
                raise NotImplementedError, 'non-gauss intensity function'
            intensityFunc = IntensityFuncMulti2D( intensityFuncClass, nSubFuncs, minWidth=minWidth[0:2] )
        
        'absolute frame pixel location does not matter for guessXVec -- just need connectivity'
        
        ijk = self.xyoToIJKLocal(self.xAll, self.yAll, self.oAll)
        ijkCOMList = []
        vCOMList   = []
        for xyoCOM in self.__xyoCOMList: # for angCOM in self.__angCOMList:
            # xyoCOM = map(float, angToXYO(*angCOM))
            'find closest pixel and assign its ijk, so that watershed does not get hosed by a dead pixel'
            'note that ijk are set up so that referenced to 0'
            dists = self.getPixelDists(*xyoCOM)
            iPxMin = num.argmin(dists)
            ijkCOM = self.xyoToIJKLocal(self.xAll[iPxMin], self.yAll[iPxMin], self.oAll[iPxMin])
            ijkCOMList.append(ijkCOM)
            vCOMList.append(self.vAll[iPxMin])
            if debug1:
                print >> fout,  'got ijkCOM at %s with v of %g, compared to mean of %g' % (str(ijkCOM), self.vAll[iPxMin], num.mean(self.vAll))
        """
        check that the peaks are sufficiently distinct;
        if we had fixed integrated intensities or the like then this might not be necessary
        """
        for iCOM, (xyoCOM_i, v_i) in enumerate(zip(self.__xyoCOMList, vCOMList)): # for iCOM, (angCOM_i, v_i) in enumerate(zip(self.__angCOMList, vCOMList)):
            for jCOM, (xyoCOM_j, v_j) in enumerate(zip(self.__xyoCOMList, vCOMList)): # for angCOM_j, v_j in zip(self.__angCOMList, vCOMList)[iCOM+1:]:
                if jCOM <= iCOM:
                    continue
                #angCOM = 0.5*(angCOM_i + angCOM_j)
                #xyoCOM = map(float, angToXYO(*angCOM))
                distIJ = self.pixelDistXYO(*(list(xyoCOM_i)+list(xyoCOM_j)))
                if distIJ < minXYODist:
                    raise UnfitableError(1, 'subspots %d and %d are only %g pixels appart' % (iCOM, jCOM, distIJ))
                if checkForValley:
                    xyoCOM = 0.5*(num.array(xyoCOM_i) + num.array(xyoCOM_j))
                    dists = self.getPixelDists(*xyoCOM)
                    iPxMin = num.argmin(dists)
                    v_mid = self.vAll[iPxMin]
                    if not (v_mid < v_i and v_mid < v_j):
                        raise UnfitableError(1, 'intensities do not have a valley: %g %g %g' % (v_i, v_mid, v_j))

        if not withOmega:
            ijk = ijk[0:2]
            ijkCOMList = [ ijkCOM[0:2] for ijkCOM in ijkCOMList ]
        
        'get xVec0 without quadrature, as guessXVec not yet coded to deal with the extra mess'
        try:
            if withOmega:
                xVec0 = intensityFunc.guessXVec(xAng, yAng, zAng, ijkCOMList, w, self.vAll, ijk, wtrshd=self.__wtrshd)
            else:
                xVec0 = intensityFunc.guessXVec(xAng, yAng,       ijkCOMList, w, self.vAll, ijk, wtrshd=self.__wtrshd)
        except UnfitableError, e:
            if debugMulti :
                print >> fout,  'message from UnfitableError : %d %s' % (e.err, e.msg)
            raise e
        except:
            (etype, eobj, etb) = sys.exc_info()
            if debugMulti :
                print >> fout,  'message from exception : %s' % (eobj.message)
            raise UnfitableError(1, 'guessXVec raised an exception : %s' % (eobj.message))
        if debugMulti:
            print >> fout,  'for func %s guess for xVec is %s' % ( str(intensityFunc), str(xVec0) )
        
        if quadr == 'auto':
            'get guess of xVec without quadrature, so that can guess appropriate quadrature'

            if debug1 : print >> fout,  'spot size and angPixelSize : %d %s ' % (len(self), str(angPixelSize))
            
            xVecWOQP = xVec0 
            fwhm = intensityFunc.getFWHM(xVecWOQP)

            if withOmega:
                nPx_FWHM = fwhm/angPixelSize
            else:
                nPx_FWHM = fwhm/angPixelSize[0:2]
            if debug1 : print >> fout,  'fwhm : '+str(fwhm)
            if debug1 : print >> fout,  'nPx_FWHM : '+str(nPx_FWHM)
            q1D_invAll = nPx_FWHM / self.nQP_FWHM 
            'given that minWidth was specified, probably do not need to check q1D_tooSmall, but do it anyway'
            q1D_tooSmall = q1D_invAll < self.minQ1DInv
            numTooSmall = num.sum(q1D_tooSmall)
            if numTooSmall:
                mesg = 'q1D_invAll : '+str(q1D_invAll)
                print >> fout,  mesg
                raise UnfitableError(1, mesg)
            q1D_all = 1.0 / q1D_invAll
            if self.uniformQP:
                q1D = max(3, int(num.round(num.max(q1D_all))) )
            else:
                q1D = num.maximum(3, num.array(num.round(q1D_all), dtype=int))
            
            if withOmega:
                if self.uniformQP:
                    qRule = '3x%d' % (q1D)
                else:
                    qRule = '%db%db%d' % tuple(q1D.tolist())
                if debug1 : print >> fout,  'from %s requested quad points : %s' % (str(q1D_all), str(qRule))
                xi, w = q3db.qLoc(qRule)
            else:
                if self.uniformQP:
                    qRule = '2x%d' % (q1D)
                else:
                    qRule = '%db%d' % tuple(q1D.tolist())
                if debug1 : print >> fout,  'from %s requested quad points : %s' % (str(q1D_all), str(qRule))
                xi, w = q2db.qLoc(qRule)
                
        else:
            'not automatic quadrature'
            if withOmega:            
                xi, w = q3db.qLoc(quadr)
            else:
                xi, w = q2db.qLoc(quadr)
        
        "intensityFunc is now set"
        if self.nPx() < intensityFunc.getNParams():
            # self.cleanFit()
            raise UnfitableError(1, 'nPx < nParams')
        
        bounds = intensityFunc.getBounds()

        """
        xi are in [0,1], so need to be centered (and scaled for omega)
        
        ... think about what this integral actually means in terms of converting intensity to measured intensity over some area on the detector ... need some geometric factors in the integral to account for angle of incidence or the like? ... would go into dvQP

        """
        xi[:,:] = xi[:,:] - 0.5
        if withOmega:
            xi[:,2] = xi[:,2] * self.delta_omega_abs
        #
        'store information so that can reconstruct quadrature rules'
        self.__setQ(w, xi, withOmega)

        xAng, yAng, zAng, wQP = self.__getQP() # get versions for use with quadrature
        #
        'set of vScale for weighting, desirable in general and for uncertainty analysis'
        vScale = self.__detectorGeom.getVScale(self.vAll + self.darkAll)
        if withOmega:
            args = (xAng, yAng, zAng, wQP, self.vAll, vScale)
        else:
            args = (xAng, yAng,       wQP, self.vAll, vScale)
        #
        if bounds is None:
            maxIter = min(300, len(xVec0)*100) # do not let it take too long if things look nasty
            x, cov_x, infodict, mesg, ier = \
                optimize.leastsq(intensityFunc, # intensityFunc.eval should work too
                                 xVec0,
                                 Dfun=intensityFunc.deval,
                                 args=args,
                                 maxfev=maxIter,
                                 full_output=1)
            message = ' message :  %s' % (mesg)
            if debug1 : print >> fout,  'leastsq did %d function calls' % (infodict['nfev'])
            if not [1,2,3,4].count(ier):
                if len(message) > 0:
                    print >> fout,  'error code %d from leastsq, message : %s' % (ier, message)
                if self.debug: self.display(title='spot %s failed to fit' % (str(self.key)))
                # self.cleanFit()
                #raise UnfitableError(2, 'error from leastsq')
                raise FitFailedError(2, 'error from leastsq')
        else:
            if debug1 : print >> fout,  'have bounds, will use lbfgsb'
            maxIter = min(300, len(xVec0)*100) # do not let it take too long if things look nasty
            cov_x = None # probably need to fix this if want to do uncertainties
            infodict = {}
            'need as scalar objective function and gradient'
            rRef = intensityFunc.eval(xVec0, *args)
            fRefInv = 1.0 / (0.5 * num.dot(rRef, rRef))
            def fThis(xCur):
                resid = intensityFunc.eval(xCur, *args)
                retval = fRefInv * 0.5 * num.dot(resid, resid) # num.linalg.norm(f)
                return retval
            def fPrimeThis(xCur): 
                resid = intensityFunc.eval(xCur, *args)
                jacob = intensityFunc.deval(xCur, *args)
                retval = fRefInv * num.dot(jacob.T, resid)
                return retval
            """
            'having trouble with slsqp, use bfgs instead'
            x, fx, its, imode, smode = optimize.fmin_slsqp(
                    fThis, xVec0, bounds=bounds, fprime=fPrimeThis, 
                    args=(), # args already wrapped up
                    full_output=True, iter=maxIter)
            if debug1 : print >> fout,  'fmin_slsqp did %d iterations' % (its)
            if imode != 0:
                if self.debug: self.display(title='spot %s failed to fit' % (str(self.key)))
                # self.cleanFit()
                #raise UnfitableError(2, 'error from leastsq')
                raise FitFailedError(imode, 'error from fmin_slsqp : '+smode)
            """
            x, fx, infodict = optimize.lbfgsb.fmin_l_bfgs_b(
                fThis, xVec0, bounds=bounds, fprime=fPrimeThis,
                pgtol=1e-5, maxfun=maxIter, disp=0, m=max(len(xVec0),20))
            if debug1 : 
                print >> fout,  'fmin_l_bfgs_b did %d iterations' % (infodict['funcalls'])
                # ...*** report on how many x are near bounds?
            if not infodict['warnflag'] == 0:
                if self.debug: self.display(title='spot %s failed to fit' % (str(self.key)))
                # self.cleanFit()
                #raise UnfitableError(2, 'error from leastsq')
                if infodict['warnflag'] == 2:
                    message = str(infodict['task'])
                else:
                    message = str(infodict['warnflag'])
                raise FitFailedError(1, 'error from fmin_l_bfgs_b : '+message)

        
        'calculate intensities at solution'
        if withOmega:
            vCalc = intensityFunc(x, xAng, yAng, zAng, wQP)
        else:
            vCalc = intensityFunc(x, xAng, yAng,       wQP)
        
        'check again that spots are sufficiently separated'
        xyoCOMList = []
        for iSubSpot in range(self.getNSplit()):
            # have not set self.__fitFunc yet, so cannot call self.xyoCOM
            angCOM = intensityFunc.getCenter(x, iSubSpot)
            if not withOmega:
                angCOM = num.hstack( (angCOM, self.oAll[0]) )
            xyoCOM = self.__detectorGeom.angToXYO(*angCOM)
            if not self.xyoIsInSpot(*xyoCOM):
                raise UnfitableError(1, 'post-fit, subspot %d at %s is not in the spot' % (iSubSpot, str(xyoCOM)) )
            xyoCOMList.append(xyoCOM)
        for iCOM, xyoCOM_i in enumerate(xyoCOMList):
            for jCOM, xyoCOM_j in enumerate(xyoCOMList):
                if jCOM <= iCOM:
                    continue
                distIJ = self.pixelDistXYO(*(list(xyoCOM_i)+list(xyoCOM_j)))
                if distIJ < minXYODist:
                    raise UnfitableError(1, 'post-fit, subspots %d and %d are only %g pixels appart' % (iCOM, jCOM, distIJ))
        
        if full_output:
            assert cov_x is not None, \
                 'singular matrix encountered -- infinite covariance in some direction'
            retval = (x, vCalc, ndim, intensityFunc, cov_x, infodict)
        else:                                      
            retval = (x, vCalc, ndim, intensityFunc)
        return retval
    def __fitSpot(self, 
                  funcType,
                  quadr, 
                  full_output, 
                  fout=sys.stdout):
        
        quadrMesgForm = 'quadrature rule %s not available, peak is perhaps too tight for fitting; pixel widths : %s; requested quad points : %s'
        
        self.__checkDG()
        xyoToAng = self.__detectorGeom.xyoToAngMap

        ndim = 3 - list(self.shape).count(1)
        withOmega = ndim == 3
        if ndim == 3:
            if funcType is 'gauss' or funcType is 'normal':
                intensityFunc = IntensityFuncGauss3D()
            elif funcType is 'gaussGenEll' or funcType is 'normalGenEll':
                intensityFunc = IntensityFuncGauss3DGenEll()
            else:
                raise NotImplementedError, 'non-gauss intensity function'
            withOmega = True
        elif ndim == 2:
            if funcType is 'gauss' or funcType is 'normal':
                intensityFunc = IntensityFuncGauss2D()
            else:
                raise NotImplementedError, 'non-gauss intensity function'
            if self.shape[2] != 1:
                raise NotImplementedError, 'coded only for no spread in omega'
            withOmega = False
        else:
            raise NotImplementedError, 'effective spot dimension is %d' % (ndim)
        
        if quadr == 'auto':
            'get guess of xVec without quadrature, so that can guess appropriate quadrature'
            xAng, yAng, zAng = xyoToAng(self.xAll, self.yAll, self.oAll)
            angPixelSize = self.angPixelSize
            if self.debug >1 : print >> fout,  'angPixelSize : '+str(angPixelSize)
            
            dropFrom3D = False; newFunc = None;
            if withOmega:
                xVecWOQP = intensityFunc.guessXVec(xAng, yAng, zAng, v=self.vAll)
                fwhm = intensityFunc.getFWHM(xVecWOQP)
                nPx_FWHM = fwhm/angPixelSize
                if self.debug >1 : print >> fout,  'fwhm : '+str(fwhm)
                if self.debug >1 : print >> fout,  'nPx_FWHM : '+str(nPx_FWHM)
                q1D_invAll = nPx_FWHM / self.nQP_FWHM 
                q1D_tooSmall = q1D_invAll < self.minQ1DInv
                numTooSmall = num.sum(q1D_tooSmall)
                if numTooSmall:
                    if numTooSmall == 1 and q1D_tooSmall[2] == True:
                        'only too small in omega'
                        dropFrom3D = True
                        newFunc = intensityFunc.get2DFunc()
                    else:
                        mesg = 'q1D_invAll : '+str(q1D_invAll)
                        print >> fout,  mesg
                        raise UnfitableError(1, mesg)
                else:
                    q1D_all = 1.0 / q1D_invAll
                    if self.uniformQP:
                        q1D = max(3, int(num.round(num.max(q1D_all))) )
                        qRule = '3x%d' % (q1D)
                    else:
                        q1D = num.maximum(3, num.array(num.round(q1D_all), dtype=int))
                        qRule = '%db%db%d' % tuple(q1D.tolist())
                    if self.debug >1 : print >> fout,  'from %s requested quad points : %s' % (str(q1D_all), str(qRule))
                    xi, w = q3db.qLoc(qRule)
            
            if dropFrom3D:
                if self.debug > 1:
                    print >> fout,  'dropping from 3D to 2D spot'
                self.flattenOmega()
                assert ndim == 3, 'internal error, ndim'
                assert newFunc is not None, 'internal error, newFunc'
                intensityFunc = newFunc
                withOmega = False
                ndim = 2
                xAng, yAng, zAng = xyoToAng(self.xAll, self.yAll, self.oAll)
                
            if not withOmega:
                xVecWOQP = intensityFunc.guessXVec(xAng, yAng, v=self.vAll)
                fwhm = intensityFunc.getFWHM(xVecWOQP)
                nPx_FWHM = fwhm/angPixelSize[0:2]
                if self.debug >1 : print >> fout,  'fwhm : '+str(fwhm)
                if self.debug >1 : print >> fout,  'nPx_FWHM : '+str(nPx_FWHM)
                iargmin = nPx_FWHM.argmin()
                if nPx_FWHM[iargmin] == 0.:
                    mesg = 'no measurable width in direction '+str(iargmin)
                    print >> fout,  mesg
                    raise UnfitableError(1, mesg)
                q1D_all = self.nQP_FWHM / nPx_FWHM
                if self.uniformQP:
                    q1D = max(3, int(num.round(num.max(q1D_all))) )
                    qRule = '2x%d' % (q1D)
                else:
                    q1D = num.maximum(3, num.array(num.round(q1D_all), dtype=int))
                    qRule = '%db%d' % tuple(q1D.tolist())
                if self.debug >1 : print >> fout,  'from %s requested quad points : %s' % (str(q1D_all), str(qRule))
                try: 
                    xi, w = q2db.qLoc(qRule)
                except:
                    mesg = quadrMesgForm % (str(qRule), str(nPx_FWHM), str(q1D_all)) 
                    print >> fout,  mesg
                    raise UnfitableError(1, mesg)
            if self.debug > 1:
                print >> fout,  'auto quadrature resulted in %d quad points' % (len(w))
        else:
            'not automatic quadrature'
            if withOmega:
                xi, w = q3db.qLoc(quadr)
            else:
                xi, w = q2db.qLoc(quadr)
        
        "intensityFunc is now set"
        if self.nPx() < intensityFunc.getNParams():
            # self.cleanFit()
            raise UnfitableError(1, 'nPx < nParams')

        """
        xi are in [0,1], so need to be centered (and scaled for omega)
        
        ... think about what this integral actually means in terms of converting intensity to measured intensity over some area on the detector ... need some geometric factors in the integral to account for angle of incidence or the like? ... would go into dvQP

        """
        xi[:,:] = xi[:,:] - 0.5
        if withOmega:
            xi[:,2] = xi[:,2] * self.delta_omega_abs
        #
        'store information so that can reconstruct quadrature rules'
        self.__setQ(w, xi, withOmega)
        #
        xAng, yAng, zAng, wQP = self.__getQP()
        #
        'set of vScale for weighting for uncertainty analysis'
        vScale = self.__detectorGeom.getVScale(self.vAll + self.darkAll)
        if withOmega:
            args = (xAng, yAng, zAng, wQP, self.vAll, vScale)
        else:
            args = (xAng, yAng,       wQP, self.vAll, vScale)
        
        # 'do not use existing fit in case it is bad' 
        # com_xyo = self.xyoCOM(useFit=False)
        # com_ang = xyoToAng(*com_xyo)
        # #
        # widths_xyo = self.shape * num.array([1., 1., self.delta_omega_abs])
        # widths_ang = \
        #     num.array(xyoToAng(*(com_xyo + 0.5 * widths_xyo))) - \
        #     num.array(xyoToAng(*(com_xyo - 0.5 * widths_xyo)))
        # if withOmega:
        #     xVec0 = intensityFunc.guessXVec(com_ang, widths_ang, 
        #                                     self.vAll.mean(), self.vAll.min(), self.vAll.max())
        # else:
        #     xVec0 = intensityFunc.guessXVec(com_ang[0:2], widths_ang[0:2], 
        #                                     self.vAll.mean(), self.vAll.min(), self.vAll.max())
        xVec0 = intensityFunc.guessXVec(*args[0:-1]) # guessXVec does not take vScale
        fwhm = intensityFunc.getFWHM(xVec0)
        if num.any( fwhm <= 0. ):
            raise UnfitableError(1, 'guess has zero width in at least one dimension : '+str(fwhm))

        intensityFunc.debug = self.debug

        if full_output:
            x, cov_x, infodict, mesg, ier = \
                optimize.leastsq(intensityFunc, # intensityFunc.eval should work too
                                 xVec0,
                                 Dfun=intensityFunc.deval,
                                 ftol=1e-12,
                                 args=args,
                                 full_output=1)
            message = ' message :  %s' % (mesg)
            assert cov_x is not None, \
                'singular matrix encountered -- infinite covariance in some direction'
        else:
            x, ier = \
                optimize.leastsq(intensityFunc, # intensityFunc.eval should work too
                                 xVec0,
                                 Dfun=intensityFunc.deval,
                                 ftol=1e-12,
                                 args=args)
            message = ''
            if not [1,2,3,4].count(ier):
                print >> fout,  'error code %d from leastsq' % (ier)
        if not [1,2,3,4].count(ier):
            if len(message) > 0:
                print >> fout,  'error code %d from leastsq, message : %s' % (ier, message)
            if self.debug: self.display(title='spot %s failed to fit' % (str(self.key)))
            # self.cleanFit()
            #raise UnfitableError(2, 'error from leastsq')
            raise FitFailedError(2, 'error from leastsq')
        
        'check that center is in the spot'
        # have not set self.__fitFunc yet, so cannot call self.xyoCOM
        angCOM = intensityFunc.getCenter(x)
        if not withOmega:
            angCOM = num.hstack( (angCOM, self.oAll[0]) )
        xyoCOM = self.__detectorGeom.angToXYO(*angCOM)
        if not self.xyoIsInSpot(*xyoCOM):
            raise UnfitableError(1, 'post-fit, center at %s is not in the spot' % (str(xyoCOM)) )
        
        'calculate intensities at solution'
        if withOmega:
            vCalc = intensityFunc(x, xAng, yAng, zAng, wQP)
        else:
            vCalc = intensityFunc(x, xAng, yAng,       wQP)
        
        if full_output:
            retval = (x, vCalc, ndim, intensityFunc, cov_x, infodict)
        else:                                      
            retval = (x, vCalc, ndim, intensityFunc)
        return retval
    def __getQP(self):
        """
        see __fitSpot
        
        do not store results because in other places we do crazy things 
        like change the oAll and do not want to have to worry about keeping these
        quadrature quantities consistent 
        """
        
        xyoToAng = self.__detectorGeom.xyoToAngMap
        
        nQP = len(self.__q_w)
        nPx = len(self.vAll)
        
        xQP = num.tile(self.xAll, (nQP,1)).T + num.tile(self.__q_xi[:,0], (nPx,1))
        yQP = num.tile(self.yAll, (nQP,1)).T + num.tile(self.__q_xi[:,1], (nPx,1))
        if self.__q_withOmega:
            oQP = num.tile(self.oAll, (nQP,1)).T + num.tile(self.__q_xi[:,2], (nPx,1))
        else:
            oQP = num.tile(self.oAll, (nQP,1)).T

        xAng, yAng, zAng, dvQP = xyoToAng(xQP, yQP, oQP, outputDV=True)
        
        #wQP = num.dot(dvQP, num.diag(w)) # correct, but do with sparse instead
        wD = sparse.lil_matrix((nQP,nQP)); wD.setdiag(self.__q_w);
        wQP = num.array(num.matrix(dvQP) * wD.tocsr())
        
        return xAng, yAng, zAng, wQP
    
    def getAngPixelSize(self):
        retval = self.__detectorGeom.getAngPixelSize(
            self.xyoCOM(useFit=False), 
            self.delta_omega_abs,
            )
        return retval
    angPixelSize = property(getAngPixelSize, None, None)

    def __setQ(self, w, xi, withOmega):
        self.__q_w  = w
        self.__q_xi = xi
        self.__q_withOmega = withOmega
        # self.__q_xAng = self.__q_yAng = self.__q_zAng = self.__q_wQP = None
        return
    def setupQuardFromFWHM(self, fwhm, fout=sys.stdout):
        """
        so far, written only for withOmega being True
        """

        if self.__fitFunc is not None:
            raise RuntimeError, 'do not want to clobber existing quadrature setup'
        
        withOmega = True
        
        angPixelSize = self.getAngPixelSize()
        nPx_FWHM = fwhm/angPixelSize
        q1D_invAll = nPx_FWHM / self.nQP_FWHM 
        q1D_tooSmall = q1D_invAll < self.minQ1DInv
        numTooSmall = num.sum(q1D_tooSmall)
        if numTooSmall:
            mesg = 'q1D_invAll : '+str(q1D_invAll)
            print >> fout,  mesg
            raise UnfitableError(1, mesg)
        
        q1D_all = 1.0 / q1D_invAll
        if self.uniformQP:
            q1D = max(3, int(num.round(num.max(q1D_all))) )
            qRule = '3x%d' % (q1D)
        else:
            q1D = num.maximum(3, num.array(num.round(q1D_all), dtype=int))
            qRule = '%db%db%d' % tuple(q1D.tolist())
        if self.debug >1 : print >> fout,  'from %s requested quad points : %s' % (str(q1D_all), str(qRule))
        xi, w = q3db.qLoc(qRule)
        
        xi[:,:] = xi[:,:] - 0.5
        if withOmega:
            xi[:,2] = xi[:,2] * self.delta_omega_abs
        
        self.__setQ(w, xi, withOmega)
        
        return

    def getVCalc(self, intensityFunc, x, noBkg=False, iSubSpot=None):
        xAng, yAng, zAng, wQP = self.__getQP()
        if iSubSpot is not None and not self.__split:
            raise RuntimeError, 'only set iSubSpot for a spot that is split'
        if iSubSpot is not None and self.__split:
            vCalc = intensityFunc(x, xAng, yAng, zAng, wQP, noBkg=noBkg, iFunc=iSubSpot)
        else:
            if self.__q_withOmega:
                vCalc = intensityFunc(x, xAng, yAng, zAng, wQP, noBkg=noBkg)
            else:
                vCalc = intensityFunc(x, xAng, yAng, wQP, noBkg=noBkg)
        return vCalc
    def getIntegratedIntensity(self, useFit=True, iSubSpot=None, getResid=True):
        if self.haveFit() and useFit:
            # assert self.__fitFunc is not None, 'have not done fit, set useFit to False'
            
            x = self.__fitX
            intensityFunc = self.__fitFunc
            vCalcNBG = self.getVCalc(intensityFunc, x, noBkg=True, iSubSpot=iSubSpot)
            integIntens = num.sum(vCalcNBG)
            
            if getResid:
                'also compute the residual for the fit'
                resid = self.__fitV - self.vAll
                res   = num.linalg.norm(resid)
                res0  = num.linalg.norm(self.vAll.astype(float))
                retval = integIntens, res, res0
            else:
                retval = integIntens
        else:
            'without having done fit, do not have a good way of doing background subtraction'
            assert iSubSpot is None,\
                'need useFit and spot.haveFit() if specify subSpot'
            retval = num.sum(self.vAll)
        return retval
    def flattenOmega(self, modifyDeltaOmega=False):
        """
        flatten spot in Omega; 
        
        unless modifyDeltaOmega is set, do not change delta_omega,
        under the assumption that the spot is being flattened if it
        does not significantly spill across multiple omega frames
        """
        assert self.finalized, 'called on spot that was not finalized'
        
        self.dataPreFlatten = self.data
        self.data = None
        self.finalize(flatten=True, modifyDeltaOmega=False)
        
        return
    
    def __xyoCOM(self):
        """estimate of center-of-mass in x, y, omega coordinates on the image stack"""
        assert self.finalized, 'called on spot that was not finalized'
        vSum = float(num.sum(self.vAll))
        com_xyo = num.array([
                float(num.dot(self.xAll, self.vAll)) / vSum,
                float(num.dot(self.yAll, self.vAll)) / vSum,
                float(num.dot(self.oAll, self.vAll)) / vSum
                ])
        return com_xyo
    def __angCOM(self):
        assert self.finalized, 'called on spot that was not finalized'
        self.__checkDG()
        vSum = float(num.sum(self.vAll))
        if vSum <= 0:
            print >> sys.stderr, "vSum %g <= 0, using only positive values" % (vSum)
            vSum = num.sum(self.vAll[self.vAll > 0])
            assert vSum > 0, 'vSum <= 0, all intensities are negative?'
        a1, a2, a3 = self.__detectorGeom.xyoToAngMap(self.xAll, self.yAll, self.oAll)
        com_ang = num.array([
                num.dot(a1, self.vAll) / vSum,
                num.dot(a2, self.vAll) / vSum,
                num.dot(a3, self.vAll) / vSum
                ])
        return com_ang
    def haveFit(self):
        retval = self.__fitFunc is not None
        return retval
    def xyoCOM(self, useFit=True, iSubSpot=None):
        """
        Get center-of-mass in x, y, omega coordinates on the image stack.
        If useFit, then return value from function fit instead of simple estimate
        """
        if self.__fitFunc is not None and useFit:
            self.__checkDG()
            if self.__split:
                assert iSubSpot is not None,\
                    'if specify useFit for split spot then also specify iSubSpot'
                angCOM = self.__fitFunc.getCenter(self.__fitX, iSubSpot)
                if self.__fitNDim == 2:
                    angCOM = num.hstack( (angCOM, self.oAll[0]) )
                # angCOM = self.__fit_angCOMList[iSubSpot]
            else:
                angCOM = self.__fitFunc.getCenter(self.__fitX)
                if self.__fitNDim == 2:
                    angCOM = num.hstack( (angCOM, self.oAll[0]) )
            retval = self.__detectorGeom.angToXYO(*angCOM)
        else:
            assert iSubSpot is None,\
                'iSubSpot specified but could not be used'
            retval = self.__xyoCOM()
        return retval
    def getNSplit(self):
        if self.__split:
            retval = len(self.__xyoCOMList)
        else:
            retval = 0
        return retval
    def angCOM(self, useFit=True, detectorGeom=None, getUncertainties=False, iSubSpot=None):
        """
        Get center-of-mass in twotheta, eta, omega coordinates.
        If useFit, then return value from function fit instead of simple estimate
        
        Could look at detectorGeom to see if it has a pVec and do something
        more elaborate if detectorGeom is already set without a pVec, but
        that is probably asking for trouble
        """
        pi = math.pi

        if self.__split:
            assert not getUncertainties, \
                'getUncertainties for split spot'
            if self.__fitFunc is not None and useFit: # self.__fit_angCOMList is not None:
                assert iSubSpot is not None,\
                    'if specify useFit for split spot then also specify iSubSpot'
                retval = self.__fitFunc.getCenter(self.__fitX, iSubSpot)
                if self.__fitNDim == 2:
                    retval = num.hstack( (retval, self.oAll[0]) )
            else:
                if iSubSpot is None:
                    retval = self.__angCOM()
                else:
                    '__xyoCOMList may just be guesses, do not want to return those numbers'
                    # retval = self.__xyoCOMList[iSubSpot]
                    raise RuntimeError, 'iSubSpot speicified, but either fit was not done successfully or useFit was not specified'
            
        else:
            
            assert iSubSpot is None, \
                'got iSubSpot for spot that is not split'
            
            if detectorGeom is not None:
                self.setDetectorGeom(detectorGeom)
            if self.__fitFunc is not None and useFit:
                retval = self.__fitFunc.getCenter(self.__fitX)
                if self.__fitNDim == 2:
                    retval = num.hstack( (retval, self.oAll[0]) )
                if getUncertainties:
                    if self.__fitUncertainties is None:
                        raise RuntimeError, 'do not have information for uncertainties'
                    unc = self.__fitUncertainties[0:self.__fitNDim]
                    if self.__fitNDim == 2:
                        unc = num.hstack( (unc, self.delta_omega_abs) )
                    # retval = (retval, unc) # do this later, so that can map more easily
            else:
                'want to convert to angular positions and do weighted average there'
                #xyoCOM = self.__xyoCOM()
                #retval = self.__detectorGeom.xyoToAngMap(*xyoCOM)
                retval = self.__angCOM()
                if getUncertainties:
                    if self.__fitUncertainties is not None:
                        assert len(self.__fitUncertainties) == 3, 'bad length in uncertainties'
                        unc = self.__fitUncertainties
                        # retval = (retval, unc) # do this later, so that can map more easily
                    else:
                        raise RuntimeError, 'missing uncertainty information'
            'map to have standard branch cut'
            retval[1::] = mapAngle(retval[1::], units='radians')
            if getUncertainties:
                retval = (retval, unc)
            
        return retval
    def getFrames(self, reader=None):
        """
        if have self.data, do not need reader
        """
        if self.data is None:
            assert not reader is None, 'need a reader if do not have self.data'
            omegas = num.unique(self.oAll)
            retval = map(reader.omegaToFrame, omegas)
        else:
            retval = self.data['iFrame']
        return retval
    def getBBox(self, reader=None):
        if not self.finalized:
            raise RuntimeError, 'needs to be finalized'
        iFrames = list(self.getFrames(reader=reader))
        iFrames.sort() # in place
        bbox = ( 
            self.xMM, 
            self.yMM,
            (num.min(iFrames), num.max(iFrames)),
            )
        return bbox
    def xyoIsInSpot(self, x, y, o, pixelDist=0.5):
        if not self.finalized:
            raise RuntimeError, 'needs to be finalized'
        if self.xMM[0] - x > pixelDist:
            return False
        if x - self.xMM[1] > pixelDist:
            return False
        if self.yMM[0] - y > pixelDist:
            return False
        if y - self.yMM[1] > pixelDist:
            return False
        if (self.oMM[0] - o)/self.delta_omega_abs > pixelDist:
            return False
        if (o - self.oMM[1])/self.delta_omega_abs > pixelDist:
            return False
        return True
    def displayFrames(self,
                      reader, # need a new one?
                      nFramesPad=0,
                      sumImg=num.maximum,
                      **kwargs
                      ):
        """
        display the frame(s) from which a spot came
        """
        localReader = reader.makeNew()
        bbox = self.getBBox(reader=reader)
        iFrames = range(*bbox[2])
        if nFramesPad > 0:
            iFrames = range(iFrames[0]-nFramesPad,iFrames[-1]+1+nFramesPad)
        thisframe = localReader(nskip=iFrames[0], nframes=len(iFrames), sumImg=sumImg)
        pw = reader.display(thisframe) 
        pw.drawBBox(bbox, style='r-')
        return pw
    def displayFlat(self,
                    vAll = None,
                    cmap = None,
                    markersize = 2,
                    **kwargs
                    ):
        """
        Flatten in Omega for display, no loss of x-y resolution
        """
        if vAll is None:
            vAll = self.vAll
        'make vAll into float data so that summing does not overrun int storage'
        vAll = num.asarray(vAll).astype(float)
        
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)

        omegas = num.unique(self.oAll)
        nFrames = len(omegas)
        xbound = (self.xAll.min(), self.xAll.max())
        ybound = (self.yAll.min(), self.yAll.max())
        nX = xbound[1]-xbound[0]+1
        nY = ybound[1]-ybound[0]+1
        if cmap is None:
            cmap = detector.getCMap(False)
        
        # use sparse.coo_matrix to do summing
        
        xNZ = self.xAll-xbound[0]
        yNZ = self.yAll-ybound[0]
        A = sparse.coo_matrix(   (vAll, (xNZ, yNZ)), 
                                 shape=(nX,nY))
        pw = plotWrap.PlotWrap(**kwargs)
        pw(A.todense(), vmin=vmin, vmax=vmax, cmap=cmap)
        Z = num.ones((nX,nY), dtype=int)
        # xNZ, yNZ = A.nonzero() # allows for badness when have done keepWithinBBox
        Z[xNZ,yNZ] = 0 
        B = sparse.coo_matrix(Z)
        pw.a.spy(B, markersize=markersize)
        return pw
    def display(self, 
                cmap = None,
                relfigsize=(1,1),
                vAll=None,
                markersize = 2,
                xyoPointsList = [([],{})],
                **kwargs
                ):
        """
        vAll, if present, is used instead of self.vAll;
        useful for results from fitting
        """
        
        angCOMList = []
        comMec = 'r'
        comSty = 'r-'
        #
        'check for uncertainty information'
        haveUnc = True
        try:
            angCOM, angCOM_unc = self.angCOM(getUncertainties=True)
            angCOMList = [(angCOM, angCOM_unc)]
        except:
            angCOMList = []
            haveUnc = False
        #
        'check for split spot with fitted centers'
        if self.getNSplit() > 0:
            angCOMList = []
            if  self.__fitFunc is not None:
                for iSubSpot in range(self.getNSplit()):
                    angCOMList.append(self.angCOM(useFit=True, iSubSpot=iSubSpot))
            else:
                comMec = 'c'
                comSty = 'c-'
                self.__checkDG()
                for xyoCOM in self.__xyoCOMList:
                    angCOM = map(float, self.__detectorGeom.xyoToAng(*xyoCOM))
                    angCOMList.append(angCOM)
        elif len(angCOMList) == 0:
            angCOMList = [self.angCOM(getUncertainties=False)]
        
        if vAll is None:
            vAll = self.vAll
        omegas = num.unique(self.oAll)
        nFrames = len(omegas)
        xbound = (self.xAll.min(), self.xAll.max())
        ybound = (self.yAll.min(), self.yAll.max())
        vmin   = vAll.min()
        vmax   = num.abs(vAll.max())
        centered = detector.getCentered(vmin,vmax)
        if centered:
            vmin = -vmax
        if vmin > 0:
            vmin = 0
        nX = xbound[1]-xbound[0]+1
        nY = ybound[1]-ybound[0]+1
        if cmap is None:
            cmap = detector.getCMap(centered)
        win = plotWrap.PlotWin(numRows=nFrames, numCols=-1, relfigsize=relfigsize, **kwargs)
        for iFrame, omega in enumerate(omegas):
            these = num.where(self.oAll == omega)
            
            xNZ = self.xAll[these]-xbound[0]
            yNZ = self.yAll[these]-ybound[0]
            A = sparse.coo_matrix(   (vAll[these], (xNZ, yNZ)), 
                                     shape=(nX,nY))
            pw = plotWrap.PlotWrap(window=win, showByDefault=False) # xbound=xbound, ybound=ybound
            pw(A.todense(), vmin=vmin, vmax=vmax, cmap=cmap)
            Z = num.ones((nX,nY), dtype=int)
            # xNZ, yNZ = A.nonzero() # allows for badness when have done keepWithinBBox
            Z[xNZ,yNZ] = 0 
            B = sparse.coo_matrix(Z)
            pw.a.spy(B, markersize=markersize)
            
            'take points as xyo instead of ang, as may want to use a detectorGeom with pVec set'
            for xyoPoints, pwKWArgs in xyoPointsList: # for angPoint, pwKWArgs in angPointList:
                for xyoPoint in xyoPoints:
                    if abs(omega - xyoPoint[2]) <= 0.5*self.delta_omega_abs: # if abs(omega - angPoint[2]) < self.delta_omega_abs:
                        # self.__checkDG()
                        # xyo_com = map(float, self.__detectorGeom.angToXYO(*angPoint))
                        # xyo_com = xyo_com - num.array([xbound[0], ybound[0], 0.])
                        xyo_com = xyoPoint - num.array([xbound[0], ybound[0], 0.])
                        if (xyo_com[0] > -0.5 and xyo_com[0] < nX-0.5) and (xyo_com[1] > -0.5 and xyo_com[1] < nY-0.5):
                            'xyoPoint is visible'
                            thesePwKWArgs = {}
                            thesePwKWArgs.update(pwKWArgs)
                            thesePwKWArgs.setdefault('marker','o')
                            thesePwKWArgs.setdefault('ls','None') 
                            thesePwKWArgs.setdefault('mec','b') 
                            thesePwKWArgs.setdefault('mfc','None')
                            pw( [xyo_com[0]], [xyo_com[1]], **thesePwKWArgs )
            for angCOM in angCOMList:
                if haveUnc:
                    angCOM, angCOM_unc = angCOM
                else:
                    'just use pixel size to have something to plot'
                    self.__checkDG()
                    angCOM_unc = self.angPixelSize
                if abs(omega - angCOM[2]) <= 0.5*self.delta_omega_abs:
                    """
                    plot on this omega frame;
                    may eventually want to somehow show distance from angCOM relative to uncertainty;
                    assume have detectorGeom
                    
                    use uncertainty in angular fit to approximate uncertainty in xyo space
                    """
                    self.__checkDG()
                    xyo_com = map(float, self.__detectorGeom.angToXYO(*angCOM))
                    xyo_min = copy.deepcopy(xyo_com)
                    xyo_max = copy.deepcopy(xyo_com)
                    angCOM_p = num.empty(3)
                    
                    for iDim in range(3):
                        for mult in (1.0, -1.0):
                            angCOM_p[:] = angCOM[:]
                            angCOM_p[iDim] = angCOM[iDim] + mult * angCOM_unc[iDim]
                            xyo_p = self.__detectorGeom.angToXYO(*angCOM_p)
                            xyo_min = map(min, xyo_min, xyo_p)
                            xyo_max = map(max, xyo_max, xyo_p)
                    xyo_com = xyo_com - num.array([xbound[0], ybound[0], 0.])
                    xyo_min = xyo_min - num.array([xbound[0], ybound[0], 0.])
                    xyo_max = xyo_max - num.array([xbound[0], ybound[0], 0.])
                    
                    pw( [xyo_com[0]], [xyo_com[1]], marker='o', ls='None', mec=comMec, mfc='None')
                    pw( [xyo_min[0], xyo_max[0], xyo_max[0], xyo_min[0], xyo_min[0]], 
                        [xyo_min[1], xyo_min[1], xyo_max[1], xyo_max[1], xyo_min[1]],
                        style=comSty ) 
            
        win.show()
        return win
    
    @staticmethod
    def cullSpots(spots, tests):
        """
        Apply a list of tests to spots and return a list with spots that fail culled
        
        spots is a list of Spot instances;
        tests is a list of tests to be applied;
        	For convenience, some tests have been defined as static methods off of the
                Spot class.
                If a test is tests has a length, then the non-first entries are 
                used as arguments.
        """
        testFuncs = []
        for testEntry in tests:
            if hasattr(testEntry,'__len__'):
                testFuncs.append( (testEntry[0], testEntry[1:]) )
            else:
                testFuncs.append( (testEntry, []) )
        testResults = num.ones(len(spots), dtype=bool)
        for iSpot, spot in enumerate(spots):
            for test, args in testFuncs:
                if not test(spot, *args):
                    testResults[iSpot] = False
                    break # no need to do other tests on this spot
        return num.array(spots)[testResults]
    @staticmethod
    def testSpotIsNotStreak(spot, width=1):
        assert spot.finalized, 'called on spot that was not finalized'
        testVal = not (spot.shape[0] <= width or spot.shape[1] <= width)
        return testVal
    @staticmethod
    def testSpotNetIntensity(spot, threshold=0):
        assert spot.finalized, 'called on spot that was not finalized'
        vSum = float(num.sum(spot.vAll))
        testVal = vSum > threshold
        return testVal
    @staticmethod
    def testSpotLenAtLeast(spot, minPx):
        assert spot.finalized, 'called on spot that was not finalized'
        testVal = spot.nPx() >= minPx
        return testVal
    @staticmethod
    def storeSpots(f, spotList, closeLocal=True):
        """
        store spots to f, which can be a filename or 
        something which behaves like a shelve instance
        
        stores a somewhat minimal set of data
        
        all spots must have the same detectorGeom
        
        consider changing to something like:
        import cPickle as pickle
        s = open(f,'w')
        pickle.dump(stuff, f)
        f.close()
        """
        if isinstance(f, str):
            shelf = shelve.open(f,'c',protocol=2)
        elif hasattr(f,'keys'):
            shelf = f
            closeLocal = False
        else:
            raise RuntimeError, 'do not know what to do with f : '+str(type(f))
        
        storedDG = False
        if shelf.has_key('nSpots'):
            print 'WARNING: working with existing shelf'
        if shelf.has_key('detectorGeom'):
            del shelf['detectorGeom'] # just in case
        shelf['nSpots'] = len(spotList)
        for iSpot, thisSpot in enumerate(spotList):
            shelf[str(iSpot)] = thisSpot.getDataMinimal()
            if thisSpot.detectorGeom is not None:
                if storedDG:
                    assert storedDG is thisSpot.detectorGeom,\
                        'detectorGeom is not the same for spots in spotList'
                else:
                    storedDG = thisSpot.detectorGeom
                    'make potentially lighter weight instance to store'
                    shelf['detectorGeom'] = thisSpot.detectorGeom.__class__(thisSpot.detectorGeom)
        if closeLocal:
            shelf.close()
            retval = None
        else:
            retval = shel
        return retval
    @staticmethod
    def loadSpots(f, minimal=True, closeLocal=True):
        """
        load spots stored with storeSpots
        
        note that is a detectorGoem was stored by storeSpots, then 
        all spots loaded here will have it attached to them
        """
        if isinstance(f, str):
            shelf = shelve.open(f)
        elif hasattr(f,'keys'):
            shelf = f
            closeLocal = False
        else:
            raise RuntimeError, 'do not know what to do with f : '+str(type(f))
        
        nSpots = shelf['nSpots'] 
        detectorGeom = None
        if shelf.has_key('detectorGeom'):
            detectorGeom = shelf['detectorGeom']
        spotKWArgs = {'detectorGeom':detectorGeom}
        spotList = []
        for iSpot in range(nSpots):
            spotDataMinimal = shelf[str(iSpot)]
            spot = Spot(*spotDataMinimal, **spotKWArgs)
            spotList.append(spot)
        
        if closeLocal:
            shelf.close()
            retval = spotList, detectorGeom
        else:
            retval = spotList, detectorGeom, shelf
        return retval

spots_MP = None
args_fitWrap_MP = None
kwargs_fitWrap_MP = None
debug_MP = None
thresMulti_MP = None
#
def doSpotThis(iSpot):
    "meant for use with multiprocessing"
    global spots_MP
    global args_fitWrap_MP
    global kwargs_fitWrap_MP
    global debug_MP
    global thresMulti_MP
    spots      = spots_MP
    args       = args_fitWrap_MP
    kwargs     = kwargs_fitWrap_MP
    debug      = debug_MP
    thresMulti = thresMulti_MP
    
    spot = spots[iSpot]

    retStr = StringIO.StringIO()
    
    kwargsThis = {}
    kwargsThis.update(kwargs)
    kwargsThis['fout'] = retStr
        
    newSpotsKeysThis = []
    newSpotsAngsThis = []
    if isinstance(spot, tuple):
        raise RuntimeError, 'should only be looping over raw spots, spot %d is %s' % (iSpot, str(spot))
    else:
        
        'figure out if the spot should be split'
        pxl      = spot.xyoToIJKLocal(spot.xAll, spot.yAll, spot.oAll) 
        shape    = [num.max(iPxL)+1 for iPxL in pxl]
        inp      = num.zeros( shape, dtype='int16')
        inp[pxl] = num.maximum(spot.vAll,0)
        try:
            wtrshd, markersIJK = getWtrShd(inp, thresMulti)
        except:
            wtrshd = None; markersIJK = [];
        splitThisSpot = len(markersIJK) > 1
        #
        if splitThisSpot:
            markersXYO = [spot.ijkLocalToXYO(*markerIJK) for markerIJK in markersIJK]
            spot.setupMulti(markersXYO, wtrshd=wtrshd) # does cleanFit
            
            tic = time.time()
            spot.fitWrap(*args, **kwargsThis)
            toc = time.time(); dt = toc - tic
            print >> retStr,     'fitWrap call for multi-spot %d took %g seconds' % (iSpot, dt)
            if debug > 1: 
                print >> sys.stderr, 'fitWrap call for multi-spot %d took %g seconds' % (iSpot, dt)
            if spot.fitHasFailed():
                if debug:
                    e = spot.exceptionOnFit()
                    if e:
                        print >> retStr,  'fit e-failed for multi-peak spot %d len %d : %d, %s' % (iSpot, len(spot), e.err, str(e.msg))
                    else:
                        print >> retStr,  'fit failed for multi-peak spot %d len %d' % (iSpot, len(spot))
                """
                in getHKLSpots this spot will no longer come up because it is 
                split but the fit has failed, so upon reclaiming no grain should grab it
                """
            else:
                # if debug > 1:
                integI, res, res0 = spot.getIntegratedIntensity(useFit=True, getResid=True)
                print >> retStr,  'fit successful for multiply-claimed spot %d ; integI, res, res0 : %g %g %g' % (iSpot, integI, res, res0)
                """
                could reassign claims in here, but that seems messy;
                better to let a grain do that if it wants so that it can keep its data structures consistent;
                and probably want to do resetClaims first
                """
                for iSubSpot in range(spot.getNSplit()):
                    newSpotsKeysThis.append( (iSpot, iSubSpot) )
                    newSpotsAngsThis.append( spot.angCOM(iSubSpot=iSubSpot) )
                # self.__claimedBy[iSpot] = self.__claimInvalid 

        else: # not splitting this spot
            tic = time.time()
            spot.fitWrap(*args, **kwargsThis)
            toc = time.time(); dt = toc - tic
            print >> retStr,     'fitWrap call for spot %d took %g seconds' % (iSpot, dt)
            if debug > 1: 
                print >> sys.stderr, 'fitWrap call for spot %d took %g seconds' % (iSpot, dt)
            'fitWrap will not necessarily raise an exception, but if it stored an exception re-raise it here'
            # self stuff done elsewhere!
            # self.__spotAngCoords[iSpot] = spot.angCOM(useFit=True)
            # self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *self.__spotAngCoords[iSpot] )
    'returning the data from spot might be much cleaner, but is a pain'
    return spot, newSpotsKeysThis, newSpotsAngsThis, retStr.getvalue()
        
class Spots(object):
    """
    Meant for holding spot data in a form useful for indexing
    """
    __keepMarks = False # careful about making this True -- can cause unexpected results, mostly it is for debugging
    # __claimInvalid = -1
    def __init__(self, planeData, data, detectorGeom, omegaMM, **kwargs):
        """
        planeData : an instance or list (for multiple phases) of the PlaneData class
        data : can any of:
        	spots : list of Spot instances, all of which must have been finalized
                spotAngCoords : spot positions in angular coordinates
                None
        detectorGeom : an instance of DetectorGeomGE or the like
        omegaMM : (min, max) of omega range, or list of such min/max tuples
        """
        
        self.__spots = None
        self.__nRawSpots = None
        self.__anyAreSplit = False
        self.debug = debugDflt

        assert hasattr(detectorGeom, 'xyoToAng'), \
            'detectorGeom, of type %s, does not seem to be the right type' % \
            (str(type(detectorGeom)))
        
        'make data into emtpy list of spots if it is None'
        if data is None:
            data = num.zeros((0,3))
            
        if hasattr(data, 'shape') and len(data.shape) == 2:
            self.__initFromAngCoords(planeData, data, detectorGeom, omegaMM, **kwargs)
        else:
            self.__initFromSpotList( planeData, data, detectorGeom, omegaMM, **kwargs)
        return
    
    def __initFromSpotList(self, planeData, spots, detectorGeom, omegaMM, **kwargs):
        
        self.__spots = num.array(spots)
        self.__nRawSpots = len(self.__spots)
        self.detectorGeom = detectorGeom
        
        'sanity checks'
        finalized = num.empty(len(spots))
        for iSpot in range(len(spots)):
            finalized[iSpot] = spots[iSpot].isFinalized()
        if not num.all(finalized):
            raise RuntimeError, 'all spots must have been finalized'

        'get angular coordinates of all spots'
        self.__spotAngCoords = num.empty([len(spots),3])
        for iSpot, spot in enumerate(spots):
            'to get angular coordinates, spots will need detector geometry'
            self.__spotAngCoords[iSpot] = spot.angCOM(detectorGeom=detectorGeom)
        
        self.__initBase(planeData, omegaMM, **kwargs)
        return
    def __initFromAngCoords(self, planeData, spotAngCoords, detectorGeom, omegaMM, **kwargs):
        
        self.detectorGeom = detectorGeom
        self.__spotAngCoords = num.array(spotAngCoords)
        assert self.__spotAngCoords.shape[1] == 3, \
            'spotsAngCorrds have wrong shape : '+str(self.__spotAngCoords.shape)
        
        self.__initBase(planeData, omegaMM, **kwargs)
        return
    def __initBase(self, planeData, omegaMM, **kwargs):
        
        if hasattr(omegaMM[0], '__len__'):
            'assume omegaMM is list of tuples already'
            omegaMM = omegaMM
        else:
            omegaMM = [omegaMM]
        mins = toFloat( zip(*omegaMM)[0], 'radians' )
        maxs = toFloat( zip(*omegaMM)[1], 'radians' )
        self.omegaMM = zip(mins, maxs)
        
        self.planeDataList = num.atleast_1d(planeData)
        
        self.__multiprocMode = False
        
        friedelToler = 0.01
        strainMag = None
        lparmsList = None
        if kwargs.has_key('friedelToler'):
            friedelToler = kwargs.pop('friedelToler')
        if kwargs.has_key('strainMag'):
            strainMag = kwargs.pop('strainMag')
        if kwargs.has_key('lparms'):
            lparmsList = num.atleast_1d(kwargs.pop('lparms'))
            assert len(lparmsList) == len(self.planeDataList), \
                'lparms wrong length'
        self.friedelToler = friedelToler
        self.strainMag = strainMag
        self.lparmsList = lparmsList
        
        if len(kwargs) > 0:
            raise RuntimeError, 'have unparsed keyword arguments with keys: ' + str(kwargs.keys())

        #self.__invalidated = None
        self.__makeDataStructures()
        
        self.resetClaims()
        
        return
        raise NotImplementedError, 'still working on implementation'
    
    def __makeDataStructures(self):
        """
        note that this gets called from multiple places
        """
        if self.__multiprocMode:
            raise RuntimeError, 'do not call splitConflicts when in multiprocMode'
        
        nSpots = len(self.__spotAngCoords)
        
        # if self.__invalidated is None:
        #     self.__invalidated = num.zeros( nSpots, dtype=bool)
        # else:
        #     oldLen = num.size( self.__invalidated )
        #     self.__invalidated.resize( nSpots )
        #     self.__invalidated[oldLen:] = False
        
        """
        get xyo coordinates of all spots
        these should get updated if spotAngCoords are updated
        """
        xyo = self.detectorGeom.angToXYO(
            self.__spotAngCoords[:,0],
            self.__spotAngCoords[:,1],
            self.__spotAngCoords[:,2] )
        
        self.__spotXYOCoords = num.array(xyo).T

        """lump into (non-unique) 2-theta groups and set up other data structures"""
        self.nTThAssoc = num.zeros( nSpots, dtype=int)
        #
        self.planeDataDict = {}
        self.hklIDListAll = []
        hklIDWithPhase = len(self.planeDataList) > 1
        for iPDE, planeDataEntry in enumerate(self.planeDataList):
            lparms = None
            if self.lparmsList is not None:
                lparms = self.lparmsList[iPDE]
            phaseID    = planeDataEntry.getPhaseID() # may return None
            tThRanges  = planeDataEntry.getTThRanges(strainMag=self.strainMag, lparms=lparms)
            hkls       = planeDataEntry.getHKLs()
            hklList    = hkls.tolist()
            dHKLInv    = dict([[tuple(hkl),iHKL] for iHKL, hkl in enumerate(hklList)])
            if phaseID is None:
                phaseID = iPDE
            if hklIDWithPhase:
                hklIDList = [ tuple([phaseID]+hkl) for hkl in hklList]
            else:
                hklIDList = [ tuple(hkl) for hkl in hklList]
            #tThAssoc = num.zeros( [len(tThRanges), nSpots], dtype='i1')
            tThAssoc = num.empty( [len(tThRanges), nSpots], dtype=bool)
            for iHKL, tThRange in enumerate(tThRanges):
                tThLo, tThHi = tThRange
                whichSpots = num.logical_and(
                    self.__spotAngCoords[:,0] > tThLo, # ... convention
                    self.__spotAngCoords[:,0] < tThHi  # ... convention
                    ) 
                # whichSpots[self.__invalidated] = False
                self.nTThAssoc += whichSpots
                # theseSpots = num.where( whichSpots ) # indices
                # theseAngCoords = self.__spotAngCoords[ theseSpots ]
                tThAssoc[iHKL, :] = whichSpots
                # ... # perhaps with azimuthal binning for faster spot finding
            d = {}
            d['phaseID']    = phaseID
            d['planeData']  = planeDataEntry
            d['tThRanges']  = tThRanges
            d['tThAssoc']   = tThAssoc
            d['phaseAssoc'] = num.any(tThAssoc, axis=0)
            d['hklInv']     = dHKLInv
            d['hklList']    = hklList
            d['hklIDList']  = hklIDList
            d['multip']     = planeDataEntry.getMultiplicity()
            self.planeDataDict[phaseID] = d
            self.hklIDListAll += hklIDList
        'now have spots associated to 2-theta bins' 
        nUnassoc    = num.sum(self.nTThAssoc == 0)
        nMultiAssoc = num.sum(self.nTThAssoc >  1)
        print 'after 2-theta binning, out of %d spots, %d are unassociated and %d are multiply-associated' % \
            (nSpots, nUnassoc, nMultiAssoc)
        
        self.friedelPair = num.empty( nSpots, dtype=int)
        self.friedelPair[:] = -1
        """
        maybe do:
        self.friedelPair[iSpot] = jSpot
        self.friedelPair[jSpot] = iSpot

        ... # do friedel pair associations, using friedelToler
        ... # and reduce to friedel pairs? based on a flag
        ... # get JVB to put in something like what is done in XFormSpots.m?

        """
        
        return
    def __len__(self):
        nSpots = len(self.__spotAngCoords)
        return nSpots
    
    # multiprocessing stuff
    def getMPM(self):
        return self.__multiprocMode
    def setMPM(self, value):
        
        if not value:
            'cleanup, could test self.__multiprocMode to see if multiprocessing mode was actually on'
            if not self.__keepMarks: 
                self.__marks = None
        else:
            assert xrdBase.haveMultiProc, \
                'do not have multiprocessing module available'
            
            'put marks in shared memory, and mostly use marks instead of claimedBy'
            self.__marks = multiprocessing.Array(ctypes.c_int, len(self))
            # self.__marks[:] = 0 # not necessary
            # NOTE: that can work on slices of marks, but not more general lists of integers as indices
            
        self.__multiprocMode = value
            
        return
    multiprocMode = property(getMPM, setMPM, None)
    
    def resetDetectorGeom(self, detectorGeom, doFits=False, fitArgs=[], fitKWArgs={}):
        """
        update detector geometry; 
        also fits all spots;
        probably only want to call this for a single-grain data set or the like;
        updates angular and cartesian coordinates, but not 2-theta associations
        or other meta-data
        """
        self.detectorGeom = detectorGeom
        if self.__spots is None:
            raise NotImplementedError, 'not sure what behavior would be desirable here'
        else:
            for iSpot, spot in enumerate(self.__spots):
                if not isinstance(spot, tuple):
                    spot.setDetectorGeom(detectorGeom, clobber=True)
                    if doFits:
                        if spot.fitHasFailed(): # spot.exceptionOnFit()
                            if self.debug:
                                print 'fit failed previously, will not try to fit again'
                        else:
                            spot.cleanFit()
                            spot.fitWrap(*fitArgs, **fitKWArgs)
                    else:
                        spot.cleanFit()
                    angCOM = spot.angCOM()
                    self.__spotAngCoords[iSpot] = angCOM
                    self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *angCOM )
                else:
                    'this is a subspot, and should come after the master spots'
                    iSpotMaster, iSubSpot = spot
                    spotMaster = self.__spots[iSpotMaster]
                    angCOM = spotMaster.angCOM(subSpot=iSubSpot) 
                    self.__spotAngCoords[iSpot] = angCOM
                    self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *angCOM )
        return
    def getAngCoords(self, indices=None):
        """
        ... make returned thing unwritable even if it is a slice
        """
        if indices is None:
            retval = self.__spotAngCoords
        else:
            retval = self.__spotAngCoords[indices]
        return retval 
    def getXYOCoords(self, indices):
        """
        ... make returned thing unwritable even if it is a slice
        """
        return self.__spotXYOCoords[indices]
    def getSpotObj(self, iSpot):
        if self.__spots is None:
            raise RuntimeError, 'do not have spot objects'
        spot = self.__spots[iSpot]
        if isinstance(spot, tuple):
            iSpotMaster, iSubSpot = spot
            spotMaster = self.__spots[iSpotMaster]
            spot = spotMaster
        return spot
    def getPixelIsInSpots(self, indices, xyo, pixelDist=0):
        """
        find spots that are within pixelDist of containing a given pixel;
        indices may be a single integer, a list of spot indices, or a boolean array;
        
        may want to use pixelDist of 1 or 2 to protect against dead pixels
        """

        single = False
        if isinstance(indices,int):
            single = True
            indices = [indices]
        if indices is None:
            indices = slice(0,len(self.__spots)) # range(len(self.__spots))
        
        if self.__spots is None:
            raise RuntimeError, 'needs Spots to have been made from spots'
        
        spots = self.__spots[indices] # should work whether indices is integers or booleans
        retval = num.ones(len(spots), dtype=bool)
        #
        'if have any spots that were split, then do not allow the parent spot to be a hit'
        disallowedList = []
        for iInd, spot in enumerate(spots):
            if isinstance(spot, tuple):
                iSpotMaster, iSubSpot = spot
                disallowedList.append(self.__spots[iSpotMaster])
        #
        for iInd, spot in enumerate(spots):
            for disallowedSpot in disallowedList:
                if spot is disallowedSpot:
                    retval[iInd] = False
                    break
            if retval[iInd] == False:
                continue
            if isinstance(spot, tuple):
                'split spot, assume came from existing fit of a sub-spot or the like'
                iSpotMaster, iSubSpot = spot
                spot = self.__spots[iSpotMaster]
            if spot.xyoIsInSpot(*xyo, pixelDist=pixelDist):
                dist = spot.getPixelDists(*xyo, asInt=True)
                retval[iInd] = num.any(dist <= pixelDist)
            else:
                retval[iInd] = False
            
        # if single: 
        #     retval = retval[0]

        return retval
    def getIntegratedIntensity(self, index, **kwargs):
        """
        wrapper that takes care of subspot details
        """
        assert self.__spots is not None,\
            'no spots data for getting integrated intensities'
        spot = self.__spots[index]
        iSubSpot = None
        if isinstance(spot, tuple):
            iSpotMaster, iSubSpot = spot
            spotMaster = self.__spots[iSpotMaster]
            spotUse = spotMaster
        else:
            spotUse = spot
        retval = spotUse.getIntegratedIntensity(iSubSpot=iSubSpot, **kwargs)
        return retval
    def fitHasFailed(self, index, subSpotOnly=False):
        "if subSpotOnly the only return true if it was a subspot that failed to fit"
        if self.__spots is None:
            retval = False
        else:
            spot = self.getSpotObj(index)
            if subSpotOnly and spot.getNSplit() > 1:
                retval = False
            else:
                retval = spot.fitHasFailed()
        return retval
    def fitSpots(self, indices, *args, **kwargs):
        """
        fit spots with given indices and update spotAngCoords;
        for indices is None, all spots are considered for fitting;
        indices may be a single integer instead of a list of spot indices
        
        set claimsBased if want to base method on spots that have been claimed
        """
        
        if self.__multiprocMode:
            raise RuntimeError, 'do not fit spots when in multiprocessing mode'
        
        uncertainties = kwargs.get('uncertainties', False)
        debug         = kwargs.get('debug', False)
        #
        claimsBased = False
        if kwargs.has_key('claimsBased'):
            claimsBased = kwargs.pop('claimsBased')
        reRaise = not claimsBased
        #
        interactive = False
        if kwargs.has_key('interactive'):
            interactive = kwargs.pop('interactive')
        
        if uncertainties and claimsBased:
            raise RuntimeError, 'have not coded uncertainties and claimsBased'

        single = False
        if isinstance(indices,int):
            single = True
            indices = [indices]
        if claimsBased:
            'must do all self.__spots'
            assert self.__spots is not None,\
                'no way to resolve multiply claimed spots without self.__spots'
            assert indices is None,\
                'specify indices of None if doing claimsBased'
            fitNonMulti = False
            if kwargs.has_key('fitNonMulti'):
                fitNonMulti = kwargs.pop('fitNonMulti')
            indices = range(self.__nRawSpots)
            # self.__spots.resize( self.__nRawSpots ) # not needed here, and want to be able to index back to master spots from current claimedBy entries
            # self.__spotAngCoords.resize( self.__nRawSpots, 3 ) # not needed here
            self.cleanMulti()
        elif indices is None:
            indices = range(len(self))

        #assert self.__spots is not None,\
        #    'called fitSpots but instance does not have spots'
        if self.__spots is None:
            assert uncertainties is False,\
                'cannot do uncertainties for spots without self.__spots'
            retval = self.__spotAngCoords[indices]
            
        else:        
            'okay, have spots to work with'
            
            retval = []
            if claimsBased:
                newSpotsKeys = []
                newSpotsAngs = []
            if debug: 
                print 'working on fitting spots ... '
            for iSpot in indices:
                spot = self.__spots[iSpot]
                if isinstance(spot, tuple):
                    raise RuntimeError, 'should only be looping over raw spots, spot %d is %s' % (iSpot, str(spot))
                    # 'split spot, assume came from existing fit of a sub-spot or the like'
                    # if debug > 1:
                    #     print 'spot %d is a subspot' % (iSpot)
                    # iSpotMaster, iSubSpot = spot
                    # spotMaster = self.__spots[iSpotMaster]
                    # spotAng = spotMaster.angCOM(subSpot=iSubSpot)
                else:
                    
                    doFit = True
                    
                    if claimsBased:
                        claimedBy = self.__claimedBy[iSpot]
                        if claimedBy is None:
                            if debug > 2:
                                print 'spot %d is not claimed' % (iSpot)
                            doFit  = False
                        elif isinstance(claimedBy,list): 
                            doFit  = False # done here, not below
                            'assume claiming objects are grains'
                            if debug > 2:
                                print 'for spot %d, will do setupMulti' % (iSpot)
                            predXYOList = []
                            for grain in claimedBy:
                                # hitReflId = num.where(grain.grainSpots['iRefl'] == iSpot)[0][0] # does not get split-spot claims right
                                hitReflId = -1
                                for iReflId, iRefl in enumerate(grain.grainSpots['iRefl']):
                                    """
                                    for now, assume a grain points to a given spot only once
                                    """
                                    if iRefl == iSpot:
                                        hitReflId = iReflId
                                        break
                                    if iRefl >= self.__nRawSpots:
                                        "iRefl points to a subspot, check the master"
                                        subSpot = self.__spots[iRefl]
                                        assert isinstance(subSpot, tuple), \
                                            'expecting subspot, but has type '+str(type(subSpot))
                                        iSpotMaster, iSubSpot = subSpot
                                        if iSpotMaster == iSpot:
                                            hitReflId = iReflId
                                            break
                                assert hitReflId >= 0, \
                                    'did not find a hit for grain %s and spot %d' % (str(grain.lineageList), iSpot)
                                predAngs = grain.grainSpots['predAngles'][hitReflId]
                                predXYO  = map(float, grain.detectorGeom.angToXYO(*predAngs))
                                predXYOList.append(predXYO)
                            if debug > 1:
                                print 'for spot %d, doing setupMulti with %s' % (iSpot, predXYOList)
                            spot.setupMulti(predXYOList) # does cleanFit
                            """
                            do not call spot.cleanFit; 
                            might want to switch to doing this if it is known that more
                            information is available, like if the intensities of the 
                            subspots are now locked down and we are going to retry
                            """
                            tic = time.time()
                            spot.fitWrap(*args, **kwargs)
                            toc = time.time(); dt = toc - tic
                            if debug > 1: 
                                print 'fitWrap call for multi-spot %d took %g seconds' % (iSpot, dt)
                            if spot.fitHasFailed():
                                if debug:
                                    e = spot.exceptionOnFit()
                                    if e:
                                        print 'fit e-failed for multi-peak spot %d len %d : %d, %s' % (iSpot, len(spot), e.err, str(e.msg))
                                    else:
                                        print 'fit failed for multi-peak spot %d len %d' % (iSpot, len(spot))
                                """
                                in getHKLSpots this spot will no longer come up because it is 
                                split but the fit has failed, so upon reclaiming no grain should grab it
                                """
                            else:
                                if debug > 1:
                                    integI, res, res0 = spot.getIntegratedIntensity(useFit=True, getResid=True)
                                    print 'fit successful for multiply-claimed spot %d ; integI, res, res0 : %g %g %g' % (iSpot, integI, res, res0)
                                """
                                could reassign claims in here, but that seems messy;
                                better to let a grain do that if it wants so that it can keep its data structures consistent;
                                and probably want to do resetClaims first
                                """
                                if interactive:
                                    doLoop = True
                                    while doLoop:
                                        userInput = raw_input('options: (d)isplay, (c)ontinue, (i)nvalidate')
                                        if   userInput == 'd':
                                            spot.display()
                                        elif userInput == 'c':
                                            doLoop = False
                                        elif userInput == 'i':
                                            spot.cleanFit()
                                            # spot.__fitFailed = FitFailedError(1, 'error from  user interaction')
                                            spot.__fitFailed = (1, 'error from  user interaction')
                                            doLoop = False
                                        else:
                                            print 'invalid choice'
                                for iSubSpot in range(spot.getNSplit()):
                                    newSpotsKeys.append( (iSpot, iSubSpot) )
                                    newSpotsAngs.append( spot.angCOM(iSubSpot=iSubSpot) )
                                # self.__claimedBy[iSpot] = self.__claimInvalid 
                        else:
                            doFit  = fitNonMulti
                            if doFit and debug > 1:
                                print 'for spot %d, will do normal fit' % (iSpot)
                    # else:
                    #     useFit = not self.__invalidated[iSpot]
                    #     doFit  = useFit
                        
                    if doFit:
                        tic = time.time()
                        spot.fitWrap(*args, **kwargs)
                        toc = time.time(); dt = toc - tic
                        if debug > 1: 
                            print 'fitWrap call for spot %d took %g seconds' % (iSpot, dt)
                        'fitWrap will not necessarily raise an exception, but if it stored an exception re-raise it here'
                        if reRaise:
                            e = spot.exceptionOnFit()
                            if e: # isinstance(spot.__fitFailed, exceptions.Exception):
                                print 'for spot %d re-raising excetiton from fitWrap' % (iSpot)
                                raise e
                    self.__spotAngCoords[iSpot] = spot.angCOM(useFit=doFit)
                    self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *self.__spotAngCoords[iSpot] )
                    spotAng = spot.angCOM(useFit=doFit, getUncertainties=uncertainties)
                retval.append(spotAng)
            
            if debug: 
                print '... done fitting spots'
            if single:
                retval = spotAng
            else:
                a = num.array(retval)
                if uncertainties:
                    retval = a[:,0], a[:,1]
                else:
                    retval = a
        
        if claimsBased:
            if debug: 
                print 'making new spots data structures'
            nNew = len(newSpotsKeys)
            """
            drop any previously split spots;
            """
            self.__anyAreSplit = False
            if nNew > 0:
                self.__spots.resize( self.__nRawSpots + nNew )
                self.__spotAngCoords.resize( self.__nRawSpots + nNew, 3 )
                for iNew, key, angs in zip( range(self.__nRawSpots, self.__nRawSpots+nNew), newSpotsKeys, newSpotsAngs ):
                    self.__spots[iNew] = key; self.__anyAreSplit = True
                    self.__spotAngCoords[iNew] = angs
                self.__makeDataStructures()
            self.resetClaims()

        return retval 
    def __fitSpotsMulti_serial(self, indices, threshold, *args, **kwargs):
        """
        like fitSpots, but consider splitting up spots with multiple peaks
        """
        
        if self.__multiprocMode:
            raise RuntimeError, 'do not fit spots when in multiprocessing mode'
        
        uncertainties = kwargs.get('uncertainties', False)
        debug         = kwargs.get('debug', False)
        #
        interactive = False
        if kwargs.has_key('interactive'):
            interactive = kwargs.pop('interactive')
        #
        if uncertainties:
            raise RuntimeError, 'have not coded uncertainties'

        'must do all self.__spots'
        assert self.__spots is not None,\
            'no way to resolve multiply claimed spots without self.__spots'
        assert indices is None,\
            'specify indices of None -- other cases not yet coded'
        
        indices = range(self.__nRawSpots)
        # self.__spots.resize( self.__nRawSpots ) # not needed here, and want to be able to index back to master spots from current claimedBy entries
        # self.__spotAngCoords.resize( self.__nRawSpots, 3 ) # not needed here
        self.cleanMulti()

        newSpotsKeys = []
        newSpotsAngs = []
        if debug: 
            print 'working on fitting spots ... '
        for iSpot in indices:
            spot = self.__spots[iSpot]
            if isinstance(spot, tuple):
                raise RuntimeError, 'should only be looping over raw spots, spot %d is %s' % (iSpot, str(spot))
            else:
                
                'figure out if the spot should be split'
                pxl      = spot.xyoToIJKLocal(spot.xAll, spot.yAll, spot.oAll) 
                shape    = [num.max(iPxL)+1 for iPxL in pxl]
                inp      = num.zeros( shape, dtype='int16')
                inp[pxl] = num.maximum(spot.vAll,0)
                try:
                    wtrshd, markersIJK = getWtrShd(inp, threshold)
                except:
                    wtrshd = None; markersIJK = [];
                splitThisSpot = len(markersIJK) > 1
                #
                if splitThisSpot:
                    markersXYO = [spot.ijkLocalToXYO(*markerIJK) for markerIJK in markersIJK]
                    spot.setupMulti(markersXYO, wtrshd=wtrshd) # does cleanFit
                    
                    tic = time.time()
                    spot.fitWrap(*args, **kwargs)
                    toc = time.time(); dt = toc - tic
                    if debug > 1: 
                        print 'fitWrap call for multi-spot %d took %g seconds' % (iSpot, dt)
                    if spot.fitHasFailed():
                        if debug:
                            e = spot.exceptionOnFit()
                            if e:
                                print 'fit e-failed for multi-peak spot %d len %d : %d, %s' % (iSpot, len(spot), e.err, str(e.msg))
                            else:
                                print 'fit failed for multi-peak spot %d len %d' % (iSpot, len(spot))
                        """
                        in getHKLSpots this spot will no longer come up because it is 
                        split but the fit has failed, so upon reclaiming no grain should grab it
                        """
                    else:
                        if debug > 1:
                            integI, res, res0 = spot.getIntegratedIntensity(useFit=True, getResid=True)
                            print 'fit successful for multiply-claimed spot %d ; integI, res, res0 : %g %g %g' % (iSpot, integI, res, res0)
                        """
                        could reassign claims in here, but that seems messy;
                        better to let a grain do that if it wants so that it can keep its data structures consistent;
                        and probably want to do resetClaims first
                        """
                        if interactive:
                            doLoop = True
                            while doLoop:
                                userInput = raw_input('options: (d)isplay, (c)ontinue, (i)nvalidate')
                                if   userInput == 'd':
                                    spot.display()
                                elif userInput == 'c':
                                    doLoop = False
                                elif userInput == 'i':
                                    spot.cleanFit()
                                    # spot.__fitFailed = FitFailedError(1, 'error from  user interaction')
                                    spot.__fitFailed = (1, 'error from  user interaction')
                                    doLoop = False
                                else:
                                    print 'invalid choice'
                        for iSubSpot in range(spot.getNSplit()):
                            newSpotsKeys.append( (iSpot, iSubSpot) )
                            newSpotsAngs.append( spot.angCOM(iSubSpot=iSubSpot) )
                        # self.__claimedBy[iSpot] = self.__claimInvalid 

                else: # not splitting this spot
                    tic = time.time()
                    spot.fitWrap(*args, **kwargs)
                    toc = time.time(); dt = toc - tic
                    if debug > 1: 
                        print 'fitWrap call for spot %d took %g seconds' % (iSpot, dt)
                    'fitWrap will not necessarily raise an exception, but if it stored an exception re-raise it here'
                    self.__spotAngCoords[iSpot] = spot.angCOM(useFit=True)
                    self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *self.__spotAngCoords[iSpot] )
        
        if debug: 
            print '... done fitting spots'
        
        if debug: 
            print 'making new spots data structures'
        nNew = len(newSpotsKeys)
        """
        drop any previously split spots;
        """
        self.__anyAreSplit = False
        if nNew > 0:
            self.__spots.resize( self.__nRawSpots + nNew )
            self.__spotAngCoords.resize( self.__nRawSpots + nNew, 3 )
            for iNew, key, angs in zip( range(self.__nRawSpots, self.__nRawSpots+nNew), newSpotsKeys, newSpotsAngs ):
                self.__spots[iNew] = key; self.__anyAreSplit = True
                self.__spotAngCoords[iNew] = angs
            self.__makeDataStructures()
        self.resetClaims()

        return 
    def fitSpotsMulti(self, indices, threshold, *args, **kwargs):
        """
        like fitSpots, but consider splitting up spots with multiple peaks
        """
        multiprocessing = xrdBase.multiprocessing # former import
        
        if self.__multiprocMode:
            raise RuntimeError, 'do not fit spots when in multiprocessing mode'
        
        uncertainties = kwargs.get('uncertainties', False)
        debug         = kwargs.get('debug', False)
        #
        multiProcMode = True
        if kwargs.has_key('multiProcMode'):
            multiProcMode = kwargs.pop('multiProcMode')
        #
        nCPUs = None
        if kwargs.has_key('nCPUs'):
            nCPUs = kwargs.pop('nCPUs')
        #
        if uncertainties:
            raise RuntimeError, 'have not coded uncertainties'

        'must do all self.__spots'
        assert self.__spots is not None,\
            'no way to resolve multiply claimed spots without self.__spots'
        assert indices is None,\
            'specify indices of None -- other cases not yet coded'
        
        indices = range(self.__nRawSpots)
        # self.__spots.resize( self.__nRawSpots ) # not needed here, and want to be able to index back to master spots from current claimedBy entries
        # self.__spotAngCoords.resize( self.__nRawSpots, 3 ) # not needed here
        self.cleanMulti()

        if debug: 
            print 'working on fitting spots ... '
        #
        global spots_MP 
        global args_fitWrap_MP 
        global kwargs_fitWrap_MP 
        global debug_MP 
        global thresMulti_MP
        spots_MP = self.__spots
        args_fitWrap_MP = args
        kwargs_fitWrap_MP = kwargs
        debug_MP = debug
        thresMulti_MP = threshold
        #
        if multiProcMode:
            multiprocessing = xrdBase.multiprocessing # former import
            nCPUs = nCPUs or xrdBase.dfltNCPU
            pool = multiprocessing.Pool(nCPUs)
            if debug : 
                print "using multiprocessing with %d cpus" % (nCPUs)
            resultsFit = pool.map(doSpotThis, indices, chunksize=10)
        else:
            resultsFit = map(doSpotThis, indices)
        #
        spots_MP = None
        args_fitWrap_MP = None
        kwargs_fitWrap_MP = None
        debug_MP = None
        thresMulti_MP = None
        if multiProcMode:
            print 'about to close pool'
            pool.close()
        #
        if debug: 
            print '... done fitting spots'
        #
        newSpotsKeys = []
        newSpotsAngs = []
        for iSpot, (spot, newSpotsKeysThis, newSpotsAngsThis, retStr) in zip(indices, resultsFit):
            print retStr
            newSpotsKeys = newSpotsKeys + newSpotsKeysThis
            newSpotsAngs = newSpotsAngs + newSpotsAngsThis
            self.__spots[iSpot] = spot # replace old spot!
            if not spot.getNSplit() > 1:
                self.__spotAngCoords[iSpot] = spot.angCOM(useFit=True)
                self.__spotXYOCoords[iSpot] = self.detectorGeom.angToXYO( *self.__spotAngCoords[iSpot] )
        
        if debug: 
            print 'making new spots data structures'
        nNew = len(newSpotsKeys)
        """
        drop any previously split spots;
        """
        self.__anyAreSplit = False
        if nNew > 0:
            self.__spots.resize( self.__nRawSpots + nNew )
            self.__spotAngCoords.resize( self.__nRawSpots + nNew, 3 )
            for iNew, key, angs in zip( range(self.__nRawSpots, self.__nRawSpots+nNew), newSpotsKeys, newSpotsAngs ):
                self.__spots[iNew] = key; self.__anyAreSplit = True
                self.__spotAngCoords[iNew] = angs
            self.__makeDataStructures()
        self.resetClaims()

        return 
    def getPlaneData(self, phaseID=None):
        if phaseID is None:
            assert len(self.planeDataDict) == 1, 'need phaseID if multiphase'
            phaseID = self.planeDataDict.keys()[0]
        planeDataThisPhase = self.planeDataDict[phaseID]
        return planeDataThisPhase['planeData']
    def getHKLSpots(self, hkl, phaseID=None, unclaimedOnly=False, disallowMasterWhenSplit=True):
        """
        get boolean array that is True for spots associated with a
        given hkl (and optionally phase)
        """
        if phaseID is None:
            assert len(self.planeDataDict) == 1, 'need phaseID if multiphase'
            phaseID = self.planeDataDict.keys()[0]
        planeDataThisPhase = self.planeDataDict[phaseID]
        if hasattr(hkl,'__len__'): 
            iHKL = planeDataThisPhase['hklInv'][tuple(hkl)]
        elif isinstance(hkl, int):
            iHKL = hkl
        else:
            raise RuntimeError, 'bad hkl : '+str(hkl)
        #theseSpots = num.where( planeDataThisPhase['tThAssoc'][iHKL, :] )
        
        theseSpots = planeDataThisPhase['tThAssoc'][iHKL, :]
        
        'take care of this in grain by calling spots.fitHasFailed'
        # 'kick out split spots that fail to fit'
        # if self.__spots is None:
        #     for iSpot in range(len(spots)):
        #         if theseSpots[iSpot]:
        #             spot = self.__spots[iSpot]
        #             if spot.__split and spot.fitHasFailed():
        #                 theseSpots[iSpotMaster] = False
        
        if self.__spots is not None and disallowMasterWhenSplit and self.__anyAreSplit:
            for iSpot in range(len(self.__spots)):
                if theseSpots[iSpot]:
                    spot = self.__spots[iSpot]
                    if isinstance(spot, tuple):
                        'this spot is allowed, so do not allow the master'
                        iSpotMaster, iSubSpot = spot
                        theseSpots[iSpotMaster] = False

        
        if unclaimedOnly:
            notClaimed = -num.array(self.__claimedBy, dtype=bool)
            theseSpots = theseSpots & notClaimed
            if self.__marks is not None:
                whereThese = num.where(theseSpots)[0]
                'theseSpots[whereThese] = self.__marks[whereThese] <= 0'
                for iSpot in whereThese:
                    theseSpots[iSpot] = self.__marks[iSpot] <= 0
        
        return theseSpots
    def getIterHKL(self, hkl, phaseID=None, unclaimedOnly=True, friedelOnly=False, iSpotLo=0, returnBothCoordTypes=False):
        """
        Returns an iterator over a given set of spots
        """
        iterator = SpotsIterator(self, hkl, phaseID, unclaimedOnly, friedelOnly, iSpotLo, returnBothCoordTypes)
        return iterator

    def getIterPhase(self, phaseID=None, unclaimedOnly=True, friedelOnly=False, iSpotLo=0, returnBothCoordTypes=False):
        """
        Returns an iterator over a given set of spots
        """
        iterator = SpotsIterator(self, None, phaseID, unclaimedOnly, friedelOnly, iSpotLo, returnBothCoordTypes)
        return iterator

    def nSpotBins(self):
        n = 0
        for phaseID, planeDataThisPhase in self.planeDataDict.iteritems():
            n += len(planeDataThisPhase['hklInv'])
        return n
    
    def checkClaims(self, indices=None):
        """
        careful: unlike claimSpots, does not check grain association;
        careful: indices should not be a boolean array
        """
        # 'check __invalidated just in case'
        if indices is None:
            claimed = num.array(self.__claimedBy, dtype=bool)
            # claimed[self.__invalidated] = True
        else:
            claimed = num.array(self.__claimedBy[indices], dtype=bool)
            # claimed[self.__invalidated[indices]] = True
        if self.__multiprocMode:
            """
            claimed = num.logical_or(
                num.array(self.__claimedBy[indices], dtype=bool),
                num.array(self.__marks[indices], dtype=bool)
                )
            """
            whereUnclaimed = num.where(-claimed)[0]
            if indices is None:
                for iSpot in whereUnclaimed:
                    claimed[iSpot] = self.__marks[iSpot] > 0
            else:
                for iInd in whereUnclaimed:
                    iSpot = indices[iInd]
                    claimed[iInd]  = self.__marks[iSpot] > 0
        return claimed
    def checkClaim(self, index):
        'careful: unlike claimSpots, does not check grain association'
        claimed = False
        if self.__multiprocMode:
            claimed = self.__marks[index] > 0
        claimed = claimed or (self.__claimedBy[index] is not None) # or self.__invalidated[index]
        return claimed
    def claimSpots(self, indices, newClaimedBy, checkOnly=False, asMaster=False):
        """
        careful: indices should not be a boolean array;
        
        returns bool of same length as indices, with True entries for indices
        that are in conflict;
        if checkOnly, then do not actaully claim the spots
        """
        
        indices = num.atleast_1d(indices)
        
        claimed = self.checkClaims(indices=indices)
        
        if self.__multiprocMode: #  and not asMaster
            """
            do marks regardless of asMaster value, in case marks are not
            otherwise going to be set; marks are important in multiprocessing as
            they get shared while claimedBy does not
            """
            if not checkOnly:
                for iSpot in indices:
                    self.__marks[iSpot] += 1
        
        if asMaster or not self.__multiprocMode:
            # 'do not check self.__marks here on purpose'
            # claimed = num.array(self.__claimedBy[indices], dtype=bool)
            for iInd, iSpot in enumerate(indices):
                addedClaim = False
                if not claimed[iInd] : 
                    if not checkOnly:
                        self.__claimedBy[iSpot] = newClaimedBy
                        addedClaim = True
                else:
                    
                    # if self.__invalidated[iSpot]:
                    #     'claimed should already be set True for invalidated spots, but reinforce that here'
                    #     claimed[iInd] = True
                    #     # if self.__claimedBy[iSpot] != self.__claimInvalid:
                    #     #     raise RuntimeError, 'this may not really be a problem, but say far it is unexpected so should figure out why it is happening: __claimedBy is not None for an invalidated spot'
                    # else:
                    claimedBy = self.__claimedBy[iSpot]
                    if isinstance(claimedBy,list): # if hasattr(claimedBy,'__len__'):
                        if claimedBy.count(newClaimedBy):
                            'no need to add to list'
                            claimed[iInd] = False 
                        else:
                            if not checkOnly:
                                self.__claimedBy[iSpot] += [newClaimedBy]
                                addedClaim = True
                    else:
                        if claimedBy is newClaimedBy:
                            'claimed only by newClaimedBy, so no conflict for newClaimedBy'
                            claimed[iInd] = False
                        else:
                            if not checkOnly:
                                self.__claimedBy[iSpot] = [claimedBy, newClaimedBy]
                                addedClaim = True
                if self.__spots is not None and not checkOnly:
                    if addedClaim:
                        spot = self.__spots[iSpot]
                        if isinstance(spot, tuple):
                            'now deal with master spot'
                            iSpotMaster, iSubSpot = spot
                            claimedBy = self.__claimedBy[iSpotMaster]
                            if claimedBy is None:
                                self.__claimedBy[iSpotMaster] = newClaimedBy
                            elif isinstance(claimedBy,list):
                                if claimedBy.count(newClaimedBy):
                                    # raise RuntimeError, 'master spot already claimed by this grain'
                                    print 'WARNING : master spot %d already claimed by this grain %s' % (iSpot, str(newClaimedBy.lineageList))
                                self.__claimedBy[iSpotMaster] += [newClaimedBy]
                            else:
                                if claimedBy is newClaimedBy:
                                    # raise RuntimeError, 'master spot already claimed by this grain'
                                    print 'WARNING : master spot %d already claimed by this grain %s' % (iSpot, str(newClaimedBy.lineageList))
                                self.__claimedBy[iSpotMaster] = [claimedBy, newClaimedBy]

        
        return claimed
    def cleanMulti(self):
        """
        clean up spots for which split fitting was done, as they might
        need to change to single spots
        """
        assert self.__spots is not None,\
            'for cleanMulti, need spot instances'
        for spot in self.__spots:
            if isinstance(spot, tuple):
                'nothing to be done'
                continue
            if spot.getNSplit() > 0:
                spot.unsetMulti() # does spot.cleanFit()
        
        return
    def resetClaims(self, inPlace=False):
        nSpots = len(self)
        
        if inPlace:
            oldLen = num.size( self.__claimedBy )
            self.__claimedBy.resize( nSpots )
            self.__claimedBy[oldLen:] = None
            if self.__keepMarks:
                raise RuntimeError, 'do not resetClaims with keepMarks set'
            if self.__multiprocMode:
                raise RuntimeError, 'do not resetClaims inPlace when in multiprocMode'
            else:
                self.__marks = None
        else:
            self.__claimedBy = num.empty( nSpots, dtype=object )
            self.__claimedBy[:] = None # may not be necessary
            if not self.__keepMarks:
                if self.__multiprocMode:
                    self.__marks[:] = num.zeros(len(self.__marks), dtype=int)
                else:
                    self.__marks = None
        
        # self.__claimedBy[self.__invalidated] = self.__claimInvalid
            
        return
    
    def reportFitFailures(self):
        assert self.__spots is not None, \
            'need list of Spot instances'
        nFailed  = 0
        nFailedM = 0
        for spot in self.__spots[0:self.__nRawSpots]:
            failed = spot.fitHasFailed()
            if failed:
                nFailed = nFailed + 1
                if spot.getNSplit() > 1:
                    nFailedM = nFailedM + 1
        return nFailed, nFailedM
    def report(self, pw=None):
        assert not self.__multiprocMode, 'do not report in multiprocmode'
        nBins = self.nSpotBins()
        nSpotsTThAll     = num.empty(nBins) # dtype=int
        nSpotsTThClaimed = num.empty(nBins) # dtype=int
        iBin = 0
        #claimed = num.array(self.__claimedBy, dtype=bool)
        claimed = self.checkClaims()
        for phaseID, planeDataThisPhase in self.planeDataDict.iteritems():
            planeDataEntry = planeDataThisPhase['planeData']
            hklList        = planeDataThisPhase['hklList']
            tThAssoc       = planeDataThisPhase['tThAssoc']
            multip         = planeDataThisPhase['multip']
            for iHKL, hkl in enumerate(hklList):
                # hklID = planeDataThisPhase['hklIDList'][iHKL]
                multipInv = 1.0 / float(multip[iHKL])
                nSpotsTThAll[iBin]     = num.sum(tThAssoc[iHKL, :]) * multipInv
                nSpotsTThClaimed[iBin] = num.sum(num.logical_and(tThAssoc[iHKL, :], claimed)) * multipInv
                iBin += 1
        if pw is not None:
            pw.clear()
            
            width = 0.8
            x = num.arange(nBins)
            bA = pw.a.bar(x, nSpotsTThAll,     width, color='0.5') # bottom=nSpotsTThClaimed
            bC = pw.a.bar(x, nSpotsTThClaimed, width, color='r')
            
            yLabel      = pw.a.set_ylabel('# spots / hkl multiplicity')
            xTicks      = pw.a.set_xticks(x+width/2.)
            xTickLables = pw.a.set_xticklabels( map(str, self.hklIDListAll), rotation='vertical' )
            title       = pw.a.set_title('spot report')
            legend      = pw.a.legend( (bC[0], bA[0]), ('claimed','total') ) # loc=(0.05,0.85)
            
            pw.show()
        
        numClaimed = num.sum(num.array(self.__claimedBy, dtype=bool))
        numMultiClaimed = num.sum(num.array(map(lambda x: isinstance(x,list), self.__claimedBy)))
        return nSpotsTThAll, nSpotsTThClaimed, numClaimed, numMultiClaimed # self.hklIDListAll

    def __markFriedel(self, iSpot, jSpot):
        if self.friedelPair[iSpot] == -1 and self.friedelPair[jSpot] == -1:
            self.friedelPair[iSpot] = jSpot
            self.friedelPair[jSpot] = iSpot
        else:
            print 'WARNING: Attempted to multiply assign friedel pair'
        return
    def findFriedelPairsHKL(self, hkl,
                            phaseID=None,
                            tthTol=None,
                            etaTol=valUnits.valWUnit('etaTol', 'angle', 0.25, 'degrees'),
                            omeTol=valUnits.valWUnit('etaTol', 'angle', 1.00, 'degrees')):

        getFriedelPair = crystallography.getFriedelPair
        angularDifference = Rotations.angularDifference

        if hasattr(etaTol, 'getVal'):
            etaTol = etaTol.getVal('radians')
        
        if hasattr(omeTol, 'getVal'):
            omeTol = omeTol.getVal('radians')
        
        # loop time
        theseSpotsB = self.getHKLSpots(hkl, phaseID) # boolean over all spots
        theseSpotsI = num.where(theseSpotsB)[0]
        for iThese, iSpot in enumerate(theseSpotsI):
            if self.friedelPair[iSpot] >= 0: continue
            angCOM_I = self.__spotAngCoords[iSpot]
            angCOM_restOfThem = self.__spotAngCoords[theseSpotsI[iThese+1:], :]
            
            # get theoretical Friedel pair of I
            omeFP_I, etaFP_I = getFriedelPair(angCOM_I[0], angCOM_I[1], angCOM_I[2],
                                              display=False,
                                              units='radians',
                                              convention='aps')
                        
            # take difference and find hits for tolerances
            fpDist = angularDifference(
                num.tile(num.r_[etaFP_I, omeFP_I], (angCOM_restOfThem.shape[0], 1)),
                angCOM_restOfThem[:, 1:] )
            
            # scale by tolerances so that rowNorm makes sense
            fpDist[:,0] = fpDist[:,0] * (1.0/etaTol)
            fpDist[:,1] = fpDist[:,1] * (1.0/omeTol)
            
            hits = num.where( num.logical_and(
                    (fpDist[:, 0] <= 1.0),
                    (fpDist[:, 1] <= 1.0)
                    ) )[0]
            
            if len(hits) > 0:           # found some
                if len(hits) > 1:       # found more than one; take closest
                    minInHits = num.argmin(mUtil.rowNorm(fpDist[hits, :]).squeeze())
                    jSpot     = theseSpotsI[iThese+1+minInHits]
                else:                   # found exactly one
                    jSpot     = theseSpotsI[iThese+1+hits[0]]
                self.__markFriedel(iSpot, jSpot)
        return
    @staticmethod
    def newActiveSpotData(spot, indexs):
        isNew = True
        assert isinstance(indexs, list), 'indexsPrev should be a list'
        d = {
            'spot'       : spot,
            'indexsPrev' : [],
            'slaved'     : None,
            'indexsNew'  : []
            }
        if isNew:
            d['indexsNew'] = indexs
        else:
            d['indexsPrev'] = indexs
        return d
    @staticmethod
    def enslaveSpotData(d, dMaster):
        """ends up, after potentially doing recursive calls, 
        slaving master of d to master of dMaster"""
        realMaster = Spots.getMaster(dMaster)
        if d is realMaster:
            'nothing to do'
            pass
        else:
            assert realMaster['slaved'] is None, 'slaving to a slave, should not happen'
            if d['slaved'] is None:
                d['slaved'] = realMaster
            else:
                Spots.enslaveSpotData(d['slaved'], realMaster)
        return
    @staticmethod
    def mergeSpotData(d):
        merged = False
        if d['slaved'] is not None:
            merged = True
            dMaster = Spots.getMaster(d)
            dMaster['spot'].merge(d['spot']) # grab spot data
            dMaster['indexsPrev'] += d['indexsPrev'] # probably not needed, but do anyway
            dMaster['indexsNew'] += d['indexsNew']
        return merged
    @staticmethod
    def ageSpotData(d):
        # d['spot'] stays unchanged
        d['indexsPrev'] = num.unique(d['indexsNew']).tolist() # get rid of redundant entries from merging
        assert d['slaved'] is None, 'should not be aging slaved spots'
        d['indexsNew'] = []
        assert len(d['indexsPrev']) > 0, 'empty indexsPrev'
        return
    @staticmethod
    def getMaster(d):
        if not d['slaved']:
            retval = d
        else:
            retval = Spots.getMaster(d['slaved'])
        return retval
    @staticmethod
    def findSpotsOmegaStack(reader, 
                            nFrames,
                            threshold, minPx, 
                            discardAtBounds=True,
                            keepWithinBBox=True,
                            overlapPixelDistance=None,
                            nframesLump=1, # probably get rid of this eventually
                            padOmega=True,
                            padSpot=True,
                            debug=False, pw=None,
                            fout=None,
                            sumImg=None):
        """
        This method does not necessarily need to hang off of the Spots class, 
        but Spots is a convenient place to put it.
        
        If nFrames == 0, read all frames.
        
        if pass overlapPixelDistance, then overlap is determined based
        on centroid distances in pixel units; with
        overlapPixelDistance being set to the tolerance for overlap;
        probably most useful when omega steps are large
        
        reader has been created by doing something like:
        	fileInfo = [('RUBY_4537.raw', 2), ('RUBY_4538.raw', 2)]
                reader = detector.ReadGE(fileInfo, subtractDark=True)
        
        if go to parallel processing, perhaps return first and last lables for doing merges
        """
        location = '  findSpotsOmegaStack'
        fout = fout or sys.stdout
        def log_l(message):
            print >> fout, location+' : '+message
        
        if overlapPixelDistance:
            assert isinstance(overlapPixelDistance,float), \
                'if specified, overlapPixelDistance should be a float'
            if debug:
                log_l(' will use pixel distance of %g for overlap determination' % (overlapPixelDistance))
        if padOmega:
            '''
            The idea of padding in omega is to include background pixels in adjacent omega frames,
            but if overlapPixelDistance is being used then these pixels may not in fact be below
            the background. Just in case other combinations of settings make this a concern, just
            set cullPadOmega to True for all cases for now
            '''
            cullPadOmega = True
        
        labelsPrev       = None
        objsPrev         = None
        comsPrev         = None
        #
        'storage for finalized spots'
        spotListF = []
        #
        """keep active spot data in a dictionary, indexed by (original iFrame, original index in iFrame)"""
        spotDictA = {} # active spots
        #
        delta_omega = reader.getDeltaOmega()
        delta_omega = delta_omega * nframesLump
        if sumImg is None:
            sumImg = nframesLump > 1

        if nFrames == 0:
            nFrames = reader.getNFrames()
        if debug:
            log_l('using reader with %d frames and delta-omega of %g' % (nFrames, delta_omega))
        #
        prevframe = None
        prevbin   = None
        prevomega = None
        darkframe = reader.getDark()
        for iFrame in range(nFrames / nframesLump): # for iFrame, omega in enumerate(omegas):
            if debug > 1:
                log_l('working on frame %d' % (iFrame))
                ticFrame = time.time()
            thisframe = reader.read(nframes=nframesLump, sumImg=sumImg) 
            omega = reader.getFrameOmega()
            
            if debug > 1:
                tic = time.time()
            'call spotFinderSingle with minPx == 1 so it does not discard any spots'
            if pw is not None:
                reader.display(thisframe, pw=pw)
            if debugFrameRefs:
                log_l('references to thisframe : %d' % (sys.getrefcount(thisframe)))
            labels, objs, spotDataList, bin = spotFinderSingle(thisframe, threshold, 1, 
                                                               keepWithinBBox, padSpot,
                                                               debug=debug>3, pw=pw, darkframe=darkframe)
            if debug > 1:
                toc = time.time()
                log_l('    spotFinderSingle took %g seconds' % ((toc - tic)))
            if debug > 1:
                log_l('len(spotDataList) : %d' % (len(spotDataList)))
                
            if debugFrameRefs:
                log_l('references to thisframe : %d (after spotFinderSingle)' % (sys.getrefcount(thisframe)))
            if overlapPixelDistance:
                'will need centers of mass'
                coms = num.empty((len(objs),2))
                for iObj in range(len(objs)):
                    index = iObj+1
                    # coms[iObj] = ndimage.center_of_mass(thisframe, labels=labels, index=index) # slow
                    coms[iObj] = getImCOM(thisframe, labels, objs, index)
            
            if labelsPrev is None:
                """starting with all new spots
                mark for deletion as touch the lower omega limit
                """ 
                if debug > 1:
                    log_l('starting all new spots')
                for spotData in spotDataList:
                    key = (iFrame, spotData['index'])
                    """
                    even if padOmega is True, do not have a previous frame with which to pad;
                    careful of this if doing multiprocessing
                    """
                    spot = Spot(key, delta_omega, spotData, omega, iFrame)
                    spot.mark(Spot.markAtOmegaLo)
                    indexPrev = spotData['index']
                    spotDictA[key] = Spots.newActiveSpotData(spot, [indexPrev])
            else:
                
                if debug > 2:
                    ticClaim = time.time()
                keysToPop = [] # spotsToPop = []
                spotClaimedBy = [None for iSpot in range(len(spotDataList))]
                #
                if debug > 2:
                    ttIndU = 0.
                    ttOver = 0.
                    ttNotO = 0.
                if debug > 1:
                    log_l('working on %d active spots' % (len(spotDictA)))
                for key, activeSpotData in spotDictA.iteritems(): 
                    activeSpot = activeSpotData['spot']
                    
                    foundOverlap = False
                    if debug > 2:
                        log_l('    indexsPrev : '+str(activeSpotData['indexsPrev']))
                    for indexPrev in activeSpotData['indexsPrev']:
                        if overlapPixelDistance:
                            """
                            this is a funny way of doing centroidal
                            overlap, but is done this way to be
                            consistent with the other method for
                            determining overlap
                            """
                            comPrev = comsPrev[indexPrev-1]
                            spotIndexsUnique = [0]
                            for iObj in range(len(objs)):
                                index = iObj+1
                                distp   = coms[iObj] - comPrev
                                dist    = num.linalg.norm(distp)
                                if dist <= overlapPixelDistance:
                                    spotIndexsUnique.append(index)
                        else:
                            if debug > 2:
                                tic = time.time()
                            spotIndexs = getValuesOnly(labels, labelsPrev, objsPrev, indexPrev)
                            spotIndexsUnique = num.unique(spotIndexs)
                            if debug > 2:
                                ttIndU += time.time() - tic
                        # overlap = ndimage.sum(labels, labels=labelsPrev, index=indexPrev) > 0
                        if debug > 2:
                            tic = time.time()
                        if debug > 2:
                            log_l('    spotIndexsUnique : '+str(spotIndexsUnique))
                        overlap = spotIndexsUnique[-1] > 0
                        if overlap:
                            foundOverlap = True
                            'get rid of 0 as necessary'
                            if spotIndexsUnique[0] == 0:
                                spotIndexsUnique = spotIndexsUnique[1:]
                            for index in spotIndexsUnique:
                                
                                iSpot = index-1
                                spotData = spotDataList[iSpot]
                                assert spotData['index'] == index, 'index mismatch'
                                activeSpot.append(spotData, omega, iFrame)
                                activeSpotData['indexsNew'] += [index]
                                
                                if spotClaimedBy[iSpot] is not None: # not unclaimed[iSpot]:
                                    Spots.enslaveSpotData(activeSpotData, spotClaimedBy[iSpot])
                                else:
                                    spotClaimedBy[iSpot] = activeSpotData
                        if debug > 2:
                            ttOver += time.time() - tic
                    'done looping over indexPrev for this active spot'
                    if debug > 2:
                        tic = time.time()
                    
                    if foundOverlap:
                        if padOmega:
                            for indexPrev in activeSpotData['indexsPrev']:
                                'previous frame onto this one'
                                spotDataPad = getSpot(thisframe, labelsPrev, objsPrev, indexPrev, 
                                                      keepWithinBBox, False, darkframe) # padSpot
                                cullSpotUsingBin(spotDataPad, bin)
                                activeSpot.append(spotDataPad, omega, iFrame)
                            for index in activeSpotData['indexsNew']:
                                'this frame onto previous one'
                                spotDataPad = getSpot(prevframe, labels, objs, index, 
                                                      keepWithinBBox, False, darkframe) # padSpot
                                cullSpotUsingBin(spotDataPad, prevbin)
                                activeSpot.append(spotDataPad, prevomega, iFrame-1)
                    else:
                        'finalize this spot -- no overlap with newly found spots'
                        if padOmega:
                            for indexPrev in activeSpotData['indexsPrev']:
                                spotDataPad = getSpot(thisframe, labelsPrev, objsPrev, indexPrev, 
                                                      keepWithinBBox, False, darkframe) # padSpot
                                if cullPadOmega:
                                    cullSpotUsingBin(spotDataPad, bin)
                                activeSpot.append(spotDataPad, omega, iFrame)
                        activeSpot.finalize()
                        spotListF.append(activeSpot)
                        keysToPop.append(key) # spotsToPop.append(spot) # spot
                    if debug > 2:
                        ttNotO += time.time() - tic
                'end for key, activeSpotData in spotDictA.iteritems()'
                if debug > 2:
                    tocClaim = time.time()
                    log_l('    claiming spots took %g seconds' % ((tocClaim - ticClaim)))
                    log_l('        for IndU : %g seconds' % (ttIndU))
                    log_l('        for Over : %g seconds' % (ttOver))
                    log_l('        for NotO : %g seconds' % (ttNotO))
                    ticMunge = time.time()
                
                nFinalize = len(keysToPop)
                'before popping, add keys for merged spots to keysToPop'
                for key, activeSpotData in spotDictA.iteritems(): 
                    merged = Spots.mergeSpotData(activeSpotData)
                    if merged:
                        keysToPop.append(key)
                nMerge = len(keysToPop) - nFinalize
                if debug:
                    log_l('finalizing %d and merging %d spots' % (nFinalize, nMerge))
                
                'now pop'
                for key in keysToPop: # for spot in spotsToPop:
                    'already put in spotListF as appropriate, just need to pop it out of spotDictA'
                    junk = spotDictA.pop(key) # spot.key
                
                'make new active spots from unclaimed spots in current frame'
                unclaimed = num.array([spotClaimedBy[iSpot] is None  
                                       for iSpot in range(len(spotDataList))],
                                      dtype=bool)
                if debug > 1:
                    log_l('creating %d new spots' % (num.sum(unclaimed)))
                for iSpotData, spotData in enumerate(spotDataList):
                    if not unclaimed[iSpotData] : continue
                    'no overlap with existing spots, make a new spot'
                    key = (iFrame, spotData['index'])
                    if padOmega:
                        'for now, keep key as it is, without changing iFrame'
                        #prevSpotData = copy.deepcopy(spotData)
                        prevSpotData = {}
                        prevSpotData.update(spotData)
                        prevSpotData['vbox']    = copyBox(prevframe, prevSpotData['obj'])
                        prevSpotData['darkbox'] = copyBox(darkframe, prevSpotData['obj'])
                        if cullPadOmega:
                            cullSpotUsingBin(prevSpotData, prevbin)
                        spot = Spot(key, delta_omega, prevSpotData, prevomega, iFrame-1)
                        spot.append(spotData, omega, iFrame)
                        pass
                    else:
                        spot = Spot(key, delta_omega, spotData, omega, iFrame)
                    indexPrev = spotData['index']
                    spotDictA[key] = Spots.newActiveSpotData(spot, [indexPrev])
                
                if debug > 2:
                    tocMunge = time.time()
                    log_l('    munging spots took %g seconds' % ((tocMunge - ticMunge)))
            'end of if for labelsPrev is None'

            if debug > 1:
                log_l('now have %d active spots' % (len(spotDictA)))
            for key, activeSpotData in spotDictA.iteritems(): 
                Spots.ageSpotData(activeSpotData)
            
            if padOmega:
                prevframe = thisframe
                prevbin   = bin
                prevomega = omega
            
            labelsPrev = labels
            objsPrev   = objs
            if overlapPixelDistance:
                comsPrev   = coms
            
            if debug > 1:
                tocFrame = time.time()
                log_l('in total, frame %d took %g seconds' % (iFrame, (tocFrame - ticFrame)))
            
            if debugFrameRefs:
                log_l('references to thisframe : %d (end of loop)' % (sys.getrefcount(thisframe)))
            del thisframe
        'done with iFrame loop'
        
        'take care of remaining active spots'
        for key, activeSpotData in spotDictA.iteritems(): 
            activeSpot = activeSpotData['spot']
            #
            'mark active spots as on an omega bound' 
            activeSpot.mark(Spot.markAtOmegaHi)
            #
            """
            even if padOmega is True, do not have a next frame with which to pad;
            careful of this if doing multiprocessing
            """
            activeSpot.finalize()
            spotListF.append(activeSpot)
        
        keepers = num.array(map(Spot.nPx, spotListF)) >= minPx # dtype=bool
        if debug > 1:
            log_l('%d of %d spots have more than %d pixels' % (num.sum(keepers), len(spotListF), minPx))
        if discardAtBounds:
            if hasattr(discardAtBounds, '__len__'):
                toDiscard = discardAtBounds
            else:
                toDiscard = Spot.atBoundMarks
            nDiscard = 0
            for iSpot in range(len(spotListF)): # spot in spotListF:
                if not keepers[iSpot] : continue
                spot = spotListF[iSpot]
                if spot.isMarked(toDiscard):
                    keepers[iSpot] = False
                    nDiscard += 1
            if debug > 1:
                log_l('discarded %d spots at omega bounds' % (nDiscard))
        
        keptSpots = []
        for iSpot in range(len(spotListF)): # spot in spotListF:
            if not keepers[iSpot] : continue
            spot = spotListF[iSpot]
            keptSpots.append(spot)
        
        return keptSpots
    def getOmegaMins(self):
        mins = zip(*self.omegaMM)[0]
        'toFloat already done'
        #retval = num.array(toFloat(mins, 'radians'))
        retval = num.array(mins)
        return retval
    def getOmegaMaxs(self):
        maxs = zip(*self.omegaMM)[1]
        'toFloat already done'
        #retval = num.array(toFloat(maxs,'radians'))
        retval = num.array(maxs)
        return retval

class SpotsIterator:
    'iterator over a given set of spots in a Spots instance, note that the iterator ignores relations among split spots'
    def __init__(self, spots, hkl, phaseID, unclaimedOnly, friedelOnly, iSpotLo, returnBothCoordTypes):
        self.spots = spots
        if phaseID is None:
            assert len(self.spots.planeDataDict) == 1, 'need phaseID if multiphase'
            phaseID = self.spots.planeDataDict.keys()[0]
        planeDataDict = self.spots.planeDataDict[phaseID]
        if hkl is None:
            iHKL = None
        else:
            iHKL = planeDataDict['hklInv'][tuple(hkl)]
        self.__iter_pnt             = iSpotLo # 0
        self.__iter_pnt_0           = iSpotLo
        self.__iter_phaseID         = phaseID
        self.__iter_iHKL            = iHKL
        self.__iter_friedelOnly     = friedelOnly
        self.__iter_unclaimedOnly   = unclaimedOnly
        self.__returnBothCoordTypes = returnBothCoordTypes
        return
    def __iter__(self):
        self.__iter_pnt = self.__iter_pnt_0
        return self
    def next(self):
        nSpots    = len(self.spots._Spots__spotAngCoords)
        planeDataDict = self.spots.planeDataDict[self.__iter_phaseID]
        if self.__iter_iHKL is None:
            tThAssoc  = planeDataDict['phaseAssoc']
        else:
            tThAssoc  = planeDataDict['tThAssoc'][self.__iter_iHKL, :]
        while True:
            if self.__iter_pnt >= nSpots:
                raise StopIteration
            
            if not tThAssoc[self.__iter_pnt]: 
                self.__iter_pnt += 1
                continue
            if self.__iter_unclaimedOnly:
                if self.spots.checkClaim(self.__iter_pnt):
                    self.__iter_pnt += 1
                    continue
            
            if self.__iter_friedelOnly:
                if self.spots.friedelPair[self.__iter_pnt] > self.__iter_pnt: 
                    """ checking > self.__iter_pnt will only be True for the
                    master spot if have done:
                    	friedelPair[iSpot] = jSpot
                        friedelPair[jSpot] = iSpot
                    to set up the data
                    """
                    iSpot = self.__iter_pnt
                    jSpot = self.spots.friedelPair[self.__iter_pnt]
                    if self.__returnBothCoordTypes:
                        raise NotImplementedError, 'both coord types and Friedel pairs'
                    retval = iSpot, jSpot, \
                             self.spots._Spots__spotAngCoords[iSpot], \
                             self.spots._Spots__spotAngCoords[jSpot]
                else:
                    self.__iter_pnt += 1
                    continue
            else:
                iSpot = self.__iter_pnt
                if self.__returnBothCoordTypes:
                    retval = iSpot, \
                        self.spots._Spots__spotAngCoords[iSpot], \
                        self.spots._Spots__spotXYOCoords[iSpot]
                else:
                    retval = iSpot, \
                        self.spots._Spots__spotAngCoords[iSpot]
            'found one if get to here'
            break
        'increment for next call'
        self.__iter_pnt += 1
        return retval

######################################################################

def testSingle(fileInfo):
    pw = plotWrap.PlotWrap()
    
    reader = detector.ReadGE(fileInfo, subtractDark=True)
    
    thisframe = reader.read()
    # thisframe = reader.read(nskip=369)
    # reader.display(thisframe, range=600, pw=pw)
    reader.display(thisframe, pw=pw)
    # detector.ReadGE.display(thisframe, range=200, pw=pw)
    labels, objs, spotDataList, bin = spotFinder.spotFinderSingle(thisframe, 20, 2, False, pw=pw, debug=True)
    
    return reader, pw, labels, objs, spotDataList, thisframe
    
def testSpotFinder(fileInfo, delta_omega, omega_low, howMuch=0.1):
    
    tic = time.time()
    
    mask = None
    nFrames = detector.ReadGE.getNFramesFromFileInfo(fileInfo)
    nFrames = int(nFrames*howMuch)
    if howMuch == 1:
        maskData = num.loadtxt('RUBY_MASK_CSC.dat', dtype=int, comments='#')
        mask = reader.getEmptyMask()
        mask[maskData[:,0],maskData[:,1]] = maskData[:,2]
        pMask = detector.ReadGE.display(mask)
        del maskData

    print 'will test with %d frames' % (nFrames)
    #omegas = num.arange(omega_low, omega_low+delta_omega*(nFrames-0.5), delta_omega)
    
    reader = detector.ReadGE(fileInfo, subtractDark=True, mask=mask)
    
    threshold = 20
    minPx = 4
    spots = spotFinder.Spots.findSpotsOmegaStack(
        reader, omega_low, delta_omega, nFrames, 
        threshold, minPx, debug=1) # pw=pw
    spotSizes = map(spotFinder.Spot.nPx, spots)
    print 'spot sizes : '+str(spotSizes)
    # locmax = num.argmax(spotSizes)
    # spot = spots[locmax]
    # win = spot.display(cmap=None)
    
    toc = time.time()
    dt = toc - tic
    print 'finding spots took %g seconds to run' % (dt)
    
    return spots

def main(argv=[]):
    
    #fileInfo = 'RUBY_4537.raw' 
    # ... need to check nempty values!
    fileInfo = [('RUBY_4537.raw', 2), ('RUBY_4538.raw', 2)]
    delta_omega = 0.25 * detector.d2r # 0.25 degrees, converted to radians
    # ... need to get correct omega_low 
    omega_low = -60.0 * detector.d2r # converted to radians
    
    howMuch = 0.1
    if len(argv) > 0:
        howMuch = argv[0]
    spots = spotFinder.testSpotFinder(fileInfo, delta_omega, omega_low, howMuch=howMuch)
    
    'cull spots and the sort by size'
    tests = [
        spotFinder.Spot.testSpotIsNotStreak,
        (spotFinder.Spot.testSpotLenAtLeast, 5)
        ]
    culledSpots = spotFinder.Spot.cullSpots(spots, tests)
    #
    sortedSpots = arrayUtil.structuredSort( map(len, culledSpots), culledSpots )

    plotSpots = False
    if plotSpots: 
        'plot each spot in turn, but first cull and then sort by size'
        
        'sortedSpots is smallest to biggest, so reverse for plotting biggest to smallest'
        for spot in sortedSpots[-1:0:-1]: # for spot in spots:
            win = spot.display(cmap=None)
            raw_input("any key to continue")
            win.destroy()

    spot = sortedSpots[-1] # biggest spot
    # win = spot.display(cmap=None)
    workDist = valUnits.valWUnit('workDist', 'length', 1.9365, 'meter')
    detectorGeom = detector.DetectorGeomGE(1024, 1024, workDist)
    fitX, vCalc = spot.fit(detectorGeom=detectorGeom)
    winData = spot.display(title='data')
    winFit  = spot.display(vAll=vCalc, title='fit')
    winDiff = spot.display(vAll=(vCalc-spot.vAll), title='difference')
    print 'angCOM estimate : '+str(spot.angCOM(useFit=False))
    print 'angCOM fit      : '+str(spot.angCOM(useFit=True))
    
    return spots

if __name__ == '__main__':
    import cProfile as profile # profile
    pr = profile.Profile()
    pr.run('spots = main(sys.argv[1:])') #  filename='main.cprof'
    pr.dump_stats('spotFinder_main.cprof')
    #import pstats
    #stats = pstats.Stats('spotFinder_main.cprof')
    print 'Profiler stats : '
    pr.print_stats(sort=1) # stats.print_stats(sort=1)
    shelf = shelve.open('spots.py_shelf')
    shelf['spots'] = spots
    shelf.close()
    
