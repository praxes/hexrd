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

import math
import copy
import os
import time
import sys

try:
    import Image
    haveImageModule = True
except:
    haveImageModule = False

haveThreading = True
try:
    import threading
except:
    haveThreading = False

import numpy as num
from scipy import sparse
from scipy.optimize import fsolve
from scipy.optimize import leastsq
import matplotlib
from matplotlib import cm, colors
from matplotlib import mlab
from matplotlib.widgets import Slider, Button, RadioButtons

from hexrd import xrd
from hexrd.xrd.xrdbase import getGaussNDParams, dataToFrame
from hexrd.xrd.rotations import mapAngle
from hexrd.xrd.rotations import rotMatOfExpMap
from hexrd.xrd.rotations import rotMatOfExpMap, arccosSafe
from hexrd.quadrature import q1db
from hexrd.quadrature import q2db
from hexrd.matrixutil import unitVector
from hexrd import valunits

havePlotWrap = True
try:
    from hexrd import plotwrap
except:
    havePlotWrap = False

d2r = piby180 = num.pi/180.
r2d = 1.0/d2r

bsize = 25000
ztol  = 1e-10

DFLT_XTOL = 1e-6
        
def angToXYIdeal(tTh, eta, workDist):
    rho = num.tan(tTh) * workDist
    x   = rho * num.cos(eta)
    y   = rho * num.sin(eta)
    return x, y

def mapAngs(eta, doMap=None):
    if doMap is None:
        doMap = num.all(num.cos(eta) < 0)
    if doMap:
        eta = mapAngle(eta, (0, 2*num.pi), units='radians')
    return eta

class FmtCoordIdeal:
    def __init__(self, planeData, workDist):
        self.planeData = planeData
        self.workDist  = workDist
        self.dataList  = []
        return
    def addDetectorData(self, iDetector, *data):
        self.dataList.append((iDetector,copy.deepcopy(data)))
    def __call__(self, x, y):
        rho = num.sqrt(x*x + y*y)
        eta = num.arctan2(y, x)
        tTh = num.arctan2(rho, self.workDist)
        for iDetector, data in self.dataList:
            xedges, yedges, h, mask = data
            ix = xedges.searchsorted(x)
            iy = yedges.searchsorted(y)
            if ix == 0 or ix == len(xedges) or iy == 0 or iy == len(yedges):
                intens = None
            elif mask[ix-1,iy-1]:
                intens = None
            else:
                intens = h[ix-1,iy-1]
            if intens is not None:
                if self.planeData is None:
                    retval = "detector %d : tth=%g eta=%g int=%g" \
                        % (iDetector, r2d*tTh, r2d*eta, intens)
                else:
                    dsp    = 0.5 * self.planeData.wavelength / num.sin(0.5*tTh)
                    HKLs   = str(self.planeData.getHKLs(asStr=True, allHKLs=True, thisTTh=tTh))
                    retval = "detector %d : d=%g tth=%g eta=%g int=%g \n HKLs=%s" \
                        % (iDetector, dsp, r2d*tTh, r2d*eta, intens, HKLs)
                return retval
        retval = 'off-detector'
        return retval

class ThreadReadFrame(threading.Thread):
    def __init__(self, img, readArgs, castArgs):
        threading.Thread.__init__(self)

        self.img   = img
        # self.dtype = dtype
        # self.count = count
        self.readArgs = readArgs
        self.castArgs = castArgs

        self.data = None
        # self.success  = None

        return

    def run(self):
        try:
            readData = num.fromfile(self.img, **self.readArgs)
            self.data = num.array(readData, **self.castArgs)
            self.success = True
        except:
            self.success = False
        return

# if readGE gets a base class, perhaps hang this there
def getCentered(vmin, vmax):
    centered = bool(-vmin > 0.4*vmax) # vmin < 0
    return centered
def getCMap(spec):
    if isinstance(spec, bool):
        if spec:
            cmap = cm.RdBu
        else:
            cmap = cm.bone
    elif isinstance(spec, str):
        if spec == 'hotAndCold':
            cdict = {'red'   : ((0.0,   0., 0.),
                                (0.457, 0., 0.),
                                (0.5,   1., 1.),
                                (0.542, 1., 1.),
                                (1.0,   1., 1.),
                                ),
                     'green' : ((0.0,   0., 0.),
                                (0.457, 1., 1.),
                                (0.5,   1., 1.),
                                (0.542, 1., 1.),
                                (1.0,   0., 0.),
                                ),
                     'blue'  : ((0.0,   1., 1.),
                                (0.457, 1., 1.),
                                (0.5,   1., 1.),
                                (0.542, 0., 0.),
                                (1.0,   0., 0.),
                                ),
                     }
            cmap = colors.LinearSegmentedColormap('hotAndCold',cdict,256)
        else:
            raise RuntimeError, 'unknown: '+str(spec)
    else:
        raise RuntimeError, 'unknown: '+str(spec)
    return cmap


class Framer2DRC(object):
    """
    Base class for readers.

    You can make an instance of this class and use it for most of the 
    things a reader would do, other than actually reading frames
    """
    def __init__(self,
                 ncols, nrows,
                 dtypeDefault='int16', dtypeRead='uint16', dtypeFloat='float64'):
        self.__ncols = ncols
        self.__nrows = nrows
        self.__frame_dtype_dflt  = dtypeDefault
        self.__frame_dtype_read  = dtypeRead
        self.__frame_dtype_float = dtypeFloat
        
        self.__nbytes_frame  = num.nbytes[dtypeRead]*nrows*ncols
        
        return
    
    def get_ncols(self):
        return self.__ncols
    ncols = property(get_ncols, None, None)

    def get_nbytesFrame(self):
        return self.__nbytes_frame
    nbytesFrame = property(get_nbytesFrame, None, None)

    def get_nrows(self):
        return self.__nrows
    nrows = property(get_nrows, None, None)
    
    def get_dtypeDefault(self):
        return self.__frame_dtype_dflt
    dtypeDefault = property(get_dtypeDefault, None, None)
    def get_dtypeRead(self):
        return self.__frame_dtype_read
    dtypeRead = property(get_dtypeRead, None, None)
    def get_dtypeFloat(self):
        return self.__frame_dtype_float 
    dtypeFloat = property(get_dtypeFloat, None, None)
    
    def getOmegaMinMax(self):
        raise NotImplementedError
    def getDeltaOmega(self):
        'needed in findSpotsOmegaStack'
        raise NotImplementedError
    def getNFrames(self):
        """
        number of total frames with real data, not number remaining
        needed in findSpotsOmegaStack
        """
        raise NotImplementedError
    def read(self, nskip=0, nframes=1, sumImg=False):
        'needed in findSpotsOmegaStack'
        raise NotImplementedError
    def getDark(self):
        'needed in findSpotsOmegaStack'
        raise NotImplementedError
    def getFrameOmega(self, iFrame=None):
        'needed in findSpotsOmegaStack'
        raise NotImplementedError


    @classmethod
    def maxVal(cls, dtypeRead):
        """
        maximum value that can be stored in the image pixel data type;
        redefine as desired
        """
        maxInt = num.iinfo(dtypeRead).max
        return maxInt

    def getEmptyMask(self):
        """convenience method for getting an emtpy mask"""
        # this used to be a class method
        mask = num.zeros([self.nrows, self.ncols], dtype=bool)
        return mask

    def getSize(self):
        retval = (self.nrows, self.ncols)
        return retval

    def frame(self, nframes=None, dtype=None, buffer=None, mask=None):
        if buffer is not None and dtype is None:
            if hasattr(buffer,'dtype'):
                dtype = buffer.dtype
        if dtype is None:
            dtype = self.__frame_dtype_dflt
        if nframes is None:
            shape = (self.nrows, self.ncols)
        else:
            assert mask is None,\
                'not coded: multiframe with mask'
            shape = (nframes, self.rows, self.ncols)
        if buffer is None:
            retval = num.zeros(shape, dtype=dtype)
        else:
            retval = num.array(buffer, dtype=dtype).reshape(shape)
        if mask is not None:
            retval = num.ma.masked_array(retval, mask, hard_mask=True, copy=False)
        return retval

    @classmethod
    def display(cls,
                thisframe,
                roi = None,
                pw  = None,
                **kwargs
                ):
        # ... interpolation method that looks like max() so that do not miss peak pixels?
        
        if roi is not None:
            dROI   = thisframe[ roi[0][0]:roi[0][1], roi[1][0]:roi[1][1] ]
        else:
            dROI = thisframe
        vmin, vmax, cmap = cls.getDisplayArgs(dROI, kwargs)

        if havePlotWrap:
            if pw is None:
                p = plotwrap.PlotWrap(**kwargs)
                kwargs = {}
            else:
                p = pw
            p(dROI, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
            # 'turn off format_coord because have not made this one report correctly'
            # p.a.format_coord = lambda x,y: ''
        # elif havePylab:
        #     assert pw is None, 'do not specify pw without plotwrap'
        #     retval = pylab.imshow(dROI, vmin=vmin, vmax=vmax, cmap=cm.bone)
        else:
            raise RuntimeError, 'no plotting pacakge available'

        retval = p
        return retval
    
    @classmethod
    def getDisplayArgs(cls, dROI, kwargs):
        range     = kwargs.pop('range',None)
        cmap      = kwargs.pop('cmap',None)
        dtypeRead = kwargs.pop('dtypeRead','uint16')
            
        roiMin = dROI.min()
        roiMax = dROI.max()
        #
        centered = getCentered(roiMin, roiMax)
        if dROI.dtype == 'bool' and range is None:
            centered = False
            vmin = 0
            vmax = 1
        elif dROI.dtype == 'float64' and \
                centered and \
                range is None:
            range = 2.0*num.max(num.abs(dROI))
            thr   = 0.0
            vmin = thr-range/2
            vmax = thr+range/2
        else:
            centered = False
            vmin, vmax = cls.getVMM(dROI, range=range, dtypeRead=dtypeRead)
        #
        if cmap is None:
            cmap = getCMap(centered)
        
        return vmin, vmax, cmap

    @classmethod
    def getVMM(cls, dROI, range=None, dtypeRead='uint16'):
        if range is None:
            range = 200.
        if hasattr(range,'__len__'):
            assert len(range) == 2, 'wrong length for value range'
            vmin = range[0]
            vmax = range[1]
        else:
            thr    = dROI.mean()
            vmin = max(0,            thr-range/2) # max(dROI.min(), thr-range/2) 
            vmax = min(cls.maxVal(dtypeRead), thr+range/2)
        return vmin, vmax

def omeToFrameRange(omega, omegas, omegaDelta):
    """
    check omega range for the frames instead of omega center;
    result can be a pair of frames if the specified omega is
    exactly on the border
    """
    retval = num.where(num.abs(omegas - omega) <= omegaDelta*0.5)[0]
    return retval

def getNFramesFromBytes(fileBytes, nbytesHeader, nbytesFrame):
    assert (fileBytes - nbytesHeader) % nbytesFrame == 0,\
        'file size not correct'
    nFrames = int((fileBytes - nbytesHeader) / nbytesFrame)
    if nFrames*nbytesFrame + nbytesHeader != fileBytes:
        raise RuntimeError, 'file size not correctly calculated'
    return nFrames

class FrameWriter(Framer2DRC):
    def __init__(self, *args, **kwargs):
        self.filename        = kwargs.pop('filename')
        self.__nbytes_header = kwargs.pop('nbytesHeader', 0)
        self.__nempty        = kwargs.pop('nempty', 0)
        self.__nempty        = kwargs.pop('nempty', 0)

        Framer2DRC.__init__(self, *args, **kwargs)
        
        self.nFrame = 0
        self.img = open(self.filename, mode='wb')
        
        # skip header for now
        self.img.seek(self.__nbytes_header, 0)
        if self.__nempty > 0:
            self.img.seek(self.nbytesFrame*self.__nempty, 1)
        
        return
    def write(self, data, doAllChecks=True):
        
        # if nskip > 0:
        #     self.img.seek(self.__nbytes_frame*nskip, 1)
        
        assert len(data.shape) == 2, 'data is not 2D'
        assert data.shape[0] == self.nrows, 'number of rows is wrong'
        assert data.shape[1] == self.ncols, 'number of rows is wrong'
        
        intType = False
        
        if   num.result_type(self.dtypeRead).kind == 'u':
            intType = True
            if data.dtype.kind == 'u':
                'all set' 
            else:
                if num.any(data < 0):
                    raise RuntimeError, 'trying to write negative data to unsigned type'
                data = data.astype(self.dtypeRead)
        elif num.result_type(self.dtypeRead).kind == 'i':
            intType = True
            data = data.astype(self.dtypeRead)
        else:
            data = data.astype(self.dtypeRead)
        
        if doAllChecks and intType:
            dataMax = data.max()
            readMax = num.iinfo(self.dtypeRead).max
            if dataMax > readMax : 
                raise RuntimeError, 'max of %g greater than supported value of %g' % (dataMax, readMax)
        
        data.tofile(self.img)
        
        return
    def __call__(self, *args, **kwargs):
        return self.write(*args, **kwargs)
    def close(self):
        self.img.close()
        return

class ReadGeneric(Framer2DRC):
    '''
    may eventually want ReadGE to inherit from this, or pull common things
    off to a base class
    '''
    def __init__(self, filename, ncols, nrows, *args, **kwargs):
        self.filename        = filename
        self.__nbytes_header = kwargs.pop('nbytes_header', 0)
        self.__nempty        = kwargs.pop('nempty', 0)
        doFlip               = kwargs.pop('doFlip', False)
        self.subtractDark    = kwargs.pop('subtractDark', False)
        
        if doFlip is not False:
            raise NotImplementedError, 'doFlip not False'
        if self.subtractDark is not False:
            raise NotImplementedError, 'subtractDark not False'

        Framer2DRC.__init__(self, ncols, nrows, **kwargs)

        self.dark = None
        self.dead = None
        self.mask = None

        self.omegaStart = None
        self.omegaDelta = None
        self.omegas = None
        #
        if len(args) == 0:
            pass
        elif len(args) == 2:
            self.omegaStart = omegaStart = args[0]
            self.omegaDelta = omegaDelta = args[1]
        else:
            raise RuntimeError, 'do not know what to do with args: '+str(args)
        self.omegas = None
        if self.omegaStart is not None:
            if hasattr(omegaStart, 'getVal'):
                omegaStart = omegaStart.getVal('radians')
            if hasattr(omegaDelta, 'getVal'):
                omegaDelta = omegaDelta.getVal('radians')
            nFramesTot = self.getNFrames()
            self.omegas = \
                num.arange(omegaStart, omegaStart+omegaDelta*(nFramesTot-0.5), omegaDelta) + \
                0.5 * omegaDelta # put omegas at mid-points of omega range for frame
            omegaEnd = omegaStart+omegaDelta*(nFramesTot)
            self.omegaMin = min(omegaStart, omegaEnd)
            self.omegaMax = max(omegaStart, omegaEnd)
            self.omegaDelta = omegaDelta
            self.omegaStart = omegaStart

        if len(kwargs) > 0:
            raise RuntimeError, 'unparsed kwargs : %s' + str(kwargs.keys())
        
        self.iFrame = -1 # counter for last global frame that was read
        
        self.img = None
        if self.filename is not None:
            self.img = open(self.filename, mode='rb')
            # skip header for now
            self.img.seek(self.__nbytes_header, 0)
            if self.__nempty > 0:
                self.img.seek(self.nbytesFrame*self.__nempty, 1)
        
        return

    def getFrameUseMask(self):
        return False
    def __flip(self, thisframe):
        return thisframe

    '''
    def read(self, nskip=0, nframes=1, sumImg=False):
        
        if not nframes == 1:
            raise NotImplementedError, 'nframes != 1'
        if not sumImg == False:
            raise NotImplementedError, 'sumImg != False'
        
        data = self.__readNext(nskip=nskip)

        self.iFrame += nskip + 1
        
        return data
    '''
    def read(self, nskip=0, nframes=1, sumImg=False):
        """
        sumImg can be set to True or to something like numpy.maximum
        """
        
        if self.img is None:
            raise RuntimeError, 'no image file open'

        'get iFrame ready for how it is used here'
        self.iFrame = num.atleast_1d(self.iFrame)[-1]
        iFrameList = []
        multiframe = nframes > 1

        nFramesInv = 1.0 / nframes
        doDarkSub = self.subtractDark # and self.dark is not None

        if doDarkSub:
            assert self.dark is not None, 'self.dark is None'

        # assign storage array
        if sumImg:
            sumImgCallable = hasattr(sumImg,'__call__')
            imgOut = self.frame(dtype=self.dtypeFloat, mask=self.dead)
        elif multiframe:
            imgOut = self.frame(nframes=nframes, dtype=self.dtypeDflt, mask=self.dead)


        # now read data frames
        for i in range(nframes):

            #data = self.__readNext(nskip=nskip)
            #thisframe = data.reshape(self.__nrows, self.__ncols)
            data = self.__readNext(nskip=nskip) # .reshape(self.__nrows, self.__ncols)
            self.iFrame += nskip + 1
            nskip=0 # all done skipping once have the first frame!
            iFrameList.append(self.iFrame)
            # dark subtraction
            if doDarkSub:
                'used to have self.dtypeFloat here, but self.dtypeDflt does the trick'
                thisframe = self.frame(buffer=data,
                                       dtype=self.dtypeDflt, mask=self.dead) - self.dark
            else:
                thisframe = self.frame(buffer=data,
                                       mask=self.dead)

            # flipping
            thisframe = self.__flip(thisframe)

            # masking (True get zeroed)
            if self.mask is not None:
                if self.getFrameUseMask():
                    thisframe[self.mask] = 0

            # assign output
            if sumImg:
                if sumImgCallable:
                    imgOut = sumImg(imgOut, thisframe)
                else:
                    imgOut = imgOut + thisframe * nFramesInv
            elif multiframe:
                imgOut[i, :, :] = thisframe[:, :]
        'end of loop over nframes'

        if sumImg:
            # imgOut = imgOut / nframes # now taken care of above
            pass
        elif not multiframe:
            imgOut = thisframe

        if multiframe:
            'make iFrame a list so that omega or whatever can be averaged appropriately'
            self.iFrame = iFrameList
        return imgOut

    def getNFrames(self, lessEmpty=True):
        fileBytes = os.stat(self.filename).st_size
        nFrames = getNFramesFromBytes(fileBytes, self.__nbytes_header, self.nbytesFrame)
        if lessEmpty:
            nFrames -= self.__nempty
        return nFrames

    def getOmegaMinMax(self):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaMin, self.omegaMax
    def getDeltaOmega(self, nframes=1):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaDelta * nframes
    def getDark(self):
        'no dark yet supported'
        return 0
    def getFrameOmega(self, iFrame=None):
        """if iFrame is none, use internal counter"""
        assert self.omegas is not None,\
            """instance does not have omega information"""
        if iFrame is None:
            iFrame = self.iFrame
        if hasattr(iFrame, '__len__'):
            'take care of case nframes>1 in last call to read'
            retval = num.mean(self.omegas[iFrame])
        else:
            retval = self.omegas[iFrame]
        return retval

    def __readNext(self, nskip=0):
        if self.img is None:
            raise RuntimeError, 'no image file open'
        
        if nskip > 0:
            self.img.seek(self.nbytesFrame*nskip, 1)
        data = num.fromfile(self.img, 
                            dtype=self.dtypeRead,
                            count=self.nrows*self.ncols)
        return data

    def getNFrames(self, lessEmpty=True):
        fileBytes = os.stat(self.filename).st_size
        nFrames = getNFramesFromBytes(fileBytes, self.__nbytes_header, self.nbytesFrame)
        if lessEmpty:
            nFrames -= self.__nempty
        return nFrames

    def getOmegaMinMax(self):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaMin, self.omegaMax
    def getDeltaOmega(self, nframes=1):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaDelta * nframes
    def getDark(self):
        'no dark yet supported'
        return 0
    def getFrameOmega(self, iFrame=None):
        """if iFrame is none, use internal counter"""
        assert self.omegas is not None,\
            """instance does not have omega information"""
        if iFrame is None:
            iFrame = self.iFrame
        if hasattr(iFrame, '__len__'):
            'take care of case nframes>1 in last call to read'
            retval = num.mean(self.omegas[iFrame])
        else:
            retval = self.omegas[iFrame]
        return retval

    def getWriter(self, filename):
        # if not self.doFlip is False:
        #     raise NotImplementedError, 'doFlip true not coded'
        new = FrameWriter(self.ncols, self.nrows, 
                          filename=filename,
                          dtypeDefault=self.dtypeDefault, 
                          dtypeRead=self.dtypeRead, 
                          dtypeFloat=self.dtypeFloat, 
                          nbytesHeader=self.__nbytes_header) 
        return new
    
class ReadGE(Framer2DRC):
    """
    Read in raw GE files; this is the class version of the foregoing functions

    NOTES

    *) The flip axis ('v'ertical) was verified on 06 March 2009 by
       JVB and UL.  This should be rechecked if the configuration of the GE
       changes or you are unsure.

    *) BE CAREFUL! nframes should be < 10 or so, or you will run out of
       memory in the namespace on a typical machine.

    *) The header is currently ignored

    *) If a dark is specified, this overrides the use of empty frames as
       background; dark can be a file name or frame

    *) In multiframe images where background subtraction is requested but no
       dark is specified, attempts to use the
       empty frame(s).	An error is returned if there are not any specified.
       If there are multiple empty frames, the average is used.

    """
    """
    It is likely that some of the methods here should be moved up to a base class
    """
    __nbytes_header    = 8192
    __idim             = 2048
    __nrows = __idim
    __ncols = __idim
    __frame_dtype_dflt = 'int16' # good for doing subtractions
    __frame_dtype_read = 'uint16'
    __frame_dtype_float = 'float64'
    __nbytes_frame     = num.nbytes[num.uint16]*__nrows*__ncols # = 2*__nrows*__ncols
    __debug = False
    __useThreading = True and haveThreading
    __location = '  ReadGE'
    __readArgs = {
        'dtype' : __frame_dtype_read, 
        'count' : __nrows*__ncols
        }
    __castArgs = {
        'dtype' : __frame_dtype_dflt
        }
    __inParmDict = {
        'omegaStart':None,
        'omegaDelta':None,
        'subtractDark':False,
        'mask':None,
        'useMask':None,
        'dark':None,
        'dead':None,
        'nDarkFrames':1,
        'doFlip':True,
        'flipArg':'v',
        }
    # 'readHeader':False
    def __init__(self,
                 fileInfo,
                 *args,
                 **kwargs):
        """
        meant for reading a series of frames from an omega sweep, with fixed delta-omega
        for each frame

        omegaStart and omegaDelta can follow fileInfo or be specified in whatever order by keyword

        fileInfo: string, (string, nempty), or list of (string, nempty) for multiple files

        for multiple files and no dark, dark is formed only from empty
        frames in the first file
        """

        Framer2DRC.__init__(self, 
                            self.__nrows, self.__ncols,
                            dtypeDefault = self.__frame_dtype_dflt, 
                            dtypeRead    = self.__frame_dtype_read,
                            dtypeFloat   = self.__frame_dtype_float,
                            )
        
        # defaults
        self.__kwPassed = {}
        for parm, val in self.__inParmDict.iteritems():
            self.__kwPassed[parm] = kwargs.has_key(parm)
            if kwargs.has_key(parm):
                val = kwargs.pop(parm)
            self.__setattr__(parm, val)
        if len(kwargs) > 0:
            raise RuntimeError, 'unparsed keyword arguments: '+str(kwargs.keys())
        if len(args) == 0:
            pass
        elif len(args) == 2:
            self.omegaStart = args[0]
            self.omegaDelta = args[1]
        else:
            raise RuntimeError, 'do not know what to do with args : '+str(args)

        # initialization
        self.omegas = None
        self.img = None
        self.th  = None
        self.fileInfo      = None
        self.fileInfoR     = None
        self.nFramesRemain = None # remaining in current file
        self.iFrame = -1 # counter for last global frame that was read

        if self.dark is not None:
            if not self.__kwPassed['subtractDark']:
                'subtractDark was not explicitly passed, set it True'
                self.subtractDark = True
            if isinstance(self.dark, str):
                darkFile = self.dark
                self.dark = ReadGE.readDark(darkFile, nframes=self.nDarkFrames)
                self.__log('got dark from %d frames in file %s' % (self.nDarkFrames, darkFile))
            elif isinstance(self.dark, num.ndarray):
                assert self.dark.size == self.__nrows * self.__ncols, \
                    'self.dark wrong size'
                self.dark.shape = (self.__nrows, self.__ncols)
                if self.dark.dtype.name == self.__frame_dtype_read:
                    'protect against unsigned-badness when subtracting'
                    self.dark = self.dark.astype(self.__frame_dtype_dflt)
                self.__log('got dark from ndarray input')
            else:
                raise RuntimeError, 'do not know what to do with dark of type : '+str(type(self.dark))
        
        if fileInfo is not None:
            self.__setupRead(fileInfo, self.subtractDark, self.mask, self.omegaStart, self.omegaDelta)

        return

    # property:  useThreading

    @property
    def useThreading(self):
        """turn threading on or off"""
        return self.__useThreading
    
    @useThreading.setter
    def useThreading(self, v):
        """Set method for useThreading"""
        self.__useThreading = haveThreading and v
        return

    @classmethod
    def display(cls, 
                thisframe,
                roi = None,
                pw  = None,
                **kwargs
                ):
        'this is a bit ugly in that it sidesteps the dtypeRead property'
        retval = Framer2DRC.display(thisframe, roi=roi, pw=pw, dtypeRead=cls.__frame_dtype_read)
        return retval
    
    @classmethod
    def readRaw(cls, fname, mode='raw', headerlen=0):
        '''
        read a raw binary file;
        if specified, headerlen is in bytes;
        does not do any flipping
        '''
        print cls
        if hasattr(cls, 'doFlip'):
            print 'has doFlip'
        img = open(fname, mode='rb')
        if headerlen > 0:
            img.seek(headerlen, 0)
        if mode == 'raw' or mode == 'avg':
            dtype = cls.__frame_dtype_read
        elif mode == 'sum':
            dtype = 'float32'
        else:
            raise RuntimeError, 'unknown mode : '+str(mode)
        thisframe = num.fromfile(img, dtype=dtype, count=cls.__nrows*cls.__ncols).reshape(cls.__nrows, cls.__ncols)
        return thisframe
    def rawRead(self, *args, **kwargs):
        '''
        wrapper around readRaw that does the same flipping as the reader instance from which it is called
        '''
        thisframe = self.__flip(self.readRaw(*args, **kwargs))
        return thisframe
    @classmethod
    def readDark(cls, darkFile, nframes=1):
        'dark subtraction is done before flipping, so do not flip when reading either'
        darkReader = ReadGE(darkFile, doFlip=False)
        dark = darkReader.read(nframes=nframes, sumImg=True).astype(cls.__frame_dtype_dflt)
        darkReader.close()
        return dark
    def makeNew(self):
        """return a clean instance for the same data files
        useful if want to start reading from the beginning"""
        inParmDict = {}
        inParmDict.update(self.__inParmDict)
        for key in self.__inParmDict.keys():
            inParmDict[key] = eval("self."+key)
        newSelf = self.__class__(self.fileInfo, **inParmDict)
        return newSelf
    def getRawReader(self, doFlip=False):
        new = self.__class__(self.fileInfo, doFlip=doFlip)
        return new

    def get_nbytes_header(self):
        return self.__nbytes_header
    nbytesHeader = property(get_nbytes_header, None, None)

    def getWriter(self, filename):
        if not self.doFlip is False:
            raise NotImplementedError, 'doFlip true not coded'
        new = FrameWriter(self.ncols, self.nrows, 
                          filename=filename,
                          dtypeDefault=self.dtypeDefault, 
                          dtypeRead=self.dtypeRead, 
                          dtypeFloat=self.dtypeFloat, 
                          nbytesHeader=self.nbytesHeader) 
        return new
    
    def __setupRead(self, fileInfo, subtractDark, mask, omegaStart, omegaDelta):

        self.fileInfo = fileInfo
        self.fileListR = self.__convertFileInfo(self.fileInfo)
        self.fileListR.reverse() # so that pop reads in order

        self.subtractDark = subtractDark
        self.mask         = mask

        if self.dead is not None:
            self.deadFlipped = self.__flip(self.dead)

        assert (omegaStart is None) == (omegaDelta is None),\
            'must provide either both or neither of omega start and delta'
        if omegaStart is not None:
            if hasattr(omegaStart, 'getVal'):
                omegaStart = omegaStart.getVal('radians')
            if hasattr(omegaDelta, 'getVal'):
                omegaDelta = omegaDelta.getVal('radians')
            nFramesTot = self.getNFrames()
            self.omegas = \
                num.arange(omegaStart, omegaStart+omegaDelta*(nFramesTot-0.5), omegaDelta) + \
                0.5 * omegaDelta # put omegas at mid-points of omega range for frame
            omegaEnd = omegaStart+omegaDelta*(nFramesTot)
            self.omegaMin = min(omegaStart, omegaEnd)
            self.omegaMax = max(omegaStart, omegaEnd)
            self.omegaDelta = omegaDelta
            self.omegaStart = omegaStart

        self.__nextFile()

        return
    def getNFrames(self):
        """number of total frames with real data, not number remaining"""
        nFramesTot = self.getNFramesFromFileInfo(self.fileInfo)
        return nFramesTot
    def getDeltaOmega(self, nframes=1):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaDelta * nframes
    def getOmegaMinMax(self):
        assert self.omegas is not None,\
            """instance does not have omega information"""
        return self.omegaMin, self.omegaMax
    def frameToOmega(self, frame):
        scalar = num.isscalar(frame)
        frames = num.asarray(frame)
        if frames.dtype == int:
            retval = self.omegas[frames]
        else:
            retval = (frames + 0.5) * self.omegaDelta + self.omegaStart
        if scalar:
            retval = num.asscalar(retval)
        return retval
    def getFrameOmega(self, iFrame=None):
        """if iFrame is none, use internal counter"""
        assert self.omegas is not None,\
            """instance does not have omega information"""
        if iFrame is None:
            iFrame = self.iFrame
        if hasattr(iFrame, '__len__'):
            'take care of case nframes>1 in last call to read'
            retval = num.mean(self.omegas[iFrame])
        else:
            retval = self.omegas[iFrame]
        return retval
    def omegaToFrameRange(self, omega):
        assert self.omegas is not None,\
            'instance does not have omega information'
        assert self.omegaDelta is not None,\
            'instance does not have omega information'
        retval = omeToFrameRange(omega, self.omegas, self.omegaDelta)
        return retval
    def omegaToFrame(self, omega, float=False):
        assert self.omegas is not None,\
            'instance does not have omega information'
        if float:
            assert omega >= self.omegaMin and omega <= self.omegaMax,\
                'omega %g is outside of the range [%g,%g] for the reader' % (omega, self.omegaMin, self.omegaMax)
            retval = (omega - self.omegaStart)/self.omegaDelta - 0.5*self.omegaDelta
        else:
            temp = num.where(self.omegas == omega)[0]
            assert len(temp) == 1, 'omega not found, or found more than once'
            retval = temp[0]
        return retval
    def getFrameUseMask(self):
        """this is an optional toggle to turn the mask on/off"""
        assert isinstance(self.iFrame, int), \
            'self.iFrame needs to be an int for calls to getFrameUseMask'
        if self.useMask is None:
            retval = True
        else:
            assert len(self.useMask) == self.getNFrames(),\
                   "len(useMask) must be %d; yours is %d" % (self.getNFrames(), len(self.useMask))
            retval = self.useMask[self.iFrame]
        return retval
    @classmethod
    def __getNFrames(cls, fileBytes):
        retval = getNFramesFromBytes(fileBytes, cls.__nbytes_header, cls.__nbytes_frame)
        return retval
    def __nextFile(self):

        # close in case already have a file going
        self.close()

        fname, nempty = self.fileListR.pop()

        # open file
        fileBytes = os.stat(fname).st_size
        self.img = open(fname, mode='rb')

        # skip header for now
        self.img.seek(self.__nbytes_header, 0)

        # figure out number of frames
        self.nFramesRemain = self.__getNFrames(fileBytes)

        if nempty > 0:  # 1 or more empty frames
            if self.dark is None:
                scale = 1.0 / nempty
                self.dark = self.frame(dtype=self.__frame_dtype_float)
                for i in range(nempty):
                    self.dark = self.dark + num.fromfile(
                        self.img, **self.__readArgs
                        ).reshape(self.__nrows, self.__ncols) * scale
                self.dark.astype(self.__frame_dtype_dflt)
                self.__log('got dark from %d empty frames in file %s' % (nempty, fname))
            else:
                self.img.seek(self.nbytesFrame*nempty, 1)
            self.nFramesRemain -= nempty

        if self.subtractDark and self.dark is None:
            raise RuntimeError, "Requested dark field subtraction, but no file or empty frames specified!"

        return
    @staticmethod
    def __convertFileInfo(fileInfo):
        if isinstance(fileInfo,str):
            fileList = [(fileInfo, 0)]
        elif hasattr(fileInfo,'__len__'):
            assert len(fileInfo) > 0, 'length zero'
            if hasattr(fileInfo[0],'__iter__'): # checking __len__ bad because has len attribute
                fileList = copy.copy(fileInfo)
            else:
                assert len(fileInfo) == 2, 'bad file info'
                fileList = [fileInfo]
        else:
            raise RuntimeError, 'do not know what to do with fileInfo '+str(fileInfo)
        # fileList.reverse()
        return fileList
    def readBBox(self, bbox, raw=True, doFlip=None):
        """
        with raw=True, read more or less raw data, with bbox = [(iLo,iHi),(jLo,jHi),(fLo,fHi)]

        careful: if raw is True, must set doFlip if want frames
        potentially flipped; can set it to a reader instance to pull
        the doFlip value from that instance
        """

        if raw:
            if hasattr(doFlip,'doFlip'):
                'probably a ReadGe instance, pull doFlip from it'
                doFlip = doFlip.doFlip
            doFlip = doFlip or False # set to False if is None
            reader = self.getRawReader(doFlip=doFlip)
        else:
            assert doFlip is None, 'do not specify doFlip if raw is True'
            reader = self.makeNew()

        nskip = bbox[2][0]
        bBox = num.array(bbox)
        sl_i = slice(*bBox[0])
        sl_j = slice(*bBox[1])
        'plenty of performance optimization might be possible here'
        if raw:
            retval = num.empty( tuple(bBox[:,1] - bBox[:,0]), dtype=self.__frame_dtype_read )
        else:
            retval = num.empty( tuple(bBox[:,1] - bBox[:,0]), dtype=self.__frame_dtype_dflt )
        for iFrame in range(retval.shape[2]):
            thisframe = reader.read(nskip=nskip)
            nskip = 0
            retval[:,:,iFrame] = copy.deepcopy(thisframe[sl_i, sl_j])
        if not raw and self.dead is not None:
            'careful: have already flipped, so need deadFlipped instead of dead here'
            mask = num.tile(self.deadFlipped[sl_i, sl_j].T, (retval.shape[2],1,1)).T
            retval = num.ma.masked_array(retval, mask, hard_mask=True, copy=False)
        return retval
    def __flip(self, thisframe):
        if self.doFlip:
            if self.flipArg == 'v':
                thisframe = thisframe[:, ::-1]
            elif self.flipArg == 'h':
                thisframe = thisframe[::-1, :]
            elif self.flipArg == 'vh' or self.flipArg == 'hv':
                thisframe = thisframe[::-1, ::-1]
            elif self.flipArg == 'cw90':
                thisframe = thisframe.T[:, ::-1]
            elif self.flipArg == 'ccw90':
                thisframe = thisframe.T[::-1, :]
            else:
                raise RuntimeError, "unrecognized flip token."
        return thisframe
    def getDark(self):
        if self.dark is None:
            retval = 0
        else:
            retval = self.dark
        return retval
    def read(self, nskip=0, nframes=1, sumImg=False):
        """
        sumImg can be set to True or to something like numpy.maximum
        """

        'get iFrame ready for how it is used here'
        self.iFrame = num.atleast_1d(self.iFrame)[-1]
        iFrameList = []
        multiframe = nframes > 1

        nFramesInv = 1.0 / nframes
        doDarkSub = self.subtractDark # and self.dark is not None

        if doDarkSub:
            assert self.dark is not None, 'self.dark is None'

        # assign storage array
        if sumImg:
            sumImgCallable = hasattr(sumImg,'__call__')
            imgOut = self.frame(dtype=self.__frame_dtype_float, mask=self.dead)
        elif multiframe:
            imgOut = self.frame(nframes=nframes, dtype=self.__frame_dtype_dflt, mask=self.dead)


        # now read data frames
        for i in range(nframes):

            #data = self.__readNext(nskip=nskip)
            #thisframe = data.reshape(self.__nrows, self.__ncols)
            data = self.__readNext(nskip=nskip) # .reshape(self.__nrows, self.__ncols)
            self.iFrame += nskip + 1
            nskip=0 # all done skipping once have the first frame!
            iFrameList.append(self.iFrame)
            # dark subtraction
            if doDarkSub:
                'used to have self.__frame_dtype_float here, but self.__frame_dtype_dflt does the trick'
                thisframe = self.frame(buffer=data,
                                       dtype=self.__frame_dtype_dflt, mask=self.dead) - self.dark
            else:
                thisframe = self.frame(buffer=data,
                                       mask=self.dead)

            # flipping
            thisframe = self.__flip(thisframe)

            # masking (True get zeroed)
            if self.mask is not None:
                if self.getFrameUseMask():
                    thisframe[self.mask] = 0

            # assign output
            if sumImg:
                if sumImgCallable:
                    imgOut = sumImg(imgOut, thisframe)
                else:
                    imgOut = imgOut + thisframe * nFramesInv
            elif multiframe:
                imgOut[i, :, :] = thisframe[:, :]
        'end of loop over nframes'

        if sumImg:
            # imgOut = imgOut / nframes # now taken care of above
            pass
        elif not multiframe:
            imgOut = thisframe

        if multiframe:
            'make iFrame a list so that omega or whatever can be averaged appropriately'
            self.iFrame = iFrameList
        return imgOut
    def __log(self, message):
        if self.__debug:
            print self.__location+' : '+message
        return
    def __thWait(self):
        if self.__useThreading:
            if self.th is not None:
                if self.th.isAlive():
                    self.__log('wait for existing thread to finish')
                    tic = time.time()
                    self.th.join()
                    toc = time.time(); dt = toc - tic;
                    self.__log('--- existing thread has finished (%g seconds)' % (dt))
                else:
                    self.__log('existing thread already finished')
            # if done:
            #    self.th = None
        return
    def __thCheck(self):
        data = None
        if self.__useThreading:
            self.__thWait()
            if self.th is not None:
                if not self.th.success:
                    raise RuntimeError, 'failed to get image data'
                data = self.th.data
                self.nFramesRemain -= 1
                self.th = None
        return data
    def __readNext(self, nskip=0):

        if self.img is None:
            raise RuntimeError, 'no image file set'

        nHave = 0
        if self.__useThreading:
            data = self.__thCheck()

        nskipThis = nskip
        if nskipThis > 0 and data is not None:
            nskipThis = nskipThis - 1
            data = None
        if data is not None : nHave = 1
        #
        while self.nFramesRemain+nHave - nskipThis < 1:
            'not enough frames left in this file'
            nskipThis = nskipThis - self.nFramesRemain
            self.nFramesRemain = 0 # = self.nFramesRemain - self.nFramesRemain
            self.__nextFile()
        if nskipThis > 0:
            # advance counter past empty frames
            self.img.seek(self.nbytesFrame*nskipThis, 1)
            self.nFramesRemain -= nskipThis

        if data is None:
            # grab current frame
            data = num.fromfile(self.img, **self.__readArgs)
            data = num.array(data, **self.__castArgs)
            self.nFramesRemain -= 1

        if self.__useThreading:
            'try to get next frame'
            if self.nFramesRemain < 1:
                try:
                    self.__nextFile()
                except:
                    'if len((self.fileListR) == 0 then at end of files'
                    assert len(self.fileListR) == 0, \
                        'problem opening next file'
                    self.img = None
            if self.img is not None:
                self.th = ThreadReadFrame(self.img, self.__readArgs, self.__castArgs)
                self.th.start()
        return data
    def __call__(self, *args, **kwargs):
        return self.read(*args, **kwargs)
    def close(self):
        # if already have a file going, close it out
        self.__thWait()
        if self.img is not None:
            self.img.close()
        return
    """
    getReadDtype function replaced by dtypeRead property
    """
    @classmethod
    def maxVal(cls):
        'maximum value that can be stored in the image pixel data type'
        # dtype = reader._ReadGE__frame_dtype
        # maxInt = num.iinfo(cls.__frame_dtype_read).max # bigger than it really is
        maxInt = 2 ** 14
        return maxInt
    @classmethod
    def getNFramesFromFileInfo(cls, fileInfo, lessEmpty=True):
        fileList = cls.__convertFileInfo(fileInfo)
        nFramesTot = 0
        for fname, nempty in fileList:
            fileBytes = os.stat(fname).st_size
            nFrames = cls.__getNFrames(fileBytes)
            if lessEmpty:
                nFrames -= nempty
            nFramesTot += nFrames
        return nFramesTot

    def indicesToMask(self, indices):
      """
      Indices can be a list of indices, as from makeIndicesTThRanges
      """
      mask = self.getEmptyMask()
      if hasattr(indices,'__len__'):
        for indThese in indices:
          mask[indThese] = True
      else:
        mask[indices] = True
      return mask

class ReadMar165(Framer2DRC):
    """
    placeholder; not yet really implemented

    """
    __frame_dtype_read = 'uint16'
    __frame_dtype_dflt = 'int16' # good for doing subtractions
    def __init__(self, mode):
        if not isinstance(mode, int) or not [1,2,4,8].count(mode):
            raise RuntimeError, 'unknown mode : '+str(mode)

        self.__mode = mode
        self.__idim = mar165IDim(mode)
        return
    def __call__(self, filename):
        if not haveImageModule:
            msg = "PIL Image module is required for this operation, "\
                "but not loaded\n"
            raise NameError(msg)

        i = Image.open(filename, mode='r')
        a = num.array(i, dtype=self.__frame_dtype_read)
        frame = num.array(a, dtype=self.__frame_dtype_dflt)
        return frame


class ReadMar165NB1(ReadMar165):
    def __init__(self, *args, **kwargs):
        ReadMar165.__init__(self, 1, *args, **kwargs)
        return
class ReadMar165NB2(ReadMar165):
    def __init__(self, *args, **kwargs):
        ReadMar165.__init__(self, 2, *args, **kwargs)
        return
class ReadMar165NB3(ReadMar165):
    def __init__(self, *args, **kwargs):
        ReadMar165.__init__(self, 3, *args, **kwargs)
        return
class ReadMar165NB4(ReadMar165):
    def __init__(self, *args, **kwargs):
        ReadMar165.__init__(self, 4, *args, **kwargs)
        return

class LineStyles:
    """
    do not want to just cycle through default plot line colors, as end up with black lines
    """
    def __init__(self, lt=None):
        self.iStyle = 0
        if lt is None:
            lt = '-'
        self.styleList = ['r'+lt,'g'+lt,'b'+lt,'c'+lt,'m'+lt,'y'+lt]
        self.nStyle = len(self.styleList)
        return
    def __call__(self):
        style = self.styleList[self.iStyle % self.nStyle]
        self.iStyle += 1
        return style

funcTypeDflt = 'pv'

class Peak1DAtLoc:
  """
  base class for 1D peak shapes at fixed location;
  fixed that is unless newCenter is passed to the __call__ method
  """
  def __init__(self, centers, xVecDflt=None):
    """
    If __init__ is called with a list, then put one peak at each location
    """
    self.setCenters(centers)
    self.bRef = num.mean(self.centers) # does not change even if centers changed later
    self.setXVecDflt(xVecDflt)
    return
  def setXVecDflt(self, xVecDflt):
      if hasattr(xVecDflt, '__len__'):
          assert len(xVecDflt) == self.getNParams(),\
              'xVecDflt wrong length'
          self.xVecDflt = copy.deepcopy(xVecDflt)
      else:
          self.xVecDflt = xVecDflt
      return
  def getNParams(self):
    raise NotImplementedError
  def guessXVec(self, xs, vals, w=None):
      raise NotImplementedError
  def getNPeaks(self):
    return len(self.centers)
  def setCenters(self, centers):
    if hasattr(centers, '__len__'):
      self.centers = centers
    else:
      self.centers = [centers]
    return
  def __call__(self, xVec, p):
      # if newCenters is not None:
      #     self.setCenters(newCenters)
      if xVec is None or len(xVec) == 0:
          assert self.xVecDflt is not None,\
              'xVec is empty and xVecDflt is None'
          retval = self.eval(self.xVecDflt, p)
      else:
          retval = self.eval(xVec, p)
      return retval
  def eval(self, xVec, p):
    """
    xVec is parameters, p is positions
    """
    raise NotImplementedError
  def d_dx(self, xVec, p):
    'derivative of call with respect to xVec'
    raise NotImplementedError
  def d_dp(self, xVec, p):
    """derivative of call with respect to p
    assuming each eval depends only on its own point!
    """
    raise NotImplementedError
  def d_dCenters(self, xVec, p):
    'derivative of call with respect to centers'
    raise NotImplementedError
  def fitFloatingCenter(self, tThVals, intensityVals,
                        xVecGuess=None, centersGuess=None,
                        weights=4, tThWidth=None,
                        fitGoodnessTol=0.5):
      '''
      Note that centers are kept as they are -- if you want to
      actually change the centers of the function you need to call
      setCenters(cFit) after calling this function
      '''

      func = self

      if isinstance(weights, int):
          'interpret weights as number of quadrature points, need tThWidth'
          assert tThWidth is not None, \
              'need tThWidht for doing internal quadrature'
          assert len(tThVals.shape) == 1,\
              'tThVals wrong shape : '+str(tThVals.shape)
          quadr = weights
          xi1d, w1d = q1db.qLoc(quadr)
          'xi are in [0,1], so need to be centered'
          xi1d = xi1d - 0.5
          'and now make the right width for tTh bin sizes'
          xi1d = xi1d * tThWidth
          tThQP = num.tile(tThVals, (len(w1d),1)).T + num.tile(xi1d[:], (len(tThVals),1))
          w2d = num.tile(w1d, (len(tThVals),1))
      elif hasattr(weights, 'shape'):
          assert len(tThVals.shape) == 2,\
              'tThVals wrong shape : '+str(tThVals.shape)
          assert len(weights.shape) == 2,\
              'weights wrong shape'
          w2d = weights
          tThQP = tThVals
      else:
          raise RuntimeError, 'do not know what to do with weights of type '+str(type(weights))

      nPeaks = func.getNPeaks()

      if xVecGuess is None:
          #xVecGuess = func.xVecDflt
          xVecGuess = func.guessXVec(tThQP, intensityVals, w=w2d)
      if centersGuess is None:
          centersGuess = func.centers
      assert nPeaks == len(centersGuess), 'failed sanity check'
      xVec0 = num.hstack( (xVecGuess, centersGuess) )

      'keep centers so that can restore them later'
      centersRef = copy.deepcopy(func.centers)

      if len(xVec0) > len(intensityVals):
          raise RuntimeError, 'more DOF than equations'

      def _f_leastsq(x):
          # retval = num.empty(nBins)
          xFit = x[:-nPeaks]
          cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
          func.setCenters(cFit)
          # evalIntensity = func(xFit, tThVals)
          funcQP = func(xFit, tThQP)
          evalIntensity = num.sum(funcQP * w2d, axis=1)
          retval = evalIntensity - intensityVals
          return retval
      def _df_leastsq(x):
          # retval = num.empty(nBins)
          xFit = x[:-nPeaks]
          cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
          func.setCenters(cFit)
          retval = num.empty((len(intensityVals),len(x)))
          # funcQP = func(xFit, tThQP)
          d_evalQP_d_x       = func.d_dx(xFit, tThQP)
          d_evalQP_d_centers = func.d_dCenters(xFit, tThQP)
          for iX in range(0,len(xFit)):
              retval[:,iX]           = num.sum(d_evalQP_d_x[:,:,iX] * w2d, axis=1)
          for iC in range(0,len(cFit)):
              retval[:,len(xFit)+iC] = num.sum(d_evalQP_d_centers[:,:,iC] * w2d, axis=1)
          return retval

      x, ier = \
          leastsq(_f_leastsq, xVec0, Dfun=_df_leastsq)
      if not [1,2,3,4].count(ier):
          print >> sys.stderr, \
              'error code %d from leastsq' % (ier)
          raise RuntimeError, 'error from leastsq'

      xFit = x[:-nPeaks]
      cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
      'return evalIntensity so caller does not have to monkey with changing centers back and forth'
      func.setCenters(cFit)
      residual = _f_leastsq(x)
      residNorm = num.linalg.norm(residual, ord=num.inf)
      fitTol = num.linalg.norm(intensityVals)*fitGoodnessTol
      if residNorm > fitTol :
          print 'fit might not be good enough : %g > %g' % (residNorm, fitTol)
          # raise RuntimeError, 'fit is not good enough : %g > %g' % (residNorm, fitTol)
      evalIntensity = residual + intensityVals # func(xFit, tThVals)
      func.setCenters(centersRef)
      return xFit, cFit, evalIntensity

class PeakPV1DAtLoc(Peak1DAtLoc):
    """
    the pseudo-Voigt:
    f = A*( n*fl + (1 - n)*fg )
    """
    def __init__(self, *args, **kwargs):
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        self.C0 = 4.0*math.log(2.0)
        self.C1 = 4.0
        self.Al = num.sqrt(self.C1)/num.pi
        self.Ag = num.sqrt(self.C0/num.pi)
        return
    def getNParams(self):
        nParams = 2 + len(self.centers) * 3
        return nParams
    def guessXVec(self, xs, vals, w=None):
        xVec = num.empty(self.getNParams())
        if len(self.centers) == 1:
            xGauss = getGaussNDParams([xs], v=vals, w=w)
            xVec[0] = xGauss[3]
            xVec[1] = 0.0e0
            iPeak = 0
            H = xGauss[1] # FWHM
            xVec[2+iPeak*3]   = xGauss[2] * H / self.Ag
            xVec[2+iPeak*3+1] = H
            xVec[2+iPeak*3+2] = 0.0 # xn
        else:
            maxV = vals.max()
            minV = vals.min()
            width = (xs.max()-xs.min())/len(self.centers)
            xVec[0] = minV
            xVec[1] = 0.0e0
            for iPeak in range(len(self.centers)):
                xVec[2+iPeak*3]   = maxV * width # A
                xVec[2+iPeak*3+1] = 0.25 * width # H
                xVec[2+iPeak*3+2] = 0.0 # xn
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            retval += A * (n*fl + (1.0 - n)*fg)
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        retval[:,0] = 1.0e0

        # dB = xVec[1]
        # retval += dB * (p - self.bRef)
        retval[:,1] = p - self.bRef

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            #####################
            # some intermediate partials
            dfl_dH = ( 2.0 * self.C1 * delx**2 * fl/(self.Al*H) - 1.0 ) * (fl / H)
            dfg_dH = ( 2.0 * self.C0 * delx**2 / H**2 - 1.0 ) * (fg / H)

            dfl_dx0 = 2.0 * self.C1 * delx * fl**2 / (self.Al * H)
            dfg_dx0 = 2.0 * self.C0 * delx * fg / H**2

            # mixing parm
            dn_dxn = 0.5 * (1.0/num.cosh(xn)**2)
            #####################

            # amplitude
            df_dA = (n*fl + (1.0 - n)*fg)

            # FWHM
            df_dH = A * ( n*dfl_dH + (1.0 - n)*dfg_dH )

            # mixing parm
            df_dxn = A * ( fl - fg ) * dn_dxn

            # assign jacobian vals
            retval[:,2+iPeak*3]   = df_dA
            retval[:,2+iPeak*3+1] = df_dH
            retval[:,2+iPeak*3+2] = df_dxn

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        dB = xVec[1]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        retval = dB

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            dfl_dp = -(2.0 * self.C1) * delx * fl**2 / (self.Al*H)
            dfg_dp = -(2.0 * self.C0) * delx * fg / H**2

            retval += A * ( n*dfl_dp + (1.0 - n)*dfg_dp )
        return retval
    def d_dCenters(self, xVec, p):
        ps = p.shape
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center).flatten()

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            dfl_dc = (2.0 * self.C1) * delx * fl**2 / (self.Al*H)
            dfg_dc = (2.0 * self.C0) * delx * fg / H**2

            retval[:,iPeak] = A * ( n*dfl_dc + (1.0 - n)*dfg_dc )

        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval

class PeakLorentzian1DAtLoc(Peak1DAtLoc):
    def __init__(self, *args, **kwargs):
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        return
    def getNParams(self):
        """2 parameters for background, 2 for intensity and width of each peak"""
        nParams = 2 + len(self.centers) * 2
        return nParams
    def guessXVec(self, xs, vals, w=None):
        # guessXVec(self, widths, mins, maxs)
        xVec = num.empty(self.getNParams())
        xVec[0] = num.min(vals) # num.mean(mins)
        xVec[1] = 0.0e0
        width = (xs.max()-xs.min())/len(self.centers)
        for iPeak in range(len(self.centers)):
            xVec[2+iPeak*2]   = num.max(vals) # maxs[iPeak]
            xVec[2+iPeak*2+1] = 0.24 * width # 0.25 * widths[iPeak]
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            dist = (p - center)
            retval += A * (0.5*w) / ( dist*dist + (0.5*w)**2  )
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        retval[:,0] = 1.0e0

        # dB = xVec[1]
        # retval += dB * (p - self.bRef)
        retval[:,1] = p - self.bRef

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            dist = (p - center)

            # retval += A * num.exp( self.expFact * ( dist * dist) )
            evalLtz = (0.5*w) / ( dist*dist + (0.5*w)**2  )
            retval[:,2+iPeak*2]   =  evalLtz
            retval[:,2+iPeak*2+1] = A * (evalLtz / w - evalLtz * evalLtz)

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        dB = xVec[1]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        retval = dB

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            dist = (p - center)
            evalLtz = (0.5*w) / ( dist*dist + (0.5*w)**2  )
            # retval += A * num.exp( self.expFact * ( dist * dist) )
            retval += -1.0 * A * 4.0 * evalLtz * evalLtz * dist * (1/w)
        return retval
    def d_dCenters(self, xVec, p):
        ps = p.shape
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            dist = (p - center).flatten()
            evalLtz = (0.5*w) / ( dist*dist + (0.5*w)**2  )
            retval[:,iPeak] = A * 4.0 * evalLtz * evalLtz * dist * (1/w)
        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval

class PeakGauss1DAtLoc(Peak1DAtLoc):

    def __init__(self, *args, **kwargs):
        linBG=True
        if kwargs.has_key('linBG'):
            linBG = kwargs.pop('linBG')
        self.linBG = linBG
        self.nPBG = 1
        if self.linBG: self.nPBG += 1
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        self.expFact  = -4.0 * math.log(2.0)
        return
    def getNParams(self):
        """2 parameters for background, 2 for intensity and width of each peak"""
        nParams = self.nPBG + len(self.centers) * 2
        return nParams
    def guessXVec(self, xs, vals, w=None):
        xVec = num.empty(self.getNParams())
        if len(self.centers) == 1:
            xGauss = getGaussNDParams([xs], v=vals, w=w)
            xVec[0] = xGauss[3]
            width = (xs.max()-xs.min())/len(self.centers)
            if self.linBG :
                xVec[1] = 0.0e0
            iPeak = 0
            xVec[self.nPBG+iPeak*2]   = xGauss[2]
            xVec[self.nPBG+iPeak*2+1] = xGauss[1]
        else:
            xVec[0] = num.min(vals)
            width = (xs.max()-xs.min())/len(self.centers)
            if self.linBG :
                xVec[1] = 0.0e0
            for iPeak in range(len(self.centers)):
                xVec[self.nPBG+iPeak*2]   = num.max(vals) # maxs[iPeak]
                xVec[self.nPBG+iPeak*2+1] = 0.25 * width # 0.25 * widths[iPeak]
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        if self.linBG:
            dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        if self.linBG:
            retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = (p - center) / w
            retval += A * num.exp( self.expFact * ( dist * dist) )
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        iXBase = 0
        retval[:,iXBase] = 1.0e0
        iXBase += 1

        if self.linBG:
            # dB = xVec[1]
            # retval += dB * (p - self.bRef)
            retval[:,iXBase] = p - self.bRef
            iXBase += 1

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = (p - center) / w

            # retval += A * num.exp( self.expFact * ( dist * dist) )
            evalExp = num.exp( self.expFact * ( dist * dist) )
            retval[:,iXBase+iPeak*2]   = evalExp
            retval[:,iXBase+iPeak*2+1] = A * evalExp * self.expFact * 2.0e0 * dist * (-1.0*dist/w)
        iXBase += len(self.centers)*2
        assert iXBase == len(xVec), 'bookkeeping error'

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        if self.linBG:
            dB = xVec[1]
            retval = dB

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = (p - center) / w
            # retval += A * num.exp( self.expFact * ( dist * dist) )
            retval += A * num.exp( self.expFact * ( dist * dist) ) * self.expFact * 2.0e0 * dist * (1/w)
        return retval
    def d_dCenters(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        #retval = num.zeros(num.hstack((p.shape,len(self.centers))))
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = ((p - center) / w).flatten()
            retval[:,iPeak] = -A * num.exp( self.expFact * ( dist * dist) ) * self.expFact * 2.0e0 * dist * (1/w)
        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval

class MultiRingBinned:
    """
    like MultiRingEval, but works with polar rebinned (or 'caked') images

    no funcXVecList to init because expectation is that always pulled from dataFrame

    should work fine whether or not corrected is True in polarRebinKWArgs;
    but note that default is changed to True

    if etaRange is not specified in polarRebinKWArgs (or is None),
    then the etaRange is calculated based on numEta so that the first
    eta bin is centered around an angle of zero

    note that this works by matching intesities for binned 1D
    functions in a least sqaures problem; one could probably instead
    form a residual on the two-theta values of the image frame
    positions for the peak centers found with independtly floating
    centers (on the output of the getTThErrors method)

    KEYWORD ARGS

    funcType = funcTypeDflt,
    refineParamsDG = True,
    refineParamsL = False,
    targetNRho = 30,
    polarRebinKWArgs = {},
    quadr = 4,
    npdivMax = 4,
    samplingFactor = 0.25,
    singleRebin = True,
    distortionFreeRefDG = False,
    log=None:  if not None, then a file-like object with a "write" method;

    """

    __debug = True
    __print = True
    __location = '  MultiRingBinned'
    def __init__(self, detectorGeom, planeData,
                 dataFrame, 
                 funcType = funcTypeDflt, 
                 refineParamsDG = True,
                 refineParamsL = False,
                 targetNRho = None,
                 polarRebinKWArgs = {},
                 quadr = 4,
                 npdivMax = 8,
                 samplingFactor = 1,
                 singleRebin = True,
                 distortionFreeRefDG = False,
                 log=None
                 ):

        if log is not None:
            self.logfile = log
        else:
            log = sys.stdout
        
        ticMethod = time.time()

        # if copyFrame:
        #     self.dataFrame = copy.deepcopy(dataFrame)
        # else:
        #     self.dataFrame = dataFrame
        #
        "defaults which want different from those for polarRebin's defaults"
        prbkw = {
            'corrected' : False,
            'ROI'       : None,
            'etaRange'  : None,
            'numEta'    : 36,
            'verbose'   : False,
            'npdiv'     : None,
            'log'       : self.logfile
            }
        "things from the user"
        prbkw.update(polarRebinKWArgs)
        "keep track of the actual values for subsequent manipulations"
        self.corrected = prbkw['corrected']
        self.ROI       = prbkw['ROI']
        self.etaRange  = prbkw['etaRange']
        self.numEta    = prbkw['numEta']
        if prbkw['etaRange'] is None:
            prbkw['etaRange'] = num.array([0.,2*num.pi]) - num.pi/self.numEta

        '''groupings of rings could change if planeData.lparms is changing, but
        need to keep fixed number of degrees of freedom so just use initial values'''
        iHKLLists = detectorGeom.makeTThRanges(planeData, cullDupl=True)
        tThs = planeData.getTTh()
        tThRanges = planeData.getTThRanges()
        iRingToHKL = reduce(lambda x, y: x+y, iHKLLists)
        iHKLToRing = -1 * num.ones(planeData.getNHKLs(), dtype=int)
        iHKLToRing[iRingToHKL] = num.arange(len(iRingToHKL))
        self.floatingCentersTTh = num.zeros([self.numEta, len(iRingToHKL)])
        self.floatingCentersIJ  = num.zeros([self.numEta, len(iRingToHKL), 2])

        'make a reference detector geom'
        self.refDG   = detectorGeom.makeNew()
        if distortionFreeRefDG:
            'make the reference detector geom be distortion free'
            self.refDG.setDParamZero()

        nParamCur = 0
        if refineParamsDG:
            nParamDG = detectorGeom.getNParams()
        else:
            nParamDG = 0
        self.xIndDG = range(nParamCur, nParamCur+nParamDG)
        nParamCur += nParamDG

        if refineParamsL:
            nParamL = len(planeData.lparms)
        else:
            nParamL = 0
        self.xIndL = range(nParamCur, nParamCur+nParamL)
        nParamCur += nParamL

        def computeBinning_l(iRingSet, npdivDflt):
            iHKLList = iHKLLists[iRingSet]
            tThMin   = num.min(tThRanges[iHKLList,0])
            tThMax   = num.max(tThRanges[iHKLList,1])
            nPeaks   = len(iHKLList)

            rhoMin     = num.tan(tThMin) * self.refDG.workDist
            rhoMax     = num.tan(tThMax) * self.refDG.workDist
            rhoPxRange = num.array([rhoMin, rhoMax])/self.refDG.pixelPitch
            rhoPxEff   = int((rhoPxRange[1] - rhoPxRange[0]) / samplingFactor)
            npdiv      = npdivDflt

            if targetNRho is None:
                numRho = rhoPxEff
            else:
                numRho = targetNRho*nPeaks

            if numRho <= rhoPxEff:
                if npdiv == None:
                    npdiv = 1
            else:
                if npdiv == None:
                    "user did not set npdiv in polarRebinKWArgs, so assume we can do what we like"
                    npdivCalc = round(float(numRho)/rhoPxEff)
                    if npdivCalc > 2*npdivMax:
                        raise RuntimeError, 'npdiv wants to be %d, substantially more than max of %d' % (npdiv, npdivMax)
                    else:
                        npdiv = npdivCalc
                else:
                    if rhoPxEff < numRho * 0.5:
                        raise RuntimeError, 'rhoPx %d for %d peaks is less than half of target' % (rhoPx, nPeaks)

            return numRho, npdiv, rhoPxRange

        polImgSingle = None
        if singleRebin:
            prbkwThis = {}
            prbkwThis.update(prbkw)

            'use innermost ring set to figure out how many rho bins are wanted across the whole range'
            iRingSet = 0
            numRhoRS0, prbkwThis['npdiv'], rhoPxRangeRS0 = computeBinning_l(iRingSet, prbkw['npdiv'])

            tThMin     = num.min(tThRanges)
            tThMax     = num.max(tThRanges)
            rhoMin     = num.tan(tThMin) * self.refDG.workDist
            rhoMax     = num.tan(tThMax) * self.refDG.workDist
            rhoPxRange = num.array([rhoMin, rhoMax])/self.refDG.pixelPitch
            'and now compute number of rho bins across all ring sets'
            numRho     = int(numRhoRS0 * (rhoPxRange[1]-rhoPxRange[0]) / (rhoPxRangeRS0[1]-rhoPxRangeRS0[0]))
            if self.__print:
                self.__log(' for ring set %d have %d peak(s) and %d rho bins' % (iRingSet, len(iHKLLists[iRingSet]), numRhoRS0))
                self.__log(' for all rings have %d peak(s) and %d rho bins' % (len(iRingToHKL), numRho))

            prbkwThis['rhoRange'] = rhoPxRange; prbkwThis['numRho'] = numRho;
            if self.__print:
                self.__log('  doing polar rebin of frame data, with npdiv %d' % (prbkwThis['npdiv']))
                self.__log('  single polar rebin call keyword arguments : %s' % (str(prbkwThis)))
            tic = time.time()
            polImgSingle = self.refDG.polarRebin(dataFrame, **prbkwThis)
            if num.any(num.isnan(polImgSingle['intensity'])):
                raise RuntimeError, 'got NaNs in rebinned image'
            toc = time.time(); dt = toc - tic;
            if self.__print:
                self.__log(' single polar rebin call took %g seconds' % (dt))

        ringDataList = [] # funcDataList = []
        nPoints = 0
        dtPRTot = 0.e0
        for iRingSet in range(len(iHKLLists)):
            iHKLList = iHKLLists[iRingSet]
            nPeaks   = len(iHKLList)
            tThMin   = num.min(tThRanges[iHKLList,0])
            tThMax   = num.max(tThRanges[iHKLList,1])

            if singleRebin:
                polImg = polImgSingle

                'tth values for figuring out where the rings fall'
                tThVals      = num.arctan2(polImg['radius'], self.refDG.workDist)

                deltaRho = polImgSingle['deltaRho']
                rhoBinIndices = num.where(num.logical_and(
                        tThVals > tThMin, tThVals < tThMax
                        ))[0]

            else:
                prbkwThis = {}
                prbkwThis.update(prbkw)

                prbkwThis['numRho'], prbkwThis['npdiv'], prbkwThis['rhoRange'] = computeBinning_l(iRingSet, prbkw['npdiv'])
                if self.__print:
                    self.__log('  doing polar rebin of frame data, with npdiv %d' % (prbkwThis['npdiv']))
                    self.__log('  polar rebin call keyword arguments : %s' % (str(prbkwThis)))
                tic = time.time()
                polImg = self.refDG.polarRebin(dataFrame, **prbkwThis)
                if num.any(num.isnan(polImg['intensity'])):
                    raise RuntimeError, 'got NaNs in rebinned image'
                toc = time.time(); dt = toc - tic; dtPRTot += dt; 
                if self.__print:
                    self.__log(' polar rebin call took %g seconds' % (dt))

                deltaRho = polImg['deltaRho']
                rhoBinIndices = range(0, prbkwThis['numRho'])

            etaVals           = num.array(polImg['azimuth'])
            rhoVals           = num.array(polImg['radius'])
            dataFrameRebinned = num.array(polImg['intensity']) # (numEta, numRho)

            xi, w = q1db.qLoc(quadr)
            'xi are in [0,1], so need to be centered'
            xi = (xi - 0.5) * deltaRho # to cover rho bin size
            nQP = len(w)

            etaList = []
            for iEta, etaThis in enumerate(etaVals):

                xVecDflt = None
                if funcType is 'pv':
                    func = PeakPV1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
                elif funcType is 'gauss' or funcType is 'normal':
                    func = PeakGauss1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
                elif funcType is 'gaussFlatBG':
                    func = PeakGauss1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt, linBG=False)
                elif funcType is 'lorentz':
                    func = PeakLorentzian1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
                else:
                    raise RuntimeError,  'unknown funcType : '+str(funcType)

                # fold quadrature in to get pixel locations
                nPx = len(rhoBinIndices)
                etaQP = etaThis * num.ones((nPx, nQP))
                rhoQP = num.tile(rhoVals[rhoBinIndices], (nQP,1)).T + num.tile(xi, (nPx,1))
                # binEtas, binRhos = num.meshgrid(self.etaVals, self.rhoVals)
                # binEtas = binEtas.T; binRhos = binRhos.T;
                # nope # iQP, jQP = self.refDG.polarToPixel(rhoQP, etaQP, ROI=self.ROI, corrected=self.corrected)
                iQP, jQP = self.refDG.polarToPixel(rhoQP, etaQP, corrected=self.corrected)
                wQP = num.tile(w, (nPx,1))

                'leastsq to do fit, with floating center'
                binnedIntensity = dataFrameRebinned[iEta,rhoBinIndices]
                tThQP, etaQP_temp = self.refDG.xyoToAng(iQP, jQP)
                if self.__print:
                    self.__log('   fitting function for iRingSet %d iEta %d' % (iRingSet, iEta))
                xFit, cFit, evalIntensity = func.fitFloatingCenter(
                    tThQP, binnedIntensity, weights=wQP)
                func.setXVecDflt(xFit)
                'keep tThc and (i,j) enters for reference fits, in case want to know how far off they were'
                iCFit, jCFit = self.refDG.angToXYO(cFit, etaThis * num.ones(nPeaks))
                self.floatingCentersTTh[iEta, iHKLToRing[iHKLList]]    = cFit
                self.floatingCentersIJ [iEta, iHKLToRing[iHKLList], 0] = iCFit
                self.floatingCentersIJ [iEta, iHKLToRing[iHKLList], 1] = jCFit

                normalization = 1.0/num.mean(dataFrameRebinned[iEta,rhoBinIndices])

                nPoints += nPx
                etaList.append((func, iQP, jQP, wQP, rhoBinIndices, normalization))
            'done with iEta'
            ringDataList.append({
                    'etaFuncData'       : etaList,
                    'etaVals'           : etaVals,
                    'rhoVals'           : rhoVals,
                    'dataFrameRebinned' : dataFrameRebinned,
                 })
        'done with iRingSet'
        if not singleRebin:
            if self.__print:
                self.__log(' all together polar rebin calls took %g seconds' % (dtPRTot))


        'store stuff'
        self.detectorGeom = detectorGeom
        self.planeData    = planeData
        self.nParamDG     = nParamDG
        self.nParamL      = nParamL
        self.ringDataList = ringDataList
        self.nParam       = nParamCur
        self.nPoints      = nPoints
        self.iHKLLists    = iHKLLists
        self.iRingToHKL   = iRingToHKL
        self.iHKLToRing   = iHKLToRing

        assert self.nParam > 0, \
            'do not have any free parameters'

        # do not yet have analytic derivatives for detector, so do not bother with Dfun
        self.useDFun = False

        tocMethod = time.time(); dt = tocMethod - ticMethod;
        if self.__print:
            self.__log(' init method took %g seconds' % (dt))


        return

    # property:  logfile

    def _get_logfile(self):
        """Get method for logfile"""
        return self._logfile

    def _set_logfile(self, v):
        """Set method for logfile"""
        self._logfile = v
        return

    logfile = property(_get_logfile, _set_logfile, None,
                                "file for log messages")

    def __log(self, message):
        """Logging facility for this class"""
        if self.logfile:
            self.logfile.write(self.__location + ': ' + message + '\n')
            pass

        return

    def guessXVec(self):
        xVec = num.zeros(self.nParam)
        if self.nParamDG > 0:
            xVec[self.xIndDG] = self.detectorGeom.getParams()
        if self.nParamL > 0:
            xVec[self.xIndL]  = self.planeData.lparms
        return xVec
    def getNParam(self):
        return self.nParam
    def deval(self, xVec):
        raise NotImplementedError, 'deval not coded yet'
        return retval
    def eval(self, xVec):
        """
        careful: this updates the settings in detectorGeom and planeData
        """

        retval = num.empty(self.nPoints)

        'first update the detector geometry'
        if self.nParamDG > 0:
            self.detectorGeom.updateParams(*xVec[self.xIndDG])

        if self.nParamL > 0:
            self.planeData.lparms = xVec[self.xIndL] # property, calls set_lparms
            tThs = self.planeData.getTTh()
            for iRingSet in range(len(self.iHKLLists)):
                iHKLList = self.iHKLLists[iRingSet]
                for iEta, etaFuncData in enumerate(self.ringDataList[iRingSet]['etaFuncData']):
                    #func, iQP, jQP, wQP, rhoBinIndices, normalization = etaFuncData
                    func = etaFuncData[0]
                    func.setCenters(tThs[(iHKLList)])

        nPoints = 0
        for iRingSet in range(len(self.iHKLLists)):
            for iEta, etaFuncData in enumerate(self.ringDataList[iRingSet]['etaFuncData']):
                func, iQP, jQP, wQP, rhoBinIndices, normalization = etaFuncData
                nPointsThis = wQP.shape[0] # len(data)

                tThQP, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)
                evalQP       = func(None, tThQP)
                evalP        = num.sum(evalQP * wQP, axis=1)
                retval[nPoints:nPoints+nPointsThis] = evalP

                nPoints += nPointsThis

        return retval

    def __call__(self, xVec, makePlots=False, etaPlotIndices=None):
        'meant for call in leastsq type algorithm, see doFit'

        retval = num.empty(self.nPoints)

        evalAllRings = self.eval(xVec)

        if makePlots:
            if isinstance(etaPlotIndices, bool):
                plotEta = num.ones(self.numEta, dtype=bool)
                plotEta[:] = etaPlotIndices
            else:
                etaPlotIndices = num.atleast_1d(etaPlotIndices)
                plotEta = num.array(
                    [num.any(etaPlotIndices == iEta) for iEta in range(self.numEta) ],
                    dtype=bool)

            plotTitlePrefix = ''
            if not num.any(plotEta):
                plotTitlePrefix='eta with worse errors'
            plotWinRadial = self.makePlotWin(sqrtIntensity=True,
                                             plotTitlePrefix=plotTitlePrefix)

        nPoints = 0
        errsByEta = num.zeros(self.numEta)
        for iRingSet, ringData in enumerate(self.ringDataList): # range(len(self.iHKLLists)):
            dataFrameRebinned = ringData['dataFrameRebinned']
            for iEta, etaFuncData in enumerate(ringData['etaFuncData']):
                func, iQP, jQP, wQP, rhoBinIndices, normalization = etaFuncData
                nPointsThis = wQP.shape[0] # len(data)

                evalThisEta  = evalAllRings[nPoints:nPoints+nPointsThis]
                dataThisEta = dataFrameRebinned[iEta, rhoBinIndices]
                diff  = evalThisEta - dataThisEta
                retval[nPoints:nPoints+nPointsThis] = normalization * diff
                errsByEta[iEta] += num.linalg.norm(retval[nPoints:nPoints+nPointsThis])

                if makePlots and plotEta[iEta]:
                    self.plotByRingEta(iRingSet, iEta,
                                       sqrtIntensity=True, win=plotWinRadial)

                nPoints += nPointsThis
            'iEta'
        'iRingSet'
        if makePlots:
            'plot maximum error eta'
            iEta = num.argmax(errsByEta)
            if not plotEta[iEta]:
                for iRingSet, ringData in enumerate(self.ringDataList):
                    self.plotByRingEta(iRingSet, iEta,
                                       sqrtIntensity=True, win=plotWinRadial)

            plotWinRadial.show()
        if self.__print:
            self.__log(' norm : '+str(num.linalg.norm(retval)))
        return retval

    def makePlotWin(self, sqrtIntensity=True, plotTitlePrefix=''):
        plotWinRadial = plotwrap.PlotWin(2,1, relfigsize=(8,3),
                                         title=plotTitlePrefix+'ring evaluation results')
        'using twinx() might be better, but not plumbed for that'
        if sqrtIntensity:
            ylabel='sqrt(intensity)'
        else:
            ylabel='intensity'
        pwREval = plotwrap.PlotWrap(window=plotWinRadial, ylabel=ylabel,
                                    showByDefault=False)
        pwRDiff = plotwrap.PlotWrap(window=plotWinRadial, ylabel='relative error',
                                    showByDefault=False, axprops={'sharex':pwREval.a})
        return plotWinRadial
    def plotByRingEta(self, iRingSet, iEta, win=None, sqrtIntensity=True, alpha=0.25):
        'may have redundant work here, but assume this is not a big deal if doing plots'

        if win is None:
            win = self.makePlotWin(sqrtIntensity=sqrtIntensity)
            retval = win
        else:
            retval = None
        aEval, pwREval = win.getAxes(0, withPW=True)
        aDiff, pwRDiff = win.getAxes(1, withPW=True)

        etaFuncData = self.ringDataList[iRingSet]['etaFuncData'][iEta]
        func, iQP, jQP, wQP, rhoBinIndices, normalization = etaFuncData

        tThQP, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)
        evalQP       = func(None, tThQP)
        evalP        = num.sum(evalQP * wQP, axis=1)
        tThCen = num.sum(tThQP * wQP, axis=1)
        width = num.mean(tThCen[1:] - tThCen[:-1])

        dataThisEta = self.ringDataList[iRingSet]['dataFrameRebinned'][iEta, rhoBinIndices]

        if sqrtIntensity:
            def mapData(v):
                return num.sqrt(num.maximum(0, v))
        else:
            def mapData(v):
                return v

        pwREval.a.bar(tThCen, mapData(dataThisEta), width=width, align='center', alpha=alpha)
        '''have not set things up to match line styles with drawRings, so make all lines black;
        might be able to fix this if switch drawRings to do something like cullDupl too
        '''
        pwREval(tThCen, mapData(evalP), style='kx-')
        scaling = dataThisEta
        scaling[num.where(scaling==0)] = 1.0
        pwRDiff(tThCen, (evalP - dataThisEta) / scaling, style='kx-')

        if retval is not None:
            retval.show()
        return retval
    def doFit(self, xVec0=None, **lsKWArgs):
        """
        lsKWArgs can have things like ftol and xtol for leastsq
        """

        if xVec0 is None:
            xVec0 = self.guessXVec()
        tic = time.time()
        if self.useDFun:
            x, ier = leastsq(self, xVec0, Dfun=self.deval, **lsKWArgs)
        else:
            x, ier = leastsq(self, xVec0, **lsKWArgs)
        toc = time.time()
        dt = toc - tic
        self.__log('leastsq took %g seconds to run' % (dt))
        if not [1,2,3,4].count(ier):
            msg = 'error code %d from leastsq' % (ier)
            self.__log(msg)
            print >> sys.stderr, msg
            raise RuntimeError, 'error from leastsq'

        # must update detector geometry with solution
        fitParams = num.atleast_1d(x)
        pList = self.detectorGeom._Detector2DRC__makePList()
        pList[num.atleast_1d(self.detectorGeom.refineFlags)] = fitParams
        self.detectorGeom._Detector2DRC__updateFromPList(pList)

        return fitParams

    def getTThErrors(self, plot=False, units='strain', outputFile=None):
        """
        convenient way of looking at the errors, though not how the
        errors are actually measured in the fitting procedure;
        get the tTh values at the image frame locations deemed to be
        the centers with the floating-center fits

        units can be:  'mm' <radius>, 'd-spacing', 'strain', 'radians' <tTh>, 'degrees' <tTh>
        """
        # baseline ("ideal") tTh values
        tThs = self.planeData.getTTh()[self.iRingToHKL]
        
        # need these to reconstruct azimuths
        #   - Are in mm (self.detectorGeom.pixelPitchUnit)
        xFl, yFl = self.detectorGeom.cartesianCoordsOfPixelIndices(
            self.floatingCentersIJ[:,:,0], self.floatingCentersIJ[:,:,1])

        # azimuthal angles from xyoToAng
        tThFloating, etaFloating = self.detectorGeom.xyoToAng(
                self.floatingCentersIJ[:,:,0], self.floatingCentersIJ[:,:,1] )

        if units == 'mm':
            rho = num.sqrt((xFl - self.detectorGeom.xc)**2 +
                           (yFl - self.detectorGeom.yc)**2)
            x0, y0 = self.detectorGeom.angToXYO(num.tile(tThs, (self.numEta, 1)),
                                                etaFloating,
                                                units=self.detectorGeom.pixelPitchUnit)
            rho0 = num.sqrt((x0 - self.detectorGeom.xc)**2 +
                            (y0 - self.detectorGeom.yc)**2)
            ylabel = '(rho - rho0)'
            errs = rho - rho0
        elif units == 'd-spacing' or units == 'strain':
            dspFloating = self.planeData.wavelength / ( 2. * num.sin(0.5 * tThFloating) )
            dsp         = self.planeData.getPlaneSpacings()
            print len(dspFloating)
            print len(dsp)
            dspErrs     = dspFloating - num.tile(dsp, (self.numEta, 1) )
            if units == 'strain':
                ylabel = '(d - d0) / d0'
                errs = dspErrs / dsp
            else:
                ylabel = 'angstroms'
                errs = dspErrs
        elif units == 'radians' or units == 'degrees':
            tThErrs = tThFloating - num.tile(tThs, (self.numEta,1) )
            if units == 'degrees':
                ylabel = 'tTh - tTh0  [deg]'
                errs = tThErrs * r2d
            else:
                ylabel = 'tTh - tTh0  [deg]'
                errs = tThErrs
        else:
            raise RuntimeError, "units keyword '%s' is not recognized" % units

        if plot:
            if hasattr(plot, '__call__'):
                'assume is a PlotWrap instance or the like'
                pwTThE = plot
            else:
                pwTThE = plotwrap.PlotWrap(title='line position errors', xlabel='eta index', ylabel=ylabel)
                retval = [errs, pwTThE]

            for iTTh in range(errs.shape[1]):
                pwTThE(num.arange(self.numEta), errs[:,iTTh], noShow=True)
            pwTThE.show()
        else:
            retval = errs    
        
        if outputFile is not None:
            lp = self.planeData.latVecOps['dparms'].reshape(6, 1)
            lp = num.dot( num.diag( num.hstack([num.ones(3,), 180*num.ones(3,)/num.pi]) ), lp )
            laueGrp = self.planeData.getLaueGroup()
            numTTh = len(tThs)
            rhoMax = 204.8
            hklIDs = self.planeData.getHKLID(self.planeData.getHKLs().T)
            hklStr = self.planeData.getHKLs(asStr=True)
            
            rho = num.sqrt((xFl - self.detectorGeom.xc)**2 +
                           (yFl - self.detectorGeom.yc)**2)
            x0, y0 = self.detectorGeom.angToXYO(num.tile(tThs, (self.numEta, 1)),
                                                etaFloating,
                                                units=self.detectorGeom.pixelPitchUnit)
            rho0 = num.sqrt((x0 - self.detectorGeom.xc)**2 +
                            (y0 - self.detectorGeom.yc)**2)
            
            # make output for nutonian
            idsOut  = num.tile(hklIDs, (2*self.numEta, 1)).T.flatten()
            etaOut  = num.vstack([etaFloating, 2*num.pi + etaFloating]).T.flatten()
            rho0Out = num.vstack([rho0, rho0]).T.flatten()/rhoMax
            rhoOut  = num.vstack([rho, rho]).T.flatten()/rhoMax
            
            dataList = ["%d,%1.6e,%1.6e,%1.6e"
                        % (idsOut[i], etaOut[i], rho0Out[i], rhoOut[i])
                        for i in range(len(idsOut))]
            
            if isinstance(outputFile, file):
                fid = outputFile
            else:
                fid = open(outputFile, 'w')
            print >> fid, '# Lattice Parameters: (%.4f, %.4f, %.4f, %.1f, %.1f, %.1f)' % (lp[0],lp[1],lp[2],lp[3],lp[4],lp[5])
            print >> fid, '# Laue Group: %s' % (laueGrp)
            print >> fid, '# Wavelength: %f' % (self.planeData.wavelength)
            print >> fid, '# \n# xc, yc = (%.3f, %.3f)' % (self.detectorGeom.xc, self.detectorGeom.yc)
            print >> fid, '# D = %.3f' % (self.detectorGeom.workDist)
            print >> fid, '# xTilt, yTilt, zTilt = (%1.3e, %1.3e, %1.3e)\n# ' % (self.detectorGeom.xTilt, self.detectorGeom.yTilt, self.detectorGeom.zTilt)
            print >> fid, '\n'.join(['# ' + str(hklIDs[i]) + ' --> ' + hklStr[i] for i in range(len(hklIDs))])
            print >> fid, '# '
            print >> fid, '# ID, azimuth, ideal, measured\n# '
            print >> fid, '\n'.join(dataList)
            fid.close()
        return retval

class MultiRingEval:
    """
    For working with data as rings, particularly for fitting detector geometry or lattice parameters.
    """
    __debug = False
    __print = True

    '''
    check_func in minpack.py calls atleast_1d on the output from Dfun, which hoses the shape output
    (and probably plenty of other things) with sparse matrices
    '''
    __dSparse = False

    def __init__(self, detectorGeom, planeData,
                 indicesList = None, iHKLLists = None,
                 dataFrame = None,
                 funcType = funcTypeDflt,
                 refineParamsDG = True,
                 refineParamsL = False,
                 funcXVecList = None,
                 copyFrame = False,
                 quadr = 3):
      """
      Mostly meant for use with DetectorGeomGE.fit

      If funcXVecList is passed, then entries in this list are used for peak
      function forms, and these peak function forms do not appear in the
      degrees of freedom

      Note that ranges for 2-thetas from planeData need to be such that
      rings are adequately covered

      Can optionally pass indicesList and iHKLLists if they are already handy

      if copyFrame is True, then data in dataFrame is copied
      """

      refineParamsXVec = funcXVecList is None

      if dataFrame is not None:
          self.setupForCall = True
          if copyFrame:
              self.dataFrame = copy.deepcopy(dataFrame)
          else:
              self.dataFrame = dataFrame
      else:
          self.setupForCall = False
          self.dataFrame = None

      self.quadr = quadr
      xi, w = q2db.qLoc(quadr)
      'xi are in [0,1], so need to be centered'
      xi = xi - 0.5
      self.xi = xi
      self.w  = w
      nQP = len(w)

      assert (indicesList is None) == (iHKLLists is None), 'inconsistent arguments for indicesList and iHKLLists'
      if indicesList is None:
          '''groupings of rings could change if planeData.lparms is changing, but
          need to keep fixed number of degrees of freedom so just use initial values'''
          indicesList, iHKLLists = detectorGeom.makeIndicesTThRanges(planeData, cullDupl=True)
      tThs = planeData.getTTh()
      tThRanges = planeData.getTThRanges()

      guessFXVL = False
      if not refineParamsXVec:
          if hasattr(funcXVecList,'__len__'):
              assert len(funcXVecList) == len(iHKLLists),\
                  'funcXVecList needs to have length '+str(len(iHKLLists))
          else:
              guessFXVL = True

      nParamCur = 0
      if refineParamsDG:
          nParamDG = detectorGeom.getNParams()
      else:
          nParamDG = 0
      self.xIndDG = range(nParamCur, nParamCur+nParamDG)
      nParamCur += nParamDG

      if refineParamsL:
          nParamL = len(planeData.lparms)
      else:
          nParamL = 0
      self.xIndL = range(nParamCur, nParamCur+nParamL)
      nParamCur += nParamL

      funcDataList = []
      nPoints = 0
      for iRingSet in range(len(iHKLLists)):
        iHKLList = iHKLLists[iRingSet]
        frameIndices  = indicesList[iRingSet]

        xVecDflt = None
        if not refineParamsXVec and not guessFXVL:
            xVecDflt = funcXVecList[iRingSet]
        if funcType is 'pv':
            func = PeakPV1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
        elif funcType is 'gauss' or funcType is 'normal':
            func = PeakGauss1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
        elif funcType is 'gaussFlatBG':
            func = PeakGauss1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt, linBG=False)
        elif funcType is 'lorentz':
            func = PeakLorentzian1DAtLoc(tThs[(iHKLList)], xVecDflt=xVecDflt)
        else:
            raise RuntimeError,  'unknown funcType : '+str(funcType)
        if not refineParamsXVec and guessFXVL:
            if self.dataFrame is None:
                'used to have something here, but not sure still need to support this'
                raise RuntimeError, 'no way to init 1D function'
            else:
                tThIJ, etaIJ = detectorGeom.xyoToAng(frameIndices[0], frameIndices[1])
                'histogram with floats so that do not overflow int16 data types'
                hist, bin_edges = num.histogram(tThIJ,
                                                weights=num.array(self.dataFrame[frameIndices], dtype=float),
                                                bins=10)
                binCenters = (bin_edges[1:] + bin_edges[:1])*0.5
                func.setXVecDflt(func.guessXVec(binCenters, hist))

        if refineParamsXVec:
            nParam   = func.getNParams()
            xIndices = range(nParamCur, nParamCur+nParam)
        else:
            nParam   = 0
            xIndices = [] # slice(-1,-1,1)
        if self.dataFrame is not None:
            #data = copy.deepcopy(dataFrame[frameIndices])
            normalization = 1.0/num.mean(self.dataFrame[frameIndices])
        else:
            #data = None
            normalization = 1.0

        nPx = len(frameIndices[0])
        iQP = num.tile(frameIndices[0], (nQP,1)).T + num.tile(xi[:,0], (nPx,1))
        jQP = num.tile(frameIndices[1], (nQP,1)).T + num.tile(xi[:,1], (nPx,1))
        'do not worry about dvQP kinds of contributions for now'
        wQP = num.tile(w, (nPx,1))
        #funcDataList.append((func, xIndices, iQP, jQP, wQP, data, indices, normalization))
        funcDataList.append((func, xIndices, iQP, jQP, wQP, frameIndices, normalization))
        nPoints += nPx
        nParamCur += nParam

      'store stuff'
      self.detectorGeom = detectorGeom
      self.planeData    = planeData
      self.nParamDG     = nParamDG
      self.nParamL      = nParamL
      self.funcDataList = funcDataList
      self.nParam       = nParamCur
      self.nPoints      = nPoints
      self.iHKLLists    = iHKLLists

      if self.nParam == 0:
          print >> sys.stderr, 'WARNING: have 0 free parameters in '+str(self.__class__.__name__)

      self.useDFun = True
      if refineParamsDG and not refineParamsL and not refineParamsXVec:
          'only doing dg params, derivatives are by finite differencing, so do not bother with Dfun'
          self.useDFun = False

      xVec = num.zeros(self.nParam)
      if self.nParam > 0:
          'go ahead and set xVecGuess too'
          if self.nParamDG > 0:
              xVec[self.xIndDG] = self.detectorGeom.getParams()
          if self.nParamL > 0:
              xVec[self.xIndL]  = self.planeData.lparms
          for iRingSet in range(len(iHKLLists)):
            iHKLList = iHKLLists[iRingSet]
            func, xIndices, iQP, jQP, wQP, frameIndices, normalization = self.funcDataList[iRingSet]
            if len(xIndices) > 0:
                if self.dataFrame is None:
                    # 'set up something in case want to make synthetic data'
                    # widths = (tThRanges[iHKLList,1] - tThRanges[iHKLList,0])/8.
                    # mins   = num.tile(1.0,   len(iHKLList))
                    # maxs   = num.tile(200.0, len(iHKLList))
                    pass
                else:
                    tThIJ, etaIJ = detectorGeom.xyoToAng(frameIndices[0], frameIndices[1])
                    'histogram with floats so that do not overflow int16 data types'
                    hist, bin_edges = num.histogram(tThIJ,
                                                    weights=num.array(self.dataFrame[frameIndices], dtype=float),
                                                    bins=10)
                    binCenters = (bin_edges[1:] + bin_edges[:1])*0.5
                    xVec[(xIndices)] = func.guessXVec(binCenters, hist)
      self.xVecGuess = xVec

      return
    def getFuncXVecList(self, xVec):
        funcXVecList = []
        for func, xIndices, iQP, jQP, wQP, frameIndices, normalization in self.funcDataList:
            if len(xIndices) > 0:
                funcXVecList.append(xVec[(xIndices)])
            else:
                funcXVecList.append(func.xVecDflt)
        return funcXVecList
    def setFuncXVecList(self, funcXVecList):
        """
        only okay if funcXVecList set on init
        """
        for iFunc, funcData in enumerate(self.funcDataList):
            func, xIndices, iQP, jQP, wQP, frameIndices, normalization = funcData
            assert len(xIndices) == 0, \
                'setting funcXVecList for instance with DOF for funcs'
            func.setXVecDflt(funcXVecList[iFunc])
        return
    def guessXVec(self):
      return self.xVecGuess
    def getNParam(self):
      return self.nParam
    def __dgfd(self, xVec, dGScalings, iQP, jQP):
        # J is 3D, b/c tThQP_ref is 2D
        assert self.nParamDG > 0,\
            'do not have detector geometry parameter degrees of freedom'
        xVecDG = xVec[self.xIndDG]
        self.detectorGeom.updateParams(*xVecDG)
        tThQP_ref, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)
        assert len(tThQP_ref.shape) == 2, \
            'have not gotten around to generalizing'
        pert = 1.0e-5
        J = num.empty(num.hstack((tThQP_ref.shape,len(xVecDG))))
        for iX, xVal in enumerate(xVecDG):
            xVecDG_p = copy.deepcopy(xVecDG)
            thisPert = pert * dGScalings[iX]
            xVecDG_p[iX] = xVal + thisPert
            self.detectorGeom.updateParams(*xVecDG_p)
            tThQP_p, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)
            J[:,:,iX] = (tThQP_p - tThQP_ref) / thisPert
        return tThQP_ref, J
    def deval(self, xVec):
      '''
      useful to pass, for example, as Dfun to leastsq;
      bit of a misnomer in that deval is the derivative of __call__, not eval
      '''

      retvalShape = (self.nPoints,len(xVec))
      if self.__dSparse:
          # retval = sparse.dok_matrix(retvalShape)
          # retval = sparse.lil_matrix(retvalShape)
          retval = sparse.coo_matrix(retvalShape)
      else:
          retval = num.zeros(retvalShape)

      if self.nParamDG > 0:
          if self.__debug: print 'xIndDG : '+str(self.xIndDG)
          self.detectorGeom.updateParams(*xVec[self.xIndDG])
          dGScalings = self.detectorGeom.getParamScalings() # for finite differencing
      else:
          dGScalings = None
      if self.nParamL > 0:
          if self.__debug: print 'xIndL : '+str(self.xIndL)
          self.planeData.lparms = xVec[self.xIndL] # property, calls set_lparms
          d_tThs_d_lparms = self.planeData.getDD_tThs_lparms()
          tThs = self.planeData.getTTh()
          d_centers_d_lparms = []
          for iRingSet in range(len(self.iHKLLists)):
              iHKLList = self.iHKLLists[iRingSet]
              func, xIndices, iQP, jQP, wQP, frameIndices, normalization = self.funcDataList[iRingSet]
              func.setCenters(tThs[(iHKLList)])
              d_centers_d_lparms.append(d_tThs_d_lparms[(iHKLList)])

      nPoints = 0
      for iRingSet in range(len(self.iHKLLists)):
        func, xIndices, iQP, jQP, wQP, frameIndices, normalization = self.funcDataList[iRingSet]
        nPointsThis = wQP.shape[0] # len(data)
        if self.__debug: print 'xIndices : '+str(xIndices)

        if self.nParamDG > 0:
            tThQP, d_tThQP_d_xVecDG = self.__dgfd(xVec, dGScalings, iQP, jQP)
        else:
            tThQP, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)

        if len(xIndices) > 0:
            #evalQP       = func(xVec[(xIndices)], tThQP)
            d_evalQP_d_x = func.d_dx(xVec[(xIndices)], tThQP)
            d_evalQP_d_p = func.d_dp(xVec[(xIndices)], tThQP)

            if self.nParamL > 0:
                d_evalQP_d_centers = func.d_dCenters(xVec[(xIndices)], tThQP)

            #eval          = normalization * num.sum(evalQP * wQP, axis=1)
            for iXInd, xInd in enumerate(xIndices):
                temp = normalization * num.sum(d_evalQP_d_x[:,:,iXInd] * wQP, axis=1)
                if self.__dSparse:
                    i = num.mgrid[nPoints:nPoints+nPointsThis,xInd:xInd+1]
                    retval = retval + sparse.coo_matrix(
                        (temp.flatten(), (i[0].flatten(), i[1].flatten())),
                        shape=retvalShape)
                else:
                    retval[nPoints:nPoints+nPointsThis,xInd] = temp
        else:
            'not len(xIndices) > 0'
            if self.nParamDG > 0:
                d_evalQP_d_p       = func.d_dp(func.xVecDflt, tThQP)
            if self.nParamL > 0:
                d_evalQP_d_centers = func.d_dCenters(func.xVecDflt, tThQP)

        for iXInd, xInd in enumerate(self.xIndDG):
            'should only end up in here if self.nParamDG > 0'
            temp = normalization * num.sum((d_evalQP_d_p * d_tThQP_d_xVecDG[:,:,iXInd]) * wQP, axis=1)
            if self.__dSparse:
                i = num.mgrid[nPoints:nPoints+nPointsThis,xInd:xInd+1]
                retval = retval + sparse.coo_matrix(
                    (temp.flatten(), (i[0].flatten(), i[1].flatten())),
                    shape=retvalShape)
            else:
                retval[nPoints:nPoints+nPointsThis,xInd] = temp
        for iXInd, xInd in enumerate(self.xIndL):
            'should only end up in here if self.nParamL > 0'
            temp = normalization * num.sum((num.dot(d_evalQP_d_centers, d_centers_d_lparms[iRingSet][:,iXInd])) * wQP, axis=1)
            if self.__dSparse:
                i = num.mgrid[nPoints:nPoints+nPointsThis,xInd:xInd+1]
                retval = retval + sparse.coo_matrix(
                    (temp.flatten(), (i[0].flatten(), i[1].flatten())),
                    shape=retvalShape)
            else:
                retval[nPoints:nPoints+nPointsThis,xInd] = temp

        nPoints += nPointsThis

      if self.__dSparse:
          retval = retval.tocsr() # or tocsr would be better?
      return retval

    def eval(self, xVec, thisframe=None):
      'if thisframe is passed, the put values on the frame'

      retval = num.empty(self.nPoints)

      'first update the detector geometry'
      if self.nParamDG > 0:
          self.detectorGeom.updateParams(*xVec[self.xIndDG])

      if self.nParamL > 0:
          self.planeData.lparms = xVec[self.xIndL] # property, calls set_lparms
          tThs = self.planeData.getTTh()
          for iRingSet in range(len(self.iHKLLists)):
              iHKLList = self.iHKLLists[iRingSet]
              func, xIndices, iQP, jQP, wQP, frameIndices, normalization = self.funcDataList[iRingSet]
              func.setCenters(tThs[(iHKLList)])

      nPoints = 0
      for func, xIndices, iQP, jQP, wQP, frameIndices, normalization in self.funcDataList:
        nPointsThis = wQP.shape[0] # len(data)

        tThQP, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)
        evalQP       = func(xVec[(xIndices)], tThQP) # if xVec ends up being [], xVecDflt gets used
        eval         = num.sum(evalQP * wQP, axis=1)
        retval[nPoints:nPoints+nPointsThis] = eval
        if thisframe is not None:
          thisframe[frameIndices] = eval

        nPoints += nPointsThis

      return retval
    def __call__(self, xVec, makePlots=False, plotTitlePrefix=''):
        'meant for call in leastsq type algorithm'

        assert self.setupForCall, 'not setup for call'

        retval = num.empty(self.nPoints)

        evalAllRings = self.eval(xVec)

        if makePlots:
            frameData = self.detectorGeom.frame()
            frameEval = self.detectorGeom.frame()
            frameDiff = self.detectorGeom.frame(dtype='float')

            plotWinRadial = plotwrap.PlotWin(2,1, relfigsize=(8,3),
                                             title=plotTitlePrefix+'ring evaluation results')
            'using twinx() might be better, but not plumbed for that'
            pwREval = plotwrap.PlotWrap(window=plotWinRadial, showByDefault=False)
            pwRDiff = plotwrap.PlotWrap(window=plotWinRadial, showByDefault=False, axprops={'sharex':pwREval.a})

        nPoints = 0
        for iFunc, funcData in enumerate(self.funcDataList):
          func, xIndices, iQP, jQP, wQP, frameIndices, normalization = funcData
          nPointsThis = wQP.shape[0] # len(data)

          data = self.dataFrame[frameIndices]

          eval  = evalAllRings[nPoints:nPoints+nPointsThis]
          diff  = eval - data

          retval[nPoints:nPoints+nPointsThis] = normalization * diff

          if makePlots:

            frameData[frameIndices] = data
            frameEval[frameIndices] = eval
            frameDiff[frameIndices] = diff

            histTThCen, histIntensity, width = self.__radialBinData(iFunc, data)
            evalIntensity = func(xVec[(xIndices)], histTThCen)
            pwREval.a.bar(histTThCen, histIntensity, width=width, align='center', alpha=0.5)
            '''have not set things up to match line styles with drawRings, so make all lines black;
            might be able to fix this if switch drawRings to do something like cullDupl too
            '''
            pwREval(histTThCen, evalIntensity, style='k-')
            pwRDiff(histTThCen, (evalIntensity - histIntensity) / histIntensity, style='k-')

          nPoints += nPointsThis
        if makePlots:

          plotWinFrames = plotwrap.PlotWin(1,3,title='ring evaluation results')
          #
          pwData  = plotwrap.PlotWrap(window=plotWinFrames,
                                      title='data; max() = %g' % (num.max(frameData))
                                      )
          axprops = {'sharex':pwData.a, 'sharey':pwData.a}
          pwEval  = plotwrap.PlotWrap(window=plotWinFrames,
                                      title='eval',
                                      axprops=axprops)
          pwDiff  = plotwrap.PlotWrap(window=plotWinFrames,
                                      title='diff; max(abs()) = %g' % (num.max(num.abs(frameDiff))),
                                      axprops=axprops)
          #
          ReadGE.display(frameData, pw=pwData)
          ReadGE.display(frameEval, pw=pwEval)
          ReadGE.display(frameDiff, pw=pwDiff, cmap=None)
          self.detectorGeom.drawRings(pwData, self.planeData, legendLoc=None)
          self.detectorGeom.drawRings(pwEval, self.planeData, legendLoc=None)
          self.detectorGeom.drawRings(pwDiff, self.planeData, legendLoc=None)

          plotWinFrames.show()
          plotWinRadial.show()
        if self.__print:
          print ' norm : '+str(num.linalg.norm(retval))
        return retval
    def doFit(self, xVec0=None, **lsKWArgs):
        """
        lsKWArgs can have things like ftol and xtol for leastsq
        """
        if xVec0 is None:
            xVec0 = self.guessXVec()
        tic = time.time()
        if self.useDFun:
            x, ier = leastsq(self, xVec0, Dfun=self.deval, **lsKWArgs)
        else:
            x, ier = leastsq(self, xVec0, **lsKWArgs)
        toc = time.time()
        dt = toc - tic
        print 'leastsq took %g seconds to run' % (dt)
        if not [1,2,3,4].count(ier):
          print >> sys.stderr, 'error code %d from leastsq' % (ier)
          raise RuntimeError, 'error from leastsq'
        return num.atleast_1d(x)
    def radialPlotData(self, dataFrame=None, plotTitlePrefix=''):
        'for simple radial plotting, useful if other things are mysteriously breaking'
        if dataFrame is None:
            dataFrame = self.dataFrame
        assert dataFrame is not None, 'need data to work with'

        pwREval = plotwrap.PlotWrap(showByDefault=False, title=plotTitlePrefix+'ring evaluation results')
        for iFunc, funcData in enumerate(self.funcDataList):
            func, xIndices, iQP, jQP, wQP, frameIndices, normalization = funcData
            histTThCen, histIntensity, width = self.__radialBinData(iFunc, dataFrame[frameIndices])
            pwREval.a.bar(histTThCen, histIntensity, width=width, align='center', alpha=0.5)
        pwREval.show()
        return pwREval
    def radialFitXVec(self, dataFrame=None, plot=False, plotTitlePrefix='', quadr1d=None):
        """
        if dataFrame is not provided, use self.dataFrame

        if quadr1d is not specified, use quadr specified in init
        """
        if dataFrame is None:
            dataFrame = self.dataFrame
        assert dataFrame is not None, 'need data to work with'

        quadr = self.quadr
        if quadr1d is not None:
            quadr = quadr1d
        funcXVecList = []
        cFitList = []
        pw = None
        if plot > 1:
            plotList = []
            plotRetval = plotList
        elif plot:
            pw = plotWinRadial = plotwrap.PlotWin(2,1, relfigsize=(8,3),
                                             title=plotTitlePrefix+'ring evaluation results')
            #pw = plotwrap.PlotWrap(title=plotTitlePrefix+'radial fit results', showByDefault=False,
            #                       figsize=(8,2))
            pwREval = plotwrap.PlotWrap(window=plotWinRadial, showByDefault=False)
            pwRDiff = plotwrap.PlotWrap(window=plotWinRadial, showByDefault=False, axprops={'sharex':pwREval.a})
            plotRetval = plotWinRadial
        hkls = self.planeData.getHKLs(asStr=True)
        for iFunc, funcData in enumerate(self.funcDataList):
            func, xIndices, iQP, jQP, wQP, frameIndices, normalization = funcData
            data = self.dataFrame[frameIndices]
            assert self.quadr, \
                'self.quadr does not evaluate to True'
            histTThCen, histIntensity, width, xFit, cFit, evalIntensity = \
                self.__radialBinData(iFunc, data, fitFuncXVec=quadr)
            funcXVecList.append(xFit)
            cFitList.append(cFit)
            if plot:
                iHKLList = self.iHKLLists[iFunc]
                if plot > 1:
                    pw = plotwrap.PlotWrap(title=plotTitlePrefix+str(hkls[(iHKLList)]))
                    pw.a.bar(histTThCen, histIntensity, width=width, align='center', alpha=0.5)
                    pw(histTThCen, evalIntensity, style='k-')
                    plotList.append(pw)
                    pw = None
                else:
                    pwREval.a.bar(histTThCen, histIntensity, width=width, align='center', alpha=0.5)
                    pwREval(histTThCen, evalIntensity, style='k-')
                    pwRDiff(histTThCen, evalIntensity-histIntensity, style='k-')
        if plot:
            if pw is not None:
                pw.show()
            retval = funcXVecList, cFitList, plotRetval
        else:
            retval = funcXVecList, cFitList
        return retval
    def __radialBinData(self, iFunc, data,
                        nBinsPerPeak = 30,
                        fitFuncXVec=False):
        """
        useful for checking on ring positions, and for fitting funcXVec

        data needs to be collected for the iFunc in question, as by doing
        self.dataFrame[frameIndices]
        """

        func, xIndices, iQP, jQP, wQP, frameIndices, normalization = self.funcDataList[iFunc]
        assert len(data) == len(frameIndices[0]),\
            'data is wrong length for iFunc'

        tThQP, etaQP = self.detectorGeom.xyoToAng(iQP, jQP)

        """
        NOTE: if wQP is made to have dv contributions, then may
        need to change the following;
        though the following is not quite consistent with the least
        squares system anyway, it is just to have something more
        to look at
        """
        #
        weights = wQP * num.tile(data, (wQP.shape[1],1)).T
        #
        'Do not do the following because need to renormalize'
        # histIntensity, histTThBins, histPatches = \
        #     pwREval.a.hist(tThQP.flatten(), weights=weights.flatten(),
        #                    bins=60,
        #                    histtype='bar', alpha=0.5)
        range = (tThQP.min(), tThQP.max())
        nBins = nBinsPerPeak*func.getNPeaks()
        width = (range[1]-range[0])/nBins
        bins  = num.arange(range[0], range[1]+width*0.5, width)
        histIW, bin_edges = num.histogram(tThQP.flatten(), range=range, bins=bins, weights=weights.flatten())
        histW,  bin_edges = num.histogram(tThQP.flatten(), range=range, bins=bins, weights=wQP.flatten())
        histTThCen   = 0.5 * (bins[:-1] + bins[1:]) # bin centers
        histIntensity = histIW / histW

        if fitFuncXVec:

            quadr = 1
            if isinstance(fitFuncXVec,int):
                quadr = int(fitFuncXVec)

            if len(xIndices) > 0:
                xVecGuess = self.xVecGuess[(xIndices)]
            else:
                xVecGuess = func.xVecDflt

            tThs = self.planeData.getTTh()
            iHKLList = self.iHKLLists[iFunc]

            xFit, cFit, evalIntensity = func.fitFloatingCenter(
                histTThCen, histIntensity, xVecGuess, tThs[(iHKLList)], weights=quadr, tThWidth=width)

        retval = [histTThCen, histIntensity, width]
        if fitFuncXVec:
            retval += xFit, cFit, evalIntensity

        return retval

class DetectorBase(object):
    """
    base class for a detector
    """

    # basis vectors
    Xl = num.vstack([1, 0, 0])                      # X in the lab frame
    Yl = num.vstack([0, 1, 0])                      # Y in the lab frame
    Zl = num.vstack([0, 0, 1])                      # Z in the lab frame

    def __init__(self, reader):
        self.__reader = reader
        self.refineFlags = self.getRefineFlagsDflt()
        return
    
    def get_reader(self):
        return self.__reader
    reader = property(get_reader, None, None)

    def frame(self, *args, **kwargs):
        retval = self.reader.frame(*args, **kwargs)
        return retval

    @classmethod
    def getRefineFlagsDflt(cls):
        raise NotImplementedError

    def setupRefinement(self, flags):
        assert len(flags) == len(self.refineFlags), 'flags wrong length, do you know what you are doing?'
        self.refineFlags = num.array(flags, dtype='bool')
        return
    def getNParams(self):
        return num.sum(self.refineFlags)
    def updateParams(self, *args):
        raise NotImplementedError
        return
    def getParams(self, allParams=False):
        retval = None
        raise NotImplementedError
        return retval
    def getParamScalings(self):
        'scalings, suitable for scaling perturbations for finite differencing'
        retval = None
        raise NotImplementedError
        return retval
    def getPVecScaling(self):
        'scaling, suitable for scaling perturbations for finite differencing'
        retval = None
        raise NotImplementedError
        return retval

class Detector2DRC(DetectorBase):
    """
    base class for 2D row-column detectors
    """

    __pixelPitchUnit = 'mm'

    tilt = num.zeros(3)   # must initialize tilt to ndarray

    chiTilt = None
    
    def __init__(self, 
                 ncols, nrows, pixelPitch, 
                 vFactorUnc, vDark, 
                 reader, *args, **kwargs):
        
        if reader is None:
            readerKWArgs = kwargs.pop('readerKWArgs', {})
            reader = newGenericReader(ncols, nrows, **readerKWArgs)
        
        """
        The following is meant to facilitate creation of generic detector types
        """
        optionalFuncs = ('getDParamDflt', 
                         'setDParamZero', 
                         'getDParamScalings', 
                         'getDParamRefineDflt', 
                         'radialDistortion') 
        for optFunc in optionalFuncs:
            if kwargs.has_key(optFunc):
                self.__setattr__(optFunc, kwargs.pop(optFunc))
        
        DetectorBase.__init__(self, reader)
        
        self.__ncols = ncols
        self.__nrows = nrows
        self.__vFactorUnc = vFactorUnc
        self.__vDark = vDark
        self.__pixelPitch = pixelPitch

        if len(args) == 0:
            self.__initWithDefault(**kwargs)
        elif hasattr(args[0], 'xyoToAng'): # this should work ok (JVB)
            self.__initFromDG(*args, **kwargs)
        elif len(args) == 1:
            if hasattr(args[0], '__len__'):
                assert len(args[0]) is 6, "your initial parameters list is the wrong length"
                xc, yc, workDist, xTilt, yTilt, zTilt = args[0]
                self.__initFromParams(xc, yc, workDist, xTilt, yTilt, zTilt, **kwargs)
            else:
                raise RuntimeError, "don't know what to do with your initialization argument "+str(args)
        else:
            self.__initFromParams(*args, **kwargs)
        self.pVecUncertainties = None

        return

    def set_ncols(self, ncols):
        raise RuntimeError, 'set of ncols not allowed'
    def get_ncols(self):
        return self.__ncols
    ncols = property(get_ncols, set_ncols, None)

    def set_pixelPitch(self, pixelPitch):
        raise RuntimeError, 'set of pixelPitch not allowed'
    def get_pixelPitch(self):
        return self.__pixelPitch
    pixelPitch = property(get_pixelPitch, set_pixelPitch, None)

    def set_pixelPitchUnit(self, pixelPitchUnit):
        raise RuntimeError, 'set of pixelPitchUnit not allowed'
    def get_pixelPitchUnit(self):
        return self.__pixelPitchUnit
    pixelPitchUnit = property(get_pixelPitchUnit, set_pixelPitchUnit, None)

    def set_nrows(self, nrows):
        raise RuntimeError, 'set of nrows not allowed'
    def get_nrows(self):
        return self.__nrows
    nrows = property(get_nrows, set_nrows, None)

    def getSize(self):
        retval = self.ncols * self.nrows
        return retval
    def __len__(self):
        return self.getSize()

    def getShape(self):
        return (self.ncols, self.nrows)
    def setShape(self):
        raise RuntimeError, 'what are you thinking'
    shape = property(getShape, setShape, None)

    def getExtent(self):
        return num.array([self.ncols, self.nrows]) * self.pixelPitch
    def setExtent(self):
        raise RuntimeError, 'what are you thinking'
    extent = property(getExtent, setExtent, None)

    # methods that a specific detector geometry needs to implement:
    def getDParamDflt(self):
        raise NotImplementedError
    def setDParamZero(self):
        raise NotImplementedError
    def getDParamScalings(self):
        raise NotImplementedError
    def getDParamRefineDflt(self):
        raise NotImplementedError
    def radialDistortion(self, xin, yin, invert=False):
        raise NotImplementedError

    def getNewReader(self, filename, *args, **kwargs):
        newR = self.reader.__class__(filename, self.ncols, self.nrows, *args, **kwargs)
        return newR

    # methods that may want to override for a specific detector geometry
    def getParamDflt(self):
        xc        = 0.5*self.ncols*self.pixelPitch
        yc        = 0.5*self.nrows*self.pixelPitch
        workDist  = 1000.0
        xTilt     = 0.0
        yTilt     = 0.0
        zTilt     = 0.0
        dparms    = self.getDParamDflt()
        return (xc, yc, workDist, xTilt, yTilt, zTilt, dparms)

    def __initFromParams(self, xc, yc, workDist, xTilt, yTilt, zTilt, distortionParams=None):
        """
        inputs:
        xc, yc : beam center, relative to pixel centers, so that
                 (xc, yc) = (0,0) puts the beam center in the
                 _center_ of the corner pixel
        workDist : working distance, as valWUnit instance, eg:
        workDist = valunits.valWUnit('workDist', 'length', 1.9365, 'meter')
        distortionParams : distortion parameters

        eventually, may want to have this be able to determine
        geometric parameters from other input data
        """

        if hasattr(xc, 'getVal'):
            self.xc        = xc.getVal(self.pixelPitchUnit)
            self.yc        = yc.getVal(self.pixelPitchUnit)
            self.workDist  = workDist.getVal(self.pixelPitchUnit)
            xTilt     = xTilt.getVal('radians')
            yTilt     = yTilt.getVal('radians')
            zTilt     = zTilt.getVal('radians')
        else:
            self.xc        = xc
            self.yc        = yc
            self.workDist  = workDist
        self.setTilt((xTilt, yTilt, zTilt))
        dParamDflt = self.getDParamDflt()
        if distortionParams is None:
            distortionParams = dParamDflt
        self.dparms    = copy.deepcopy(distortionParams)
        assert len(self.dparms) == len(dParamDflt), 'wrong length for distortion parameters'

        self.__initBase()

        return

    def getRefineFlagsDflt(cls):
        # for refinement
        # --------------->       (   xc,    yc,     D,    xt,    yt,    zt)
        retval = num.atleast_1d(
            tuple( [True,  True,  True,  True,  True, False] +
                   list(cls.getDParamRefineDflt())
                   ) )
        return retval
    def __makePList(self, flatten=True):
        pList = [self.xc, self.yc, self.workDist, self.xTilt, self.yTilt, self.zTilt, self.dparms]
        if flatten:
            pList = num.hstack(pList)
        return pList
    def __updateFromPList(self, plist):
        self.xc, self.yc, self.workDist, self.xTilt, self.yTilt, self.zTilt = plist[0:6]
        if hasattr(plist[6],'__len__'):
            assert len(plist) >= 6, 'plist of wrong length : '+str(plist)
            self.dparms = plist[6]
        else:
            self.dparms = plist[6:]
        assert len(self.dparms) == len(self.getDParamDflt()), \
            'dparms wrong length : '+str(self.dparms)
        return plist
    def updateParams(self, *args): # xc, yc, workDist, xTilt, yTilt, zTilt, *dparms
        if len(args) == 1:
            args = num.atleast_1d(args[0])
        assert len(args) == self.getNParams(), 'wrong length for distortion parameters'
        plist = self.__makePList()
        plist[num.array(self.refineFlags)] = args
        self.__updateFromPList(plist)
        return
    def getParams(self, allParams=False):
        plist = num.array(self.__makePList())
        if allParams:
            retval = plist
        else:
            idx = num.array(self.refineFlags)
            retval = plist[idx]
        return retval
    def getParamScalings(self):
        'scalings, suitable for scaling perturbations for finite differencing'
        plistScalings = num.array([
                self.pixelPitch*self.__nrows, self.pixelPitch*self.__ncols, # self.xc, self.yc
                1000.0,                 # self.workDist
                1.0, 1.0, 1.0,          # self.xTilt, self.yTilt, self.zTilt
                ] + list(self.getDParamScalings()) )
        return plistScalings[num.array(self.refineFlags)]
    def getPVecScaling(self):
        return 0.01
    def __initBase(self):
        
        ' no precession or chiTilt by default init'
        self.pVec    = None
        self.chiTilt = None
        
        self.refineFlags = self.getRefineFlagsDflt()

        'stuff for drawing and fitting rings'
        self.withRanges = True
        self.asMasked = False
        self.withTThWidth = True
        self.nEta = 180
        self.xFitRings = None
        # self.fitRingsFunc = None # potential for circular references if keep this
        self.lineArgs = {
            'antialiased' : True, # False is faster, but uglier
            'linewidth' : 0.5
            }
        return

    def __initWithDefault(self):

        self.__updateFromPList(self.getParamDflt())
        self.__initBase()
        
        return
    def __initFromDG(self, detectorGeomOther, pVec=None, chiTilt=None):
        """...not sure we'll keep chiTilt here as is"""
        self.__initFromParams( *detectorGeomOther.__makePList(flatten=False) )
        self.refineFlags = copy.deepcopy(detectorGeomOther.refineFlags)
        self.pVec = copy.deepcopy(pVec)
        self.chiTilt = copy.deepcopy(chiTilt)
        return

    def setTilt(self, tilt):
        self.tilt = num.array(tilt)

        # tilt angles
        gX = self.tilt[0]
        gY = self.tilt[1]
        gZ = self.tilt[2]

        # rotation 1: gX about Xl
        ROTX = rotMatOfExpMap(gX * self.Xl)

        # the transformed Yl axis (Yd)
        Yd = num.dot(ROTX, self.Yl)

        # rotation 2: gY about Yd
        ROTY = rotMatOfExpMap(gY * Yd)

        # the transformed Zl axis (Zd)
        Zd = num.dot(num.dot(ROTY, ROTX), self.Zl)

        # rotation 3: gZ about Zd
        ROTZ = rotMatOfExpMap(gZ * Zd)

        # change of basis matrix from hatPrime to Prime
        self.ROT_l2d = num.dot(ROTZ, num.dot(ROTY, ROTX))

        # tilted X-Y plane normal
        self.N = Zd

        return
    def getXTilt(self):
        return self.tilt[0]
    def setXTilt(self, xTilt):
        self.tilt[0] = xTilt
        self.setTilt(self.tilt)
        return
    xTilt = property(getXTilt, setXTilt, None)
    #
    def getYTilt(self):
        return self.tilt[1]
    def setYTilt(self, yTilt):
        self.tilt[1] = yTilt
        self.setTilt(self.tilt)
        return
    yTilt = property(getYTilt, setYTilt, None)
    #
    def getZTilt(self):
        return self.tilt[2]
    def setZTilt(self, zTilt):
        self.tilt[2] = zTilt
        self.setTilt(self.tilt)
        return
    zTilt = property(getZTilt, setZTilt, None)

    def getVDark(self):
        return self.__vDark
    
    def getVScale(self, vThese):
        """
        get scale factors for use in uncertainty quantification
        """
        vScale      = self.__vFactorUnc * (vThese)
        vScaleFloor = 0.1 * (vScale[vScale > 0].mean())
        assert vScaleFloor > 0, \
            'vScaleFloor <= 0'
        vScale[vScale < vScaleFloor] = vScaleFloor
        vScale = 1.0 / vScale
        return vScale

    def getAngPixelSize(self, xyo, delta_omega):
        'get pixel size in angular coordinates at a given cartesian coordinate position'
        xyoCorners = num.tile(xyo, (8,1)).T + \
            num.array([[[[i,j,k] for i in [-0.5,0.5]] for j in [-0.5,0.5]] for k in [-0.5*delta_omega, 0.5*delta_omega]]).reshape(8,3).T
        angCorners = num.array(self.xyoToAngMap(*xyoCorners))
        unc = num.max(angCorners, axis=1) - num.min(angCorners, axis=1)
        return unc

    def cartesianCoordsOfPixelIndices(self, row, col, ROI=None):
        """
        converts [i, j] pixel array indices to cartesian spatial coords
        where the lower left corner of the image is (0, 0)

        Output units are in the pixel pitch units (see self.pixelPitch)

        Will optionally take the upper left-hand corner (min row, min col) of
        a ROI when dealing with subregions on the detector as in when zooming in on
        diffraction spots...

        *) explicitly enforce this to be self-consistent with radial distortion correction, etc...
        """

        # properly offset in case
        if ROI is not None:
            assert len(ROI) is 2, 'wrong length for ROI; should be 2 integers representing the UL corner'
            row = row + ROI[0]
            col = col + ROI[1]

        xout = self.pixelPitch*(col + 0.5)
        yout = self.pixelPitch*(self.__nrows - (row + 0.5))

        return xout, yout

    def pixelIndicesOfCartesianCoords(self, x, y, ROI=None):
        """
        converts [i, j] pixel array indices to cartesian spatial coords
        where the lower left corner of the image is (0, 0)

        Output units are in the pixel pitch units (see self.pixelPitch)

        Will optionally take the upper left-hand corner (min row, min col) of
        a ROI when dealing with subregions on the detector as in when zooming in on
        diffraction spots...

        *) explicitly enforce this to be self-consistent with radial distortion correction, etc...
        """

        row = (self.__nrows - 0.5) - y/self.pixelPitch
        col = x/self.pixelPitch - 0.5

        # properly offset in case
        if ROI is not None:
            assert len(ROI) is 2, 'wrong length for ROI; should be 2 integers representing the UL corner'
            row = row - ROI[0]
            col = col - ROI[1]

        return row, col

    #
    # Real geometric stuff below -- proceed at own risk
    #
    def angToXYO_V(self, tth, eta_l, *args, **kwargs):
        """
        opposite of xyoToAng
        """

        outputDV              = False
        outputForGve          = False
        rhoRange              = ()
        toPixelUnits          = True
        applyRadialDistortion = True

        outputShape = num.shape(tth)
        tth         = num.asarray(tth).flatten()
        eta_l       = num.asarray(eta_l).flatten()

        numPts = len(tth)                # ... no check for y0 or omega
        ome = ()
        if len(args) is not 0:
            ome = num.atleast_1d(args[0])

        nzeros = num.zeros(numPts)

        # argument handling
        kwarglen = len(kwargs)
        if kwarglen > 0:
            argkeys = kwargs.keys()
            for i in range(kwarglen):
                if argkeys[i] is 'outputGve':
                    outputForGve = kwargs[argkeys[i]]
                elif argkeys[i] is 'outputDV':
                    outputDV = kwargs[argkeys[i]]
                elif argkeys[i] is 'units':
                    if kwargs[argkeys[i]] is 'pixels':
                        toPixelUnits = True
                    elif kwargs[argkeys[i]] is self.pixelPitchUnit:
                        toPixelUnits = False
                    else:
                        raise RuntimeError, 'Output units \'%s\'not understood!' % (str(kwargs[argkeys[i]]))
                elif argkeys[i] is 'rhoRange':
                    tthRange = kwargs[argkeys[i]]
                    assert len(tthRange) is 2, 'Radial range should have length 2'
                elif argkeys[i] is 'rdist':
                    if not isinstance(kwargs[argkeys[i]], bool):
                        raise RuntimeError, 'Expecting boolean for rdist kewyord argument; got' \
                              + str(applyRadialDistortion)
                    else:
                        applyRadialDistortion = kwargs[argkeys[i]]
                else:
                    raise RuntimeError, 'Unrecognized keyword argument: ' + str(argkeys[i])

        # make center-based cartesian coord's
        #   - SHOULD BE IN PIXEL PITCH UNITS (MM)
        #   ... maybe add hard check on this in future
        #
        if self.pVec is None:
            XC = nzeros
            YC = nzeros

            D  = num.tile(self.workDist, numPts)
        else:
            assert len(ome) == numPts, 'with precession, omega argument consistent with ' \
                   + 'x and y (or i and j) is required'
            # if here, we have a precession vector and must deal with it
            #
            #   - ome is taken as a CCW (+) rotation of the SAMPLE FRAME about Y
            #   - when the BASIS is transformed by R, vector comp's must transform by R'
            R_s2l = rotMatOfExpMap( num.tile(ome, (3, 1)) * num.tile(self.Yl, (1, numPts)) )

            # array of rotated precession vector components
            if not hasattr(self.pVec, '__len__'):
                raise RuntimeError, 'pVec must be array-like'

            self.pVec = num.asarray(self.pVec)

            grainCM_l = num.dot(R_s2l, self.pVec.reshape(3, 1))
            if grainCM_l.ndim == 3:
                grainCM_l = grainCM_l.squeeze().T

            XC = grainCM_l[0, :]
            YC = grainCM_l[1, :]

            D  = self.workDist + grainCM_l[2, :] # now array of D's

        # make radii
        rho_l = D * num.tan(tth)

        #
        # ------- ASSIGN POINT COORD'S AND FORM ROTATION
        #
        # origins of the scattering (P1) and lab (P2) frames
        #   - the common frame for arithmatic is the scattering frame
        P1 = num.vstack([XC, YC, D])
        P2 = num.zeros((3, numPts))

        # tilt calculations moved into setTilt

        # Convert to cartesian coord's in lab frame
        P3 = num.vstack( [ rho_l * num.cos(eta_l) + XC, rho_l * num.sin(eta_l) + YC, nzeros ] )

        #
        # ------- SOLVE FOR RAY-PLANE INTERSECTION
        #
        u = num.tile( num.dot(self.N.T, (P2 - P1)) / num.dot(self.N.T, P3 - P1), (3, 1) )

        P4_l = P1 + u * (P3 - P1)

        P4_d = num.dot(self.ROT_l2d.T, P4_l)

        if applyRadialDistortion:
            X_d, Y_d = self.radialDistortion(P4_d[0, :], P4_d[1, :], invert=True)
            # note that the Z comps should all be zeros anyhow
            P4_d = num.vstack( [X_d, Y_d, nzeros] )

        if len(rhoRange) is 2:
            rhoMin = min(rhoRange)
            rhoMax = max(rhoRange)
            #
            minidx = P4_d[0, :] >= rhoMin
            maxidx = P4_d[0, :] <= rhoMax
            #
            arein  = minidx and maxidx
            #
            P4_d = P4_d[:, arein]

        # full comps in ref cartesian frame on image (origin in lower left)
        P4_i  = P4_d + num.tile( [self.xc, self.yc, 0], (numPts, 1) ).T

        if toPixelUnits:
            xout, yout = self.pixelIndicesOfCartesianCoords(P4_i[0, :], P4_i[1, :], ROI=None)
        else:
            xout, yout = (P4_i[0, :], P4_i[1, :])

        # assign the return value here
        retval = [xout.reshape(outputShape), yout.reshape(outputShape)]
        if len(ome) > 0:
            retval.append(ome.reshape(outputShape))
        if outputDV:
            dv = num.ones(outputShape)
            retval.append(dv)
        if outputForGve:
            risoeCOB = num.dot(rotMatOfExpMap( 0.5*num.pi*num.c_[0,0,1].T ),
                               rotMatOfExpMap( 0.5*num.pi*num.c_[0,1,0].T ))
            risoeLabGvec = num.dot(risoeCOB.T, num.dot(self.ROT_l2d, P4_d)) \
                           + num.vstack([D, num.zeros((2, numPts))])
            retval = 1e3 * risoeLabGvec # these have to be in microns (Soeren)
        return retval

    def xyoToAngMap(self, x0, y0, *args, **kwargs):
        """
        eta by default is in [-pi,pi]
        if all data are in the left quadrants, remap eta into [0,2*pi]
        """
        doMap = kwargs.pop('doMap',None)
        angs = self.xyoToAng(x0, y0, *args, **kwargs)
        angs[1] = mapAngs(angs[1], doMap=doMap)
        return angs
    def makeNew(self, *args, **kwargs):
        kwargs.setdefault('getDParamDflt', self.getDParamDflt)
        kwargs.setdefault('setDParamZero', self.setDParamZero)
        kwargs.setdefault('getDParamScalings', self.getDParamScalings)
        kwargs.setdefault('getDParamRefineDflt', self.getDParamRefineDflt)
        kwargs.setdefault('radialDistortion', self.radialDistortion)
        newDG = self.__class__(self.__ncols, self.__nrows, self.__pixelPitch, 
                               self.__vFactorUnc, self.__vDark, self.reader,
                               self,
                               *args, **kwargs)
        return newDG

    def angToXYO(self, x0, y0, *args, **kwargs):
        """convert Cartesian to polar

        uses blocking to call vectorized version
        """
        #
        # Block the data if a 1D array, otherwise just call the old one.
        #
        haveOme = False
        wantDV  = False
        inShape = None
        try:
            ndim = x0.ndim
            if ndim == 2:
                inShape = x0.shape
                x0 = x0.flatten()
                y0 = y0.flatten()
            elif ndim == 0:
                return self.angToXYO_V(x0, y0, *args, **kwargs)
            lenx0 = len(x0)
        except:
            # case of nonarray as arg
            return self.angToXYO_V(x0, y0, *args, **kwargs)

        # don't forget about the pass-through ome arg
        if len(args) is not 0:
            haveOme = True
            ome = num.atleast_1d(args[0])
            if ome.ndim >= 2:
                ome = ome.flatten()
                pass
            pass

        # need this in case something asks for the dV values...
        if kwargs.has_key('outputDV'):
            wantDV = kwargs['outputDV']

        extraRows = sum([haveOme, wantDV])

        sofar = 0; tmpRetv = num.zeros((2+extraRows, lenx0))
        while lenx0 > sofar:
            #
            #  Find out how many to do
            #
            somany   = min(lenx0 - sofar, bsize)
            newsofar = sofar + somany
            x0i = x0[sofar:newsofar]
            y0i = y0[sofar:newsofar]
            if haveOme:
                rvi = self.angToXYO_V(x0i, y0i, ome[sofar:newsofar], **kwargs)
            else:
                rvi = self.angToXYO_V(x0i, y0i, **kwargs)
            for j in range(len(rvi)):
                tmpRetv[j, sofar:newsofar] = rvi[j]
                pass
            sofar = newsofar
            pass

        # ... inShape should be set properly from above
        return [tmpRetv[i, :].reshape(inShape) for i in range(tmpRetv.shape[0])]

    def xyoToAng(self, x0, y0, *args, **kwargs):
        """convert Cartesian to polar

        uses blocking to call vectorized version
        """
        #
        # Block the data if a 1D array, otherwise just call the old one.
        #
        haveOme = False
        wantDV  = False
        inShape = None
        try:
            ndim = x0.ndim
            if ndim == 2:
                inShape = x0.shape
                x0 = x0.flatten()
                y0 = y0.flatten()
            elif ndim == 0:
                return self.xyoToAng_V(x0, y0, *args, **kwargs)
            lenx0 = len(x0)
        except:
            # case of nonarray as arg
            return self.xyoToAng_V(x0, y0, *args, **kwargs)

        # don't forget about the pass-through ome arg
        if len(args) is not 0:
            haveOme = True
            ome = num.atleast_1d(args[0])
            if ome.ndim >= 2:
                ome = ome.flatten()
                pass
            pass

        # need this in case something asks for the dV values...
        if kwargs.has_key('outputDV'):
            wantDV = kwargs['outputDV']

        extraRows = sum([haveOme, wantDV])

        sofar = 0; tmpRetv = num.zeros((2+extraRows, lenx0))
        while lenx0 > sofar:
            #
            #  Find out how many to do
            #
            somany   = min(lenx0 - sofar, bsize)
            newsofar = sofar + somany
            x0i = x0[sofar:newsofar]
            y0i = y0[sofar:newsofar]
            if haveOme:
                rvi = self.xyoToAng_V(x0i, y0i, ome[sofar:newsofar], **kwargs)
            else:
                rvi = self.xyoToAng_V(x0i, y0i, **kwargs)
            for j in range(len(rvi)):
                tmpRetv[j, sofar:newsofar] = rvi[j]
                pass
            sofar = newsofar
            pass

        # ... inShape should be set properly from above
        return [tmpRetv[i, :].reshape(inShape) for i in range(tmpRetv.shape[0])]

    def xyoToAng_V(self, x0, y0, *args, **kwargs):
        """
        Convert radial spectra obtained from polar
        rebinned powder diffraction images to angular spectra.

        USAGE:
        mappedData = XFormRadialSpectra(t, D, tilt, xydata, azim, tthRange, radDistFuncHandle, radDistArgs)

        INPUTS:

        1) t is 2 x 1 (double), the origin translation.
        The convention is from `true' to `estimated' centers.
        2) D is 1 x 1 (double), the sample-to-detector distance in mm.
        3) gammaYprime is 1 x 1 (double), the angle between the `ideal'
        and `experimental' X-axes (horizontal).
        4) gammaXhatPrime is 1 x 1 (double), the angle between the
        `ideal' and `experimental' Y-axes (vertical).
        5) xydata is 1 x n (cell), the cell array of data
        6) azim
        7) tthRange
        8) radDistFuncHandle

        OUTPUT:

        1) mappedData is 1 x n (cell), the cell array of mapped radial
        data corresponding to the input `xydata'.

        """

        outputDV              = False
        tthRange              = ()
        inputPixelUnits       = True
        applyRadialDistortion = True

        outputShape = num.shape(x0)
        x0          = num.asarray(x0).flatten()
        y0          = num.asarray(y0).flatten()

        numPts = len(x0)                # ... no check for y0 or omega
        ome = ()
        if len(args) is not 0:
            ome = num.atleast_1d(args[0])

        nzeros = num.zeros(numPts)

        # argument handling
        kwarglen = len(kwargs)
        if kwarglen > 0:
            argkeys = kwargs.keys()
            for i in range(kwarglen):
                if argkeys[i] is 'outputDV':
                    outputDV = kwargs[argkeys[i]]
                elif argkeys[i] is 'units':
                    if kwargs[argkeys[i]] is 'pixels':
                        inputPixelUnits = True
                    elif kwargs[argkeys[i]] is self.pixelPitchUnit:
                        inputPixelUnits = False
                    else:
                        raise RuntimeError, 'Input units \'%s\' not understood!' % (str(kwargs[argkeys[i]]))
                elif argkeys[i] is 'tthRange':
                    tthRange = kwargs[argkeys[i]]
                    assert len(tthRange) is 2, 'Two-theta range should have length 2'
                elif argkeys[i] is 'rdist':
                    if not isinstance(kwargs[argkeys[i]], bool):
                        raise RuntimeError, 'Expecting boolean for rdist kewyord argument; got' \
                              + str(applyRadialDistortion)
                    else:
                        applyRadialDistortion = kwargs[argkeys[i]]
                else:
                    raise RuntimeError, 'Unrecognized keyword argument: ' + str(argkeys[i])

        # make center-based cartesian coord's
        #   - SHOULD BE IN PIXEL PITCH UNITS (MM)
        #   ... maybe add hard check on this in future
        #   - x0, y0 are now written in the cartesian coords
        #     where [0, 0] is the lower left corner of detector
        if inputPixelUnits:
            x0, y0 = self.cartesianCoordsOfPixelIndices(x0, y0)

        if self.pVec is None:
            X_d = x0 - self.xc              # is 1-d!
            Y_d = y0 - self.yc              # is 1-d!

            XC = nzeros
            YC = nzeros

            D = num.tile(self.workDist, numPts)
        else:
            assert len(ome) == numPts, 'with precession, omega argument consistent with ' \
                   + 'x and y (or i and j) is required'
            # if here, we have a precession vector and must deeal with it
            #
            #   - ome is taken as a CCW (+) rotation of the SAMPLE FRAME about Y
            #   - when the BASIS is transformed by R, vector comp's must transform by R'
            R_s2l = rotMatOfExpMap( num.tile(ome, (3, 1)) * num.tile(self.Yl, (1, numPts)) )

            if not hasattr(self.pVec, '__len__'):
                raise RuntimeError, 'pVec must be array-like'

            self.pVec = num.asarray(self.pVec)

            # array of rotated precession vector components
            grainCM_l = num.dot(R_s2l, self.pVec.reshape(3, 1))
            if grainCM_l.ndim == 3:
                grainCM_l = grainCM_l.squeeze().T

            # precession-corrected polar detector coord's
            # X_d = x0 - (self.xc + grainCM_l[0, :]) # is 1-d!
            # Y_d = y0 - (self.yc + grainCM_l[1, :]) # is 1-d!
            X_d = x0 - self.xc          # is 1-d!
            Y_d = y0 - self.yc          # is 1-d!

            XC = grainCM_l[0, :]
            YC = grainCM_l[1, :]

            D = self.workDist + grainCM_l[2, :] # now array of D's

        if applyRadialDistortion:
            # apply distortion
            X_d, Y_d = self.radialDistortion(X_d, Y_d, invert=False)

        #
        # ------- ASSIGN POINT COORD'S AND FORM ROTATION
        #
        # origins of the scattering (P1) and lab (P2) frames
        #   - the common frame for arithmatic is the scattering frame
        P1 = num.vstack([XC, YC, D])
        # P2 = num.vstack([XC, YC, nzeros])
        P2 = num.zeros((3, numPts))

        # tilt calculations moved into setTilt

        # full 3-d components in tilted the detector frame
        P4_d = num.vstack( (X_d, Y_d, nzeros) )

        # rotate components into the lab frame
        P4_l = num.dot(self.ROT_l2d, P4_d)

        # apply translation to get equations of diffracted rays in lab frame
        rays = P4_l - P1

        # solve for P3 coord's in lab frame
        u = num.tile( num.dot(self.N.T, (P2 - P1)) / num.dot(self.N.T, rays), (3, 1) )


        P3 = P1 + u * rays

        # X-Y components of P3 in lab frame
        X_l = P3[0, :] - XC
        Y_l = P3[1, :] - YC

        # polar coords in lab frame
        rho_l = num.sqrt(X_l*X_l + Y_l*Y_l)
        eta_l = num.arctan2(Y_l, X_l)

        # get two-theta from dot products with lab-frame beam direction
        dotProds = num.dot(-self.Zl.T, unitVector(rays)).squeeze()

        # two-theta
        measTTH = arccosSafe(dotProds)

        # transform data
        tmpData = num.vstack( [measTTH, eta_l] )

        if len(tthRange) is 2:
            tthMin = min(tthRange)
            tthMax = max(tthRange)

            minidx = r2d*tmpData[0, :] >= tthMin
            maxidx = r2d*tmpData[0, :] <= tthMax

            arein  = minidx and maxidx

            tmpData = tmpData[:, arein]

        retval = [tmpData[0, :].reshape(outputShape), tmpData[1, :].reshape(outputShape)]
        if len(ome) > 0:
            retval.append(ome.reshape(outputShape))

        if outputDV:
            dv = num.ones(outputShape)
            retval.append(dv)

        return retval

    def makeMaskTThRanges(self, planeData):
        """
        Mask in the sense that reader with the mask will exclude all else
        """
        indicesList, iHKLLists = self.makeIndicesTThRanges(planeData)
        mask = -self.reader.indicesToMask(indicesList)
        return mask
    def xyoToAngAll(self):
        """
        get angular positions of all pixels
        """
        jVals = num.tile(num.arange(self.__ncols),(self.__nrows,1))
        iVals = jVals.T
        twoTheta, eta = self.xyoToAng(iVals, jVals)
        return twoTheta, eta
    def xyoToAngCorners(self):
        """
        get angular positions of corner pixels
        """
        iVals = num.array([0, self.__nrows-1, 0, self.__nrows-1])
        jVals = num.array([0, 0, self.__ncols-1, self.__ncols-1])
        twoTheta, eta = self.xyoToAng(iVals, jVals)
        return twoTheta, eta
    def angOnDetector(self, tTh, eta, *args):
        '''
        note: returns a scalar if tTh and eta have single entries
        '''
        i, j = self.angToXYO(tTh, eta)
        retval = num.logical_and (
            num.logical_and( i >= 0, i <= self.nrows-1 ),
            num.logical_and( j >= 0, j <= self.ncols-1 ) )
        retval = num.atleast_1d(retval)
        return retval
    def makeTThRanges(self, planeData, cullDupl=False):
        tThs      = planeData.getTTh()
        tThRanges = planeData.getTThRanges()

        nonoverlapNexts = num.hstack((tThRanges[:-1,1] < tThRanges[1:,0], True))
        iHKLLists = []
        hklsCur = []
        for iHKL, nonoverlapNext in enumerate(nonoverlapNexts):
          if not nonoverlapNext:
            if cullDupl and abs(tThs[iHKL] - tThs[iHKL+1]) < 1e-8:
              'do not append, is a duplicate'
              pass
            else:
              hklsCur.append(iHKL)
          else:
            hklsCur.append(iHKL)
            iHKLLists.append(hklsCur)
            hklsCur = []

        return iHKLLists
    def makeIndicesTThRanges(self, planeData, cullDupl=False):
        """
        return a list of indices for sets of overlaping two-theta ranges;
        to plot, can do something like:
        	mask = self.reader.getEmptyMask()
          mask[indices] = True

        With cullDupl set true, eliminate HKLs with duplicate 2-thetas
        """
        tThs      = planeData.getTTh()
        tThRanges = planeData.getTThRanges()

        iHKLLists = self.makeTThRanges(planeData, cullDupl=cullDupl)

        indicesList = []
        twoTheta, eta = self.xyoToAngAll()
        for iHKLList in iHKLLists:
            'for some reason, this does not work well when made into a one-liner'
            #indices = num.where(num.logical_and(twoTheta > tThRanges[iHKLList[0],0], twoTheta < tThRanges[iHKLList[-1],1]))
            b1 = twoTheta > tThRanges[iHKLList[0],0]
            b2 = twoTheta < tThRanges[iHKLList[-1],1]
            b = num.logical_and(b1,b2)
            indices = num.where(b)
            indicesList.append(indices)
        return indicesList, iHKLLists
    def getTThMax(self, func=num.min):
        x0_max = self.ncols-1;   y0_max = self.nrows-1;
        x0_mid = 0.5 * x0_max;   y0_mid = 0.5 * y0_max;
        x0_test = num.array([     0., x0_max, x0_mid, x0_mid, 0., x0_max,     0., x0_max ])
        y0_test = num.array([ y0_mid, y0_mid,     0., y0_max, 0., y0_max, y0_max,     0. ])
        tTh_test, eta_test = self.xyoToAng(x0_test, y0_test)
        return func(tTh_test)

    def getRings(self, planeData, ranges=False):
        """
        Return a list of rings for the given hkls

        Already filters on the exclusions.
        """
        rList = []
        nEta = self.nEta
        dEta = 2*math.pi/nEta

        etaRing = num.arange(0., 2.0*math.pi+dEta/2., dEta)
        nEta    = len(etaRing) # why this?

        excl   = num.array(planeData.exclusions)

        # grab full tThs and tile if ranges are desired
        tThs = planeData.getTTh()
        if ranges:
            tThs = planeData.getTThRanges().flatten()

        # grab the relevant hkls and loop
        for i in range(len(tThs)):
            tThRing = num.tile(tThs[i], nEta)
            r = self.angToXYO(tThRing, etaRing) # omegaRing
            rList.append(r)
            pass
        return rList

    def getPRBOverlay(self, polarRebinKWArgs):
        """
        Return plottable coordinates of rebinning sector.

        Takes in dictionary of keyword args for polarRebin

        for etas, right now assumes stopEta > startEta, CCW
        """
        startEta = polarRebinKWArgs['etaRange'][0]
        stopEta  = polarRebinKWArgs['etaRange'][1]
        numEta   = polarRebinKWArgs['numEta']

        startRho = polarRebinKWArgs['rhoRange'][0]*self.pixelPitch
        stopRho  = polarRebinKWArgs['rhoRange'][1]*self.pixelPitch

        nrows = self.nrows   # total number of rows in the full image
        ncols = self.ncols   # total number of columns in the full image

        nEtaRing = round(self.nEta * abs(stopEta - startEta) / (2*num.pi))
        dEta     = abs(stopEta - startEta) / nEtaRing

        # DEBUGGING # print 'nEtaRing: %d' %(nEtaRing)
        # DEBUGGING # print 'dEta: %f' %(dEta * 180. / num.pi)

        # this is the generic ring segment
        etaRing = num.arange(startEta, stopEta + 0.5*dEta, dEta)
        # print etaRing, startRho, stopRho

        # arc segments as [rho, eta] pairs (vstacked)
        innerArc_pol = num.vstack([startRho * num.ones(len(etaRing)), etaRing])
        outerArc_pol = num.vstack([stopRho  * num.ones(len(etaRing)), etaRing])
        innerArc_pixI, innerArc_pixJ = self.polarToPixel(innerArc_pol[0, :],
                                                         innerArc_pol[1, :],
                                                         corrected=polarRebinKWArgs['corrected'])
        outerArc_pixI, outerArc_pixJ = self.polarToPixel(outerArc_pol[0, :],
                                                         outerArc_pol[1, :],
                                                         corrected=polarRebinKWArgs['corrected'])
        # print innerArc_pol, outerArc_pol
        del(innerArc_pol)
        del(outerArc_pol)

        sector_dEta = abs(stopEta - startEta) / numEta
        sectorEtas = num.arange(startEta, stopEta + 0.5*sector_dEta, sector_dEta).tolist()
        # sector edges
        if numEta == 1:
            edge_pixI, edge_pixJ = self.polarToPixel([startRho,  stopRho, startRho, stopRho],
                                                     [startEta, startEta,  stopEta, stopEta],
                                                     corrected=polarRebinKWArgs['corrected'])
            edge_pixL = zip(edge_pixI.reshape(2,2).tolist(),
                            edge_pixJ.reshape(2,2).tolist())
        else:
            edge_pixI = []
            edge_pixJ = []
            for isec in range(numEta+1):
                tEta  = sectorEtas[isec]
                tmpI, tmpJ = self.polarToPixel([startRho,  stopRho],
                                               [    tEta,     tEta],
                                               corrected=polarRebinKWArgs['corrected'])
                edge_pixI = edge_pixI + tmpI.tolist()
                edge_pixJ = edge_pixJ + tmpJ.tolist()
                pass
            edge_pixI = num.array(edge_pixI).reshape(numEta+1, 2)
            edge_pixJ = num.array(edge_pixJ).reshape(numEta+1, 2)
            edge_pixL = zip(edge_pixI.tolist(),
                            edge_pixJ.tolist())
            pass

        # DEBUGGING # import pdb; pdb.set_trace()
        # pixI = innerArc_pixI.tolist() + outerArc_pixI.tolist() + edge_pixI
        # pixJ = innerArc_pixJ.tolist() + outerArc_pixJ.tolist() + edge_pixJ
        retval = [edge_pixL[i] for i in range(len(edge_pixL))]
        retval.append([innerArc_pixI.tolist(), innerArc_pixJ.tolist()])
        retval.append([outerArc_pixI.tolist(), outerArc_pixJ.tolist()])
        return retval

    def angToXYOBBox(self, *args, **kwargs):
        """
        given either angBBox or angCOM (angular center) and angPM (+-values), compute the bounding box on the image frame
        
        if forSlice=True, then returned bbox is appropriate for use in array slicing
        """
        
        units    = kwargs.setdefault('units', 'pixels')
        #
        # reader = kwargs.get('reader', None)
        reader = None
        if kwargs.has_key('reader'):
            reader = kwargs.pop('reader')
        #
        omegas = None
        if kwargs.has_key('omegas'):
            omegas = kwargs.pop('omegas')
        #
        forSlice = True
        if kwargs.has_key('forSlice'):
            forSlice = kwargs.pop('forSlice')
        slicePad = 0
        if forSlice:
            slicePad = 1
        
        if len(args) == 1:
            angBBox = args[0]
        elif len(args) == 2:
            angCOM, angPM = args
            angBBox = (
                (angCOM[0]-angPM[0], angCOM[0]+angPM[0]),
                (angCOM[1]-angPM[1], angCOM[1]+angPM[1]),
                (angCOM[2]-angPM[2], angCOM[2]+angPM[2]),
                )
        else:
            raise RuntimeError, 'specify either angBBox or angCOM, angPM'

        'along eta, try multiple points if the spread is wide enough so that do not mis apex of arc'
        delAng = 0.1 # about 6 degrees
        nTest = max(int(math.ceil((angBBox[1][1]-angBBox[1][0]) / delAng)), 1) + 1
        etas = num.hstack( (num.linspace(angBBox[1][0], angBBox[1][1], nTest), angBBox[1][0], angBBox[1][1] ) )
        tths = num.hstack( (num.tile(angBBox[0][1],                    nTest), angBBox[0][0], angBBox[0][0] ) )
        x, y = self.angToXYO(tths, etas, **kwargs)

        xyoBBox = [
            ( x.min(), x.max() ),
            ( y.min(), y.max() ),
            angBBox[2],
            ]
        if units == 'pixels':
            'make into integers'
            xyoBBox[0] = ( max( int(math.floor(xyoBBox[0][0])), 0), 
                           min( int(math.floor(xyoBBox[0][1])), self.nrows-1)+slicePad,
                           )
            xyoBBox[1] = ( max( int(math.floor(xyoBBox[1][0])), 0), 
                           min( int(math.floor(xyoBBox[1][1])), self.ncols-1)+slicePad,
                           )
        if reader is not None:
            'convert bounding box from omegas to frames'
            xyoBBox[2] = ( num.hstack( (reader.omegaToFrameRange(xyoBBox[2][0]), 0) )[0],
                           num.hstack( (reader.getNFrames()-1, reader.omegaToFrameRange(xyoBBox[2][1]) ) )[-1] + slicePad,
                           )
        elif omegas is not None:
            'convert bounding box from omegas to frames'
            omegaDelta = num.mean(omegas[1:]-omegas[:-1])
            nFrames = len(omegas)
            xyoBBox[2] = ( 
                num.hstack( (omeToFrameRange(xyoBBox[2][0], omegas, omegaDelta), 0) )[0],
                num.hstack( (nFrames-1, omeToFrameRange(xyoBBox[2][1], omegas, omegaDelta) ) )[-1] + slicePad,
                )

        return xyoBBox

    def drawRings(self, drawOn, planeData, withRanges=False, legendLoc=(0.05,0.5), legendMaxLen=10,
                  ideal=None, range=None, lineType=None, lineWidth=1.0):
        """
        If drawOn is a PlotWrap instance, draw on the existing instance,
        otherwise pass drawOn to display and return the resulting
        PlotWrap instance

        planeData.exclusions can be used to work with a subset of rings;

        set legendLoc to None or False to skip making the legend

        removes any existing lines in the axes

        if pass ideal, then display rings on an ideal detector with the working distance taken from the value of the ideal argument
        """
        nEta = self.nEta
        dEta = 2*math.pi/nEta

        if ideal is None:
            workDist = self.workDist
            angToXY = self.angToXYO
        else:
            workDist = ideal
            angToXY = lambda tTh_l, eta_l: angToXYIdeal(tTh_l, eta_l, workDist)

        if isinstance(drawOn, plotwrap.PlotWrap):
            pw = drawOn
            retval = None
        else:
            if ideal is None:
                pw = self.display(drawOn, planeData=planeData, range=range)
                retval = pw
            else:
                pw = self.displayIdeal(drawOn, planeData=planeData, workDist=workDist, range=range)
                retval = pw

        'get rid of any existing lines and legends'
        pw.removeLines()

        tThs      = planeData.getTTh()
        hkls      = planeData.getHKLs(asStr=True)

        if lineType is not None:
            lineStyles = LineStyles(lt=lineType)
        else:
            lineStyles = LineStyles()

        etaRing   = num.arange(0., 2.0*math.pi+dEta/2., dEta)
        nEta      = len(etaRing)
        # omegaRing = num.zeros([nEta]) # 'omega is arbitrary'
        pw.showByDefault = False
        linesForLegend = []
        for tTh in tThs:
            tThRing = num.tile(tTh, nEta)
            x,y = angToXY(tThRing, etaRing)
            pw(x,y, style=lineStyles(), **self.lineArgs)
            linesForLegend.append(pw.a.get_lines()[-1])
        if withRanges:
          tThRanges = planeData.getTThRanges()
          for tTh in tThRanges.flatten():
            tThRing = num.tile(tTh, nEta)
            x,y = angToXY(tThRing, etaRing)
            pw(x,y, ls='-', color='1.0', **self.lineArgs)
        if legendLoc:
            'remove of legends not supported'
            # for legend in pw.f.legends:
            #     legend.remove()
            #
            #pw.a.legend(linesForLegend, hkls, loc=legendLoc)
            if legendMaxLen and len(linesForLegend) > legendMaxLen:
                linesForLegend = linesForLegend[0:legendMaxLen]
                hkls           = hkls[0:legendMaxLen]
                hkls[-1]       = "..."
            pw.f.legend(linesForLegend, hkls, legendLoc,
                        handlelength=1, handletextpad=0.4, borderaxespad=1, labelspacing=0.1,
                        borderpad=0.1, numpoints=1)
        pw.showByDefault = True
        pw.show()
        return retval
    def renderIdeal(self, thisframe, nlump=None,
                    workDist=None):
        """
        render the frame on an ideal detector plane;
        returns interpolated frame data zi on a regular grid xi, yi;
        suitable for use with pcolormesh(xim, yim, zi), with xim, yim = num.meshgrid(xi, yi);
        note that pcolormesh is used instead of pcolor because zi may be a masked array
        """

        nlump = nlump or 4
        workDist = workDist or self.workDist
        assert nlump >= 2,\
            'due to histogram2d based method, this make not work well for nlump < 2'

        nx = float(self.shape[0])/float(nlump)+1
        ny = float(self.shape[1])/float(nlump)+1

        # 'figure out the range for plotting, assuming the corners bound things well enough'
        # tthBox, etaBox = self.xyoToAngCorners()

        'create data on the ideal plane, ends up being irregularly spaced'
        tTh, eta = self.xyoToAngAll()
        x, y = angToXYIdeal(tTh, eta, workDist)
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        nx = int(nx * (xmax-xmin)/self.extent[0])
        ny = int(ny * (ymax-ymin)/self.extent[1])

        'Generate a regular grid to interpolate the data.'
        xi = num.linspace(xmin, xmax, nx)
        yi = num.linspace(ymin, ymax, ny)

        'interpolate'
        # have tried other things, like interpolate.interp2d and mlab.griddata without much luck
        # xim, yim = num.meshgrid(xi, yi)
        # zi = mlab.griddata(x.flatten(),y.flatten(),thisframe.flatten(),xim,yim)
        h, xedges, yedges = num.histogram2d(
            x.flatten(), y.flatten(), bins=[xi,yi], weights=thisframe.flatten())
        c, xedges, yedges = num.histogram2d(x.flatten(), y.flatten(), bins=[xi,yi])
        mask = c == 0
        h[-mask] = h[-mask] / c[-mask]

        return h, xedges, yedges, mask

    def displayIdeal(self, thisframe, planeData=None,
                     workDist=None, nlump=None,
                     **kwargs):
        """
        render and display frame on ideal detector plane;
        if workDist is not specified, then use self.workDist
        """

        h, xedges, yedges, mask = self.renderIdeal(thisframe, workDist=workDist, nlump=nlump)
        workDist = workDist or self.workDist

        if kwargs.has_key('pw'):
            pw = kwargs.pop('pw')
        else:
            pwKWArgs = plotwrap.PlotWrap.popKeyArgs(kwargs)
            pw = plotwrap.PlotWrap(**pwKWArgs)
        retval = pw
        
        vmin, vmax, cmap = self.reader.getDisplayArgs(h, kwargs)
        pw.a.set_axis_bgcolor('white')
        cmap.set_under(color='white', alpha=0.0)
        norm = matplotlib.colors.Normalize(clip=False, vmin=vmin, vmax=vmax)
        if True:
            '''
            hack around to make for better plotting;
            needed due to tolerancing in cmap for when a color is under?
            '''
            h[h < 0] =  vmin
            h[mask]  = -vmax
        pw.a.pcolor(xedges, yedges, h.T, cmap=cmap, norm=norm, **kwargs)
        pw.a.set_aspect('equal')
        pw.a.set_autoscale_on(False)
        pw.show()

        fmtCoord = FmtCoordIdeal(planeData, workDist)
        fmtCoord.addDetectorData(1, xedges, yedges, h, mask)
        pw.a.format_coord = fmtCoord

        return retval

    def display(self, thisframe, planeData=None, **kwargs):
        """
        wrap reader display method;
        display coordinates as 2-theta and eta given that self knows how to do this

        if pass planeData, then it is used to list HKLs overlapping the given 2-theta location

        ...*** option for drawing lab-frame glyph
        """
        pw = self.reader.display(thisframe, **kwargs)
        # pw.a is same as pw.win.getAxes(0)
        if planeData is None:
            def fmtCoord(x,y):
                tTh, eta = num.array(self.xyoToAng(y, x))
                return "tth=%g eta=%g int=%g" \
                    % (r2d*tTh, r2d*eta, thisframe[round(y), round(x)])
            pw.a.format_coord = fmtCoord
        else:
            def fmtCoord(x,y):
                tTh, eta = num.array(self.xyoToAng(y, x))
                cartX, cartY = self.cartesianCoordsOfPixelIndices(y, x)
                cx = (cartX - self.xc)/self.pixelPitch
                cy = (cartY - self.yc)/self.pixelPitch
                rho      = num.sqrt(cx*cx + cy*cy)
                dsp      = 0.5 * planeData.wavelength / num.sin(0.5*tTh)
                HKLs     = str(planeData.getHKLs(asStr=True, allHKLs=True, thisTTh=tTh))
                return "d=%g tth=%g eta=%g int=%g \n HKLs=%s" \
                       % (dsp, r2d*tTh, r2d*eta, thisframe[round(y), round(x)], HKLs)
            pw.a.format_coord = fmtCoord
        return pw

    def drawRingsGUI(self, thisframe, planeData, displayKWArgs={}, sliderRangeFactor=1.0, funcType=funcTypeDflt):
        """
        a simple GUI
        """
        doDragging = True
        sliderYArea = 0.18+0.03
        buttonsXArea = 0.15
        'pass planeData to display so that HKLs for the 2-theta at the cursor position are shown'
        pw = self.display(thisframe, figsize=(7,7), planeData=planeData, **displayKWArgs)
        pw.a.set_position([buttonsXArea, sliderYArea, 1.0-buttonsXArea, 1.0-sliderYArea], which='original')

        delxy = min(self.__nrows, self.__ncols) * 0.1 * sliderRangeFactor * self.pixelPitch
        axcolor = 'lightgoldenrodyellow'
        #
        rect_cur   = [0.125, sliderYArea-0.03,  0.75, 0.02]
        rect_dy    = -0.03
        #
        ax_xc    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_xc     = Slider(ax_xc, 'xc', self.xc-delxy, self.xc+delxy, valinit=self.xc, dragging=doDragging)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        ax_yc    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_yc     = Slider(ax_yc, 'yc', self.yc-delxy, self.yc+delxy, valinit=self.yc, dragging=doDragging)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        ref      = self.workDist
        delwd    = 2.0*sliderRangeFactor
        ax_wd    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_wd     = Slider(ax_wd, 'D ', ref/delwd, ref*delwd, valinit=ref, dragging=doDragging)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        deltilt  = max(0.2*sliderRangeFactor, abs(self.xTilt)*1.2, abs(self.yTilt)*1.2)
        #
        ref      = 0.
        ax_xt    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_xt     = Slider(ax_xt, 'xt', ref-deltilt, ref+deltilt, valinit=ref, dragging=doDragging)
        s_xt.set_val(self.xTilt)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        ref      = 0.
        ax_yt    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_yt     = Slider(ax_yt, 'yt', ref-deltilt, ref+deltilt, valinit=ref, dragging=doDragging)
        s_yt.set_val(self.yTilt)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        if planeData.tThWidth is None:
            self.withTThWidth = False
            ref = planeData.strainMag
        else:
            self.withTThWidth = True
            ref = planeData.tThWidth
            if hasattr(ref, 'getVal'):
                ref = ref.getVal('radians')
        deltw    = 10.0*sliderRangeFactor
        ax_tw    = pw.f.add_axes(rect_cur, axisbg=axcolor)
        s_tw     = Slider(ax_tw, 'width', ref/deltw, ref*deltw, valinit=ref, dragging=doDragging)
        rect_cur[1] = rect_cur[1] + rect_dy
        #
        ax_reset = pw.f.add_axes([0.03, 0.9, buttonsXArea-0.02, 0.04])
        b_reset  = Button(ax_reset, 'Reset')
        #
        ax_range = pw.f.add_axes([0.03, 0.8, buttonsXArea-0.02, 0.09])
        b_range  = RadioButtons(ax_range, ['Off','Ranges'], active=int(self.withRanges), activecolor='blue')
        #
        ax_stw   = pw.f.add_axes([0.03, 0.7, buttonsXArea-0.02, 0.09])
        b_stw    = RadioButtons(ax_stw, ['Off','tThWidth'], active=int(self.withTThWidth), activecolor='blue')
        #
        self.asMasked = False
        ax_mask  = pw.f.add_axes([0.03, 0.65, buttonsXArea-0.02, 0.04])
        #b_mask   = RadioButtons(ax_mask, ['Off','Mask'], active=int(self.asMasked), activecolor='blue')
        b_mask   = Button(ax_mask, 'Mask')
        #
        ax_fit   = pw.f.add_axes([0.03, 0.6, buttonsXArea-0.02, 0.04])
        b_fit    = Button(ax_fit, 'Fit')
        #
        legendLoc = (0.03, sliderYArea+0.02) # 'center left'

        def undo_mask():
            self.asMasked = False
            self.display(thisframe,  pw=pw, **displayKWArgs)
            self.drawRings(pw, planeData,
                           withRanges=self.withRanges, legendLoc=None)
        def update_self(val):
            if self.asMasked:
              undo_mask()
            self.xc = s_xc.val
            self.yc = s_yc.val
            self.workDist = s_wd.val
            self.xTilt = s_xt.val
            self.yTilt = s_yt.val
            self.drawRings(pw, planeData,
                           withRanges=self.withRanges, legendLoc=None)
            pw.show()
        def update_range(val):
            if val == 'Off':
              self.withRanges = False
            else:
              self.withRanges = True
            self.drawRings(pw, planeData,
                           withRanges=self.withRanges, legendLoc=None)
            pw.show()
        def update_stw(val):
            if val == 'Off':
              self.withTThWidth = False
            else:
              self.withTThWidth = True
            'call update_tw to redo everything needed'
            update_tw(None)
        def do_mask(event):
            self.asMasked = True
            mask = -self.makeMaskTThRanges(planeData)
            maskedframe = copy.deepcopy(thisframe)
            maskedframe[mask] = 0
            self.display(maskedframe, pw=pw, **displayKWArgs)
            self.drawRings(pw, planeData,
                           withRanges=self.withRanges, legendLoc=None)
            pw.show()
        def update_tw(val):
            if self.withTThWidth:
                planeData.tThWidth  = s_tw.val
            else:
                planeData.strainMag = s_tw.val
            if self.asMasked:
              'call do_mask so that masking is redone'
              do_mask(None)
            self.drawRings(pw, planeData,
                           withRanges=self.withRanges, legendLoc=None)
            pw.show()
        def do_fit(event):
            """ ... consider only allowing fit if b_range is 'Ranges', so that know
            the user at least could have done the sanity check to see that rings are
            covered;
            or could code 'smart' adjustment of strainMag -- increase until stuff that
            is being added looks like background
            """
            self.fitRings(thisframe, planeData, funcType=funcType)
            '''calling set_val on the sliders causes trouble;
            inside of set_val they call their observers, so update_self
            which do not want called before slide values are set;
            but the following is not ideal either in that it sets the
            values in the sliders without using a method on them'''
            s_xc.val = self.xc # s_xc.set_val(self.xc)
            s_yc.val = self.yc # s_yc.set_val(self.yc)
            s_wd.val = self.workDist # s_wd.set_val(self.workDist)
            s_xt.val = self.xTilt # s_xt.set_val(self.xTilt)
            s_yt.val = self.yTilt # s_yt.set_val(self.yTilt)
            update_self(None)
        def do_reset(event):
            if self.asMasked:
              undo_mask()
            s_xc.reset()
            s_yc.reset()
            s_wd.reset()
            s_xt.reset()
            s_yt.reset()
            s_tw.reset()

        s_xc.on_changed(update_self)
        s_yc.on_changed(update_self)
        s_wd.on_changed(update_self)
        s_xt.on_changed(update_self)
        s_yt.on_changed(update_self)
        s_tw.on_changed(update_tw)
        b_range.on_clicked(update_range)
        b_stw.on_clicked(update_stw)
        b_mask.on_clicked(do_mask)
        b_reset.on_clicked(do_reset)
        b_fit.on_clicked(do_fit)

        'only draw legend the first time as legend.remove does not currently work'
        self.drawRings(pw, planeData,
                       withRanges=self.withRanges, legendLoc=legendLoc)

        pw.show()
        matplotlib.interactive(True)
        return pw
    def clean(self):
        # self.fitRingsFunc = None
        self.xFitRings = None
        return
    def fitRings(self, thisframe, planeData, xVec0=None,
                 funcType=funcTypeDflt, quadr=1, makePlots=False):

      # 'free up memory'
      # self.fitRingsFunc = None

      func = MultiRingEval(self, planeData, dataFrame=thisframe,
                           funcType=funcType, quadr=quadr)

      if xVec0 is None:
          if self.xFitRings is not None and len(self.xFitRings) == func.getNParam():
              'use previous fit as starting point'
              xVec0 = copy.deepcopy(self.xFitRings)
          else:
              xVec0 = func.guessXVec()
      self.xFitRings = None
      
      x = func.doFit(xtol=DFLT_XTOL)
      
      self.xFitRings = x
      # self.fitRingsFunc = func
      'call func one more time to make sure that parameters are set to values from x solution'
      if makePlots:
          funcEval = func(x, makePlots=makePlots,
                          plotTitlePrefix='(auxiliary, not of prime importance for fits!) ') # self.updateParams(x[range(func.nParamBase)])
      print 'fit detector parameters : ' + str(self.getParams()) + '\n'

      return

    def pixelToPolar(self, rowInd, colInd, corrected=False, startEta=None):

        # WTF?! # if ROI is None:
        # WTF?! #     ROI = [0, 0, self.nrows, self.ncols]  # integer pixels indices
        # WTF?! #
        # WTF?! # rowInd = ROI[0] + num.arange(ROI[2])
        # WTF?! # colInd = ROI[1] + num.arange(ROI[3])

        pixelGrid = num.meshgrid( num.atleast_1d(rowInd),
                                  num.atleast_1d(colInd) )
        pixelIs = pixelGrid[0].T.flatten()
        pixelJs = pixelGrid[1].T.flatten()

        # convert to cartesian frame
        #     - ouput is in self.pixelPitchUnit
        #     - origin is LOWER LEFT CORNER
        # WTF ?! # x0, y0 = self.cartesianCoordsOfPixelIndices(pixelIs, pixelJs, ROI=ROI[:2])
        x0, y0 = self.cartesianCoordsOfPixelIndices(pixelIs, pixelJs)


        if corrected:

            # do conversion to tTh
            tTh, eta = self.xyoToAng(x0, y0, units=self.pixelPitchUnit)

            # get rho in ideal frame
            rho = self.workDist * num.tan(tTh.flatten())
            eta = eta.flatten()
            if startEta is not None:
                eta = mapAngle(eta, [startEta, 2*num.pi + startEta], units='radians')


            x = rho * num.cos(eta)
            y = rho * num.sin(eta)
        else:
            # move to center
            x = (x0 - self.xc).flatten()
            y = (y0 - self.yc).flatten()

            # convert to polar
            #   - map eta into specified monotonic range
            rho = num.sqrt( x*x + y*y )
            eta = num.arctan2(y, x)
            if startEta is not None:
                eta = mapAngle(eta, [startEta, 2*num.pi + startEta], units='radians')

        return rho, eta, x, y
    def polarToPixel(self, rho, eta, corrected=False):
        # WTF ?! # if not ROI is None:
        # WTF ?! #     ROI = ROI[:2]
        if corrected:
            tth = num.arctan2(rho, self.workDist)
            x, y = self.angToXYO(tth, eta, units=self.pixelPitchUnit)
        else:
            x = rho * num.cos(eta) + self.xc
            y = rho * num.sin(eta) + self.yc
        # WTF ?! # pixelI, pixelJ = self.pixelIndicesOfCartesianCoords(x, y, ROI=ROI)
        pixelI, pixelJ = self.pixelIndicesOfCartesianCoords(x, y)
        return pixelI, pixelJ

    def polarRebin(self, thisFrame,
                   npdiv=2,
                   rhoRange=[100, 1000],
                   numRho=1200,
                   etaRange=num.pi*num.r_[-5, 355]/180.,
                   numEta=36,
                   ROI=None,
                   corrected=False,
                   verbose=True,
                   log=None
                   ):
        """
        Caking algorithm

        INPUTS

        thisFrame
        npdiv=2, pixel subdivision (n x n) to determine bin membership
        rhoRange=[100, 1000] - radial range in pixels
        numRho=1200 - number of radial bins
        etaRange=num.pi*num.r_[-5, 355]/180. -- range of eta
        numEta=36 - number of eta subdivisions
        ROI=None - region of interest (four vector)
        corrected=False - uses 2-theta instead of rho
        verbose=True,

        """

        startEta = etaRange[0]
        stopEta  = etaRange[1]

        startRho = rhoRange[0]*self.pixelPitch
        stopRho  = rhoRange[1]*self.pixelPitch

        nrows = self.nrows   # total number of rows in the full image
        ncols = self.ncols   # total number of columns in the full image

        subPixArea = 1/float(npdiv)**2 # areal rescaling for subpixel intensities

        # import pdb;pdb.set_trace()

        if ROI is None:
            ROI = [0, 0, nrows, ncols]  # integer pixels indices

        # MASTER COORDINATES
        #   - in pixel indices, UPPER LEFT PIXEL is [0, 0] --> (row, col)
        #   - in fractional pixels, UPPER LEFT CORNER is [-0.5, -0.5] --> (row, col)
        #   - in cartesian frame, the LOWER LEFT CORNER is [0, 0] --> (col, row)

        rowInd = ROI[0] + num.arange(ROI[2])
        colInd = ROI[1] + num.arange(ROI[3])

        # ROI data in proper shape
        roiData = thisFrame[num.ix_(rowInd, colInd)].flatten()

        rho, eta, x, y = self.pixelToPolar(rowInd, colInd, corrected=corrected, startEta=startEta)


        # MAKE POLAR BIN CENTER ARRAY
        deltaEta = (stopEta - startEta) / numEta
        deltaRho = (stopRho - startRho) / numRho

        rowEta = startEta + deltaEta * ( num.arange(numEta) + 0.5 )
        colRho = startRho + deltaRho * ( num.arange(numRho) + 0.5 )

        # initialize output dictionary
        polImg = {}
        polImg['radius']    = colRho
        polImg['azimuth']   = rowEta
        polImg['intensity'] = num.empty( (numEta, numRho) )
        polImg['deltaRho']  = deltaRho


        if verbose:
            msg = "INFO: Masking pixels\n"
            if log:
                log.write(msg)
            else:
                print msg
                pass

        rhoI = startRho - 1
        rhoF = stopRho + 1
        inAnnulus = num.where( (rho >= rhoI) & (rho <= rhoF) )[0]
        for i in range(numEta):
            if verbose:
                msg = "INFO: Processing sector %d of %d\n" % (i+1, numEta)
                if log:
                    log.write(msg)
                else:
                    print msg
                    pass

            # import pdb;pdb.set_trace()
            etaI = rowEta[i] - 0.5*deltaEta
            etaF = rowEta[i] + 0.5*deltaEta

            tmpEta = eta[ inAnnulus ]
            inSector = num.where( (tmpEta >= etaI) & (tmpEta <= etaF) )[0]

            nptsIn = len(inSector)

            tmpX = x[ inAnnulus[inSector] ]
            tmpY = y[ inAnnulus[inSector] ]
            tmpI = roiData[ inAnnulus[inSector] ]

            # import pdb;pdb.set_trace()
            # subdivide pixels
            #   - note that these are in fractional pixel coordinates (centered)
            #   - must convert to working units (see 'self.pixelPitchUnits')
            subCrds    = (num.arange(npdiv) + 0.5) / npdiv

            intX, intY = num.meshgrid(subCrds, subCrds)

            intX = num.tile(intX.flatten(), (nptsIn, 1)).T
            intY = num.tile(intY.flatten(), (nptsIn, 1)).T

            # expand coords using pixel subdivision
            tmpX = num.tile(tmpX, (npdiv**2, 1)) + (intX - 0.5)*self.pixelPitch
            tmpY = num.tile(tmpY, (npdiv**2, 1)) + (intY - 0.5)*self.pixelPitch
            tmpI = num.tile(tmpI, (npdiv**2, 1)) / subPixArea

            # import pdb;pdb.set_trace()
            # if corrected:
            #     # do conversion to tTh instead
            #     tempTTh, tmpEta = self.xyoToAng(tmpX+self.xc, tmpY+self.yc, units=self.pixelPitchUnit)
            #
            #     # can get rho in ideal frame
            #     tmpRho = self.workDist * num.tan(tempTTh)
            #     tmpEta = mapAngle(tmpEta, [startEta, 2*num.pi + startEta], units='radians')
            # else:
            # convert to polar
            #   - map eta into specified monotonic range
            tmpRho = num.sqrt( tmpX*tmpX + tmpY*tmpY )
            tmpEta = mapAngle(num.arctan2(tmpY, tmpX), [startEta, 2*num.pi + startEta], units='radians')

            inSector2 = ( (tmpRho >= startRho) & (tmpRho <= stopRho) ) \
                        & ( (tmpEta >= etaI) & (tmpEta <= etaF) )

            # import pdb;pdb.set_trace()
            tmpRho = tmpRho[inSector2]
            tmpI   = tmpI[inSector2]

            binId = num.floor( ( tmpRho - startRho ) / deltaRho )
            nSubpixelsIn = len(binId)

            # import pdb;pdb.set_trace()
            tmpI  = sparse.csc_matrix( \
                ( tmpI, (binId, num.arange(nSubpixelsIn)) ), shape=(numRho, nSubpixelsIn) )
            binId = sparse.csc_matrix( \
                ( num.ones(nSubpixelsIn), (binId, num.arange(nSubpixelsIn)) ), shape=(numRho, nSubpixelsIn) )

            # Normalized contribution to the ith sector's radial bins
            binIdSum = binId.sum(1)
            if num.any(binIdSum <= 0):
                import string
                raise RuntimeError, 'got binId sum of '+string.join(num.array(binIdSum).flatten().astype(str), ',')
            polImg['intensity'][i, :] = (tmpI.sum(1) / binIdSum).T
            
        return polImg


def mar165IDim(mode):
    if not isinstance(mode, int) or not [1,2,4,8].count(mode):
        raise RuntimeError, 'unknown mode : '+str(mode)
    idim = 4096 / mode
    return idim

class DetectorGeomMar165(Detector2DRC):
    __vfu = 0.2 # made up
    __vdk = 1800 # made up
    def __init__(self, *args, **kwargs):

        mode = 1
        if kwargs.has_key('mode'):
            mode = kwargs.pop('mode')
        readerClass = eval('ReadMar165NB%d' % (mode))
        idim = mar165IDim(mode)
        nrows = ncols = idim
        pixelPitch = 165.0 / idim # mm
        reader = readerClass()
        
        self.mode = mode
        
        Detector2DRC.__init__(self, 
                              nrows, ncols, pixelPitch, 
                              self.__vfu, self.__vdk,
                              reader,
                              *args, **kwargs)
        return

    def getDParamDflt(self):    
        return []
    def setDParamZero(self):    
        return
    def getDParamScalings(self):    
        return []
    def getDParamRefineDflt(self):
        return []
    #
    def radialDistortion(self, xin, yin, invert=False):
        'no distortion correction'
        return xin, yin

class DetectorGeomGE(Detector2DRC):
    """
    handle geometry of GE detector, such as geometric and radial distortion corrections;
    x and y are in pixels, as is rho;
    pixels are numbered from (0,0);
    """
    
    __vfu            = 0.2 # made up
    __vdk            = 1800 # made up
    # 200 x 200 micron pixels
    __pixelPitch     = 0.2      # in mm
    __idim           = ReadGE._ReadGE__idim
    __nrows          = ReadGE._ReadGE__nrows
    __ncols          = ReadGE._ReadGE__ncols
    __dParamDflt     = [   0.0,       0.0,       0.0,      2.0,      2.0,      2.0]
    __dParamZero     = [   0.0,       0.0,       0.0,      2.0,      2.0,      2.0]
    __dParamScalings = [   1.0,       1.0,       1.0,      1.0,      1.0,      1.0]
    __dParamRefineDflt = (True,      True,      True,    False,    False,    False)
    
    def __init__(self, *args, **kwargs):
        
        kwargs.setdefault('getDParamDflt',       self.getDParamDflt)
        kwargs.setdefault('setDParamZero',       self.setDParamZero)
        kwargs.setdefault('getDParamScalings',   self.getDParamScalings)
        kwargs.setdefault('getDParamRefineDflt', self.getDParamRefineDflt)
        kwargs.setdefault('radialDistortion',    self.radialDistortion)

        reader = kwargs.pop('reader', None)
        if reader is None:
            readerKWArgs = kwargs.pop('readerKWArgs', {})
            reader = ReadGE(None, **readerKWArgs)
        
        Detector2DRC.__init__(self, 
                              self.__nrows, self.__ncols, self.__pixelPitch, 
                              self.__vfu, self.__vdk,
                              reader,
                              *args, **kwargs)
        return

    def getNewReader(self, filename, *args, **kwargs):
        # newR = self.reader.__class__(filename, self.ncols, self.nrows, *args, **kwargs)
        newR = ReadGE(filename, *args, **kwargs)
        return newR
    def makeNew(self, *args, **kwargs):
        kwargs.setdefault('reader',self.reader)
        newDG = self.__class__(self, *args, **kwargs)
        return newDG

    def getDParamDflt(self):
        return self.__dParamDflt
    def setDParamZero(self):
        self.dparm = self.__dParamZero
        return
    def getDParamScalings(self):
        return self.__dParamScalings
    def getDParamRefineDflt(self):
        return self.__dParamRefineDflt
    def radialDistortion(self, xin, yin, invert=False):
        """
        Apply radial distortion to polar coordinates on GE detector

        xin, yin are 1D arrays or scalars, assumed to be relative to self.xc, self.yc
        Units are [mm, radians].  This is the power-law based function of Bernier.

        Available Keyword Arguments :

        invert = True or >False< :: apply inverse warping
        """
        if self.dparms[0] == 0 and self.dparms[1] == 0 and self.dparms[2] == 0:
            xout = xin
            yout = yin
        else:
            # canonical max radius based on perfectly centered beam
            #   - 204.8 in mm or 1024 in pixel indices
            rhoMax = self.__idim * self.__pixelPitch / 2
            
            # make points relative to detector center
            x0 = (xin + self.xc) - rhoMax
            y0 = (yin + self.yc) - rhoMax
            
            # detector relative polar coordinates
            #   - this is the radius that gets rescaled
            rho0 = num.sqrt( x0*x0 + y0*y0 )
            eta0 = num.arctan2( y0, x0 )
                    
            if invert:
                # in here must do nonlinear solve for distortion
                # must loop to call fsolve individually for each point
                rho0   = num.atleast_1d(rho0)
                rShape = rho0.shape
                rho0   = num.atleast_1d(rho0).flatten()
                rhoOut = num.zeros(len(rho0), dtype=float)
                
                eta0   = num.atleast_1d(eta0).flatten()
                
                rhoSclFuncInv = lambda ri, ni, ro, rx, p: \
                    (p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) + \
                     p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) + \
                     p[2]*(ri/rx)**p[5] + 1)*ri - ro 
                
                rhoSclFIprime = lambda ri, ni, ro, rx, p: \
                    p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) * (p[3] + 1) + \
                    p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) * (p[4] + 1) + \
                    p[2]*(ri/rx)**p[5] * (p[5] + 1) + 1
                
                for iRho in range(len(rho0)):
                    rhoOut[iRho] = fsolve(rhoSclFuncInv, rho0[iRho], 
                                          fprime=rhoSclFIprime, 
                                          args=(eta0[iRho], rho0[iRho], rhoMax, self.dparms) )
                    pass
                
                rhoOut = rhoOut.reshape(rShape)            
            else:
                # usual case: calculate scaling to take you from image to detector plane
                # 1 + p[0]*(ri/rx)**p[2] * num.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
                rhoSclFunc = lambda ri, rx=rhoMax, p=self.dparms, ni=eta0: \
                             p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) + \
                             p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) + \
                             p[2]*(ri/rx)**p[5] + 1
                
                rhoOut = num.squeeze( rho0 * rhoSclFunc(rho0) )
                pass
            
            xout = rhoOut * num.cos(eta0) + rhoMax - self.xc
            yout = rhoOut * num.sin(eta0) + rhoMax - self.yc
        
        return xout, yout

class DetectorGeomFrelon(Detector2DRC):
    """
    handle geometry of GE detector, such as geometric and radial distortion corrections;
    x and y are in pixels, as is rho;
    pixels are numbered from (0,0);
    """
    
    # 50 X 50 micron pixels
    __pixelPitch     = 0.05      # in mm
    __idim           = ReadGE._ReadGE__idim
    __nrows          = ReadGE._ReadGE__nrows
    __ncols          = ReadGE._ReadGE__ncols
    __dParamDflt     = [   0.0,      0.0,     0.0,      2.0,      2.0,      2.0]
    __dParamZero     = [   0.0,      0.0,     0.0,      2.0,      2.0,      2.0]
    __dParamScalings = [   1.0,      1.0,     1.0,      1.0,      1.0,      1.0]
    __dParamRefineDflt = (True,     True,    True,    False,    False,    False)
    
    def __init__(self, *args, **kwargs):
        
        Detector2DRC.__init__(self, 
                              self.__nrows, self.__ncols, self.__pixelPitch,
                              ReadGE,
                              *args, **kwargs)
        return
    
    def getDParamDflt(self):
        return self.__dParamDflt
    def setDParamZero(self):
        self.dparm = self.__dParamZero
        return 
    def getDParamScalings(self):    
        return self.__dParamScalings
    def getDParamRefineDflt(self):
        return self.__dParamRefineDflt
    def radialDistortion(self, xin, yin, invert=False):
        """    
        Apply radial distortion to polar coordinates on GE detector
        
        xin, yin are 1D arrays or scalars, assumed to be relative to self.xc, self.yc
        Units are [mm, radians].  This is the power-law based function of Bernier.
        
        (p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) + \
         p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) + \
         p[2]*(ri/rx)**p[5] + 1)*ri
         
         1 + \
         p[0]*(ri/rx)**p[2] * num.cos(p[4] * ni) + \
         p[1]*(ri/rx)**p[3]
         
        Available Keyword Arguments :
        
        invert = True or >False< :: apply inverse warping
        """
        if self.dparms[0] == 0 and self.dparms[1] == 0 and self.dparms[2] == 0:
            xout = xin
            yout = yin
        else:
            # canonical max radius based on perfectly centered beam
            #   - 204.8 in mm or 1024 in pixel indices
            rhoMax = self.__idim * self.__pixelPitch / 2
            
            # make points relative to detector center
            x0 = (xin + self.xc) - rhoMax
            y0 = (yin + self.yc) - rhoMax
            
            # detector relative polar coordinates
            #   - this is the radius that gets rescaled
            rho0 = num.sqrt( x0*x0 + y0*y0 )
            eta0 = num.arctan2( y0, x0 )
                    
            if invert:
                # --> in here must do nonlinear solve for distortion
                # must loop to call fsolve individually for each point
                rho0   = num.atleast_1d(rho0)
                rShape = rho0.shape
                rho0   = num.atleast_1d(rho0).flatten()
                rhoOut = num.zeros(len(rho0), dtype=float)
                
                eta0   = num.atleast_1d(eta0).flatten()
                
                rhoSclFunc = lambda ri, ni, ro, p=self.dparms, rx=rhoMax: \
                    (p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) + \
                     p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) + \
                     p[2]*(ri/rx)**p[5] + 1)*ri - ro 
                
                for iRho in range(len(rho0)):
                    rhoOut[iRho] = fsolve(rhoSclFunc, rho0[iRho], args=(eta0[iRho], rho0[iRho]))
                    pass
                rhoOut = rhoOut.reshape(rShape)            
            else:
                # usual case: calculate scaling to take you from image to detector plane
                # 1 + p[0]*(ri/rx)**p[2] * num.cos(p[4] * ni) + p[1]*(ri/rx)**p[3]
                rhoSclFunc = lambda ri, p=self.dparms, rx=rhoMax, ni=eta0: \
                             p[0]*(ri/rx)**p[3] * num.cos(2.0 * ni) + \
                             p[1]*(ri/rx)**p[4] * num.cos(4.0 * ni) + \
                             p[2]*(ri/rx)**p[5] + 1
                
                rhoOut = num.squeeze( rho0 * rhoSclFunc(rho0) )
                
            xout = rhoOut * num.cos(eta0) + rhoMax - self.xc
            yout = rhoOut * num.sin(eta0) + rhoMax - self.yc
        return xout, yout

class DetectorGeomQuadGE(DetectorBase):
    """
    No global parameters -- all detector parameters hang off of the sub-detectors;
    although some data are stored off of this higher level class for convenience
    """

    __inParmDict = {
        'quadAngle'   : -53.0*(num.pi/180.), # 53 degrees, converted to radians
        'quadPad'     : 20., # in-plane distance between active surfaces of detectors; in pixelPitchUnit
        'quadShift'   : 100., # in-plane shift of detectors to allow for central opening; in pixelPitchUnit
        'quadOffsets' : 'hydra',
        }

    # order by diagram numbering
    __quadHydraLQuad    = [0,3,1,2] # draw ones in front first
    __quadHydraPush     = [0.,98.,0.,98.] # offset in working distance to overlap detector frames; in pixelPitchUnit
    __quadHydraIOffsets = [ [0,0], [1,0], [0,-1], [1,-1] ]
    __quadHydraIPad     = [ [0,0], [1,0], [0,-1], [1,-1] ]
    __quadHydraIXShift  = [ 0,  0,  1,  1 ]
    __quadHydraIYShift  = [ 0,  1,  0,  1 ]

    def __init__(self, *args, **kwargs):
        
        'pass ReadGE instance as the reader for now; perhaps make a ReadQuadGE class later if it turns out to be needed'
        reader = ReadGE(None)
        DetectorBase.__init__(self, reader)
        
        # defaults and kwargs parsing        
        for parm, val in self.__inParmDict.iteritems():
            if kwargs.has_key(parm):
                val = kwargs.pop(parm)
            self.__setattr__(parm, val)

        dgDummy = DetectorGeomGE()

        'cleanup after auto-parsing of keyword args'
        self.quadAngle = valunits.toFloat(self.quadAngle, 'radians')
        self.quadPad   = valunits.toFloat(self.quadPad,   dgDummy.pixelPitchUnit)
        self.quadShift = valunits.toFloat(self.quadShift, dgDummy.pixelPitchUnit)
        # self.quadPush  = valunits.toFloat(self.quadPush,  dgDummy.pixelPitchUnit)

        iRefQuad = 0

        if len(args) == 1:
            if hasattr(args[0], 'xyoToAng'):
                dgIn = args[0]
                dRef = dgIn.makeNew(**kwargs)
                dRef.zTilt = self.quadAngle
            else:
                raise RuntimeError, 'do not know how to parse arguments'
            self.detectors = [dRef]
            self.detectors.append(dRef.makeNew())
            self.detectors.append(dRef.makeNew())
            self.detectors.append(dRef.makeNew())
        elif len(args) == 4:
            'assume these whould be used more or less as they are, so do not set zTilt'
            self.detectors = []
            for iQuad in range(4): # not self.lQuad, because doing append below!
                dgIn = args[iQuad]
                assert hasattr(dgIn, 'xyoToAng'),\
                    'argument does not appear to be a detector geometry'
                self.detectors.append(dgIn.makeNew(**kwargs))
        else:
            raise RuntimeError, 'do not know how to parse arguments'
        dRef = self.detectors[iRefQuad]

        if self.quadOffsets == 'hydra':
            'do not do quadAngle rotation in centers, as rotation happens _after_ translation'
            self.lQuad = self.__quadHydraLQuad
            # assert len(self.quadPush) == 4,\
            #     'quadPush wrong length, should be 4 and is : '+str(self.quadPush)
            translations = [None, None, None, None]
            for iQuad in self.lQuad:
                x0 = dRef.extent[0] * self.__quadHydraIOffsets[iQuad][0] + \
                    self.quadPad * self.__quadHydraIPad[iQuad][0] + \
                    self.quadShift * self.__quadHydraIXShift[iQuad]
                y0 = dRef.extent[1] * self.__quadHydraIOffsets[iQuad][1] + \
                    self.quadPad * self.__quadHydraIPad[iQuad][1] + \
                    self.quadShift * self.__quadHydraIYShift[iQuad]
                translations[iQuad] = [
                        x0, # x0 * num.cos(self.quadAngle) - y0 * num.sin(self.quadAngle),
                        y0, # x0 * num.sin(self.quadAngle) + y0 * num.cos(self.quadAngle),
                        self.__quadHydraPush[iQuad],
                        ]
            self.quadOffsets = num.array(translations)
            dRef.xc =  dRef.extent[0]*1.0 + (self.quadShift - self.quadPad)*0.5 + self.quadPad
            dRef.yc =  dRef.extent[1]*0.0 + (self.quadShift - self.quadPad)*0.5
        else:
            'might want to order based on workDist offsets'
            self.lQuad = [0,1,2,3]
            self.quadOffsets = num.atleast_2d(self.quadOffsets)
            assert self.quadOffsets.shape == (4,3), \
                'quadOffsets wrong shape, should be (4,3) and is : '+str(self.quadOffsets)

        self.setCentersFromRef(iRefQuad=iRefQuad)

        return

    @classmethod
    def getRefineFlagsDflt(cls):
        '''
        no parameters to refine for this detector;
        call fitProcedureA for a procedure to refine individual sub-detectors
        '''
        retval = []
        return retval

    def getTThMax(self):
        '''
        min over sub-detectors, where for each sub-detector max two-theta is evaluated as
        the max over points checked in getTThMax for the sub-detector
        '''
        tThMaxList = [ self.detectors[iQuad].getTThMax(func=num.max) for iQuad in self.lQuad ]
        tThMax = num.min(num.array(tThMaxList))
        return tThMax

    def setCentersFromRef(self, iRefQuad=0):
        'this assumes all of the tilts are the same'
        lQuad = copy.copy(self.lQuad)
        lQuad.remove(iRefQuad)
        dRef = self.detectors[iRefQuad]
        for iQuad in lQuad:
            dg = self.detectors[iQuad]
            dg.xc       = dRef.xc       + ( self.quadOffsets[iRefQuad][0] - self.quadOffsets[iQuad][0] )
            dg.yc       = dRef.yc       + ( self.quadOffsets[iRefQuad][1] - self.quadOffsets[iQuad][1] )
            dg.workDist = dRef.workDist + ( self.quadOffsets[iRefQuad][2] - self.quadOffsets[iQuad][2] )
        return

    def setQuadOffsets(self, iRefQuad=0):
        'this assumes all of the tilts are the same'
        lQuad = copy.copy(self.lQuad)
        lQuad.remove(iRefQuad)
        dRef = self.detectors[iRefQuad]
        self.quadOffsets[iRefQuad][:] = num.zeros(3)
        for iQuad in lQuad:
            dg = self.detectors[iQuad]
            self.quadOffsets[iQuad][0] = self.quadOffsets[iRefQuad][0]  + dRef.xc        - dg.xc
            self.quadOffsets[iQuad][1] = self.quadOffsets[iRefQuad][1]  + dRef.yc        - dg.yc
            self.quadOffsets[iQuad][2] = self.quadOffsets[iRefQuad][2]  + dRef.workDist  - dg.workDist
        return

    def drawRings(self, drawOn, planeData, workDist=None, **kwargs):
        """
        assumes ideal geometry, as from displayIdeal
        """
        workDist = workDist or self.detectors[0].workDist
        dg = self.detectors[0]
        dg.drawRings(drawOn, planeData, ideal=workDist, **kwargs)
        return
    def displaySub(self, iQuad, thisframe, planeData=None, **kwargs):
        """
        convenience for displaying a sub-detector
        ...*** need to code labAxesGlyph support in display
        """
        dg = self.detectors[iQuad]
        retval = dg.display(thisframe, planeData=planeData, labAxesGlyph=True, **kwargs)
        return retval

    def displayIdeal(self, framesQuad, planeData=None, workDist=None, nlump=None,
                     doFmtCoord=True, **kwargs):
        """
        display all sub-detectors on an idealized detector plane

        If matplotlib gets around to enabling the transform argument
        to imshow, that might be a much faster approach than what is
        currently done here, although what is done here is nice in
        that it takes account of all of the distortions, not just the
        in-plane rotation. The idea would be that the in-plane
        rotation would be, by far, the biggest effect.  import
        matplotlib.transforms as mtransforms
        tr = mtransforms.Affine2D()
        tr.rotate(self.zTilt)
        imshow( , transform=tr)
        """
        pwKWArgs = plotwrap.PlotWrap.popKeyArgs(kwargs)
        pw = plotwrap.PlotWrap(**pwKWArgs)
        pw.a.set_axis_bgcolor('white')
        pw.a.set_aspect('equal')
        retval = pw

        'unless specified, use working distance of 0th sub-detector'
        workDist = workDist or self.detectors[0].workDist

        nlump = nlump or 4

        cmap = None; vmin = None; vmax = None; norm = None;

        if doFmtCoord:
            fmtCoord = FmtCoordIdeal(planeData, workDist)
        for iQuad in self.lQuad:
            dg = self.detectors[iQuad]
            thisframe = framesQuad[iQuad]

            h, xedges, yedges, mask = dg.renderIdeal(thisframe, workDist=workDist, nlump=nlump)

            if cmap is None:
                vmin, vmax, cmap = dg.reader.getDisplayArgs(h, kwargs)
                cmap.set_under(color='white', alpha=0.0)
                norm = matplotlib.colors.Normalize(clip=False, vmin=vmin, vmax=vmax)
            if True:
                '''
                hack around to make for better plotting;
                needed due to tolerancing in cmap for when a color is under?
                '''
                h[h < 0] =  vmin + (vmax-vmin)*0.001
                h[mask]  = -vmax

            pw.a.pcolor(xedges, yedges, h.T, cmap=cmap, norm=norm, **kwargs)

            if doFmtCoord:
                iDetector = iQuad+1 # people like things numbered from 1
                fmtCoord.addDetectorData(iDetector, xedges, yedges, h, mask)

            del thisframe

        if doFmtCoord:
            pw.a.format_coord = fmtCoord
        pw.a.set_autoscale_on(False)
        pw.show()
        return retval
    
    def fitProcedureA(self, planeData, framesQuad, iRefQuad=0, 
                      funcType=funcTypeDflt, funcXVecList = None, quadr=1, 
                      doGUI=0, doMRingPlot=False):
        """
        A procedure for fitting the set of detectors;
        do not need to click 'fit' button in GUI -- done inside the procedure.

        Watch out -- MultiRingEval instances a memory hogs, especially while creating Jacobian matrices!

        If want to just refine detector geometry and not the functional forms for the rings,
        pass funcXVecList as True or as something like a list of arrays from MultiRingEval.getFuncXVecList()
        """
        
        assert len(framesQuad) == 4,\
            'need len(framesQuad) to be 4'

        if funcXVecList:
            if hasattr(funcXVecList,'__len__'):
                funcXVecListList = funcXVecList
            else:
                funcXVecListList = [funcXVecList for iQuad in self.lQuad]
        else:
                funcXVecListList = [None for iQuad in self.lQuad]
        assert len(funcXVecListList) == len(self.lQuad), \
            'funcXVecListList wrong length'

        lQuad = copy.copy(self.lQuad)
        lQuad.remove(iRefQuad)
        #dRef = self.detectors[iRefQuad]
        # refPos = num.array([dRef.xc, dRef.yc, dRef.workDist])
        'set quadOffsets from current centers, gets used in setCentersFromRef below'
        self.setQuadOffsets(iRefQuad)

        iQuad = iRefQuad
        dg = self.detectors[iQuad]
        frame = dataToFrame(framesQuad[iQuad])
        if doGUI:
            'offer GUI to let the user twiddle parameters'
            displayKWArgs={'winTitle':'detector %d'%(iQuad+1)}
            pw    = dg.drawRingsGUI(frame, planeData, displayKWArgs=displayKWArgs)
            raw_input('Enter to continue (will destroy GUI and proceed with fit)')
            pw.destroy()
            del pw
        #
        'now go ahead and do fit'
        mRing  = MultiRingEval(dg, planeData, frame, funcType=funcType,
                               funcXVecList=funcXVecListList[iQuad])
        xFit   = mRing.doFit()
        if doGUI > 1 or doMingPlot:
            mRing(xFit, makePlots=True)
            if doGUI > 1:
                raw_input('Enter to continue (close mRing windows manually)')
        # params = dg.getParams(allParams=True)
        del mRing, frame
        #
        'update centers and working distances of others to correspond to results from fit of iRefQuad'
        self.setCentersFromRef(iRefQuad=iRefQuad)
        #
        'now fit everyone else'
        for iQuad in lQuad:
            dg = self.detectors[iQuad]
            frame = dataToFrame(framesQuad[iQuad])
            if doGUI > 1:
                'offer GUI to let the user twiddle parameters'
                displayKWArgs={'winTitle':'detector %d'%(iQuad+1)}
                pw    = dg.drawRingsGUI(frame, planeData, displayKWArgs=displayKWArgs)
                raw_input('Enter to continue (will destroy GUI and proceed with fit)')
                pw.destroy()
                del pw
            #
            'now go ahead and do fit'
            mRing  = MultiRingEval(dg, planeData, frame, funcType=funcType,
                                   funcXVecList=funcXVecListList[iQuad])
            xFit   = mRing.doFit()
            if doGUI > 1 or doMRingPlot:
                mRing(xFit, makePlots=True)
                if doGUI > 1:
                    raw_input('Enter to continue (close mRing windows manually)')
            # params = dg.getParams(allParams=True)
            del mRing, frame

        self.setQuadOffsets(iRefQuad)
        return

def getOmegaMMReaderList(readerList, overall=False):
    """
    get omega min/max information from a list of readers
    """
    retval = []
    for reader in num.atleast_1d(readerList):
        omegaMin, omegaMax = reader.getOmegaMinMax()
        retval.append((omegaMin,omegaMax))
    if overall:
        retval = (min(zip(*retval)[0]), max(zip(*retval)[1]))
    return retval

# ============================== Utility functions for instantiating detectors
#
def detectorList():
    return ["ge", "mar165", "generic"]

def newDetector(detectorType, *args, **kwargs):
    """Return a detector of the requested type

    INPUTS

    detectorType - a string in the detector type list [see detectorList()]
    """
    if len(args) == 0:
        if kwargs.has_key('gParms'):
            args = kwargs.pop('gParms')
        if kwargs.has_key('dParms'):
            dp = kwargs.pop('dParms')
            if dp:
                kwargs['distortionParams'] = dp
                pass
            pass
        pass
    
    dt = detectorType.lower()
    if dt == "ge":
        d = DetectorGeomGE(*args, **kwargs)
    elif dt == "mar165":
        d = DetectorGeomMar165(*args, **kwargs)
    elif dt == "generic":
        ncols = kwargs.pop('ncols')
        nrows = kwargs.pop('nrows')
        pixelPitch = kwargs.pop('pixelPitch')
        d = newGenericDetector(ncols, nrows, pixelPitch, *args, **kwargs)
    else:
        #d = None
        emsg = 'No such detector:  "%s"' % dt
        raise ValueError(emsg)

        pass


    return d

def newGenericReader(ncols, nrows, *args, **kwargs):
    '''
    currently just returns a Framer2DRC
    '''
    
    # retval = Framer2DRC(ncols, nrows, **kwargs)
    filename = kwargs.pop('filename', None)
    retval = ReadGeneric(filename, ncols, nrows, *args, **kwargs)
    
    return retval

def newGenericDetector(ncols, nrows, pixelPitch, *args, **kwargs):
    """
    If reader is passed as None, then a generic reader is created

    Keyword Arguments:
	vFactorUnc
	vDark
	reader
	readerKWArgs
	getDParamDflt
	setDParamZero
	getDParamScalings
	getDParamRefineDflt
	radialDistortion
        
    If *args is an existing detector geometry, then
    additional keyword arguments may include:
	pVec
    
    If *args is (xc, yc, workDist, xTilt, yTilt, zTilt) detector parameters, then 
    additional keyword arguments may include:
	distortionParams

    """
    
    vFactorUnc   = kwargs.pop('vFactorUnc',0.2)
    vDark        = kwargs.pop('vDark',1800)
    reader       = kwargs.pop('reader',None)
    readerKWArgs = kwargs.pop('readerKWArgs',{})

    'default functions corresponding to no distortion'
    def getDParamDflt_dflt():    
        return []
    def setDParamZero_dflt():    
        return
    def getDParamScalings_dflt():    
        return []
    def getDParamRefineDflt_dflt():
        return []
    def radialDistortion_dflt(xin, yin, invert=False):
        'no distortion correction'
        return xin, yin
    
    getDParamDflt        = kwargs.setdefault('getDParamDflt', getDParamDflt_dflt)
    setDParamZero        = kwargs.setdefault('setDParamZero', setDParamZero_dflt)
    getDParamScalings    = kwargs.setdefault('getDParamScalings', getDParamScalings_dflt)
    getDParamRefineDflt  = kwargs.setdefault('getDParamRefineDflt', getDParamRefineDflt_dflt)
    radialDistortion     = kwargs.setdefault('radialDistortion', radialDistortion_dflt)

    if reader is None:
        reader = newGenericReader(ncols, nrows, **readerKWArgs)

    detector = Detector2DRC(ncols, nrows, pixelPitch, 
                            vFactorUnc, vDark,
                            reader, *args, **kwargs)

    return detector
