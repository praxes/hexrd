#! /usr/bin/env python
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
import mdef

import numpy as num

try:
    import multiprocessing
    haveMultiProc = True
    dfltNCPU = int(multiprocessing.cpu_count()/2)
    if dfltNCPU > 3: 
        dfltNCPU = dfltNCPU/2
except:
    haveMultiProc = False

def dataToFrame(data, sumImg=True):
    """
    utility function to allow flexibility in input
    
    data can be:
    (*) an instance of ReadGE or the like, which is already set up, in which
        case all frames are used and flattened
    (*) a frame
    """
    if hasattr(data, 'getNFrames'):
        reader = data.makeNew()
        nFrames = reader.getNFrames()
        frame = reader.read(nframes=nFrames, sumImg=sumImg)
    elif hasattr(data, 'shape') and len(data.shape) == 2:
        'assume the data is a frame'
        frame = data
    else:
        raise RuntimeError, 'do not know what to do with data : '+str(type(data))
    return frame

def getGaussNDParams(xList, w=None, v=None):
    from math import sqrt, log
    
    nDim = len(xList)
    
    xVec = num.empty(nDim+nDim+2)
    
    if w is not None:
         assert xList[0].shape == w.shape,\
             'w not the same shape as other arguments; %s != %s' % (str(x.shape), str(w.shape))
    
    if v is None:
        bg   = 0.0e0
        vNbg = num.ones(x.shape)
    else:
        bg   = num.min(v)
        vNbg = v - bg
    if w is None:
        vwNbg = vNbg
    else:
        if len(w.shape) == 2:
            nQP = w.shape[1]
            vwNbg = num.tile(vNbg, (nQP,1)).T * w
        elif len(w.shape) == 1:
            vwNbg = vNbg * w
        else:
            raise NotImplementedError,\
                'shape of w has length %d' % (len(shape(w)))
    vwSum = float(num.sum(vwNbg))
    
    if vwSum <= 0:
        'just swag something'
        # raise RuntimeError, 'vwSum <= 0'
        vwNbg = num.ones_like( vwNbg )
        vwSum = float(num.sum(vwNbg))

    com = []
    for iDim, xThis in enumerate(xList):
        mu = num.sum(xThis * vwNbg) / vwSum
        com.append(mu)
    #
    for iDim, xThis in enumerate(xList):
        diffs = xThis - com[iDim]
        sigma = sqrt(num.sum(vwNbg*diffs*diffs)/vwSum)
        xVec[iDim]      = com[iDim]
        xVec[iDim+nDim] = sqrt(8.0 * log(2.0)) * sigma # FWHM
    xVec[nDim+nDim]   = vwNbg.max()
    xVec[nDim+nDim+1] = bg
    
    return xVec

