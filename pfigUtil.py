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
import arrayUtil
from arrayUtil import getMem
from arrayUtil import num # import numpy as num

def sph2n(coords_sph):
    '''convert from chi/eta spherical coordinates to normal vectors;
    can use with coords from femODF.FemHemisphere
    '''
    
    assert coords_sph.shape[0] == 2, 'shape[0] not 2'
    
    z   = num.cos(coords_sph[0,:])
    rho = num.sin(coords_sph[0,:])
    x = rho * num.cos(coords_sph[1,:])
    y = rho * num.sin(coords_sph[1,:])

    # nVectors = arrayUtil.toArray(num.vstack((x,y,z)))
    nVectors = getMem((3,coords_sph.shape[1]))
    nVectors[0,:] = x
    nVectors[1,:] = y
    nVectors[2,:] = z
    
    return nVectors

def n2sph(nVectors):
    from math import atan2, sqrt
    assert nVectors.shape[0] == 3, 'shape[0] not 3'
    sph = arrayUtil.getMem((2,nVectors.shape[1]))
    for iVector in range(nVectors.shape[1]):
        x = nVectors[0,iVector]
        y = nVectors[1,iVector]
        z = nVectors[2,iVector]
        r = sqrt(x*x+y*y)
        sph[0,iVector] = atan2(r,z)
        sph[1,iVector] = atan2(y,x)
    return sph

def n2eap(nVectors, flip=True):
    """
    unit vectors to equal-area projection
    """
    import copy
    nVecs = copy.deepcopy(nVectors)
    nPnts = nVectors.shape[1]
    retval = getMem((2,nPnts))
    
    belowEquator = num.where(nVecs[2,:] < -1e-4)
    if flip:
        nVecs[:,belowEquator] = -nVecs[:,belowEquator]
    else:
        'put on the equator'
        nVecs[2,belowEquator] = 0.
        # num.apply_along_axis(num.linalg.norm, 0, nVecs[0:2,belowEquator])
        norms = num.sqrt( nVecs[0,belowEquator]*nVecs[0,belowEquator] + nVecs[1,belowEquator]*nVecs[1,belowEquator] )
        nVecs[0:2,belowEquator] = nVecs[0:2,belowEquator] / norms
    
    r2    = nVecs[0,:]*nVecs[0,:] + nVecs[1,:]*nVecs[1,:]
    r2pos = num.where(r2 > 0.)
    
    n31 = 1.0 - num.abs(nVecs[2,:])
    
    den        = num.ones(nPnts)
    den[r2pos] = num.sqrt(2.0*r2[r2pos])
    
    dist_np = num.zeros(nPnts)
    dist_np = num.sqrt(r2 + n31*n31) / den
    
    retval = num.vstack((nVecs[0,:] * dist_np, nVecs[1,:] * dist_np)) # .T
    
    return retval

def renderEAProj(nVecs, vals, n, patch=False, sum=False, nzByContrib=True, northernOnly=False):
    """
    render an equal-area projects of general pole values, on an n-by-n grid;
    if sum=True, then sum contributions, otherwise average them;
    returns a masked array
    """
    
    if northernOnly:
        bNVecsN = nVecs[2,:] >= 0
        nVecs = nVecs[:,bNVecsN]
        vals  = vals.flatten()[bNVecsN]
        eaProj = n2eap(nVecs, flip=False)
    else:
        vals  = vals.flatten()
        eaProj = n2eap(nVecs, flip=True)
    
    eaProjX, eaProjY = eaProj[0,:], eaProj[1,:]
    if patch:
        xyd = 0.5 * max(eaProjX.max() - eaProjX.min(), eaProjY.max() - eaProjY.min())
        xMid = 0.5 * (eaProjX.max() + eaProjX.min())
        yMid = 0.5 * (eaProjY.max() + eaProjY.min())
        rng = ((xMid-xyd,xMid+xyd),(yMid-xyd,yMid+xyd))
    else:
        rng = ((-1.,1.),(-1.,1.))
    H, xedges, yedges = num.histogram2d(
        eaProjX, eaProjY,
        weights=vals,
        bins=n, range=rng)
    nH, xedges, yedges = num.histogram2d(
        eaProjX, eaProjY, 
        bins=n, range=rng)
    nonz = nH > 0
    if not sum:
        H[nonz] = H[nonz] / nH[nonz]
    if not nzByContrib:
        nonz = H > 0
    retval = num.ma.masked_array(H, mask=-nonz, fill_value=0., hard_mask=False)
    return retval

def fromSouthern(nVecs, invert):
    if invert:
        'invert through origin'
        retval = -nVecs
    else:
        'rotate about vertical axis in plane of pole figure'
        retval = num.vstack((-nVecs[0,:], nVecs[1,:], -nVecs[2,:]))
    return retval

def drawLines(pw, pointLists=[], 
              netStyle = None, netNDiv = 12, netAlpha = 0.5, rMat=None, 
              southern=False, invertFromSouthern=True,
              origin=(0.,0.), r=1.0):
    import matrixUtils
    import orientations as ors
    
    x0, y0 = origin
    
    lines = pw.a.get_lines()
    for line in lines:
        line.remove()
    
    ringAngs = num.linspace(0,num.pi*2.,181)
    if netStyle:
        arcAngsM  = num.linspace(-num.pi/2.,num.pi/2.,91)
        arcAngsL  = num.linspace(0, num.pi, 181)
        rMat_yDiv = ors.RotInv(-num.pi/netNDiv, [0,1,0]).toMatrix()
    pw(r*num.cos(ringAngs)+x0,r*num.sin(ringAngs)+y0,style='k-')
    
    if netStyle:
        nVMerid = num.vstack( ( num.cos(arcAngsM), num.sin(arcAngsM), num.zeros(len(arcAngsM)) ) )
        for iDiv in range(netNDiv-1):
            nVMerid = num.dot(rMat_yDiv, nVMerid)
            eaProj = n2eap(nVMerid, flip=False)
            pw(r*eaProj[0,:]+x0, r*eaProj[1,:]+y0, alpha=netAlpha, style=netStyle)
        for latAng in num.pi/netNDiv*num.arange(1,netNDiv):
            polZ = num.cos(latAng)
            polR = num.sin(latAng)
            nVLat = num.vstack( ( polR*num.cos(arcAngsL), polZ*num.ones(len(arcAngsL)), polR*num.sin(arcAngsL)) )
            eaProj = n2eap(nVLat, flip=False)
            pw(r*eaProj[0,:]+x0, r*eaProj[1,:]+y0, alpha=netAlpha, style=netStyle)
    'done with drawing net'

    for points, pwKWArgs in pointLists:
        nVecs  = matrixUtils.unitVector(points)
        if rMat is not None:
            'rotate as did elsewhere'
            nVecs = num.dot(rMat, nVecs)
        bNVecsS = nVecs[2,:] < 0
        if southern:
            nVecsS = nVecs[:,bNVecsS]
            nVecsS = fromSouthern(nVecsS, invertFromSouthern)
            eaProj = n2eap(nVecsS, flip=False)
            pw(r*eaProj[0,:]+x0, r*eaProj[1,:]+y0, **pwKWArgs)
        else:
            nVecsN = nVecs[:,num.logical_not(bNVecsS)]
            eaProj = n2eap(nVecsN, flip=False)
            pw(r*eaProj[0,:]+x0, r*eaProj[1,:]+x0, **pwKWArgs)
    'done with pointLists'
    
    return
    
