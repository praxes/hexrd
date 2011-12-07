#!/usr/bin/env python
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

from femODFUtil import pfig as pfigPkg

import xrdBase
if xrdBase.haveMultiProc:
    from xrdBase import multiprocessing

polePoints_MP = None
nPerChunk_MP = None
mRodr_MP = None
applySym_MP = None
cV_MP = None
#
def polePathIntegThis(iLo):
    from femODFUtil.mRodrUtil import convertFemODFCoo
    
    global polePoints_MP
    global nPerChunk_MP
    global mRodr_MP
    global applySym_MP
    global cV_MP
    polePoints = polePoints_MP
    nPerChunk  = nPerChunk_MP 
    mRodr      = mRodr_MP     
    applySym   = applySym_MP
    cV         = cV_MP        
    
    iHi = min(iLo+nPerChunk, polePoints.shape[1])
    pp = polePoints[:,iLo:iHi]

    'getPoleProjection is where the work sits'
    odfPfProj = mRodr.getPoleProjection(cV, pp, applySym)
    
    retval = convertFemODFCoo(odfPfProj)
    return retval
#
def polePathInteg(
    cV,
    polePoints, 
    mRodr, 
    applySym=True,
    nPerChunk=None,
    nCPUs=None,
    doMultiProc=True,
    verbose=False,
    ):
    from scipy import sparse
    
    multiProcMode = xrdBase.haveMultiProc and doMultiProc
    
    if multiProcMode:
        from xrdBase import multiprocessing
        nCPUs = nCPUs or xrdBase.dfltNCPU
    else:
        nCPUs = 1

    nPerChunk = nPerChunk or max(int(polePoints.shape[1]/(10*nCPUs)),1)
    iLoList = range(0, polePoints.shape[1], nPerChunk)
    if verbose: 
        print 'nCPUs, nPerChunk : %g %g' % (nCPUs, nPerChunk)
    
    global polePoints_MP
    global nPerChunk_MP
    global mRodr_MP
    global applySym_MP
    global cV_MP
    polePoints_MP = polePoints
    nPerChunk_MP  = nPerChunk
    mRodr_MP      = mRodr
    applySym_MP   = applySym
    cV_MP         = cV

    if multiProcMode:
        pool = multiprocessing.Pool(nCPUs)
        results = pool.map(polePathIntegThis, iLoList, chunksize=1)
    else:
        results = map(polePathIntegThis, iLoList)
    retval = (sparse.vstack(results)).tocsr()
    
    polePoints_MP = None
    nPerChunk_MP = None
    mRodr_MP = None
    applySym_MP = None
    cV_MP = None
    if multiProcMode:
        pool.close()
    
    return retval
