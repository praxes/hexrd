#!/usr/bin/env python
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
