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

import arrayUtil
from arrayUtil import getMem
from arrayUtil import num
import pfigUtil

from femODFUtil.pfig import makePfigDict as makePfigDict_pfig
from femODFUtil.pfig import calcPfigProjOps as calcPfigProjOps_pfig
from femODFUtil.pfig import calcPfigProj as calcPfigProj_pfig
from femODFUtil.pfig import addH1ToSystem, buildRHS, loadSysBlocks, storeSysBlocks

def makePfigDict(hkl, **kwargs):
    d = makePfigDict_pfig(hkl, **kwargs)
    cV = d['crystalVector']
    n  = cV / num.linalg.norm(cV)
    d['nnProj'] = num.outer(n,n)
    return d

def getEIScaling(projOpLists, pfigsP, xODF):
    from pfig import getNRByBlock

    # projOpLists    = systemBlocks['projOpLists']
    nRByBlock = getNRByBlock(pfigsP)
    lbByBlock = num.cumsum(num.hstack( (0, nRByBlock[:-1]) ))
    
    pEI_scaling = num.zeros(num.sum(nRByBlock))
    
    for projOpList, pfigP, iR, nR in zip(projOpLists, pfigsP, lbByBlock, nRByBlock):
        pI  = pfigP['poleVals']
        pA  = num.zeros(nR, dtype=pEI_scaling.dtype)
        contrib = num.zeros( nR, dtype=pEI_scaling.dtype)
        for projOp in projOpList:
            pA = pA + projOp * xODF
        pA_NZ = pA > 0.
        contrib[pA_NZ]  = pI[pA_NZ] / pA[pA_NZ]
        pEI_scaling[iR:iR+nR] = contrib
    return pEI_scaling

def calcPfigProj(systemBlocks, pfigsP, xODF, xE):
    from pfig import getNRByBlock
    from tens import symmToVecds
    
    projOpLists    = systemBlocks['projOpLists']
    cVLists        = systemBlocks['cVLists'] 
    nRByBlock      = getNRByBlock(pfigsP)
    lbByBlock      = num.cumsum(num.hstack( (0, nRByBlock[:-1]) ))
    pEI_scaling    = getEIScaling(projOpLists, pfigsP, xODF)

    vals = []
    for projOpList, cVList, iR, nR in zip(projOpLists, cVLists, lbByBlock, nRByBlock):
        pEIThis = num.zeros(nR)
        for projOp, cV in zip(projOpList, cVList):
            nCV  = cV / num.linalg.norm(cV)
            nnProjVecds = symmToVecds(num.outer(nCV,nCV))
            for iSvec in range(6):
                if num.abs(nnProjVecds[iSvec]) < 1.0e-16:
                    'know contribution will be negligible, so skip it'
                    continue
                contrib = nnProjVecds[iSvec] * ( pEI_scaling[iR:iR+nR] * ( projOp * ( xODF * xE[iSvec,:] ) ) )
                # pEI[iR:iR+nR] = pEI[iR:iR+nR] + contrib
                pEIThis = pEIThis + contrib
        vals.append(pEIThis)
    return vals

def calcPfigProjOps(pfigs, mRodr, **kwargs):
    retval = calcPfigProjOps_pfig(pfigs, mRodr, splitHKLVariants=True, **kwargs)
    return retval

def invert(systemBlocksE, pfigsE, xODF, pfigsP, verbose=True, solverTol=1e-7, maxfun=None):
    """
    pfigsP and pfigsE should have scaleFactor values included (multiplied into these, instead of divided into calcPfigProj results!)
    
    unlike the invert function for the ODF itself, this function does not employ any constraints; eventually one might change this to constrain the average stress state to be uniaxial or the like
    
    strain on the boundary is, like the ODF, force to zero
    
    returns and array with shape (6, xODF.size)
    """
    from scipy import sparse
    from scipy.sparse import linalg
    from scipy.sparse.linalg import LinearOperator, cg
    from tens import symmToVecds
    import copy
    import time
    
    systemBlocks = systemBlocksE
    
    rhs, nRByBlock = buildRHS(pfigsE, getN=True)

    h1Factor   = None
    h1         = None
    if systemBlocks.has_key('h1Factor'):
        h1Factor   = systemBlocks['h1Factor']
        h1         = systemBlocks['H1']
    projOpLists    = systemBlocks['projOpLists']
    cVLists        = systemBlocks['cVLists'] 
    nProjOps       = len(projOpLists)
    lbByBlock      = num.cumsum(num.hstack( (0, nRByBlock[:-1]) ))

    M = rhs.size
    N = projOpLists[0][0].shape[1]
    if systemBlocks['irnp_master'] is not None:
        'check elsewhere that irnp_master == mRodr.getNNodesReduced()-1'
        N = N - 1
        xODFN = xODF[:-1]
        if h1Factor is not None:
            h1 = h1[:-1,:-1]
    else:
        xODFN = xODF
    Nby6 = N
    N = N * 6 # for six components of strain
    lbByN6 = num.arange(0, N, Nby6)
    maxfun = maxfun or 3*N

    if verbose:
        print 'strain pole projection system size is : %d by %d' % (M,N)

    x0 = num.zeros(N)
    
    # could save some work by multiplying through with diagA
    diagA = sparse.spdiags(xODFN, 0, xODFN.size, xODFN.size)
    
    # factors for intensity rescaling
    pEI_scaling = getEIScaling(projOpLists, pfigsP, xODF)
    
    def matt_x_vec(rVec):
        retval = num.zeros_like(x0)
        
        for projOpList, cVList, iR, nR in zip(projOpLists, cVLists, lbByBlock, nRByBlock):
            for projOp, cV in zip(projOpList, cVList):
                nCV  = cV / num.linalg.norm(cV)
                nnProjVecds = symmToVecds(num.outer(nCV,nCV))
                if systemBlocks['irnp_master'] is not None:
                    projOpThis = projOp[:,:-1]
                else:
                    projOpThis = projOp
                for iSvec, iX in enumerate(lbByN6):
                    # note that diagA.T = diagA
                    contrib = nnProjVecds[iSvec] * (
                        diagA * (projOpThis.T * (pEI_scaling[iR:iR+nR] * rVec[iR:iR+nR]))
                        )
                    retval[iX:iX+Nby6] = retval[iX:iX+Nby6] + contrib
                    
        return retval
    
    def mat_x_vec(xVec, withOpT=True, withH1=True):
        xE = xVec.reshape((6, Nby6))
        if systemBlocks['irnp_master'] is not None:
            # pad x with a zero for each strain component
            xE = num.hstack( (xE, num.zeros((6,1)) ) )
        
        pEI = num.zeros_like(rhs)
        # pfigDataXS = calcPfigProj_pfig(sysBlocksP, xODF, applyScaleFactors=False)
        for projOpList, cVList, iR, nR in zip(projOpLists, cVLists, lbByBlock, nRByBlock):
            for projOp, cV in zip(projOpList, cVList):
                nCV  = cV / num.linalg.norm(cV)
                nnProjVecds = symmToVecds(num.outer(nCV,nCV))
                for iSvec in range(6):
                    if num.abs(nnProjVecds[iSvec]) < 1.0e-16:
                        'know contribution will be negligible, so skip it'
                        continue
                    contrib = nnProjVecds[iSvec] * ( pEI_scaling[iR:iR+nR] * ( projOp * ( xODF * xE[iSvec,:] ) ) )
                    pEI[iR:iR+nR] = pEI[iR:iR+nR] + contrib
        if withOpT:
            retval = matt_x_vec(pEI)
            if withH1 and h1Factor is not None:
                for iSvec, iX in enumerate(lbByN6):
                    retval[iX:iX+Nby6] = retval[iX:iX+Nby6] + h1Factor * ( h1 * xE[iSvec,:-1] )
        else:
            retval = pEI
        return retval
    
    A = LinearOperator( (N,N), mat_x_vec, dtype=x0.dtype)
    b = matt_x_vec(rhs)

    if verbose:
        tic = time.time()
        print 'about to run cg ...'
    x, info = cg(A,b,x0=x0,tol=solverTol,maxiter=maxfun)
    if verbose:
        toc = time.time()
        print '... done with CG, took %g seconds' % (toc-tic)

    if not info == 0:
        raise RuntimeError, 'info == '+str(info)
    xE = x.reshape((6, Nby6))
    if systemBlocks['irnp_master'] is not None:
        xE = num.hstack( (xE, num.zeros((6,1)) ) )
    retval = xE

    if verbose:
        resid = mat_x_vec(x, withOpT=False, withH1=False) - rhs
        import math
        relativeErr = math.sqrt(num.dot(resid,resid)) / math.sqrt(num.dot(rhs, rhs))
        print 'relative error : '+str(relativeErr)

    return retval
