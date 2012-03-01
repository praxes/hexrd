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

import arrayUtil
from arrayUtil import getMem
from arrayUtil import num
import pfigUtil

# sysBlockKeys = ['cV', 'hkls', 'scaleFactors', 'hkls_unique', 'projOps', 'h1Factor', 'H1', 'nnr', 'irnp_master', 'hkls_unique', 'hkl_index', 'h1Factor', 'H1']
sysBlockIndirectKeys = ['projOps','projOpLists']

import scipy
haveBFGSDisp = tuple(map(int,scipy.version.version.split('.')[0:2])) >= (0,9)

def calcPfigProj(systemBlocks, x, applyScaleFactors=True):
    """
    if systemBlocks has scaleFactors and applyScaleFactors is True, 
    then divide by scale factors for comparison to _unscaled_ pole figure data
    """
    vals = []
    if systemBlocks.has_key('scaleFactors') and applyScaleFactors:
        systemBlocks['scaleFactors'] = num.array(systemBlocks['scaleFactors'])
        scaleFactors     = systemBlocks['scaleFactors'][systemBlocks['hkl_index']]
        for projOp, scaleFactor in zip( systemBlocks['projOps'], scaleFactors ): 
            vals.append((1./scaleFactor) * (projOp * x))
    else:
        for projOp in systemBlocks['projOps']: 
            vals.append(projOp * x)
    return vals


def makePfigDict(hkl, crystalVector=None, 
                 nVecs=None, 
                 nVecsQP=None, wQP=None,
                 pVals=None, pfigName=None):
    '''
    hkl can be more or less anything, but must be unique to a given
    hkl type if scale factors are to be shared among that type
    properly
    '''
    import arrayUtil
    if crystalVector is None:
        hklFrom = hkl
    else:
        hklFrom = crystalVector
    assert len(hklFrom) == 3, 'length of crystalVector needs to be 3'
    pfigName = pfigName or str(hkl)
    d = {
        'hkl' : hkl,
        'crystalVector' : arrayUtil.toArray(num.array(hklFrom,dtype=float)),
        'pfigName' : pfigName,
        }
    lenV = lenP = 0
    if nVecs is not None:
        assert nVecsQP is None, 'specify only one of nVecs and nVecsQP'
        assert nVecs.shape[0] == 3, 'nVecs should be 3xN in shape'
        d['polePoints'] = nVecs
        lenP = d['polePoints'].shape[1]
    if nVecsQP is not None:
        assert nVecs is None, 'specify only one of nVecs and nVecsQP'
        for nVecsThis in nVecsQP:
            assert nVecsThis.shape[0] == 3, 'nVecs should be 3xN in shape'
        d['polePointsQP'] = nVecsQP
        d['polePoints'] = None
        assert wQP is not None, 'if have nVecsQP, also need wQP'
        d['wQP'] = wQP
        lenP = d['polePointsQP'][0].shape[1]
    if pVals is not None:
        d['poleVals']   = pVals.flatten()
        lenV = len(d['poleVals'])
    
    'sanity checks'
    if lenV != lenP:
        raise RuntimeError, 'mismatch in number of values, points : %d %d' % (lenV, lenP)
    
    return d

def getNVecs(pf):
    assert pf.has_key('polePoints'), 'need point points'
    if pf['polePoints'] is None:
        ppList = pf['polePointsQP']
        wList  = pf['wQP']
        pp = num.zeros_like(pf['polePointsQP'][0])
        for ppThis, wThis in zip(ppList, wList):
            pp += wThis * ppThis
    else:
        pp = pf['polePoints']
    return pp
    
def calcPfigProjOps(pfigs, mRodr, polePoints=None, h1Factor=None, verbose=True,
                    nCPUs=None, nPerChunk=None, # multiprocessing stuff
                    splitHKLVariants=False,
                    ):
    '''
    If polePoints not given, pull out of pfigs
    
    splitHKLVariants=True is meant to be useful when doing strain pole figures
    '''
    from scipy import sparse
    import orientations as ors
    from femODFUtil.mRodrUtil import convertFemODFCoo
    import xrdBase
    useMultiproc = xrdBase.haveMultiProc and (nCPUs is None or nCPUs > 0)
    if useMultiproc:
        if verbose:
            print 'will attempt to use multiprocessing'
        from femODFUtil import pfigMultiProc
        mRodr.setVerbosity(0)

    nnr = mRodr.getNNodesReduced()
    irnp1_master = mRodr.getMasterRNP() # result numbers from 1

    systemBlocks = {}
    systemBlocks['cV'] = []
    systemBlocks['hkls'] = []
    systemBlocks['hkls_unique'] = []
    if splitHKLVariants:
        systemBlocks['projOpLists'] = []
        systemBlocks['cVLists'] = []
    else:
        systemBlocks['projOps'] = []
    systemBlocks['h1Factor'] = None
    systemBlocks['H1'] = None
    systemBlocks['residWeights'] = None
    systemBlocks['nnr'] = nnr
    if irnp1_master > 0:
        irnp_master = irnp1_master - 1
        systemBlocks['irnp_master'] = irnp_master
        assert irnp_master == nnr-1, 'master node is not the highest-numbered reduced node'
    else:
        'no master node'
        systemBlocks['irnp_master'] = None
    
    symmGroupString = mRodr.getSymmGroup()
    symm = ors.SymmGroup(symmGroupString)
    qArray = symm.getQArray()
    
    for pf in pfigs:
        cVMaster = pf['crystalVector']
        hkl = pf['hkl']
        pfigName = pf['pfigName']
        if polePoints is None:
            assert pf.has_key('polePoints'), 'need point points'
            if pf['polePoints'] is None:
                ppList = pf['polePointsQP']
                wList  = pf['wQP']
                nameList = ['qp%d'%(iQP) for iQP in range(len(wList))]
            else:
                ppList = [pf['polePoints']]
                wList  = [None]
                nameList = ['']
        else:
            ppList = [polePoints]
            wList  = [None]
            nameList = ['']
        
        lPM = cvecPMInSym(cVMaster, qArray)
        systemBlocks['cV'].append(cVMaster)
        systemBlocks['hkls'].append(hkl)
        if systemBlocks['hkls_unique'].count(hkl) == 0:
            systemBlocks['hkls_unique'].append(hkl)

        if splitHKLVariants:
            projOpList = []
            from XRD import Symmetry
            cVList = (Symmetry.applySym(cVMaster.reshape(3,1), qArray.T, csFlag=False, cullPM=True)).T
            # len(cVList) is number of specific HKLs
            applySym = False
            pmFactor = 1.0 # because not applying symmetries inside getPoleProjection
            lPM = False # because not applying symmetries inside getPoleProjection
        else:
            cVList = [cVMaster]
            applySym = True
            pmFactor = 0.5
            if verbose:
                print 'for %s, lPM is %s' % (pfigName, str(lPM))
                
        for cV in cVList:
            for pp, name, w in zip(ppList, nameList, wList):
                projTotal = None
                
                if verbose:
                    print 'working on %s %s : (+) %s' % (pfigName, name, str(cV))
                if useMultiproc:
                    Pp = pfigMultiProc.polePathInteg(cV, pp, mRodr, applySym=applySym, nCPUs=nCPUs, nPerChunk=nPerChunk)
                else:
                    odfPfProj = mRodr.getPoleProjection(cV, pp, applySym)
                    Pp = convertFemODFCoo(odfPfProj)
                #
                if not lPM:
                    if verbose:
                        print 'working on %s %s : (-) %s' % (pfigName, name, str(cV))
                    if useMultiproc:
                        Pn = pfigMultiProc.polePathInteg(-cV, pp, mRodr, applySym=applySym, nCPUs=nCPUs, nPerChunk=nPerChunk)
                    else:
                        odfPfProj = mRodr.getPoleProjection(-cV, pp, applySym)
                        Pn = convertFemODFCoo(odfPfProj)
                    #
                    P = pmFactor * (Pp + Pn)
                else:
                    P = Pp
                #
                if w is not None:
                    P = w * P
                if projTotal is None:
                    projTotal = P
                else:
                    projTotal = projTotal + P
            'end of quadrature loop'
            if splitHKLVariants:
                projOpList.append(projTotal)
        'end of cV loop'
        if splitHKLVariants:
            systemBlocks['projOpLists'].append(projOpList)
            systemBlocks['cVLists'].append(cVList)
        else:
            systemBlocks['projOps'].append(projTotal) 
    'end of loop over pfigs'
    
    'constructure hkl_index for use with scaleFactors_hkl or whatever'
    hkls_unique = systemBlocks['hkls_unique']
    hkls        = systemBlocks['hkls']
    hkl_index   = num.array([hkls_unique.index(hkl) for hkl in hkls])
    systemBlocks['hkl_index'] = num.array(hkl_index)
    
    'H1 if requested'
    if h1Factor: 
        addH1ToSystem(systemBlocks, mRodr, h1Factor)
    
    return systemBlocks

def getNPfigPatches(systemBlocks):
    retval = len(systemBlocks['hkl_index'])
    return retval

def setResidWeights(systemBlocks, weights):
    nPfigs = getNPfigPatches(systemBlocks)
    assert len(weights) == nPfigs, \
        'expect length of weights to equal the number of pole figure patches'
    systemBlocks['residWeights'] = num.asarray(weights)
    return

def addH1ToSystem(systemBlocks, mRodr, h1Factor):
    from scipy import sparse
    from femODFUtil.mRodrUtil import convertFemODFCoo
    
    systemBlocks['h1Factor'] = h1Factor
    
    if systemBlocks['H1'] is not None:
        'assume H1 is correct'
        pass
    else:
        opH1 = mRodr.getH1()
        H1 = convertFemODFCoo(opH1)
        # coo supports duplicate entries, but better to use csr or csc for subsequent use;
        # note that entries are summed as desired when converting to csr or csc;
        # nnz may drop substantially when covert from coo to csr or csc
        systemBlocks['H1'] = H1.tocsr()
    
    return

def getNRByBlock(pfigs):
    rhsBlocks = []
    for pf in pfigs:
        pv = pf['poleVals']
        rhsBlocks.append(pv)
    'the following is not necessary given how the solve is done'
    # if systemBlocks['h1Factor'] is not None:
    #     rhsBlocks.append(num.zeros(systemBlocks['nnr']))
    nRByBlock = zip(*map(num.shape, rhsBlocks))[0]
    retval = nRByBlock
    return retval

def buildRHS(pfigs, getN=False):
    rhsBlocks = []
    for pf in pfigs:
        pv = pf['poleVals']
        rhsBlocks.append(pv)
    'the following is not necessary given how the solve is done'
    # if systemBlocks['h1Factor'] is not None:
    #     rhsBlocks.append(num.zeros(systemBlocks['nnr']))
    rhs = num.hstack(rhsBlocks)
    retval = num.array(rhs,dtype='float64')
    if getN:
        nRByBlock = zip(*map(num.shape, rhsBlocks))[0]
        retval = (retval, nRByBlock)
    return retval

def getUniformX(systemBlocks, val):
    mThis, n = systemBlocks['projOps'][0].shape
    
    x = val * num.ones(n, dtype=type(val))
    if systemBlocks['irnp_master'] is not None:
        x[-1] = 0
    
    return x

def invert(systemBlocks, rhs, boundConstrained=True, verbose=True, scaleIntensities=False, 
           solverTol=1.0e-7, maxfun=None):
    from scipy import sparse
    from scipy.sparse import linalg
    from scipy.sparse.linalg import LinearOperator, cg
    import copy
    
    if hasattr(rhs, 'pop'): # type(rhs) == list
        'assume needs to be stacked'
        rhs = num.array(num.hstack(rhs), dtype='float64')
    
    h1Factor   = None
    h1         = None
    if systemBlocks.has_key('h1Factor'):
        h1Factor   = systemBlocks['h1Factor']
        h1         = systemBlocks['H1']
    projOps    = systemBlocks['projOps']
    nProjOps   = len(projOps)
    nRByBlock  = zip(*map(num.shape, projOps))[0]
    lbByBlock  = num.cumsum(num.hstack( (0, nRByBlock[:-1]) ))

    if systemBlocks['residWeights'] is not None:
        residWeights = systemBlocks['residWeights']
        if verbose:
            print 'will use residual scaling factors : '+str(residWeights)
        assert len(residWeights) == len(nRByBlock),\
            'residWeights are wrong length'
        # expandedResidWeights = num.hstack(map( lambda s,n: num.tile(s, n), residWeights, nRByBlock ))
        rhs = copy.copy(rhs)
        projOpsNew = []
        for iBlock, (residWeight, projOp) in enumerate(zip(residWeights, projOps)): 
            iR   = lbByBlock[iBlock]
            nR   = nRByBlock[iBlock]
            # projOp[iR:iR+nR,:] = residWeight * projOp[iR:iR+nR,:]
            projOpsNew.append( projOp * residWeight )
            rhs[iR:iR+nR] = residWeight * rhs[iR:iR+nR]
        projOps = projOpsNew
        
    projOpFull = (sparse.vstack(projOps)).tocsr()

    if systemBlocks['irnp_master'] is not None:
        'check elsewhere that irnp_master == nnr-1'
        projOp = projOpFull[:,:-1]
        if h1Factor is not None:
            h1 = h1[:-1,:-1]
    else:
        projOp = projOpFull
    
    M,N = projOp.shape
    maxfun = maxfun or 3*N
    
    if verbose:
        print 'pole projection system size is : %d by %d' % (M,N)
    
    if h1Factor is not None:
        def matvec(x):
            return projOp.T * (projOp * x) + h1Factor * (h1 * x)
    else:
        def matvec(x):
            return projOp.T * (projOp * x)
    A = LinearOperator( (N,N), matvec, dtype=projOp.dtype)

    if systemBlocks.has_key('scaleFactors'):
        systemBlocks['scaleFactors'] = num.array(systemBlocks['scaleFactors'])
        scaleFactorsHKLRef = systemBlocks['scaleFactors']
        scaleFactorsRef = scaleFactorsHKLRef[systemBlocks['hkl_index']]
        expandedSFS = num.hstack(map( lambda s,n: num.tile(s, n), scaleFactorsRef, nRByBlock ))
        rhsScaled = expandedSFS * rhs
        b = projOp.T * rhsScaled
    else:
        scaleFactorsHKLRef = num.tile(1., nProjOps)
        rhsScaled = rhs
        b = projOp.T * rhsScaled
    
    if verbose:
        print 'using CG to solve the unconstrained problem ...'
    xUnBnd, info = cg(A,b,tol=solverTol,maxiter=maxfun)
    if verbose:
        print '... done with CG'
    
    if not info == 0:
        raise RuntimeError, 'info == '+str(info)
    if systemBlocks['irnp_master'] is not None:
        x = num.hstack((xUnBnd, 0.))
    else:
        x = xUnBnd
    
    if verbose:
        print 'for problem without bound constraints or solving for scalings: '
        print '    min, max, and mean of solution : %g %g %g' % (x.min(), x.max(), x.mean())
        resid = projOpFull * x - rhsScaled
        import math
        relativeErr = math.sqrt(num.dot(resid,resid)) / math.sqrt(num.dot(rhsScaled, rhsScaled))
        print '    relative error : '+str(relativeErr)
    
    if not boundConstrained and not scaleIntensities:
        return x
    
    'solve using unconstrained problem without variable scaleFactors as a starting point'
    from scipy import optimize 

    nSFDOF = None
    if scaleIntensities:
        """
        could condense the scale factor degrees of freedom out of the system, 
        but there should not be that many extra variables due to the scale factors 
        and leaving them as DOF may keep things cleaner
        """
        
        # nSFDOF = nProjOps-1
        nSFDOF = len(systemBlocks['hkls_unique'])-1
        
        def evalResid(xFull, *args):
            x     = xFull[:-nSFDOF]
            sfDOF = xFull[-nSFDOF:]
            scaleFactorsHKLMod = num.hstack( (1., sfDOF) )
            scaleFactorsMod = scaleFactorsHKLMod[systemBlocks['hkl_index']]
            expandedSFDOF = num.hstack(map( lambda s,n: num.tile(s, n), scaleFactorsMod, nRByBlock ))
            resid = projOp * x - (expandedSFDOF * rhsScaled)
            return resid
        if h1Factor is not None:
            def func(xFull, *args):
                resid = evalResid(xFull, *args)
                x     = xFull[:-nSFDOF]
                retval = 0.5 * num.dot(resid, resid) + (h1Factor * 0.5) * num.dot(x, h1 * x)
                return retval
            def fprime(xFull, *args): 
                resid = evalResid(xFull, *args)
                x     = xFull[:-nSFDOF]
                retvalX = projOp.T * resid + h1Factor * (h1 * x) # assums h1 is symmetric
                retvalW = num.zeros(nSFDOF)
                for iBlock, hklIndex in enumerate(systemBlocks['hkl_index']):
                    if hklIndex == 0:
                        'first hkl type does not have a scaling DOF'
                        continue
                    iSFDOF = hklIndex-1
                    iR   = lbByBlock[iBlock]
                    nR   = nRByBlock[iBlock]
                    contrib = -num.dot( resid[iR:iR+nR], rhsScaled[iR:iR+nR] )
                    retvalW[iSFDOF] += contrib
                retval = num.hstack( (retvalX, retvalW) )
                return retval
        else:
            def func(xFull, *args):
                resid = evalResid(xFull, *args)
                retval = 0.5 * num.dot(resid, resid)
                return retval
            def fprime(xFull, *args): 
                resid = evalResid(xFull, *args)
                retvalX = projOp.T * resid
                retvalW = num.zeros(nSFDOF)
                for iBlock, hklIndex in enumerate(systemBlocks['hkl_index']):
                    if hklIndex == 0:
                        'first hkl type does not have a scaling DOF'
                        continue
                    iSFDOF = hklIndex-1
                    iR   = lbByBlock[iBlock]
                    nR   = nRByBlock[iBlock]
                    contrib = -num.dot( resid[iR:iR+nR], rhsScaled[iR:iR+nR] )
                    retvalW[iSFDOF] += contrib
                retval = num.hstack( (retvalX, retvalW) )
                return retval
            
    else:
        def evalResid(x, *args):
            resid = projOp * x - rhsScaled
            return resid
        if h1Factor is not None:
            def func(x, *args):
                resid = evalResid(x, *args)
                retval = 0.5 * num.dot(resid, resid) + (h1Factor * 0.5) * num.dot(x, h1 * x)
                return retval
            def fprime(x, *args): 
                resid = evalResid(x, *args)
                retval = projOp.T * resid + h1Factor * (h1 * x) # assums h1 is symmetric
                return retval
        else:
            def func(x, *args):
                resid = evalResid(x, *args)
                retval = 0.5 * num.dot(resid, resid)
                return retval
            def fprime(x, *args): 
                resid = evalResid(x, *args)
                retval = projOp.T * resid
                return retval
    
    x0 = copy.copy(xUnBnd)
    minVal0 = 1.0e-8 * xUnBnd.max()
    x0[num.where(x0 < minVal0)] = minVal0
    x0Ref = num.zeros(len(x0))
    if scaleIntensities:
        x0    = num.hstack( (x0,    num.tile(1., nSFDOF)) )
        x0Ref = num.hstack( (x0Ref, num.tile(1., nSFDOF)) )
    # pgtol default is 1.0e-5
    # pgtol = solverTol / func( x0Ref )
    pgtol = solverTol
    eps = num.finfo(float).eps # machine precision
    factr = solverTol / eps
    if verbose:
        print 'using pgtol of %g' % (pgtol)
        print 'using factr of %g' % (factr)
    bfgsKWArgs = {}
    if haveBFGSDisp:
        disp = 0
        if verbose: disp = 2
        bfgsKWArgs['disp'] = disp
    if boundConstrained:
        bounds = [ (0., None) for i in range(len(xUnBnd)) ]
        if scaleIntensities:
            bounds += [ (0., None) for i in range(nSFDOF) ]
    else:
        bounds = None
    xBnd, fVal, d = optimize.lbfgsb.fmin_l_bfgs_b(func, x0, bounds=bounds, fprime=fprime, 
                                                  pgtol=pgtol, factr=factr, maxfun=maxfun, **bfgsKWArgs)
    if not d['warnflag'] == 0:
        import sys
        print >> sys.stderr, 'trouble in lbfgsb : ' + str(d['task'])
    
    resid = None
    if verbose:
        resid    = evalResid(xBnd)
        x0       = num.zeros_like(xBnd)
        x0[:-nSFDOF] = xBnd[:-nSFDOF]
        residRef = evalResid(x0)
    
    if scaleIntensities:
        'pull scale factor degrees of freedom off the end of xBnd'
        sfDOF = xBnd[-nSFDOF:]
        xBnd = xBnd[:-nSFDOF]
        scaleFactorsHKLMod = scaleFactorsHKLRef * num.hstack( (1., sfDOF) )
        scaleFactorsMod = scaleFactorsHKLMod[systemBlocks['hkl_index']]
        # 'store in systemBlocks as well as returning'
        # systemBlocks['scaleFactors'] = scaleFactorsHKLMod
        'check that satisfy analytic form of weight factor to some tolerance'
        ww = num.zeros(nSFDOF)
        bb = num.zeros(nSFDOF)
        for iBlock, hklIndex in enumerate(systemBlocks['hkl_index']):
            if hklIndex == 0:
                'first hkl type does not have a scaling DOF'
                continue
            iSFDOF = hklIndex-1
            sfDOFThis = sfDOF[iSFDOF]
            iR   = lbByBlock[iBlock]
            nR   = nRByBlock[iBlock]
            projOpThis = projOp[iR:iR+nR,:]
            rhsThisBlock = rhsScaled[iR:iR+nR]
            bb[iSFDOF] += num.dot( rhsThisBlock, rhsThisBlock )
            ww[iSFDOF] += num.dot( rhsThisBlock , projOpThis * xBnd )
        for iSFDOF in range(nSFDOF):
            hklIndex = iSFDOF+1
            sfDOFThis = sfDOF[iSFDOF]
            if bb[iSFDOF] > 0.:
                wThis = ww[iSFDOF] / bb[iSFDOF]
                print 'analytic versus solved weight : %g vs %g' % (wThis, sfDOFThis)
                if (wThis - sfDOFThis) > min(solverTol*1e2,1e-3) * wThis:
                    import sys
                    print >> sys.stderr, 'error in weight factor solution seems too high'
    
    if systemBlocks['irnp_master'] is not None:
        x = num.hstack((xBnd, 0.))
    else:
        x = xBnd
    if scaleIntensities:
        retval = (x, scaleFactorsHKLMod)
    else:
        retval = x
    
    if verbose:
        print 'for problem boundConstrained %s and scaleIntensities %s : ' % (boundConstrained, scaleIntensities)
        print '    min, max, and mean of solution : %g %g %g' % (x.min(), x.max(), x.mean())
        if scaleIntensities:
            print '    scale Factors : %s' % (str(scaleFactorsHKLMod))
        import math
        # relativeErr = math.sqrt(num.dot(resid,resid)) / math.sqrt(num.dot(rhsScaled, rhsScaled))
        relativeErr = math.sqrt(num.dot(resid,resid)) / math.sqrt(num.dot(residRef, residRef))
        print '    relative error : '+str(relativeErr)
    
    return retval
    

def getAutoMesh(angRes):
    
    import femODF.FemHemisphere
    mPfig = femODF.FemHemisphere.FemHemisphere()
    mPfig.autoMesh(angRes)
    
    coords = mPfig.getCoords()
    from femODFUtil import pfig
    nVectors = pfigUtil.sph2n(coords)
    
    return (mPfig, nVectors)

def getAutoMeshInv(angRes, symmGroupString):
    #import femODF.FemMesh
    import femODF.FemHemisphere
    import femODF.ElemType
    from math import atan2, sqrt, sin, cos
    
    elemType = femODF.ElemType.stdT3
    quadRule = arrayUtil.toArray([22004, 11002])
    
    #mInvPfig = femODF.FemMesh.FemMesh()
    mInvPfig = femODF.FemHemisphere.FemHemisphere()
    
    if symmGroupString == 'cub':
        'based on liverne_input/invpfig/cub/mkmesh.awk'

        ang_max = atan2(sqrt(2.),1.);
        nela = int(ang_max/angRes)+1;
        dang = ang_max/nela;
        
        # make nodes
        nnodes=0;
        ang=0.;
        coordList = []
        for i in range(nela+1):
            ldist = sin(ang)/cos(ang);
            x = ldist/sqrt(2);
        
            y = 0.; dely = 1.;
            if (i>0): dely = x/i;
            for j in range(i+1):
                # nnodes += 1
                # coordList.append([x, y, 1.]);
                coordList.append([x, y]);
                nnodes += 1
                y += dely;
        
            ang += dang;

        coords = arrayUtil.toArray(coordList).T

        # make elements; connectivity numbered from 1;
        nelms = 0;
        in_base_l = 0; in_base_r = 1;
        connList = []
        for i in range(1,nela+1):
            for j in range(1,i):
                nelms += 1
                connList.append([
                        in_base_l+(j-1)+1,
                        in_base_r+(j-1)+1,
                        in_base_r+(j-1)+2
                        ])
                nelms += 1
                connList.append([
                        in_base_r+(j-1)+2,
                        in_base_l+(j-1)+2,
                        in_base_l+(j-1)+1
                        ])
        
            nelms += 1
            connList.append([
                    in_base_l+(i-1)+1,
                    in_base_r+(i-1)+1,
                    in_base_r+(i-1)+2
                    ])
        
            in_base_l += i;
            in_base_r += i+1;

        conn = arrayUtil.toArray(connList).T

    elif symmGroupString == 'hex':
        'depends on covera'
        'see liverne_input/invpfig/hex/mkmake.awk'
        raise RuntimeError, 'getAutoMeshInv not coded for '+symmGroupString
    else:
        raise RuntimeError, 'getAutoMeshInv not coded for '+symmGroupString
    
    # make unit vectors
    zs = num.ones(coords.shape[1])
    vectors = num.vstack((coords, zs)).T
    from numpy import linalg
    mags = num.array(map(linalg.norm, vectors))
    nVectors = arrayUtil.toArray([vectors[iV]/mags[iV] for iV in range(len(vectors))]).T
    
    # replace coords, which are funny things sitting on the z=1 plane,
    # with their spherical coordinate versions
    coords = pfigUtil.n2sph(nVectors)
    
    mInvPfig.init(elemType, quadRule, conn, coords)
    
    return (mInvPfig, nVectors)

def plotPfigDA(quatsData, 
               symmGroupString, # 'cub'
               crystalVectors, 
               dataFilePrefix=None,
               weightsData=None,
               names=None,
               angTol=0.1, angRes=0.01, 
               plotList=None,
               **args):
    '''All-in-one utility;
    '''
    import sys
    import dxWrap
    
    nCV = crystalVectors.shape[1]

    mPfig, nVectors = getAutoMesh(angRes)
    
    import femODF.DiscreteAggregate
    dagg = femODF.DiscreteAggregate.DiscreteAggregate()
    if weightsData is not None:
        dagg.setOrientationsWeighted(quatsData, weightsData)
    else:
        dagg.setOrientations(quatsData)
  
    pFigVals = dagg.calculatePoleValues(crystalVectors, nVectors, symmGroupString, angTol)
    
    imageFileNames = dxWrap.plotPfigs(pFigVals.T, mPfig, plotList, dataFilePrefix, names=names, **args)

    return imageFileNames
    
def genNvecSph(*args, **keyArgs):
    '''theta is down from pole, phi is around azimuth;
    both can be made using numpy arange
    '''
    import numpy as num
    import math
    
    deg = True
    if keyArgs.has_key('deg'):
        deg = keyArgs.pop('deg')
    if len(keyArgs) > 0:
        raise RuntimeError, 'unprocessed keyword args : '+str(keyArgs)
    
    degToRad = math.pi / 180.
    if len(args) == 1:
        dTheta = args[0]
        if deg: dTheta = dTheta * degToRad
        thetaVals = num.arange(0., math.pi/2.0+dTheta/2.0, dTheta)
        phiVals   = num.arange(0., 2.0*math.pi, dTheta) 
    elif len(args) == 2:
        thetaVals = args[0]
        phiVals   = args[1]
        if deg:
            thetaVals = thetaVals * degToRad
            phiVals   = phiVals   * degToRed
    else:
        raise RuntimeError, 'do not know what to do with args : '+str(args)
    
    nAz = len(phiVals)
    cosPhi = num.cos(phiVals)
    sinPhi = num.sin(phiVals)
    polePoints = getMem((3,len(thetaVals)*nAz))
    iBase = 0
    for theta in thetaVals:
        z = math.cos(theta)
        r = math.sin(theta)
        polePoints[2,iBase:iBase+nAz] = z
        for phi in phiVals:
            polePoints[0,iBase:iBase+nAz] = r * cosPhi
            polePoints[1,iBase:iBase+nAz] = r * sinPhi
        iBase += nAz
    return (polePoints, len(phiVals), len(thetaVals))

def cvecPMInSym(vec, qArray):
    '''
    this is a bit of a mongrel routine in that it
    uses stuff from XRD;
    
    qArray needs to be shape (*,8)
    '''
    from XRD import Symmetry
    cvec = vec.reshape(3,1)
    cvecs1 = Symmetry.applySym(cvec, qArray.T, csFlag=False, cullPM=False)
    cvecs2 = Symmetry.applySym(cvec, qArray.T, csFlag=True,  cullPM=False)
    retval = cvecs2.shape[1] == cvecs1.shape[1]
    return retval

def main(argv):
    operation = argv[0]
    if operation == 'pfigDA':
        quatsFile = argv[1]
        pfigDataFile = argv[2]
        prefix = argv[3]
        angTol = 0.1
        angRes = 0.06
        if len(argv)>4: angTol = float(argv[4])
        if len(argv)>5: angRes = float(argv[5])
        
        import os, sys
        dataFilePath = os.path.split(pfigDataFile)[0]
        pfigDataFileName = os.path.split(pfigDataFile)[1]
        sys.path.append(dataFilePath)
        
        import fileUtil
        quatsData = num.array(fileUtil.readFloatData(quatsFile)).T
        
        exec('import %s as pfigData' % (os.path.splitext(pfigDataFileName)[0]))
        crystalVectors   = num.array([cV for inFix,cV in pfigData.poleInfo]).T
        names            = [prefix+inFix for inFix,cV in pfigData.poleInfo]
        
        import pfig
        pfig.plotPfigDA(quatsData,
                        pfigData.symmGroupString,
                        crystalVectors, 
                        names=names,
                        angTol=angTol, angRes=angRes)
    else:
        raise RuntimeError, 'bad operation '+str(operation)
    return

def loadSysBlocks(prefix):
    """
    returns None if the data are not found
    """
    import shelve
    try:
        s = shelve.open(prefix+'.py_shelf', 'r')
    except:
        return None
    
    sysBlocks = {}
    keys = s['_keys'] # s.keys()
    for key in keys:
        print 'reading '+str(key)
        if sysBlockIndirectKeys.count(key) > 0:
            if key == 'projOps':
                from scipy.io import loadmat, savemat
                d = loadmat(prefix+"_"+key+".mat", appendmat=False)
                projOps = []
                hkls = s['hkls']
                for iProjOp in range(len(hkls)): # cannot do len(d.keys()) because it has other stuff
                    projOps.append(d[str(iProjOp)])
                value = projOps
            elif key == 'projOpLists':
                from scipy.io import loadmat, savemat
                d = loadmat(prefix+"_"+key+".mat", appendmat=False)
                projOpLists = []
                cVLists = s['cVLists']
                for iProjOpList, cVList in enumerate(cVLists):
                    projOpList = []
                    for iProjOp, cV in enumerate(cVList):
                        projOpList.append(d['%d_%d' % (iProjOpList, iProjOp)])
                    projOpLists.append(projOpList)
                value = projOpLists
            else:
                raise RuntimeError, 'do not know what to do with key '+str(key)
            sysBlocks[key] = value
        else:
            sysBlocks[key] = s[key]
    s.close()
    return sysBlocks

def storeSysBlocks(prefix, sysBlocks):
    ''' implemented mostly because having trouble with shelve and
    scipy.sparse objects
    '''
    import shelve, copy
    s = shelve.open(prefix+'.py_shelf', 'c', protocol=2)
    
    keys = copy.deepcopy(sysBlocks.keys())
    s['_keys'] = keys
    
    for key, value in sysBlocks.iteritems():
        # if sysBlockKeys.count(key) == 0:
        #     raise RuntimeError, 'unknown key: '+str(key)
        print 'writing '+str(key)
        if sysBlockIndirectKeys.count(key) > 0:
            if key == 'projOps':
                from scipy.io import loadmat, savemat
                # d = dict(enumerate(value)) # 'does not work because savemeat seems to want string-valued keys'
                d = dict(zip(map(str, range(len(value))), value))
                savemat(prefix+"_"+key+".mat", d, appendmat=False)
            elif key == 'projOpLists':
                from scipy.io import loadmat, savemat
                projOpLists = value
                'loadmat does not work with matrices with no entries, so make sure have at least one entry'
                for projOpList in projOpLists:
                    for projOp in projOpList:
                        if projOp.nnz == 0:
                            projOp[0,0] = 0.
                d = dict( 
                    reduce( 
                        lambda x,y:x+y, 
                        [[ ('%d_%d' % (iProjOpList,iProjOp), projOp) for iProjOp, projOp in enumerate(projOpList)] for iProjOpList, projOpList in enumerate(projOpLists)] 
                        ) )
                savemat(prefix+"_"+key+".mat", d, appendmat=False, oned_as='row')
            else:
                raise RuntimeError, 'do not know what to do with key '+str(key)
                # for iThing, thing in enumerate(value):
                #     thisKey = key+'_%d'%(iThing)
                #     print 'writing '+str(thisKey)
                #     s[thisKey] = thing
            s[key] = True
        else:
            s[key] = value
        #     s[key] = value
        # elif key == 'H1':
        #     H1 = sysBlocks[key]
        #     dok = H1.todok()
        #     keys = dok.keys()
        #     vals = dok.values()
        #     ...
        # elif key == 'projOps':
        #     ...
        # else:
        #     raise RuntimeError, 'unknown key: '+str(key)
    s.close()
    return

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])


