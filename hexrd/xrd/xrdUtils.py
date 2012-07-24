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
import copy
import time
import math
from math import pi
import shelve

import numpy as num
from scipy import sparse
from scipy.linalg import svd
from scipy import ndimage

import matplotlib
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import cm, colors
import matplotlib.collections as collections

from hexrd import plotWrap
from hexrd import tens
from hexrd import matrixutil
from hexrd import pfigUtil
from hexrd.valUnits import toFloat
import hexrd.orientations as ors
from hexrd.xrd import crystallography
from hexrd.xrd.crystallography import latticeParameters, latticeVectors
from hexrd.xrd import detector
from hexrd.xrd.detector import Reader
from hexrd.xrd.detector import getCMap
from hexrd.xrd import xrdBase
from hexrd.xrd.xrdBase import dataToFrame
from hexrd.xrd.xrdBase import multiprocessing
from hexrd.xrd import Rotations as rot
from hexrd.xrd.Rotations import mapAngle
from hexrd.xrd import spotFinder

'quadr1d of 8 is probably overkill, but we are in 1d so it is inexpensive'
quadr1dDflt = 8
#
debugDflt = False

d2r = piby180 = num.pi/180.
r2d = 1.0/d2r

class FormatEtaOme:
    'for plotting data as a matrix, with ijAsXY=True'
    def __init__(self, etas, omes, A, T=False, debug=False):
        self.etas  = etas
        self.omes  = omes
        self.A     = A
        self.debug = debug
        self.T     = T
        return
    def __call__(self,x,y):
        if self.T:
            iX = int(round(y))
            iY = int(round(x))
        else:
            iX = int(round(x))
            iY = int(round(y))
        eta = r2d*(self.etas[iY])
        ome = r2d*(self.omes[iX])
        val = self.A[iY,iX]
        retval = 'eta=%g; ome=%g; val=%g' % \
            (eta, ome, val)
        if self.debug:
            retval += "; x_image=%g y_image=%g" % (x, y)
        return retval

def fitLParm(data, detectorGeom, planeData,
             funcType = detector.funcTypeDflt,
             lParm0 = None,
             funcXVecList = None,
             quadr1d = quadr1dDflt,
             debug = debugDflt,
             ):
    """
    fit lattice parameters to data; using dataToFrame to map the data to a frame

    input planeData is not changed

    """
    frame = dataToFrame(data)

    if funcXVecList is None:
        funcXVecList = True

    'copy planeData, so that the instance passed in does not getting changed'
    #pd = crystallography.PlaneData(planeData.hkls, planeData.getParams())
    pd = copy.deepcopy(planeData)
    if lParm0 is not None:
        pd.lparms = lParm0

    mRing = detector.MultiRingEval(detectorGeom, pd, dataFrame=frame,
                                   funcType=funcType,
                                   refineParamsDG=False,
                                   refineParamsL=True,
                                   funcXVecList=funcXVecList)
    if not hasattr(funcXVecList, '__len__'):
        'still need to get funcXVecList; want something better than default guess'
        funcXVecList, cFitList, pw = \
            mRing.radialFitXVec(plot=debug,
                                      plotTitlePrefix='Before lparm fit: ',
                                      quadr1d=quadr1d)
        mRing.setFuncXVecList(funcXVecList)
    if debug:
        print 'number of parameters in geom-only fit : %d' % (mRing.getNParam())

    x = mRing.doFit()
    'return refined lattice parameters'
    retval = pd.lparms

    if debug:
        funcEval = mRing(x, makePlots=True)

    return retval

def fitDG(data, detectorGeom, planeData,
          funcType = detector.funcTypeDflt,
          funcXVecList = None,
          quadr1d = quadr1dDflt,
          debug = debugDflt,
          ):
    """
    fit detector geometry parameters to data; using dataToFrame to map the data to a frame

    pass funcXVecList as True or as something like detectorGeom.fitRingsFunc.getFuncXVecList()
    if want to just refine detector geometry and not the functional forms for the rings

    input detectorGeom is used to guess parameters and is not modified -- a new
    detectorGeom is returned
    """

    frame = dataToFrame(data)

    'do the following to keep ring functional forms from being refined'
    #if funcXVecList is None:
    #    funcXVecList = True

    # make a new instance to modify, so that detectorGeom is not changed
    dg = detectorGeom.makeNew()

    mRing = detector.MultiRingEval(dg, planeData, dataFrame=frame,
                                   funcType=funcType,
                                   refineParamsDG=True,
                                   refineParamsL=False,
                                   funcXVecList=funcXVecList)
    if funcXVecList and not hasattr(funcXVecList, '__len__'):
        'still need to get funcXVecList; want something better than default guess'
        funcXVecList, cFitList, pw = \
            mRing.radialFitXVec(plot=debug,
                                      plotTitlePrefix='Before lparm fit: ',
                                      quadr1d=quadr1d)
        mRing.setFuncXVecList(funcXVecList)
    if debug:
        print 'number of parameters in fit : %d' % (mRing.getNParam())

    x = mRing.doFit()
    if funcXVecList:
        'return just detector geom'
        retval = dg
    else:
        retval = dg, mRing.getFuncXVecList()

    if debug:
        funcEval = mRing(x, makePlots=True)

    return retval

def fitDGX(data, detectorGeom, planeData,
          funcType = detector.funcTypeDflt,
          quadr1d = quadr1dDflt,
          debug = debugDflt,
          nGlIter = 2,
          xFuncs = None,
          xDG    = None,
          ):
    """
    fit detector geometry parameters to data; using dataToFrame to map the data to a frame

    uses a procedure that might end up being more robust than fitDG

    input detectorGeom is used to guess parameters and is not modified -- a new
    detectorGeom is returned
    """

    'compute these up front so that they can be shared'
    indicesList, iHKLLists = detectorGeom.makeIndicesTThRanges(planeData, cullDupl=True)

    frame = dataToFrame(data)

    # make a new instance to modify, so that detectorGeom is not changed
    dg = detectorGeom.makeNew()

    'mRing for just the DOF in the ring functions'
    mRingFuncs = detector.MultiRingEval(dg, planeData, dataFrame=frame,
                                        indicesList=indicesList, iHKLLists=iHKLLists,
                                        funcType=funcType,
                                        refineParamsDG=False,
                                        refineParamsL=False,
                                        funcXVecList=None)
    if debug:
        print 'number of parameters in ring-funciton fit : %d' % (mRingFuncs.getNParam())

    'mRing for just the DG DOF'
    mRingDG = detector.MultiRingEval(dg, planeData, dataFrame=frame,
                                        indicesList=indicesList, iHKLLists=iHKLLists,
                                        funcType=funcType,
                                        refineParamsDG=True,
                                        refineParamsL=False,
                                        funcXVecList=mRingFuncs.getFuncXVecList(mRingFuncs.xVecGuess))
    if debug:
        print 'number of parameters in dg fit : %d' % (mRingDG.getNParam())

    for iGlIter in range(nGlIter):

        if debug:
            print 'doing mRingFuncs fit'
        xFuncs = mRingFuncs.doFit(xVec0=xFuncs)
        if debug:
            print '    got : '+str(xFuncs)
        mRingDG.setFuncXVecList(mRingFuncs.getFuncXVecList(xFuncs))

        if debug:
            print 'doing mRingDG fit'
        xDG    = mRingDG.doFit(xVec0=xDG)
        if debug:
            print '    got : '+str(xDG)
        mRingFuncs.detectorGeom = mRingDG.detectorGeom

    dg.clean() # free up memory

    if debug:
        funcEval = mRingDG(xDG, makePlots=True)

    return dg, xFuncs, xDG


def textureToSpots(texture,
                   planeData, detectorGeom, omeMin, omeMax,
                   pVecs = None,
                   ):
    """
    take texture as returned from pyMps and make spots
    """

    stretches = None

    quats = texture['quats']
    nGrains = quats.shape[1]
    if texture.has_key('Vs'):
        stretches = texture['Vs'].T

    #
    fMatRef = planeData.latVecOps['F']
    if stretches is not None:
        bMats = num.zeros([nGrains,3,3])
        for iGrain, stretch in enumerate(stretches):
            bMats[iGrain,:,:] = stretchToLV(tens.svecToSymm(stretch), fMatRef)['B']
    else:
        bMats = None
    rMats = rot.rotMatOfQuat(quats)
    #
    spotAngs = makeSynthSpots(rMats, pVecs, bMats,
                              planeData, detectorGeom, omeMin, omeMax)

    return spotAngs

def makeSynthSpots(rMats, pVecs, bMats, planeData, detectorGeom, omeMin, omeMax, hklList=None, beamSize=None):
    """
    make synthetic spots
    """
    if beamSize is not None:
        "the beam size shoud by [width (X), hieght (Y)]"
        omeAxis = num.c_[0, 1, 0].T     # ... may need to grab this from some global to be same!
        pass
    nGrains = rMats.shape[0]
    assert rMats.shape[1] == 3 and rMats.shape[2] == 3,\
        'rMats is wrong shape'

    omeMin = toFloat(num.atleast_1d(omeMin), 'radians')
    omeMax = toFloat(num.atleast_1d(omeMax), 'radians')

    spotAngs = []
    completeGrains = num.ones(nGrains, dtype='bool')
    for i in range(nGrains):

        if bMats is None:
            bMat = None
        else:
            bMat = bMats[i, :, :]
        # full set of symmetric HKLs subject to exclusions
        if hklList is None:
            qVec, qAng0, qAng1 = \
                planeData.makeAllScatteringVectors(rMats[i, :, :], bMat=bMat)
        else:
            qVec, qAng0, qAng1 = \
                planeData.makeTheseScatteringVectors(hklList, rMats[i, :, :], bMat=bMat)

        validAng0 = num.zeros(qVec.shape[1], dtype='bool')
        validAng1 = num.zeros(qVec.shape[1], dtype='bool')
        for j in range(len(omeMin)):
            validAng0 = validAng0 | (
                ( qAng0[2, :] >= omeMin[j] ) & ( qAng0[2, :] <= omeMax[j] ) )
            validAng1 = validAng1 | (
                ( qAng1[2, :] >= omeMin[j] ) & ( qAng1[2, :] <= omeMax[j] ) )

        tmpAngs = num.hstack( [
                qAng0[:, validAng0], qAng1[:, validAng1] ] )

        if pVecs is not None:
            tmpDG = detectorGeom.makeNew(pVec=pVecs[i, :])

            if beamSize is not None:
                numGrainSpots = tmpAngs.shape[1] # number of spots in this grain
                inBeam = num.zeros(numGrainSpots, dtype='bool')
                for k in range(numGrainSpots):
                    comCrds = num.dot(rot.rotMatOfExpMap(tmpAngs[2, k]*omeAxis), pVecs[i, :])
                    if (comCrds[0] >= -0.5*beamSize[0] and comCrds[0] <= 0.5*beamSize[0]) and \
                       (comCrds[1] >= -0.5*beamSize[1] and comCrds[1] <= 0.5*beamSize[1]):
                        inBeam[k] = True
                    else:
                        completeGrains[i] = False
                        pass
                    pass
                if num.any(inBeam):
                    tmpX, tmpY, tmpO = tmpDG.angToXYO(tmpAngs[0, inBeam],
                                                      tmpAngs[1, inBeam],
                                                      tmpAngs[2, inBeam])
                else:
                    continue
            else:
                tmpX, tmpY, tmpO = tmpDG.angToXYO(tmpAngs[0, :],
                                                  tmpAngs[1, :],
                                                  tmpAngs[2, :])
                pass
            tmpTTh, tmpEta, tmpOme = detectorGeom.xyoToAng(tmpX, tmpY, tmpO)
            tmpAngs = num.vstack([tmpTTh, tmpEta, tmpOme])

        spotAngs += list(tmpAngs.T)
        # inefficient:
        # if i == 0:
        #     spotAngs = tmpAngs
        # else:
        #     spotAngs = num.hstack( [
        #         spotAngs, tmpAngs ] )

    retval = num.array(spotAngs).T
    return retval, completeGrains

markerListDflt = [
    'D','o','p','s','v','x', #
    '+','*',
    '<','>',
    '1','2','3','4',
    '^',
    ]
# 'h', 'd', 'H' # hard to distinguish
# ',' , '.' # too hard to see
# '_' # odd

def makePathVariantPoles(rMatRef, fromPhase,
                         pathList,
                         planeDataDict,
                         hklID,
                         ):

    qVecList = []
    for pathName, rMatTransfList, phaseIDTo, phaseIDFrom, pRefSvecList in pathList:
        if not phaseIDFrom == fromPhase:
            'this path not of interest'
            continue
        planeData = planeDataDict[phaseIDTo]

        for iRMat, rMatTransf in enumerate(rMatTransfList):
            rMat = num.dot(rMatRef, rMatTransf)
            qVec, qAng0, qAng1 = \
                planeData.makeTheseScatteringVectors(hklID, rMat)
            qVecList.append(qVec)

    return qVecList

def displayPathVariants(data, rMatRef, fromPhase,
                        pathList,
                        planeDataDict,
                        detectorGeom, omeMin, omeMax,
                        phaseForDfltPD=None,
                        markerList = markerListDflt,
                        hklList = None,
                        color=None,
                        pointKWArgs={},
                        hklIDs=None, pw=None):
    '''
    '''
    pw, local = plotWrap.argToPW(pw)
    if local and data is not None:
        frame = dataToFrame(data, sumImg=num.maximum)
        planeData=None
        if phaseForDfltPD is not None:
            planeData = planeDataDict[phaseForDfltPD]
        detectorGeom.display(frame, planeData, pw=pw)

    col = color
    if col is None:
        col = 'r'
    iMarker = 0
    for pathName, rMatTransfList, phaseIDTo, phaseIDFrom, pRefSvecList in pathList:
        if not phaseIDFrom == fromPhase:
            'this path not of interest'
            continue
        planeData = planeDataDict[phaseIDTo]
        if hasattr(color,'keys'):
            col = color[phaseIDTo]

        for iRMat, rMatTransf in enumerate(rMatTransfList):
            rMat    = num.dot(rMatRef, rMatTransf)
            pVecs = None
            bMats = None
            spotAngs = makeSynthSpots(rMat.reshape([1,3,3]),
                                      pVecs, bMats, planeData, detectorGeom,
                                      omeMin, omeMax,
                                      hklList=hklList,
                                      ) # 3,nSpots
            if spotAngs.size > 0:
                x, y, o = detectorGeom.angToXYO(spotAngs[0,:], spotAngs[1,:], spotAngs[2,:])
                #
                if iMarker == len(markerList):
                    raise RuntimeError, 'not enough markers in list'
                marker = markerList[iMarker]
                kwArgs = {'ls':'None', 'mec':col, 'mfc':'None', 'ms':3.}
                kwArgs.update(pointKWArgs)
                kwArgs.update({'marker':marker})
                pw(x, y, **kwArgs) # style=
            iMarker += 1

    return pw

def makeRMatList(pathList, fromPhase):
    rMatList = []
    for pathName, rMatTransfList, phaseIDTo, phaseIDFrom, pRefSvecList in pathList:
        if not phaseIDFrom == fromPhase:
            'this path not of interest'
            continue
        for iRMat, rMatTransf in enumerate(rMatTransfList):
            rMatList.append( num.dot(rMat, rMatTransf) )
    return rMatList

def findGrainsNewSpots(grainList0, spots, doFitting=True, minCompleteness=None, pathList=[],
                       eosDict=None, refPDDict=None,
                       tK=300.,
                       debug=True,
                       indepFitPVec=False, # fit precession independently
                       findByPixelDist=1, # meant to deal better with overlapped spots
                       maxIterRefit=3,
                       pVecTol=None,
                       ):
    """
    see if grains in grainList0 show up in spots;
    meant to be useful for taking existing grains from a load step
    and looking for them in a new load step;
    returns a new list of the same length, with None wherever a
    grain near the existing one was not found, and a list of grains for
    each path in pathList; grains in
    grainList0 are not modified;
    spots should be a Spots instance
    """
    # this import here because of circular reference
    from hexrd.xrd import grain

    newGrainList = [None for iGrain in range(len(grainList))]

    interactive = debug > 1

    'check that have both eosDict and refPDDict'
    assert (eosDict is None) == (refPDDict is None), \
        'need both or neither of eosDict and refPPDict'

    for iGrain, oldGrain in enumerate(grainList):
        if oldGrain is None:
            'do this check in case not culling lost grains so that numbering stays the same'
            continue
        newGrain = oldGrain.newGrain(spots,
                                     lineage = (iGrain),
                                     findByPixelDist=findByPixelDist)
        if debug:
            print 'completeness of candidate grain is : %g' % \
                (newGrain.completeness)
        if minCompleteness is not None:
            """
            # make this default to something like 0.5?
            # better to check this after fit?
            """
            if newGrain.completeness < minCompleteness:
                newGrain = None
        if newGrain:
            if debug:
                print 'doing claims for new grain %d' % (iGrain)
            newGrain.claimSpots()
        newGrainList[iGrain] = newGrain

    grainLists = [newGrainList]

    "now work on pathList"
    for pathName, rMatTransfList, phaseIDTo, phaseIDFrom, pRefSvecList in pathList:
        pathGrainList = []
        phaseTrans = False
        if phaseIDTo is not phaseIDFrom:
            phaseTrans = True
            if eosDict is not None:
                eosFr = eosDict[phaseIDFrom]
                eosTo = eosDict[phaseIDTo]
                pdFr  = refPDDict[phaseIDFrom]
                pdTo  = refPDDict[phaseIDTo]
        assert minCompleteness is not None,\
            'if specify paths, also specify minCompleteness'

        """
        use newGrainList here in case did fitting, and in case lost anyone,
        but still call them oldGrain
        """
        for iGrain, oldGrain in enumerate(newGrainList): # grainList
            if oldGrain is None:
                continue
            lp = None
            if oldGrain.phaseID:
                assert phaseIDFrom is not None,\
                    'grains have phases specified, need to provide a non-None "from" phase in the path'
                if oldGrain.phaseID is not phaseIDFrom:
                    "this grain is not of the right phase for the given path"
                    continue
            if phaseTrans:
                assert oldGrain.phaseID, 'for phase change need grains to have a phaseID!'

                vFr      = oldGrain.vol / pdFr.latVecOps['vol']
                p,   eFr = eosFr.pe_of_vt(vFr, tK)
                vTo, eTo = eosTo.ve_of_pt(p, tK)
                if debug:
                    print 'pressure in %s phase parent grain determined to be : %e' % \
                        (oldGrain.phaseID, p)
                try:
                    lpTag = pdTo.getLatticeType()
                    lpAtP = eosTo.lparms_of_ve(vTo, eTo)
                    # lp needs to be for 'triclinic', so need to convert
                    lp = crystallography.getDparms(lpAtP, lpTag)
                except NotImplementedError:
                    print >> sys.stderr, "WARNING: lparms_of_ve not implemented"
                    # lp = None # already done above
            for iRMat, rMatTransf in enumerate(rMatTransfList):
                lineage = (iGrain, pathName, iRMat)
                newGrain = oldGrain.newGrain(spots,
                                             phaseID=phaseIDTo,
                                             rMatTransf=rMatTransf,
                                             lineage=lineage,
                                             claimingSpots=False,
                                             lp=lp)
                if debug:
                    print 'for iRMat %d completeness is : %g' % \
                        (iRMat, newGrain.completeness)
                if newGrain.completeness < minCompleteness:
                    newGrain = None
                if newGrain is not None:
                    if debug:
                        print 'doing claims for new grain %s' % (str(lineage))
                    newGrain.claimSpots()
                    pathGrainList.append(newGrain)
        'end of newGrainList loop'
        grainLists += [pathGrainList]
    'end of pathList loop'

    if doFitting:

        lenSpotsPrev = len(spots)

        for iRefitLoop in range(maxIterRefit):
            'allow some looping'

            if debug:
                print 'about to call fitSpots with claimsBased=True, currently have %d spots ...' % (len(spots))
            spots.fitSpots(None, claimsBased=True, funcType='gaussGenEll', interactive=interactive)
            if debug:
                print ' ... back from fitSpots, now have %d spots' % (len(spots))
            spots.resetClaims()

            nConfTot = 0
            for iList, grainList in enumerate(grainLists):
                for iGrain, grain in enumerate(grainList):

                    """
                    redo findMatches now that have looked at spots in more detail
                    do not need doFit=True because did fits above
                    """
                    grain.findMatches(findByPixelDist=findByPixelDist, updateSelf=True)

                    # 'ignore claims so that can check out spots that may be overlapped'
                    # validSpotIdx, hitReflId = grain.getValidSpotIdx(ignoreClaims=True)

                    d_resRel = {}
                    d_integI = {}
                    d_spotsl = {}
                    validSpotIdx, hitReflId = grain.getValidSpotIdx(ignoreClaims=True) # ignore claims so that can check out spots that may be overlapped
                    iHKLs = grain.grainSpots['iHKL'][hitReflId]
                    for spotIdx, iHKL in zip(validSpotIdx, iHKLs):
                        spot = spots[spotIdx]
                        integIntens, res, res0  = spots.getIntegratedIntensity(spotIdx, useFit=True)
                        resRelList = d_resRel.setdefault( iHKL, [])
                        integIList = d_integI.setdefault( iHKL, [])
                        spotsList  = d_spotsl.setdefault( iHKL, [])
                        resRelList.append( (res,res0) )
                        integIList.append(integList)
                        spotsList.append(spot)
                    for iHKL in d_resRel.keys():
                        resRelList = d_resRel[iHKL]
                        res, res0 = zip(*resRelList)
                        resRel = num.array(res)/num.array(res0)
                        rrMean = resRel.mean()
                        rrMax  = num.max(resRel - rrMean)
                        rrStd  = resRel.std()
                        integI = d_integI[iHKL]
                        iiMean = integI.mean()
                        iiMax  = num.max(num.abs(integI - iiMean))
                        iiStd  = integI.std()
                        if rrMax > rrStd * 1.5 or iiMax > iiStd * 1.5:
                            if interactive or debug:
                                print 'for iHKL %d, rr and ii numbers are : %g %g %g ; %g %g %g ' % (
                                    iHKL,
                                    rrMax, rrMean, rrStd,
                                    iiMax, iiMean, iiStd,
                                    )
                            if interactive:
                                doLoop = True
                                while doLoop:
                                    userInput = raw_input('options: (d)isplay, (c)ontinue, (i:#)nvalidate spot #')
                                    char = userInput.split(':')[0]
                                    if   char == 'd':
                                        for iSpot, spot in enumerate(d_spotsl[iHKL]):
                                            spot.display(title='spot %d for iHKL %d' % (iSpot, iHKL))
                                    elif char == 'c':
                                        doLoop = False
                                    elif char == 'i':
                                        iSpot = userInput.split(':')[1]
                                        spot = d_spotsl[iHKL][iSpot]
                                        self.cleanFit()
                                        self.__fitFailed = FitFailedError(1, 'error from  user interaction')
                                        doLoop = True
                                    else:
                                        print 'invalid choice'

                    # ...*** may want to add logic to grain so that it knows when a repeated fit is necessary -- would have to have it ask spots about whether or not they have changed? might get a bit nasty
                    if indepFitPVec:
                        'fit precession and deformation gradient separately'
                        grain.fitPrecession()
                        grain.fit(fitPVec=False)
                    else:
                        grain.fit()
                    'redo findMatches now that did fit'
                    grain.findMatches(findByPixelDist=findByPixelDist, updateSelf=True)

                    if len(grain.lineageList) > 0:
                        lineage = grain.lineageList[-1]
                        prevGrain = grainList[lineage[0]]
                        if pVecTol is not None:
                            pVecDist = num.linalg.norm(grain.pVec - prevGrain.pVec)
                            if pVecDist > pVecTol:
                                print 'grain %d in list %d : have pVec mismatch with previous: %s versus %s' % (
                                    iGrain, iList, str(grain.pVec), str(prevGrain.pVec),
                                    )
                        if len(lineage) == 1:
                            'normal parentage'
                            validSpotIdx_new, hitReflId_new = grain.getValidSpotIdx(ignoreClaims=True)
                            validSpotIdx_prv, hitReflId_prv = prevgrain.getValidSpotIdx(ignoreClaims=True)
                            lostSpots = num.setdiff1d(hitReflId_prv, hitReflId_new, assume_unique=True)
                            newSpots  = num.setdiff1d(hitReflId_new, hitReflId_prv, assume_unique=True)
                            nLost = len(lostSpots); nNew = len(newSpots);
                            if nLost > 0 or nNew > 0:
                                print 'grain %d in list %d : have %d lost and %d new spots; new spots : ' % (
                                    iGrain, iList, nLost, nNew, str(newSpots),
                                    )
                        elif len(lineage) == 3:
                            'path parentage; nothing more to do'
                            pass
                        else:
                            raise RuntimeError, 'bad lineage : '+str(lineage)

                    grain.claimSpots()
                    #
                    nAll, nHit, nConf, nFail, nMiss = spot.getSpotCounts()
                    nConfTot += nConf

                'end of grainList loop'
            'end of grainLists loop'

            if len(spots) == lenSpotsPrev and nConfTot == 0:
                'assume there is no need to go through the loop again'
                if debug:
                    print 'exiting refit loop on iteration number %d because len(spots) did not change' % (iRefitLoop+1)
                break
            else:
                if debug:
                    print 'on iteration %d, len(spots) was %d and is now %d, and grains have %d conflicts' % (
                        iRefitLoop+1, lenSpotsPrev, len(spots), nConfTot,
                        )

        'end of iRefitLoop'

    if len(pathList) == 0:
        assert len(grainLists) == 1, 'expecting len(grainLists) == 1'
        grainLists = grainLists[0]
    return grainLists

def stretchToLV(V, fMat):
    """
    from stretch V in crystal frame and fMat, compute
    new lattice parameters; fMat can be the 'F' from
    lparms.latticeVectors

    V = B + I, where B is the Biot strain
    """
    # Finit = latticeVectors(initLP, 'triclinic' )['F']
    Fnew  = num.dot(V, fMat)
    lp    = latticeParameters(Fnew)
    lv    = latticeVectors(lp, tag='triclinic', radians=True)

    return lv

def calculateBiotStrain(initLP, finalLP, display=False):
    """
    """
    Finit = latticeVectors(initLP, 'triclinic' )['F']
    Fref  = latticeVectors(finalLP, 'triclinic' )['F']

    mF = num.zeros((3, 3))

    # diagonals are straight stretches
    mF[0, 0] = Fref[0, 0] / Finit[0, 0]
    mF[1, 1] = Fref[1, 1] / Finit[1, 1]
    mF[2, 2] = Fref[2, 2] / Finit[2, 2]

    # upper triangle
    mF[0, 1] = (Fref[0, 1] - mF[0, 0]*Finit[0, 1]) / Finit[1, 1]
    mF[0, 2] = (Fref[0, 2] - mF[0, 0]*Finit[0, 2] - mF[0, 1]*Finit[1, 2]) / Finit[2, 2]
    mF[1, 2] = (Fref[1, 2] - mF[1, 1]*Finit[1, 2]) / Finit[2, 2]

    u, s, v = svd(mF)
    biotStrain = num.dot(v.T, num.dot(num.diag(s), v)) - num.eye(3)
    if display:
        print "refined biot strain tensor: \n" + \
              "%1.3e, %1.3e, %1.3e\n" % (biotStrain[0, 0], biotStrain[0, 1], biotStrain[0, 2]) + \
              "%1.3e, %1.3e, %1.3e\n" % (biotStrain[1, 0], biotStrain[1, 1], biotStrain[1, 2]) + \
              "%1.3e, %1.3e, %1.3e\n" % (biotStrain[2, 0], biotStrain[2, 1], biotStrain[2, 2])
    return biotStrain, u, s, v

def makeMNConn(m,n,tri=True):
    """
    m and n are number of edges, so 1 more than number of zones
    """
    c1 = num.zeros( (m - 1)*(n - 1) , dtype=int)
    c2 = num.zeros( (m - 1)*(n - 1) , dtype=int)
    c3 = num.zeros( (m - 1)*(n - 1) , dtype=int)
    c4 = num.zeros( (m - 1)*(n - 1) , dtype=int)
    ii = 0
    for i in num.arange(1, m):
        for j in num.arange(1, n):
            c1[ii] = j + n*(i - 1)
            c2[ii] = j + n*(i - 1) + 1
            c3[ii] = j + n*i + 1
            c4[ii] = j + n*i

            ii += 1
            pass
        pass
    #
    'make quadCon be 0 based'
    quadCon = num.vstack([c1, c2, c3, c4]) - 1
    if tri:
        triCon1 = quadCon[:3, :]
        triCon2 = quadCon[[0, 2, 3], :]
        triCon  = num.hstack((triCon1, triCon2))
        retval = triCon
    else:
        retval = quadCon
    return retval

def makeNVecs(tth, etaIn, omeIn, asGrid=False):
    if asGrid:
        eta, ome = num.meshgrid(etaIn, omeIn)
    else:
        eta, ome = etaIn, omeIn
    qVecs = makeMeasuredScatteringVectors(
        num.asarray(tth).flatten(),
        num.asarray(eta).flatten(),
        num.asarray(ome).flatten())
    # qVecs is 3xN
    nVecs  = matrixutil.unitVector(qVecs)
    return nVecs

def omeEtaGridToNVecs(tTh, omegas, etas):
    etaGrid, omeGrid = num.meshgrid(etas, omegas)
    tThGrid = num.tile(tTh, (1, omeGrid.size) )
    nVecs = makeNVecs(tThGrid, etaGrid, omeGrid)
    return nVecs

rMat_y90 = ors.RotInv(num.pi/2.0, [0,1,0]).toMatrix()
def roty90(v):
    v_r = num.dot(rMat_y90, v)
    return v_r

class OmeEtaPfig(object):
    """
    works with data from with CollapseOmeEta
    """
    def __init__(self,
                 omeEdges, etaEdges,
                 cmap=None,
                 vMM=None,
                 doRot = False,
                 invertFromSouthern = True,
                 netStyle = None,
                 netNDiv = 12,
                 netAlpha = 0.5,
                 pointLists = [],
                 drawColorbar = True,
                 ):
        self.cmap = cmap
        if self.cmap is None:
            self.cmap = getCMap(False)

        if vMM is not None:
            self.vmin = vMM[0]; self.vmax = vMM[1];
        else:
            self.vmin = 0; self.vmax = None
            pass

        if isinstance(omeEdges, list): # hasattr(omeEdges,'pop'):
            self.omeEdges = omeEdges
        else:
            self.omeEdges = [omeEdges]
        if isinstance(etaEdges, list): # hasattr(etaEdges,'pop'):
            self.etaEdges = etaEdges
        else:
            self.etaEdges = [etaEdges]

        self.doRot = doRot
        self.netStyle = netStyle
        self.netNDiv = netNDiv
        self.netAlpha = netAlpha
        self.invertFromSouthern = invertFromSouthern
        self.pointLists = pointLists
        self.drawColorbar = drawColorbar

        self.antialiased = False
        self.pList = []
        self.colList = []
        self.cbar = None
        self.p = None

        return

    def __nVecsPatches(self, nP):
        'note: must test for bool first, as bools look like ints'
        if isinstance(nP,bool):
            retval = True
        elif isinstance(nP,int):
            retval = False
        else:
            raise RuntimeError, 'do not know what to do with nP of '+str(nP)
        return retval
    def __display(self, omeEdges, etaEdges, data, nVecs, nP, opacity, rangeVV_w):

        if self.allNorthern:
            """all points above the equator, so safe to use QuadMesh
            without worry of things being seen from the back side
            """
            showedges = False
            if not self.__nVecsPatches(nP):
                'render to a regular grid'
                vals = data.flatten()
                pfigR = pfigUtil.renderEAProj(nVecs, vals, nP)
                if opacity is not None:
                    norm = matplotlib.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
                    wData = pfigUtil.renderEAProj(nVecs, opacity.flatten(), nP)
                    vMM_w = Reader.getVMM(wData, range=rangeVV_w)
                    norm_w = matplotlib.colors.Normalize(vmin=vMM_w[0], vmax=vMM_w[1])
                    # pfigR.T here instead of ijAsXY=True
                    cData = self.cmap(norm(pfigR.T,clip=True))
                    cData[:,:,3] = norm_w(wData.T, clip=True)
                    im = self.pList[0](cData, ijAsXY=False)
                    self.colList.append(im)
                    'for colorbar, make a mappable so that the range shown is correct'
                    mappable = cm.ScalarMappable(cmap=self.cmap, norm=norm)
                    mappable.set_array(vals)
                    forColorBar = mappable
                else:
                    im = self.pList[0](pfigR, cmap=self.cmap, vmin=self.vmin, vmax=self.vmax, ijAsXY=True)
                    self.colList.append(im)
            else:
                'draw as a collection'
                if opacity is not None:
                    assert opacity.dtype == bool, 'for non-rendered pole figure, make opacities boolean'
                # if opacity is not None:
                #     raise RuntimeError, 'not coded: opacity for non-rendered pole figure, specify integer-valued nP'
                eaProj = pfigUtil.n2eap(nVecs, flip=False)
                if opacity is None:
                    col = collections.QuadMesh(len(etaEdges)-1,
                                               len(omeEdges)-1,
                                               eaProj,
                                               showedges,
                                               antialiased=self.antialiased)
                    col.set_array(data) # vals
                    # vals = (data[:,:].flatten()-vmin)/(vmax-vmin)
                    # vals = num.minimum(num.maximum(data[:,:].flatten(), vmin), vmax)
                    #
                else:
                    conn = makeMNConn(len(omeEdges), len(etaEdges), tri=False)
                    nQuads = conn.shape[1]
                    vals = data.flatten()
                    vertsN = []
                    valsN = []
                    for iQuad, quadConn in enumerate(conn.T):
                        if (opacity is not None) and (not opacity[iQuad]):
                            'quad is not to be drawn'
                            continue
                        nVecsPatch = nVecs[(slice(3), quadConn)] # 3x4
                        eaProjPatch = pfigUtil.n2eap(nVecsPatch, flip=False)
                        #verts[iQuad, :, :] = eaProjPatch.T
                        vertsN.append(eaProjPatch.T)
                        valsN.append(vals[iQuad])
                    colN = collections.PolyCollection(num.array(vertsN),
                                                      edgecolors='none',
                                                      antialiased=self.antialiased)
                    colN.set_array(num.array(valsN))
                    col = colN
                col.set_cmap(self.cmap)
                col.set_clim(vmin=self.vmin, vmax=self.vmax)
                #
                self.colList.append(col)
                self.pList[0].a.add_collection(col)

        else:
            """
            take a bit more care, due to patches ending up on the
            back side
            """
            if not self.__nVecsPatches(nP):
                'render to a regular grid'
                vals = data.flatten()
                northern = nVecs[2,:] >= 0
                southern = nVecs[2,:] <= 0
                nVecsN = nVecs[(slice(3), northern)]
                nVecsS = nVecs[(slice(3), southern)]
                if self.invertFromSouthern:
                    'invert through origin'
                    nVecsS = -nVecsS
                else:
                    'rotate about vertical axis in plane of pole figure'
                    nVecsS = num.vstack((-nVecsS[0,:], nVecsS[1,:], -nVecsS[2,:]))
                #
                if opacity is not None:

                    opacity = opacity.flatten()
                    norm = matplotlib.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
                    wDataN = pfigUtil.renderEAProj(nVecsN, opacity[northern], nP)
                    wDataS = pfigUtil.renderEAProj(nVecsS, opacity[southern], nP)
                    vMM_w = Reader.getVMM(num.vstack((wDataN,wDataS)), range=rangeVV_w)
                    norm_w = matplotlib.colors.Normalize(vmin=vMM_w[0], vmax=vMM_w[1])
                    #
                    pfigR = pfigUtil.renderEAProj(nVecsN, vals[northern], nPw)
                    cData = self.cmap(norm(pfigR.T,clip=True))
                    cData[:,:,3] = norm_w(wDataN.T, clip=True)
                    im = self.pList[0](cData, ijAsXY=False)
                    self.colList.append(im)
                    #
                    pfigR = pfigUtil.renderEAProj(nVecsS, vals[southern], nP)
                    cData = self.cmap(norm(pfigR.T,clip=True))
                    cData[:,:,3] = norm_w(wDataS.T, clip=True)
                    im = self.pList[1](cData, ijAsXY=False)
                    self.colList.append(im)
                    #
                    'for colorbar, make a mappable so that the range shown is correct'
                    mappable = cm.ScalarMappable(cmap=self.cmap, norm=norm)
                    mappable.set_array(vals)
                    forColorBar = mappable

                else:
                    pfigR = pfigUtil.renderEAProj(nVecsN, vals[northern], nP)
                    im = self.pList[0](pfigR, cmap=self.cmap, vmin=self.vmin, vmax=self.vmax, ijAsXY=True)
                    self.colList.append(im)
                    #
                    pfigR = pfigUtil.renderEAProj(nVecsS, vals[southern], nP)
                    im = self.pList[1](pfigR, cmap=self.cmap, vmin=self.vmin, vmax=self.vmax, ijAsXY=True)
                    self.colList.append(im)
            else:
                'draw as a collection'
                if opacity is not None:
                    assert opacity.dtype == bool, 'for non-rendered pole figure, make opacities boolean'
                # if opacity is not None:
                #     raise RuntimeError, 'not coded: opacity for non-rendered pole figure, specify integer-valued nP'
                conn = makeMNConn(len(omeEdges), len(etaEdges), tri=False)
                nQuads = conn.shape[1]
                #nVecsPatches = num.empty([nQuads, 4, 3])
                #verts = num.empty([nQuads, 4, 2])
                #vals = num.minimum(num.maximum(data[:,:].flatten(), vmin), vmax)# handled with set_clim
                vals = data.flatten()
                vertsN = []
                vertsS = []
                valsN = []
                valsS = []
                for iQuad, quadConn in enumerate(conn.T):
                    if (opacity is not None) and (not opacity[iQuad]):
                        'quad is not to be drawn'
                        continue
                    nVecsPatch = nVecs[(slice(3), quadConn)] # 3x4
                    if num.all(nVecsPatch[2,:] < 0):
                        'put in northern hemisphere'
                        if self.invertFromSouthern:
                            'invert through origin'
                            nVecsPatch = -nVecsPatch
                        else:
                            'rotate about vertical axis in plane of pole figure'
                            nVecsPatch = num.vstack((-nVecsPatch[0,:], nVecsPatch[1,:], -nVecsPatch[2,:]))
                        eaProjPatch = pfigUtil.n2eap(nVecsPatch, flip=False)
                        #verts[iQuad, :, :] = eaProjPatch.T
                        vertsS.append(eaProjPatch.T)
                        valsS.append(vals[iQuad])
                    else:
                        eaProjPatch = pfigUtil.n2eap(nVecsPatch, flip=False)
                        #verts[iQuad, :, :] = eaProjPatch.T
                        vertsN.append(eaProjPatch.T)
                        valsN.append(vals[iQuad])

                #col = collections.PolyCollection(verts, edgecolors=None, antialiased=antialiased)
                colN = collections.PolyCollection(num.array(vertsN),
                                                  edgecolors='none',
                                                  antialiased=self.antialiased)
                colN.set_clim(vmin=self.vmin, vmax=self.vmax)
                colN.set_array(num.array(valsN))
                colN.set_cmap(self.cmap)
                self.colList.append(colN)
                #
                colS = collections.PolyCollection(num.array(vertsS),
                                                  edgecolors='none',
                                                  antialiased=self.antialiased)
                colS.set_clim(vmin=self.vmin, vmax=self.vmax)
                colS.set_array(num.array(valsS))
                colS.set_cmap(self.cmap)
                self.colList.append(colS)

                self.pList[0].a.add_collection(colN)
                self.pList[1].a.add_collection(colS)
            'end of if for nP'
        'endif for all in upper hemisphere'
        return

    def display(self,
                data,
                tTh,
                nP=True,
                opacity=[None],
                rangeVV_w = None,
                winKWArgs={}
                ):
        """
        if want to do interpolated results, matplotlib.delaunay
        has some triangulation and interpolation stuff, but have
        not sorted that out yet
        """
        # bkgCirc = matplotlib.patches.Circle(xy=(0,0), radius=1, fill=True,
        #                                     facecolor='0.5',
        #                                     edgecolor='black')
        # pw.a.add_patch(bkgCirc)

        'clean up any old stuff'
        for p in self.pList:
            p.destroy()
        self.pList = []

        if self.doRot is True:
            'rotate to get "x" into "z" so that projection along 3 direction works well'
            rMat = rMat_y90
        elif hasattr(self.doRot, 'shape'):
            'assume doRot is the rMat to use'
            rMat = self.doRot
            assert len(rMat.shape) == 2, 'bad doRot'
            assert rMat.shape[0] == 3 and rMat.shape[1] == 3, 'bad doRot'
        elif self.doRot is False:
            rMat = None
        else:
            raise RuntimeError, 'do not know what to do with doRot of '+str(self.doRot)

        forColorBar = None

        if isinstance(data, list): # hasattr(data,'pop'):
            dataList = data
        else:
            dataList = [data]
        assert len(dataList) == len(self.omeEdges), \
            'data and omeEdges lists do not have same length : %d versus %d' % (len(dataList), len(self.omeEdges))
        #
        if isinstance(opacity, list): # hasattr(opacity,'pop'):
            opacityList = opacity
        else:
            opacityList = [opacity]
            assert len(opacityList) == len(self.omeEdges), \
            'opacity and omeEdges lists do not have same length : %d versus %d' % (len(opacityList), len(self.omeEdges))

        if self.vmax is None:
            self.vmax = num.hstack(dataList).max()

        nVecsList = []
        anySouthern = False
        for omeEdges, etaEdges in zip(self.omeEdges, self.etaEdges):
            if self.__nVecsPatches(nP):
                nVecs = omeEtaGridToNVecs(tTh, omeEdges, etaEdges)
            else:
                omegas = 0.5*(omeEdges[1:] + omeEdges[:-1])
                etas = 0.5*(etaEdges[1:] + etaEdges[:-1])
                nVecs = omeEtaGridToNVecs(tTh, omegas, etas)
            if rMat is not None:
                nVecs = num.dot(rMat, nVecs)
            anySouthern = anySouthern or num.any(nVecs[2,:] < -1e-5)
            nVecsList.append(nVecs)

        self.allNorthern = not anySouthern
        'need a window, make a new one'
        pw = None
        if self.drawColorbar:
            relfigsize = (6,6.3)
            if self.allNorthern:
                axesRectList = [
                    (0.15,0.15,0.6,0.6),
                    (0.85,0.15,0.05,0.75),
                    ]
            else:
                axesRectList = [
                    (0.03,0.1,0.4,0.8),
                    (0.43,0.1,0.4,0.8),
                    (0.85,0.15,0.05,0.75),
                    ]
        else:
            relfigsize = (6,6)
            if self.allNorthern:
                axesRectList = [(0.15,0.15,0.6,0.6)]
            else:
                axesRectList = [
                    (0.075,0.1,0.4,0.8),
                    (0.525,0.1,0.4,0.8),
                    ]
        kwArgs = {}
        kwArgs.update(winKWArgs)
        kwArgs.setdefault('relfigsize',relfigsize)
        kwArgs.setdefault('noAutoAxesList',True)
        win = plotWrap.PlotWin(1,2, **kwArgs)
        self.pList.append(plotWrap.PlotWrap(window=win,
                                        axesRect=axesRectList[0], # title='northern hemisphere'
                                        showByDefault=False, makeWinByDefault=True,
                                        )
                      )
        if not self.allNorthern:
            self.pList.append(plotWrap.PlotWrap(window=win,
                                                axesRect=axesRectList[1], # title='southern hemisphere'
                                                showByDefault=False, makeWinByDefault=True,
                                                )
                              )
        self.pList.append(win)
        self.p = win

        'now have pList'
        for pw in self.pList:
            if hasattr(pw,'getAxes'):
                'plotWin type thing instead of plotWrap type thing'
                continue
            pw.a.set_aspect('equal')
            pw.a.format_coord = lambda x,y:''
            pw.a.axis('off')
            pw.a.set_autoscale_on(False)
            pw.a.set_xlim(-1.05,1.05)
            pw.a.set_ylim(-1.05,1.05)
        'done with fun stuff for self.pList'

        self.colList = []

        for omeEdges, etaEdges, data, opacity, nVecs in zip(self.omeEdges, self.etaEdges, dataList, opacityList, nVecsList):
            self.__display(omeEdges, etaEdges, data, nVecs, nP, opacity, rangeVV_w)

        self.drawLines(nP, rMat=rMat)

        if self.drawColorbar:
            aCur = self.p.f.add_axes(axesRectList[-1])
            if forColorBar is None:
                forColorBar = self.colList[0]
            self.cbar = self.p.f.colorbar(forColorBar, cax=aCur)
        self.p.show()

        return
    def destroy(self):
        for p in self.pList:
            p.destroy()
        return
    def clearLines(self):
        'get rid of any existing lines and legends'
        for pw in self.pList:
            if hasattr(pw,'getAxes'):
                'plotWin type thing instead of plotWrap type thing'
                continue
            lines = pw.a.get_lines()
            for line in lines:
                line.remove()
        return
    def drawLines(self, nP, rMat=None, pointLists=None):
        """if specify a pointLists argument then it replaces self.pointLists"""
        if pointLists is not None:
            self.pointLists = pointLists

        # for pw in self.pList:
        #     if hasattr(pw,'getAxes'):
        #         'plotWin type thing instead of plotWrap type thing'
        #         continue

        auxKWArgs = {}
        if not self.__nVecsPatches(nP):
            r = 0.5 * float(nP)
            auxKWArgs = {'origin':(r,r), 'r':r}
        pw = self.pList[0]
        pfigUtil.drawLines(pw, self.pointLists,
                           rMat=rMat,
                           netStyle=self.netStyle, netNDiv=self.netNDiv, netAlpha=self.netAlpha,
                           **auxKWArgs)
        if not self.allNorthern:
            pw = self.pList[1]
            pfigUtil.drawLines(pw, self.pointLists,
                               rMat=rMat,
                               netStyle=self.netStyle, netNDiv=self.netNDiv, netAlpha=self.netAlpha,
                               southern=True, invertFromSouthern=self.invertFromSouthern,
                               **auxKWArgs)
        return

    def show(self):
        self.p.show()
        return
    def setProps(self, vMM=None,
                 cmap=None, pointLists=None, netStyle=None,
                 title=None,
                 doShow=True,
                 kwArgs={},
                 ):
        didChange = False
        if vMM is not None:
            self.vmin = vMM[0] or self.vmin
            self.vmax = vMM[1] or self.vmax
            for col in self.colList:
                col.set_clim(vmin=self.vmin, vmax=self.vmax)
            didChange = True
        if title is not None:
            self.p.f.suptitle(title, **kwArgs)
            didChange = True
        if cmap is not None:
            self.cmap = cmap
            for col in self.colList:
                col.set_cmap(self.cmap)
            didChange = True
        if pointLists is not None or netStyle is not None:
            if pointLists is not None:
                self.pointLists = pointLists
            if netStyle is not None:
                self.netStyle = netStyle
            self.drawLines()
            didChange = True
        if doShow and didChange:
            self.p.show()
        return

    def save(self, *args, **kwargs):
        self.p.save(*args, **kwargs)

class CollapseOmeEta(object):
    """
    Can pass a mask to use in addition to whatever the readers are already
    set up to do; with frames set zero where mask is True

    If pass pw: after plot is made, can do calls on self.pw to plot other things
    on top of the image

    while this takes a list of readers, it is not yet really coded for disjoint omega ranges
    """
    def __init__(self,
                 readerList,
                 planeData, hklList,
                 detectorGeom,
                 strainMag=None,
                 nEtaBins = 480,
                 mask=None,
                 threshold=None,
                 applyLorenz=True,
                 nframesLump=1,
                 debug=debugDflt,
                 computeMeanTwoTheta=False,
                 ):

        location = '  '+self.__class__.__name__+'.__init__'

        self.planeData = planeData

        if not isinstance(readerList, list): # hasattr(readerList,'pop'):
            readerList = [readerList]
        self.readerList = readerList

        self.iHKLList = num.atleast_1d(planeData.getHKLID(hklList))

        tThRanges = planeData.getTThRanges(strainMag=strainMag)

        'set up data'
        twoTheta, eta = detectorGeom.xyoToAngAll()
        # omegaMin, omegaMax = detector.getOmegaMMReaderList(readerList, overall=True) # retain information about omega range of each reader
        # omegaMM = detector.getOmegaMMReaderList(self.readerList) # can use this is need to know omega ranges for each reader

        'the following need to be consistent'
        etaDelta      = 2.*num.pi/nEtaBins
        eta0          = -num.pi
        self.etaEdges = num.linspace(-num.pi, num.pi, nEtaBins+1)
        self.etas = 0.5*(self.etaEdges[1:] + self.etaEdges[:-1])

        'count frames'
        nFramesTot = 0
        for reader in readerList:
            nFramesTot += reader.getNFrames() / nframesLump

        self.dataStore = num.empty([len(self.iHKLList), nFramesTot, nEtaBins])
        if computeMeanTwoTheta:
            self.storeTTh = num.empty([len(self.iHKLList), nFramesTot, nEtaBins])
        else:
            self.storeTTh = None
        self.omegas   = num.empty(nFramesTot)
        self.omeEdges = num.empty(nFramesTot+1)

        indicesList = []
        if computeMeanTwoTheta:
            'compute and store iEtaBinHKL'
            iEtaBinHKL = []
        for iData, iHKL in enumerate(self.iHKLList):
            indices = num.where(num.logical_and(
                    twoTheta > tThRanges[iHKL,0],
                    twoTheta < tThRanges[iHKL,1]
                    ))
            indicesList.append(indices)
            if computeMeanTwoTheta:
                etasThese = eta[indices]
                iEtaBin = (etasThese - eta0) / etaDelta
                'protect against numerical silliness'
                iEtaBin[num.where(iEtaBin) < 0         ] = 0
                iEtaBin[num.where(iEtaBin) > nEtaBins-1] = nEtaBins-1
                iEtaBinHKL.append(iEtaBin)


        if applyLorenz:
            lorenzFrame = getLorenz(detectorGeom)
        iFrameTot = 0
        for reader in readerList:
            delta_omega = reader.getDeltaOmega(nframes=nframesLump)
            sumImg = nframesLump > 1
            nFrames = reader.getNFrames()
            #
            for iFrame in range(reader.getNFrames() / nframesLump): # for iFrame, omega in enumerate(omegas):
                if debug:
                    print location+' : working on frame %d' % (iFrame)
                    tic = time.time()
                thisframe = reader.read(nframes=nframesLump, sumImg=sumImg)
                omega = reader.getFrameOmega()
                self.omegas[iFrameTot] = omega
                self.omeEdges[iFrameTot] = omega-delta_omega*0.5

                if not mask is None:
                    thisframe[mask] = 0

                if threshold is not None:
                    thisframe[thisframe < threshold] = 0

                for iData, iHKL in enumerate(self.iHKLList):
                    indices = indicesList[iData]
                    'make valsThese into float so that do not overrun max int value'
                    # valsThese = num.array(thisframe[indices], dtype=float) # problematic for masked arrays!
                    valsThese = thisframe[indices].astype(float)
                    if applyLorenz:
                        valsThese = valsThese / lorenzFrame[indices]
                    if num.ma.is_masked(valsThese):
                        keepThese = num.logical_not(valsThese.mask)
                        valsThese = valsThese.compressed()
                    else:
                        keepThese = None
                    if computeMeanTwoTheta:
                        ''' use scipy sparse matrices because
                        coo_matrix allows multiple values to be
                        additively assigned to an entry; surely there
                        is a better (and still easy) way to do this,
                        but I have not found it yet
                        '''
                        tThsThese = twoTheta[indices]
                        iEtaBin = iEtaBinHKL[iData]
                        if keepThese is not None:
                            iEtaBin   = iEtaBin[keepThese]
                            tThsThese = tThsThese[keepThese]
                        iRZ = num.zeros( valsThese.size, dtype=int )

                        integIntens_s = sparse.coo_matrix(
                            (valsThese, (iRZ, iEtaBin) ),
                            shape=(1,nEtaBins)
                            )
                        integIntens = integIntens_s.todense()
                        #
                        weightedTTh_s = sparse.coo_matrix(
                            (valsThese*tThsThese, (iRZ, iEtaBin) ),
                            shape=(1,nEtaBins)
                            )
                        weightedTTh = weightedTTh_s.todense()
                        #
                        nzIntegIntens = num.where(integIntens > 0.)
                        averageTTh = num.zeros( weightedTTh.shape )
                        averageTTh[nzIntegIntens] = weightedTTh[nzIntegIntens] / integIntens[nzIntegIntens]

                        self.dataStore[iData, iFrameTot, :] = integIntens[:]
                        self.storeTTh [iData, iFrameTot, :] = averageTTh[:]

                    else:
                        etasThese = eta[indices]
                        if keepThese is not None:
                            etasThese = etasThese[keepThese]
                        hist, bin_edges = num.histogram(etasThese, bins=self.etaEdges, weights=valsThese)
                        self.dataStore[iData, iFrameTot, :] = hist[:]
                'done with iHKLList'

                if debug:
                    toc = time.time()
                    elapsed = toc - tic
                    print location+' : frame %d took %g seconds' % (iFrameTot, elapsed)
                iFrameTot += 1
            'done with frames'

        'done with readerList'
        self.omeEdges[nFramesTot] = omega+delta_omega*0.5

        return
    def getPfig(self, **kwargs):
        retval = OmeEtaPfig(self.omeEdges, self.etaEdges, **kwargs)
        return retval
    def display(self,
                iData = 0,
                rangeVV = None,
                cmap = None,
                pw  = None,
                debug = debugDflt,
                pfig = False,
                comment = '',
                tTh = False,
                rangeVV_w = None,
                winKWArgs = {},
                pfigKWArgs = {}
                ):
        """
        iData is index into hklList passed on initialization, NOT into
        planeData

        pointList, if specified is a list of tuples of the form (3xn
        ndarray, style) where style gets passed to plotWrap

        can set tTh to 'withOpacity' if the 2-theta plot should be
        made with opacity (alpha) constructed from the intensity
        weights; so far this only works for pfig=False

        """

        if tTh:
            assert self.storeTTh is not None, 'do not have tThs stored'
            data = self.storeTTh [iData,:,:]
        else:
            data = self.dataStore[iData,:,:]

        # vmin, vmax = detector.Reader.getVMM(data)
        vMM = self.readerList[0].getVMM(data, range=rangeVV)

        hklStr = self.planeData.getHKLs(asStr=True)[self.iHKLList[iData]]
        winTitle = 'hkl: '+str(hklStr)+' '+comment
        if pfig:

            assert pw is None, \
                'do not support input pw for pfig plots'
            #retval = OmeEtaPfig(self.omeEdges, self.etaEdges, vMM=vMM, **pfigKWArgs)
            retval = self.getPfig(vMM=vMM, cmap=cmap, **pfigKWArgs)
            tThNominal = self.planeData.getTTh()[self.iHKLList[iData]]
            kwArgs = {}
            kwArgs.update(winKWArgs)
            kwArgs.setdefault('title', winTitle)
            if tTh == 'withOpacity':
                retval.display(data, tThNominal, nP=pfig,
                               opacity=self.dataStore[iData,:,:],
                               rangeVV_w=rangeVV_w,
                               winKWArgs=kwArgs)
            else:
                retval.display(data, tThNominal, nP=pfig, winKWArgs=kwArgs)

        else: # pfig
            if pw is None:
                pw = plotWrap.PlotWrap(
                    winTitle=winTitle,
                    showByDefault=False, makeWinByDefault=True,
                    **winKWArgs
                    )

                pw.a.format_coord = FormatEtaOme(0.5*(self.etaEdges[1:]+self.etaEdges[:-1]),
                                                 self.omegas,
                                                 data)
            retval = pw
            if cmap is None:
                cmap = getCMap(True)
            if tTh == 'withOpacity':
                norm = matplotlib.colors.Normalize(vmin=vMM[0], vmax=vMM[1])
                # data.T here instead of ijAsXY=True
                cData = cmap(norm(data.T,clip=True))

                wData = self.dataStore[iData,:,:]
                vMM_w = self.readerList[0].getVMM(wData, range=rangeVV_w)
                norm_w = matplotlib.colors.Normalize(vmin=vMM_w[0], vmax=vMM_w[1])

                cData[:,:,3] = norm_w(wData.T, clip=True)

                pw(cData, ijAsXY=False)
                'for colorbar, make a mappable so that the range shown is correct'
                mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
                mappable.set_array(data)
                pw.colorbar(thing=mappable)

            else:
                pw(data, vmin=vMM[0], vmax=vMM[1], cmap=cmap, ijAsXY=True)
            pw.show()
        'endif pfig'

        return retval
    def getNVecs(self, iData, edges=False):
        tThNominal = self.planeData.getTTh()[self.iHKLList[iData]]
        if edges:
            nVecs = omeEtaGridToNVecs(tThNominal, self.omeEdges, self.etaEdges)
        else:
            nVecs = omeEtaGridToNVecs(tThNominal, self.omegas, self.etas)
        return nVecs
    def getData(self, iData):
        return self.dataStore[iData,:,:]

def makeMeasuredScatteringVectors(tTh, eta, ome, convention='LLNL', frame='sample'):
    """
    Construct scattering vectors from angles (2theta, eta, ome)
    will do LLNL/APS and Risoe frames, sample or lab.

    for risoe frame geomtery, see:

    http://sourceforge.net/apps/trac/fable/attachment/wiki/WikiStart/Geometry_version_1.0.8.pdf
    """
    ome = num.atleast_1d(ome)

    rIndex = None
    if not hasattr(convention, 'lower'):
        raise RuntimeError, "something is amiss with your keywork argument"
    else:
        if convention.lower() == 'risoe':
            # must convert eta coordinate
            eta = 0.5*num.pi - eta

            'ome axis is lab Z'
            # raxis = num.c_[0, 0, 1].T
            rIndex = 2

            # the unit qVec
            retVal = num.vstack( [
                -num.sin(0.5*tTh),
                -num.cos(0.5*tTh)*num.sin(eta),
                 num.cos(0.5*tTh)*num.cos(eta) ] )

        elif convention.lower() == 'llnl':
            'ome axis is lab Y'
            # raxis = num.c_[0, 1, 0].T
            rIndex = 1

            # the unit qVec
            retVal = num.vstack( [
                num.cos(eta)*num.cos(0.5*tTh),
                num.sin(eta)*num.cos(0.5*tTh),
                num.sin(0.5*tTh) ] )
        else:
            raise RuntimeError, "unrecognized token '%s'" % (convention)

        # convert to sample frame from lab frame (as measured)
        if frame.lower() == 'sample':
            ## for i in range(retVal.shape[1]):
            ##     retVal[:, i] = num.dot(R.rotMatOfExpMap(ome[i]*raxis).T, retVal[:, i])
            cOme = num.cos(ome)
            sOme = num.sin(ome)
            if rIndex == 1:
                retVal = num.vstack( [
                        cOme * retVal[0,:] - sOme * retVal[2,:], # Note: -sOme due to transpose
                        retVal[1,:],
                        sOme * retVal[0,:] + cOme * retVal[2,:], # Note: +sOme due to transpose
                        ] )
            elif rIndex == 2:
                retVal = num.vstack( [
                        cOme  * retVal[0,:] + sOme * retVal[1,:], # Note: +sOme due to transpose
                        -sOme * retVal[0,:] + cOme * retVal[1,:], # Note: -sOme due to transpose
                        retVal[2,:],
                        ])
            else:
                raise RuntimeError, "unrecognized rIndex %d" % (rIndex)
    return retVal

######################################################################
#
# fun with multiprocessing

reader_MP = None
func_MP = None
nPerChunk_MP = None
debug_MP = True
def doItMP(nSkip):
    location = '  doItMP'

    global reader_MP, func_MP, debug_MP
    print 'in doItMP, nPerChunk nSkip : %g %g' % (nPerChunk_MP, nSkip)
    tic = time.time()
    readerLocal = reader_MP.makeNew()
    toc = time.time(); dt = toc - tic; tic = toc
    print location+' : %g seconds to make new reader' % (dt)
    retval = []
    nLocalFrame = min(nPerChunk_MP, readerLocal.getNFrames()-nSkip)
    for iLocalFrame in range(nLocalFrame):
        frame = readerLocal.read(nskip=nSkip)
        toc = time.time(); dt = toc - tic; tic = toc
        print location+' : %g seconds to read frame' % (dt)
        retval.append(func_MP(frame))
        toc = time.time(); dt = toc - tic; tic = toc
        print location+' : %g seconds to process frame' % (dt)
        nSkip = 0 # only need to skip for first read
    return retval
def readFrameStack_multiproc(reader, func, nPerChunk=None, nCPUs=None, debug=False):
    """
    read the full frame stack using multiprocessing, applying fund
    to each frame to obtain the result;

    use makeNew for each chunk; if reader.dark is a shared memory
    array then it remains so
    """
    assert xrdBase.haveMultiProc, \
        'do not have multiprocessing available'

    nFrames = reader.getNFrames()
    nCPUs = nCPUs or xrdBase.dfltNCPU
    nPerChunk = nPerChunk or max(int(nFrames/(4*nCPUs)),1)

    global reader_MP, func_MP, nPerChunk_MP, debug_MP
    reader_MP = reader
    func_MP = func
    nPerChunk_MP = nPerChunk
    debug_MP = debug
    'if reader has a dark frame, make sure it is in shared memory'
    if reader.dark is not None:
        #dark_temp  = detector.ReadGE.readDark(fileInfo[0][0], nframes=fileInfo[0][1])
        dark_temp   = reader.dark
        dark_shmem  = multiprocessing.RawArray(ctypes.c_double, dark_temp.size)
        dark        = num.frombuffer(dark_shmem, dtype=num.float64, count=dark_temp.size)
        dark.shape  = dark_temp.shape
        dark[:,:]   = dark_temp[:,:]
        reader.dark = dark
        del dark_temp

    nSkipList = range(0,nFrames,nPerChunk)

    pool = multiprocessing.Pool(nCPUs)

    if debug:
        print 'nSkipList : '+str(nSkipList)
    tic = time.time(); results_chunked = pool.map(doItMP, nSkipList, chunksize=1); toc = time.time()
    dt = toc - tic
    print 'elapsed time running under multiprocessing pool : %g' % (dt)
    results = reduce(lambda x,y: x+y, results_chunked, [])

    pool.close()

    reader_MP = None
    func_MP = None
    nPerChunk_MP = None
    return results

def thresholdStackDisplay(data, threshold, cmap=None, pw=None,
                          detectorGeom=None, planeData=None,
                          displayKWArgs={}, drawRingsKWArgs={}
                          ):
    '''
    passes sumImg=num.maximum to dataToFrame so that if data is a
    reader then frame ends up being the maximum over the image stack
    '''
    thisframe = dataToFrame(data, sumImg=num.maximum)
    if cmap is None:
        cmap = cm.bone
        cmap.set_over('red')
    if detectorGeom is None:
        retval = detector.Reader.display(
            thisframe, cmap=cmap, range=(0,threshold), pw=pw,
            **displayKWArgs)
    else:
        retval = detectorGeom.display(
            thisframe, planeData=planeData, cmap=cmap, range=(0,threshold), pw=pw,
            **displayKWArgs)
        detectorGeom.drawRings(retval, planeData, **drawRingsKWArgs)
    return retval

class GrainPoles:
    def __init__(self, omeEtaPfigs):

        self.omeEtaPfigs = omeEtaPfigs

        self.refPointsLists = []
        for omeEtaPfig in self.omeEtaPfigs:
            self.refPointsLists.append(copy.deepcopy(omeEtaPfig.pointLists))
            # refPointsList = []
            # for refPoints, kwArgs in omeEtaPfig.pointLists:
            #     refPointsList.append(copy.deepcopy(refPoints))
            # refPointsLists.append(refPointsList)

        self.__skew = [0.,0.,0.]

        return
    def __call__(self, skew):
        self.__skew = num.copy(skew)
        mat = self.getMat()
        for omeEtaPfig, refPointsList in zip(self.omeEtaPfigs, self.refPointsLists):
            # omeEtaPfig.clearLines() 'not necessary, drawLines does this'
            pointLists = []
            for refPoints, kwArgs in refPointsList:
                points    = num.dot(mat, refPoints)
                pointLists.append((points, kwArgs))
            # print 'pointLists : '+str(pointLists)
            omeEtaPfig.drawLines(pointLists=pointLists)
            # print 'pointLists : '+str(omeEtaPfig.pointLists)
            omeEtaPfig.show()
        return
    def getMat(self):
        mat  = ors.RotInv(self.__skew).toMatrix()
        return mat

def grainPolesGUI(omeEtaPfigs):
    """
    GUI with sliders for rotating a grain's spots

    execfile('examples/feSynthSpotsPfig.py')
    gui = xrdUtils.grainPolesGUI([pwSB])
    """
    win = plotWrap.PlotWin(-1,
                           title='rotation sliders',
                           relfigsize=(8,2),
                           )
    ax1 = win.getNextAxes(rect=[0.1, 0.1, 0.8, 0.2])
    ax2 = win.getNextAxes(rect=[0.1, 0.4, 0.8, 0.2])
    ax3 = win.getNextAxes(rect=[0.1, 0.7, 0.8, 0.2])
    s1  = Slider(ax1, 'r1', -pi, pi, 0., dragging=False)
    s2  = Slider(ax2, 'r2', -pi, pi, 0., dragging=False)
    s3  = Slider(ax3, 'r3', -pi, pi, 0., dragging=False)

    grainPoles = GrainPoles(omeEtaPfigs)

    def update_self(val):
        skew = [s1.val, s2.val, s3.val]
        grainPoles(skew)
        print 'angles are now : %24.16e %24.16e %24.16e' % \
            (skew[0], skew[1], skew[2], )
        # print 'rotation matrix is: '+str(mat)
        win.show()
    s1.on_changed(update_self)
    s2.on_changed(update_self)
    s3.on_changed(update_self)

    win.show()
    matplotlib.interactive(True)
    return win, grainPoles

class MultiSlopeFunc:
    """
    function that transitions smoothly from slope of m1 to slope of m2 in the vacinity of xcrit, with the smoothness of the transition dictated by a specified power;
    for large values of power, might eventually want to put in code to protect against numerical overflow
    """
    def __init__(self, m1, m2, xcrit, power=0.2):
        self.m1 = m1
        self.m2 = m2
        self.power = power
        self.power_inv = 1.0/power
        self.xcrit = xcrit
        return
    def _f1(self,x):
        retval = x * self.m1
        return retval
    def _f2(self,x):
        if x > self.xcrit:
            retval = (x - self.xcrit) * self.m2
        else:
            retval = 0.
        return retval
    def __call__(self,x):
        xa = math.fabs(x)
        f1 = self._f1(xa)
        f2 = self._f2(xa)
        retval = (f1**self.power_inv + f2**self.power_inv)**self.power
        if x < 0:
            retval = -retval
        return retval

class MultiSlopeFuncSmooth:
    """
    function that transitions smoothly from slope of m1 to slope of m2 in the vacinity of xcrit, with the smoothness of the transition dictated by a specified power;
    """
    def __init__(self, m1, w1, p1, m2, xcrit, power=0.2):
        self.m1 = m1
        self.w1 = w1
        self.p1 = p1
        self.m2 = m2
        self.power = power
        self.power_inv = 1.0/power
        self.xcrit = xcrit
        return
    def _f1(self,x):
        if x > self.xcrit:
            retval = ((x - self.xcrit)/self.w1)  ** self.p1
        else:
            retval = 0.
        retval += x * self.m1
        return retval
    def _f2(self,x):
        retval = x * self.m2
        return retval
    def __call__(self,x):
        xa = math.fabs(x)
        f1 = self._f1(xa)
        f2 = self._f2(xa)
        if f1 == 0 or f2 == 0:
            return 0.0
        retval = (f1**(-self.power_inv) + f2**(-self.power_inv))**(-self.power)
        if x < 0:
            retval = -retval
        return retval

def darkFromStack(reader, nFrames=0, nChunk=4,
                  medianSize=None, medianRange=(-15,15),
                  cutMinFactor=None,
                  checkIntensityResolution=False):
    """
    If nFrames == 0, read all frames.

    If medianSize is specified then a median filter of the given size is used to find
    dead pixels, with pixels outside of medianRange marked as dead.
    """
    nChop = 1
    nRow, nCol = reader.getSize()
    if nFrames == 0:
        nFrames = reader.getNFrames()
    if nChunk is not None:
        nChop = max(1, int(nFrames / nChunk))

    # would be perhaps more efficient to just read what is wanted each time
    darkFrame  = reader.frame()
    stdFrame   = num.empty(darkFrame.shape,dtype=float)
    if checkIntensityResolution:
        meanDFrame  = reader.frame() # num.empty(darkFrame.shape,dtype=float)
    deadPixels = num.zeros(darkFrame.shape,dtype=bool)
    iRU = 0
    for iChop in range(nChop):
        iRL = iRU
        iRU = int(nRow*(float(iChop+1)/float(nChop)))
        thisReader = reader.getRawReader()
        allThisChop = num.empty((iRU-iRL,nCol,nFrames), dtype=darkFrame.dtype)
        for iFrame in range(nFrames): # for iFrame, omega in enumerate(omegas):
            print >> sys.stdout, '+',
            thisframe = thisReader.read()
            allThisChop[:,:,iFrame] = thisframe[iRL:iRU,:]
        print >> sys.stdout, ''
        thisReader.close()
        if cutMinFactor is not None:
            mins = num.min(allThisChop, axis=2)
            # mask = allThisChop > num.tile(mins.flatten()*cutMinFactor, (allThisChop.shape[2],1)).T.reshape(allThisChop.shape)
            allThisChopMedian = num.empty((iRU-iRL,nCol), dtype=darkFrame.dtype)
            for iThis in range(allThisChop.shape[0]):
                for jThis in range(allThisChop.shape[1]):
                    pStack = allThisChop[iThis,jThis,:]
                    allThisChopMedian[iThis,jThis] = num.median(pStack[pStack < mins[iThis,jThis] * cutMinFactor])
            darkFrame[iRL:iRU,:] = allThisChopMedian[:,:]
        else:
            darkFrame[iRL:iRU,:] = num.median(allThisChop, axis=2)
        std = num.std(allThisChop.astype('float'),axis=2)
        stdFrame[iRL:iRU,:] = std
        deadPixels[iRL:iRU,:] = std == 0.
        if checkIntensityResolution:
            # num.unique does not support an axis argument
            for iRowChop in range(allThisChop.shape[0]):
                iRow = iRowChop + iRL
                for iCol in range(allThisChop.shape[1]):
                    pUnique = num.unique(allThisChop[iRowChop,iCol,:])
                    if len(pUnique) > 1:
                        # meanDFrame[iRow,iCol] = num.mean(pUnique[1:] - pUnique[:-1])
                        meanDFrame[iRow,iCol] = num.min(pUnique[1:] - pUnique[:-1])
                    else:
                        meanDFrame[iRow,iCol] = 0
            # perhaps mark as dead if meanDFrame > 1.75
        # whereDead = num.where(std == 0.)
    print >> sys.stdout, 'number of dead pixels with 0 std : %d' % (num.sum(deadPixels))

    if medianSize is not None:
        darkFrameMF = ndimage.median_filter(darkFrame, size=medianSize)
        darkDiff = darkFrame - darkFrameMF
        print >> sys.stdout, 'number of dead pixels with below median filter by %g : %d' % \
            (medianRange[0], num.sum(darkDiff < medianRange[0]))
        print >> sys.stdout, 'number of dead pixels with above median filter by %g : %d' % \
            (medianRange[1], num.sum(darkDiff > medianRange[1]))
        deadPixels = num.logical_or(deadPixels, num.logical_or(
                darkDiff < medianRange[0], darkDiff > medianRange[1]
                ) )

    if checkIntensityResolution:
        deadPixels = num.logical_or(deadPixels, meanDFrame > checkIntensityResolution)

    if checkIntensityResolution:
        retval = darkFrame, deadPixels, stdFrame, meanDFrame
    else:
        retval = darkFrame, deadPixels, stdFrame
    return retval

def tryFromShelf(shelfFileName, thingName):
    """
    try to pull the thing from the shelf and return None if it does not work
    """
    try:
        shelf = shelve.open(shelfFileName, 'r')
        retval = shelf[thingName]
        shelf.close()
    except:
        retval = None
    return retval

def putInShelf(shelfFileName, thingName, thing):
    shelf = shelve.open(shelfFileName, 'c')
    shelf[thingName] = thing
    shelf.close()
    return

def pullFromStack(reader, detectorGeom, tThMM, angWidth, angCen,
                  threshold=20,
                  distInAng=False, padSpot=True,
                  mask3D=None, exitOnFail=False
                  ):
    """
    angWidth is distance from angCen, so actually half-width
    """
    omeStep = reader.getDeltaOmega()
    etaStep = getEtaResolution(detectorGeom, angCen[0])
    dr_Eta, dr_Ome, dA = drEtaOme(angCen, etaStep, omeStep)

    tThMin, tThMax = tThMM
    tthPM = 0.5 * (tThMax - tThMin)
    etaPM = angWidth*etaStep/dr_Eta
    omePM = angWidth*omeStep/dr_Ome
    angPM = [tthPM, etaPM, omePM]
    print 'angWidth %g mapped to eta ome widths of %g %g' % (angWidth, etaPM, omePM)

    'use xyf instead of xyo to reinforce that is in frame instead of omega coordinates'
    xyfCen  = num.array(detectorGeom.angToXYO(*angCen, units='pixels')) # omega is in angles
    xyfCen[2] = reader.omegaToFrame(xyfCen[2], float=True)

    xyfBBox = detectorGeom.angToXYOBBox(angCen, angPM, units='pixels', reader=reader) # omega is in frame -- passed the reader
    xyfBBox_0 = num.array([sl[0] for sl in xyfBBox])

    pixelData = reader.readBBox(xyfBBox, raw=False)
    #
    # pw = detector.Reader.display(num.max(pixelData,axis=2), range=(0,400))

    if mask3D is not None:
        maskThis = mask3D[slice(*xyfBBox[0]), slice(*xyfBBox[1]), slice(*xyfBBox[2])]
        if hasattr(pixelData, 'mask'):
            mask = num.logical_or( pixelData.mask, maskThis )
        else:
            mask = maskThis
        pixelData = num.ma.masked_array(pixelData, mask=mask, hard_mask=True, copy=False)

    # ndimage.watershed_ift may be worth trying at some point
    bin = spotFinder.getBin(pixelData, threshold, padSpot)
    labelStructure = ndimage.generate_binary_structure(pixelData.ndim,1)
    labels, numSpots = ndimage.label(bin, labelStructure)
    # labels may include masked pixels, but left them in on purpose so that spots would not get sliced by rows of bad pixels
    # masked pixels should get excluded by getIndices as appropriate
    objs = ndimage.find_objects(labels)
    if distInAng:
        'find spot based on distance in angular space'
        com  = num.empty([numSpots,3])
        vSum = num.empty(numSpots)
        if numSpots <= 0:
            iSpot = -1
        else:
            for iSpot in range(numSpots):
                index = iSpot+1
                'this com is with respect to pixelData box, so add in xyfBBox_0'
                comThis, vSumThis = spotFinder.getImCOM(pixelData, labels, objs, index, floor=0, getVSum=True)
                vSum[iSpot] = vSumThis
                com[iSpot,:] = comThis + xyfBBox_0
            comTTh, comEta, comOme = detectorGeom.xyoToAng(com[:,0], com[:,1], reader.frameToOmega(com[:,2]))
            dTTh = (comTTh - angCen[0])/tthPM
            dEta = (comEta - angCen[1])/etaPM
            dOme = (comOme - angCen[2])/omePM
            dSq = dTTh*dTTh + dEta*dEta + dOme*dOme
            iSpot = dSq.argmin()
    else:
        'get spot in which xyfCen falls'
        # okay if this lands in a masked pixel -- still associated with labeled spot
        index = labels[tuple(num.round(xyfCen).astype(int) - xyfBBox_0)]
        if index < 1:
            if exitOnFail:
                raise RuntimeError, 'no spot found at '+str(xyfCen)
            else:
                iSpot = -1
        else:
            iSpot = index-1

    return iSpot, labels, objs, pixelData, xyfBBox

def spotFromStack(reader, detectorGeom, tThMM, angWidth, angCen, threshold,
                  fullBackground=False,
                  asFrame=False, exitOnFail=True,
                  distInAng=False,
                  debug=True,
                  ):
    '''
    if asFrame, then omegas come out as frame indices; note that
    angCen should still be specified as omega, not a frame index
    '''
    iSpot, labels, objs, pixelData, xyfBBox = pullFromStack(reader, detectorGeom, tThMM, angWidth, angCen,
                                                            threshold=threshold,
                                                            distInAng=distInAng, exitOnFail=exitOnFail)
    if iSpot < 0:
        retval = None
    else:
        xyfBBox_0 = num.array([sl[0] for sl in xyfBBox])

        # grab pixel data
        index = iSpot+1
        obj = objs[iSpot]
        # indices = list(num.where(labels[obj] == index))
        indices = list(spotFinder.getIndices(pixelData, labels, objs[index-1], index)) # convert tuple to list so that can modify it
        for iX in range(len(indices)):
            indices[iX] = indices[iX] + obj[iX].start
        xAll = []
        for iX, xL in enumerate(indices):
            xAll.append( xL + xyfBBox_0[iX] )
        if not asFrame:
            xAll[2] = reader.frameToOmega(xAll[2]) # convert from frame to omega
        vAll = pixelData[indices]
        # note that num.sum(vAll) may not agree with vSumThis above due to floor=0 setting
        if hasattr(vAll, 'mask'):
            numMasked = num.sum(vAll.mask)
            if numMasked > 0:
                raise RuntimeError, 'got %d masked pixels' % (numMasked)
            vAll =vAll.shrink_mask()
        data = xAll
        data.append(vAll)
        retval = data

        if fullBackground:
            'add in pixels that are inside angWidth and are below background'
            # indicesBkg = list(spotFinder.getIndices(pixelData, labels, None, 0)) # convert tuple to list so that can modify it
            indicesBkg = spotFinder.getIndices(pixelData, labels, None, 0)
            'no need to add in obj[iX].start -- did full pixelData'
            xAllBkg = []
            for iX, xL in enumerate(indicesBkg):
                xAllBkg.append( xL + xyfBBox_0[iX] )
            if not asFrame:
                xAllBkg[2] = reader.frameToOmega(xAllBkg[2]) # convert from frame to omega
            'could check to see which pixels are within angWidth of angCen, but they should more or less all be'
            vAllBkg = pixelData[indicesBkg]
            # note that num.sum(vAll) may not agree with vSumThis above due to floor=0 setting
            if hasattr(vAllBkg, 'mask'):
                numMasked = num.sum(vAllBkg.mask)
                if numMasked > 0:
                    raise RuntimeError, 'got %d masked pixels' % (numMasked)
                vAllBkg = vAllBkg.shrink_mask()
            if debug:
                print 'have %d foreground and %d background spots' % (len(vAll), len(vAllBkg))
            dataBkg = xAllBkg
            dataBkg.append(vAllBkg)
            retval = map(num.hstack, zip(data, dataBkg))
            pass
        pass
    # spot = spotFinder.Spot(key, reader.getDeltaOmega(), data=data)
    return retval

def collapse(vAll, eta, ome,
             nOme=None, nEta=None,
             weightedList=[], averagedList=[], auxList=[],
             debug=False,
             ):
    """
    Returns a sparse matrix, with zeros where no intensity falls;

    pass nEta and nOme to control the number of eta and omega bins,
    otherwise they are determined automatically, with the omega
    binning assuming that omegas fall on frames at regular intervals
    with no gaps

    for each entry in weightedList, also returns that data collapsed and
    weighted by vAll; similarly for averagedList
    """

    'make vAll into float data so that summing does not overrun int storage'
    vAll = num.asarray(vAll).astype(float)

    ome = num.asarray(ome)
    eta = num.asarray(eta)

    if nOme is None:
        omegas = num.unique(ome)
        dOme = num.min(omegas[1:]-omegas[:-1])
        nOme = len(omegas) # assume omega frames are adjacent -- some data on each frame
        omeLo = omegas.min()-0.5*dOme
    else:
        if not hasattr(nOme,'__len__'): # or len(nOme) == 1:
            dOme = (ome.max() - ome.min())/max(1., nOme-1.)
            omeLo = ome.min()-0.5*dOme
        elif len(nOme) == 3:
            nOme, omeMin, omeMax = nOme
            omeLo = omeMin
            if nOme is None:
                omegas = num.unique(ome)
                dOme = num.min(omegas[1:]-omegas[:-1])
                nOme = int(round((omeMax - omeMin)/dOme))
            else:
                dOme = (omeMax - omeMin)/nOme
        else:
            raise RuntimeError, 'nOme is funny'
    if debug:
        print 'ome n, d, lo : %d %g %g' % (nOme, dOme, omeLo)

    if nEta is None:
        nEta = int(0.5*num.sqrt(eta.size/nOme))
        dEta = (eta.max() - eta.min())/max(nEta-1., 1.)
        etaLo = eta.min()-0.5*dEta
    else:
        if not hasattr(nEta,'__len__'):
            dEta = (eta.max() - eta.min())/max(nEta-1., 1.)
            etaLo = eta.min()-0.5*dEta
        elif len(nEta) == 3:
            nEta, etaMin, etaMax = nEta
            assert nEta is not None, 'nEta is None'
            etaLo = etaMin
            dEta = (etaMax - etaMin)/nEta
        else:
            raise RuntimeError, 'nEta is funny'
    if debug:
        print 'eta n, d, lo : %d %g %g' % (nEta, dEta, etaLo)


    # use sparse.coo_matrix to do summing
    # could instead use num.histogram2d

    nX, nY = nEta, nOme
    xI = ((eta - etaLo)/dEta).astype(int)
    yJ = ((ome - omeLo)/dOme).astype(int)
    a = sparse.coo_matrix( (vAll, (xI, yJ)),  shape=(nX,nY) ) # overlapped entries get summed
    a = a.tocsr() # to sum duplicate entries

    etas = num.linspace(etaLo+0.5*dEta, etaLo+(nEta-1)*dEta+0.5*dEta, nEta)
    omes = num.linspace(omeLo+0.5*dOme, omeLo+(nOme-1)*dOme+0.5*dOme, nOme)

    retval = [a, etas, omes]

    for thisData in weightedList:
        thisWeighted = vAll * thisData
        b = sparse.coo_matrix( (thisWeighted, (xI, yJ)),  shape=(nX,nY) ) # overlapped entries get summed
        b = b.tocsr() # to sum duplicate entries
        retval.append(b)
    if len(averagedList)>0 :
        xNZ, yNZ = a.nonzero() # not the same as (xI,yJ) due to overlaps
    for thisData in averagedList:
        thisWeighted = vAll * thisData
        b = sparse.coo_matrix( (thisWeighted, (xI, yJ)),  shape=(nX,nY) ) # overlapped entries get summed
        b = b.tocsr() # to sum duplicate entries
        for x, y in zip(xNZ,yNZ):
            axy = a[x,y]
            if axy > 0:
                b[x,y] = b[x,y] / axy
            else:
                b[x,y] = 0.
        retval.append(b)
    for vAll, eta, ome in auxList:
        xI = ((eta - etaLo)/dEta).astype(int)
        yJ = ((ome - omeLo)/dOme).astype(int)
        b = sparse.coo_matrix( (vAll, (xI, yJ)),  shape=(nX,nY) ) # overlapped entries get summed
        b = b.tocsr()
        retval.append(b)

    return tuple(retval)

def displaySparse(a,
                  vmin = None, vmax = None,
                  cmap = None,
                  fmt = None,
                  markersize = 3,
                  colorUnder = None,
                  ijNZ = None,
                  **kwargs):

    cmap   = cmap or detector.getCMap(False)

    pw = kwargs.pop('pw',None)
    if pw is None:
        pwKWArgs = plotWrap.PlotWrap.popKeyArgs(kwargs)
        pw = plotWrap.PlotWrap(**pwKWArgs)
    retval = pw
    #
    if len(kwargs) > 0:
        raise RuntimeError, 'unparsed keyword arguments : '+str(kwargs.keys())

    # xNZ, yNZ, vNZ = sparse.extract.find(a)
    if ijNZ is None:
        xNZ, yNZ = a.nonzero()
    else:
        xNZ, yNZ = ijNZ

    if sparse.issparse(a):
        A = a.todense()
    else:
        A = a
    if colorUnder is not None:
        vmin = vmin or 0
        B = - 10. * max(num.max(num.abs(A)), vmax) * num.ones(A.shape, dtype=A.dtype)
        B[(xNZ,yNZ)] = num.array(A[(xNZ,yNZ)]).reshape(len(xNZ))
        cmap.set_under(colorUnder)
        pw(B, vmin=vmin, vmax=vmax, cmap=cmap)
    else:
        pw(A, vmin=vmin, vmax=vmax, cmap=cmap)
        Z = num.ones(A.shape, dtype=int)
        Z[xNZ,yNZ] = 0
        B = sparse.coo_matrix(Z)
        pw.a.spy(B, markersize=markersize)

    if fmt is not None:
        fmt.A = A
        pw.a.format_coord = fmt

    pw.show()

    return retval

def displayEtaOme(eta, ome, vAll,
                  nEta = None, nOme = None,
                  **kwargs):

    vCollapsed, etas, omes = collapse(vAll, eta, ome, nEta=nEta, nOme=nOme)
    fmt = FormatEtaOme(etas, omes, None)
    pw = displaySparse(vCollapsed, fmt=fmt, **kwargs)
    retval = pw

    return retval

def getLorenz(detectorGeom, *args):
    if len(args) == 2:
        tTh, eta = args
    elif len(args) == 0:
        tTh, eta = detectorGeom.xyoToAngAll()
    else:
        raise RuntimeError, 'pass either detectorGeom only for full frame or (detectorGeom, tth, eta)'
    theta  = num.asarray(0.5 * tTh)
    eta    = num.asarray(eta)
    sTh    = num.sin(theta)
    sTTh   = num.sin(tTh)
    cTh    = num.cos(theta)
    cTTh   = num.cos(tTh)
    cEt    = num.cos(eta)
    lorenz = ( 1.0 + cTTh*cTTh ) / (2.0 * sTTh) / num.sqrt( sTh*sTh + cEt*cEt * cTh*cTh )
    return lorenz

def getEtaResolution(detectorGeom, tTh):
    "angular resolution along eta"
    radius = num.linalg.norm(
        num.array(detectorGeom.angToXYO(tTh, 0.)) - num.array(detectorGeom.angToXYO(0., 0.))
        )
    'have 2*pi radians over 2*pi*radius pixels'
    etaRes = 1./radius
    return etaRes

def getTThResolution(detectorGeom, tTh):
    detectorGeom.angToXYO
    x0, y0 = detectorGeom.angToXYO(tTh, 0., units='pixels')
    tThP, eta = detectorGeom.xyoToAng(x0, y0+1., units='pixels')
    tThRes = tThP - tTh
    return tThRes

def drEtaOme(angCen, dEta, dOme):
    "compute true angular changes and dA corresponding to steps in eta and omega"
    tThNominal, eta, ome = angCen
    nVecsPatch = makeNVecs( # 3xN
        num.tile(tThNominal, 3),
        num.array([eta, eta+dEta, eta]),
        num.array([ome, ome, ome+dOme]),
        )
    dr_Eta = num.arccos(num.dot(nVecsPatch[:,0], nVecsPatch[:,1])) # * r, but r==1
    dr_Ome = num.arccos(num.dot(nVecsPatch[:,0], nVecsPatch[:,2])) # * r, but r==1
    # dA_cen = dr_Eta * dr_Ome # neglects cross-product factor
    dA_cen = num.linalg.norm(num.cross((nVecsPatch[:,1]-nVecsPatch[:,0]), (nVecsPatch[:,2]-nVecsPatch[:,0])))
    return dr_Eta, dr_Ome, dA_cen

def omeEtaGridToDA(tThNominal, etaEdges, omeEdges):
    """
    get grid patch areas, in the sense of solid angle (pole figure) coverage;
    """
    tTh = num.tile(tThNominal, len(etaEdges)*len(omeEdges))
    omeGrid, etaGrid = num.meshgrid(omeEdges, etaEdges)
    nVecs = makeNVecs(tTh, etaGrid, omeGrid) # 3xN
    nVecs = nVecs.reshape((3,  len(etaEdges), len(omeEdges)))
    nEta = len(etaEdges)-1
    nOme = len(omeEdges)-1
    dA = num.zeros((nEta, nOme))
    for iEta in range(nEta):
        for iOme in range(nOme):
            n0 = nVecs[:,iEta  ,iOme  ]
            n1 = nVecs[:,iEta+1,iOme  ]
            n2 = nVecs[:,iEta  ,iOme+1]
            # 'the following neglects the cross product factor'
            # dr_Eta = num.arccos(num.dot(n0, n1)) # * r, but r==1
            # dr_Ome = num.arccos(num.dot(n0, n2)) # * r, but r==1
            # dA[iEta, iOme] = dr_Eta * dr_Ome
            dA[iEta, iOme] = num.linalg.norm( num.cross( n1-n0, n2-n0 ) )
    return dA

def bboxIntersect3D(b1, b2):
    '''
    0-tolerance bounding box intersection in 3d
    '''
    (b1xmin, b1xmax), (b1ymin, b1ymax), (b1zmin,b1zmax) = b1
    (b2xmin, b2xmax), (b2ymin, b2ymax), (b2zmin,b2zmax) = b2

    if ((((b1xmin >= b2xmin) and (b1xmin <= b2xmax)) or
         ((b1xmax >= b2xmin) and (b1xmax <= b2xmax)) or
         ((b2xmin >= b1xmin) and (b2xmin <= b1xmax)) or
         ((b2xmax >= b1xmin) and (b2xmax <= b1xmax)))
        and
        (((b1ymin >= b2ymin) and (b1ymin <= b2ymax)) or
         ((b1ymax >= b2ymin) and (b1ymax <= b2ymax)) or
         ((b2ymin >= b1ymin) and (b2ymin <= b1ymax)) or
         ((b2ymax >= b1ymin) and (b2ymax <= b1ymax)))
        and
        (((b1zmin >= b2zmin) and (b1zmin <= b2zmax)) or
         ((b1zmax >= b2zmin) and (b1zmax <= b2zmax)) or
         ((b2zmin >= b1zmin) and (b2zmin <= b1zmax)) or
         ((b2zmax >= b1zmin) and (b2zmax <= b1zmax)))):
        retval = True
    else:
        retval = False
    return retval

def pfigFromSpots(spots, iHKL, phaseID=None,
                  nOme=None, nEta=None,
                  tThTol=None, etaExcl=0.1,
                  plot=False, plotPrefix=None, debug=False):
    """
    probably want to have collected spots without discarding those at boundaries (discardAtBounds=False)

    depending on the context, may need to have done something like iHKL = hklIDs.index(hklID)

    if nEta is negative, it is treated as the target lumping of pixels

    etaExcl is in radians -- etas within this range of +-pi/2 degrees are left out;
    can set to None to turn this behavior off

    can use tThTol to tighten down the two-theta tolerance
    """
    nEta = nEta or -3
    plot = plot or plotPrefix is not None

    planeData = spots.getPlaneData(phaseID=phaseID)
    hklTuple = tuple(planeData.getHKLs()[iHKL])
    tThNominal = planeData.getTTh()[iHKL]
    crystalVector = planeData.hklDataList[iHKL]['latPlnNrmls'][:,0]

    lorenz = getLorenz(spots.detectorGeom)
    tThAll, etaAll = spots.detectorGeom.xyoToAngAll()

    tThList, etaList, omeList, valList = [], [], [], []
    #
    for iSpot, angCoords in spots.getIterHKL(hklTuple, phaseID=phaseID):
        spot = spots.getSpotObj(iSpot)
        if tThTol is not None:
            relErr = abs(angCoords[0] - tThNominal) / tThTol
            keepThis = relErr < 1.
            if debug:
                print 'keep is %6s for spot %4d of size %9d intensity %9g due to tThTol, tTh relative error is %g' % (
                    str(keepThis), iSpot, len(spot), spot.vAll.sum(), relErr,
                    )
            if not keepThis:
                continue
        # print 'contributing spot '+str(iSpot)
        key, deltaOmega, (xAll, yAll, oAll, vAll) = spot.getDataMinimal()
        tThList.append(tThAll[xAll,yAll])
        etaList.append(etaAll[xAll,yAll])
        omeList.append(oAll)
        # 'apply Lorenz correction'
        vAllByL = vAll / lorenz[xAll,yAll]
        valList.append(vAllByL)

    if etaExcl is None:
        etaRanges = [( (-num.pi, num.pi) , (-num.pi, num.pi) )]
    else:
        e1 = 0.5*num.pi - etaExcl
        e2 = 0.5*num.pi + etaExcl
        e3 = 0.5*num.pi + num.pi - etaExcl
        etaRanges = [
            ( (-e1, e1), (-num.pi, num.pi) ),
            ( ( e2, e3), (0., 2.0*num.pi) ),
            ]

    tThList = num.hstack(tThList)
    etaList = num.hstack(etaList)
    omeList = num.hstack(omeList)
    valList = num.hstack(valList)

    valLow = min(1.0e-8, 1.0e-12*valList.max())
    valList[valList < valLow] = valLow

    retval = []
    for iRange, (etaRange, etaMM) in enumerate(etaRanges):

        if nEta < 0:
            radius = 1.0 / getEtaResolution(spots.detectorGeom, tThNominal)
            nEta   = int(((etaRange[1]-etaRange[0]) * radius) / abs(nEta))

        etaTemp = mapAngle(etaList, etaMM)
        indices = num.logical_and( etaTemp > etaRange[0], etaTemp < etaRange[1] )
        intensVals, etasCen, omesCen = \
            collapse(valList[indices],
                     etaTemp[indices],
                     omeList[indices],
                     nOme = ( nOme, spots.getOmegaMins().min(), spots.getOmegaMaxs().max() ),
                     nEta = ( nEta, etaRange[0], etaRange[1] ),
                     debug=debug
                     # averagedList=[tThList[indices]],
                     )

        if plot:
            fmt = FormatEtaOme(etasCen, omesCen, intensVals.todense())
            #pw = plotWrap.PlotWrap()
            #pw(intensVals.todense())
            pw = displaySparse(intensVals, fmt=fmt, colorUnder=(0.5,0.5,0.5), cmap=cm.jet)
            if plotPrefix is not None:
                pw.save(filename=str(plotPrefix)+'_omeEtaIntensVals_iHKL_%d_%d.png'%(iHKL,iRange))
                pw.destroy()

        dEta = num.mean(etasCen[1:]-etasCen[:-1])
        dOme = num.mean(omesCen[1:]-omesCen[:-1])
        etaEdges = num.hstack((etasCen[0] - 0.5*dEta, etasCen + 0.5*dEta))
        omeEdges = num.hstack((omesCen[0] - 0.5*dOme, omesCen + 0.5*dOme))
        dA = omeEtaGridToDA(tThNominal, etaEdges, omeEdges)
        'covert to dense and then (from matrix) to array, so that flatten works properly'
        pVals = num.array(intensVals.todense() / dA) # all missing points converted to zero
        nVecs = makeNVecs(
            num.tile( tThNominal, pVals.size ),
            etasCen, omesCen, asGrid=True)
        if debug:
            print 'making pfig dict with pVals nVecs shape %s and pVals size %s' % (str(nVecs.shape), str(pVals.size))
        retval.append( (
                omeEdges,
                etaEdges,
                intensVals,
                ) )

    return retval
