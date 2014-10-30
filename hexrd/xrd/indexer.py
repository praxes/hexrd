#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Joel Bernier <bernier2@llnl.gov> and others.
# LLNL-CODE-529294.
# All rights reserved.
#
# This file is part of HExrd. For details on dowloading the source,
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
import ctypes
import tempfile
import glob
import logging
import time
import pdb

import numpy as num
num.seterr(invalid='ignore')

import hexrd.matrixutil as mUtil

from hexrd.xrd.grain     import Grain, makeMeasuredScatteringVectors
from hexrd.xrd.rotations import \
     discreteFiber, mapAngle, \
     quatOfRotMat, quatProductMatrix, \
     rotMatOfExpMap, rotMatOfQuat
from hexrd.xrd.symmetry  import toFundamentalRegion
from hexrd.xrd           import xrdbase

from hexrd.xrd import transforms      as xf
from hexrd.xrd import transforms_CAPI as xfcapi

if xrdbase.haveMultiProc:
    multiprocessing = xrdbase.multiprocessing # formerly import


logger = logging.getLogger(__name__)


# module vars
piby2 = num.pi * 0.5
r2d = 180. / num.pi
d2r = num.pi / 180.

Xl = num.c_[1, 0, 0].T
Yl = num.c_[0, 1, 0].T
Zl = num.c_[0, 0, 1].T

fableSampCOB = num.dot( rotMatOfExpMap(piby2*Zl),
                        rotMatOfExpMap(piby2*Yl) )

# global vars for multiproc paintGrid method
planeData_MP  = None
omeMin_MP     = None
omeMax_MP     = None
etaMin_MP     = None
etaMax_MP     = None
hklList_MP    = None
hklIDs_MP     = None
etaOmeMaps_MP = None
bMat_MP       = None
threshold_MP  = None

class GrainSpotter:
    """
    Interface to grain spotter, which must be in the user's path
    """
    __execName = 'grainspotter'
    def __init__(self):
        self.__tempFNameList = []

        if (os.system('which '+self.__execName) != 0):
            print >> sys.stderr, "need %s to be in the path" % (self.__execName)
            raise RuntimeError, "unrecoverable error"

        return

    def __call__(self, spotsArray, **kwargs):
        """
        writes gve and ini files to system, calls grainspotter, parses results.

        A word on spacegroup numbers: it appears that grainspotter is using the 'VolA' tag for calls to SgInfo
        """
        location = self.__class__.__name__
        tic = time.time()

        phaseID   = None
        gVecFName = 'tmp'

        kwarglen = len(kwargs)
        if kwarglen > 0:
            argkeys = kwargs.keys()
            for i in range(kwarglen):
                if argkeys[i] == 'phaseID':
                    phaseID = kwargs[argkeys[i]]
                elif argkeys[i] == 'filename':
                    gVecFName = kwargs[argkeys[i]]
                    pass
                pass
            pass

        planeData = spotsArray.getPlaneData(phaseID=phaseID)
        U0        = planeData.latVecOps['U0']
        symTag    = planeData.getLaueGroup()

        writeGVE(spotsArray, gVecFName, **kwargs)

        toc = time.time()
        print 'in %s, setup took %g' % (location, toc-tic)
        tic = time.time()

        # tempFNameStdout = tempfile.mktemp()
        # self.__tempFNameList.append(tempFNameStdout)
        # tempFNameStdout = 'tmp.out'
        # grainSpotterCmd = '%s %s > %s' % (self.__execName, gVecFName, tempFNameStdout)
        grainSpotterCmd = '%s %s' % (self.__execName, gVecFName+'.ini')
        os.system(grainSpotterCmd)
        toc = time.time()
        print 'in %s, execution took %g' % (location, toc-tic)
        tic = time.time()

        # add output files to cleanup list
        # self.__tempFNameList += glob.glob(gVecFName+'.*')

        # collect data from gff file'
        gffFile = gVecFName+'.gff'
        gffData = num.loadtxt(gffFile)
        if gffData.ndim == 1:
            gffData = gffData.reshape(1, len(gffData))
        gffData_U = gffData[:,6:6+9]

        # process for output
        retval = convertUToRotMat(gffData_U, U0, symTag=symTag)

        toc = time.time()
        print 'in %s, post-processing took %g' % (location, toc-tic)
        tic = time.time()

        return retval

    def __del__(self):
        self.cleanup()
        return
    def cleanup(self):
        for fname in self.__tempFNameList:
            os.remove(fname)
        return

def convertUToRotMat(Urows, U0, symTag='Oh', display=False):
    """
    Takes GrainSpotter gff ouput in rows

    U11 U12 U13 U21 U22 U23 U13 U23 U33

    and takes it into the hexrd/APS frame of reference

    Urows comes from grainspotter's gff output
    U0 comes from xrd.crystallography.latticeVectors.U0
    """

    numU, testDim = Urows.shape
    assert testDim == 9, "Your input must have 9 columns; yours has %d" % (testDim)

    qin  = quatOfRotMat(Urows.reshape(numU, 3, 3))
    qout = num.dot( quatProductMatrix( quatOfRotMat(fableSampCOB), mult='left' ), \
                    num.dot( quatProductMatrix( quatOfRotMat(U0.T), mult='right'),  \
                             qin ).squeeze() ).squeeze()
    if qout.ndim == 1:
        qout = toFundamentalRegion(qout.reshape(4, 1), crysSym=symTag, sampSym=None)
    else:
        qout = toFundamentalRegion(qout, crysSym=symTag, sampSym=None)
    if display:
        print "quaternions in (Fable convention):"
        print qin.T
        print "quaternions out (hexrd convention, symmetrically reduced)"
        print qout.T
        pass
    Uout = rotMatOfQuat(qout)
    return Uout

def convertRotMatToFableU(rMats, U0=num.eye(3), symTag='Oh', display=False):
    """
    Makes GrainSpotter gff ouput

    U11 U12 U13 U21 U22 U23 U13 U23 U33

    and takes it into the hexrd/APS frame of reference

    Urows comes from grainspotter's gff output
    U0 comes from xrd.crystallography.latticeVectors.U0
    """
    numU = num.shape(num.atleast_3d(rMats))[0]

    qin  = quatOfRotMat(num.atleast_3d(rMats))
    qout = num.dot( quatProductMatrix( quatOfRotMat(fableSampCOB.T), mult='left' ), \
                    num.dot( quatProductMatrix( quatOfRotMat(U0), mult='right'),  \
                             qin ).squeeze() ).squeeze()
    if qout.ndim == 1:
        qout = toFundamentalRegion(qout.reshape(4, 1), crysSym=symTag, sampSym=None)
    else:
        qout = toFundamentalRegion(qout, crysSym=symTag, sampSym=None)
    if display:
        print "quaternions in (hexrd convention):"
        print qin.T
        print "quaternions out (Fable convention, symmetrically reduced)"
        print qout.T
        pass
    Uout = rotMatOfQuat(qout)
    return Uout

######################################################################

"""
things for doing fiberSearch with multiprocessing;
multiprocessing has a hard time pickling a function defined in the local scope of another function,
so stuck putting the function out here;
"""
debugMultiproc = 0
if xrdbase.haveMultiProc:
    foundFlagShared = multiprocessing.Value(ctypes.c_bool)
    foundFlagShared.value = False
multiProcMode_MP   = None
spotsArray_MP      = None
candidate_MP       = None
dspTol_MP          = None
minCompleteness_MP = None
doRefinement_MP    = None
nStdDev_MP         = None
def testThisQ(thisQ):
    """
    NOTES:
    (*) doFit is not done here -- in multiprocessing, that would end
        up happening on a remote process and then different processes
        would have different data, unless spotsArray were made to be
        fancier

    (*) kludge stuff so that this function is outside of fiberSearch
    """
    global multiProcMode_MP
    global spotsArray_MP
    global candidate_MP
    global dspTol_MP
    global minCompleteness_MP
    global doRefinement_MP
    global nStdDev_MP
    # assign locals
    multiProcMode   = multiProcMode_MP
    spotsArray      = spotsArray_MP
    candidate       = candidate_MP
    dspTol          = dspTol_MP
    minCompleteness = minCompleteness_MP
    doRefinement    = doRefinement_MP
    nStdDev         = nStdDev_MP
    nSigmas = 2                         # ... make this a settable option?
    if multiProcMode:
        global foundFlagShared

    foundGrainData = None
    #print "testing %d of %d"% (iR+1, numTrials)
    thisRMat = rotMatOfQuat(thisQ)

    ppfx = ''
    if multiProcMode:
        proc = multiprocessing.current_process()
        ppfx = str(proc.name)+' : '
        if multiProcMode and foundFlagShared.value:
            """
            map causes this function to be applied to all trial orientations,
            but skip evaluations after an acceptable grain has been found
            """
            if debugMultiproc > 1:
                print ppfx+'skipping on '+str(thisQ)
            return foundGrainData
        else:
            if debugMultiproc > 1:
                print ppfx+'working on  '+str(thisQ)
            pass
    candidate.findMatches(rMat=thisRMat,
                          strainMag=dspTol,
                          claimingSpots=False,
                          testClaims=True,
                          updateSelf=True)
    if debugMultiproc > 1:
        print ppfx+' for '+str(thisQ)+' got completeness : '+str(candidate.completeness)
    if candidate.completeness >= minCompleteness:
        ## attempt to filter out 'junk' spots here by performing full
        ## refinement before claiming
        fineEtaTol = candidate.etaTol
        fineOmeTol = candidate.omeTol
        if doRefinement:
            if multiProcMode and foundFlagShared.value:
                'some other process beat this one to it'
                return foundGrainData
            print ppfx+"testing candidate q = [%1.2e, %1.2e, %1.2e, %1.2e]" %tuple(thisQ)
            # not needed # candidate.fitPrecession(display=False)
            ## first fit
            candidate.fit(display=False)
            ## auto-tolerace based on statistics of current matches
            validRefls = candidate.grainSpots['iRefl'] > 0
            fineEtaTol = nStdDev * num.std(candidate.grainSpots['diffAngles'][validRefls, 1])
            fineOmeTol = nStdDev * num.std(candidate.grainSpots['diffAngles'][validRefls, 2])
            ## next fits with finer tolerances
            for iLoop in range(3):
                candidate.findMatches(etaTol=fineEtaTol,
                                      omeTol=fineOmeTol,
                                      claimingSpots=False,
                                      testClaims=True,
                                      updateSelf=True)
                # not needed # candidate.fitPrecession(display=False)
                candidate.fit(display=False)
            if candidate.completeness < minCompleteness:
                print ppfx+"candidate failed"
                return foundGrainData
            if multiProcMode and foundFlagShared.value:
                'some other process beat this one to it'
                return foundGrainData
            # not needed # candidate.fitPrecession(display=False)
            # not needed? # candidate.fit(display=False)
            # not needed? # candidate.findMatches(etaTol=fineEtaTol,
            # not needed? #                       omeTol=fineOmeTol,
            # not needed? #                       claimingSpots=False,
            # not needed? #                       testClaims=True,
            # not needed? #                       updateSelf=True)
        else:
            ## at least do precession correction
            candidate.fitPrecession(display=False)
            candidate.findMatches(rMat=thisRMat,
                                  strainMag=dspTol,
                                  claimingSpots=False,
                                  testClaims=True,
                                  updateSelf=True)
            fineEtaTol = candidate.etaTol
            fineOmeTol = candidate.omeTol
            if candidate.completeness < minCompleteness:
                print ppfx+"candidate failed"
                return foundGrainData
            if multiProcMode and foundFlagShared.value:
                'some other process beat this one to it'
                return foundGrainData
        if multiProcMode:
            foundFlagShared.value = True
        # # newGrain uses current candidate.rMat
        # # do not do claims here -- those handled outside of this call
        # foundGrain = candidate.newGrain(
        #     spotsArray, claimingSpots=False,
        #     omeTol=fineOmeTol,
        #     etaTol=fineEtaTol)
        # if multiProcMode:
        #     foundGrain.strip()
        cInfo = quatOfRotMat(candidate.rMat).flatten().tolist()
        cInfo.append(candidate.completeness)
        print ppfx+"Grain found at q = [%1.2e, %1.2e, %1.2e, %1.2e] with completeness %g" \
              % tuple(cInfo)
        foundGrainData = candidate.getGrainData()
        'tolerances not actually set in candidate, so set them manually'
        foundGrainData['omeTol'] = fineOmeTol
        foundGrainData['etaTol'] = fineEtaTol

    return foundGrainData
#
def fiberSearch(spotsArray, hklList,
                iPhase=0,
                nsteps=120,
                minCompleteness=0.60,
                minPctClaimed=0.95,
                preserveClaims=False,
                friedelOnly=True,
                dspTol=None,
                etaTol=0.025,
                omeTol=0.025,
                etaTolF=0.00225,
                omeTolF=0.00875,
                nStdDev=2,
                quitAfter=None,
                doRefinement=True,
                debug=True,
                doMultiProc=True,
                nCPUs=None,
                outputGrainList=False
                ):
    """
    This indexer finds grains by performing 1-d searches along the fibers under the
    valid spots associated with each reflection order specified in hklList.  The set
    of spots used to generate the candidate orientations may be restricted to Friedel
    pairs only.

    hklList *must* have length > 0;
    Dach hkl entry in hklList *must* be a tuple, not a list

    the output is a concatenated list of orientation matrices ((n, 3, 3) numpy.ndarray).
    """

    assert hasattr(hklList, '__len__'), "the HKL list must have length, and len(hklList) > 0."

    nHKLs = len(hklList)
    grainList = []
    nGrains = 0
    planeData = spotsArray.getPlaneData(iPhase)
    csym = planeData.getLaueGroup()
    bMat = planeData.latVecOps['B']
    if dspTol is None:
        dspTol = planeData.strainMag

    centroSymRefl = planeData.getCentroSymHKLs()

    candidate = Grain(spotsArray, rMat=None,
                      etaTol=etaTol, omeTol=omeTol)
    multiProcMode = xrdbase.haveMultiProc and doMultiProc
    #
    global foundFlagShared
    global multiProcMode_MP
    global spotsArray_MP
    global candidate_MP
    global dspTol_MP
    global minCompleteness_MP
    global doRefinement_MP
    global nStdDev_MP
    multiProcMode_MP   = multiProcMode
    spotsArray_MP      = spotsArray
    candidate_MP       = candidate
    dspTol_MP          = dspTol
    minCompleteness_MP = minCompleteness
    doRefinement_MP    = doRefinement
    nStdDev_MP         = nStdDev
    """
    set up for shared memory multiprocessing
    """
    if multiProcMode:
        nCPUs = nCPUs or xrdbase.dfltNCPU
        spotsArray.multiprocMode = True
        pool = multiprocessing.Pool(nCPUs)

    """
    HKL ITERATOR
    """
    if isinstance(quitAfter, dict):
        n_hkls_to_search = quitAfter['nHKLs']
    else:
        n_hkls_to_search = nHKLs

    if isinstance(quitAfter, int):
        quit_after_ngrains = quitAfter
    else:
        quit_after_ngrains = 0

    numTotal = len(spotsArray)
    pctClaimed = 0.
    time_to_quit = False
    tic = time.time()

    for iHKL in range(n_hkls_to_search):
        print "\n#####################\nProcessing hkl %d of %d\n" % (iHKL+1, nHKLs)
        thisHKLID = planeData.getHKLID(hklList[iHKL])
        thisRingSpots0   = spotsArray.getHKLSpots(thisHKLID)
        thisRingSpots0W  = num.where(thisRingSpots0)[0]
        unclaimedOfThese = -spotsArray.checkClaims(indices=thisRingSpots0W)
        thisRingSpots    = copy.deepcopy(thisRingSpots0)
        thisRingSpots[thisRingSpots0W] = unclaimedOfThese
        if friedelOnly:
            # first, find Friedel Pairs
            spotsArray.findFriedelPairsHKL(hklList[iHKL],
                                           etaTol=etaTolF,
                                           omeTol=omeTolF)
            spotsIteratorI = spotsArray.getIterHKL(hklList[iHKL], unclaimedOnly=True, friedelOnly=True)
            # make some stuff for counters
            maxSpots = 0.5*(sum(thisRingSpots) - sum(spotsArray.friedelPair[thisRingSpots] == -1))
        else:
            spotsIteratorI = spotsArray.getIterHKL(hklList[iHKL], unclaimedOnly=True, friedelOnly=False)
            maxSpots = sum(thisRingSpots)
        """
        SPOT ITERATOR
          - this is where we iterate over all 'valid' spots for the current HKL as
            subject to the conditions of claims and ID as a friedel pair (when requested)
        """
        for iRefl, stuff in enumerate(spotsIteratorI):
            unclaimedOfThese = -spotsArray.checkClaims(indices=thisRingSpots0W)
            thisRingSpots    = copy.deepcopy(thisRingSpots0)
            thisRingSpots[thisRingSpots0W] = unclaimedOfThese
            if friedelOnly:
                iSpot, jSpot, angs_I, angs_J = stuff

                Gplus  = makeMeasuredScatteringVectors(*angs_I)
                Gminus = makeMeasuredScatteringVectors(*angs_J)

                Gvec = 0.5*(Gplus - Gminus)
                maxSpots = 0.5*(sum(thisRingSpots) - sum(spotsArray.friedelPair[thisRingSpots] == -1))
            else:
                iSpot, angs_I = stuff
                Gvec  = makeMeasuredScatteringVectors(*angs_I)
                maxSpots = sum(thisRingSpots)
            print "\nProcessing reflection %d (spot %d), %d remain unclaimed\n" % (iRefl+1, iSpot, maxSpots)
            if multiProcMode and debugMultiproc > 1:
                marks = spotsArray._Spots__marks[:]
                print 'marks : '+str(marks)
            # make the fiber;
            qfib = discreteFiber(hklList[iHKL], Gvec,
                                 B=bMat,
                                 ndiv=nsteps,
                                 invert=False,
                                 csym=csym, ssym=None)[0]
            # if +/- hkl aren't in the symmetry group, need '-' fiber
            if not centroSymRefl[thisHKLID]:
                minusHKL = -num.r_[hklList[iHKL]]
                qfibM = discreteFiber(minusHKL, Gvec,
                                      B=bMat,
                                      ndiv=nsteps,
                                      invert=False,
                                      csym=csym, ssym=None)[0]
                qfib = num.hstack([qfib, qfibM])
                pass
            # cull out duplicate orientations
            qfib = mUtil.uniqueVectors(qfib, tol=1e-4)
            numTrials = qfib.shape[1]
            """
            THIS IS THE BIGGIE; THE LOOP OVER THE DISCRETE ORIENTATIONS IN THE CURRENT FIBER
            """
            if multiProcMode:
                foundFlagShared.value = False
                qfibList = map(num.array, qfib.T.tolist())
                #if debugMultiproc:
                #    print 'qfibList : '+str(qfibList)
                results = num.array(pool.map(testThisQ, qfibList, chunksize=1))
                trialGrains = results[num.where(num.array(results, dtype=bool))]
                # for trialGrain in trialGrains:
                #     trialGrain.restore(candidate)
            else:
                trialGrains = []
                for iR in range(numTrials):
                    foundGrainData = testThisQ(qfib[:, iR])
                    if foundGrainData is not None:
                        trialGrains.append(foundGrainData)
                        break
            'end of if multiProcMode'

            if len(trialGrains) == 0:
                print "No grain found containing spot %d\n" % (iSpot)
                # import pdb;pdb.set_trace()
            else:
                asMaster = multiProcMode
                'sort based on completeness'
                trialGrainCompletenesses = [tgd['completeness'] for tgd in trialGrains]
                order = num.argsort(trialGrainCompletenesses)[-1::-1]
                for iTrialGrain in order:
                    foundGrainData = trialGrains[iTrialGrain]
                    foundGrain = Grain(spotsArray, grainData=foundGrainData, claimingSpots=False)
                    'check completeness before accepting, especially important for multiproc'
                    foundGrain.checkClaims() # updates completeness
                    if debugMultiproc:
                        print 'final completeness of candidate is %g' % (foundGrain.completeness)
                    if foundGrain.completeness >= minCompleteness:
                        conflicts = foundGrain.claimSpots(asMaster=asMaster)
                        numConfl = num.sum(conflicts)
                        if numConfl > 0:
                            'tried to claim %d spots that are already claimed' % (numConfl)
                        grainList.append(foundGrain)
                        nGrains += 1
                numUnClaimed = num.sum(-spotsArray.checkClaims())
                numClaimed = numTotal - numUnClaimed
                pctClaimed = num.float(numClaimed) / numTotal
                print "Found %d grains so far, %f%% claimed" % (nGrains,100*pctClaimed)

                time_to_quit = (pctClaimed > minPctClaimed) or\
                  ((quit_after_ngrains > 0) and (nGrains >= quit_after_ngrains))
                if time_to_quit:
                    break
        'end of iRefl loop'

        if time_to_quit:
            break

    'end of iHKL loop'
    rMats = num.empty((len(grainList), 3, 3))
    for i in range(len(grainList)):
        rMats[i, :, :] = grainList[i].rMat

    if outputGrainList:
        retval = (rMats, grainList)
    else:
        retval = rMats

    if not preserveClaims:
        spotsArray.resetClaims()
    toc = time.time()
    print 'fiberSearch execution took %g seconds' % (toc-tic)

    if multiProcMode:
        pool.close()
        spotsArray.multiprocMode = False
        foundFlagShared.value = False
    # global foundFlagShared
    # global multiProcMode_MP
    # global spotsArray_MP
    # global candidate_MP
    # global dspTol_MP
    # global minCompleteness_MP
    # global doRefinement_MP
    multiProcMode_MP = None
    spotsArray_MP = None
    candidate_MP = None
    dspTol_MP = None
    minCompleteness_MP = None
    doRefinement_MP = None

    return retval

def pgRefine(x, etaOmeMaps, omegaRange, threshold):
    phi = sum(x*x)
    if phi < 1e-7:
        q = [num.r_[1.,0.,0.,0.],]
    else:
        phi = num.sqrt(phi)
        n = (1. / phi) * x.flatten()
        cphi2 = num.cos(0.5*phi)
        sphi2 = num.sin(0.5*phi)
        q = [num.r_[cphi2, sphi2*n[0], sphi2*n[1], sphi2*n[2]],]
    c = paintGrid(q, etaOmeMaps, threshold=threshold, bMat=None, omegaRange=omegaRange, etaRange=None, debug=False, progressBar=False)
    f = abs(1. - c)
    return f

def paintGrid(quats, etaOmeMaps,
              threshold=None, bMat=None,
              omegaRange=None, etaRange=None,
              omeTol=d2r, etaTol=d2r,
              omePeriod=(-num.pi, num.pi),
              progressBar=False, doMultiProc=False,
              nCPUs=None, debug=False):
    """
    do a direct search of omega-eta maps to paint each orientation in
    quats with a completeness

    bMat is in CRYSTAL frame

    etaOmeMaps is instance of xrd.xrdutil.CollapseOmeEta

    omegaRange=([-num.pi/3., num.pi/3.],) for example

    *) lifted mainly from grain.py

    *) self.etaGrid, self.omeGrid = num.meshgrid(self.etaEdges, self.omeEdges)
       this means that ETA VARIES FASTEST!

    ...make a new function that gets called by grain to do the g-vec angle computation?
    """

    quats = num.atleast_2d(quats)
    if quats.size == 4:
        quats = quats.reshape(4, 1)

    planeData = etaOmeMaps.planeData

    hklIDs    = num.r_[etaOmeMaps.iHKLList]
    hklList   = num.atleast_2d(planeData.hkls[:, hklIDs].T).tolist()
    nHKLS     = len(hklIDs)

    numEtas   = len(etaOmeMaps.etaEdges) - 1
    numOmes   = len(etaOmeMaps.omeEdges) - 1

    if threshold is None:
        threshold = num.zeros(nHKLS)
        for i in range(nHKLS):
            threshold[i] = num.mean( num.r_[num.mean(etaOmeMaps.dataStore[i]),
                                            num.median(etaOmeMaps.dataStore[i]) ] )
    elif threshold is not None and not hasattr(threshold, '__len__'):
        threshold = threshold * num.ones(nHKLS)
    elif hasattr(threshold, '__len__'):
        if len(threshold) != nHKLS:
            raise RuntimeError, "threshold list is wrong length!"
        else:
            print "INFO: using list of threshold values"
    else:
        raise RuntimeError, "unknown threshold option.  should be a list of numbers or None"
    if bMat is None:
        bMat = planeData.latVecOps['B']

    """
    index munging here -- look away

    order of ome-eta map arrays is (i, j) --> (ome, eta)
    i.e. eta varies fastest.
    """
    # mapIndices = num.indices([numEtas, numOmes])
    # etaIndices = mapIndices[0].flatten()
    # omeIndices = mapIndices[1].T.flatten()
    # etaIndices = num.tile(range(numEtas), (numOmes))
    # omeIndices = num.tile(range(numOmes), (numEtas))
    # j_eta, i_ome = np.meshgrid(range(numEtas), range(numOmes))
    # etaIndices = j_eta.flatten()
    # omeIndices = i_ome.flatten()
    etaIndices = num.r_[range(numEtas)]
    omeIndices = num.r_[range(numOmes)]

    omeMin = None
    omeMax = None
    if omegaRange is None:              # this NEEDS TO BE FIXED!
        omeMin = [num.min(etaOmeMaps.omeEdges),]
        omeMax = [num.max(etaOmeMaps.omeEdges),]
    else:
        omeMin = [omegaRange[i][0] for i in range(len(omegaRange))]
        omeMax = [omegaRange[i][1] for i in range(len(omegaRange))]
    etaMin = None
    etaMax = None
    if etaRange is not None:
        etaMin = [etaRange[i][0] for i in range(len(etaRange))]
        etaMax = [etaRange[i][1] for i in range(len(etaRange))]

    # obselete # # make list of rMats from input quats
    # obselete # rMatsList = [rotMatOfQuat(quats[:, i]) for i in range(quats.shape[1])]

    multiProcMode = xrdbase.haveMultiProc and doMultiProc

    if multiProcMode:
        nCPUs = nCPUs or xrdbase.dfltNCPU
        chunksize = min(quats.shape[1] // nCPUs, 10)
        logger.info(
            "using multiprocessing with %d processes and a chunk size of %d",
            nCPUs, chunksize
            )
    else:
        logger.info("running in serial mode")
        nCPUs = 1

    # assign the globals for paintGridThis
    global symHKLs_MP, wavelength_MP
    global omeMin_MP, omeMax_MP, omeTol_MP, omePeriod_MP
    global etaMin_MP, etaMax_MP, etaTol_MP
    global omeIndices_MP, etaIndices_MP
    global omeEdges_MP, etaEdges_MP
    global hklList_MP, hklIDs_MP
    global etaOmeMaps_MP
    global bMat_MP
    global threshold_MP
    symHKLs_MP    = planeData.getSymHKLs()
    wavelength_MP = planeData.wavelength
    hklIDs_MP     = hklIDs
    hklList_MP    = hklList
    omeMin_MP     = omeMin
    omeMax_MP     = omeMax
    omeTol_MP     = omeTol
    omeIndices_MP = omeIndices
    omePeriod_MP  = omePeriod
    omeEdges_MP   = etaOmeMaps.omeEdges
    etaMin_MP     = etaMin
    etaMax_MP     = etaMax
    etaTol_MP     = etaTol
    etaIndices_MP = etaIndices
    etaEdges_MP   = etaOmeMaps.etaEdges
    etaOmeMaps_MP = etaOmeMaps.dataStore
    bMat_MP       = bMat
    threshold_MP  = threshold

    # do the mapping
    start = time.time()                      # time this
    retval = None
    if multiProcMode:
        pool = multiprocessing.Pool(nCPUs)
        retval = pool.map(paintGridThis, quats.T, chunksize=chunksize)
    else:
        retval = map(paintGridThis, quats.T)
    elapsed = (time.time() - start)
    logger.info("paintGrid took %.3f seconds", elapsed)

    symHKLs_MP    = None
    wavelength_MP = None
    hklIDs_MP     = None
    hklList_MP    = None
    omeMin_MP     = None
    omeMax_MP     = None
    omeIndices_MP = None
    omePeriod_MP  = None
    omeEdges_MP   = None
    etaMin_MP     = None
    etaMax_MP     = None
    etaIndices_MP = None
    etaEdges_MP   = None
    etaOmeMaps_MP = None
    bMat_MP       = None
    threshold_MP  = None

    if multiProcMode:
        pool.close()

    return retval

def paintGridThis(quat):
    """
    """
    # pull local vars from globals set in paintGrid
    global symHKLs_MP, wavelength_MP
    global omeMin_MP, omeMax_MP, omeTol_MP, omePeriod_MP
    global etaMin_MP, etaMax_MP, etaTol_MP
    global omeIndices_MP, etaIndices_MP
    global omeEdges_MP, etaEdges_MP
    global hklList_MP, hklIDs_MP
    global etaOmeMaps_MP
    global bMat_MP
    global threshold_MP
    symHKLs    = symHKLs_MP
    wavelength = wavelength_MP
    hklIDs     = hklIDs_MP
    hklList    = hklList_MP
    omeMin     = omeMin_MP
    omeMax     = omeMax_MP
    omeTol     = omeTol_MP
    omePeriod  = omePeriod_MP
    omeIndices = omeIndices_MP
    omeEdges   = omeEdges_MP
    etaMin     = etaMin_MP
    etaMax     = etaMax_MP
    etaTol     = etaTol_MP
    etaIndices = etaIndices_MP
    etaEdges   = etaEdges_MP
    etaOmeMaps = etaOmeMaps_MP
    bMat       = bMat_MP
    threshold  = threshold_MP

    # need this for proper index generation

    omegas = [omeEdges[0] + (i+0.5)*(omeEdges[1] - omeEdges[0]) for i in range(len(omeEdges) - 1)]
    etas   = [etaEdges[0] + (i+0.5)*(etaEdges[1] - etaEdges[0]) for i in range(len(etaEdges) - 1)]

    delOmeSign = num.sign(omegas[1] - omegas[0])

    del_ome = abs(omegas[1] - omegas[0])
    del_eta = abs(etas[1] - etas[0])

    dpix_ome = round(omeTol / del_ome)
    dpix_eta = round(etaTol / del_eta)

    debug = False
    if debug:
        print "using ome, eta dilitations of (%d, %d) pixels" % (dpix_ome, dpix_eta)

    nHKLs = len(hklIDs)

    rMat = rotMatOfQuat(quat)

    nPredRefl = 0
    nMeasRefl = 0
    reflInfoList = []
    dummySpotInfo = num.nan * num.ones(3)

    hklCounterP = 0                 # running count of excpected (predicted) HKLs
    hklCounterM = 0                 # running count of "hit" HKLs
    for iHKL in range(nHKLs):
        # select and C-ify symmetric HKLs
        these_hkls = num.array(symHKLs[hklIDs[iHKL]].T, dtype=float, order='C')

        # oscillation angle arrays
        oangs   = xfcapi.oscillAnglesOfHKLs(these_hkls, 0., rMat, bMat, wavelength)
        angList = num.vstack(oangs)
        if not num.all(num.isnan(angList)):
            idx = -num.isnan(angList[:, 0])
            angList = angList[idx, :]
            angList[:, 1] = xf.mapAngle(angList[:, 1])
            angList[:, 2] = xf.mapAngle(angList[:, 2], omePeriod)

            if omeMin is None:
                omeMin = [-num.pi, ]
                omeMax = [ num.pi, ]
            if etaMin is None:
                etaMin = [-num.pi, ]
                etaMax = [ num.pi, ]

            angMask = num.logical_and(
                xf.validateAngleRanges(angList[:, 1], etaMin, etaMax),
                xf.validateAngleRanges(angList[:, 2], omeMin, omeMax))

            allAngs_m = angList[angMask, :]

            # not output # # duplicate HKLs
            # not output # allHKLs_m = num.vstack([these_hkls, these_hkls])[angMask, :]

            culledTTh  = allAngs_m[:, 0]
            culledEta  = allAngs_m[:, 1]
            culledOme  = allAngs_m[:, 2]
            # not output # culledHKLs = allHKLs_m.T

            nThisPredRefl = len(culledTTh)
            hklCounterP += nThisPredRefl
            for iTTh in range(nThisPredRefl):
                culledEtaIdx = num.where(etaEdges - culledEta[iTTh] > 0)[0]
                if len(culledEtaIdx) > 0:
                    culledEtaIdx = culledEtaIdx[0] - 1
                    if culledEtaIdx < 0:
                        culledEtaIdx = None
                else:
                    culledEtaIdx = None
                culledOmeIdx = num.where(omeEdges - culledOme[iTTh] > 0)[0]
                if len(culledOmeIdx) > 0:
                    if delOmeSign > 0:
                        culledOmeIdx = culledOmeIdx[0] - 1
                    else:
                        culledOmeIdx = culledOmeIdx[-1]
                    if culledOmeIdx < 0:
                        culledOmeIdx = None
                else:
                    culledOmeIdx = None

                if culledEtaIdx is not None and culledOmeIdx is not None:
                    if dpix_ome > 0 or dpix_eta > 0:
                        i_dil, j_dil = num.meshgrid(num.arange(-dpix_ome, dpix_ome + 1),
                                                    num.arange(-dpix_eta, dpix_eta + 1))
                        i_sup = omeIndices[culledOmeIdx] + num.array([i_dil.flatten()], dtype=int)
                        j_sup = etaIndices[culledEtaIdx] + num.array([j_dil.flatten()], dtype=int)
                        # catch shit that falls off detector...
                        # ...maybe make this fancy enough to wrap at 2pi?
                        i_max, j_max = etaOmeMaps[iHKL].shape
                        idx_mask = num.logical_and(num.logical_and(i_sup >= 0, i_sup < i_max),
                                                   num.logical_and(j_sup >= 0, j_sup < j_max))
                        pixelVal = etaOmeMaps[iHKL][i_sup[idx_mask], j_sup[idx_mask]]
                    else:
                        pixelVal = etaOmeMaps[iHKL][omeIndices[culledOmeIdx], etaIndices[culledEtaIdx] ]
                    isHit = num.any(pixelVal >= threshold[iHKL])
                    if isHit:
                        hklCounterM += 1
                        pass
                    pass
                # disabled # if debug:
                # disabled #     print "hkl %d -->\t%d\t%d\t%d\t" % (iHKL, culledHKLs[0, iTTh], culledHKLs[1, iTTh], culledHKLs[2, iTTh]) + \
                # disabled #           "isHit=%d\tpixel value: %g\n" % (isHit, pixelVal) + \
                # disabled #           "eta index: %d,%d\tetaP: %g\n" % (culledEtaIdx, etaIndices[culledEtaIdx], r2d*culledEta[iTTh]) + \
                # disabled #           "ome index: %d,%d\tomeP: %g\n" % (culledOmeIdx, omeIndices[culledOmeIdx], r2d*culledOme[iTTh])
                # disabled #     pass
                pass # close conditional on valid reflections
            pass # close loop on signed reflections
        pass # close loop on HKL
    if hklCounterP == 0:
        retval = 0.
    else:
        retval = float(hklCounterM) / float(hklCounterP)
    return retval

def progress_bar(progress):
    """
    Prints a progress bar to stdout.

    Inputs:
        progress - a float between 0. and 1.

    Example:
        >> progress_bar(0.7)
            |===================================               |
    """
    progressConditionStr = "ERROR: The argument of function NeuroTools.plotting.progress_bar(...) must be a float between 0. and 1.!"
    assert (type(progress) == float) and (progress >= 0.) and (progress <= 1.), progressConditionStr
    length = 50
    filled = int(round(length*progress))
    print "|" + "=" * filled + " " * (length-filled) + "|\r",
    sys.stdout.flush()

def writeGVE(spotsArray, fileroot, **kwargs):
    """
    write Fable gve file from Spots object

    fileroot is the root string used to write the gve and ini files

    Outputs:

    No return value, but writes the following files:

    <fileroot>.gve
    <fileroot>_grainSpotter.ini (points to --> <fileroot>_grainSpotter.log)

    Keyword arguments:

    Mainly for GrainSpotter .ini file, but some are needed for gve files

    keyword        default              definitions
    -----------------------------------------------------------------------------------------
    'sgNum':       <225>
    'phaseID':     <None>
    'cellString':  <F>
    'omeRange':    <-60, 60, 120, 240>  the oscillation range(s) **currently pulls from spots
    'deltaOme':    <0.25, 0.25>         the oscillation delta(s) **currently pulls from spots
    'minMeas':     <24>
    'minCompl':    <0.7>
    'minUniqn':    <0.5>
    'uncertainty': <[0.10, 0.25, .50]>  the min [tTh, eta, ome] uncertainties in degrees
    'eulStep':     <2>
    'nSigmas':     <2>
    'minFracG':    <0.90>
    'nTrials':     <100000>
    'positionfit': <True>

    Notes:

    *) The omeRange is currently pulled from the spotsArray input; the kwarg has no effect
       as of now.  Will change this to 'override' the spots info if the user, say, wants to
       pare down the range.

    *) There is no etaRange argument yet, but presumably GrainSpotter knows how to deal
       with this.  Pending feature...
    """
    # check on fileroot
    assert isinstance(fileroot, str)

    # keyword argument processing
    phaseID     = None
    sgNum       = 225
    cellString  = 'P'
    omeRange    = num.r_[-60, 60]   # in DEGREES
    deltaOme    = 0.25              # in DEGREES
    minMeas     = 24
    minCompl    = 0.7
    minUniqn    = 0.5
    uncertainty = [0.10, 0.25, .50] # in DEGREES
    eulStep     = 2                 # in DEGREES
    nSigmas     = 2
    minFracG    = 0.90
    numTrials   = 100000
    positionFit = True

    kwarglen = len(kwargs)
    if kwarglen > 0:
        argkeys = kwargs.keys()
        for i in range(kwarglen):
            if argkeys[i] == 'sgNum':
                sgNum = kwargs[argkeys[i]]
            elif argkeys[i] == 'phaseID':
                phaseID = kwargs[argkeys[i]]
            elif argkeys[i] == 'cellString':
                cellString = kwargs[argkeys[i]]
            elif argkeys[i] == 'omeRange':
                omeRange = kwargs[argkeys[i]]
            elif argkeys[i] == 'deltaOme':
                deltaOme = kwargs[argkeys[i]]
            elif argkeys[i] == 'minMeas':
                minMeas = kwargs[argkeys[i]]
            elif argkeys[i] == 'minCompl':
                minCompl = kwargs[argkeys[i]]
            elif argkeys[i] == 'minUniqn':
                minUniqn = kwargs[argkeys[i]]
            elif argkeys[i] == 'uncertainty':
                uncertainty = kwargs[argkeys[i]]
            elif argkeys[i] == 'eulStep':
                eulStep = kwargs[argkeys[i]]
            elif argkeys[i] == 'nSigmas':
                nSigmas = kwargs[argkeys[i]]
            elif argkeys[i] == 'minFracG':
                minFracG = kwargs[argkeys[i]]
            elif argkeys[i] == 'nTrials':
                numTrials = kwargs[argkeys[i]]
            elif argkeys[i] == 'positionfit':
                positionFit = kwargs[argkeys[i]]
            else:
                raise RuntimeError, "Unrecognized keyword argument '%s'" % (argkeys[i])

    # grab some detector geometry parameters for gve file header
    mmPerPixel = float(spotsArray.detectorGeom.pixelPitch) # ...these are still hard-coded to be square
    nrows_p = spotsArray.detectorGeom.nrows - 1
    ncols_p = spotsArray.detectorGeom.ncols - 1

    row_p, col_p = spotsArray.detectorGeom.pixelIndicesOfCartesianCoords(spotsArray.detectorGeom.xc,
                                                                         spotsArray.detectorGeom.yc)
    yc_p = ncols_p - col_p
    zc_p = nrows_p - row_p

    wd_mu = spotsArray.detectorGeom.workDist * 1e3 # in microns (Soeren)

    osc_axis = num.dot(fableSampCOB.T, Yl).flatten()

    # start grabbing stuff from planeData
    planeData = spotsArray.getPlaneData(phaseID=phaseID)
    cellp   = planeData.latVecOps['dparms']
    U0      = planeData.latVecOps['U0']
    wlen    = planeData.wavelength
    dsp     = planeData.getPlaneSpacings()
    fHKLs   = planeData.getSymHKLs()
    tThRng  = planeData.getTThRanges()
    symTag  = planeData.getLaueGroup()

    tThMin, tThMax = (r2d*tThRng.min(), r2d*tThRng.max()) # single range should be ok since entering hkls
    etaMin, etaMax = (0, 360)   # not sure when this will ever *NOT* be the case, so setting it

    omeMin = spotsArray.getOmegaMins()
    omeMax = spotsArray.getOmegaMaxs()

    omeRangeString = ''
    for iOme in range(len(omeMin)):
        if hasattr(omeMin[iOme], 'getVal'):
            omeRangeString += 'omegarange %g %g\n' % (omeMin[iOme].getVal('degrees'), omeMax[iOme].getVal('degrees'))
        else:
            omeRangeString += 'omegarange %g %g\n' % (omeMin[iOme] * r2d, omeMax[iOme] * r2d)

    # convert angles
    cellp[3:] = r2d*cellp[3:]

    # make the theoretical hkls string
    gvecHKLString = ''
    for i in range(len(dsp)):
        for j in range(fHKLs[i].shape[1]):
            gvecHKLString += '%1.8f %d %d %d\n' % (1/dsp[i], fHKLs[i][0, j], fHKLs[i][1, j], fHKLs[i][2, j])

    # now for the measured data section
    # xr yr zr xc yc ds eta omega
    gvecString = ''
    spotsIter = spotsArray.getIterPhase(phaseID, returnBothCoordTypes=True)
    for iSpot, angCOM, xyoCOM in spotsIter:
        sR, sC, sOme     = xyoCOM                          # detector coords
        sTTh, sEta, sOme = angCOM                          # angular coords (radians)
        sDsp = wlen / 2. / num.sin(0.5*sTTh)               # dspacing

        # get raw y, z (Fable frame)
        yraw = ncols_p - sC
        zraw = nrows_p - sR

        # convert eta to fable frame
        rEta = mapAngle(90. - r2d*sEta, [0, 360], units='degrees')

        # make mesaured G vector components in fable frame
        mGvec = makeMeasuredScatteringVectors(sTTh, sEta, sOme, convention='fable', frame='sample')

        # full Gvec components in fable lab frame (for grainspotter position fit)
        gveXYZ = spotsArray.detectorGeom.angToXYO(sTTh, sEta, sOme, outputGve=True)

        # no 4*pi
        mGvec = mGvec / sDsp

        # make contribution
        gvecString += '%1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %d %1.8f %1.8f %1.8f\n' \
                      % (mGvec[0], mGvec[1], mGvec[2], \
                         sR, sC,  \
                         1/sDsp, rEta, r2d*sOme, \
                         iSpot, \
                         gveXYZ[0, :], gveXYZ[1, :], gveXYZ[2, :])
        pass

    # write gve file for grainspotter
    fid = open(fileroot+'.gve', 'w')
    print >> fid, '%1.8f %1.8f %1.8f %1.8f %1.8f %1.8f ' % tuple(cellp)        + \
          cellString + '\n'                                                    + \
          '# wavelength = %1.8f\n'                       % (wlen)              + \
          '# wedge = 0.000000\n'                                               + \
          '# axis = %d %d %d\n'                          % tuple(osc_axis)     + \
          '# cell__a %1.4f\n'                            %(cellp[0])           + \
          '# cell__b %1.4f\n'                            %(cellp[1])           + \
          '# cell__c %1.4f\n'                            %(cellp[2])           + \
          '# cell_alpha %1.4f\n'                         %(cellp[3])           + \
          '# cell_beta  %1.4f\n'                         %(cellp[4])           + \
          '# cell_gamma %1.4f\n'                         %(cellp[5])           + \
          '# cell_lattice_[P,A,B,C,I,F,R] %s\n'          %(cellString)         + \
          '# chi 0.0\n'                                                        + \
          '# distance %.4f\n'                            %(wd_mu)              + \
          '# fit_tolerance 0.5\n'                                              + \
          '# o11  1\n'                                                         + \
          '# o12  0\n'                                                         + \
          '# o21  0\n'                                                         + \
          '# o22 -1\n'                                                         + \
          '# omegasign %1.1f\n'                          %(num.sign(deltaOme)) + \
          '# t_x 0\n'                                                          + \
          '# t_y 0\n'                                                          + \
          '# t_z 0\n'                                                          + \
          '# tilt_x 0.000000\n'                                                + \
          '# tilt_y 0.000000\n'                                                + \
          '# tilt_z 0.000000\n'                                                + \
          '# y_center %.6f\n'                            %(yc_p)               + \
          '# y_size %.6f\n'                              %(mmPerPixel*1.e3)    + \
          '# z_center %.6f\n'                            %(zc_p)               + \
          '# z_size %.6f\n'                              %(mmPerPixel*1.e3)    + \
          '# ds h k l\n'                                                       + \
          gvecHKLString                                                        + \
          '# xr yr zr xc yc ds eta omega\n'                                    + \
          gvecString
    fid.close()

    ###############################################################
    # GrainSpotter ini parameters
    #
    # fileroot = tempfile.mktemp()
    if positionFit:
        positionString = 'positionfit'
    else:
        positionString = '!positionfit'

    if numTrials == 0:
        randomString = '!random\n'
    else:
        randomString = 'random %g\n' % (numTrials)

    fid = open(fileroot+'_grainSpotter.ini', 'w')
    # self.__tempFNameList.append(fileroot)
    print >> fid, \
          'spacegroup %d\n' % (sgNum) + \
          'tthrange %g %g\n' % (tThMin, tThMax) + \
          'etarange %g %g\n' % (etaMin, etaMax) + \
          'domega %g\n' % (deltaOme) + \
          omeRangeString + \
          'filespecs %s.gve %s_grainSpotter.log\n' % (fileroot, fileroot) + \
          'cuts %d %g %g\n' % (minMeas, minCompl, minUniqn) + \
          'eulerstep %g\n' % (eulStep)+ \
          'uncertainties %g %g %g\n' % (uncertainty[0], uncertainty[1], uncertainty[2]) + \
          'nsigmas %d\n' % (nSigmas) + \
          'minfracg %g\n' % (minFracG) + \
          randomString + \
          positionString + '\n'
    fid.close()
    return
