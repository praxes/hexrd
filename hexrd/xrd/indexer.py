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
import ctypes
import tempfile
import glob
import time

import numpy as num

import hexrd.matrixutil as mUtil
import hexrd.xrd.grain
from hexrd.xrd.grain import makeMeasuredScatteringVectors
import hexrd.xrd.rotations
from hexrd.xrd.rotations import mapAngle
from hexrd.xrd.symmetry import toFundamentalRegion
from hexrd.xrd import xrdBase

if xrdBase.haveMultiProc:
    multiprocessing = xrdBase.multiprocessing # formerly import

# module vars
piby2 = num.pi * 0.5
r2d = 180. / num.pi
d2r = num.pi / 180.

Xl = num.c_[1, 0, 0].T
Yl = num.c_[0, 1, 0].T
Zl = num.c_[0, 0, 1].T

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
        A word on spacegroup numbers: it appears that grainspotter is using the 'VolA' tag for calls to SgInfo
        """
        location = self.__class__.__name__
        tic = time.time()

        # keyword argument processing
        phaseID     = None
        gVecFName   = 'tmpGve'
        sgNum       = 225
        cellString  = 'F'
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
        positionFit = False

        kwarglen = len(kwargs)
        if kwarglen > 0:
            argkeys = kwargs.keys()
            for i in range(kwarglen):
                if argkeys[i] == 'sgNum':
                    sgNum = kwargs[argkeys[i]]
                elif argkeys[i] == 'phaseID':
                    phaseID = kwargs[argkeys[i]]
                elif argkeys[i] == 'gVecFName':
                    gVecFName = kwargs[argkeys[i]]
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

        # cleanup stuff from any previous run
        self.cleanup()

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
        ii = 0
        spotsIter = spotsArray.getIterPhase(phaseID, returnBothCoordTypes=True)
        for iSpot, angCOM, xyoCOM in spotsIter:
            sX, sY, sOme     = xyoCOM                          # detector coords
            sTTh, sEta, sOme = angCOM                          # angular coords (radians)
            sDsp = wlen / 2. / num.sin(0.5*sTTh)               # dspacing

            # convert eta to risoe frame
            rEta = mapAngle(90. - r2d*sEta, [0, 360], units='degrees')

            # make mesaured G vector components in risoe frame
            mGvec = makeMeasuredScatteringVectors(sTTh, sEta, sOme, convention='risoe', frame='sample')
            # mQvec = makeMeasuredScatteringVectors(sTTh, sEta, sOme, convention='llnl', frame='lab')
            gveXYZ = spotsArray.detectorGeom.angToXYO(sTTh, sEta, sOme, outputGve=True)

            mGvec = mGvec / sDsp
            # mQvec = 4*num.pi*num.sin(0.5*sTTh)*mQvec/wlen

            # make contribution
            gvecString += '%1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %d %1.8f %1.8f %1.8f\n' \
                          % (mGvec[0], mGvec[1], mGvec[2], \
                             sX, sY, \
                             1/sDsp, rEta, r2d*sOme, \
                             iSpot, \
                             gveXYZ[0, :], gveXYZ[1, :], gveXYZ[2, :])

            # advance counter
            ii += 1

        # write gve file for grainspotter
        f = open(gVecFName+'.gve', 'w')
        print >> f, '%1.8f %1.8f %1.8f %1.8f %1.8f %1.8f ' % tuple(cellp) + \
              cellString + '\n' + \
              '# wavelength = %1.8f\n' % (wlen) + \
              '# wedge = 0.000000\n# ds h k l\n' + \
              gvecHKLString + \
              '# xr yr zr xc yc ds eta omega\n' + \
              gvecString
        f.close()

        ###############################################################
        # GrainSpotter ini parameters
        #
        # tempFNameIn = tempfile.mktemp()
        if positionFit:
            positionString = 'positionfit'
        else:
            positionString = '!positionfit'

        if numTrials == 0:
            randomString = '!random\n'
        else:
            randomString = 'random %g\n' % (numTrials)

        tempFNameIn = 'tmpIni'
        f = open(tempFNameIn, 'w')
        # self.__tempFNameList.append(tempFNameIn)
        print >> f, \
              'spacegroup %d\n' % (sgNum) + \
              'tthrange %g %g\n' % (tThMin, tThMax) + \
              'etarange %g %g\n' % (etaMin, etaMax) + \
              'domega %g\n' % (deltaOme) + \
              omeRangeString + \
              'filespecs %s.gve %s.log\n' % (gVecFName, tempFNameIn) + \
              'cuts %d %g %g\n' % (minMeas, minCompl, minUniqn) + \
              'eulerstep %g\n' % (eulStep)+ \
              'uncertainties %g %g %g\n' % (uncertainty[0], uncertainty[1], uncertainty[2]) + \
              'nsigmas %d\n' % (nSigmas) + \
              'minfracg %g\n' % (minFracG) + \
              randomString + \
              positionString + '\n'
        f.close()

        toc = time.time()
        print 'in %s, setup took %g' % (location, toc-tic)
        tic = time.time()

        # tempFNameStdout = tempfile.mktemp()
        # self.__tempFNameList.append(tempFNameStdout)
        tempFNameStdout = 'tmp.out'
        # grainSpotterCmd = '%s %s > %s' % (self.__execName, tempFNameIn, tempFNameStdout)
        grainSpotterCmd = '%s %s' % (self.__execName, tempFNameIn)
        os.system(grainSpotterCmd)
        toc = time.time()
        print 'in %s, execution took %g' % (location, toc-tic)
        tic = time.time()

        # add output files to cleanup list
        # self.__tempFNameList += glob.glob(tempFNameIn+'.*')

        # collect data from gff file'
        gffFile = tempFNameIn+'.gff'
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

def convertUToRotMat(Urows, U0, symTag='Oh'):
    """
    Takes GrainSpotter gff ouput

    U11 U12 U13 U21 U22 U23 U13 U23 U33

    and takes it into the LLNL/APS frame of reference

    Urows comes from grainspotter's gff output
    U0 comes from xrd.crystallography.latticeVectors.U0
    """
    R = hexrd.xrd.rotations

    numU, testDim = Urows.shape
    assert testDim == 9, "Your input must have 9 columns; yours has %d" % (testDim)

    Rsamp = num.dot( R.rotMatOfExpMap(piby2*Zl), R.rotMatOfExpMap(piby2*Yl) )
    qin  = R.quatOfRotMat(Urows.reshape(numU, 3, 3))
    print "quaternions in (Risoe convention):"
    print qin.T
    qout = num.dot( R.quatProductMatrix( R.quatOfRotMat(Rsamp), mult='left' ), \
                    num.dot( R.quatProductMatrix( R.quatOfRotMat(U0.T), mult='right'),  \
                             qin ).squeeze() ).squeeze()
    if qout.ndim == 1:
        qout = toFundamentalRegion(qout.reshape(4, 1), crysSym=symTag, sampSym=None)
    else:
        qout = toFundamentalRegion(qout, crysSym=symTag, sampSym=None)
    print "quaternions out (LLNL convention, symmetrically reduced)"
    print qout.T
    Uout = R.rotMatOfQuat(qout)
    return Uout

def convertRotMatToRisoeU(rMats, U0, symTag='Oh'):
    """
    Makes GrainSpotter gff ouput

    U11 U12 U13 U21 U22 U23 U13 U23 U33

    and takes it into the LLNL/APS frame of reference

    Urows comes from grainspotter's gff output
    U0 comes from xrd.crystallography.latticeVectors.U0
    """
    R = hexrd.xrd.rotations # formerly import

    numU = num.shape(num.atleast_3d(rMats))[0]

    Rsamp = num.dot( R.rotMatOfExpMap(piby2*Zl), R.rotMatOfExpMap(piby2*Yl) )
    qin  = R.quatOfRotMat(num.atleast_3d(rMats))
    print "quaternions in (LLNL convention):"
    print qin.T
    qout = num.dot( R.quatProductMatrix( R.quatOfRotMat(Rsamp.T), mult='left' ), \
                    num.dot( R.quatProductMatrix( R.quatOfRotMat(U0), mult='right'),  \
                             qin ).squeeze() ).squeeze()
    if qout.ndim == 1:
        qout = toFundamentalRegion(qout.reshape(4, 1), crysSym=symTag, sampSym=None)
    else:
        qout = toFundamentalRegion(qout, crysSym=symTag, sampSym=None)
    print "quaternions out (Risoe convention, symmetrically reduced)"
    print qout.T
    Uout = R.rotMatOfQuat(qout)
    return Uout

######################################################################

"""
things for doing fiberSearch with multiprocessing;
multiprocessing has a hard time pickling a function defined in the local scope of another function,
so stuck putting the function out here;
"""
debugMultiproc = 1
if xrdBase.haveMultiProc:
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
    (*) doFit is not done here -- in multiprocessing, that would end up happening on a remote process
    and then different processes would have different data, unless spotsArray were made to be fancier
    """
    G = hexrd.xrd.grain
    R = hexrd.xrd.rotations

    """
    kludge stuff so that this function is outside of fiberSearch
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
    thisRMat = R.rotMatOfQuat(thisQ)

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
            print ppfx+"testing candidate..."+str(thisQ)
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
        cInfo = R.quatOfRotMat(candidate.rMat).flatten().tolist()
        cInfo.append(candidate.completeness)
        print ppfx+"Grain found at q = [%g  %g  %g  %g] with completeness %g" \
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
    G = hexrd.xrd.grain
    R = hexrd.xrd.rotations
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

    candidate = G.Grain(spotsArray, rMat=None,
                        etaTol=etaTol, omeTol=omeTol)
    multiProcMode = xrdBase.haveMultiProc and doMultiProc
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
        nCPUs = nCPUs or xrdBase.dfltNCPU
        spotsArray.multiprocMode = True
        pool = multiprocessing.Pool(nCPUs)

    """
    HKL ITERATOR
    """
    numTotal = len(spotsArray)
    pctClaimed = 0.
    tic = time.time()
    for iHKL in range(nHKLs):
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

                Gplus  = G.makeMeasuredScatteringVectors(*angs_I)
                Gminus = G.makeMeasuredScatteringVectors(*angs_J)

                Gvec = 0.5*(Gplus - Gminus)
                maxSpots = 0.5*(sum(thisRingSpots) - sum(spotsArray.friedelPair[thisRingSpots] == -1))
            else:
                iSpot, angs_I = stuff
                Gvec  = G.makeMeasuredScatteringVectors(*angs_I)
                maxSpots = sum(thisRingSpots)
            print "\nProcessing reflection %d (spot %d), %d remain unclaimed\n" % (iRefl+1, iSpot, maxSpots)
            if multiProcMode and debugMultiproc > 1:
                marks = spotsArray._Spots__marks[:]
                print 'marks : '+str(marks)
            # make the fiber;
            qfib = R.discreteFiber(hklList[iHKL], Gvec,
                                   B=bMat,
                                   ndiv=nsteps,
                                   invert=False,
                                   csym=csym, ssym=None)[0]
            # if +/- hkl aren't in the symmetry group, need '-' fiber
            if not centroSymRefl[thisHKLID]:
                minusHKL = -num.r_[hklList[iHKL]]
                qfibM = R.discreteFiber(minusHKL, Gvec,
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
                print "No grain found containing reflection %d\n" % (iSpot)
                # import pdb;pdb.set_trace()
            else:
                asMaster = multiProcMode
                'sort based on completeness'
                trialGrainCompletenesses = [tgd['completeness'] for tgd in trialGrains]
                order = num.argsort(trialGrainCompletenesses)[-1::-1]
                for iTrialGrain in order:
                    foundGrainData = trialGrains[iTrialGrain]
                    foundGrain = G.Grain(spotsArray, grainData=foundGrainData, claimingSpots=False)
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
                if pctClaimed > minPctClaimed:
                    break
                if quitAfter is not None:
                    if nGrains >= quitAfter:
                        break
        'end of iRefl loop'
        if pctClaimed > minPctClaimed:
            break
    'end of iHKL loop'
    rMats = num.empty((len(grainList), 3, 3))
    for i in range(len(grainList)):
        rMats[i, :, :] = grainList[i].rMat

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

    return rMats

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

def paintGrid(quats, etaOmeMaps, threshold=None, bMat=None, omegaRange=None, etaRange=None, progressBar=False, doMultiProc=False, nCPUs=None, debug=False):
    """
    do a direct search of omega-eta maps to paint each orientation in quats with a completeness

    bMat is in CRYSTAL frame

    etaOmeMaps is instance of xrd.xrdUtils.CollapseOmeEta

    omegaRange=([-num.pi/3., num.pi/3.],) for example

    *) lifted mainly from grain.py

    *) self.etaGrid, self.omeGrid = num.meshgrid(self.etaEdges, self.omeEdges)
       this means that ETA VARIES FASTEST!

    ...make a new function that gets called by grain to do the g-vec angle computation?
    """
    rot = hexrd.xrd.rotations

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

    'index munging here -- look away'
    # mapIndices = num.indices([numEtas, numOmes])
    # etaIndices = mapIndices[0].flatten()
    # omeIndices = mapIndices[1].T.flatten()
    etaIndices = num.tile(range(numEtas), (numOmes))
    omeIndices = num.tile(range(numOmes), (numEtas))

    omeMin = None
    omeMax = None
    if omegaRange is not None:
        omeMin = [omegaRange[i][0] for i in range(len(omegaRange))]
        omeMax = [omegaRange[i][1] for i in range(len(omegaRange))]
    etaMin = None
    etaMax = None
    if etaRange is not None:
        etaMin = [etagaRange[i][0] for i in range(len(etagaRange))]
        etaMax = [etagaRange[i][1] for i in range(len(etagaRange))]

    # make list of rMats from input quats
    rMatsList = [rot.rotMatOfQuat(quats[:, i]) for i in range(quats.shape[1])]

    multiProcMode = xrdBase.haveMultiProc and doMultiProc

    if multiProcMode:
        nCPUs = nCPUs or xrdBase.dfltNCPU
        print "INFO: using multiprocessing with %d processes\n" % (nCPUs)
    else:
        print "INFO: running in serial mode\n"
        nCPUs = 1

    # assign the globals for paintGridThis
    global planeData_MP
    global omeMin_MP, omeMax_MP
    global etaMin_MP, etaMax_MP
    global omeIndices_MP, etaIndices_MP
    global hklList_MP, hklIDs_MP
    global etaOmeMaps_MP
    global bMat_MP
    global threshold_MP
    planeData_MP  = planeData
    hklIDs_MP     = hklIDs
    hklList_MP    = hklList
    omeMin_MP     = omeMin
    omeMax_MP     = omeMax
    omeIndices_MP = omeIndices
    etaMin_MP     = etaMin
    etaMax_MP     = etaMax
    etaIndices_MP = etaIndices
    etaOmeMaps_MP = etaOmeMaps
    bMat_MP       = bMat
    threshold_MP  = threshold

    # do the mapping
    retval = None
    if multiProcMode:
        pool = multiprocessing.Pool(nCPUs)
        retval = pool.map(paintGridThis, rMatsList)
    else:
        retval = map(paintGridThis, rMatsList)

    planedata_mp  = None
    hklIDs_MP     = None
    hklList_MP    = None
    omeMin_MP     = None
    omeMax_MP     = None
    omeIndices_MP = None
    etaMin_MP     = None
    etaMax_MP     = None
    etaIndices_MP = None
    etaOmeMaps_MP = None
    bMat_MP       = None
    threshold_MP  = None

    if multiProcMode:
        pool.close()

    return retval

def paintGridThis(rMat):
    """
    """
    # pull local vars from globals set in paintGrid
    global planeData_MP
    global omeMin_MP, omeMax_MP
    global etaMin_MP, etaMax_MP
    global omeIndices_MP, etaIndices_MP
    global hklList_MP, hklIDs_MP
    global etaOmeMaps_MP
    global bMat_MP
    global threshold_MP
    planeData  = planeData_MP
    hklIDs     = hklIDs_MP
    hklList    = hklList_MP
    omeMin     = omeMin_MP
    omeMax     = omeMax_MP
    omeIndices = omeIndices_MP
    etaMin     = etaMin_MP
    etaMax     = etaMax_MP
    etaIndices = etaIndices_MP
    etaOmeMaps = etaOmeMaps_MP
    bMat       = bMat_MP
    threshold  = threshold_MP

    debug = False

    nHKLS     = len(hklIDs)

    nPredRefl = 0
    nMeasRefl = 0
    reflInfoList = []
    dummySpotInfo = num.nan * num.ones(3)

    hklCounterP = 0                 # running count of excpected (predicted) HKLs
    hklCounterM = 0                 # running count of "hit" HKLs
    for iHKL in range(nHKLS):
        # for control of tolerancing
        symHKLs   = planeData.getSymHKLs()[hklIDs[iHKL]]
        tThRanges = planeData.getTThRanges()[hklIDs[iHKL]]

        # make all theoretical scattering vectors
        predQvec, predQAng0, predQAng1 = \
                  planeData.makeTheseScatteringVectors([hklList[iHKL]], rMat, bMat=bMat)

        # work with generated spots for iHKL
        # allPredAng = zip(predQAng0[iHKL], predQAng1[iHKL])
        allPredAng = zip(predQAng0, predQAng1)

        # filter using omega range
        if omeMin is not None:
            reflInRangeOme0 = num.zeros(allPredAng[2][0].shape, dtype=bool)
            reflInRangeOme1 = num.zeros(allPredAng[2][1].shape, dtype=bool)
            for iOme in range(len(omeMin)):
                reflInRangeOme0 = reflInRangeOme0 | (
                    (allPredAng[2][0] >= omeMin[iOme]) &
                    (allPredAng[2][0] <= omeMax[iOme]) )
                reflInRangeOme1 = reflInRangeOme1 | (
                    (allPredAng[2][1] >= omeMin[iOme]) &
                    (allPredAng[2][1] <= omeMax[iOme]) )
                pass
        else:
            reflInRangeOme0 = num.ones(allPredAng[2][0].shape, dtype=bool)
            reflInRangeOme1 = num.ones(allPredAng[2][1].shape, dtype=bool)
            pass

        if etaMin is not None:
            reflInRangeEta0 = num.zeros(allPredAng[2][0].shape, dtype=bool)
            reflInRangeEta1 = num.zeros(allPredAng[2][1].shape, dtype=bool)
            for iEta in range(len(etaMin)):
                reflInRangeEta0 = reflInRangeEta0 | (
                    (allPredAng[1][0] >= etaMin[iEta]) &
                    (allPredAng[1][0] <= etaMax[iEta]) )
                reflInRangeEta1 = reflInRangeEta1 | (
                    (allPredAng[1][1] >= etaMin[iEta]) &
                    (allPredAng[1][1] <= etaMax[iEta]) )
                pass
        else:
            reflInRangeEta0 = num.ones(allPredAng[2][0].shape, dtype=bool)
            reflInRangeEta1 = num.ones(allPredAng[2][1].shape, dtype=bool)
            pass

        reflInRange0 = reflInRangeOme0 & reflInRangeEta0
        reflInRange1 = reflInRangeOme1 & reflInRangeEta1

        # get culled angle and hkl lists for predicted spots
        culledTTh = num.r_[ allPredAng[0][0][reflInRange0], allPredAng[0][1][reflInRange1] ]
        culledEta = num.r_[ allPredAng[1][0][reflInRange0], allPredAng[1][1][reflInRange1] ]
        culledOme = num.r_[ allPredAng[2][0][reflInRange0], allPredAng[2][1][reflInRange1] ]

        culledHKLs = num.hstack( [
            symHKLs[:, reflInRange0],
            symHKLs[:, reflInRange1] ] )

        culledQvec = num.c_[ predQvec[:, reflInRange0], predQvec[:, reflInRange1] ]

        nThisPredRefl = len(culledTTh)
        hklCounterP += nThisPredRefl
        for iTTh in range(nThisPredRefl):
            culledEtaIdx = num.where(etaOmeMaps.etaEdges - culledEta[iTTh] > 0)[0]
            if len(culledEtaIdx) > 0:
                culledEtaIdx = culledEtaIdx[0] - 1
                if culledEtaIdx < 0:
                    culledEtaIdx = None
            else:
                culledEtaIdx = None
            culledOmeIdx = num.where(etaOmeMaps.omeEdges - culledOme[iTTh] > 0)[0]
            if len(culledOmeIdx) > 0:
                culledOmeIdx = culledOmeIdx[0] - 1
                if culledOmeIdx < 0:
                    culledOmeIdx = None
            else:
                culledOmeIdx = None

            if culledEtaIdx is not None and culledOmeIdx is not None:
                pixelVal = etaOmeMaps.dataStore[iHKL][omeIndices[culledOmeIdx], etaIndices[culledEtaIdx] ]
                isHit = pixelVal >= threshold[iHKL]
                if isHit:
                    hklCounterM += 1
                    pass
                pass
            if debug:
                print "hkl %d -->\t%d\t%d\t%d\t" % (iHKL, culledHKLs[0, iTTh], culledHKLs[1, iTTh], culledHKLs[2, iTTh]) + \
                      "isHit=%d\tpixel value: %g\n" % (isHit, pixelVal) + \
                      "eta index: %d,%d\tetaP: %g\n" % (culledEtaIdx, etaIndices[culledEtaIdx], r2d*culledEta[iTTh]) + \
                      "ome index: %d,%d\tomeP: %g\n" % (culledOmeIdx, omeIndices[culledOmeIdx], r2d*culledOme[iTTh])
                pass
            pass # close loop on signed reflections
        pass # close loop on HKL
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
