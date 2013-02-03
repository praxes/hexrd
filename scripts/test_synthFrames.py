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

usage = '''\
For running tests using synthetic frames

Examples:

test_synthFrames.py --test=makeSynthStack --test=stackToSpots --test=index

test_synthFrames.py --test=makeSynthStack --test=stackToSpots --num-grains=25

test_synthFrames.py --test=makeSynthStack --test=stackToSpots --num-grains=2 --args-in-degrees=True --ome-min=60 --ome-max=-60 --n-ome=120

test_synthFrames.py --detector-geom='generic;[];{"ncols":512,"nrows":512,"pixelPitch":0.8}' --fwhm='(0.002, 0.01, 0.01)' --test=makeSynthStack --test=stackToSpots --test=index
'''

import os
import time

import numpy as num
import optparse

from hexrd import valunits
from hexrd import orientations as ors

from hexrd.xrd            import grain, xrdutil, detector, spotfinder, indexer
from hexrd.xrd.spotfinder import Spots
from hexrd.xrd.material   import Material, loadMaterialList

# from examples.rubyNIST import planeData
#
## from examples import fe
## planeData = fe.getAlphaPD()
mat_list = loadMaterialList("../hexrd/data/all_materials.cfg")

planeData = mat_list[2].planeData
planeData.exclusions = num.zeros( len(planeData.exclusions), dtype=bool )
planeData.tThMax = 0.17

USAGE = '%s [options]' % __file__ + '\n\n' + usage

'defaults'
degToRad = num.pi/180.
omeMin =  -90.*degToRad
omeMax =   90.*degToRad
nOme   =   360

synthStackName = 'synth.stack'
synthSpotsName = 'synth.spots'

debug = 1

def createParser(usage=USAGE):
    """
    Returns:  optparse.OptionParser object
    
    Creates an OptionParser for the command line options.
    """
    parser = optparse.OptionParser(usage=USAGE) # , version=version

    parser.set_defaults(tests=[])
    parser.add_option (
       "-t", "--test", action = "append", dest = "tests", metavar = "TEST",
       help = "add TEST to list of tests to run [default=%default]")

    parser.set_defaults(inputsInDeg=False)
    parser.add_option(
        "--args-in-degrees", action = "store_true", dest = "inputsInDeg",
        help = "take command-line arguments in degrees [default=%default]")

    parser.set_defaults(doMultiProc=False)
    parser.add_option(
        "--multiproc", action = "store_true", dest = "doMultiProc",
        help = "run with multiprocessing [default=%default]")

    parser.set_defaults(numGrains=1)
    parser.add_option(
        "-g", "--num-grains", action = "store", type = "int",
        dest = "numGrains", 
        help = "number of grains to use [default=%default]")

    parser.set_defaults(detectorGeomStr='ge;[];{}')
    parser.add_option (
        "-d", "--detector-geom", action = "store", type = "string",
        dest = "detectorGeomStr", 
        help = "detector [default=%default]")

    parser.set_defaults(fwhmStr='(0.005, 0.01, 0.01)')
    parser.add_option (
        "--fwhm", action = "store", type = "string",
        dest = "fwhmStr", 
        help = "fwhm [default=%default]")

    parser.set_defaults(omeMin=None)
    parser.add_option (
        "--ome-min", action = "store", type = "float",
        dest = "omeMin", 
        help = "minimum omega (not center of bin) [default=%default]")

    parser.set_defaults(omeMax=None)
    parser.add_option (
        "--ome-max", action = "store", type = "float",
        dest = "omeMax", 
        help = "maximum omega (not center of bin) [default=%default]")

    parser.set_defaults(nOme=None)
    parser.add_option (
        "-n", "--n-ome", action = "store", type = "int",
        dest = "nOme", 
        help = "number of omega frames [default=%default]")

    parser.set_defaults(randomSeed=0)
    parser.add_option (
        "--seed", action = "store", type = "int",
        dest = "randomSeed", 
        help = "random number seed [default=%default]")
    
    return parser

######################################################################

def makeRandGrains(numGrains, omegaMM, detectorGeom):
    
    spotsDummy = Spots(planeData, None, detectorGeom, omegaMM)

    if numGrains == 1:
        rMats = num.atleast_3d(num.eye(3)).T
    else:
        quats = num.atleast_2d(ors.Quat.getRandQuat(n=numGrains)).T
        rMats = ors.quatToMat(quats)
        
    gList = []
    for rMat in rMats:
        g = grain.Grain(spotsDummy, rMat=rMat)
        gList.append(g)
        
    return gList

def getSpotsFromGrains(gList, fwhm, a):
    '''
    see also makeSynthSpots in xrdutil for another approach 
    
    eventually, would want to fold in structure factors and the like;
    perhaps scaling factors as is done for pole figures to detector angular space
    '''

    spotParamsList = []

    if hasattr(a,'__len__'):
        'assume fwhm and a specified on a grain-by-grain basis'
        for g, fwhmThis, aThis in zip(gList, fwhm, a):
            for angs in g.grainSpots['predAngles'] :
                spotParamsList.append( [angs, fwhmThis, aThis] )
    else:
        for g in gList:
            for angs in g.grainSpots['predAngles'] :
                spotParamsList.append( [angs, fwhm, a] )
    
    return spotParamsList

######################################################################

if __name__ == '__main__':
    import sys
    # main(sys.argv[1:])
    argv = sys.argv[1:]
    # def main(argv):
    
    parser = createParser()
    (opts, args) = parser.parse_args(argv)
    
    ######################################################################
    
    angConv = 1.0
    if opts.inputsInDeg:
        angConv = degToRad
    #
    if opts.omeMin is not None:
        omeMin = opts.omeMin * angConv
    if opts.omeMax is not None:
        omeMax = opts.omeMax * angConv
    #
    if opts.nOme is not None:
        nOme = opts.nOme

    omegaMM = (omeMin, omeMax)
    deltaOme = (omeMax-omeMin)/nOme # may be negative
    omegas = num.linspace(omeMin + 0.5*deltaOme, omeMax - 0.5*deltaOme, nOme)

    num.random.seed(seed=opts.randomSeed)
    
    ######################################################################
    
    # ... add generic detector capability, so that can have smaller number of pixels -- detector.newGenericDetector
    # ... add detector geom parameters
    #detectorGeom = detector.newDetector('ge')
    dgType, argsStr, kwargsStr = opts.detectorGeomStr.split(';')
    args = eval(argsStr); kwargs = eval(kwargsStr);
    kwargs.setdefault('readerKWArgs',{"doFlip":False})
    detectorGeom = detector.newDetector(dgType, *args, **kwargs)
    detectorGeom.workDist = 2000.
    tThMax = detectorGeom.getTThMax(func=num.max)
    print tThMax/degToRad
    planeData.tThMax = tThMax

    ######################################################################
    
    if opts.tests.count("makeSynthStack") > 0:
        tic = time.time()
        
        gList = makeRandGrains(opts.numGrains, omegaMM, detectorGeom)
        '''    
        spotParamsList = [
                [(0.5*tThMax,0.5*num.pi,omegas[1]), (0.002, 0.05, 0.05), 1e2], # center, fwhm, A 
                [(0.5*tThMax,0.0*num.pi,omegas[0]), (0.002, 0.05, 0.05), 1e2],
                [(0.5*tThMax,1.0*num.pi,omegas[0]), (0.002, 0.05, 0.05), 1e2],
                ]
        '''
        
        fwhm = eval(opts.fwhmStr)
        a = 1e4 # num.iinfo('uint16').max = 65535
        if debug > 1:
            for g in gList:
                print g.grainSpots['predAngles']
        spotParamsList = getSpotsFromGrains(gList, fwhm, a)
    
        if os.path.exists(synthStackName):
            os.remove(synthStackName)
    
        # writer = detector.ReadGE(None, omeMin, deltaOme, doFlip=False).getWriter(synthStackName)
        writer = detectorGeom.reader.getWriter(synthStackName)
        stack = xrdutil.makeSynthFrames(spotParamsList, detectorGeom, omegas, 
                                        cutoffMult=3.0, # to keep things less computationally expensive
                                        output=writer, debug=debug)
        writer.close()
        '''
        # if a stack of sparse matrices is returned, can do things like:
        frame = stack[-1].todense()
        pw = detectorGeom.display(frame)
        detectorGeom.drawRings(pw, planeData)
        '''
        
        '''
        # readerTmp = reader.makeNew()
        readerTmp = detectorGeom.getNewReader(synthStackName, omeMin, deltaOme, doFlip=False, subtractDark=False)
        flatFrame = readerTmp.read(nframes=nOme, sumImg=num.maximum)
        pw = detectorGeom.display(flatFrame)
        '''
        
        toc = time.time()
        print '[time] took %g seconds to make the stack' % (toc-tic)
        
    # ... add noise -- normal with mean and stddev
    # ... add background -- normal with mean and stddev
    # ... add rings
    
    ######################################################################
    
    if opts.tests.count("stackToSpots") > 0:
        'now pull the spots back out'

        if not os.path.exists(synthStackName):
            raise RuntimeError, 'no stack, try running makeSynthStack'
    
        tic = time.time()

        # reader = detector.ReadGE(synthStackName, omeMin, deltaOme, doFlip=False, subtractDark=False)
        reader = detectorGeom.getNewReader(synthStackName, omeMin, deltaOme, doFlip=False, subtractDark=False)
        assert reader.getNFrames() == nOme, 'not the expected number of frames'
        
        threshold = 20
        minPx = 4
        spots = spotfinder.Spots.findSpotsOmegaStack(
            reader, 0, 
            threshold, minPx, debug=True, 
            overlapPixelDistance=None,
            discardAtBounds=True)
        print 'found %d spots in the image stack' % (len(spots))
        toc = time.time()
        print '[time] took %g seconds to run findSpotsOmegaStack' % (toc-tic)
        
        if os.path.exists(synthSpotsName):
            os.remove(synthSpotsName)
    
        import glob
        fnames = glob.glob(synthSpotsName+'*') 
        for filename in fnames:
            os.remove(filename)
        spotfinder.Spot.storeSpots(synthSpotsName, spots)

    ######################################################################

    if opts.tests.count("index") > 0:
        'index and fit'

        'use glob because filename can end up with a suffix'
        import glob
        fnames = glob.glob(synthSpotsName+'*') 
        if not len(fnames) > 0 : # if not os.path.exists(synthSpotsName):
            raise RuntimeError, 'no spots, try running stackToSpots'

        spots, dgTemp = spotfinder.Spot.loadSpots(synthSpotsName)
        
        tThMax = detectorGeom.getTThMax(func=num.min)
        planeData.tThMax = tThMax
        
        spotsForFit = spotfinder.Spots(planeData, spots, detectorGeom, omegaMM)
    
        if debug > 1:
            foundSpotAngs = []
            for spot in spots:
                foundSpotAngs.append(spot.angCOM())
            foundSpotAngs = num.array(foundSpotAngs)
            print foundSpotAngs
            '''
            from hexrd import plotwrap
            pwA = plotwrap.PlotWrap()
            pwA(predAngles[:,0],    predAngles[:,1],    style='r+')
            pwA(foundSpotAngs[:,0], foundSpotAngs[:,1], style='gx')
            # alignment is good, but definitely missing somex
            pwB = plotwrap.PlotWrap()
            pwB(predAngles[:,0],    predAngles[:,2],    style='r+')
            pwB(foundSpotAngs[:,0], foundSpotAngs[:,2], style='gx')
            '''
    
        
        nSpotsTThAll, nSpotsTThClaimed, numClaimed, numMultiClaimed = spotsForFit.report(scaleByMultip=False)
        # iHKL = num.argmax(nSpotsTThAll)
        iHKL = num.where(nSpotsTThAll == num.unique(nSpotsTThAll)[1])[0][0]
        hklListForIndex = [ tuple(planeData.getHKLs()[iHKL]), ]
        
        etaTol = valunits.valWUnit('etaTol', 'angle', 0.5, 'degrees')
        omeTol = valunits.valWUnit('omeTol', 'angle', 0.5, 'degrees')
        rMat = indexer.fiberSearch(spotsForFit, 
                                   hklListForIndex,
                                   iPhase=planeData.phaseID,
                                   nsteps=720,
                                   minCompleteness=0.70,
                                   minPctClaimed=0.75,
                                   friedelOnly=False,
                                   dspTol=None,
                                   etaTol=etaTol,
                                   omeTol=omeTol,
                                   debug=True,
                                   doMultiProc=opts.doMultiProc,
                                   doRefinement=False)
        print 'found %d grains' % (rMat.shape[0])
    
        etaTol    = valunits.valWUnit("etaTol","ANGLE",0.0625,"degrees")
        omeTol    = valunits.valWUnit("omeTol","ANGLE",0.2500,"degrees")
        strainMag = 1e-3
        firstGrain = grain.Grain(spotsForFit,
                                 rMat=rMat[0, :, :],
                                 etaTol=etaTol,
                                 omeTol=omeTol)
        
        firstGrain.fit()

    ######################################################################
    # add more tests here as see fit
    
    # return
