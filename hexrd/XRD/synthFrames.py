"""
to do:
(*) make ReadGE take None for fileInfo -- 
	nFramesRemain = 0
	__convertFileInfo() to take None, and return []
	__nextFile call is skipped
	(*) or make a stripped-down reader class that does not actually do reading? call it a Frame class?

data:
    	centers (3)
        FWHMs   (3)
        scaling  (1)
... need to be able to set scaling to correspond to a given integrated intensity? -- for general normal distrib in 3D? 
... or just scale integrated intensity in vCalc.sum() once have it?!

"""

import detector
import spotFinder
import numpy as num

gParms = ( 2.06233014e+02,   
           2.08046997e+02,   
           1.05065037e+03,
           1.06710838e-03,  
          -2.93290955e-03,  
          -3.80295099e-04 )

dParms = (-7.29999194e-05,  
          -4.93935260e-04,   
           4.00000000e+00,
           2.00000000e+00 )

detectorGeom = detector.DetectorGeomGE(gParms, distortionParams=dParms)

tThMax = detectorGeom.getTThMax()

omegas = (num.linspace(1.,2.,2)-0.5-1.)*num.pi/180.

# inputs
spotParamsList = [
    [(0.5*tThMax,0.5*num.pi,omegas[1]), (0.002, 0.05, 0.05), 1e2],
    [(0.5*tThMax,0.0*num.pi,omegas[0]), (0.002, 0.05, 0.05), 1e2], # center, fwhm, A 
    [(0.5*tThMax,1.0*num.pi,omegas[0]), (0.002, 0.05, 0.05), 1e2],
    ]
intensityFunc = spotFinder.IntensityFuncGauss3D()
asSparse = True
outputPrefix = None
cutoffMult = 4.0 # multiplies FWHMs to get extent of coverage
debug=10

########################################
# def makeSynthFrames(spotParamsList, intensityFunc, detectorGeom, omegas, asSparse=True, outputPrefix=None, cutoffMult=4.0, debug=1)

import numpy as num
import spotFinder
haveSparse = True
try:
    from scipy import sparse
except:
    haveSparse = False

nFrames = len(omegas)

'might eventually want to add a check for constant delta-omega'
omegaDelta = num.mean(omegas[1:]-omegas[:-1])

# nParams = intensityFunc.getNParams(noBkg=False) # not needed

spotList = []
#
for iSpot, spotParams in enumerate(spotParamsList):
    xVec = intensityFunc.constructParams(*spotParams)

    'bbox from center and widths'
    'do not worry about excessive pixel coverage for spots at eta around 45 degrees and the like?'
    angCen = intensityFunc.getCenter(xVec)
    fwhm   = intensityFunc.getFWHM(xVec)
    angPM  = fwhm * cutoffMult
    xyfBBox = detectorGeom.angToXYOBBox(angCen, angPM, units='pixels', omegas=omegas)
    # xyfBBox_0 = num.array([sl[0] for sl in xyfBBox])
    # stack[slice(*xyfBBox[0]), slice(*xyfBBox[1]), slice(*xyfBBox[2])]
    
    'make spot instance, set up for just a single omega frame slice'
    xM, yM = num.meshgrid(num.arange(*xyfBBox[0]), num.arange(*xyfBBox[1]))
    xAll = xM.flatten(); yAll = yM.flatten();
    data = ( xAll, yAll, num.zeros_like(xAll), num.ones_like(xAll) )
    spot = spotFinder.Spot(iSpot, omegaDelta, data=data, detectorGeom=detectorGeom)
    spot.setupQuardFromFWHM(fwhm)
    
    spotList.append( (spot, xyfBBox, xVec) )
    if debug: print 'created spot %d'%(iSpot)

if outputPrefix is None:
    if asSparse:
        assert haveSparse,\
            'do not have sparse module to use'
        stack = [] # sparse.coo_matrix( (detectorGeom.nrows, detectorGeom.ncols), float )
    else:
        'this can eat up a big chunk of memory'
        stack = num.zeros( (nFrames, detectorGeom.nrows, detectorGeom.ncols), dtype=float)
else:
    if asSparse:
        raise NotImplementedError, 'have not implemented output to file in sparse format'
    stack = None
#
for iFrame, omega in enumerate(omegas):
    
    if outputPrefix is None:
        if asSparse:
            vThese = []
            xThese = []
            yThese = []
        else:
            frame = stack[iFrame,:,:]
    else:
        frame = detectorGeom.frame()
    
    for iSpot, (spot, xyfBBox, xVec) in enumerate(spotList):
        
        if iFrame < xyfBBox[2][0] or iFrame >= xyfBBox[2][1]:
            '>= for upper end because is meant to be used as a slice'
            if debug>2: print 'spot %d not on frame %d' % (iSpot, iFrame)
            continue
        
        'calculate intensities at the omega for this frame'
        spot.oAll[:] = omega
        vCalc = spot.getVCalc(intensityFunc, xVec, noBkg=True) 
        if debug>2:print vCalc.max()
        
        'put intensity on frames'
        if asSparse:
            vThese.append(vCalc)
            xThese.append(spot.xAll)
            yThese.append(spot.yAll)
        else:
            frame[spot.xAll, spot.yAll] += vCalc
    'done with spot loop'
    
    if outputPrefix is None:
        if asSparse:
            if len(vThese) > 0:
                vThese = num.hstack(vThese)
                xThese = num.hstack(xThese)
                yThese = num.hstack(yThese)
                if debug>1: print 'frame %d will have up to %d nonzeros' % (iFrame,len(vThese))
                frame = sparse.coo_matrix( ( vThese , (xThese,yThese) ) , 
                                           shape=(detectorGeom.nrows, detectorGeom.ncols) )
                frame = frame.tocsr() # note that coo->csr sums entries as desired for overlapped spots
            else:
                'just make an empty matrix'
                frame = sparse.csr_matrix( (detectorGeom.nrows, detectorGeom.ncols) )
            stack.append(frame)
        else:
            'nothing to do here'
            pass
    else:
        frame.tofile(outputPrefix+'_%04d.dat'%iFrame)

    if debug: print 'created frame %d'%(iFrame)

'done with loop over frames'
    
# return stack # may be None    

a = stack[-1].todense()
detectorGeom.display(a, range=(0,a.max()))
