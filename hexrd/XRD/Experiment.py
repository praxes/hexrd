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
#
######################################################################
## TOP-LEVEL MODULES AND SOME GLOBALS
##
"""Module for wrapping the main functionality of the XRD package.

The Experiment class is the primary interface.  Other classes 
are helpers.
"""
import sys, os, copy
import cPickle

import numpy

from hexrd.XRD import detector
from hexrd.XRD import spotFinder
from hexrd.XRD import crystallography
from hexrd.XRD import grain as G
from hexrd.XRD import indexer

from hexrd.XRD.hydra    import Hydra
from hexrd.XRD.Material import Material, loadMaterialList


from hexrd import valUnits
from hexrd.XRD import xrdUtils

#
#  Defaults (will eventually make to a config file)
#
HERE = os.path.dirname(__file__)
toMatFile    = os.path.join(HERE, '..', 'Data', 'materials.cfg')
DFLT_MATFILE = os.path.normpath(toMatFile) # check whether it exists
matfileOK    = os.access(DFLT_MATFILE, os.F_OK)
if not matfileOK:  # use relative path
    DFLT_MATFILE = os.path.join('Data', 'materials.cfg')
    pass
matfileOK = os.access(DFLT_MATFILE, os.F_OK)
if not matfileOK:  # set to null
    DFLT_MATFILE = ''
    pass

print 'default matfile:  ', DFLT_MATFILE.join(2*['"'])
print 'current dir:  ', os.getcwd()
#
#
__all__ = ['Experiment',
           'FitModes', 'ImageModes',
           'ReaderInput', 'CalibrationInput', 'PolarRebinOpts',
           'saveExp', 'loadExp']
#
# ---------------------------------------------------CLASS:  FitModes
#
class FitModes(object):
    """Indicators for single-frame or multiframe data files"""
    #
    DIRECT    = 0
    MULTIRING = 1
    #
    DEFAULT = MULTIRING
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:FitModes
# ---------------------------------------------------CLASS:  ImageModes
#
class ImageModes(object):
    """Indicators for single-frame or multiframe data files"""
    #
    SINGLE_FRAME = 0
    MULTI_FRAME  = 1
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  ImageModes
# ---------------------------------------------------CLASS:  Experiment
#
class Experiment(object):
    """Experiment class:  wrapper for XRD functionality"""
    def __init__(self, cfgFile, matFile):
        """Constructor for Experiment

        INPUTS
        
        cfgFile -- name of the config file to use for initialization; 
                   an empty string indicates that default values for 
                   options are used
        matFile -- name of the materials data file; a real file name
                   is required here
"""
        #
        #  Reader inputs and info
        #
        self.__active_rdr   = ReaderInput()
        self.__savedReaders = [self.__active_rdr]
        #
        self.__active_img = None
        self.__curFrame = 0
        self.__numFrame = 0
        #
        #  Load material lists
        #
        self.matList = loadMaterialList(matFile)
        self.activeMaterial = 0
        #
        #  Detector and calibration information.  
        #
        self._detInfo  = DetectorInfo()
        self._calInput = CalibrationInput(self.matList[0])
        #
        #  Add hydra interface
        #
        self._hydra = Hydra()

        return

    def __str__(self):
        """Description"""
        s = 'Experiment instance'
        #
        return s
    #
    # ============================== API
    #
    #                     ========== Properties
    #
    # property:  detector

    @property
    def detector(self):
        """(read only) detector"""
        return self._detInfo.detector

    @property
    def refineFlags(self):
        """(read only) refinement flags for calibration"""
        return self._detInfo.refineFlags

    def saveDetector(self, fname):
        """Save the detector information to a file
        
        INPUTS
        fname -- the name of the file to save in
"""
        f = open(fname, 'w')
        self._detInfo.mrbImages = [] # remove images before saving
        cPickle.dump(self._detInfo, f)
        f.close()
        
        return

    def loadDetector(self, fname):
        """Load the detector information from a file

        INPUTS
        fname -- the name of the file to load from
"""
        # should check the loaded file here
        f = open(fname, 'r')
        self._detInfo = cPickle.load(f)
        f.close()
        
        return

    #
    # ==================== Calibration Input
    #
    # property:  calInput

    @property
    def calInput(self):
        """(get only) Calibration input instance"""
        return self._calInput
    
    # ==================== Hydra
    #
    # property:  hydra

    @property
    def hydra(self):
        """(read only) hydra image class"""
        return self._hydra

    # ==================== Materials
    #
    # property:  activeMaterial

    def _get_activeMaterial(self):
        """Get method for activeMaterial"""
        return self._active_mat

    def _set_activeMaterial(self, v):
        """Set method for activeMaterial"""
        if isinstance(v, int):
            self._active_mat = self.matList[v]
        else:
            # v is a string
            self._active_mat = self.matDict[v]
            pass

        return

    _amdoc = r"""Active Material

    Can be set by number (index in material list) or by name. 

    On output, it is always a material instance.
"""
    activeMaterial = property(_get_activeMaterial, _set_activeMaterial, None,
                              _amdoc)

    # property:  matList

    def _get_matList(self):
        """Get method for matList"""
        return self._matList

    def _set_matList(self, v):
        """Set method for matList"""
        self._matList = v
        self.activeMaterial = 0
        # On initialization, this is called before calInput exists
        print 'calibrant name: ', self.matList[0]
        try:
            self.calInput.calMat = self.matList[0]
        except:
            pass

        return

    matList = property(_get_matList, _set_matList, None,
                                "List of materials")

    @property
    def matNames(self):
        """(read only) List of material names"""
        return [m.name for m in self.matList]

    @property
    def matDict(self):
        """(read only) Dictionary mapping material names to material"""
        return dict(zip(self.matNames, self.matList))

    def newMaterial(self):
        """Create a new material and add it to the list"""
        self._active_mat = Material()

        # find name not already in list
        n  = self._active_mat.name
        self._active_mat.name = newName(n, self.matNames)
        # 
        self._matList.append(self.activeMaterial)

        return

    def loadMaterialList(self, fname):
        """Load the pickled material list from a file

        INPUTS
        fname -- the name of the file to load from
"""
        # should check the loaded file here
        f = open(fname, 'r')
        self.matList = cPickle.load(f)
        f.close()
        
        return
    #
    # ==================== Readers
    #
    # property:  activeReader

    def _get_activeReader(self):
        """Get method for activeReader

        Reader is set by using index in reader list or by name.
"""
        return self.__active_rdr

    def _set_activeReader(self, v):
        """Set method for activeReader"""
        if isinstance(v, int):
            self.__active_rdr = self.__savedReaders[v]
        else:
            # v is a string
            for r in self.__savedReaders:
                if r.name == v:
                    self.__active_rdr = r
                    return
                pass
            pass

        return

    _ardoc = r"""Active Material

    Can be set by number (index in saved readers list) or by name. 

    On output, it is always a ReaderInput instance.
"""

    activeReader = property(_get_activeReader, _set_activeReader,
                            _ardoc)

    def saveReaderList(self, fname):
        """Save the reader list to a file
        
        INPUTS
        fname -- the name of the file to save in
"""
        f = open(fname, 'w')
        cPickle.dump(self.__savedReaders, f)
        f.close()
        
        return

    def loadReaderList(self, fname):
        """Load the reader list from a file

        INPUTS
        fname -- the name of the file to load from
"""
        # should check the loaded file here
        f = open(fname, 'r')
        self.__savedReaders = cPickle.load(f)
        self.activeReader = 0
        f.close()
        
        return

    def newReader(self):
        """Add new reader to the list and make it active

        Changes name if necessary.
"""
        self.__active_rdr   = ReaderInput()
        # find name not already in list
        n  = self.__active_rdr.name
        nl = [r.name for r in self.__savedReaders]
        self.__active_rdr.name = newName(n, nl)
        # 
        self.__savedReaders.append(self.__active_rdr)
        
        return

    def getSavedReader(self, which):
        """Get a specified reader"""
        if isinstance(which, int):
            return self.__savedReaders[v]
        else:
            # which is a string
            for r in self.__savedReaders:
                if r.name == which:
                    return r
                pass
            pass
        return r

    @property
    def savedReaders(self):
        """Return list of saved readers"""
        return self.__savedReaders

    @property
    def readerNames(self):
        """Return list of saved readers"""
        return [r.name for r in self.__savedReaders]
    #
    # ==================== Image Info
    #
    @property
    def curFrameNumber(self):
        """Current frame number"""
        return self.__curFrame

    @property
    def numFramesTotal(self):
        """Number of frames available for reading"""
        return self.__numFrame

    @property
    def activeImage(self):
        """Active image"""
        return self.__active_img
    #
    # ==================== Calibration
    #
    # property:  calInput

    @property
    def calInput(self):
        """(read only) Calibration input data"""
        return self._calInput
    #
    #                     ========== Public Methods
    #
    def readerListAddCurrent(self):
        """Add current list to list of saved readers"""
        
        return

    def readImage(self, frameNum=1):
        """Read and return an image

        DESCRIPTION
        
        This reads an image according to the active reader
        specification, saving it in the activeImage attribute.
"""
        #
        # Now read the current frame
        #
        aggMode = self.activeReader.aggModeOp
        nrFrame = self.activeReader.getNumberOfFrames() # number of reader frames
        if aggMode:
            rdFrames = nrFrame
            self.__numFrame = 1
        else:
            rdFrames = 1
            self.__numFrame = nrFrame
            pass
        #
        #  If current frame is 0, no reader has yet been instantiated.
        #
        haveReader = (self.__curFrame > 0)
        #
        #  Check frameNum
        #
        if (frameNum > nrFrame) or (frameNum < 1):
            msg = 'frame number out of range: requesting frame %d (max = %d)' % (frameNum, nrFrame)
            raise ValueError(msg)

        #if (frameNum == self.__curFrame): return
        # NOTE:  instantiate new reader even when requested frame is current
        # frame because reader properties may have changed
        
        if haveReader and (frameNum > self.__curFrame):
            nskip = frameNum - self.__curFrame - 1
            self.__active_img = self.__active_reader.read(nframes= rdFrames, 
                                                          nskip  = nskip,
                                                          sumImg = aggMode)
        else:
            # instantiate new reader
            self.__active_reader = self.activeReader.makeReader()
            nskip = frameNum - 1
            self.__active_img = self.__active_reader.read(nframes= rdFrames, 
                                                          nskip  = nskip,
                                                          sumImg = aggMode)
            
            pass
        
        self.__curFrame = frameNum
        print 'frame:  (exp) %d, (rdr) %s' % (self.__curFrame, 
                                              str(self.__active_reader.iFrame))

        return 

    def calibrate(self, log=None):
        """Calibrate the detector

        Currently, uses polar rebin only.
"""
        try:
            self._detInfo.calibrate(self.calInput,
                                    self.activeReader, 
                                    self.activeMaterial, log=log)
        except Exception as e:
            if log:
                log.write(str(e) + '\n')
                raise
            else:
                raise
            pass
        
        if log:
            log.write('done')
        
        return
        #
        # ==================== Polar Rebinning (Caking)
        #
        def polarRebin(self, opts):
            """Rebin the image according to certain parameters

            opts -- an instance of PolarRebinOpts
"""
            
            img_info = det.polarRebin(self.activeImage, opts.kwArgs)

            return img_info
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  Experiment
#
# ================================================== Utility Functions
#
def newName(name, nlist):
    """return a name not in the list, but based on name input"""
    if name not in nlist: return name

    i=1
    while i:
        name_i = '%s %d' % (name, i)
        if name_i not in nlist:
            break

        i += 1
        pass

    return name_i

def saveExp(e, f):
    """save experiment to file"""
    print 'saving file:  %s' % f
    fobj = open(f, 'w')
    cPickle.dump(e, fobj)
    fobj.close()
    print 'save succeeded'

    return

def loadExp(inpFile, matFile=DFLT_MATFILE):
    """Load an experiment from a config file or from a saved exp file

    inpFile -- the name of either the config file or the saved exp file;
               empty string means start new experiment
    matFile -- name of the materials file 
"""
    #
    if not matFile:
        print >> sys.stderr, 'no material file found'
        sys.exit(1)
        pass
        
    root, ext = os.path.splitext(inpFile)
    #
    if ext == '.cfg' or not inpFile:
        #  Instantiate from config file
        exp = Experiment(inpFile, matFile)
    elif ext == '.exp':
        #  Load existing experiment
        try:
            print 'loading saved file:  %s' % inpFile
            f = open(inpFile, 'r')
            exp = cPickle.load(f)
            f.close()
            print '... load succeeded'
        except:
            print '... load failed ... please check your data file'
            raise
            sys.exit()
            pass
    else:
        #  Not recognized
        print 'file is neither .cfg or .exp', inpFile
        sys.exit(1)
        pass

    return exp

# ---------------------------------------------------CLASS:  geReaderInput
#
class ReaderInput(object):
    """ReaderInput

This class is for holding input required to 
instantiate a reader object.  Currently, only
GE reader is supported.    
"""
    # Class data
    DFLT_NAME = 'unnamed reader'
    #
    DARK_MODES = (DARK_MODE_NONE,
                  DARK_MODE_FILE,
                  DARK_MODE_EMPTY,
                  DARK_MODE_FRAME) = range(4)
    #
    AGG_MODES = (AGG_FUN_NONE, AGG_FUN_SUM, AGG_FUN_MAX, AGG_FUN_MIN) = range(4)
    AGG_DICT = {
        AGG_FUN_NONE: False, 
        AGG_FUN_SUM : True,  # or (1/nframes)*numpy.add
        AGG_FUN_MAX : numpy.maximum, 
        AGG_FUN_MIN : numpy.minimum
        }
    #
    FLIP_MODES = (FLIP_NONE, FLIP_VERT, FLIP_HORIZ, FLIP_180, FLIP_M90, FLIP_P90) = range(6)
    FLIP_STRS  = ('',        'v',       'h',        'hv',     'cw90',   'ccw90')
    FLIP_DICT  = dict(zip(FLIP_MODES, FLIP_STRS))
    #
    def __init__(self, name='reader', desc='no description'):
        """Constructor for ReaderInput

        INPUT
        name -- [optional] (str) name
        desc -- [optional] (str) description

        NOTES
        * currently only GE reader is supported
"""
        self.name = name
        self.desc = desc
        # Image Mode
        # . imageDict contains tuple of (#empty, o-min, o-max, o-del)
        self.imageMode = ImageModes.SINGLE_FRAME
        self.imageDir = ''
        self.imageNames = []
        self.imageNameD = dict()
        # Dark file
        self.darkMode = ReaderInput.DARK_MODE_NONE
	self.darkDir  = ''
        self.darkName = ''
        # File aggregation
        self.aggFun = ReaderInput.AGG_FUN_NONE
        # Flip options
        self.flipMode = ReaderInput.FLIP_NONE
        #
        return

    def _check(self):
        """Check that input is ok for making a reader instance
"""
        # * Empty frames = 0 for single frame mode
        return
    #
    # ============================== API
    #
    #                     ========== Properties
    #
    # property:  hasImages

    @property
    def hasImages(self):
        """(get only) true if list of images has been set """
        return len(self.imageNames) > 0

    # property:  imageNames

    def _get_imageNames(self):
        """Get method for imageNames"""
        return self._imageNames

    def _set_imageNames(self, v):
        """Set method for imageNames"""
        self._imageNames = v
        # Set up dictionary for each file
        #     value is:  (numEmpty, o-min, o-max, o-del)
        for fn in v:
            self.imageNameD[fn] = (0, None, None, None)
            pass

        return

    imageNames = property(_get_imageNames, _set_imageNames,
                                "List of file image file names")
    # property:  "darkFile"

    def _get_darkFile(self):
        """Get method for darkFile"""
        if self.darkMode == ReaderInput.DARK_MODE_NONE:
            s = '<no dark subtraction>'
        elif self.darkMode == ReaderInput.DARK_MODE_FILE:
            s = os.path.join(self.darkDir, self.darkName)
        else:
            s = '<using empty frames>'
        return s

    darkFile = property(_get_darkFile, None, None,
                                "Full pathname of dark file" )

    @property
    def aggMode(self):
        """Mode identifier for frame aggregation"""
        return self.aggFun

    # property:  aggModeOp

    @property
    def aggModeOp(self):
        """(read only) option to pass to GE reader instances for aggregation mode"""
        return ReaderInput.AGG_DICT[self.aggFun]
    #
    #                     ========== Public Methods
    #  Omega Info
    #
    def setOmegaInfo(self, imgName, omin, omax, odel):
        """Set omega info for the specified image"""
        if imgName not in self.imageNameD:
            raise KeyError('image not in image list')

        info = self.imageNameD[imgName]
        self.imageNameD[imgName] = (info[0], omin, omax, odel)

        return
    #
    def getNumberOfFrames(self):
        """Return number of frames available in data files"""
        try:
            r = self.makeReader()
            n = r.getNFrames()
        except Exception:
            raise
            n = 0
            pass

        return n

    def makeReader(self):
        """Return a reader instance based on self
"""
        # check validity of input
        self._check()
        #
        # Set up image names in right format
        #
        fullPath = lambda fn: os.path.join(self.imageDir, fn)
        numEmpty = lambda fn: self.imageNameD[fn][0]
        imgInfo = [(fullPath(f), numEmpty(f)) for f in self.imageNames]
        #
        # Check for omega info
        #
        nfile = len(imgInfo)
        dinfo = [self.imageNameD[f] for f in self.imageNames]
        omin = dinfo[0][1]
        if omin is not None:
            odel = dinfo[nfile - 1][3]
            print "omega min and delta: ", omin, odel
            omargs = (float(omin), float(odel))
        else:
            omargs = ()
            pass
        print 'omargs:  ', omargs
        #
        # Dark file
        #
        subDark = not (self.darkMode == ReaderInput.DARK_MODE_NONE)
        if (self.darkMode == ReaderInput.DARK_MODE_FILE):
            drkFile = os.path.join(self.darkDir, self.darkName) 
        else:
            drkFile = None
            pass
        #
        # Flip options
        #
        doFlip  = not (self.flipMode == ReaderInput.FLIP_NONE)
        flipArg = ReaderInput.FLIP_DICT[self.flipMode]
        #
        # Make the reader
        #
        print 'reader:  \n', imgInfo, subDark, drkFile, doFlip, flipArg
        r = detector.ReadGE(imgInfo, *omargs,
                            subtractDark = subDark,
                            dark         = drkFile,
                            doFlip       = doFlip,
                            flipArg      = flipArg)

        return r
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  geReaderInput
# ---------------------------------------------------CLASS:  CalibrationInput
#
class CalibrationInput(object):
    """CalibrationInput"""
    def __init__(self, mat):
        """Constructor for CalibrationInput"""
        #
        self.numRho = 20 # for multiring binned image
        self.numEta = 36
        
        self.corrected = False
        #
	self.calMat = mat

        return
    #
    # ============================== API
    #
    # property:  fitType

    def _get_fitType(self):
        """Get method for fitType"""
        if not hasattr(self, '_fitType'):
            self._fitType = FitModes.DEFAULT
            pass
        
        return self._fitType

    def _set_fitType(self, v):
        """Set method for fitType"""
        self._fitType = v
        return

    fitType = property(_get_fitType, _set_fitType, None,
                                "fit type:  direct or caked")
    # property:  calMat

    def _get_calMat(self):
        """Get method for calMat"""
        return self._calMat

    def _set_calMat(self, v):
        """Set method for calMat"""
        self._calMat = v
        print 'setting calMat:  \n', str(self.calData)
        return

    calMat = property(_get_calMat, _set_calMat, None,
                                "Calibration material (calibrant)")

    # property:  calData

    @property
    def calData(self):
        """(get only) Lattice parameter data for calibrant

        This provides a deepcopy with wavelength, strain magnitude and
        two-theta width set.
"""
        return self.calMat.planeData

    # property:  cakeArgs

    @property
    def cakeArgs(self):
        """(get only) Keyword arguments for polar rebinning"""
        
        return {'verbose':True, 'numEta': self.numEta, 'corrected':self.corrected}


    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  CalibrationInput
# ---------------------------------------------------CLASS:  DetectorInfo
#
class DetectorInfo(object):
    """detectorInfo"""
    # refinement-specific things
    #  ---------------> (   xc,    yc,     D,    xt,    yt,    zt,   
    #                      dp1,   dp2,   dp3,   dp4,   dp5,   dp6)
    DFLT_REFINE_FLAGS = ( True,  True,  True,  True,  True, False, 
                          True,  True,  True, False, False, False)
    def __init__(self):
        """Constructor for detectorInfo"""
        #
	self.detector  = detector.newDetector('ge')
        self.mrbImages = []
        self.fitParams = []
        #
        # Refinement variables
        #
        self.refineFlags = list(DetectorInfo.DFLT_REFINE_FLAGS)
        #
        self.calInp = self.rdrInp = None
        #
        return
    #
    # ============================== API
    #
    def calibrate(self, calInp, rdrInp, mat, log=None):
        """Calibrate this detector using specified reader and options"""
        #
        #  Save inputs for the record
        #
        self.calInp = calInp
        self.rdrInp = rdrInp
        self.mat    = mat
        #
        self.detector.setupRefinement(tuple(self.refineFlags))
        #
        #  Set calibrant properties
        #
        calDat = mat.planeData
        #
        #  Generate images
        #
        self.fitParams = []
        self.mrbImages = []
        #
        reader = rdrInp.makeReader()
        nf = reader.getNFrames()
        for i in range(nf):
            msg = '*** fitting frame %d of %d\n' %(i+1, nf) + \
                  '*** detector x, y, D = (%g, %g, %g)\n' \
                  %(self.detector.xc, self.detector.yc, self.detector.workDist) + \
                  '*** material name:\t%s\n' %(mat.name) + \
                  '*** target num rho:\t%d\n' %(calInp.numRho) + \
                  'using the following HKLs:\n' \
                  + mat.planeData.hkls.__repr__()
            #
            #  Reset tthmax here in case detector geometry pushes
            #  a ring off the detector.
            #
            # NO!!! # calDat.tThMax = self.detector.getTThMax() 
            
            cFrame = reader()
            if log:
                log.write(msg)
            else:
                print msg
                pass
            
            if calInp.fitType == FitModes.MULTIRING:
                mrb = detector.MultiRingBinned(self.detector, calDat, cFrame, 
                                               targetNRho=calInp.numRho,
                                               polarRebinKWArgs=calInp.cakeArgs,
                                               log=log)
                tmp = mrb.doFit()
                self.mrbImages.append(mrb)
            else:
                print '... using direct fit mode'
                
                self.detector.fitRings(cFrame, calDat)
                tmp = self.detector.xFitRings
                pass
            
            self.fitParams.append(self.detector.getParams(allParams=True))

            print 'fit parameters(%d):\n' % i, self.fitParams[-1]

            pass

        # mean and std dev of geometric parameters
        meanParams = numpy.mean(numpy.array(self.fitParams), 0).tolist()
        stdvParams = numpy.std(numpy.array(self.fitParams), 0).tolist()


        # make detector object from mean of phi = 0-180 scans
        self.detector._Detector2DRC__updateFromPList(meanParams)
        
        return
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  DetectorInfo
# ---------------------------------------------------CLASS:  PolarRebinOpts
#
class PolarRebinOpts(object):
    """Options for polar rebinning"""
    #
    cakeMethods = [CAKE_IMG, CAKE_RNG, CAKE_SPH] = ['standard', 'multiring', 'spherical']
    
    def __init__(self):
        """Constructor for PolarRebinOpts

        This routine sets default values for caking options.

        The following attributes (with initial values) can be modified directly. 
            etaMin  =    0
            etaMax  =  360
            rhoMin  =  100
            rhoMax  = 1000
            numEta  =   36
            numRho  =  500
            correct = True
"""
        #
	self.type = cakeMethods[0]
        #
        #  Standard (whole image) rebinning
        #
        etaMin  =    0
        etaMax  =  360
        rhoMin  =  100
        rhoMax  = 1000
        numEta  =   36
        numRho  =  500
        correct = True
        #
        return
    #
    # ============================== API
    #
    # property:  kwArgs

    @property
    def kwArgs(self):
        """(get only) Return keyword args to pass to polarRebin"""
        kwa = dict()
        #
        if self.type == CAKE_IMG:
            kwa = {
                'etaRange' : [self.etaMin, self.etaMax],
                'numEta'   : self.numEta,
                'rhoRange' : [self.rhoMin, self.rhoMax],
                'numRho'   : self.numRho,
                'corrected': self.correct
                }
            pass
        
        return kwa

    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  PolarRebinOpts
