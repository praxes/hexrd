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
#
######################################################################
## TOP-LEVEL MODULES AND SOME GLOBALS
##
"""Module for wrapping the main functionality of the xrd package.

The Experiment class is the primary interface.  Other classes
are helpers.
"""
import sys, os, copy
import cPickle
import numpy

from scipy.linalg          import inv
from scipy.linalg.matfuncs import logm
from scipy                 import optimize

from hexrd import matrixutil
from hexrd import valunits

from hexrd.xrd import detector
from hexrd.xrd import grain      as G
from hexrd.xrd import indexer    
from hexrd.xrd import rotations  as ROT
from hexrd.xrd import spotfinder as SPT
from hexrd.xrd import xrdutil

from hexrd.xrd.hydra    import Hydra
from hexrd.xrd.material import Material, loadMaterialList

from math import pi
r2d = 180. / pi
d2r = pi / 180.

#
#  Defaults (will eventually make to a config file)
#
HERE = os.path.dirname(__file__)
toMatFile    = os.path.join(HERE, '..', 'data', 'materials.cfg')
DFLT_MATFILE = os.path.normpath(toMatFile) # check whether it exists
matfileOK    = os.access(DFLT_MATFILE, os.F_OK)
if not matfileOK:  # use relative path
    DFLT_MATFILE = os.path.join('data', 'materials.cfg')
    pass
matfileOK = os.access(DFLT_MATFILE, os.F_OK)
if not matfileOK:  # set to null
    DFLT_MATFILE = ''
    pass
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
    """Wrapper for xrd functionality"""
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
        #  Spots information.
        #
        self._spotOpts = SpotOptions()
        self._spots = []
        self._spots_ind = []
        self._spot_readers = []
        #
        #  Index Information
        #
        self._index_opts = IndexOptions()
        self._fitRMats = []
        #
        #  Image Lists
        #
        self._img_list = []
        self._img_names = []
        #
        #  Add hydra interface (not really functional)
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
    # ==================== Image List
    #
    # property:  active_img

    @property
    def active_img(self):
        """Current image"""
        return self.__active_img

    @active_img.setter
    def active_img(self, v):
        """Set active image from list by number or name"""
        if isinstance(v, int):
            self.__active_img = self.img_list[v]
        else:
            # v is a string
            istar = 0
            for i in range(len(self.img_list)):
                if v == self.img_names[i]:
                    istar = i
                    break
            
            self.__active_img = self.img_list[istar]
            pass
        
    # property:  img_list

    @property
    def img_list(self):
        """(get only) List of saved images (get only)"""
        if not hasattr(self, '_img_list'):
            self._img_list = []
            self._img_names = []
            
        return self._img_list

    def add_to_img_list(self, name):
        """Append the active image to the image list"""
        # NEED TO: check names for duplication
        self.img_names.append(name)
        self.img_list.append(self.active_img)
        
        return

    # property:  img_names

    @property
    def img_names(self):
        """(get-only) List of names for saved images"""
        if not hasattr(self, '_img_names'):
            self._img_names = []
            
        return self._img_names
#
    # ==================== Indexing
    #
    # property:  fitRMats

    @property
    def fitRMats(self):
        """(get-only) Rotation matrices from indexing"""
        if not hasattr(self, '_fitRMats'):
            self._fitRMats = []
        return self._fitRMats
    
    # property:  index_opts
    
    def refine_grains(self, 
                      minCompl,
                      nSubIter=3,
                      doFit=False,
                      etaTol=valunits.valWUnit('etaTol', 'angle', 1.0, 'degrees'), 
                      omeTol=valunits.valWUnit('etaTol', 'angle', 1.0, 'degrees'), 
                      fineDspTol=5.0e-3, 
                      fineEtaTol=valunits.valWUnit('etaTol', 'angle', 0.5, 'degrees'), 
                      fineOmeTol=valunits.valWUnit('etaTol', 'angle', 0.5, 'degrees')):
        """
        refine a grain list
        """
        # refine grains formally using a multi-pass refinement
        nGrains    = self.rMats.shape[0]
        grainList = []
        for iG in range(nGrains):
            indexer.progress_bar(float(iG) / nGrains)
            grain = G.Grain(self.spots_for_indexing,
                                rMat=self.rMats[iG, :, :], 
                                etaTol=etaTol,
                                omeTol=omeTol,
                                claimingSpots=False)
            if grain.completeness > minCompl:
                for i in range(nSubIter):
                    grain.fit()
                    s1, s2, s3 = grain.findMatches(etaTol=etaTol, omeTol=omeTol, strainMag=fineDspTol,
                                                   updateSelf=True, claimingSpots=False, doFit=doFit, 
                                                   testClaims=True)
                if grain.completeness > minCompl:                    
                    grainList.append(grain)
                    pass
                pass
            pass
        self.grainList = grainList
        self._fitRMats = numpy.array([self.grainList[i].rMat for i in range(len(grainList))])
        return
    
    def saveRMats(self, f):
        """save rMats to npy file"""
        numpy.save(f, self.rMats)
        return
    
    def dump_grainList(self, f):
        """dump grainList to cPickle"""
        if isinstance(f, file):
            fid = f
        elif isinstance(f, str) or isinstance(f, unicode):
            fid = open(f, 'w')
            pass
        cPickle.dump(self.grainList, fid)
        fid.close()
        return
    
    def export_grainList(self, f, 
                         dspTol=5.0e-3,
                         etaTol=valunits.valWUnit('etaTol', 'angle', 0.5, 'degrees'), 
                         omeTol=valunits.valWUnit('etaTol', 'angle', 0.5, 'degrees'),
                         doFit=False, 
                         sort=True):
        """
        export method for grainList
        """
        if isinstance(f, file):
            fid = f
        elif isinstance(f, str) or isinstance(f, unicode):
            fid = open(f, 'w')
            pass
        
        if sort:
            loop_idx = numpy.argsort([self.grainList[i].completeness 
                                      for i in range(len(self.grainList))])[::-1]
        else:
            loop_idx = range(len(self.grainList))
            pass
        
        for iG in loop_idx:
            #
            # this grain
            grain = self.grainList[iG]
            # 
            # useful locals
            q   = ROT.quatOfRotMat(grain.rMat)
            R   = grain.rMat
            V   = grain.vMat
            FnT = inv(numpy.dot(V, R)).T
            E   = logm(V)
            Es  = logm(grain.uMat)
            lp  = grain.latticeParameters
            p   = grain.detectorGeom.pVec
            #
            # the output
            print >> fid, '\n#####################\n#### grain %d\n' % (iG) + \
                  '\n#    orientation:\n#\n' + \
                  '#    q = [%1.6e, %1.6e, %1.6e, %1.6e]\n#\n' % (q[0], q[1], q[2], q[3]) + \
                  '#    R = [[%1.3e, %1.3e, %1.3e],\n' % (R[0, 0], R[0, 1], R[0, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e],\n' % (R[1, 0], R[1, 1], R[1, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e]]\n' % (R[2, 0], R[2, 1], R[2, 2]) + \
                  '#\n#    left stretch tensor:\n#\n' + \
                  '#    V = [[%1.3e, %1.3e, %1.3e],\n' % (V[0, 0], V[0, 1], V[0, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e],\n' % (V[1, 0], V[1, 1], V[1, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e]]\n' % (V[2, 0], V[2, 1], V[2, 2]) + \
                  '#\n#    logarithmic strain tensor (log(V) --> sample frame):\n#\n' + \
                  '#    E = [[%1.3e, %1.3e, %1.3e],\n' % (E[0, 0], E[0, 1], E[0, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e],\n' % (E[1, 0], E[1, 1], E[1, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e]]\n' % (E[2, 0], E[2, 1], E[2, 2]) + \
                  '#\n#    logarithmic strain tensor (log(U) --> crystal frame):\n#\n' + \
                  '#    E = [[%1.3e, %1.3e, %1.3e],\n' % (Es[0, 0], Es[0, 1], Es[0, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e],\n' % (Es[1, 0], Es[1, 1], Es[1, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e]]\n' % (Es[2, 0], Es[2, 1], Es[2, 2]) + \
                  '#\n#    F^-T ( hkl --> (Xs, Ys, Zs), reciprocal lattice to sample frame ):\n#\n' + \
                  '#        [[%1.3e, %1.3e, %1.3e],\n' % (FnT[0, 0], FnT[0, 1], FnT[0, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e],\n' % (FnT[1, 0], FnT[1, 1], FnT[1, 2]) + \
                  '#         [%1.3e, %1.3e, %1.3e]]\n' % (FnT[2, 0], FnT[2, 1], FnT[2, 2]) + \
                  '#\n#    lattice parameters:\n#\n' + \
                  '#    %g, %g, %g, %g, %g, %g\n' % tuple(numpy.hstack([lp[:3], r2d*numpy.r_[lp[3:]]])) + \
                  '#\n#    COM coordinates (Xs, Ys, Zs):\n' +\
                  '#    p = (%g, %g, %g)\n' % (p[0], p[1], p[2]) + \
                  '#\n#    reflection table:'
            s1, s2, s3 = grain.findMatches(etaTol=etaTol, omeTol=omeTol, strainMag=dspTol,
                                           updateSelf=True, claimingSpots=True, doFit=doFit, filename=fid)
            print >> fid, '#\n#    final completeness for grain %d: %g%%\n' % (iG, grain.completeness*100) + \
                  '#####################\n'
            pass
        
        fid.close()
        
        return
    
    def simulateGrain(self, 
                      rMat=numpy.eye(3),
                      vMat=numpy.r_[1., 1., 1., 0., 0., 0.],
                      planeData=None,
                      detector=None,
                      omegaRanges=[(-pi, pi),], 
                      output=None):
        """
        Simulate a grain with choice of active material
        """ 
        if planeData is None:
            planeData = self.activeMaterial.planeData
        if detector is None:
            detector = self.detector
        dummySpots = SPT.Spots(planeData, None, detector, omegaRanges)
        sg = G.Grain(dummySpots, rMat=rMat, vMat=vMat)
        if output is not None:
            if isinstance(output, file):
                fid = output
            elif isinstance(output, str):
                fid = open(output, 'w')
            else:
                raise RuntimeError, "output must be a file object or string"
            sg.findMatches(filename=output)
        return sg
    
    @property
    def index_opts(self):
        """(get-only) Options for indexing"""
        if not hasattr(self, '_index_opts'):
            self._index_opts = IndexOptions()
        return self._index_opts

    def _run_grainspotter(self):
        return
    
    def _run_fiber_search(self):

        iopts = self.index_opts
        fsearch = indexer.fiberSearch
        myspots = self.spots_for_indexing
        print 'my spots: ', myspots
        retval = fsearch(myspots, iopts.fsHKLs,
                         nsteps=iopts.nsteps,
                         minCompleteness=iopts.minCompleteness,
                         minPctClaimed=iopts.minPctClaimed,
                         preserveClaims=iopts.preserveClaims,
                         friedelOnly=iopts.friedelOnly,
                         dspTol=iopts.dspTol,
                         etaTol=iopts.etaTol * d2r,
                         omeTol=iopts.omeTol * d2r,
                         doRefinement=iopts.doRefinement,
                         doMultiProc=iopts.doMultiProc,
                         nCPUs=iopts.nCPUs,
                         quitAfter=iopts.quitAfter,
                         outputGrainList=True)
        iopts._fitRMats = retval[0] # HUH?!
        self.rMats = retval[0]
        self.grainList = retval[1]
        return
            
    def run_indexer(self):
        """Run indexer"""
        iopts = self.index_opts

        if iopts.index_method == iopts.IND_FIBER:
            self._run_fiber_search()
        else:
            self._run_grainspotter()
        
        return
    #
    # ==================== Spots
    #
    def clear_spots(self):
        """Reset the list of spots"""
        self._spots = []
        self._spots_ind = []
        
        return
    # property:  spot_readers

    def saveRawSpots(self, fname):
        """Save the detector information to a file
        
        INPUTS
        fname -- the name of the file to save in
        """
        f = open(fname, 'w')
        cPickle.dump(self._spots, f)
        f.close()
        
        return

    def loadRawSpots(self, fname):
        """Load the detector information from a file

        INPUTS
        fname -- the name of the file to load from

        """
        # should check the loaded file here
        f = open(fname, 'r')
        self._spots = cPickle.load(f)
        f.close()
        
        return

    @property
    def spot_readers(self):
        """(get-only) list of readers used to generate spots"""
        if not hasattr(self, '_spot_readers'):
            self._spot_readers = []
        return self._spot_readers
    
    # property:  spots_for_indexing

    @property
    def spots_for_indexing(self):
        """(get-only) spots associated with rings"""
        if not hasattr(self, '_spots_ind'):
            self._spots_ind = []
        return self._spots_ind
    
    # property:  raw_spots

    @property
    def raw_spots(self):
        """(get-only) spots from image before culling and association with rings"""
        if not hasattr(self, '_spots'):
            self._spots = []
        return self._spots
    
    # property:  spotOpts

    @property
    def spotOpts(self):
        """(get-only) spot finding options"""
        # in case of loading an old pickle without spot options
        if not hasattr(self, '_spotOpts'):
            self._spotOpts = SpotOptions()
        return self._spotOpts

    def find_raw_spots(self):
        """find spots using current reader and options"""
        findSOS = SPT.Spots.findSpotsOmegaStack
        opts = self.spotOpts
        #
        newspots = findSOS(self.activeReader.makeReader(), 
                           opts.nframes,
                           opts.thresh, 
                           opts.minPix, 
                           discardAtBounds=opts.discardAtBounds,
                           keepWithinBBox=opts.keepWithinBBox,
                           overlapPixelDistance=opts.overlap,
                           nframesLump=opts.nflump,
                           padOmega=opts.padOmega,
                           padSpot=opts.padSpot,
                           )
        self._spots += newspots
        self._spot_readers.append(self.activeReader.name)
        
        return

    def get_spots_ind(self):
        """Select spots for indexing"""
        # cull spots that have integrated intensity <= 0
        spots_get_II = SPT.Spot.getIntegratedIntensity
        integIntensity = numpy.array(map(spots_get_II, self.raw_spots))
        rspots = numpy.array(self.raw_spots)
        culledSpots = rspots[integIntensity > 0.]

        pd = self.activeMaterial.planeData
        readerInpList = [self.getSavedReader(rn) for rn in self._spot_readers]
        readerList = [ri.makeReader() for ri in readerInpList]
        ominfo = [reader.getOmegaMinMax() for reader in readerList]
        self._spots_ind = SPT.Spots(pd, culledSpots, self.detector, ominfo)        
        return
    #
    # ==================== Detector
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

    def newDetector(self, gp, dp):
        """Create a new detector with given geometry and distortion parameters

        *gp* - initial geometric parameters
        *dp* - initial distortion parameters

        """
        self._detInfo  = DetectorInfo(gParms=gp, dParms=dp)

        return
    
    def saveDetector(self, fname):
        """Save the detector information to a file

        INPUTS
        fname -- the name of the file to save in
        """
        f = open(fname, 'w')
        # self._detInfo.mrbImages = [] # remove images before saving
        # cPickle.dump(self._detInfo, f)
        det_class = self.detector.__class__
        det_plist = self.detector._Detector2DRC__makePList()
        det_rflag = self.detector.refineFlags
        print >> f, "# DETECTOR PARAMETERS"
        print >> f, "# \n# %s\n#" % (det_class)
        for i in range(len(det_plist)):
            print >> f, "%1.8e\t%d" % (det_plist[i], det_rflag[i])
        f.close()

        return

    def loadDetector(self, fname):
        """Load the detector information from a file

        INPUTS
        fname -- the name of the file to load from

        """
        # should check the loaded file here
        f = open(fname, 'r')
        # 
        lines = f.readlines()
        # self._detInfo = cPickle.load(f)
        det_class_str = None
        for i in range(len(lines)):
            if 'class' in lines[i]:
                det_class_str = lines[i]
        f.seek(0)
        if det_class_str is None:
            raise RuntimeError, "detector class label not recongined in file!"
        else:
            plist_rflags = numpy.loadtxt(f)
            plist = plist_rflags[:, 0]
            rflag = numpy.array(plist_rflags[:, 1], dtype=bool)
            
            exec_str = "DC = detector." + det_class_str.split('.')[-1].split("'")[0]
            exec(exec_str)
            
            gp = plist[:6].tolist()
            if len(plist[6:]) == 0:
                dp = None
            else:
                dp = plist[6:].tolist()
            self._detInfo  = DetectorInfo(gParms=gp, dParms=dp)
            self.detector.setupRefinement(rflag)
            self._detInfo.refineFlags = rflag
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

    def clear_reader(self):
        """Close current reader"""
        self.__active_reader = None
        self.__curFrame = 0
        return

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
    def activeImage(self): # to be removed (use active_img instead)
        """Active image"""
        return self.active_img
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
            msg = 'frame number out of range: requesting frame %d (max = %d)' \
                  % (frameNum, nrFrame)
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
                  DARK_MODE_ARRAY,
                  DARK_MODE_EMPTY,
                  DARK_MODE_FRAME) = range(5)
    #
    AGG_MODES = (AGG_FUN_NONE, AGG_FUN_SUM, AGG_FUN_MAX, AGG_FUN_MIN) = range(4)
    AGG_DICT = {
        AGG_FUN_NONE: False,
        AGG_FUN_SUM : True,  # or (1/nframes)*numpy.add
        AGG_FUN_MAX : numpy.maximum,
        AGG_FUN_MIN : numpy.minimum
        }
    #
    FLIP_MODES = (FLIP_NONE, FLIP_VERT, FLIP_HORIZ, FLIP_180, FLIP_M90, FLIP_P90) \
                 = range(6)
    FLIP_STRS  = ('',        'v',       'h',        'hv',     'cw90',   'ccw90')
    FLIP_DICT  = dict(zip(FLIP_MODES, FLIP_STRS))
    #
    RC = detector.ReadGE        # HARD CODED DETECTOR CHOICE!!!
    
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
        elif self.darkMode == ReaderInput.DARK_MODE_FILE or self.darkMode == ReaderInput.DARK_MODE_ARRAY:
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
        
        ref_reader = self.RC(imgInfo)
        #
        # Check for omega info
        #
        nfile = len(imgInfo)
        dinfo = [self.imageNameD[f] for f in self.imageNames]
        omin = dinfo[0][1]
        if omin is not None:
            odel = dinfo[nfile - 1][3]
            print "omega min and delta: ", omin, odel
            omargs = (valunits.valWUnit('omin', 'angle', float(omin), 'degrees'), 
                      valunits.valWUnit('odel', 'angle', float(odel), 'degrees'))
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
        elif (self.darkMode == ReaderInput.DARK_MODE_ARRAY):
            drkFileName = os.path.join(self.darkDir, self.darkName) 
            drkFile     = ref_reader.frame(
                buffer=numpy.fromfile(drkFileName, 
                                      dtype=ref_reader.dtypeRead
                                      )
                )
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
        r = self.RC(imgInfo, *omargs,
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

        return {'verbose':True, 
                'numEta': self.numEta, 
                'etaRange':numpy.array([0., 2. * numpy.pi]),
                'corrected':self.corrected}
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  CalibrationInput
# ---------------------------------------------------CLASS:  DetectorInfo
#
class DetectorInfo(object):
    """Class for detector and associated data"""
    # refinement-specific things
    #  ---------------> (   xc,    yc,     D,    xt,    yt,    zt,
    #                      dp1,   dp2,   dp3,   dp4,   dp5,   dp6)
    DFLT_REFINE_FLAGS = ( True,  True,  True,  True,  True, False,
                          True,  True,  True, False, False, False)
    def __init__(self, gParms=[], dParms=[]):
        """Constructor for detectorInfo"""
        #
	if gParms:              # maintain this keyword for compatability with Don's usage
            self.detector  = detector.newDetector('ge', gParms=gParms, dParms=dParms) 
        else:
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
        """
        Constructor for PolarRebinOpts

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
# ---------------------------------------------------CLASS:  SpotOptions
#
class SpotOptions(object):
    """Manage options available for spot finding and analysis

    Mainly, this manages the keyword options to the findSpotsOmegaStack()
    static method in the Spots class.
    """
    #
    # call to findSpotsOmegaStack for reference:
    #
    ### def findSpotsOmegaStack(reader, 
    ###                         nFrames,
    ###                         threshold, minPx, 
    ###                         discardAtBounds=True,
    ###                         keepWithinBBox=True,
    ###                         overlapPixelDistance=None,  # float, if specified
    ###                         nframesLump=1,              # probably get rid of this eventually
    ###                         padOmega=True,
    ###                         padSpot=True,
    ###                         debug=False, pw=None):
    
    def __init__(self):
        """SpotOptions Constructor"""
        #
        self.nframes = 0   # means use all
        self.thresh = 1000 # need reasonable initial value
	self.minPix = 4    # reasonable initial value
        self.discardAtBounds = True
        self.keepWithinBBox = True
        self.overlap = None
        self.nflump = 1
        self.padOmega = True
        self.padSpot = True
        #
        # Keep debug=False, pw=None
        #
        return
    #
    # ============================== API
    #

    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  SpotOptions
# ---------------------------------------------------CLASS:  indexOptions
#
class IndexOptions(object):
    """indexOptions"""
    #
    # Class data
    #
    INDEX_CHOICES = ['Fiber Search', 'GrainSpotter']
    INDEX_CHOICE_IDS = [IND_FIBER, IND_GSPOT] = range(2)
    
    def __init__(self):
        """Constructor for indexOptions"""
        self.index_method = self.IND_FIBER
        
        self.fsHKLs=[]
        self.preserveClaims=False
        self.friedelOnly=False
        self.doRefinement=True
        self.doMultiProc=True
        self.etaTol=0.25
        self.omeTol=0.50
        self.minCompleteness=0.67
        self.minPctClaimed=0.70
        self.nsteps=360
        self.nCPUs=None
        self.dspTol=None
        self.quitAfter=None,

        return
    #
    pass
	
# -----------------------------------------------END CLASS:  indexOptions
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
    fobj = open(f, 'w')
    e.clear_reader() # close open files inside exp
    cPickle.dump(e, fobj)
    fobj.close()

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

## test utilities
def refineDetector(grainList, scl=None, gtol=1.0e-6):
    """
    """
    if scl is None:
        scl = numpy.r_[0.005, 200., 1000.]

    # need to grab initial guess for xc, zTilt
    # use first grain by default (they all have the same parameters)
    xc      = grainList[0].detectorGeom.xc
    zTilt   = grainList[0].detectorGeom.zTilt
    chiTilt = grainList[0].detectorGeom.chiTilt
    if chiTilt is None:
        chiTilt = 0.
    x0 = scl * numpy.r_[xc, zTilt, chiTilt]

    # call to optimization routine
    xopt = optimize.fmin_bfgs(objFunc, x0, args=(grainList, scl), gtol=gtol)
    
    # recall objective to set detector geometries properly with solution
    objFunc(xopt, grainList, scl)
    
    return xopt / scl

def objFunc(x, grainList, scl):
    """
    """
    x = x / scl                         # remove scaling
    xc      = x[0]                        # beam x-center
    zTilt   = x[1]                        # zTilt --> inclination of oscill. axis on detector
    chiTilt = x[2]                        # zTilt --> inclination of oscill. axis on detector
    
    for i in range(len(grainList)):
        grainList[i].detectorGeom.xc      = xc
        grainList[i].detectorGeom.zTilt   = zTilt
        grainList[i].detectorGeom.chiTilt = chiTilt
    # need a fresh detector object to hand to spots
    # use first grain by default (any will do)
    tmpDG      = grainList[0].detectorGeom.makeNew()
    tmpDG.pVec = None
    
    # reset the detector used by all spots
    # each grain currently carried th
    # ...PRIME CANDIDATE FOR OPTIMIZATION/CLEANUP/GENERALIZATION...
    # ...perhaps loop over only the spots used by the grains in grainList?...
    # grainList[0].spots.resetDetectorGeom(tmpDG)
    
    # strip out quantities to hand off to the fit objective fuction to get residual contribution
    resd = []
    for i in range(len(grainList)):
        spotIDs = grainList[i].grainSpots['iRefl']
        spotIDs = spotIDs[spotIDs >= 0]
        for j in range(len(spotIDs)):
            spot = grainList[i].spots._Spots__spots[spotIDs[j]]
            spot.setDetectorGeom(tmpDG, clobber=True)
            angCOM = spot.angCOM()
            grainList[i].spots._Spots__spotAngCoords[spotIDs[j]] = angCOM
            grainList[i].spots._Spots__spotXYOCoords[spotIDs[j]] = grainList[i].spots.detectorGeom.angToXYO( *angCOM )
            pass
        # refit grains to new detector -- mainly for fixing pVecs
        grainList[i].updateGVecs()
        grainList[i].fit(display=False)
        
        angAxs = ROT.angleAxisOfRotMat(grainList[i].rMat)
        biotT  = matrixutil.symmToVecMV(grainList[i].vMat - numpy.eye(3))
        pVec   = grainList[i].detectorGeom.pVec
        
        x = numpy.vstack([angAxs[0]*angAxs[1], biotT.reshape(6, 1), pVec.reshape(3, 1)])
        resd.append(grainList[i]._fitF_objFunc(x))
        pass
    resd = numpy.hstack(resd).flatten()
    return sum(resd**2)
