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
# -*-python-*-
import re
import copy
from math import pi
import warnings

import numpy as num
from scipy import constants as C

from hexrd.matrixutil import *
from hexrd import valunits
from hexrd.valunits import toFloat

# module vars
r2d = 180./pi
d2r = pi/180.

dUnit = 'angstrom'

outputDegrees = False
outputDegrees_bak = outputDegrees

def hklToStr(x):
    return re.sub('\[|\]|\(|\)','',str(x))

def tempSetOutputDegrees(val):
    global outputDegrees, outputDegrees_bak
    outputDegrees_bak = outputDegrees
    outputDegrees = val
    return
def revertOutputDegrees():
    global outputDegrees, outputDegrees_bak
    outputDegrees = outputDegrees_bak
    return

def cosineXform(a, b, c):
    """
    Spherical trig transform to take alpha, beta, gamma to expressions
    for cos(alpha*).  See ref below.

    [1] R. J. Neustadt, F. W. Cagle, Jr., and J. Waser, ``Vector algebra and
        the relations between direct and reciprocal lattice quantities''. Acta
        Cryst. (1968), A24, 247--248

    """
    cosar = ( num.cos(b)*num.cos(c) - num.cos(a) ) / ( num.sin(b)*num.sin(c) )
    sinar = num.sqrt(1 - cosar**2)
    return cosar, sinar

def processWavelength(arg):
    """
    Convert an energy value to a wavelength.  If argument has units of length or energy,
    will convert to globally specified unit type for wavelength (dUnit).  If argument
    is a scalar, assumed input units are keV.
    """
    if hasattr(arg, 'getVal'):
        if arg.isLength():
            retval = arg.getVal(dUnit)
        elif arg.isEnergy():
            try:
                speed  = C.c
                planck = C.h
            except:
                raise NotImplementedError, 'scipy does not have constants'
                # speed  = ...
                # planck = ...
            e = arg.getVal('J')
            retval = valunits.valWUnit('wavelength', 'length', planck*speed/e, 'm').getVal(dUnit)
        else:
            raise RuntimeError, 'do not know what to do with '+str(arg)
    else:
        keV2J = 1.e3*C.e
        e = keV2J * arg
        retval = valunits.valWUnit('wavelength', 'length', C.h*C.c/e, 'm').getVal(dUnit)

    return retval

def latticeParameters(lvec):
    """
    Generates direct and reciprocal lattice vector components in a
    crystal-relative RHON basis, X. The convention for fixing X to the
    lattice is such that a || x1 and c* || x3, where a and c* are
    direct and reciprocal lattice vectors, respectively.
    """
    lnorm = num.sqrt(num.sum(lvec**2, 0))

    a = lnorm[0]
    b = lnorm[1]
    c = lnorm[2]

    ahat = lvec[:, 0]/a
    bhat = lvec[:, 1]/b
    chat = lvec[:, 2]/c

    gama = num.arccos(num.dot(ahat, bhat))
    beta = num.arccos(num.dot(ahat, chat))
    alfa = num.arccos(num.dot(bhat, chat))
    if outputDegrees:
        gama = r2d*gama
        beta = r2d*beta
        alfa = r2d*alfa

    return [a, b, c, alfa, beta, gama]

def latticePlanes(hkls, lparms, ltype='cubic', wavelength=1.54059292, strainMag=None):
    """
    Generates lattice plane data in the direct lattice for a given set
    of Miller indices.  Vector components are written in the
    crystal-relative RHON basis, X. The convention for fixing X to the
    lattice is such that a || x1 and c* || x3, where a and c* are
    direct and reciprocal lattice vectors, respectively.

    USAGE:

    planeInfo = latticePlanes(hkls, lparms, **kwargs)

    INPUTS:

    1) hkls (3 x n float ndarray) is the array of Miller indices for
       the planes of interest.  The vectors are assumed to be
       concatenated along the 1-axis (horizontal).

    2) lparms (1 x m float list) is the array of lattice parameters,
       where m depends on the symmetry group (see below).

    3) The following optional keyword arguments are recognized:

       *) ltype=(string) is a string representing the symmetry type of
          the implied Laue group.  The 11 available choices are shown
          below.  The default value is 'cubic'. Note that each group
          expects a lattice parameter array of the indicated length
          and order.

          latticeType      lparms
          -----------      ------------
          'cubic'          a
          'hexagonal'      a, c
          'trigonal'       a, c
          'rhombohedral'   a, alpha (in degrees)
          'tetragonal'     a, c
          'orthorhombic'   a, b, c
          'monoclinic'     a, b, c, beta (in degrees)
          'triclinic'      a, b, c, alpha, beta, gamma (in degrees)

       *) wavelength=(float) is a value represented the wavelength in
          Angstroms to calculate bragg angles for.  The default value
          is for Cu K-alpha radiation (1.54059292 Angstrom)

    OUTPUTS:

    1) planeInfo is a dictionary containing the following keys/items:

       normals   (3, n) double array    array of the components to the
                                        unit normals for each {hkl} in
                                        X (horizontally concatenated)

       dspacings (n,  ) double array    array of the d-spacings for
                                        each {hkl}

       2thetas   (n,  ) double array    array of the Bragg angles for
                                        each {hkl} relative to the
                                        specified wavelength

    NOTES:

    *) This function is effectively a wrapper to 'latticeVectors'.
       See 'help(latticeVectors)' for additional info.

    *) Lattice plane d-spacings are calculated from the reciprocal
       lattice vectors specified by {hkl} as shown in Appendix 1 of
       [1].

    REFERENCES:

    [1] B. D. Cullity, ``Elements of X-Ray Diffraction, 2
        ed.''. Addison-Wesley Publishing Company, Inc., 1978. ISBN
        0-201-01174-3

    """
    location = 'latticePlanes'

    assert hkls.shape[0] is 3, "hkls aren't column vectors in call to '%s'!" % (location)

    tag = ltype
    wlen = wavelength

    # get B
    L = latticeVectors(lparms, tag)

    # get G-vectors -- reciprocal vectors in crystal frame
    G = num.dot(L['B'], hkls)

    # magnitudes
    d = 1 / num.sqrt(num.sum(G**2, 0))

    angConv = 1.
    if outputDegrees:
        angConv = r2d
    # two thetas
    tth = angConv * 2 * num.arcsin(wlen / 2 / d)

    p = {'normals':unitVector(G),
         'dspacings':d,
         'tThetas':tth}

    if strainMag is not None:
        p['tThetasLo'] = angConv * 2 * num.arcsin(wlen / 2 / (d*(1.+strainMag)))
        p['tThetasHi'] = angConv * 2 * num.arcsin(wlen / 2 / (d*(1.-strainMag)))

    return p

def latticeVectors(lparms, tag='cubic', radians=False, debug=False):
    """
    Generates direct and reciprocal lattice vector components in a
    crystal-relative RHON basis, X. The convention for fixing X to the
    lattice is such that a || x1 and c* || x3, where a and c* are
    direct and reciprocal lattice vectors, respectively.

    USAGE:

    lattice = LatticeVectors(lparms, <symmTag>)

    INPUTS:

    1) lparms (1 x n float list) is the array of lattice parameters,
       where n depends on the symmetry group (see below).

    2) symTag (string) is a case-insensitive string representing the
       symmetry type of the implied Laue group.  The 11 available choices
       are shown below.  The default value is 'cubic'. Note that each
       group expects a lattice parameter array of the indicated length
       and order.

       latticeType      lparms
       -----------      ------------
       'cubic'          a
       'hexagonal'      a, c
       'trigonal'       a, c
       'rhombohedral'   a, alpha (in degrees)
       'tetragonal'     a, c
       'orthorhombic'   a, b, c
       'monoclinic'     a, b, c, beta (in degrees)
       'triclinic'      a, b, c, alpha, beta, gamma (in degrees)

    OUTPUTS:

    1) lattice is a dictionary containing the following keys/items:

       F         (3, 3) double array    transformation matrix taking
                                        componenents in the direct
                                        lattice (i.e. {uvw}) to the
                                        reference, X

       B         (3, 3) double array    transformation matrix taking
                                        componenents in the reciprocal
                                        lattice (i.e. {hkl}) to X

       BR        (3, 3) double array    transformation matrix taking
                                        componenents in the reciprocal
                                        lattice to the Risoe reference
                                        frame (see notes)

       U0        (3, 3) double array    transformation matrix
                                        (orthogonal) taking
                                        componenents in the
                                        Risoe reference frame to X

       vol       double                 the unit cell volume


       dparms    (6, ) double list      the direct lattice parameters:
                                        [a b c alpha beta gamma]

       rparms    (6, ) double list      the reciprocal lattice
                                        parameters:
                                        [a* b* c* alpha* beta* gamma*]

    NOTES:

    *) The conventions used for assigning a RHON basis,
       X -> {x1, x2, x3}, to each point group are consistent with
       those published in Appendix B of [1]. Namely: a || x1 and
       c* || x3.  This differs from the convention chosen by the Risoe
       group, where a* || x1 and c || x3 [2].

    *) The unit cell angles are defined as follows:
       alpha=acos(b'*c/|b||c|), beta=acos(c'*a/|c||a|), and
       gamma=acos(a'*b/|a||b|).

    *) The reciprocal lattice vectors are calculated using the
       crystallographic convention, where the prefactor of 2*pi is
       omitted. In this convention, the reciprocal lattice volume is
       1/V.

    *) Several relations from [3] were employed in the component
       calculations.

    REFERENCES:

    [1] J. F. Nye, ``Physical Properties of Crystals: Their
        Representation by Tensors and Matrices''. Oxford University
        Press, 1985. ISBN 0198511655

    [2] E. M. Lauridsen, S. Schmidt, R. M. Suter, and H. F. Poulsen,
        ``Tracking: a method for structural characterization of grains
        in powders or polycrystals''. J. Appl. Cryst. (2001). 34,
        744--750

    [3] R. J. Neustadt, F. W. Cagle, Jr., and J. Waser, ``Vector
        algebra and the relations between direct and reciprocal
        lattice quantities''. Acta Cryst. (1968), A24, 247--248


    """

    # build index for sorting out lattice parameters
    lattStrings = [
        'cubic'       ,
        'hexagonal'   ,
        'trigonal'    ,
        'rhombohedral',
        'tetragonal'  ,
        'orthorhombic',
        'monoclinic'  ,
        'triclinic'
        ]

    if radians:
        angConv = 1.
    else:
        angConv = pi/180. # degToRad
    deg90  = pi/2.
    deg120 = 2.*pi/3.
    #
    if tag == lattStrings[0]:   # cubic
        cellparms = num.r_[num.tile(lparms[0], (3,)), deg90*num.ones((3,))]
    elif tag == lattStrings[1] or tag == lattStrings[2]: # hexagonal | trigonal (hex indices)
        cellparms = num.r_[lparms[0], lparms[0], lparms[1], deg90, deg90, deg120]
    elif tag == lattStrings[3]: # rhombohedral
        cellparms = num.r_[num.tile(lparms[0], (3,)), num.tile(angConv*lparms[1], (3,))]
    elif tag == lattStrings[4]: # tetragonal
        cellparms = num.r_[lparms[0], lparms[0], lparms[1], deg90, deg90, deg90]
    elif tag == lattStrings[5]: # orthorhombic
        cellparms = num.r_[lparms[0], lparms[1], lparms[2], deg90, deg90, deg90]
    elif tag == lattStrings[6]: # monoclinic
        cellparms = num.r_[lparms[0], lparms[1], lparms[2], deg90, angConv*lparms[3], deg90]
    elif tag == lattStrings[7]: # triclinic
        cellparms = lparms
    else:
        raise RuntimeError('lattice tag \'%s\' is not recognized' % (tag))

    if debug:
        print str(cellparms[0:3]) + ' ' + str(r2d*cellparms[3:6])
    alfa = cellparms[3]
    beta = cellparms[4]
    gama = cellparms[5]

    cosalfar, sinalfar = cosineXform(alfa, beta, gama)

    a = cellparms[0]*num.r_[1, 0, 0]
    b = cellparms[1]*num.r_[num.cos(gama), num.sin(gama), 0]
    c = cellparms[2]*num.r_[num.cos(beta), -cosalfar*num.sin(beta), sinalfar*num.sin(beta)]

    ad = num.sqrt(sum(a**2))
    bd = num.sqrt(sum(b**2))
    cd = num.sqrt(sum(c**2))

    # Cell volume
    V = num.dot(a, num.cross(b, c))

    # F takes components in the direct lattice to X
    F = num.c_[a, b, c]

    # Reciprocal lattice vectors
    astar = num.cross(b, c)/V
    bstar = num.cross(c, a)/V
    cstar = num.cross(a, b)/V

    # and parameters
    ar = num.sqrt(sum(astar**2))
    br = num.sqrt(sum(bstar**2))
    cr = num.sqrt(sum(cstar**2))

    alfar = num.arccos(num.dot(bstar, cstar)/br/cr)
    betar = num.arccos(num.dot(cstar, astar)/cr/ar)
    gamar = num.arccos(num.dot(astar, bstar)/ar/br)

    # B takes components in the reciprocal lattice to X
    B = num.c_[astar, bstar, cstar]

    cosalfar2, sinalfar2 = cosineXform(alfar, betar, gamar)

    arisoe = ar*num.r_[1, 0, 0]
    brisoe = br*num.r_[num.cos(gamar), num.sin(gamar), 0]
    crisoe = cr*num.r_[num.cos(betar), -cosalfar2*num.sin(betar), sinalfar2*num.sin(betar)]

    BR = num.c_[arisoe, brisoe, crisoe]
    U0 = num.dot(B, num.linalg.inv(BR))
    if outputDegrees:
        dparms = num.r_[ad, bd, cd, r2d*num.r_[ alfa,  beta,  gama]]
        rparms = num.r_[ar, br, cr, r2d*num.r_[alfar, betar, gamar]]
    else:
        dparms = num.r_[ad, bd, cd, num.r_[ alfa,  beta,  gama]]
        rparms = num.r_[ar, br, cr, num.r_[alfar, betar, gamar]]

    L = {'F':F,
         'B':B,
         'BR':BR,
         'U0':U0,
         'vol':V,

         'dparms':dparms,
         'rparms':rparms}

    return L

def hexagonalIndicesFromRhombohedral(hkl):
    """
    converts rhombohedral hkl to hexagonal indices
    """
    HKL = num.zeros((3, hkl.shape[1]), dtype='int')

    HKL[0, :] = hkl[0, :] - hkl[1, :]
    HKL[1, :] = hkl[1, :] - hkl[2, :]
    HKL[2, :] = hkl[0, :] + hkl[1, :] + hkl[2, :]

    return HKL

def rhombohedralIndicesFromHexagonal(HKL):
    """
    converts hexagonal hkl to rhombohedral indices
    """
    hkl = num.zeros((3, HKL.shape[1]), dtype='int')

    hkl[0, :] = 2 * HKL[0, :] +     HKL[1, :] + HKL[2, :]
    hkl[1, :] =    -HKL[0, :] +     HKL[1, :] + HKL[2, :]
    hkl[2, :] =    -HKL[0, :] - 2 * HKL[1, :] + HKL[2, :]

    hkl = hkl / 3.
    return hkl

def rhombohedralParametersFromHexagonal(a_h, c_h):
    """
    converts hexagonal lattice parameters (a, c) to rhombohedral
    lattice parameters (a, alpha)
    """
    a_r    = num.sqrt(3 * a_h**2 + c_h**2) / 3.
    alfa_r = 2 * num.arcsin( 3. / (2 * num.sqrt(3 + (c_h / a_h)**2)))
    if outputDegrees:
        alfa_r = r2d * alfa_r
    return a_r, alfa_r

class PlaneData(object):
    """
    Careful with ordering:
    Outputs are ordered by the 2-theta for the hkl unless you get self.__hkls directly,
    and this order can change with changes in lattice parameters (lparms);
    setting and getting exclusions works on the current hkl ordering, not the original ordering (in self.__hkls), but exclusions are stored in the original ordering in case the hkl ordering does change with lattice parameters

    if not None, tThWidth takes priority over strainMag in setting two-theta ranges;
    changing strainMag automatically turns off tThWidth
    """

    def __init__(self,
                 hkls,
                 *args,
                 **kwargs):
        import symmetry as S

        self.phaseID = None
        self.__doTThSort = True
        self.__exclusions = None
        self.__tThMax = None
        #
        if len(args) == 4:
            lparms, laueGroup, wavelength, strainMag = args
            tThWidth = None
            self.__wavelength = processWavelength(wavelength)
            self.__lparms = self.__parseLParms(lparms)
        elif len(args) == 1 and hasattr(args[0],'getParams'):
            other = args[0]
            lparms, laueGroup, wavelength, strainMag, tThWidth = other.getParams()
            self.__wavelength = wavelength
            self.__lparms = lparms
            self.phaseID = other.phaseID
            self.__doTThSort = other.__doTThSort
            self.__exclusions = other.__exclusions
            self.__tThMax = other.__tThMax
            if hkls is None:
                hkls = other.__hkls
        else:
            raise NotImplementedError, 'args : '+str(args)

        self.__laueGroup = laueGroup
        self.__qsym = S.quatOfLaueGroup(self.__laueGroup)
        self.__hkls = copy.deepcopy(hkls)
        self.__strainMag = strainMag
        self.tThWidth = tThWidth

        # ... need to implement tThMin too
        if kwargs.has_key('phaseID'):
            self.phaseID = kwargs.pop('phaseID')
        if kwargs.has_key('doTThSort'):
            self.__doTThSort = kwargs.pop('doTThSort')
        if kwargs.has_key('exclusions'):
            self.__exclusions = kwargs.pop('exclusions')
        if kwargs.has_key('tThMax'):
            self.__tThMax = toFloat(kwargs.pop('tThMax'), 'radians')
        if kwargs.has_key('tThWidth'):
            self.tThWidth = kwargs.pop('tThWidth')
        if len(kwargs) > 0:
            raise RuntimeError, 'have unparsed keyword arguments with keys: '+str(kwargs.keys())

        self.__calc()

        return

    def __calc(self):
        import symmetry as S
        symmGroup = S.ltypeOfLaueGroup(self.__laueGroup)
        latPlaneData, latVecOps, hklDataList = PlaneData.makePlaneData(
            self.__hkls, self.__lparms, self.__qsym, symmGroup, self.__strainMag, self.wavelength)
        'sort by tTheta'
        tThs = num.array([hklDataList[iHKL]['tTheta'] for iHKL in range(len(hklDataList))])
        if self.__doTThSort:
            self.tThSort = num.argsort(tThs) # sorted hkl -> __hkl
            self.tThSortInv = num.empty(len(hklDataList), dtype=int) # __hkl -> sorted hkl
            self.tThSortInv[self.tThSort] = num.arange(len(hklDataList))
            self.hklDataList = [hklDataList[iHKL] for iHKL in self.tThSort]
        else:
            self.tThSort = num.arange(len(hklDataList))
            self.tThSortInv = num.arange(len(hklDataList))
            self.hklDataList = hklDataList
        self.__latVecOps = latVecOps
        self.nHKLs = len(self.getHKLs())
        return

    def __str__(self):
        s  = '========== plane data ==========\n'
        s += 'lattice parameters:\n   ' + str(self.lparms) +'\n'
        s += 'two theta width: (%s)\n'  % str(self.tThWidth)
        s += 'strain magnitude: (%s)\n' % str(self.strainMag)
        s += 'beam energy (%s)\n' % str(self.wavelength)
        s += 'hkls: (%d)\n' % self.nHKLs
        s += str(self.getHKLs())
        return s

    def getNHKLs(self):
        return self.nHKLs
    def getPhaseID(self):
        'may return None if not set'
        return self.phaseID

    def getParams(self):
        return (self.__lparms, self.__laueGroup, self.__wavelength, self.__strainMag, self.tThWidth)

    def getNhklRef(self):
        'does not use exclusions or the like'
        retval = len(self.hklDataList)
        return retval
    def get_hkls(self):
        """do not do return self.__hkls, as everywhere else hkls are returned in 2-theta order;
        transpose is to comply with lparm convention"""
        return self.getHKLs().T
    def set_hkls(self, hkls):
        raise RuntimeError,\
            'for now, not allowing hkls to be reset'
        #self.__exclusions = None
        #self.__hkls = hkls
        #self.__calc()
        return
    hkls = property(get_hkls, set_hkls, None)

    def get_tThMax(self):
        return self.__tThMax
    def set_tThMax(self, tThMax):
        self.__tThMax = toFloat(tThMax, 'radians')
        # self.__calc() # no need to redo calc for tThMax
        return
    tThMax = property(get_tThMax, set_tThMax, None)

    def get_exclusions(self):
        import copy
        retval = num.zeros(self.getNhklRef(), dtype=bool)
        if self.__exclusions is not None:
            'report in current hkl ordering'
            retval[:] = self.__exclusions[self.tThSortInv]
        if self.__tThMax is not None:
            for iHKLr, hklData in enumerate(self.hklDataList):
                if hklData['tTheta'] > self.__tThMax:
                    retval[iHKLr] = True
        return retval
    def set_exclusions(self, exclusions):
        excl = num.zeros(len(self.hklDataList), dtype=bool)
        if exclusions is not None:
            exclusions = num.atleast_1d(exclusions)
            if len(exclusions) == len(self.hklDataList):
                assert exclusions.dtype == 'bool', 'exclusions should be bool if full length'
                'convert from current hkl ordering to __hkl ordering'
                excl[:] = exclusions[self.tThSort]
            else:
                if len(exclusions.shape) == 1:
                    'treat exclusions as indices'
                    excl[self.tThSort[exclusions]] = True
                elif len(exclusions.shape) == 2:
                    raise NotImplementedError, 'have not yet coded treating exclusions as two-theta ranges'
                else:
                    raise RuntimeError, \
                        'do not now what to do with exclusions with shape '+str(exclusions.shape)
        self.__exclusions = excl
        self.nHKLs = num.sum(num.logical_not(self.__exclusions))
        return
    exclusions = property(get_exclusions, set_exclusions, None)

    def get_lparms(self):
        return self.__lparms
    def __parseLParms(self, lparms):
        lparmsDUnit = []
        for lparmThis in lparms:
            if hasattr(lparmThis, 'getVal'):
                if lparmThis.isLength():
                    lparmsDUnit.append(lparmThis.getVal(dUnit))
                elif lparmThis.isAngle():
                    'plumbing set up to default to degrees for lattice parameters'
                    lparmsDUnit.append(lparmThis.getVal('degrees'))
                else:
                    raise RuntimeError, 'do not know what to do with '+str(lparmThis)
            else:
                lparmsDUnit.append(lparmThis)
        return lparmsDUnit
    def set_lparms(self, lparms):
        self.__lparms = self.__parseLParms(lparms)
        self.__calc()
        return
    lparms = property(get_lparms, set_lparms, None)

    def get_strainMag(self):
        return self.__strainMag
    def set_strainMag(self, strainMag):
        self.__strainMag = strainMag
        self.tThWidth = None
        self.__calc()
        return
    strainMag = property(get_strainMag, set_strainMag, None)

    def get_wavelength(self):
        return self.__wavelength
    def set_wavelength(self, wavelength):
        self.__wavelength = processWavelength(wavelength)
        self.__calc()
        return
    wavelength = property(get_wavelength, set_wavelength, None)

    @staticmethod
    def makePlaneData(hkls, lparms, qsym, symmGroup, strainMag, wavelength): # spots
        """
        hkls       : need to work with crystallography.latticePlanes
        lparms     : need to work with crystallography.latticePlanes
        laueGroup  : see symmetry module
        wavelength : wavelength
        strainMag  : swag of strian magnitudes
        """
        import crystallography
        crystallography.tempSetOutputDegrees(False)
        import symmetry as S

        latPlaneData = crystallography.latticePlanes(hkls, lparms, ltype=symmGroup,
                                           strainMag=strainMag,
                                           wavelength=wavelength)

        latVecOps = latticeVectors(lparms, symmGroup)

        hklDataList = []
        for iHKL in range(len(hkls.T)): # need transpose because of convention for hkls ordering

            # JVB # latVec = latPlaneData['normals'][:,iHKL]
            # JVB # # ... if not spots, may be able to work with a subset of these
            # JVB # latPlnNrmlList = S.applySym(num.c_[latVec], qsym, csFlag=True, cullPM=False)

            # returns UN-NORMALIZED lattice plane normals
            latPlnNrmls = S.applySym(
                num.dot(latVecOps['B'], hkls[:,iHKL].reshape(3, 1)),
                qsym,
                csFlag=True,
                cullPM=False)

            # check for +/- in symmetry group
            latPlnNrmlsM = S.applySym(
                num.dot(latVecOps['B'], hkls[:,iHKL].reshape(3, 1)),
                qsym,
                csFlag=False,
                cullPM=False)

            csRefl = latPlnNrmls.shape[1] == latPlnNrmlsM.shape[1]

            # added this so that I retain the actual symmetric integer hkls as well
            symHKLs = num.array( num.round( num.dot(latVecOps['F'].T, latPlnNrmls) ), dtype='int' )

            hklDataList.append({
                'hklID'       : iHKL                           ,
                'hkl'         : hkls[:,iHKL]                   ,
                'tTheta'      : latPlaneData['tThetas'][iHKL]  ,
                'dSpacings'   : latPlaneData['dspacings'][iHKL],
                'tThetaLo'    : latPlaneData['tThetasLo'][iHKL],
                'tThetaHi'    : latPlaneData['tThetasHi'][iHKL],
                'latPlnNrmls' : unitVector(latPlnNrmls)        ,
                'symHKLs'     : symHKLs,
                'centrosym'   : csRefl
                })

        crystallography.revertOutputDegrees()
        return latPlaneData, latVecOps, hklDataList
    def getLatticeType(self):
        """This is the lattice type"""
        import symmetry as S
        return S.ltypeOfLaueGroup(self.__laueGroup)
    def getLaueGroup(self):
        """This is the Schoenflies tag"""
        return self.__laueGroup
    def getQSym(self):
        #import symmetry as S
        return self.__qsym # S.quatOfLaueGroup(self.__laueGroup)
    def getPlaneSpacings(self):
        """
        gets plane spacings
        """
        dspacings = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue
            dspacings.append(hklData['dSpacings'])
        return dspacings
    def getPlaneNormals(self):
        """
        gets both +(hkl) and -(hkl) normals
        """
        plnNrmls = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue
            plnNrmls.append(hklData['latPlnNrmls'])
        return plnNrmls

    def getLatticeOperators(self):
        """
        gets lattice vector operators as a new (deepcopy)
        """
        import copy
        return copy.deepcopy(self.__latVecOps)
    def setLatticeOperators(self, val):
        raise RuntimeError, 'do not set latVecOps directly, change other things instead'
    latVecOps = property(getLatticeOperators, setLatticeOperators, None)

    def __thisHKL(self,iHKLr):
        retval = True
        hklData   = self.hklDataList[iHKLr]
        if self.__exclusions is not None:
            if self.__exclusions[self.tThSortInv[iHKLr]]:
                retval = False
        if self.__tThMax is not None:
            # if hklData['tThetaHi'] > self.__tThMax: continue
            if hklData['tTheta'] > self.__tThMax:
                retval = False
        return retval
    def __getTThRange(self,iHKLr):
        hklData = self.hklDataList[iHKLr]
        if self.tThWidth is not None: # tThHi-tThLo < self.tThWidth
            tTh   = hklData['tTheta']
            tThHi = tTh + self.tThWidth * 0.5
            tThLo = tTh - self.tThWidth * 0.5
        else:
            tThHi = hklData['tThetaHi']
            tThLo = hklData['tThetaLo']
        return (tThLo, tThHi)

    def getTThRanges(self, strainMag=None, lparms=None):
        """Return 2-theta ranges for included hkls

        return array is n x 2
        """
        if lparms is None:
            tThRanges = []
            for iHKLr, hklData in enumerate(self.hklDataList):
                if not self.__thisHKL(iHKLr): continue
                #tThRanges.append([hklData['tThetaLo'], hklData['tThetaHi']])
                if strainMag is None:
                    tThRanges.append(self.__getTThRange(iHKLr))
                else:
                    hklData = self.hklDataList[iHKLr]
                    d = hklData['dSpacings']
                    tThLo = 2.0 * num.arcsin(self.__wavelength / 2.0 / (d*(1.+strainMag)))
                    tThHi = 2.0 * num.arcsin(self.__wavelength / 2.0 / (d*(1.-strainMag)))
                    tThRanges.append((tThLo, tThHi))
        else:
            new = self.__class__(self.__hkls, self)
            new.lparms = lparms
            tThRanges = new.getTThRanges(strainMag=strainMag)
        return num.array(tThRanges)
    def makeNew(self):
        new = self.__class__(None, self)
        return new
    def getTTh(self, lparms=None):
        if lparms is None:
            tTh = []
            for iHKLr, hklData in enumerate(self.hklDataList):
                if not self.__thisHKL(iHKLr): continue
                tTh.append(hklData['tTheta'])
        else:
            new = self.makeNew()
            new.lparms = lparms
            tTh = new.getTTh()
        return num.array(tTh)
    def getDD_tThs_lparms(self):
        """
        derivatives of tThs with respect to lattice parameters;
        have not yet done coding for analytic derivatives, just wimp out and finite difference
        """
        import copy

        pert = 1.0e-5 # assume they are all around unity
        pertInv = 1.0/pert

        lparmsRef = copy.deepcopy(self.__lparms)
        tThRef    = self.getTTh()
        ddtTh     = num.empty((len(tThRef), len(lparmsRef)))

        for iLparm in range(len(lparmsRef)):
            self.__lparms = copy.deepcopy(lparmsRef)
            self.__lparms[iLparm] += pert
            self.__calc()

            iTTh = 0
            for iHKLr, hklData in enumerate(self.hklDataList):
                if not self.__thisHKL(iHKLr): continue
                ddtTh[iTTh, iLparm] = (hklData['tTheta'] - tThRef[iTTh]) * pertInv
                iTTh += 1

        'restore'
        self.__lparms = lparmsRef
        self.__calc()

        return ddtTh

    def getMultiplicity(self):          # ... JVB: is this incorrect?
        multip = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue
            multip.append( hklData['symHKLs'].shape[1] )
        return num.array(multip)

    def getHKLID(self, hkl):
        'can call on a single hkl or list of hkls'
        if hasattr(hkl,'__setitem__'): # tuple does not have __setitem__
            if hasattr(hkl,'shape'):
                'if is ndarray, assume is 3xN'
                retval = map(self.__getHKLID, hkl.T)
            else:
                retval = map(self.__getHKLID, hkl)
        else:
            retval = self.__getHKLID(hkl)
        return retval
    def __getHKLID(self, hkl):
        """
        for hkl that is a tuple, return externally visible hkl index
        """
        if isinstance(hkl,int):
            retval = hkl
        else:
            hklList = self.getHKLs().tolist()
            dHKLInv   = dict([[tuple(hklThis),iHKL] for iHKL, hklThis in enumerate(hklList)])
            retval = dHKLInv[tuple(hkl)]
        return retval
    def getHKLs(self, asStr=False, thisTTh=None, allHKLs=False):
        """
        if pass thisTTh, then only return hkls overlapping the specified 2-theta;
        if set allHKLs to true, the ignore exlcusions, tThMax, etc
        """
        hkls = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not allHKLs:
                if not self.__thisHKL(iHKLr):
                    continue
            if thisTTh is not None:
                tThLo, tThHi = self.__getTThRange(iHKLr)
                if thisTTh < tThHi and thisTTh > tThLo:
                    hkls.append(hklData['hkl'])
            else:
                hkls.append(hklData['hkl'])
        if asStr:
            retval = map(hklToStr, num.array(hkls))
        else:
            retval = num.array(hkls)
        return retval
    def getSymHKLs(self, asStr=False, indices=None):
        """
        new function that returns all symmetric hkls
        """
        retval = []
        iRetval = 0
        if indices is not None:
            indB = num.zeros(self.nHKLs,dtype=bool)
            indB[num.array(indices)] = True
        else:
            indB = num.ones(self.nHKLs,dtype=bool)
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue
            if indB[iRetval]:
                hkls = hklData['symHKLs']
                if asStr:
                    myStr = lambda x: re.sub('\[|\]|\(|\)','',str(x))
                    retval.append(map(myStr, num.array(hkls).T))
                else:
                    retval.append(num.array(hkls))
            iRetval += 1
        return retval
    def getCentroSymHKLs(self):
        retval = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue
            retval.append( hklData['centrosym'] )
        return retval
    def makeRecipVectors(self,**RUVorF):
        """
        input keywords: specify enough to construct a deformation gradient. Accepts R, U, V, F. If no arguments are given defaults to identity.
        uses F = RU = VR, and F^{-T}gI = rI to construct recip vectors.
        makeRecipVectors(R = [R])
        makeRecipVectors(U = [U])
        makeRecipVectors(V = [V])
        makeRecipVectors(R = [R], U = [U])
        makeRecipVectors(R = [R], U = [U])
        makeRecipVectors(F = [F])
        returns 3x~200 array with columns as the recip Vectors in the sample frame
        """

        if RUVorF.has_key('F'):
            F = RUVorF['F']
        else:
            if RUVorF.has_key('R'):
                R = RUVorF['R']
            else:
                R = num.eye(3)

            if RUVorF.has_key('U'):
                U = RUVorF['U']
                F = num.dot(R,U)
            elif RUVorF.has_key('V'):
                V = RUVorF['V']
                F = num.dot(V,R)
            else:
                U = num.eye(3)
                F = num.dot(R,U)

        from scipy.linalg import inv
        rMat = num.eye(3)
        Qvec, predQAng0,predQAng1 = self.__makeScatteringVectors(rMat, bMat=None)

        Qvecs = []
        FinvT = inv(F).T

        n_cols = 0
        for Qblock in Qvec:
            rows,cols = Qblock.shape
            n_cols+=cols
        rIs = num.zeros([rows,n_cols])
        col = 0
        for Qblock in Qvec:
            ncols = Qblock.shape[1]
            rIs_ = num.dot(FinvT,Qblock)
            rIs[0:rows,col:col + Qblock.shape[1]] = rIs_[:,:]
            col += Qblock.shape[1]

        return rIs

    def makeTheseScatteringVectors(self, hklList, rMat, bMat=None, wavelength=None, chiTilt=None):
        iHKLList = num.atleast_1d(self.getHKLID(hklList))
        fHKLs = num.hstack(self.getSymHKLs(indices=iHKLList))
        if bMat is None:
            bMat = self.__latVecOps['B']
        if wavelength is None:
            wavelength = self.__wavelength
        retval = PlaneData.makeScatteringVectors(fHKLs, rMat, bMat, wavelength, chiTilt=chiTilt)
        return retval
    def makeAllScatteringVectors(self, rMat, bMat=None, wavelength=None, chiTilt=None):
        fHKLs = num.hstack(self.getSymHKLs())
        if bMat is None:
            bMat = self.__latVecOps['B']
        if wavelength is None:
            wavelength = self.__wavelength
        retval = PlaneData.makeScatteringVectors(fHKLs, rMat, bMat, wavelength, chiTilt=chiTilt)
        return retval
    @staticmethod
    def makeScatteringVectors(hkls, rMat, bMat, wavelength, chiTilt=None):
        """
        modeled after QFromU.m
        """
        from Rotations import angleAxisOfRotMat, rotMatOfExpMap, arccosSafe, mapAngle

        # basis vectors
        Xl = num.vstack([1, 0, 0])          # X in the lab frame
        Yl = num.vstack([0, 1, 0])          # Y in the lab frame
        Zl = num.vstack([0, 0, 1])          # Z in the lab frame

        chiRotMat = num.eye(3)
        if chiTilt is not None:
            chiRotMat = rotMatOfExpMap(chiTilt*Yl)

        # projection operators
        Pxy = num.eye(3) - num.dot(Zl, Zl.T)    # xy-plane
        Pxz = num.eye(3) - num.dot(Yl, Yl.T)    # xz-plane

        zTol = 1.0e-7                       # zero tolerance for checking vectors

        Qs_vec  = []
        Qs_ang0 = []
        Qs_ang1 = []

        # grab invariants of rMat
        phi, n = angleAxisOfRotMat(rMat)

        # these are the scattering vectors in the CRYSTAL FRAME
        Qc = num.dot( bMat, hkls )

        dummy, nRefl = Qc.shape
        assert dummy == 3, "Looks like something is wrong with your lattice plane normals son!"

        Qc_mag = columnNorm(Qc)
        Qc_hat = unitVector(Qc)

        # construct bragg angle
        dSpacing = 1. / Qc_mag
        tht      = num.arcsin( wavelength / 2. / dSpacing  )

        # move to sample frame
        Qs_hat = num.dot(chiRotMat, num.dot(rMat.squeeze(), Qc_hat))

        # first find Q's that can satisfy the bragg condition
        dtplus  = num.dot( Yl.T, Qs_hat)
        dtminus = num.dot(-Yl.T, Qs_hat)

        keepers = ( (dtplus  <= (num.cos(tht) - zTol)) & (dtplus  >= 0.0) ) | \
                  ( (dtminus <= (num.cos(tht) - zTol)) & (dtminus >= 0.0) )

        bangOn  = ( (dtplus  > (num.cos(tht) - zTol)) & (dtplus  <= (num.cos(tht) + zTol)) ) | \
                  ( (dtminus > (num.cos(tht) - zTol)) & (dtminus <= (num.cos(tht) + zTol)) )

        if nRefl == 1:
            keepers = keepers.reshape(1) # these guys can diffract
            bangOn  = bangOn.reshape(1)  # these guys are tangent to the cone
        else:
            keepers = keepers.squeeze()  # these guys can diffract
            bangOn  = bangOn.squeeze()   # these guys are tangent to the cone
            pass

        #############################################################################
        # MUST CALCULATE CORRESPONDING ANGULAR COORDINATES
        #
        # After some algebra, the equation can be reduced to something of the form:
        #
        #                       a*sin(x) + b*cos(x) = c
        #
        # which has two unique solutions UNLESS Qs happens to be tangent to the cone.
        # In that case, x is a double root.
        cphi = num.cos(phi)
        sphi = num.sin(phi)
        if chiTilt is None:
            cchi = 1.
            schi = 0.
        else:
            cchi = num.cos(chiTilt)
            schi = num.sin(chiTilt)
            pass
        nchi = num.c_[0, cchi, schi].T
        Pchi = num.eye(3) - num.dot(nchi, nchi.T)

        # a = cchi * ( n[0]**2   * (cphi - 1) - cphi      ) * Qc_hat[0, :] \
        #          + ( n[0]*n[1] * (cphi - 1) + sphi*n[2] ) * Qc_hat[1, :] \
        #          + ( n[0]*n[2] * (cphi - 1) - sphi*n[1] ) * Qc_hat[2, :]
        #
        # b = cchi * ( n[0]*n[2] * (1 - cphi) - sphi*n[1] ) * Qc_hat[0, :] \
        #          + ( n[1]*n[2] * (1 - cphi) + sphi*n[0] ) * Qc_hat[1, :] \
        #          + ( n[2]**2   * (1 - cphi) + cphi      ) * Qc_hat[2, :]
        #
        # c = num.sin(tht) \
        #         - schi * ( ( (1 - cphi)*n[0]*n[1] + sphi*n[2] ) * Qc_hat[0, :] \
        #                  + ( (1 - cphi)*n[1]**2   + cphi      ) * Qc_hat[1, :] \
        #                  + ( (1 - cphi)*n[1]*n[2] - sphi*n[0] ) * Qc_hat[2, :] )

        a = -cchi * Qs_hat[0, :]
        b =  cchi * Qs_hat[2, :]
        c =  num.sin(tht) - schi * Qs_hat[1, :]

        onesNRefl = num.ones(nRefl)

        ome0 = num.nan * onesNRefl
        ome1 = num.nan * onesNRefl
        eta0 = num.nan * onesNRefl
        eta1 = num.nan * onesNRefl
        for iRefl in range(nRefl):
            if keepers[iRefl]:
                abDist = num.sqrt(a[iRefl]*a[iRefl] + b[iRefl]*b[iRefl])
                abAngl = num.arctan2(b[iRefl], a[iRefl])
                cdAngl = num.arcsin( c[iRefl] / abDist )

                ome0[iRefl] = float( mapAngle( cdAngl - abAngl, units='radians' ) )
                ome1[iRefl] = float( mapAngle( num.pi - cdAngl - abAngl, units='radians' ) )
            elif bangOn[iRefl]:
                # ... needs checking
                if chiTilt is not None and chiTilt !=0:
                    qxz_p = unitVector( num.dot(Pchi, Qs_hat[:,iRefl]).reshape(3, 1) )
                    Zl_p  = unitVector( num.dot(Pchi, Zl).reshape(3, 1) )
                    if abs(qxz_p[0] > zTol):
                        tmp = -qxz_p[0] / abs(qxz_p[0]) * arccosSafe(qxz_p[-1])
                    else:
                        tmp = arccosSafe(qxz_p[-1])
                else:
                    qxz = unitVector( Qs_hat[[0,2],iRefl].reshape(2, 1) )
                    if abs(qxz[0]) > zTol:
                        tmp = -qxz[0] / abs(qxz[0]) * arccosSafe(qxz[1])
                    else:
                        tmp = arccosSafe(qxz[1])
                        pass
                    pass

                tmp = float( mapAngle( tmp, units='radians' ) )

                # trick ome range filter into selecting just one...
                ome0[iRefl] = tmp
                ome1[iRefl] = tmp + 2*num.pi
                pass # close conditional on bangOn reflections
            if not num.isnan(ome0[iRefl]):
                # MUST NORMALIZE projected vectors to use arccosSafe!
                cOme0 = num.cos(ome0[iRefl])
                sOme0 = num.sin(ome0[iRefl])
                cOme1 = num.cos(ome1[iRefl])
                sOme1 = num.sin(ome1[iRefl])

                # get X-Y projections in lab frame
                qxy0 = unitVector( num.vstack([
                    cOme0*Qs_hat[0,iRefl] + sOme0*Qs_hat[2,iRefl],
                    schi*sOme0*Qs_hat[0,iRefl] + cchi*Qs_hat[1,iRefl] - schi*cOme0*Qs_hat[2,iRefl]]
                                              ) )
                qxy1 = unitVector( num.vstack([
                    cOme1*Qs_hat[0,iRefl] + sOme1*Qs_hat[2,iRefl],
                    schi*sOme1*Qs_hat[0,iRefl] + cchi*Qs_hat[1,iRefl] - schi*cOme1*Qs_hat[2,iRefl]]
                                              ) )
                # eta0
                if abs(qxy0[1]) > zTol:
                    eta0[iRefl] = qxy0[1] / abs(qxy0[1]) * arccosSafe(qxy0[0])
                else:
                    eta0[iRefl] = arccosSafe(qxy0[0])
                    pass
                # eta1
                if abs(qxy1[1]) > zTol:
                    eta1[iRefl] = qxy1[1] / abs(qxy1[1]) * arccosSafe(qxy1[0])
                else:
                    eta1[iRefl] = arccosSafe(qxy1[0])
                    pass
                pass # close conditional on ome0
            pass # close loop on nrefl
        Qs_vec  = num.tile(1./dSpacing, (3, 1)) * Qs_hat
        Qs_ang0 = num.vstack([2*tht, eta0, ome0])
        Qs_ang1 = num.vstack([2*tht, eta1, ome1])

        return Qs_vec, Qs_ang0, Qs_ang1

    def __makeScatteringVectors(self, rMat, bMat=None, chiTilt=None):
        """
        modeled after QFromU.m
        """
        from Rotations import angleAxisOfRotMat, rotMatOfExpMap, arccosSafe, mapAngle

        if bMat is None:
            bMat = self.__latVecOps['B']

        Qs_vec  = []
        Qs_ang0 = []
        Qs_ang1 = []
        for iHKLr, hklData in enumerate(self.hklDataList):
            if not self.__thisHKL(iHKLr): continue

            thisQs, thisAng0, thisAng1 = \
                    PlaneData.makeScatteringVectors(hklData['symHKLs'], rMat, bMat,
                                                    self.__wavelength, chiTilt=chiTilt)

            Qs_vec.append( thisQs )
            Qs_ang0.append( thisAng0 )
            Qs_ang1.append( thisAng1 )

        return Qs_vec, Qs_ang0, Qs_ang1

def getFriedelPair(tth0, eta0, *ome0, **kwargs):
    """
    Get the diffractometer angular coordinates in degrees for
    the Friedel pair of a given reflection (min angular distance).

    AUTHORS:

    J. V. Bernier -- 10 Nov 2009

    USAGE:

    ome1, eta1 = getFriedelPair(tth0, eta0, *ome0,
                                display=False,
                                units='degrees',
                                convention='aps')

    INPUTS:

    1) tth0 is a list (or ndarray) of 1 or n the bragg angles (2theta) for
       the n reflections (tiled to match eta0 if only 1 is given).

    2) eta0 is a list (or ndarray) of 1 or n azimuthal coordinates for the n
       reflections  (tiled to match tth0 if only 1 is given).

    3) ome0 is a list (or ndarray) of 1 or n reference oscillation
       angles for the n reflections (denoted omega in [1]).  This argument
       is optional.

    4) Keyword arguments may be one of the following:

    Keyword             Values|{default}        Action
    --------------      --------------          --------------
    'display'           True|{False}            toggles display info to cmd line
    'units'             'radians'|{'degrees'}   sets units for input angles
    'convention'        'risoe'|{'aps'}         sets conventions defining
                                                the angles (see below)
    'chiTilt'           None                    the inclination (about Xlab) of
                                                the oscillation axis

    OUTPUTS:

    1) ome1 contains the oscialltion angle coordinates of the
       Friedel pairs associated with the n input reflections, relative to ome0
       (i.e. ome1 = <result> + ome0).  Output is in DEGREES!

    2) eta1 contains the azimuthal coordinates of the Friedel
       pairs associated with the n input reflections.  Output units are
       controlled via the module variable 'outputDegrees'

    NOTES:

    JVB) The ouputs ome1, eta1 are written using the selected convention, but the
         units are alway degrees.  May change this to work with Nathan's global...

    JVB) In the 'risoe' convention [1], {XYZ} form a RHON basis where X is
         downstream, Z is vertical, and eta is CCW with +Z defining eta = 0.

    JVB) In the 'aps' convention [2], {XYZ} form a RHON basis where Z is upstream,
         Y is vertical, and eta is CCW with +X defining eta = 0.

    REFERENCES:

    [1] E. M. Lauridsen, S. Schmidt, R. M. Suter, and H. F. Poulsen,
        ``Tracking: a method for structural characterization of grains in
        powders or polycrystals''. J. Appl. Cryst. (2001). 34, 744--750

    [2] J. V. Bernier, M. P. Miller, J. -S. Park, and U. Lienert,
        ``Quantitative Stress Analysis of Recrystallized OFHC Cu Subject
        to Deformed In Situ'', J. Eng. Mater. Technol. (2008). 130.
        DOI:10.1115/1.2870234
    """
    from Rotations import arccosSafe, mapAngle, rotMatOfExpMap

    dispFlag  = False
    risoeFlag = False
    chiTilt   = None
    c1        = 1.
    c2        = pi/180.
    zTol      = 1.e-7

    # cast to arrays (in case they aren't)
    if num.isscalar(eta0):
        eta0 = [eta0]

    if num.isscalar(tth0):
        tth0 = [tth0]

    if num.isscalar(ome0):
        ome0 = [ome0]

    eta0 = num.asarray(eta0)
    tth0 = num.asarray(tth0)
    ome0 = num.asarray(ome0)

    if eta0.ndim != 1:
        raise RuntimeError, 'your azimutal input was not 1-D, so I do not know what you expect me to do'

    npts = len(eta0)

    if tth0.ndim != 1:
        raise RuntimeError, 'your Bragg angle input was not 1-D, so I do not know what you expect me to do'
    else:
        if len(tth0) != npts:
            if len(tth0) == 1:
                tth0 = tth0*num.ones(npts)
            elif npts == 1:
                npts = len(tth0)
                eta0 = eta0*num.ones(npts)
            else:
                raise RuntimeError, 'the azimuthal and Bragg angle inputs are inconsistent'

    if len(ome0) == 0:
        ome0 = num.zeros(npts)                  # dummy ome0
    elif len(ome0) == 1 and npts > 1:
        ome0 = ome0*num.ones(npts)
    else:
        if len(ome0) != npts:
            raise RuntimeError('your oscialltion angle input is inconsistent; ' \
                               + 'it has length %d while it should be %d' % (len(ome0), npts) )

    ome1 = num.nan*num.ones(npts)
    eta1 = num.nan*num.ones(npts)

    # keyword args processing
    kwarglen = len(kwargs)
    if kwarglen > 0:
        argkeys = kwargs.keys()
        for i in range(kwarglen):
            if argkeys[i] == 'display':
                dispFlag = kwargs[argkeys[i]]
            elif argkeys[i] == 'convention':
                if kwargs[argkeys[i]].lower() == 'risoe':
                    risoeFlag = True
            elif argkeys[i] == 'units':
                if kwargs[argkeys[i]] == 'radians':
                    c1 = 180./pi
                    c2 = 1.
            elif argkeys[i] == 'chiTilt':
                if kwargs[argkeys[i]] is not None:
                    chiTilt = kwargs[argkeys[i]]

    # a little talkback...
    if dispFlag:
        if risoeFlag:
            print '\nUsing Risoe angle convention\n'
        else:
            print '\nUsing image-based angle convention\n'

    # mapped eta input
    #   - in DEGREES, thanks to c1
    eta0 = mapAngle(c1*eta0, [-180, 180], units='degrees')
    if risoeFlag:
        eta0  = 90 - eta0

    # must put args into RADIANS
    #   - eta0 is in DEGREES,
    #   - the others are in whatever was entered, hence c2
    eta0  = d2r*eta0
    tht0  = c2*tth0/2
    if chiTilt is not None:
        chiTilt = c2*chiTilt

    # ---------------------
    # SYSTEM SOLVE
    #
    #
    # cos(chi)cos(eta)cos(theta)sin(x) - cos(chi)sin(theta)cos(x) = sin(theta) - sin(chi)sin(eta)cos(theta)
    #
    #
    # Identity: a sin x + b cos x = sqrt(a**2 + b**2) sin (x + alfa)
    #
    #       /
    #       |      atan(b/a) for a > 0
    # alfa <
    #       | pi + atan(b/a) for a < 0
    #       \
    #
    # => sin (x + alfa) = c / sqrt(a**2 + b**2)
    #
    # must use both branches for sin(x) = n: x = u (+ 2k*pi) | x = pi - u (+ 2k*pi)
    #
    cEta0 = num.cos(eta0)
    sEta0 = num.sin(eta0)
    cTht0 = num.cos(tht0)
    sTht0 = num.sin(tht0)
    if chiTilt is None or chiTilt == 0:
        cChi  = 1.
        cChi2 = 1.
        sChi  = 0.
        sChi2 = 0.
        s2Chi = 0.
    else:
        cChi  = num.cos(chiTilt)
        cChi2 = cChi**2
        sChi  = num.sin(chiTilt)
        sChi2 = sChi**2
        s2Chi = num.sin(2*chiTilt)
        pass

    a = cChi*cEta0*cTht0
    b = (0.5*s2Chi*sEta0*cTht0 - cChi2*sTht0)
    c = (1 + sChi2)*sTht0 + 0.5*s2Chi*sEta0*cTht0

    agt0 = a >= 0
    alt0 = a <  0

    ome1[agt0] =      num.arcsin( c[agt0] / num.sqrt(a[agt0]**2 + b[agt0]**2) ) - num.arctan2(b[agt0], a[agt0])
    ome1[alt0] = pi - num.arcsin( c[alt0] / num.sqrt(a[alt0]**2 + b[alt0]**2) ) - num.arctan2(b[alt0], a[alt0])

    # if chiTilt isn't 0 (None) then the solution for eta1 is not trivial
    if chiTilt is not None and chiTilt != 0:
        cOme1 = num.cos(ome1)
        sOme1 = num.sin(ome1)
        qxy = unitVector( -num.vstack(
            [ cTht0*cEta0*cOme1 - cTht0*sEta0*sOme1*sChi + sTht0*sOme1*cChi,
              cTht0*cEta0*sOme1*sChi + cTht0*sEta0*(cChi2 + cOme1*sChi2) + sTht0*(1 - cOme1)*sChi*cChi ] ) )

        eta1 = qxy[1, :] / abs(qxy[1, :]) * arccosSafe(qxy[0, :]) * r2d

        if risoeFlag:
            eta1  =  90 - eta1
            pass
    else:
        eta1 = r2d*eta0
        if risoeFlag:
            eta1  = 270 - eta1
        else:
            eta1  = 180 + eta1
            pass
        pass

    # everybody back to DEGREES!
    #     - ome1 is in RADIANS here
    #     - convert and put into [-180, 180]
    ome1 = mapAngle( mapAngle(r2d*ome1, [-180, 180], units='degrees') + c1*ome0, [-180, 180], units='degrees')

    # DEBUGGING # qfs = num.vstack([num.cos(d2r*eta1)*cTht0, num.sin(d2r*eta1)*cTht0, sTht0])
    # DEBUGGING # rsd = num.dot(qfs.T, num.c_[0,0,1].T) - sTht0
    # DEBUGGING # print rsd, num.min(rsd[-num.isnan(rsd)])

    # put eta1 in [-180, 180]
    eta1 = mapAngle(eta1, [-180, 180], units='degrees')

    if not outputDegrees:
        ome1 = d2r * ome1
        eta1 = d2r * eta1

    return ome1, eta1

def getDparms(lp, lpTag, radians=True):
    """
    Utility routine for getting dparms, that is the lattice parameters without symmetry -- 'triclinic'
    """
    latVecOps = latticeVectors(lp, tag=lpTag, radians=radians)
    return latVecOps['dparms']

# def computeRLVectors():
#     """
#     """
#
#     return

# def findFriedelPairs(ome, eta, tth, tol=(1.0, 0.2, 0.05), units='degrees'):
#     """
#     A wrapper for getFriedelPair to operate on scattering vector data
#     """
#     tolType = 'angular'
#     if num.isscalar(tol) or len(tol) == 1:
#         tolType = 'cosine'
#
#     c1        = 1.
#     c2        = pi/180.
#     zTol      = 1.e-7
#
#     # cast to arrays (in case they aren't)
#     if num.isscalar(eta0):
#         eta0 = [eta0]
#
#     if num.isscalar(tth0):
#         tth0 = [tth0]
#
#     if num.isscalar(ome0):
#         ome0 = [ome0]
#
#     eta0 = num.asarray(eta0)
#     tth0 = num.asarray(tth0)
#     ome0 = num.asarray(ome0)
#
#     npts = len(eta)
#
# Fpairs = [];
# if numPts > 1
#
#     % search for Friedel pairs
#     [ome_f, eta_f] = GetFriedelPair(eta, tth, ome);
#
#     ome_f = ome_f/r2d;
#     eta_f = eta_f/r2d;
#     th    = th/r2d;
#
#     R_l2s_f = RMatOfQuat(QuatOfAngleAxis(-ome_f, [0 1 0]'));
#
#     % idealized friedel pairs in LAB FRAME
#     Qhat_f = [...
#         cos(eta_f).*cos(th), ...
#         sin(eta_f).*cos(th), ...
#         sin(th)]';
#
#     % calculated friedel pairs in sample frame
#     Qf = reshape(SparseOfMatArray(R_l2s_f)*Qhat_f(:), [3, numPts]);
#
#     jj = 1;
#     skip = [];
#     for j = 1:numPts
#         if ~ismember(j, skip)
#
#             fdots = Qf(:, j)' * UnitVector(Q);
#             hits = find(fdots >= dTol);
#
#             if ~isempty(hits)
#                 if length(hits) > 1
#                     [best, ibest] = max(fdots(hits));
#
#                     Fpairs(jj, :) = [j, hits(ibest), best, 1];
#                     skip = cat(1, hits(ibest), skip);
#                 else
#                     Fpairs(jj, :) = [j, hits, fdots(hits), 0];
#                     skip = cat(1, hits, skip);
#                 end
#                 jj = jj + 1;
#             end
#         end
#     end
# end
#

"""
UNIT TESTS
"""
# from az31 import planeData as mgData
# for friedel pairs...  a start
# rMat = rot.rotMatOfQuat(mUtil.unitVector(num.random.randn(4, 1)))
# qv = mgData._PlaneData__makeScatteringVectors(rMat, chiTilt=0.25)
# xtal.getFriedelPair(qv[1][0][0, 0], qv[1][0][1, 0], qv[1][0][2, 0], chiTilt=0.25, units='radians')
# xtal.getFriedelPair(qv[2][0][0, 0], qv[2][0][1, 0], qv[2][0][2, 0], chiTilt=0.25, units='radians')
# qv[1][0][:, [0, 3]], qv[2][0][:, [0, 3]]
# rMat = rot.rotMatOfExpMap(0.5*mgData.getTTh()[0]*mUtil.unitVector(num.c_[1,0,1].T))
# qv = mgData._PlaneData__makeScatteringVectors(rMat, chiTilt=0.)
# qv[1][0][:, 1], qv[2][0][:, 1]
