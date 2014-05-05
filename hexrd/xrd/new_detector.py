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

import os, sys

import numpy as np

from hexrd import matrixutil as mutil

from hexrd.paramters import *

from hexrd.xrd import fitting
from hexrd.xrd import distortion      as dfunc
from hexrd.xrd import rotations       as rot
from hexrd.xrd import transforms_CAPI as xf

from hexrd.xrd.crystallography import processWavelength

class planarDetector():
    """
    base class for 2D row-column flat panel detector
    """
    __pixelPitchUnit = 'mm'
    
    def __init__(self, 
                 row_col_dims, pixelPitch_rc, 
                 tilt=tilt_DFLT, tVec=tVec_d_DFLT,
                 bVec=bVec_DFLT, eVec=eVec_DFLT,
                 rMat_s=rMat_s_DFLT, tVec_s=tVec_s_DFLT,
                 tVec_c=tVec_c_DFLT,
                 distortion=None):
        
        self.__nrows = row_col_dims[0]
        self.__ncols = row_col_dims[1]
        
        self.__pixelPitch_r = pixelPitch_rc[0]
        self.__pixelPitch_c = pixelPitch_rc[1]
        
        self.__rdim = nrows * pixelPitch_rc[0]
        self.__cdim = ncols * pixelPitch_rc[1]
        
        self.__tVec_LL = np.c_[-0.5 * self.__cdim, -0.5 * self.__rdim, 0.].T 
        self.__tVec_UL = np.c_[-0.5 * self.__cdim,  0.5 * self.__rdim, 0.].T 
        self.__tVec_UR = np.c_[ 0.5 * self.__cdim,  0.5 * self.__rdim, 0.].T 
        self.__tVec_LR = np.c_[ 0.5 * self.__cdim, -0.5 * self.__rdim, 0.].T

        self.__i_cen = np.arange(nrows)
        self.__j_cen = np.arange(ncols)

        self.__x_cen = self.__pixelPitch_c * self.__j_cen + 0.5 * (self.__pixelPitch_c - self.__cdim)
        self.__y_cen = self.__pixelPitch_r * self.__i_cen + 0.5 * (self.__pixelPitch_r - self.__rdim)
                
        # set:
        # self.__tilt, self.__rMat, self.__nVec
        self.set_tilt(tilt)

        self.__bVec = bVec
        self.__eVec = eVec
        self.__tVec = tVec

        self.rMat_s = rMat_s
        self.tVec_s = tVec_s
        self.tVec_c = tVec_c
        
        self.distortion = distortion

        return
    
    # properties for physical size of rectangular detector
    def set_nrows(self, nrows):
        raise RuntimeError, 'set of nrows not allowed'
    def get_nrows(self):
        return self.__nrows
    nrows = property(get_nrows, set_nrows, None)

    def set_ncols(self, ncols):
        raise RuntimeError, 'set of ncols not allowed'
    def get_ncols(self):
        return self.__ncols
    ncols = property(get_ncols, set_ncols, None)

    def set_pixelPitch_r(self, pixelPitch_r):
        raise RuntimeError, 'set of pixelPitch not allowed'
    def get_pixelPitch_r(self):
        return self.__pixelPitch_r
    pixelPitch_r = property(get_pixelPitch_r, set_pixelPitch_r, None)

    def set_pixelPitch_c(self, pixelPitch_c):
        raise RuntimeError, 'set of pixelPitch not allowed'
    def get_pixelPitch_c(self):
        return self.__pixelPitch_c
    pixelPitch_c = property(get_pixelPitch_c, set_pixelPitch_c, None)

    def set_rdim(self, rdim):
        raise RuntimeError, 'set of nrows not allowed'
    def get_rdim(self):
        return self.__rdim
    rdim = property(get_rdim, set_rdim, None)

    def set_cdim(self, cdim):
        raise RuntimeError, 'set of nrows not allowed'
    def get_cdim(self):
        return self.__cdim
    cdim = property(get_cdim, set_cdim, None)
    
    # want tvec to be a vector
    def get_tVec(self):
        return self.__tVec
    def set_tVec(self, tVec):
        assert len(tVec) == 3, 'tvec must be a list or array with length 3'
        self.__tVec = np.array(tVec).reshape(3, 1)
        return
    tVec = property(get_tVec, set_tVec, None)
    
    # tx
    def get_tx(self):
        return self.tVec[0]
    def set_tx(self, tx):
        self.__tVec[0] = tx
        return
    tx = property(get_tx, set_tx, None)
    
    # ty
    def get_ty(self):
        return self.tVec[1]
    def set_ty(self, ty):
        self.__tVec[1] = ty
        return
    ty = property(get_ty, set_ty, None)
    
    # tz
    def get_tz(self):
        return self.tVec[2]
    def set_tz(self, tz):
        self.__tVec[2] = tz
        return
    tz = property(get_tz, set_tz, None)
    
    def get_rMat(self):
        return self.__rMat
    def set_rMat(self):
        raise RuntimeError, 'set of rMat is not allowed; computed from tilt'
        return
    rMat = property(get_rMat, set_rMat, None)
    
    def get_normal(self):
        return self.__nVec
    def set_normal(self):
        raise RuntimeError, 'set of normal is not allowed; computed from tilt'
        return
    nVec = property(get_normal, set_normal, None)
    
    # setting tilt updates detector normal and COB matrix
    def get_tilt(self):
        return self.__tilt
    def set_tilt(self, tilt):
        """
        default geometry defined in transforms...
        detector normal is local Z
        """
        self.__rMat = xf.makeDetectorRotMat(tilt)
        self.__nVec = np.dot(self.__rMat, Z_ref)
        self.__tilt = tilt
        return
    tilt = property(get_tilt, set_tilt, None)
    
    # x tilt
    def get_xTilt(self):
        return self.tilt[0]
    def set_xTilt(self, xTilt):
        self.set_tilt([xTilt, self.tilt[1], self.tilt[2]])
        return
    xTilt = property(get_xTilt, set_xTilt, None)
    
    # y tilt
    def get_yTilt(self):
        return self.tilt[0]
    def set_yTilt(self, yTilt):
        self.set_tilt([self.tilt[0], yTilt, self.tilt[2]])
        return
    yTilt = property(get_yTilt, set_yTilt, None)
    
    # z tilt
    def get_zTilt(self):
        return self.tilt[0]
    def set_zTilt(self, zTilt):
        self.set_tilt([self.tilt[0], self.tilt[1], zTilt])
        return
    zTilt = property(get_zTilt, set_zTilt, None)

    # have to carry the beam vector around
    def get_bVec(self):
        return self.__bVec
    def set_bVec(self, bVec):
        assert len(bVec) == 3, 'beam vector must be a list or array with length 3'
        self.__bVec = np.array(bVec).reshape(3, 1)
        return
    bVec = property(get_bVec, set_bVec, None)

    # have to cary the eta ref around too
    def get_eVec(self):
        return self.__eVec
    def set_eVec(self, eVec):
        assert len(eVec) == 3, 'eta ref vector must be a list or array with length 3'
        self.__eVec = np.array(eVec).reshape(3, 1)
        return
    eVec = property(get_eVec, set_eVec, None)
    
    # beam position (where applicable)
    def get_beamXY(self):
        """
        returns the in-plane coordinates of the beam in the cartesian detector 
        frame {Xd, Yd, Zd}.  NaNs if no intersection.
        """
        beamPos = np.nan * np.r_[1, 1]
        
        # point coordinate arrays IN LAB FRAME
        #     *) P0 = [0, 0, 0] not needed
        b  = self.bVec
        P1 = self.tVec
        
        bDotN = np.dot(b.T, self.nVec).flatten()
        if np.logical_and(abs(bDotN) > sqrt_epsf, np.sign(bDotN) == -1):
            u = np.dot(self.nVec.T, P1) / bDotN
            
            P2_l = u*b
            P2_d = np.dot(self.rMat.T, P2_l - self.tVec)
            beamPos = P2_d[:2].flatten()
        return beamPos
    def set_beamXY(self):
        raise RuntimeError, "cannot set beam position; calculated quantity"
        return
    beamXY = property(get_beamXY, set_beamXY, None)
    
    def get_pixel_xy(self):
        return np.meshgrid(self.__x_cen, self.__y_cen)
    def set_pixel_xy(self):
        raise RuntimeError, "cannot set pixel xy coords"
        return
    pixelXY = property(get_pixel_xy, set_pixel_xy, None)

    def get_pixel_angles(self):
        gvecs = self.angs_of_xy(self.pixelXY)
        return (tth, eta)
    def set_pixel_angles(self):
        raise RuntimeError, "cannot set pixel angles; calculated quantity"
        return
    pixelAngles = property(get_pixel_angles, set_pixel_angles, None)
    
    def getAngPixelSize(self, xy_det):
        """
        get pixel size in angular coordinates
        at a given cartesian coordinate position

        xy_det.shape = (2, n)
        """
        ij_pix = self.ij_of_xy(xy_det)

        return unc

    def ij_of_xy(self, xy_det, pixels=False):
        """
        Convert vstacked array or list of [x,y] points in the center-based 
        cartesian frame {Xd, Yd, Zd} to (i, j) edge-based indices
        
        i is the row index, measured from the upper-left corner
        j is the col index, measured from the upper-left corner
        
        if pixels=True, then (i,j) are integer pixel indices.  
        else (i,j) are continuous coords
        """
        xy_det = np.atleast_2d(xy_det)
        
        npts   = len(xy_det)
        
        tmpJI = xy_det - np.tile(self.__tVec_UL[:2].flatten(), (npts, 1))
        I = -tmpJI[:, 1] / self.pixelPitch_r
        J =  tmpJI[:, 0] / self.pixelPitch_c
        
        ij_det = np.vstack([I, J]).T        
        if pixels:
            ij_det = np.floor(ij_det)
        return ij_det
    
    def xy_of_ij(self, ij_det):
        """
        Convert vstacked array or list of [i,j] pixel indices
        (or UL corner-based points) and convert to (x,y) in the 
        cartesian frame {Xd, Yd, Zd}

        cartesianCoordsOfPixelIndices
        """
        ij_det = np.atleast_2d(ij_det)
        
        npts   = len(ij_det)
        
        X =  ij_det[:, 1] * self.pixelPitch_c + self.__tVec_UL[0]
        Y = -ij_det[:, 0] * self.pixelPitch_r + self.__tVec_UL[1]
        return np.vstack([X, Y]).T
    
    def angs_of_xy(self, xy_det, *args):
        """
        vstacked (tth, eta) pairs of cartesian coordinates
        *) wrapper for transforms.detectorXYtoGvec
        
        args: ome, chi
        """
        xy_det = np.atleast_2d(xy_det) # so len() gives right number even if one point
        npts   = len(xy_det)
        if len(args) > 0:
            ome = np.atleast_1d(args[0]).flatten()
            assert len(ome) == npts, "ome must be the same length as xy_det!"
            chi = args[1]
            angles = np.zeros((npts, 2))
            gVecs  = np.zeros((npts, 3))
            for i in range(npts):
                rMat_s = xf.makeOscillRotMat([chi, ome[i]]) 
                tmp = xf.detectorXYToGvec(xy_det,
                                          self.rMat, rMat_s,
                                          self.tVec, self.tVec_s, self.tVec_c,
                                          distortion=self.distortion,
                                          beamVec=self.bVec,
                                          etaVec=self.eVec)
                angles[i, :] = tmp[0]
                gVecs[i, :]  = tmp[1]
                pass
        else:
            angles, gVecs = xf.detectorXYToGvec(xy_det,
                                                self.rMat, self.rMat_s,
                                                self.tVec, self.tVec_s, self.tVec_c,
                                                distortion=self.distortion,
                                                beamVec=self.bVec,
                                                etaVec=self.eVec)
            pass
        
        # filter out NaNs
        not_there = np.isnan(angles[:, 0])
        tTh = angles[-not_there, 0]
        eta = angles[-not_there, 1]
        
        return tTh, eta
    
    def xy_of_angs(self, angs):
        """
        Cartesion coordinates of vstacked (tth, eta) pairs
        *) wrapper for transforms.anglesToGVec
        """
        gVec_l = xf.anglesToGVec(angs, self.bVec, self.eVec,
                                 rMat_s=None, rMat_c=None)
        xy = xf.gvecToDetectorXY(gVec_l,
                                 self.rMat, I3, I3, 
                                 self.tVec, self.tVec_s, self.tVec_c,
                                 distortion=self.distortion,
                                 beamVec=self.bVec, 
                                 etaVec=self.eVec)
        return xy
        
    def simulateLauePattern(self, planeData, minEnergy=5, maxEnergy=25, rMat_s=np.eye(3), rMat=None, vInv=None, doGnomonic=False):
        
        multipleEnergyRanges = False
        if hasattr(maxEnergy, '__len__'):
            assert len(maxEnergy) == len(minEnergy), 'energy cutoff ranges must have the same length'
            multipleEnergyRanges = True; lmin = []; lmax = []
            for i in range(len(maxEnergy)):
                lmin.append(processWavelength(maxEnergy[i]))
                lmax.append(processWavelength(minEnergy[i]))                
        else:
            lmin = processWavelength(maxEnergy)
            lmax = processWavelength(minEnergy)
        
        gvec_c = planeData.getPlaneNormals()
        hkls   = planeData.getSymHKLs() 
        dsp    = planeData.getPlaneSpacings()
        
        if rMat is None:
            rMat = [np.eye(3),]
        if vInv is None:
            vInv = [vInv_ref,]
        
        # rMult     = planeData.getMultiplicity()
        # nHKLs_tot = rMult.sum()
        # rMask     = np.ones(nHKLs_tot, dtype=bool)
        
        retval = []
        for iG in range(len(rMat)):
            tmp = {'detXY':[], 'gnoXY':[], 'angles':[], 'dspacing':[], 'hkl':[], 'energy':[]}
            for iHKL in range(planeData.nHKLs):                
                # stretch them: V^(-1) * R * Gc
                gvec_s_str = mutil.unitVector(
                    np.dot( vInv[iG], np.dot( rMat[iG], gvec_c[iHKL] ) ) )
                gvec_c_str = np.dot(rMat[iG].T, gvec_s_str)
                gvec_l_str = np.dot(rMat_s, gvec_s_str)
                #
                # dpts  = self.gVecToDet(gvec_c_str, rMat=rMat[iG], rMat_s=rMat_s)
                # gpts  = self.gVecToDet(gvec_c_str, rMat=rMat[iG], rMat_s=rMat_s, doGnomonic=True)
                dpts  = self.gVecToDet(gvec_c_str, rMat=rMat[iG], rMat_s=rMat_s)
                gpts  = self.gVecToDet(gvec_c_str, rMat=rMat[iG], rMat_s=rMat_s, doGnomonic=True)
                canIntersect = -np.isnan(dpts[0, :])
                npts_in = sum(canIntersect)
                if np.any(canIntersect):
                    dpts = dpts[:, canIntersect].reshape(3, npts_in)
                    dhkl = hkls[iHKL][:, canIntersect].reshape(3, npts_in)

                    gvl_hat = gvec_l_str[:, canIntersect].reshape(3, npts_in)
                    gvl_xy  = gvec_l_str[:2, canIntersect].reshape(2, npts_in)
                    
                    # dot with the beam
                    dotWbeam = np.dot(Z_ref.T, gvl_hat).flatten()

                    # angles
                    theta = piby2 - rot.arccosSafe( dotWbeam )
                    wlen  = 2*dsp[iHKL]*np.sin(theta)
                    eta   = np.arccos(gvl_xy[0, :])
                    
                    # find on spatial extent of detector
                    # for corner # xTest = np.logical_and(dpts[0, :] > 0, dpts[0, :] < self.cdim)
                    # for corner # yTest = np.logical_and(dpts[1, :] > 0, dpts[1, :] < self.rdim)
                    xTest = np.logical_and(dpts[0, :] > -0.5 * self.cdim, dpts[0, :] < 0.5 * self.cdim)
                    yTest = np.logical_and(dpts[1, :] > -0.5 * self.rdim, dpts[1, :] < 0.5 * self.rdim)
                    
                    onDetector  = np.logical_and(xTest, yTest)
                    if multipleEnergyRanges:
                        validEnergy = np.zeros(len(wlen), dtype=bool)
                        for i in range(len(lmin)):
                            validEnergy = validEnergy | np.logical_and(wlen >= lmin[i], wlen <= lmax[i])
                            pass
                    else:
                        validEnergy = np.logical_and(wlen >= lmin, wlen <= lmax)
                        pass
                    
                    keepers = np.logical_and(onDetector, validEnergy)

                    dsp_this = 1. / mutil.columnNorm(np.dot(planeData.latVecOps['B'], dhkl[:, keepers]))
                    
                    tmp['detXY'].append(dpts[:2, keepers])
                    tmp['gnoXY'].append(gpts[:2, keepers])
                    tmp['angles'].append(np.vstack([2*theta[keepers]*r2d, eta[keepers]*r2d]))
                    tmp['hkl'].append(dhkl[:, keepers])
                    tmp['dspacing'].append(dsp_this)
                    tmp['energy'].append(processWavelength(wlen[keepers]))
                else:
                    tmp['detXY'].append(np.empty((2, 0)))
                    tmp['gnoXY'].append(np.empty((2, 0)))
                    tmp['angles'].append(np.empty((2, 0)))
                    tmp['hkl'].append(np.empty((3, 0)))
                    tmp['dspacing'].append(np.empty(0))
                    tmp['energy'].append(np.empty(0))
                    pass
                pass
            retval.append(tmp)
        return retval
    

