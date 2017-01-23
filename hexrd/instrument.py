# -*- coding: utf-8 -*-
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
"""
Created on Fri Dec  9 13:05:27 2016

@author: bernier2
"""

import numpy as np

from hexrd import gridutil as gutil
from hexrd import matrixutil as mutil
from hexrd.xrd.transforms_CAPI import anglesToGVec, \
                                      detectorXYToGvec, \
                                      gvecToDetectorXY, \
                                      makeDetectorRotMat, \
                                      mapAngle
from hexrd.xrd import xrdutil
from hexrd import constants as cnst

from hexrd.xrd.distortion import GE_41RT # BAD, VERY BAD!!! FIX!!!

beam_energy_DFLT = 65.351
beam_vec_DFLT = cnst.beam_vec

eta_vec_DFLT = cnst.eta_vec

panel_id_DFLT = "generic"
nrows_DFLT = 2048
ncols_DFLT = 2048
pixel_size_DFLT = (0.2, 0.2)

tilt_angles_DFLT = np.zeros(3)
t_vec_d_DFLT = np.r_[0., 0., -1000.]

chi_DFLT = 0.
t_vec_s_DFLT = np.zeros(3)

def calc_beam_vec(azim, pola):
    """
    Calculate unit beam propagation vector from 
    spherical coordinate spec in DEGREES

    ...MAY CHANGE; THIS IS ALSO LOCATED IN XRDUTIL!
    """
    tht = np.radians(azim)
    phi = np.radians(pola)
    bv = np.r_[
        np.sin(phi)*np.cos(tht),
        np.cos(phi),
        np.sin(phi)*np.sin(tht)]
    return -bv

class HEDMInstrument(object):
    """
    * Distortion needs to be moved to a class with registry; tuple unworkable
    * where should reference eta be defined? currently set to default config
    """
    def __init__(self, instrument_config=None,
                 image_series=None,
                 instrument_name="instrument",
                 ):
        self._id = instrument_name

        if instrument_config is None:
            self._num_panels = 1
            self._beam_energy = beam_energy_DFLT
            self._beam_vector = beam_vec_DFLT
            
            self._eta_vector = eta_vec_DFLT
            
            self._detectors = {
                panel_id_DFLT:PlanarDetector(
                    rows=nrows_DFLT, cols=ncols_DFLT,
                    pixel_size=pixel_size_DFLT,
                    tvec=t_vec_d_DFLT,
                    tilt=tilt_angles_DFLT,
                    bvec=self._beam_vector,
                    evec=self._eta_vector,
                    distortion=None),
                }
                        
            self._t_vec_s = t_vec_s_DFLT
            self._chi = chi_DFLT
        else:
            self._num_panels = len(instrument_config['detectors'])
            self._beam_energy = instrument_config['beam']['energy'] # keV
            self._beam_vector = calc_beam_vec(
                instrument_config['beam']['vector']['azimuth'],
                instrument_config['beam']['vector']['polar_angle'],
                )
            cnst.eta_vec
            # now build detector dict
            detector_ids = instrument_config['detectors'].keys()
            pixel_info = [instrument_config['detectors'][i]['pixels'] for i in detector_ids]
            affine_info = [instrument_config['detectors'][i]['transform'] for i in detector_ids]        
            distortion =  []
            for i in detector_ids:
                try:
                    distortion.append(
                        instrument_config['detectors'][i]['distortion']
                        )
                except KeyError:
                    distortion.append(None)
            det_list = []
            for pix, xform, dist in zip(pixel_info, affine_info, distortion):
                # HARD CODED GE DISTORTION !!! FIX
                dist_tuple = None
                if dist is not None: dist_tuple = (GE_41RT, dist['parameters'])
                
                det_list.append(
                    PlanarDetector(
                        rows=pix['rows'], cols=pix['columns'],
                        pixel_size=pix['size'],
                        tvec=xform['t_vec_d'],
                        tilt=xform['tilt_angles'],
                        bvec=self._beam_vector,
                        evec=cnst.eta_vec,
                        distortion=dist_tuple)
                    )
                pass
            self._detectors = dict(zip(detector_ids, det_list))

            self._t_vec_s = instrument_config['oscillation_stage']['t_vec_s']
            self._chi = instrument_config['oscillation_stage']['chi']            
            
        return
    
    # properties for physical size of rectangular detector
    @property
    def id(self):
        return self._id
    @property
    def num_panels(self):
        return self._num_panels
    @property
    def detectors(self):
        return self._detectors

    @property
    def tvec(self):
        return self._t_vec_s
    @tvec.setter
    def tvec(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3, 'input must have length = 3'
        self._t_vec_s = x

    @property
    def chi(self):
        return self._chi
    @chi.setter
    def chi(self, x):
        self._chi = float(x)

    @property
    def beam_energy(self):
        return self._beam_energy
    @beam_energy.setter
    def beam_energy(self, x):
        self._beam_energy = float(x)
    
    @property
    def beam_vector(self):
        return self._beam_vector
    @beam_vector.setter
    def beam_vector(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3 and sum(x*x) > 1-cnst.sqrt_epsf, \
            'input must have length = 3 and have unit magnitude'
        self._beam_vector = x
        # ...maybe change dictionary item behavior for 3.x compatibility?
        for detector_id in self._detectors:
            panel = self._detectors[detector_id]
            panel.bvec = self._beam_vector

    @property
    def eta_vector(self):
        return self._eta_vector
    @eta_vector.setter
    def eta_vector(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3 and sum(x*x) > 1-cnst.sqrt_epsf, \
            'input must have length = 3 and have unit magnitude'
        self._eta_vector = x
        # ...maybe change dictionary item behavior for 3.x compatibility?
        for detector_id in self._detectors:
            panel = self._detectors[detector_id]
            panel.evec = self._eta_vector

    # methods
    pass # end class: HEDMInstrument

class PlanarDetector(object):
    """
    base class for 2D row-column detector
    """

    __pixelPitchUnit = 'mm'
    __delta_eta = np.radians(10.)

    def __init__(self,
                 rows=2048, cols=2048,
                 pixel_size=(0.2, 0.2),
                 tvec=np.r_[0., 0., -1000.],
                 tilt=cnst.zeros_3,
                 bvec=cnst.beam_vec,
                 evec=cnst.eta_vec,
                 distortion=None):
        """
        """
        self._rows = rows
        self._cols = cols

        self.pixel_size_row = pixel_size[0]
        self.pixel_size_col = pixel_size[1]

        self._tvec = np.array(tvec).flatten()
        self._tilt = tilt

        self._bvec = np.array(bvec).flatten()
        self._evec = np.array(evec).flatten()

        self._distortion = distortion
        return

    # properties for physical size of rectangular detector
    @property
    def rows(self):
        return self._rows
    @rows.setter
    def rows(self, x):
        assert isinstance(x, int)
        self._rows = x
    @property
    def cols(self):
        return self._cols
    @cols.setter
    def cols(self, x):
        assert isinstance(x, int)
        self._cols = x

    @property
    def row_dim(self):
        return self.rows * self.pixel_size_row
    @property
    def col_dim(self):
        return self.cols * self.pixel_size_col

    @property
    def row_pixel_vec(self):
        return self.pixel_size_row*(0.5*(self.rows-1)-np.arange(self.rows))
    @property
    def row_edge_vec(self):
        return self.pixel_size_row*(0.5*self.rows-np.arange(self.rows+1))
    @property
    def col_pixel_vec(self):
        return self.pixel_size_col*(np.arange(self.cols)-0.5*(self.cols-1))
    @property
    def col_edge_vec(self):
        return self.pixel_size_col*(np.arange(self.cols+1)-0.5*self.cols)

    @property
    def corner_ul(self):
        return np.r_[-0.5 * self.col_dim,  0.5 * self.row_dim]
    @property
    def corner_ll(self):
        return np.r_[-0.5 * self.col_dim, -0.5 * self.row_dim]
    @property
    def corner_lr(self):
        return np.r_[ 0.5 * self.col_dim, -0.5 * self.row_dim]
    @property
    def corner_ur(self):
        return np.r_[ 0.5 * self.col_dim,  0.5 * self.row_dim]

    @property
    def tvec(self):
        return self._tvec
    @tvec.setter
    def tvec(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3, 'input must have length = 3'
        self._tvec = x

    @property
    def tilt(self):
        return self._tilt
    @tilt.setter
    def tilt(self, x):
        assert len(x) == 3, 'input must have length = 3'
        self._tilt = np.array(x).squeeze()

    @property
    def bvec(self):
        return self._bvec
    @bvec.setter
    def bvec(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3 and sum(x*x) > 1-cnst.sqrt_epsf, \
            'input must have length = 3 and have unit magnitude'
        self._bvec = x

    @property
    def evec(self):
        return self._evec
    @evec.setter
    def evec(self, x):
        x = np.array(x).flatten()
        assert len(x) == 3 and sum(x*x) > 1-cnst.sqrt_epsf, \
            'input must have length = 3 and have unit magnitude'
        self._evec = x

    @property
    def distortion(self):
        return self._distortion
    @distortion.setter
    def distortion(self, x):
        """
        Probably should make distortion a class...
        """
        assert len(x) == 2 and hasattr(x[0], '__call__'), \
            'distortion must be a tuple: (<func>, params)'
        self._distortion = x

    @property
    def rmat(self):
        return makeDetectorRotMat(self.tilt)

    @property
    def normal(self):
        return self.rmat[:, 2]

    @property
    def beam_position(self):
        """
        returns the coordinates of the beam in the cartesian detector
        frame {Xd, Yd, Zd}.  NaNs if no intersection.
        """
        output = np.nan * np.ones(2)
        b_dot_n = np.dot(self.bvec, self.normal)
        if np.logical_and(
            abs(b_dot_n) > cnst.sqrt_epsf,
            np.sign(b_dot_n) == -1
            ):
            u = np.dot(self.normal, self.tvec) / b_dot_n
            p2_l = u*self.bvec
            p2_d = np.dot(self.rmat.T, p2_l - self.tvec)
            output = p2_d[:2]
        return output

    @property
    def pixel_coords(self):
        pix_i, pix_j = np.meshgrid(
            self.row_pixel_vec, self.col_pixel_vec,
            indexing='ij')
        return pix_i, pix_j

    @property
    def pixel_angles(self):
        pix_i, pix_j = self.pixel_coords
        xy = np.ascontiguousarray(
            np.vstack([
                pix_j.flatten(), pix_i.flatten()
                ]).T
            )
        angs, g_vec = detectorXYToGvec(
            xy, self.rmat, cnst.identity_3x3,
            self.tvec, cnst.zeros_3, cnst.zeros_3,
            beamVec=self.bvec, etaVec=self.evec)
        del(g_vec)
        tth = angs[0].reshape(self.rows, self.cols)
        eta = angs[1].reshape(self.rows, self.cols)
        return tth, eta


    """
    ##################### METHODS
    """
    def cartToPixel(self, xy_det, pixels=False):
        """
        Convert vstacked array or list of [x,y] points in the center-based
        cartesian frame {Xd, Yd, Zd} to (i, j) edge-based indices

        i is the row index, measured from the upper-left corner
        j is the col index, measured from the upper-left corner

        if pixels=True, then (i,j) are integer pixel indices.
        else (i,j) are continuous coords
        """
        xy_det = np.atleast_2d(xy_det)

        npts = len(xy_det)

        tmp_ji = xy_det - np.tile(self.corner_ul, (npts, 1))
        i_pix = -tmp_ji[:, 1] / self.pixel_size_row - 0.5
        j_pix =  tmp_ji[:, 0] / self.pixel_size_col - 0.5

        ij_det = np.vstack([i_pix, j_pix]).T
        if pixels:
            ij_det = np.array(np.round(ij_det), dtype=int)
        return ij_det

    def pixelToCart(self, ij_det):
        """
        Convert vstacked array or list of [i,j] pixel indices
        (or UL corner-based points) and convert to (x,y) in the
        cartesian frame {Xd, Yd, Zd}
        """
        ij_det = np.atleast_2d(ij_det)

        x = (ij_det[:, 1] + 0.5)*self.pixel_size_col + self.corner_ll[0]
        y = (self.rows - ij_det[:, 0] - 0.5)*self.pixel_size_row + self.corner_ll[1]
        return np.vstack([x, y]).T

    def angularPixelSize(self, xy, rMat_s=None, tVec_s=None, tVec_c=None):
        """
        Wraps xrdutil.angularPixelSize
        """
        # munge kwargs
        if rMat_s is None: rMat_s = cnst.identity_3x3
        if tVec_s is None: tVec_s = cnst.zeros_3x1
        if tVec_c is None: tVec_c = cnst.zeros_3x1

        # call function
        ang_ps = xrdutil.angularPixelSize(
            xy, (self.pixel_size_row, self.pixel_size_col),
            self.rmat, rMat_s,
            self.tvec, tVec_s, tVec_c,
            distortion=self.distortion, 
            beamVec=self.bvec, etaVec=self.evec)
        return ang_ps

    def clip_to_panel(self, xy, buffer_edges=False):
        """
        """
        xy = np.atleast_2d(xy)
        xlim = 0.5*self.col_dim
        if buffer_edges:
            xlim -= 0.5*self.pixel_size_col
        ylim = 0.5*self.row_dim
        if buffer_edges:
            ylim -= 0.5*self.pixel_size_row
        on_panel_x = np.logical_and(xy[:, 0] >= -xlim, xy[:, 0] <= xlim)
        on_panel_y = np.logical_and(xy[:, 1] >= -ylim, xy[:, 1] <= ylim)
        on_panel = np.where(np.logical_and(on_panel_x, on_panel_y))[0]
        return xy[on_panel, :], on_panel

    def interpolate_bilinear(self, xy, img, pad_with_nans=True):
        """
        """
        is_2d = img.ndim == 2
        right_shape = img.shape[0] == self.rows and img.shape[1] == self.cols
        assert is_2d and right_shape, \
          "input image must be 2-d with shape (%d, %d)" %(self.rows, self.cols)

        # initialize output with nans
        if pad_with_nans:
            int_xy = np.nan*np.ones(len(xy))
        else:
            int_xy = np.zeros(len(xy))
        
        # clip away points too close to or off the edges of the detector
        xy_clip, on_panel = self.clip_to_panel(xy, buffer_edges=True)

        # grab fractional pixel indices of clipped points
        ij_frac = self.cartToPixel(xy_clip)

        # get floors/ceils from array of pixel _centers_
        i_floor = gutil.cellIndices(self.row_pixel_vec, xy_clip[:, 1])
        j_floor = gutil.cellIndices(self.col_pixel_vec, xy_clip[:, 0])
        i_ceil = i_floor + 1
        j_ceil = j_floor + 1

        # first interpolate at top/bottom rows
        row_floor_int = \
            (j_ceil - ij_frac[:, 1])*img[i_floor, j_floor] \
            + (ij_frac[:, 1] - j_floor)*img[i_floor, j_ceil]
        row_ceil_int = \
            (j_ceil - ij_frac[:, 1])*img[i_ceil, j_floor] \
            + (ij_frac[:, 1] - j_floor)*img[i_ceil, j_ceil]

        # next interpolate across cols
        int_vals = \
            (i_ceil - ij_frac[:, 0])*row_floor_int \
            + (ij_frac[:, 0] - i_floor)*row_ceil_int
        int_xy[on_panel] = int_vals
        return int_xy

    def make_powder_rings(
        self, pd, merge_hkls=False, delta_eta=None, eta_period=None,
        rmat_s=cnst.identity_3x3,  tvec_s=cnst.zeros_3,
        tvec_c=cnst.zeros_3
        ):
        """
        """
        
        # for generating rings
        if delta_eta is None: delta_eta=self.__delta_eta
        if eta_period is None: eta_period = (-np.pi, np.pi)
        
        neta = int(360./float(delta_eta))
        eta = mapAngle(
            np.radians(delta_eta*np.linspace(0, neta-1, num=neta)) + eta_period[0], 
            eta_period
        )
        
        if merge_hkls:
            tth_idx, tth_ranges = pd.getMergedRanges()
            tth = [0.5*sum(i) for i in tth_ranges]
        else:
            tth = pd.getTTh()
        angs = [np.vstack([i*np.ones(neta), eta, np.zeros(neta)]) for i in tth]

        # need xy coords and pixel sizes
        valid_ang = []
        valid_xy = []
        for i_ring in range(len(angs)):
            these_angs = angs[i_ring].T
            gVec_ring_l = anglesToGVec(these_angs, bHat_l=self.bvec)
            xydet_ring = gvecToDetectorXY(
                gVec_ring_l,
                self.rmat, rmat_s, cnst.identity_3x3,
                self.tvec, tvec_s, tvec_c)
            #
            xydet_ring, on_panel = self.clip_to_panel(xydet_ring)
            #
            valid_ang.append(these_angs[on_panel, :2])
            valid_xy.append(xydet_ring)
            pass
        return valid_ang, valid_xy

    def map_to_plane(self, pts, rmat, tvec):
        """
        map detctor points to specified plane

        by convention      
        
        n * (u*pts_l - tvec) = 0
        
        [pts]_l = rmat*[pts]_m + tvec
        """
        
        # arg munging
        pts = np.atleast_2d(pts); npts = len(pts)

        # map plane normal & translation vector, LAB FRAME
        nvec_map_lab = rmat[:, 2].reshape(3, 1) 
        tvec_map_lab = np.atleast_2d(tvec).reshape(3, 1)
        tvec_d_lab = np.atleast_2d(self.tvec).reshape(3, 1)
        
        # put pts as 3-d in panel CS and transform to 3-d lab coords
        pts_det = np.hstack([pts, np.zeros((npts, 1))])
        pts_lab = np.dot(self.rmat, pts_det.T) + tvec_d_lab
        
        # scaling along pts vectors to hit map plane
        u = np.dot(nvec_map_lab.T, tvec_map_lab) \
            / np.dot(nvec_map_lab.T, pts_lab)
        
        # pts on map plane, in LAB FRAME
        pts_map_lab = np.tile(u, (3, 1)) * pts_lab
        
        return np.dot(rmat.T, pts_map_lab - tvec_map_lab)[:2, :].T
        
        
        
