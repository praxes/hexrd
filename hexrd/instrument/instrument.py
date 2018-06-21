# -*- coding: utf-8 -*-
# =============================================================================
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
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
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
# =============================================================================
"""
Created on Fri Dec  9 13:05:27 2016

@author: bernier2
"""
from __future__ import print_function

import os

import yaml

import h5py

import numpy as np

from scipy import ndimage
from scipy.linalg.matfuncs import logm

from hexrd.gridutil import make_tolerance_grid
from hexrd.xrd.transforms_CAPI import anglesToGVec, \
                                      detectorXYToGvec, \
                                      gvecToDetectorXY, \
                                      makeDetectorRotMat, \
                                      makeOscillRotMat, \
                                      makeRotMatOfExpMap, \
                                      mapAngle, \
                                      oscillAnglesOfHKLs, \
                                      rowNorm, \
                                      validateAngleRanges
from hexrd.xrd import xrdutil
from hexrd import constants as ct

# FIXME: distortion kludge
from hexrd.xrd.distortion import GE_41RT  # BAD, VERY BAD!!!


from skimage.draw import polygon

# =============================================================================
# PARAMETERS
# =============================================================================

# [wavelength, chi, tvec_s, expmap_c, tec_c], len is 11
instr_param_flags_DFLT = np.array(
    [0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1],
    dtype=bool)
panel_param_flags_DFLT = np.array(
    [1, 1, 1, 1, 1, 1],
    dtype=bool)

# =============================================================================
# UTILITY METHODS
# =============================================================================


def migrate_instrument_config(instrument_config):
    """utility function to generate old instrument config dictionary"""
    cfg_list = []
    for detector_id in instrument_config['detectors']:
        cfg_list.append(
            dict(
                detector=instrument_config['detectors'][detector_id],
                oscillation_stage=instrument_config['oscillation_stage'],
            )
        )
    return cfg_list


def angle_in_range(angle, ranges, ccw=True, units='degrees'):
    """
    Return the index of the first wedge the angle is found in

    WARNING: always clockwise; assumes wedges are not overlapping
    """
    tau = 360.
    if units.lower() == 'radians':
        tau = 2*np.pi
    w = np.nan
    for i, wedge in enumerate(ranges):
        amin = wedge[0]
        amax = wedge[1]
        check = amin + np.mod(angle - amin, tau)
        if check < amax:
            w = i
            break
    return w


# ???: move to gridutil?
def centers_of_edge_vec(edges):
    assert np.r_[edges].ndim == 1, "edges must be 1-d"
    return np.average(np.vstack([edges[:-1], edges[1:]]), axis=0)


# =============================================================================
# CLASSES
# =============================================================================


class HEDMInstrument(object):
    """
    * Distortion needs to be moved to a class with registry; tuple unworkable
    * where should reference eta be defined? currently set to default config
    """
    def __init__(self, beam, detectors, oscillation_stage,
                 image_series=None, eta_vector=None,
                 name="instrument"):
        self._id = name
        self._beam = beam
        self._detectors = detectors
        self._oscillation_stage = oscillation_stage

        self._eta_vector = ct.eta_vec if eta_vector is None else eta_vector
        self._num_panels = len(detectors)
        self._tvec = oscillation_stage.tvec
        self._chi = oscillation_stage.chi

        self._param_flags = np.hstack(
            [instr_param_flags_DFLT,
             np.tile(panel_param_flags_DFLT, self._num_panels)]
        )
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
    def detector_parameters(self):
        pdict = {}
        for key, panel in self.detectors.iteritems():
            pdict[key] = panel.config_dict(self.chi, self.tvec)
        return pdict

    # Oscillation stage
    @property
    def oscillation_stage(self):
        return self._oscillation_stage

    @property
    def tvec(self):
        return self._tvec

    @property
    def chi(self):
        return self._chi

    # Beam properties
    @property
    def beam(self):
        return self._beam

    @property
    def beam_energy(self):
        return self.beam.energy

    @property
    def beam_wavelength(self):
        return self.beam.wavelength

    @property
    def beam_vector(self):
        return self.beam.vector

    # Miscellaneous parameters
    @property
    def eta_vector(self):
        return self._eta_vector

    @property
    def param_flags(self):
        return self._param_flags

    @param_flags.setter
    def param_flags(self, x):
        x = np.array(x, dtype=bool).flatten()
        assert len(x) == 11 + 6*self.num_panels, \
            "length of parameter list must be %d; you gave %d" \
            % (len(self._param_flags), len(x))
        self._param_flags = x

    # =========================================================================
    # METHODS
    # =========================================================================

    def calibration_params(self, expmap_c, tvec_c):
        plist = np.zeros(11 + 6*self.num_panels)

        plist[0] = self.beam_wavelength
        plist[1] = self.chi
        plist[2], plist[3], plist[4] = self.tvec
        plist[5], plist[6], plist[7] = expmap_c
        plist[8], plist[9], plist[10] = tvec_c

        ii = 11
        for panel in self.detectors.itervalues():
            plist[ii:ii + 6] = np.hstack([
                panel.tilt.flatten(),
                panel.tvec.flatten(),
            ])
            ii += 6

        # FIXME: FML!!!
        # this assumes old style distiortion = (func, params)
        retval = plist
        for panel in self.detectors.itervalues():
            if panel.distortion is not None:
                retval = np.hstack([retval, panel.distortion[1]])
        return retval

    def write_config(self, filename=None, calibration_dict={}):
        """ WRITE OUT YAML FILE """
        # initialize output dictionary

        par_dict = {}

        azim, pola = calc_angles_from_beam_vec(self.beam_vector)
        beam = dict(
            energy=self.beam_energy,
            vector=dict(
                azimuth=azim,
                polar_angle=pola,
            )
        )
        par_dict['beam'] = beam

        if calibration_dict:
            par_dict['calibration_crystal'] = calibration_dict

        ostage = dict(
            chi=self.chi,
            t_vec_s=self.tvec.tolist()
        )
        par_dict['oscillation_stage'] = ostage

        det_names = self.detectors.keys()
        det_dict = dict.fromkeys(det_names)
        for det_name in det_names:
            panel = self.detectors[det_name]
            pdict = panel.config_dict(self.chi, self.tvec)
            det_dict[det_name] = pdict['detector']
        par_dict['detectors'] = det_dict
        if filename is not None:
            with open(filename, 'w') as f:
                yaml.dump(par_dict, stream=f)
        return par_dict

    def extract_polar_maps(self, plane_data, imgser_dict,
                           active_hkls=None, threshold=None,
                           tth_tol=None, eta_tol=0.25):
        """
        Quick and dirty way to histogram angular patch data for make
        pole figures suitable for fiber generation

        TODO: streamline projection code
        TODO: normalization
        """
        if tth_tol is not None:
            plane_data.tThWidth = np.radians(tth_tol)
        else:
            tth_tol = np.degrees(plane_data.tThWidth)

        tth_ranges = plane_data.getTThRanges()
        if active_hkls is not None:
            assert hasattr(active_hkls, '__len__'), \
                "active_hkls must be an iterable with __len__"
            tth_ranges = tth_ranges[active_hkls]

        # need this for making eta ranges
        eta_tol_vec = 0.5*np.radians([-eta_tol, eta_tol])

        ring_maps_panel = dict.fromkeys(self.detectors)
        for i_d, det_key in enumerate(self.detectors):
            print("working on detector '%s'..." % det_key)

            # grab panel
            panel = self.detectors[det_key]
            # native_area = panel.pixel_area  # pixel ref area

            # make rings clipped to panel
            pow_angs, pow_xys, eta_idx, full_etas = panel.make_powder_rings(
                plane_data,
                merge_hkls=False, delta_eta=eta_tol,
                output_etas=True)

            ptth, peta = panel.pixel_angles
            ring_maps = []
            for i_r, tthr in enumerate(tth_ranges):
                print("working on ring %d..." % i_r)
                rtth_idx = np.where(
                    np.logical_and(ptth >= tthr[0], ptth <= tthr[1])
                )
                etas = pow_angs[i_r][:, 1]
                netas = len(etas)
                eta_ranges = np.tile(etas, (2, 1)).T \
                    + np.tile(eta_tol_vec, (netas, 1))
                ring_map = []
                for i_e, etar in enumerate(eta_ranges):
                    # WARNING: assuming start/stop
                    emin = np.r_[etar[0]]
                    emax = np.r_[etar[1]]
                    reta_idx = np.where(
                        validateAngleRanges(peta[rtth_idx], emin, emax)
                    )
                    ijs = (rtth_idx[0][reta_idx],
                           rtth_idx[1][reta_idx])
                    ring_map.append(ijs)
                    pass

                try:
                    omegas = imgser_dict[det_key].metadata['omega']
                except(KeyError):
                    msg = "imageseries for '%s' has no omega info" % det_key
                    raise RuntimeError(msg)
                nrows_ome = len(omegas)
                ncols_eta = len(full_etas)
                this_map = np.nan*np.ones((nrows_ome, ncols_eta))
                for i_row, image in enumerate(imgser_dict[det_key]):
                    psum = np.zeros(len(ring_map))
                    for i_k, k in enumerate(ring_map):
                        pdata = image[k[0], k[1]]
                        if threshold:
                            pdata[pdata <= threshold] = 0
                        psum[i_k] = np.average(pdata)
                    this_map[i_row, eta_idx[i_r]] = psum
                    # this_map[i_row, eta_idx[i_r]] = [
                    #         np.sum(image[k[0], k[1]]) for k in ring_map
                    #     ]
                ring_maps.append(this_map)
                pass
            ring_maps_panel[det_key] = ring_maps
        return ring_maps_panel, full_etas

    def extract_line_positions(self, plane_data, imgser_dict,
                               tth_tol=None, eta_tol=1., npdiv=2,
                               collapse_eta=True, collapse_tth=False,
                               do_interpolation=True):
        """
        TODO: handle wedge boundaries

        FIXME: must handle merged ranges!!!
        """
        if tth_tol is None:
            tth_tol = np.degrees(plane_data.tThWidth)
        tol_vec = 0.5*np.radians(
            [-tth_tol, -eta_tol,
             -tth_tol,  eta_tol,
             tth_tol,  eta_tol,
             tth_tol, -eta_tol])
        #
        # pbar = ProgressBar(
        #     widgets=[Bar('>'), ' ', ETA(), ' ', ReverseBar('<')],
        #     maxval=self.num_panels,
        # ).start()
        #
        panel_data = dict.fromkeys(self.detectors)
        for i_det, detector_id in enumerate(self.detectors):
            print("working on detector '%s'..." % detector_id)
            # pbar.update(i_det + 1)
            # grab panel
            panel = self.detectors[detector_id]
            instr_cfg = panel.config_dict(self.chi, self.tvec)
            native_area = panel.pixel_area  # pixel ref area
            images = imgser_dict[detector_id]
            if images.ndim == 2:
                n_images = 1
                images = np.tile(images, (1, 1, 1))
            elif images.ndim == 3:
                n_images = len(images)
            else:
                raise RuntimeError("images must be 2- or 3-d")

            # make rings
            pow_angs, pow_xys = panel.make_powder_rings(
                plane_data, merge_hkls=True, delta_eta=eta_tol)
            n_rings = len(pow_angs)

            ring_data = []
            for i_ring in range(n_rings):
                print("working on ring %d..." % i_ring)
                these_angs = pow_angs[i_ring]

                # make sure no one falls off...
                npts = len(these_angs)
                patch_vertices = (np.tile(these_angs, (1, 4)) +
                                  np.tile(tol_vec, (npts, 1))
                                  ).reshape(4*npts, 2)

                # find points that fall on the panel
                # WARNING: ignoring effect of crystal tvec
                det_xy, rMat_s = xrdutil._project_on_detector_plane(
                    np.hstack([patch_vertices, np.zeros((4*npts, 1))]),
                    panel.rmat, ct.identity_3x3, self.chi,
                    panel.tvec, ct.zeros_3, self.tvec,
                    panel.distortion)
                tmp_xy, on_panel = panel.clip_to_panel(det_xy)

                # all vertices must be on...
                patch_is_on = np.all(on_panel.reshape(npts, 4), axis=1)

                # reflection angles (voxel centers) and
                # pixel size in (tth, eta)
                ang_centers = these_angs[patch_is_on]
                ang_pixel_size = panel.angularPixelSize(tmp_xy[::4, :])

                # make the tth,eta patches for interpolation
                patches = xrdutil.make_reflection_patches(
                    instr_cfg, ang_centers, ang_pixel_size,
                    tth_tol=tth_tol, eta_tol=eta_tol,
                    distortion=panel.distortion,
                    npdiv=npdiv, quiet=True,
                    beamVec=self.beam_vector)

                # loop over patches
                # FIXME: fix initialization
                if collapse_tth:
                    patch_data = np.zeros((len(ang_centers), n_images))
                else:
                    patch_data = []
                for i_p, patch in enumerate(patches):
                    # strip relevant objects out of current patch
                    vtx_angs, vtx_xy, conn, areas, xy_eval, ijs = patch
                    if collapse_tth:
                        ang_data = (vtx_angs[0][0, [0, -1]],
                                    vtx_angs[1][[0, -1], 0])
                    else:
                        ang_data = (vtx_angs[0][0, :],
                                    ang_centers[i_p][-1])
                    prows, pcols = areas.shape
                    area_fac = areas/float(native_area)
                    # need to reshape eval pts for interpolation
                    xy_eval = np.vstack([
                        xy_eval[0].flatten(),
                        xy_eval[1].flatten()]).T

                    # interpolate
                    if not collapse_tth:
                        ims_data = []
                    for j_p, image in enumerate(images):
                        # catch interpolation type
                        if do_interpolation:
                            tmp = panel.interpolate_bilinear(
                                    xy_eval,
                                    image,
                                ).reshape(prows, pcols)*area_fac
                        else:
                            tmp = image[ijs[0], ijs[1]]*area_fac

                        # catch collapsing options
                        if collapse_tth:
                            patch_data[i_p, j_p] = np.sum(tmp)
                            # ims_data.append(np.sum(tmp))
                        else:
                            if collapse_eta:
                                ims_data.append(np.sum(tmp, axis=0))
                            else:
                                ims_data.append(tmp)
                        pass  # close image loop
                    if not collapse_tth:
                        patch_data.append((ang_data, ims_data))
                    pass  # close patch loop
                ring_data.append(patch_data)
                pass  # close ring loop
            panel_data[detector_id] = ring_data
            pass  # close panel loop
            # pbar.finish()
        return panel_data

    def simulate_laue_pattern(self, crystal_data,
                              minEnergy=5., maxEnergy=35.,
                              rmat_s=None, grain_params=None):
        """
        TODO: revisit output; dict, or concatenated list?
        """
        results = dict.fromkeys(self.detectors)
        for det_key, panel in self.detectors.iteritems():
            results[det_key] = panel.simulate_laue_pattern(
                crystal_data,
                minEnergy=minEnergy, maxEnergy=maxEnergy,
                rmat_s=rmat_s, tvec_s=self.tvec,
                grain_params=grain_params,
                beam_vec=self.beam_vector)
        return results

    def simulate_rotation_series(self, plane_data, grain_param_list,
                                 eta_ranges=[(-np.pi, np.pi), ],
                                 ome_ranges=[(-np.pi, np.pi), ],
                                 ome_period=(-np.pi, np.pi),
                                 wavelength=None):
        """
        TODO: revisit output; dict, or concatenated list?
        """
        results = dict.fromkeys(self.detectors)
        for det_key, panel in self.detectors.iteritems():
            results[det_key] = panel.simulate_rotation_series(
                plane_data, grain_param_list,
                eta_ranges=eta_ranges,
                ome_ranges=ome_ranges,
                ome_period=ome_period,
                chi=self.chi, tVec_s=self.tvec,
                wavelength=wavelength)
        return results

    def pull_spots(self, plane_data, grain_params,
                   imgser_dict,
                   tth_tol=0.25, eta_tol=1., ome_tol=1.,
                   npdiv=2, threshold=10,
                   eta_ranges=[(-np.pi, np.pi), ],
                   ome_period=(-np.pi, np.pi),
                   dirname='results', filename=None, output_format='text',
                   save_spot_list=False,
                   quiet=True, check_only=False,
                   interp='nearest'):
        """
        Exctract reflection info from a rotation series encoded as an
        OmegaImageseries object
        """

        # grain parameters
        rMat_c = makeRotMatOfExpMap(grain_params[:3])
        tVec_c = grain_params[3:6]

        # grab omega ranges from first imageseries
        #
        # WARNING: all imageseries AND all wedges within are assumed to have
        # the same omega values; put in a check that they are all the same???
        oims0 = imgser_dict[imgser_dict.keys()[0]]
        ome_ranges = [np.radians([i['ostart'], i['ostop']])
                      for i in oims0.omegawedges.wedges]

        # delta omega in DEGREES grabbed from first imageseries in the dict
        delta_ome = oims0.omega[0, 1] - oims0.omega[0, 0]

        # make omega grid for frame expansion around reference frame
        # in DEGREES
        ndiv_ome, ome_del = make_tolerance_grid(
            delta_ome, ome_tol, 1, adjust_window=True,
        )

        # generate structuring element for connected component labeling
        if ndiv_ome == 1:
            label_struct = ndimage.generate_binary_structure(2, 2)
        else:
            label_struct = ndimage.generate_binary_structure(3, 3)

        # simulate rotation series
        sim_results = self.simulate_rotation_series(
            plane_data, [grain_params, ],
            eta_ranges=eta_ranges,
            ome_ranges=ome_ranges,
            ome_period=ome_period)

        # patch vertex generator (global for instrument)
        tol_vec = 0.5*np.radians(
            [-tth_tol, -eta_tol,
             -tth_tol,  eta_tol,
             tth_tol,  eta_tol,
             tth_tol, -eta_tol])

        # prepare output if requested
        if filename is not None and output_format.lower() == 'hdf5':
            this_filename = os.path.join(dirname, filename)
            writer = GrainDataWriter_h5(
                os.path.join(dirname, filename),
                self.write_config(), grain_params)

        # =====================================================================
        # LOOP OVER PANELS
        # =====================================================================
        iRefl = 0
        compl = []
        output = dict.fromkeys(self.detectors)
        for detector_id in self.detectors:
            # initialize text-based output writer
            if filename is not None and output_format.lower() == 'text':
                output_dir = os.path.join(
                    dirname, detector_id
                    )
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                this_filename = os.path.join(
                    output_dir, filename
                )
                writer = PatchDataWriter(this_filename)

            # grab panel
            panel = self.detectors[detector_id]
            instr_cfg = panel.config_dict(self.chi, self.tvec)
            native_area = panel.pixel_area  # pixel ref area

            # pull out the OmegaImageSeries for this panel from input dict
            ome_imgser = imgser_dict[detector_id]

            # extract simulation results
            sim_results_p = sim_results[detector_id]
            hkl_ids = sim_results_p[0][0]
            hkls_p = sim_results_p[1][0]
            ang_centers = sim_results_p[2][0]
            xy_centers = sim_results_p[3][0]
            ang_pixel_size = sim_results_p[4][0]

            # now verify that full patch falls on detector...
            # ???: strictly necessary?
            #
            # patch vertex array from sim
            nangs = len(ang_centers)
            patch_vertices = (
                np.tile(ang_centers[:, :2], (1, 4)) +
                np.tile(tol_vec, (nangs, 1))
            ).reshape(4*nangs, 2)
            ome_dupl = np.tile(
                ang_centers[:, 2], (4, 1)
            ).T.reshape(len(patch_vertices), 1)

            # find vertices that all fall on the panel
            det_xy, _ = xrdutil._project_on_detector_plane(
                np.hstack([patch_vertices, ome_dupl]),
                panel.rmat, rMat_c, self.chi,
                panel.tvec, tVec_c, self.tvec,
                panel.distortion)
            _, on_panel = panel.clip_to_panel(det_xy, buffer_edges=True)

            # all vertices must be on...
            patch_is_on = np.all(on_panel.reshape(nangs, 4), axis=1)
            patch_xys = det_xy.reshape(nangs, 4, 2)[patch_is_on]

            # re-filter...
            hkl_ids = hkl_ids[patch_is_on]
            hkls_p = hkls_p[patch_is_on, :]
            ang_centers = ang_centers[patch_is_on, :]
            xy_centers = xy_centers[patch_is_on, :]
            ang_pixel_size = ang_pixel_size[patch_is_on, :]

            # TODO: add polygon testing right here!
            # done <JVB 06/21/16>
            if check_only:
                patch_output = []
                for i_pt, angs in enumerate(ang_centers):
                    # the evaluation omegas;
                    # expand about the central value using tol vector
                    ome_eval = np.degrees(angs[2]) + ome_del

                    # ...vectorize the omega_to_frame function to avoid loop?
                    frame_indices = [
                        ome_imgser.omega_to_frame(ome)[0] for ome in ome_eval
                    ]
                    if -1 in frame_indices:
                        if not quiet:
                            msg = """
                            window for (%d%d%d) falls outside omega range
                            """ % tuple(hkls_p[i_pt, :])
                            print(msg)
                        continue
                    else:
                        these_vertices = patch_xys[i_pt]
                        ijs = panel.cartToPixel(these_vertices)
                        ii, jj = polygon(ijs[:, 0], ijs[:, 1])
                        contains_signal = False
                        for i_frame in frame_indices:
                            contains_signal = contains_signal or np.any(
                                ome_imgser[i_frame][ii, jj] > threshold
                            )
                        compl.append(contains_signal)
                        patch_output.append((ii, jj, frame_indices))
            else:
                # make the tth,eta patches for interpolation
                patches = xrdutil.make_reflection_patches(
                    instr_cfg, ang_centers[:, :2], ang_pixel_size,
                    omega=ang_centers[:, 2],
                    tth_tol=tth_tol, eta_tol=eta_tol,
                    rMat_c=rMat_c, tVec_c=tVec_c,
                    distortion=panel.distortion,
                    npdiv=npdiv, quiet=True,
                    beamVec=self.beam_vector)

                # GRAND LOOP over reflections for this panel
                patch_output = []
                for i_pt, patch in enumerate(patches):

                    # strip relevant objects out of current patch
                    vtx_angs, vtx_xy, conn, areas, xy_eval, ijs = patch

                    prows, pcols = areas.shape
                    nrm_fac = areas/float(native_area)
                    nrm_fac = nrm_fac / np.min(nrm_fac)

                    # grab hkl info
                    hkl = hkls_p[i_pt, :]
                    hkl_id = hkl_ids[i_pt]

                    # edge arrays
                    tth_edges = vtx_angs[0][0, :]
                    delta_tth = tth_edges[1] - tth_edges[0]
                    eta_edges = vtx_angs[1][:, 0]
                    delta_eta = eta_edges[1] - eta_edges[0]

                    # need to reshape eval pts for interpolation
                    xy_eval = np.vstack([xy_eval[0].flatten(),
                                         xy_eval[1].flatten()]).T

                    # the evaluation omegas;
                    # expand about the central value using tol vector
                    ome_eval = np.degrees(ang_centers[i_pt, 2]) + ome_del

                    # ???: vectorize the omega_to_frame function to avoid loop?
                    frame_indices = [
                        ome_imgser.omega_to_frame(ome)[0] for ome in ome_eval
                    ]

                    if -1 in frame_indices:
                        if not quiet:
                            msg = """
                            window for (%d%d%d) falls outside omega range
                            """ % tuple(hkl)
                            print(msg)
                        continue
                    else:
                        # initialize spot data parameters
                        # !!! maybe change these to nan to not fuck up writer
                        peak_id = -999
                        sum_int = None
                        max_int = None
                        meas_angs = None
                        meas_xy = None

                        # quick check for intensity
                        contains_signal = False
                        patch_data_raw = []
                        for i_frame in frame_indices:
                            tmp = ome_imgser[i_frame][ijs[0], ijs[1]]
                            contains_signal = contains_signal or np.any(
                                tmp > threshold
                            )
                            patch_data_raw.append(tmp)
                            pass
                        patch_data_raw = np.stack(patch_data_raw, axis=0)
                        compl.append(contains_signal)

                        if contains_signal:
                            # initialize patch data array for intensities
                            if interp.lower() == 'bilinear':
                                patch_data = np.zeros(
                                    (len(frame_indices), prows, pcols))
                                for i, i_frame in enumerate(frame_indices):
                                    patch_data[i] = \
                                        panel.interpolate_bilinear(
                                            xy_eval,
                                            ome_imgser[i_frame],
                                            pad_with_nans=False
                                        ).reshape(prows, pcols)  # * nrm_fac
                            elif interp.lower() == 'nearest':
                                patch_data = patch_data_raw  # * nrm_fac
                            else:
                                msg = "interpolation option " + \
                                    "'%s' not understood"
                                raise(RuntimeError, msg % interp)

                            # now have interpolated patch data...
                            labels, num_peaks = ndimage.label(
                                patch_data > threshold, structure=label_struct
                            )
                            slabels = np.arange(1, num_peaks + 1)

                            if num_peaks > 0:
                                peak_id = iRefl
                                coms = np.array(
                                    ndimage.center_of_mass(
                                        patch_data,
                                        labels=labels,
                                        index=slabels
                                    )
                                )
                                if num_peaks > 1:
                                    center = np.r_[patch_data.shape]*0.5
                                    center_t = np.tile(center, (num_peaks, 1))
                                    com_diff = coms - center_t
                                    closest_peak_idx = np.argmin(
                                        np.sum(com_diff**2, axis=1)
                                    )
                                else:
                                    closest_peak_idx = 0
                                    pass  # end multipeak conditional
                                coms = coms[closest_peak_idx]
                                # meas_omes = \
                                #     ome_edges[0] + (0.5 + coms[0])*delta_ome
                                meas_omes = \
                                    ome_eval[0] + coms[0]*delta_ome
                                meas_angs = np.hstack(
                                    [tth_edges[0] + (0.5 + coms[2])*delta_tth,
                                     eta_edges[0] + (0.5 + coms[1])*delta_eta,
                                     mapAngle(
                                         np.radians(meas_omes), ome_period
                                         )
                                     ]
                                )

                                # intensities
                                #   - summed is 'integrated' over interpolated
                                #     data
                                #   - max is max of raw input data
                                sum_int = np.sum(
                                    patch_data[
                                        labels == slabels[closest_peak_idx]
                                    ]
                                )
                                max_int = np.max(
                                    patch_data_raw[
                                        labels == slabels[closest_peak_idx]
                                    ]
                                )
                                # ???: Should this only use labeled pixels?
                                # Those are segmented from interpolated data,
                                # not raw; likely ok in most cases.

                                # need MEASURED xy coords
                                gvec_c = anglesToGVec(
                                    meas_angs,
                                    chi=self.chi,
                                    rMat_c=rMat_c,
                                    bHat_l=self.beam_vector)
                                rMat_s = makeOscillRotMat(
                                    [self.chi, meas_angs[2]]
                                )
                                meas_xy = gvecToDetectorXY(
                                    gvec_c,
                                    panel.rmat, rMat_s, rMat_c,
                                    panel.tvec, self.tvec, tVec_c,
                                    beamVec=self.beam_vector)
                                if panel.distortion is not None:
                                    # FIXME: distortion handling
                                    meas_xy = panel.distortion[0](
                                        np.atleast_2d(meas_xy),
                                        panel.distortion[1],
                                        invert=True).flatten()
                                    pass
                                # FIXME: why is this suddenly necessary???
                                meas_xy = meas_xy.squeeze()
                                pass  # end num_peaks > 0
                        else:
                            patch_data = patch_data_raw
                            pass  # end contains_signal
                        # write output
                        if filename is not None:
                            if output_format.lower() == 'text':
                                writer.dump_patch(
                                    peak_id, hkl_id, hkl, sum_int, max_int,
                                    ang_centers[i_pt], meas_angs,
                                    xy_centers[i_pt], meas_xy)
                            elif output_format.lower() == 'hdf5':
                                xyc_arr = xy_eval.reshape(
                                    prows, pcols, 2
                                ).transpose(2, 0, 1)
                                writer.dump_patch(
                                    detector_id, iRefl, peak_id, hkl_id, hkl,
                                    tth_edges, eta_edges, np.radians(ome_eval),
                                    xyc_arr, ijs, frame_indices, patch_data,
                                    ang_centers[i_pt], xy_centers[i_pt],
                                    meas_angs, meas_xy)
                            pass  # end conditional on write output
                        pass  # end conditional on check only
                        patch_output.append([
                                peak_id, hkl_id, hkl, sum_int, max_int,
                                ang_centers[i_pt], meas_angs, meas_xy,
                                ])
                        iRefl += 1
                    pass  # end patch conditional
                pass  # end patch loop
            output[detector_id] = patch_output
            if filename is not None and output_format.lower() == 'text':
                writer.close()
            pass  # end detector loop
        if filename is not None and output_format.lower() == 'hdf5':
            writer.close()
        return compl, output

    """def fit_grain(self, grain_params, data_dir='results'):"""

    pass  # end class: HEDMInstrument
