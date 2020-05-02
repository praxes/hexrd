#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 09:37:05 2018

@author: ken38
""
#original findorientations first
Created on Wed Mar 22 19:04:10 2017

@author: bernier2

"""
#%%
# MUST BE RUN FROM FOLDER WHERE ETA OMEGA MAPS LIVE
# THIS SCRIPT IS DESIGNED TO BE USED ONCE ALL ETA OME MAPS HAVE BEEN FORMED AND ONLY FOR SUBSEQUENT CLOUD SEARCHES.
# IF ETA OMEGA MAP IS NOT OKAY THEN PLEASE RUN WITH ANOTHER SCRIPT TO GENERATE - OR - CHANGE CLOBBER_MAPS TO TRUE
# it is better to not change clobber maps to true but simply generate a new series of eta_ome maps and change yml directory

#%%
from __future__ import print_function

import time
import logging

import os

import glob

import multiprocessing

import numpy as np

from scipy import ndimage

import timeit

import argparse

try:
    import dill as cpl
except(ImportError):
    import cPickle as cpl

import yaml

from hexrd import constants as cnst
from hexrd import config
from hexrd import imageseries
from hexrd.imageseries.omega import OmegaImageSeries
from hexrd import instrument
from hexrd.findorientations import \
    generate_orientation_fibers, \
    run_cluster
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import indexer
from matplotlib import pyplot as plt
from hexrd.xrd.xrdutil import EtaOmeMaps

from hexrd.xrd          import rotations  as rot

logger = logging.getLogger(__name__)

# just require scikit-learn?
have_sklearn = False
try:
    import sklearn
    vstring = sklearn.__version__.split('.')
    if vstring[0] == '0' and int(vstring[1]) >= 14:
        from sklearn.cluster import dbscan
        from sklearn.metrics.pairwise import pairwise_distances
        have_sklearn = True
except ImportError:
    pass

# plane data
def load_pdata(cpkl, key):
    with file(cpkl, "r") as matf:
        mat_list = cpl.load(matf)
    return dict(zip([i.name for i in mat_list], mat_list))[key].planeData

# images
def load_images(yml):
    return imageseries.open(yml, format="frame-cache", style="npz")

# instrument
def load_instrument(yml):
    with file(yml, 'r') as f:
        icfg = yaml.load(f)
    return instrument.HEDMInstrument(instrument_config=icfg)

"""
Created on Fri Dec  9 13:05:27 2016

@author: bernier2
"""

import os

import yaml

import h5py

import numpy as np

from scipy import ndimage
from scipy.linalg.matfuncs import logm

from hexrd.gridutil import cellIndices, make_tolerance_grid
from hexrd import matrixutil as mutil
from hexrd.valunits import valWUnit
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
from hexrd.xrd.crystallography import PlaneData
from hexrd import constants as ct

# from hexrd.utils.progressbar import ProgressBar, Bar, ETA, ReverseBar

# FIXME: distortion kludge
from hexrd.xrd.distortion import GE_41RT  # BAD, VERY BAD!!!

from skimage.draw import polygon

#%%
if __name__ == '__main__':

    Run preprocessor

    parser = argparse.ArgumentParser(
        description="batchfndorigrains")

    parser.add_argument('yaml_file',
                        help="yaml file", type=str)
    parser.add_argument('sample_name',
                        help="sample name", type=str)
    parser.add_argument('scan_num',
                        help="scan num to fit grains", type=int)

    args = parser.parse_args()
    cfg_filename = args.yaml_file
    scan_number = args.scan_num
    samp_name = args.sample_name
# %%
#----- The following parameters are passed through as arguements if running from command line ------
#cfg_filename = 'ti7-05-cloud-tvecs.yml'
#samp_name = 'ti7-05'
#scan_number = 68
#inital_or_final = 'final'#'initial'

# =============================================================================
# START USER INPUT
# =============================================================================

#NEW GOE VARIABLES --- USER SHOULD EDIT ONCE -- same for all loadsteps
start_scan = 11 #must match GOE
misorientation_bnd = 3.0 #must match GOE
misorientation_spacing = 0.25 #must match GOE
id_analysisname = 'ti7-11' # must match GOE builder

#location of master grains.out
dir_string = '/nfs/chess/user/ken38/Ti7_project/ti7-11-1percent/'
load_step_zero_dir = dir_string + 'ti7-11-scan-%d' % start_scan
save_folder = 'saved_GOEs_centered/'

#%%
#location and name of npz file output
npz_save_dir = dir_string + save_folder + inital_or_final + '/'
# make output directory if doesn't exist
if not os.path.exists(npz_save_dir):
    os.mkdir(npz_save_dir)
#%%
# ------------------------------------------------------------------------------
#cfg file -- currently ignores image_series block

data_dir = dir_string
fc_stem = "%s_%s_%%s*.npz" % (samp_name, scan_number)

make_max_frames = False
use_direct_search = False

# for clustering neighborhood
# FIXME
min_samples = 2

# maps options
clobber_maps = False
show_maps = False

#%% one grain only
grain_out = '/grains.out'
load_data_zero = np.loadtxt(load_step_zero_dir + grain_out)
grain_id = load_data_zero[145:,0]

#%% LOAD YML FILE
cfg = config.open(cfg_filename)[0]

#analysis_id = '%s_%s' % (
#    cfg.analysis_name.strip().replace(' ', '-'),
#    cfg.material.active.strip().replace(' ', '-'),
#    )

active_hkls = cfg.find_orientations.orientation_maps.active_hkls
if active_hkls == 'all':
    active_hkls = None

max_tth = cfg.fit_grains.tth_max
if max_tth:
    if type(cfg.fit_grains.tth_max) != bool:
        max_tth = np.degrees(float(max_tth))
else:
    max_tth = None

# load plane data
plane_data = load_pdata(cfg.material.definitions, cfg.material.active)
plane_data.tThMax = max_tth

# load instrument
instr = load_instrument(cfg.instrument.parameters)
det_keys = instr.detectors.keys()

# !!! panel buffer setting is global and assumes same typ of panel!
for det_key in det_keys:
    instr.detectors[det_key].panel_buffer = \
        np.array(cfg.fit_grains.panel_buffer)

# grab eta ranges
eta_ranges = cfg.find_orientations.eta.range

# for indexing
build_map_threshold = cfg.find_orientations.orientation_maps.threshold

on_map_threshold = cfg.find_orientations.threshold
fiber_ndiv = cfg.find_orientations.seed_search.fiber_ndiv
fiber_seeds = cfg.find_orientations.seed_search.hkl_seeds

tth_tol = np.degrees(plane_data.tThWidth)
eta_tol = cfg.find_orientations.eta.tolerance
ome_tol = cfg.find_orientations.omega.tolerance
# omega period...
# QUESTION: necessary???
ome_period = np.radians(cfg.find_orientations.omega.period)

npdiv = cfg.fit_grains.npdiv

compl_thresh = cfg.find_orientations.clustering.completeness
cl_radius = cfg.find_orientations.clustering.radius

# %

imsd = dict.fromkeys(det_keys)
for det_key in det_keys:
    fc_file = sorted(
        glob.glob(
            os.path.join(
                data_dir,
                fc_stem % det_key.lower()
            )
        )
    )
    if len(fc_file) != 1:
        raise(RuntimeError, 'cache file not found, or multiple found')
    else:
        ims = load_images(fc_file[0])
        imsd[det_key] = OmegaImageSeries(ims)


if make_max_frames:
    max_frames_output_name = os.path.join(
        data_dir,
        "%s_%d-maxframes.hdf5" % (samp_name, scan_number)
    )

    if os.path.exists(max_frames_output_name):
        os.remove(max_frames_output_name)

    max_frames = dict.fromkeys(det_keys)
    for det_key in det_keys:
        max_frames[det_key] = imageseries.stats.max(imsd[det_key])

    ims_out = imageseries.open(
            None, 'array',
            data=np.array([max_frames[i] for i in max_frames]),
            meta={'panels': max_frames.keys()}
        )
    imageseries.write(
            ims_out, max_frames_output_name,
            'hdf5', path='/imageseries'
        )

#%%
class GenerateEtaOmeMaps(object):
    """
    eta-ome map class derived from new image_series and YAML config

    ...for now...

    must provide:

    self.dataStore
    self.planeData
    self.iHKLList
    self.etaEdges # IN RADIANS
    self.omeEdges # IN RADIANS
    self.etas     # IN RADIANS
    self.omegas   # IN RADIANS

    """
    def __init__(self, grain, image_series_dict, instrument, plane_data,
                 eta_step=0.25, threshold=None,
                 ome_period=(0, 360)):
        """
        image_series must be OmegaImageSeries class
        instrument_params must be a dict (loaded from yaml spec)
        active_hkls must be a list (required for now)
        """
        grain_params = np.squeeze(load_data_zero[grain,:])

        analysis_id = id_analysisname + '-grain-%d' % grain + '-%s' % initial_or_final

        self._planeData = plane_data

        # ???: change name of iHKLList?
        # ???: can we change the behavior of iHKLList?
        active_hkls = [0,1,2,3,4]

        if active_hkls is None:
            n_rings = len(plane_data.getTTh())
            self._iHKLList = range(n_rings)
        else:
            self._iHKLList = active_hkls
            n_rings = len(active_hkls)

        # ???: need to pass a threshold?
        eta_mapping, etas = instrument.extract_polar_maps_grain(
            plane_data, image_series_dict, grain_params,
            active_hkls=active_hkls, threshold=threshold,
            tth_tol=None, eta_tol=eta_step)


        # grab a det key
        # WARNING: this process assumes that the imageseries for all panels
        # have the same length and omegas
        det_key = eta_mapping.keys()[0]
        data_store = []
        for i_ring in range(n_rings):
            full_map = np.zeros_like(eta_mapping[det_key][i_ring])
            nan_mask_full = np.zeros(
                (len(eta_mapping), full_map.shape[0], full_map.shape[1])
            )
            i_p = 0
            for det_key, eta_map in eta_mapping.iteritems():
                nan_mask = ~np.isnan(eta_map[i_ring])
                nan_mask_full[i_p] = nan_mask
                full_map[nan_mask] += eta_map[i_ring][nan_mask]
                i_p += 1
            re_nan_these = np.sum(nan_mask_full, axis=0) == 0
            full_map[re_nan_these] = np.nan
            data_store.append(full_map)
        self._dataStore = data_store

        # handle omegas
        omegas_array = image_series_dict[det_key].metadata['omega']
        self._omegas = mapAngle(
            np.radians(np.average(omegas_array, axis=1)),
            np.radians(ome_period)
        )
        self._omeEdges = mapAngle(
            np.radians(np.r_[omegas_array[:, 0], omegas_array[-1, 1]]),
            np.radians(ome_period)
        )

        # !!! must avoid the case where omeEdges[0] = omeEdges[-1] for the
        # indexer to work properly
        if abs(self._omeEdges[0] - self._omeEdges[-1]) <= ct.sqrt_epsf:
            # !!! SIGNED delta ome
            del_ome = np.radians(omegas_array[0, 1] - omegas_array[0, 0])
            self._omeEdges[-1] = self._omeEdges[-2] + del_ome

        # handle etas
        # WARNING: unlinke the omegas in imageseries metadata,
        # these are in RADIANS and represent bin centers
        self._etas = etas
        self._etaEdges = np.r_[
            etas - 0.5*np.radians(eta_step),
            etas[-1] + 0.5*np.radians(eta_step)]

        self.save(npz_save_dir + analysis_id + "_maps.npz")

    @property
    def dataStore(self):
        return self._dataStore

    @property
    def planeData(self):
        return self._planeData

    @property
    def iHKLList(self):
        return np.atleast_1d(self._iHKLList).flatten()

    @property
    def etaEdges(self):
        return self._etaEdges

    @property
    def omeEdges(self):
        return self._omeEdges

    @property
    def etas(self):
        return self._etas

    @property
    def omegas(self):
        return self._omegas

    def save(self, filename):
        """
        self.dataStore
        self.planeData
        self.iHKLList
        self.etaEdges
        self.omeEdges
        self.etas
        self.omegas
        """
        args = np.array(self.planeData.getParams())[:4]
        args[2] = valWUnit('wavelength', 'length', args[2], 'angstrom')
        hkls = self.planeData.hkls
        save_dict = {'dataStore': self.dataStore,
                     'etas': self.etas,
                     'etaEdges': self.etaEdges,
                     'iHKLList': self.iHKLList,
                     'omegas': self.omegas,
                     'omeEdges': self.omeEdges,
                     'planeData_args': args,
                     'planeData_hkls': hkls}
        np.savez_compressed(filename, **save_dict)
        return
    pass  # end of class: GenerateEtaOmeMaps

#%%
from multiprocessing import Pool
from functools import partial

num_processors=24
grain_id= list(np.array(grain_id).astype('int').T)
#%%
print('building eta_ome maps using multiprocessing...')

pool = Pool(processes=num_processors)

#active hkls hardcoded to [0,1,2,3,4]
eta_ome_partial = partial(GenerateEtaOmeMaps, image_series_dict=imsd, instrument=instr, plane_data=plane_data, threshold=build_map_threshold, ome_period=cfg.find_orientations.omega.period)

eta_ome = pool.map(eta_ome_partial, grain_id, chunksize=1)
pool.close()
