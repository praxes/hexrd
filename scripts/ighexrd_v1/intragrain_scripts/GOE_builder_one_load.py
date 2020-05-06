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

#%%
if __name__ == '__main__':
    #
    #  Run preprocessor
    #
    parser = argparse.ArgumentParser(
        description="batchfndorigrains")

    parser.add_argument('grain_num',
                        help="grain_num", type=int)
    parser.add_argument('yaml_file',
                        help="yaml file", type=str)
    parser.add_argument('sample_name',
                        help="sample name", type=str)
    parser.add_argument('initial_or_final',
	                help="intial or final", type=str)
    parser.add_argument('scan_num',
	                help="scan num to fit grains", type=int)


    args = parser.parse_args()
    cfg_filename = args.yaml_file
    scan_number = args.scan_num
    samp_name = args.sample_name
    grain = args.grain_num
    initial_or_final = args.initial_or_final

# %%
# =============================================================================
# START USER INPUT
# =============================================================================
#----- The following parameters are passed through as arguements if running from command line ------

#cfg_filename = 'GOE_ti705.yml'
#samp_name = 'ti7-05'
#scan_number = 68
#initial_or_final = final
#grain_= 0
#----------------------------------------------------------------------------------------------------

#NEW GOE VARIABLES --- USER SHOULD EDIT ONCE -- same for all loadsteps

scan_to_center_GOE = 11 #name of scan - typically taken to be the initial zero load step
misorientation_bnd = 3.0 #in degrees
misorientation_spacing = 0.20 #in degrees
id_analysisname = 'ti7-11'

#location of master grains.out
dir_string = '/nfs/chess/user/ken38/Ti7_project/ti7-11-1percent/'
master_scan = dir_string + 'ti7-11-scan-%d' % scan_to_center_GOE

#location and name of npz file output
npz_string = 'grain_%d' %grain + '_goe_map_data_%s.npz' % initial_or_final
goe_path = '/nfs/chess/user/ken38/Ti7_project/ti7-11-1percent/GOE/'
npz_save_dir = goe_path + '%s_goe/' % initial_or_final #for npz file
GOE_directory = goe_path + '%s_goe/' % initial_or_final #for grains.out file

analysis_id = goe_path + id_analysisname + '-grain-%d' % grain + '-%s' % initial_or_final

# make output directory if doesn't exist
if not os.path.exists(npz_save_dir):
    os.mkdir(npz_save_dir)

# %%
# =============================================================================
# END USER INPUT
# =============================================================================
# ------------------------------------------------------------------------------
#cfg file -- currently ignores image_series block

data_dir = os.getcwd()
fc_stem = "%s_%s_%%s*.npz" % (samp_name, scan_number)

make_max_frames = False
use_direct_search = False

# for clustering neighborhood
# FIXME
min_samples = 2

# maps options
clobber_maps = False
show_maps = False

# =============================================================================
# END USER INPUT
# =============================================================================
# %%
cfg = config.open(cfg_filename)[0]

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

# %%

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
# %%

maps_fname = analysis_id + "_maps.npz"
if os.path.exists(maps_fname) and not clobber_maps:
    eta_ome = EtaOmeMaps(maps_fname)
else:
    print("INFO:\tbuilding eta_ome maps")
    start = timeit.default_timer()

    # make eta_ome maps
    eta_ome = instrument.GenerateEtaOmeMaps(
        imsd, instr, plane_data,
        active_hkls=active_hkls, threshold=build_map_threshold,
        ome_period=cfg.find_orientations.omega.period)

    print("INFO:\t\t...took %f seconds" % (timeit.default_timer() - start))

    # save them
    eta_ome.save(maps_fname)
#%%
# =============================================================================
# BOX TEST POINT GENERATION
# =============================================================================
# =============================================================================
# Set up multiprocessing from yml
# =============================================================================

ncpus = cfg.multiprocessing
#Reload original data always to start from master grains.out

exp_maps = np.zeros([1,3])
grain_id = np.zeros([1,1])

grain_out = '/grains.out'
load_data_master = np.loadtxt(master_scan + grain_out)
exp_map1 = load_data_master[grain,3:6]
grain_id1 = load_data_master[grain,0]

exp_maps[0,:] = exp_map1
grain_id[0,:] = grain_id1

mis_amt=misorientation_bnd*np.pi/180
spacing=misorientation_spacing*np.pi/180

ori_pts = np.arange(-mis_amt, (mis_amt+(spacing*0.999)), spacing)
num_ori_grid_pts=ori_pts.shape[0]**3
num_oris = exp_maps.shape[0]

Xs0, Ys0, Zs0 = np.meshgrid(ori_pts, ori_pts, ori_pts)
grid0 = np.vstack([Xs0.flatten(), Ys0.flatten(), Zs0.flatten()]).T

exp_maps_expanded=np.zeros([num_ori_grid_pts*num_oris,3])



for ii in np.arange(num_oris):
    pts_to_use=np.arange(num_ori_grid_pts) + ii*num_ori_grid_pts
    exp_maps_expanded[pts_to_use,:] =grid0 + np.r_[exp_maps[ii,:]]

exp_maps=exp_maps_expanded

rMat_c = rot.quatOfExpMap(exp_maps.T)

qfib=rMat_c
print("INFO: will test %d quaternions using %d processes"
      % (qfib.shape[1], ncpus))

# %%
# =============================================================================
# ORIENTATION SCORING
# =============================================================================

if use_direct_search:
    def test_orientation_FF_init(params):
        global paramMP
        paramMP = params

    def test_orientation_FF_reduced(quat):
        """
        input parameters are [
        plane_data, instrument, imgser_dict,
        tth_tol, eta_tol, ome_tol, npdiv, threshold
        ]
        """
        plane_data = paramMP['plane_data']
        instrument = paramMP['instrument']
        imgser_dict = paramMP['imgser_dict']
        tth_tol = paramMP['tth_tol']
        eta_tol = paramMP['eta_tol']
        ome_tol = paramMP['ome_tol']
        npdiv = paramMP['npdiv']
        threshold = paramMP['threshold']

        phi = 2*np.arccos(quat[0])
        n = xfcapi.unitRowVector(quat[1:])
        grain_params = np.hstack([
            phi*n, cnst.zeros_3, cnst.identity_6x1,
        ])

        compl, scrap = instrument.pull_spots(
            plane_data, grain_params, imgser_dict,
            tth_tol=tth_tol, eta_tol=eta_tol, ome_tol=ome_tol,
            npdiv=npdiv, threshold=threshold,
            eta_ranges=np.radians(cfg.find_orientations.eta.range),
            ome_period=(-np.pi, np.pi),
            check_only=True)

        return sum(compl)/float(len(compl))

    params = dict(
            plane_data=plane_data,
            instrument=instr,
            imgser_dict=imsd,
            tth_tol=tth_tol,
            eta_tol=eta_tol,
            ome_tol=ome_tol,
            npdiv=npdiv,
            threshold=cfg.fit_grains.threshold)

    print("INFO:\tusing direct seach")
    pool = multiprocessing.Pool(ncpus, test_orientation_FF_init, (params, ))
    completeness = pool.map(test_orientation_FF_reduced, qfib.T)
    pool.close()
else:
    print("INFO:\tusing map search with paintGrid on %d processes"
          % ncpus)
    start = timeit.default_timer()

    completeness = indexer.paintGrid(
        qfib,
        eta_ome,
        etaRange=np.radians(cfg.find_orientations.eta.range),
        omeTol=np.radians(cfg.find_orientations.omega.tolerance),
        etaTol=np.radians(cfg.find_orientations.eta.tolerance),
        omePeriod=np.radians(cfg.find_orientations.omega.period),
        threshold=on_map_threshold,
        doMultiProc=ncpus > 1,
        nCPUs=ncpus
       )


    print("INFO:\t\t...took %f seconds" % (timeit.default_timer() - start))
completeness = np.array(completeness)

# %%
# =============================================================================
# SAVE AS NPZ IN NEW FOLDER
# =============================================================================

goe_box_quat = np.zeros([1,4, len(pts_to_use)])
goe_box_con = np.zeros([1,len(pts_to_use)])

#for grain in range (0,len(grain_id)) :
goe_box_quat[0,:,:] = qfib[:,:]
goe_box_con[0,:] = completeness[:]

np.savez(npz_save_dir + npz_string,goe_box_con=goe_box_con,goe_box_quat=goe_box_quat,Xs0=Xs0,Ys0=Ys0,Zs0=Zs0)

#%%#==============================================================================
#GRAINS.OUT #currently used for nf - will eliminate
#==============================================================================

#if not os.path.exists(cfg.analysis_dir):
#    os.makedirs(cfg.analysis_dir)

print("INFO:writing misorientation clouds to grain_id_#.out files" )

gw = instrument.GrainDataWriter(os.path.join(GOE_directory, 'grain_id_%i.out') % grain )
grain_params_list = []

for gid, q in enumerate(goe_box_quat[0,:,:].T):
	phi = 2*np.arccos(q[0])
	n = xfcapi.unitRowVector(q[1:])
	grain_params = np.hstack([phi*n, cnst.zeros_3, cnst.identity_6x1])
	com = goe_box_con[0,gid]
	gw.dump_grain(grain, com, 0., grain_params)
	grain_params_list.append(grain_params)
gw.close()
