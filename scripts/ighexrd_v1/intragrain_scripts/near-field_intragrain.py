#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:48:29 2018

@author: ken38
"""

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import os
import copy

from hexrd import valunits
from hexrd.grainmap import nfutil

from hexrd.xrd import rotations  as rot
from hexrd.xrd import symmetry   as sym
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import xrdutil
from hexrd.xrd.transforms_CAPI import anglesToGVec, \
                                      makeRotMatOfExpMap, makeDetectorRotMat, makeOscillRotMat, \
                                      gvecToDetectorXY, detectorXYToGvec
import numba
import argparse
import contextlib
import multiprocessing
import tempfile
import shutil

import yaml
import cPickle as cpl

#==============================================================================
# %% INPUT FILES: Location of layer data - nf map and nf images
#==============================================================================

#location of specific layer of near-field .NPZ array (initial)
file_dir='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/'
npz_file_stem = 'initial_nf_uniform_diffvol_1'

#location of nearfield images
data_folder='/nfs/chess/raw/2018-1/f2/miller-774-1/ti7-11/7/nf/' #layer 1
img_start=31602#for 0.25 degree/steps and 5 s exposure, end up with 6 junk frames up front, this is the 7th
num_imgs=1440
img_nums=np.arange(img_start,img_start+num_imgs,1)

#location of detector file (.yml)
det_file='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/retiga.yml'
#location of material file (.cpl)
mat_file='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/materials.cpl'


#grain_id.out file generated per grain with misorientation orientations for guesses
grain_id_out_dir='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/GOE/initial_goe/'

grain_out_file='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-scan-14/grains.out'

#==============================================================================
# %% OUTPUT FILES: Location to save new .npz
#==============================================================================

output_dir = '/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/NEAR-FIELD/initial_volume_1/'

# make output directory if doesn't exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

#==============================================================================
# %% USER INPUT - X-RAY DATA experiment and analysis parameters (copy from nf)
#==============================================================================
x_ray_energy=61.332 #keV from experiment

#name of the material in materials.cpl file
mat_name='ti7al'

#keep at zero, cannot delete right now
misorientation_bnd=0.0 #degrees
misorientation_spacing=0.1 #degrees

beam_stop_width=0.6#mm, assumed to be in the center of the detector

ome_range_deg=[(0.,359.75)] #degrees

max_tth=-1. #degrees, if a negative number is input, all peaks that will hit the detector are calculated

#image processing
num_for_dark=250#num images to use for median data
threshold=6. #set to 7 for initial. Currently using dark image 'min'
num_erosions=2 #num iterations of images erosion, don't mess with unless you know what you're doing
num_dilations=3 #num iterations of images erosion, don't mess with unless you know what you're doing
ome_dilation_iter=1 #num iterations of 3d image stack dilations, don't mess with unless you know what you're doing

chunk_size=500#chunksize for multiprocessing, don't mess with unless you know what you're doing

cross_sectional_dim=1.35 #cross sectional to reconstruct (should be at least 20%-30% over sample width)
#voxel spacing for the near field reconstruction
voxel_spacing = 0.005 #in mm
##vertical (y) reconstruction voxel bounds in mm
v_bnds=[-0.085,0.085]
#v_bnds=[-0.,0.]

#==============================================================================
# %% Set threshold values for misorientation data (ff - clouds)
#==============================================================================

#thresholds for grains in reconstructions
comp_thresh=0.7 #only use orientations with completnesses ABOVE this threshold
chi2_thresh=1.0 #is not used; make sure value > 0

#==============================================================================
# %% LOAD GRAIN AND EXPERIMENT DATA
#==============================================================================

experiment, nf_to_ff_id_map  = nfutil.gen_trial_exp_data(grain_out_file,det_file,mat_file, x_ray_energy, mat_name, max_tth, comp_thresh, chi2_thresh, misorientation_bnd, \
                       misorientation_spacing,ome_range_deg, num_imgs, beam_stop_width)

#==============================================================================
# %% NEAR FIELD - MAKE MEDIAN DARK
#==============================================================================
print '>>>>>>>>>>>>>>>>>>loading images>>>>>>>>>>>>>>>>>>'
dark=nfutil.gen_nf_dark(data_folder,img_nums,num_for_dark,experiment.nrows,experiment.ncols,dark_type='median',num_digits=6)

#==============================================================================
# %% NEAR FIELD - LOAD IMAGE DATA AND PROCESS
#==============================================================================

image_stack=gen_nf_cleaned_image_stack(data_folder,img_nums,dark,ome_dilation_iter,threshold,experiment.nrows,experiment.ncols,num_digits=6)#,grey_bnds=(5,5),gaussian=4.5)

#==============================================================================
# %% Load STITCHED-DATA from .npz
#==============================================================================

print ('>>>>>>>>>>>>>>>>>>loading nf map>>>>>>>>>>>>>>>>>>')
hold = np.load(file_dir + npz_file_stem + '_grain_map_data.npz')

grain_map = hold['grain_map']
confidence_map = hold['confidence_map']
Xs = hold['Xs']
Ys = hold['Ys']
Zs = hold['Zs']
ori_list = hold['ori_list']
id_remap = hold['id_remap']

all_grains_in_layer = np.unique(grain_map)

#%%

print ('.......multi process..........')
#new gen_trial_data definition
progress_handler = nfutil.progressbar_progress_observer()
save_handler=nfutil.forgetful_result_handler()

controller = nfutil.ProcessController(save_handler, progress_handler,
                                   ncpus=44, chunk_size=chunk_size)

multiprocessing_start_method = 'fork' if hasattr(os, 'fork') else 'spawn'

#%%
mis_all = np.zeros(grain_map.shape)
confidence_index_new = np.copy(confidence_map)
grain_map_new = np.copy(grain_map)
compiled_map = np.copy(grain_map.astype('float'))

#%%
import scipy.ndimage.morphology as morphology

#%%
for grain in range(1,all_grains_in_layer.shape[0]):
    #GRAIN ID -- USER INPUT
    grain_id = all_grains_in_layer[grain]

    print ('>>>>>>>>>>>>>>>>iteration %d>>>>>>>>>>>>>>>' % grain)
    print ('>>>>>>>>>>>>>>>>>>grain %d>>>>>>>>>>>>>>>>>' % grain_id)

    grain_id_out_file= grain_id_out_dir + 'grain_id_%s.out' % (grain_id)

    #LOAD EXP_MAPS CONFIDENCE TO DECIDE THRESHOLDS
    ori_out = np.loadtxt(grain_id_out_file)
    ori_data = ori_out[:,3:6]
    ori_comp=ori_out[:,1]

    comp_thresh=np.amax(ori_comp)*0.8 #only use orientations with completnesses ABOVE this threshold
    chi2_thresh=1.0 #is not used; make sure value > 0

    #GENERATES GRAIN MASK FROM GRAIN_MAP
    grains_plot_binary_0 = np.copy(grain_map)

    grains_plot_binary_0[grains_plot_binary_0 < grain_id] = 0
    grains_plot_binary_0[grains_plot_binary_0 > grain_id] = 0
    grains_plot_binary_0[grains_plot_binary_0 == grain_id] = 1

    grains_plot_binary = morphology.binary_dilation(grains_plot_binary_0,iterations=5).astype('int')

    # GRAIN MASK PROCESSING - CREATE GRAIN TEST GRID
    test=np.where(grains_plot_binary)

    test_crd_grain=np.zeros([len(test[0]),3])

    for ii in np.arange(len(test[0])):
        test_crd_grain[ii,0]=Xs[test[0][ii],test[1][ii],test[2][ii]]
        test_crd_grain[ii,1]=Ys[test[0][ii],test[1][ii],test[2][ii]]
        test_crd_grain[ii,2]=Zs[test[0][ii],test[1][ii],test[2][ii]]

    print ('----number of test coordinates = %d -----' % test_crd_grain.shape[0])

    if test_crd_grain.shape[0] == 0:
        pass

    else:
        # LOAD GRAIN AND EXPERIMENT DATA
        experiment_g, nf_to_ff_id_map_g = gen_trial_exp_data(grain_id_out_file,det_file,mat_file, x_ray_energy, mat_name, max_tth, comp_thresh, chi2_thresh, \
                              ome_range_deg, num_imgs, beam_stop_width)

        print ('----number of orientations = %d -----' % len(experiment_g.exp_maps))

        if len(experiment_g.exp_maps) < 10:
            pass

        else:

            # INSTANTIATE CONTROLLER - RUN BLOCK NO EDITING
            progress_handler = nfutil.progressbar_progress_observer()
            save_handler=nfutil.forgetful_result_handler()

            controller = nfutil.ProcessController(save_handler, progress_handler,
                                           ncpus=44, chunk_size=chunk_size)

            multiprocessing_start_method = 'fork' if hasattr(os, 'fork') else 'spawn'
            global _multiprocessing_start_method
            _multiprocessing_start_method = 'fork'
            #==============================================================================
            #  TEST ORIENTATIONS WITHIN GRAIN
            #==============================================================================
            try:

                raw_confidence_mis=nfutil.test_orientations(image_stack, experiment_g, test_crd_grain,
                              controller,multiprocessing_start_method)

                #==============================================================================
                # PUT DATA BACK INTO MESH
                #==============================================================================
                #full mesh test_crds
                print(' >>>>>>>>>>>>>>>>>>>>>putting data back in mesh>>>>>>>>>>>>>>>>>>>>>')
                test_crd_all = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T

                raw_confidence_new = np.empty([experiment_g.n_grains,len(test_crd_all[:,0])])

                for iii in range(0,test_crd_grain.shape[0]):
                    grain_location = np.where((test_crd_all == test_crd_grain[iii,:]).all(axis=1))
                    raw_confidence_new[:,grain_location[0][0]] = raw_confidence_mis[:,iii]

                #==============================================================================
                #  MASK TEST COORDINATES
                #==============================================================================

                print('Compiling Confidence Map...')
                confidence_map_g=np.max(raw_confidence_new,axis=0).reshape(Xs.shape)
                grain_map_g=np.argmax(raw_confidence_new,axis=0).reshape(Xs.shape)
                #id_remap
                max_orientation_no=np.max(grain_map_g)
                grain_map_copy=np.copy(grain_map_g)
                print('Remapping grain ids to ff...')
                for ii in np.arange(max_orientation_no):
                    this_orientation=np.where(grain_map_g==ii)
                    grain_map_copy[this_orientation]=nf_to_ff_id_map_g[ii]
                    grain_map_g=grain_map_copy

                #==============================================================================
                #   MISORIENTATION
                #==============================================================================
                print('calculate misorientation')

                tmp_data_avg=np.loadtxt(grain_out_file)
                tmp_data_grain=np.loadtxt(grain_id_out_file)

                id_avg=tmp_data_avg[:,0]
                ori_avg=tmp_data_avg[:,3:6]
                id_mis=tmp_data_grain[:,0]
                ori_mis=tmp_data_grain[:,3:6]

                mis = np.zeros(grain_map_g.shape)

                q_mor = rot.quatOfExpMap(ori_mis.T)
                q_avg = rot.quatOfExpMap(ori_avg.T)

                material_file_loc = mat_file # hexrd material file in cpickle format
                mat_name='ti7al'

                mat_list = cpl.load(open(material_file_loc, 'r'))
                mat_idx = np.where([mat_list[i].name == mat_name for i in range(len(mat_list))])[0]

                # grab plane data, and useful things hanging off of it
                pd = mat_list[mat_idx[0]].planeData
                qsyms=sym.quatOfLaueGroup(pd.getLaueGroup())

                for w in range(0,len(test[0])):
                    q2 = np.atleast_2d(q_mor[:,grain_map_g[test[0][w],test[1][w],test[2][w]]]).T
                    q1 = np.atleast_2d(q_avg[:,grain_id]).T
                    mis[test[0][w],test[1][w],test[2][w]] = rot.misorientation(q1,q2)[0]*180./np.pi

            #    plt.imshow(confidence_map_g[20,:,:],cmap='gray')
            #    plt.hold(True)
            #    plt.imshow(mis[20,:,:], alpha = 0.5)
            #    plt.colorbar()
                #==============================================================================
                # put back into the master mesh
                #==============================================================================
                print ('put in master mesh')
                #fresh array
                empty_map_full = np.zeros(compiled_map.shape)
                #confidence_index_local = np.copy(empty_map_full)
                #grain_map_local = np.copy(empty_map_full)
                mis_local = np.copy(empty_map_full)
                compiled_map_local = np.copy(empty_map_full)
                #place in fresh array
                #confidence_index_local[test] = confidence_map_g[test]
                #grain_map_local[test] = grain_map_g[test]
                mis_local[test] = mis[test]
                compiled_map_local[test] = grain_map[test].astype('float') + grain_map_g[test].astype('float')/100000.

                print ('saving as npz')

                save_string = 'grain_%d_nf' % grain_id
                #np.savez('/nfs/chess/user/ken38/Ti7_project/nf_data_all/ti7-05-nf/grain-by-grain-nf/'+ save_string, combined_map = local_compiled_map[test], test_crds = test, local_grain_map = grain_map_g[test], local_confidence_map = confidence_map_g[test])
                np.savez(output_dir + save_string, compiled_map_local = compiled_map_local, test_crds = test, local_grain_map = grain_map_g, local_confidence_map = confidence_map_g, mis_local = mis_local)
                print ('loop end')

            except:
                pass
    #==============================================================================
    #% save grain mesh
    #==============================================================================
#%
#np.savez(output_ext,grain_map_new = grain_map_new,confidence_index_new = confidence_index_new,compiled_map = compiled_map,mis_all = mis_all,Xs=Xs,Ys=Ys,Zs=Zs, grain_map = grain_map, confidence_map = confidence_map)
