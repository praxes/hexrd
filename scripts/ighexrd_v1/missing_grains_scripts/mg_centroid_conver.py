#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 14:42:02 2020

@author: ken38
"""

import numpy as np


l=5
output_dir='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/'
output_stem = 'ti7-11-initial_average_vol_%d' % l

#global information
voxel_spacing = 0.005# in mm
overlap_in_microns = 0.020
overlap_amt = overlap_in_microns/voxel_spacing
num_diffvols = 5

#diffraction_volume information
layer = 0
hold = np.load(output_dir+output_stem+'_grain_map_data.npz')

#==============================================================================
# %% Calculate Additional Parameters from Input (run without changing)
#==============================================================================
grain_map = hold['grain_map']
confidence_map = hold['confidence_map']
Xs_1 = hold['Xs']
Ys_1 = hold['Ys']
Zs_1 = hold['Zs']
ori_list = hold['ori_list']
id_remap = hold['id_remap']

px = grain_map.shape[1]
stack = grain_map.shape[0]

stack0= stack-overlap_amt
half_bnds = ((stack0*voxel_spacing*num_diffvols)+(overlap_amt*voxel_spacing))/2
v_bnds = [-half_bnds,half_bnds]

save_file = 'centroid_diff_%i_miss_grain.npy' #FILE NAME SAVED FROM MISSING_GRAIN_CENTROIDS.PY SCRIPT

layer_1_coordinates = np.load(save_file % 1).astype('int')
layer_2_coordinates = np.load(save_file % 2).astype('int')
layer_3_coordinates = np.load(save_file % 3).astype('int')
layer_4_coordinates = np.load(save_file % 4).astype('int')
layer_5_coordinates = np.load(save_file % 5).astype('int')
#%%%

test_crds_1 = np.zeros([layer_1_coordinates.shape[0],3])
for ii in np.arange(layer_1_coordinates.shape[0]):
    test_crds_1[ii,0]=Xs_1[layer_1_coordinates[ii,0],layer_1_coordinates[ii,1],layer_1_coordinates[ii,2]]
    test_crds_1[ii,1]=Ys_1[layer_1_coordinates[ii,0],layer_1_coordinates[ii,1],layer_1_coordinates[ii,2]]
    test_crds_1[ii,2]=Zs_1[layer_1_coordinates[ii,0],layer_1_coordinates[ii,1],layer_1_coordinates[ii,2]]

test_crds_2 = np.zeros([layer_2_coordinates.shape[0],3])
for ii in np.arange(layer_2_coordinates.shape[0]):
    test_crds_2[ii,0]=Xs_1[layer_2_coordinates[ii,0],layer_2_coordinates[ii,1],layer_2_coordinates[ii,2]]
    test_crds_2[ii,1]=Ys_1[layer_2_coordinates[ii,0],layer_2_coordinates[ii,1],layer_2_coordinates[ii,2]]
    test_crds_2[ii,2]=Zs_1[layer_2_coordinates[ii,0],layer_2_coordinates[ii,1],layer_2_coordinates[ii,2]]

test_crds_3 = np.zeros([layer_3_coordinates.shape[0],3])
for ii in np.arange(layer_3_coordinates.shape[0]):
    test_crds_3[ii,0]=Xs_1[layer_3_coordinates[ii,0],layer_3_coordinates[ii,1],layer_3_coordinates[ii,2]]
    test_crds_3[ii,1]=Ys_1[layer_3_coordinates[ii,0],layer_3_coordinates[ii,1],layer_3_coordinates[ii,2]]
    test_crds_3[ii,2]=Zs_1[layer_3_coordinates[ii,0],layer_3_coordinates[ii,1],layer_3_coordinates[ii,2]]

test_crds_4 = np.zeros([layer_4_coordinates.shape[0],3])
for ii in np.arange(layer_4_coordinates.shape[0]):
    test_crds_4[ii,0]=Xs_1[layer_4_coordinates[ii,0],layer_4_coordinates[ii,1],layer_4_coordinates[ii,2]]
    test_crds_4[ii,1]=Ys_1[layer_4_coordinates[ii,0],layer_4_coordinates[ii,1],layer_4_coordinates[ii,2]]
    test_crds_4[ii,2]=Zs_1[layer_4_coordinates[ii,0],layer_4_coordinates[ii,1],layer_4_coordinates[ii,2]]
    
test_crds_5 = np.zeros([layer_5_coordinates.shape[0],3])
for ii in np.arange(layer_5_coordinates.shape[0]):
    test_crds_5[ii,0]=Xs_1[layer_5_coordinates[ii,0],layer_5_coordinates[ii,1],layer_5_coordinates[ii,2]]
    test_crds_5[ii,1]=Ys_1[layer_5_coordinates[ii,0],layer_5_coordinates[ii,1],layer_5_coordinates[ii,2]]
    test_crds_5[ii,2]=Zs_1[layer_5_coordinates[ii,0],layer_5_coordinates[ii,1],layer_5_coordinates[ii,2]]
    
np.save('missing_coords_vol_1.npy',test_crds_1)
np.save('missing_coords_vol_2.npy',test_crds_2)
np.save('missing_coords_vol_3.npy',test_crds_3)
np.save('missing_coords_vol_4.npy',test_crds_4)
np.save('missing_coords_vol_5.npy',test_crds_5)

