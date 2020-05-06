#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:13:57 2020

@author: ken38
"""

import numpy as np
import h5py
from hexrd.grainmap import nfutil
import matplotlib.pyplot as plt
#==============================================================================
# %% Input Parameters and Load DATA from .npz (FOR GENERATING OTHER PARAMETERS)
#==============================================================================
l=5 #diffraction volume - each of mine were labeled in the output stem
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
Xs = hold['Xs']
Ys = hold['Ys']
Zs = hold['Zs']
ori_list = hold['ori_list']
id_remap = hold['id_remap']

px = grain_map.shape[1]
stack = grain_map.shape[0]

stack0= stack-overlap_amt
half_bnds = ((stack0*voxel_spacing*num_diffvols)+(overlap_amt*voxel_spacing))/2
v_bnds = [-half_bnds,half_bnds]

#%% IDENTIFY LOW CONFIDENCE REGION
conf_threshold_high = 0.45 #please note, I had very poor quality data in this example so you will likely want a much higher number than this. 
conf_threshold_low = 0.0

low_conf = np.logical_and(confidence_map<conf_threshold_high, confidence_map>conf_threshold_low)

#%% CHECK THAT YOUR THRESHOLD IS IDENTIFYING AREAS YOU WANT
layer_no=30

plt.figure('area check')
plt.imshow(confidence_map[layer_no,:,:])
plt.hold('on')
plt.imshow(low_conf[layer_no,:,:], alpha = 0.2)

#%% ADDITIONAL CHECK - LOOK AT HOW THE AREAS WILL SEGMENT AND MAKE SURE IT DOESNT RECOGNIZE AS ONE CONTIGUOUS BLOB
layer_no=11

from skimage import measure
from skimage import filters

all_labels = measure.label(low_conf[0:32,:,:])
blob_labels = measure.label(low_conf[0:32,:,:], background = 0)

plt.figure('labels',figsize=(9, 3.5))
plt.subplot(131)
plt.imshow(low_conf[layer_no,:,:], cmap='gray')
plt.axis('off')
plt.subplot(132)
plt.imshow(all_labels[layer_no,:,:], cmap='nipy_spectral')
plt.axis('off')
plt.subplot(133)
plt.imshow(blob_labels[layer_no,:,:], cmap='nipy_spectral')
plt.axis('off')

plt.tight_layout()
plt.show()

#%% CREATE CENTROID MAP OF LOW CONFIDENCE REGION

from scipy import ndimage

blob_labels_2 = ndimage.label(low_conf[0:32,:,:])[0]
centroids_2 = ndimage.measurements.center_of_mass(low_conf[0:32,:,:], blob_labels_2, np.unique(blob_labels_2))

centroid_point_map = np.zeros(np.shape(confidence_map))
centroid_new = np.empty([0,3])
for i in range(1,len(centroids_2)):
    where=len(np.where(blob_labels_2==i)[0])
    if where>10:
        print i
        centroid_new = np.append(centroid_new,np.reshape(np.array(centroids_2[i]),[1,3]),axis=0)
        centroid_point_map[np.rint(centroids_2[i][0]).astype('int'),np.rint(centroids_2[i][1]).astype('int'), np.rint(centroids_2[i][2]).astype('int')] = 10

#%% CAN CHECK THE CENTROIDS ARE IN LOCATIONS EXPECTED. 
        
layer_no=10
plt.figure('area check')
plt.imshow(confidence_map[layer_no,:,:])
plt.hold('on')
plt.imshow(low_conf[layer_no,:,:], alpha = 0.2)
plt.imshow(centroid_point_map[layer_no,:,:], alpha = 0.5)

#%% SAVE FILE OF CENTROIDS

save_file = 'centroid_diff_%d_miss_grain.npy' % l
np.save(save_file,centroid_new)

#%%

