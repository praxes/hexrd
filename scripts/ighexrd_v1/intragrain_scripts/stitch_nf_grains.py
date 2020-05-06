#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 13:23:50 2019

@author: ken38
"""

#%%
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

#LOCATION OF NF DATA IN AVERAGE MAPS
file_dir='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/'
npz_file_stem = 'initial_nf_uniform_diffvol_1'

#LOCATION OF INDIVIDUAL NF DATA NPZ FROM GOEs
npz_location = '/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/NEAR-FIELD/'

GOE_loc='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/GOE/initial_goe'
tomo_mask_file='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/tomo_mask.npy'

#SAVE DIRECTORY
save_dir = npz_location
#SAVE NAME
save_name =  'initial_nf_intragrain_diffvol_1.npz'

#%% LOAD NF DATA GENERATED FROM GRAIN.OUT
hold = np.load(file_dir + npz_file_stem + '_grain_map_data.npz')

grain_map = hold['grain_map']
confidence_map = hold['confidence_map']
Xs = hold['Xs']
Ys = hold['Ys']
Zs = hold['Zs']
ori_list = hold['ori_list']
id_remap = hold['id_remap']

all_grains_in_layer = np.unique(grain_map)

empty_map = np.zeros(grain_map.shape)
full_map_confidence = np.copy(empty_map)
full_map_combined = np.copy(empty_map)
full_misorientation = np.copy(empty_map)
full_grain_map = np.copy(empty_map)

#%%
for i in range(1,len(all_grains_in_layer)):

    grain = all_grains_in_layer[i]

    try:
        local = np.load(npz_location + 'grain_%d_nf.npz' % grain)
        local_test_crds = local['test_crds']
        local_grain_map = local['local_grain_map']
        local_confidence_map = local['local_confidence_map']
        local_combined_map = local['compiled_map_local']
        local_mis = local['mis_local']
        for ii in range (0,local_test_crds.shape[1]):
            if full_map_confidence[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]] > local_confidence_map[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]]:
                pass
            else :
                full_map_confidence[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]] = local_confidence_map[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]]
                full_map_combined[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]] = local_combined_map[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]]
                full_misorientation[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]] = local_mis[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]]
                full_grain_map[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]] = local_grain_map[local_test_crds[0][ii],local_test_crds[1][ii],local_test_crds[2][ii]]
    except:
        print grain


#%% plot
plt.close('all')

layer_to_plot = 30

plt.figure()
plt.imshow(confidence_map[layer_to_plot,:,:], vmax=1)
plt.figure()
plt.imshow(full_map_confidence[layer_to_plot,:,:], vmax=1)

plt.figure()
plt.imshow(full_map_confidence[layer_to_plot,:,:]-confidence_map[layer_to_plot,:,:], vmin=-0.5, vmax=0.5)

#%%
tomo_mask = np.load(tomo_mask_file)

#%%
layer_to_plot = 15
plt.close('all')
mask = np.where(tomo_mask == False)

mask_confidence = np.copy(full_map_confidence)
mask_confidence[:,mask[0],mask[1]]=-.001

mask_grain = np.copy(np.floor(full_map_combined))
mask_grain[:,mask[0],mask[1]] =-1

mask_misorientation=np.copy(full_misorientation)
mask_misorientation[:,mask[0],mask[1]]=-.001

mask_combined_index=np.copy(full_map_combined)
mask_combined_index[:,mask[0],mask[1]]=-1

mask_local_id = np.copy(full_grain_map)
mask_local_id[:,mask[0],mask[1]] = -1

#%%
plt.close('all')
#plt.figure()
#plt.imshow(mask_confidence[layer_to_plot,:,:], vmax=1, cmap='gray')
#plt.figure()
#plt.imshow(mask_grain[layer_to_plot,:,:])
#plt.figure()
#plt.imshow(mask_grain[layer_to_plot,:,:])#, vmax=1)#, cmap='gray')
plt.imshow(mask_misorientation[layer_to_plot,:,:], vmax=3, alpha = 1)
plt.colorbar()
#plt.figure()
#plt.imshow(mask_combined_index[layer_to_plot,:,:])

#%%
grain_avg_orientation_list = np.zeros([np.unique(mask_grain).shape[0], 4])
mask_combin_id = np.zeros([mask_confidence.shape[0],mask_confidence.shape[1],mask_confidence.shape[2]])
all_grains_in_mask = np.unique(mask_grain)

for iii in range(0,np.unique(mask_grain).shape[0]):
    grain = all_grains_in_mask[iii]
    if grain != -1:
        where = np.where(mask_grain==grain)
        mode_subgrain_id = stats.mode(mask_local_id[where].astype('int'))
        grain_id_out = np.loadtxt(GOE_loc+'/grain_id_%d.out' % grain)
        exp_map_mode = grain_id_out[mode_subgrain_id[0],3:6]
        grain_avg_orientation_list[iii] = [grain, exp_map_mode[0][0], exp_map_mode[0][1], exp_map_mode[0][2]]

#%%SAVE NEW NPZ
np.savez(save_dir+save_name, mask_confidence = mask_confidence,mask_grain = mask_grain, mask_misorientation = mask_misorientation, mask_combined_index=mask_combined_index, grain_avg_orientation_list=grain_avg_orientation_list)
