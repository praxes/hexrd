#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 11:19:45 2019

@author: ken38
"""

import numpy as np
#import h5py
from hexrd.grainmap import nfutil
#==============================================================================
# %% Load DATA from .npz (ONLY FOR RELOADING DATA)
#==============================================================================

voxel_spacing=0.005#in mm

px = 271
stack = 34
layers = 4
overlap_amt = 4
stack0= 30

half_bnds = ((stack0*voxel_spacing*layers)+(overlap_amt*voxel_spacing))/2
v_bnds = [-half_bnds,half_bnds]

#%%
#full_stitch = np.empty([stack*layers,px,px])

full_stitch = np.empty([0,px,px])
full_stitch_con = full_stitch
full_stitch_ori = full_stitch
full_stitch_mis = full_stitch
full_stitch_combo = full_stitch
#full_stitch_conori = full_stitch

npz_save_dir = '/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/NEAR-FIELD/'
npz_string = 'ti7-11-stitched-initial_intragrain'

#%%
output_dir='/nfs/chess/aux/user/ken38/Ti7_project/ti7-11-1percent/NEAR-FIELD/'

for i in range(layers,0,-1):
#np.savez(save_dir+save_name, mask_confidence = mask_confidence,mask_grain = mask_grain, mask_misorientation = mask_misorientation, mask_combined_index=mask_combined_index, grain_avg_orientation_list=grain_avg_orientation_list)
    output_stem='initial_nf_intragrain_diffvol_%d' % i
    hold = np.load(output_dir+output_stem+'.npz')
    grain_map = hold['mask_grain']
    misorientation_map = hold['mask_misorientation']
    confidence_map = hold['mask_confidence']
    #grain_map_ori= hold['grain_map']
    #confidence_map_ori= hold['confidence_map']
    combined_index = hold['mask_combined_index']

    if i == layers:
        mismap_red = misorientation_map
        grain_map_red = grain_map
        con_red = confidence_map
        combo_red = combined_index
        #oricon_red = confidence_map_ori

    if i < layers:
        for ii in range(0,overlap_amt):
            layer_idx = full_stitch_ori.shape[0]-1-ii
            for iii in range(0,grain_map.shape[0]):
                for iv in range(0,grain_map.shape[1]):
                    if full_stitch_con[layer_idx,iii,iv] < confidence_map[overlap_amt-1-ii,iii,iv]:
                        pass
                    else:
                        full_stitch_con[layer_idx,iii,iv] = confidence_map[overlap_amt-1-ii,iii,iv]
                        full_stitch_ori[layer_idx,iii,iv] = grain_map[overlap_amt-1-ii,iii,iv]
                        full_stitch_mis[layer_idx,iii,iv] = misorientation_map[overlap_amt-1-ii,iii,iv]
                        full_stitch_combo[layer_idx,iii,iv] = combined_index[overlap_amt-1-ii,iii,iv]


        grain_map_red = grain_map[overlap_amt:,:,:]
        con_red = confidence_map[overlap_amt:,:,:]
        combo_red = combined_index[overlap_amt:,:,:]
        mismap_red = misorientation_map[overlap_amt:,:,:]
        #oricon_red = confidence_map_ori[overlap_amt:,:,:]

    full_stitch_ori = np.append(full_stitch_ori,grain_map_red,axis=0)
    full_stitch_con = np.append(full_stitch_con,con_red,axis=0)
    full_stitch_mis = np.append(full_stitch_mis,mismap_red,axis=0)
    full_stitch_combo = np.append(full_stitch_combo,combo_red,axis=0)
    #full_stitch_conori = np.append(full_stitch_conori, oricon_red, axis=0)
#%%

test_crds, n_crds, Xs, Ys, Zs = nfutil.gen_nf_test_grid_tomo(grain_map.shape[1], grain_map.shape[2], v_bnds, voxel_spacing)

#==============================================================================
# %% SAVE PROCESSED GRAIN MAP DATA
#==============================================================================


#nfutil.save_nf_data(output_dir,output_stem,full_stitch_ori,full_stitch_con,Xs,Ys,Zs, full_stitch_mis, full_stitch_conori)  #ori_list,id_remap=id_remap)
np.savez(npz_save_dir + npz_string,grain_map=full_stitch_ori,confidence_map=full_stitch_con, misorientation_map=full_stitch_mis,combined_index=full_stitch_combo, Xs=Xs,Ys=Ys,Zs=Zs)

#==============================================================================
# %% SAVE DATA AS .H5
#==============================================================================

hf=h5py.File(output_dir+npz_string+'_data.h5', 'w')#('data.h5','w')
#save_dir+save_stem+'_data.h5'

g1=hf.create_group('group1')
g1.create_dataset('grain_map', data=full_stitch_ori)
g1.create_dataset('confidence_map', data=full_stitch_con)
#g1.create_dataset('original_confidence_map', data=full_stitch_conori)
g1.create_dataset('misorientation_map', data=full_stitch_mis)
g1.create_dataset('combined', data = full_stitch_combo)
g1.create_dataset('Xs', data=Xs)
g1.create_dataset('Ys', data=Ys)
g1.create_dataset('Zs', data=Zs)


hf.close()
