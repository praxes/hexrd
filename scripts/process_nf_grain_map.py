#%% Necessary Dependencies


import numpy as np

import matplotlib.pyplot as plt

import multiprocessing as mp

import os

from hexrd.grainmap import nfutil
from hexrd.grainmap import tomoutil
from hexrd.grainmap import vtkutil




#==============================================================================
# %% FILES TO LOAD -CAN BE EDITED
#==============================================================================
#These files are attached, retiga.yml is a detector configuration file
#The near field detector was already calibrated

#A materials file, is a cPickle file which contains material information like lattice
#parameters necessary for the reconstruction

det_file='/####/retiga.yml'

mat_file='/####/materials.cpl' 

#==============================================================================
# %% OUTPUT INFO -CAN BE EDITED
#==============================================================================


output_dir='/####/'
output_stem='####'

#==============================================================================
# %% NEAR FIELD DATA FILES -CAN BE EDITED
#==============================================================================

#These are the near field data files used for the reconstruction, a grains.out file
#from the far field analaysis is used as orientation guess for the grid that will 
#be used for the near field reconstruction
grain_out_file='/####/grains.out'


#Locations of near field images
data_folder='/###/'

img_start=##
num_imgs=1441

img_nums=np.arange(img_start,img_start+num_imgs,1)


#==============================================================================
# %% TOMOGRAPHY DATA FILES -CAN BE EDITED
#==============================================================================


#Locations of tomography bright field images
tbf_data_folder='/####/'

tbf_img_start=##
tbf_num_imgs=20

#Locations of tomography images
tomo_data_folder='/####/'

tomo_img_start=##
tomo_num_imgs=720

#==============================================================================
# %% USER OPTIONS -CAN BE EDITED
#==============================================================================

x_ray_energy=41.991 #keV

#name of the material for the reconstruction
mat_name='az31'

#reconstruction with misorientation included, for many grains, this will quickly
#make the reconstruction size unmanagable
misorientation_bnd=0.0 #degrees 
misorientation_spacing=0.25 #degress

beam_stop_width=0.55#mm, assumed to be in the center of the detector


ome_range_deg=[(0.,360.)] #degrees 


max_tth=-1. #degrees, if a negative number is input, all peaks that will hit the detector are calculated

#image processing
num_for_dark=250#num images to use for median data
threshold=3.
num_erosions=3 #num iterations of images erosion, don't mess with unless you know what you're doing
num_dilations=2 #num iterations of images erosion, don't mess with unless you know what you're doing
ome_dilation_iter=1 #num iterations of 3d image stack dilations, don't mess with unless you know what you're doing

chunk_size=500#chunksize for multiprocessing, don't mess with unless you know what you're doing

#thresholds for grains in reconstructions
comp_thresh=0.5 #only use orientations from grains with completnesses ABOVE this threshold
chi2_thresh=0.05 #only use orientations from grains BELOW this chi^2


#tomography options
layer_row=1024 # row of layer to use to find the cross sectional specimen shape

#Don't change these unless you know what you are doing, this will close small holes
#and remove noise
recon_thresh=0.00006#usually varies between 0.0001 and 0.0005
noise_obj_size=5000
min_hole_size=5000


cross_sectional_dim=1.35 #cross sectional to reconstruct (should be at least 20%-30% over sample width)
#voxel spacing for the near field reconstruction
voxel_spacing=0.005#in mm
##vertical (y) reconstruction voxel bounds in mm
v_bnds=[-0.005,0.005]




#==============================================================================
# %% LOAD GRAIN DATA
#==============================================================================

experiment = nfutil.gen_trial_exp_data(grain_out_file,det_file,mat_file, x_ray_energy, mat_name, max_tth, comp_thresh, chi2_thresh, misorientation_bnd, \
                       misorientation_spacing,ome_range_deg, num_imgs, beam_stop_width)

#==============================================================================
# %% TOMO PROCESSING - GENERATE BRIGHT FIELD
#==============================================================================

tbf=tomoutil.gen_bright_field(tbf_data_folder,tbf_img_start,tbf_num_imgs,experiment.nrows,experiment.ncols)

#==============================================================================
# %% TOMO PROCESSING - BUILD RADIOGRAPHS
#==============================================================================


rad_stack=tomoutil.gen_attenuation_rads(tomo_data_folder,tbf,tomo_img_start,tomo_num_imgs,experiment.nrows,experiment.ncols)
    

#==============================================================================
# %% TOMO PROCESSING - INVERT SINOGRAM
#==============================================================================

reconstruction_fbp=tomoutil.tomo_reconstruct_layer(rad_stack,cross_sectional_dim,layer_row=layer_row,\
                                                   start_tomo_ang=ome_range_deg[0][0],end_tomo_ang=ome_range_deg[0][1],\
                                                   tomo_num_imgs=tomo_num_imgs, center=experiment.detector_params[3])

#==============================================================================
# %% TOMO PROCESSING - VIEW RAW FILTERED BACK PROJECTION
#==============================================================================

plt.close('all')
plt.imshow(reconstruction_fbp,vmin=0.75e-4,vmax=2e-4)
#Use this image to view the raw reconstruction, estimate threshold levels. and
#figure out if the rotation axis position needs to be corrected


#==============================================================================
# %% TOMO PROCESSING - CLEAN TOMO RECONSTRUCTION
#==============================================================================

binary_recon=tomoutil.threshold_and_clean_tomo_layer(reconstruction_fbp,recon_thresh, noise_obj_size,min_hole_size)


#==============================================================================
# %%  TOMO PROCESSING - RESAMPLE TOMO RECONSTRUCTION
#==============================================================================


tomo_mask=tomoutil.crop_and_rebin_tomo_layer(binary_recon,recon_thresh,voxel_spacing,experiment.pixel_size[0],cross_sectional_dim)

#==============================================================================
# %%  TOMO PROCESSING - VIEW TOMO_MASK FOR SAMPLE BOUNDS
#==============================================================================
plt.close('all')
plt.imshow(tomo_mask,interpolation='none')

#==============================================================================
# %%  TOMO PROCESSING - CONSTRUCT DATA GRID
#==============================================================================

test_crds, n_crds, Xs, Ys, Zs = nfutil.gen_nf_test_grid_tomo(tomo_mask.shape[1], tomo_mask.shape[0], v_bnds, voxel_spacing)

#==============================================================================
# %% NEAR FIELD - MAKE MEDIAN DARK
#==============================================================================
dark=nfutil.gen_nf_dark(data_folder,img_nums,num_for_dark,experiment.nrows,experiment.ncols)


#==============================================================================
# %% NEAR FIELD - LOAD IMAGE DATA AND PROCESS
#==============================================================================

image_stack=nfutil.gen_nf_image_stack(data_folder,img_nums,dark,num_erosions,num_dilations,ome_dilation_iter,threshold,experiment.nrows,experiment.ncols)


#==============================================================================
# %% VIEW IMAGES FOR DEBUGGING TO LOOK AT IMAGE PROCESSING PARAMETERS
#==============================================================================
plt.close('all')
img_to_view=0
plt.imshow(image_stack[img_to_view,:,:],interpolation='none')


#==============================================================================
# %% INSTANTIATE CONTROLLER - RUN BLOCK NO EDITING
#==============================================================================



progress_handler = nfutil.progressbar_progress_observer()
save_handler=nfutil.forgetful_result_handler()
    
controller = nfutil.ProcessController(save_handler, progress_handler,
                               ncpus=mp.cpu_count(), chunk_size=chunk_size)

multiprocessing_start_method = 'fork' if hasattr(os, 'fork') else 'spawn'


#==============================================================================
# %% TEST ORIENTATIONS - RUN BLOCK NO EDITING
#==============================================================================


raw_confidence=nfutil.test_orientations(image_stack, experiment, test_crds,
                  controller,multiprocessing_start_method)
    
    

#==============================================================================
# %% POST PROCESS W WHEN TOMOGRAPHY HAS BEEN USED
#==============================================================================

grain_map, confidence_map = nfutil.process_raw_confidence(raw_confidence,Xs.shape,tomo_mask=tomo_mask)



#==============================================================================
# %% SAVE RAW CONFIDENCE FILES 
#============================================================================

#This will be a very big file, don't save it if you don't need it
nfutil.save_raw_confidence(output_dir,output_stem,raw_confidence)


#==============================================================================
# %% SAVE PROCESSED GRAIN MAP DATA
#==============================================================================

nfutil.save_nf_data(output_dir,output_stem,grain_map,confidence_map,Xs,Ys,Zs,experiment.exp_maps)

#%%

layer_no=0
nfutil.plot_ori_map(grain_map, confidence_map, experiment.exp_maps, layer_no)


#==============================================================================
# %% SAVE DATA AS VTK
#==============================================================================

vtkutil.output_grain_map_vtk(output_dir,[output_stem],output_stem,0.1)
