# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details on downloading the source,
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
#%% User Input
###############################################################################

#File Locations
grains_file='###' #HEXRD Grains.out file
matl_file='###' #HEXRD Materials cpl
active_matl='###'
cfg_file='###' #HEXRD cfg yml
instr_file='###' #HEXRD instrument yml


#Output Image Location
output_location='###'
output_name='test'#Frame Cache Name

#Script Options
det_psf_fwhm=2.
cts_per_event=1000.
delta_ome = 0.25
min_I=5.


###############################################################################
#%%

import os


import cPickle as cpl

import numpy as np
import scipy as sp

import yaml

from hexrd import config

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd import transforms_CAPI as xfcapi


from hexrd.gridutil import cellIndices


from hexrd.xrd.xrdutil import simulateGVecs

#%% File Loading
cfg = config.open(cfg_file)[0] # NOTE: always a list of cfg objects

instr_cfg = yaml.load(open(instr_file, 'r'))

mat_list = cpl.load(open(matl_file, 'r'))

grain_params_list=np.loadtxt(grains_file)

#%% Extract Quantities from Loaded data

#Instrument Info
# stack into array for later
detector_params = np.hstack([
    instr_cfg['detector']['transform']['tilt_angles'],
    instr_cfg['detector']['transform']['t_vec_d'],
    instr_cfg['oscillation_stage']['chi'],
    instr_cfg['oscillation_stage']['t_vec_s'],
    ])

# pixel pitches
pixel_pitch = cfg.instrument.detector.pixels.size

# panel dimensions calculated from pixel pitches
row_dim = pixel_pitch[0]*cfg.instrument.detector.pixels.rows
col_dim = pixel_pitch[1]*cfg.instrument.detector.pixels.columns

# panel is ( (xmin, ymin), (xmax, ymax) )
panel_dims = (
    (-0.5*col_dim, -0.5*row_dim),
    ( 0.5*col_dim,  0.5*row_dim),
              )

detector_x_edges = np.arange(cfg.instrument.detector.pixels.columns+1)*pixel_pitch[1] + panel_dims[0][0]
detector_y_edges = np.arange(cfg.instrument.detector.pixels.rows+1)*pixel_pitch[0] + panel_dims[0][1]          

# UGH! hard-coded distortion...  still needs fixing when detector is rewritten
if instr_cfg['detector']['distortion']['function_name'] == 'GE_41RT':
    distortion = (dFuncs.GE_41RT,
                  instr_cfg['detector']['distortion']['parameters'],
                  )
else:
    distortion = None

#Image Info
nrows = int((panel_dims[1][1] - panel_dims[0][1]) / float(pixel_pitch[0]))
ncols = int((panel_dims[1][0] - panel_dims[0][0]) / float(pixel_pitch[1]))
row_edges = (np.arange(nrows+1)*pixel_pitch[0] + panel_dims[0][1])[::-1]
col_edges = np.arange(ncols+1)*pixel_pitch[1] + panel_dims[0][0]


nframes = int(360./float(delta_ome))
ome_edges = np.arange(nframes + 1)*delta_ome - 180.


#extract transform objects; rotations and translations
# detector first, rotation, then translation
#  - rotation takes comps from det frame to lab
rMat_d = xfcapi.makeDetectorRotMat(detector_params[:3])
tVec_d = np.r_[detector_params[3:6]]

# rotation stage (omega)
#  - chi is ccw tilt about lab X; rMat_s is omega dependent
#  - takes comps in sample to lab frame
chi = detector_params[6]
tVec_s = np.zeros((3,1))

# crystal; this will be a list of things, computed from quaternions
#  - trivial case here...
rMat_c = np.eye(3)
tVec_c = np.zeros((3,1))



#Material Info
mat_name = cfg.material.active # str that is the material name in database

# need to find the index of the active material
mat_idx = np.where([mat_list[i].name == mat_name for i in range(len(mat_list))])[0]

# grab plane data, and useful things hanging off of it
plane_data = mat_list[mat_idx].planeData
plane_data.tThMax=np.radians(20)
plane_data.set_exclusions(np.zeros(len(plane_data.exclusions), dtype=bool))


#%% Filters For Point Spread
def make_gaussian_filter(size,fwhm):
    sigma=fwhm/(2.*np.sqrt(2.*np.log(2.)))
#    size=[5,5]
#    sigma=1.    
    gaussFilter=np.zeros(size)
    cenRow=size[0]/2.
    cenCol=size[1]/2.
    
    pixRowCens=np.arange(size[0])+0.5
    pixColCens=np.arange(size[1])+0.5
    
    y=cenRow-pixRowCens
    x=pixColCens-cenCol
    
    xv, yv = np.meshgrid(x, y, sparse=False)
    
    r=np.sqrt(xv**2.+yv**2.)
    gaussFilter=np.exp(-r**2./(2*sigma**2))
    gaussFilter=gaussFilter/gaussFilter.sum()


    return gaussFilter

def make_lorentzian_filter(size,fwhm):
    
    gamma=fwhm/2.  
    
    lorentzianFilter=np.zeros(size)
    cenRow=size[0]/2.
    cenCol=size[1]/2.
    
    pixRowCens=np.arange(size[0])+0.5
    pixColCens=np.arange(size[1])+0.5
    
    y=cenRow-pixRowCens
    x=pixColCens-cenCol
    
    xv, yv = np.meshgrid(x, y, sparse=False)
    
    r=np.sqrt(xv**2.+yv**2.)
    lorentzianFilter=gamma**2 / ((r)**2 + gamma**2)
    lorentzianFilter=lorentzianFilter/lorentzianFilter.sum()


    return lorentzianFilter

#%%
#Calculate Intercepts for diffraction events from grains

pixel_data = []

for ii in np.arange(grain_params_list.shape[0]):
    print "processing grain %d..." %ii
    
    simg = simulateGVecs(plane_data, detector_params, grain_params_list[ii,3:15],distortion=None)
    
    valid_ids, valid_hkl, valid_ang, valid_xy, ang_ps = simg
    
    #ax.plot(valid_xy[:, 0], valid_xy[:, 1], 'b.', ms=2)
    this_frame = sp.sparse.coo_matrix((nrows, ncols), np.uint16)
    frame_indices = cellIndices(ome_edges, np.degrees(valid_ang[:, 2]))
    i_row = cellIndices(row_edges, valid_xy[:, 1])
    j_col = cellIndices(col_edges, valid_xy[:, 0])
    pixel_data.append(np.vstack([i_row, j_col, frame_indices]))
   

pixd = np.hstack(pixel_data)



frame_cache_data=[sp.sparse.coo_matrix([2048,2048],dtype='uint16')]*nframes

filter_size=np.round(det_psf_fwhm*5)
if filter_size % 2 == 0:
    filter_size+=1
    
psf_filter=make_gaussian_filter([filter_size,filter_size],det_psf_fwhm)

#Make pad four fast fourier tranform
filterPad=np.zeros((nrows, ncols), dtype=float)
filterPad[:psf_filter.shape[0],:psf_filter.shape[1]]=psf_filter
filterPadTransform=np.fft.fft2(filterPad)


#Build images and apply point spread
for i in np.arange(nframes):
    print "processing frame %d of %d" % (i,nframes)    
    
    this_frame = np.zeros((nrows, ncols), dtype=float)
    these_ij = pixd[:2, pixd[2, :] == i]
    
    this_frame[these_ij[0], these_ij[1]] += cts_per_event
    
    
    this_frame_transform=np.fft.fft2(this_frame)
    this_frame_convolved=np.real(np.fft.ifft2(this_frame_transform*filterPadTransform))  
    tmp=np.where(this_frame_convolved<min_I)
    this_frame_convolved[tmp]=0.
    frame_cache_data[i]=sp.sparse.coo_matrix(this_frame_convolved,dtype='uint16')


#%% Save Frame Cache Data
frame_cache=[None]*2
frame_cache[0]=frame_cache_data
frame_cache[1]=(0.,delta_ome*np.pi/180.)   
    
np.savez_compressed(os.path.join(output_location,output_name+'.npz'),frame_cache=frame_cache)

