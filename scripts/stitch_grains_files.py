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
hexrd_script_directory='###'

material_file_loc='###' # hexrd material file in cpickle format
mat_name='###'

grain_file_locs=['###',\
    '###']#Can be more than 2

output_data=True
output_dir='###'


#Position and Misorientation differences to merge grains
dist=0.05 #mm
misorientation=1. #degrees
completeness_diff=0.1 #if two grains are matched, the completenesses are checked,
#if the differences in completion are within completeness_diff, the grain values are averaged,
#if not, the data from the grain with higher completion is kept and the other data is discarded


#Offsets, these can be input as arguments
#Each dataset can have positional or rotation offsets
#Position offsets are in mm, 3 x n matrix where n is the number of grains.out files being stitched
#Rotation offsets exponential maps

#3 x 2 examples
#pos_offset=np.array([[0.,0.],[0.,0.],[0.,0.]])
#rot_offset=np.array([[0.,0.],[0.,0.],[0.,0.]])
pos_offset=None
rot_offset=None

###############################################################################

#%% #Import Modules
import os

import numpy as np

import copy

import cPickle as cpl

from hexrd              import matrixutil as mutil
from hexrd.xrd          import rotations  as rot
from hexrd.xrd          import symmetry   as sym



#%%

def remove_duplicate_grains(grain_data,qsyms,dist_thresh=0.01,misorient_thresh=0.1,comp_diff=0.1):
    total_grains=grain_data.shape[0]

    all_comp=grain_data[:,1]
    grain_quats=rot.quatOfExpMap(grain_data[:,3:6].T)
    dup_list=np.array([])    
    
    print 'Removing duplicate grains'
    for i in np.arange(total_grains-1):
        cur_pos=grain_data[i,3:6]
        other_pos=grain_data[(i+1):,3:6]
        xdist=cur_pos[0]-other_pos[:,0]
        ydist=cur_pos[1]-other_pos[:,1]
        zdist=cur_pos[2]-other_pos[:,2]
        
        dist=np.sqrt(xdist**2.+ydist**2.+zdist**2.)
        
        if np.min(dist)<dist_thresh:
            
            q1=sym.toFundamentalRegion(np.atleast_2d(grain_quats[:,i]).T,crysSym=qsyms)
            q2=sym.toFundamentalRegion(np.atleast_2d(grain_quats[:,np.argmin(dist)+i+1]).T,crysSym=qsyms)
            misorientation=rot.misorientation(q1,q2)[0]*180./np.pi

            if misorientation < misorient_thresh:             
                dup_list=np.append(dup_list,np.argmin(dist)+i+1)
                if (np.abs(all_comp[i] - all_comp[np.argmin(dist)+i+1]) <=comp_diff):
                    grain_data[i,:]=(grain_data[i,:]+grain_data[np.argmin(dist)+i+1,:])/2.
                elif ( (all_comp[np.argmin(dist)+i+1]-all_comp[i]) >0):
                    grain_data[i,:]=grain_data[np.argmin(dist)+i+1,:]
    
    grain_data=np.delete(grain_data,dup_list,axis=0)

    print 'Removed %d Grains' % (len(dup_list))

    grain_data[:,0]=np.arange(grain_data.shape[0])

    return grain_data




#%%


def assemble_grain_data(grain_data_list,pos_offset=None,rotation_offset=None):
    num_grain_files=len(grain_data_list)    
    
    num_grains_list=[None]*num_grain_files
    
    for i in np.arange(num_grain_files):
        num_grains_list[i]=grain_data_list[i].shape[0]
        
    num_grains=np.sum(num_grains_list)  
    
    grain_data=np.zeros([num_grains,grain_data_list[0].shape[1]])    
    
    for i in np.arange(num_grain_files):
        
        tmp=copy.copy(grain_data_list[i])
        
        if pos_offset is not None:
            pos_tile=np.tile(pos_offset[:,i],[num_grains_list[i],1])
            tmp[:,6:9]=tmp[:,6:9]+pos_tile
        #Needs Testing    
        if rotation_offset is not None:  
            rot_tile=np.tile(np.atleast_2d(rotation_offset[:,i]).T,[1,num_grains_list[i]])
            quat_tile=rot.quatOfExpMap(rot_tile)
            grain_quats=rot.quatOfExpMap(tmp[:,3:6].T)
            new_quats=rot.quatProduct(grain_quats,quat_tile)
            
            sinang = mutil.columnNorm(new_quats[1:,:])
            ang=2.*np.arcsin(sinang)
            axis   = mutil.unitVector(new_quats[1:,:])
            tmp[:,3:6]=np.tile(np.atleast_2d(ang).T,[1,3])*axis.T

            
        grain_data[int(np.sum(num_grains_list[:i])):int(np.sum(num_grains_list[:(i+1)])),:]=tmp
    
    
    grain_data[:,0]=np.arange(num_grains)    
    return grain_data



#%% Load data

mat_list = cpl.load(open(material_file_loc, 'r'))
mat_idx = np.where([mat_list[i].name == mat_name for i in range(len(mat_list))])[0]

# grab plane data, and useful things hanging off of it
pd = mat_list[mat_idx[0]].planeData
qsyms=sym.quatOfLaueGroup(pd.getLaueGroup())


num_grain_files=len(grain_file_locs)

grain_data_list=[None]*num_grain_files

for i in np.arange(num_grain_files):
    
    
    grain_data_list[i]=np.loadtxt(os.path.join(grain_file_locs[i],'grains.out'))

grain_data=assemble_grain_data(grain_data_list,pos_offset,rot_offset)

grain_data=remove_duplicate_grains(grain_data,qsyms,dist,misorientation,completeness_diff)


#%% Write data

if output_data:
    f = open(os.path.join(output_dir, 'grains.out'), 'w')
    
    header_items = (
        'grain ID', 'completeness', 'chi2',
        'xi[0]', 'xi[1]', 'xi[2]', 'tVec_c[0]', 'tVec_c[1]', 'tVec_c[2]',
        'vInv_s[0]', 'vInv_s[1]', 'vInv_s[2]', 'vInv_s[4]*sqrt(2)',
        'vInv_s[5]*sqrt(2)', 'vInv_s[6]*sqrt(2)', 'ln(V[0,0])',
        'ln(V[1,1])', 'ln(V[2,2])', 'ln(V[1,2])', 'ln(V[0,2])', 'ln(V[0,1])',
        )
    len_items = []
    for i in header_items[1:]:
        temp = len(i)
        len_items.append(temp if temp > 19 else 19) # for %19.12g
    fmtstr = '#%13s  ' + '  '.join(['%%%ds' % i for i in len_items]) + '\n'
    f.write(fmtstr % header_items)
    for i in grain_data.shape[0]:
        res_items = (
            grain_data[i,0], grain_data[i,1], grain_data[i,2], grain_data[i,3], grain_data[i,4], grain_data[i,5],
            grain_data[i,6], grain_data[i,7], grain_data[i,8], grain_data[i,9],
            grain_data[i,10], grain_data[i,11], grain_data[i,12], grain_data[i,13],
            grain_data[i,14], grain_data[i,15], grain_data[i,16], grain_data[i,17], grain_data[i,18],
            grain_data[i,19], grain_data[i,20],
            )
        fmtstr = (
            '%14d  ' + '  '.join(['%%%d.12g' % i for i in len_items]) + '\n'
            )
        f.write(fmtstr % res_items)
        
    f.close()
