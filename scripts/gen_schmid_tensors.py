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
#%%
import sys

import argparse

import numpy as np
import cPickle as cpl


from hexrd              import matrixutil as mutil
from hexrd.xrd          import symmetry   as sym

#%%   
def gen_schmid_tensors(pd,uvw,hkl):        

    # slip plane directions    
    slipdir  = mutil.unitVector( np.dot( pd.latVecOps['F'], uvw) ) #  2 -1 -1  0
    slipdir_sym  = sym.applySym(slipdir, pd.getQSym(), csFlag=False, cullPM=True, tol=1e-08)
    
    # slip plane plane normals
    n_plane = mutil.unitVector( np.dot( pd.latVecOps['B'], hkl ) )
    n_plane_sym = sym.applySym(n_plane, pd.getQSym(), csFlag=False, cullPM=True, tol=1e-08)

    
    num_slip_plane= n_plane_sym.shape[1]
    
    num_slip_sys=0
    for i in range(num_slip_plane):
        planeID = np.where(abs(np.dot(n_plane_sym[:, i],slipdir_sym)) < 1.e-8)[0]
        num_slip_sys +=planeID.shape[0]
        
    T= np.zeros((num_slip_sys, 3, 3))
    counter=0
        #
    for i in range(num_slip_plane):
        planeID = np.where(abs(np.dot(n_plane_sym[:, i],slipdir_sym)) < 1.e-8)[0]
        for j in np.arange(planeID.shape[0]):    
            T[counter, :, :] = np.dot(slipdir_sym[:, planeID[j]].reshape(3, 1), n_plane_sym[:, i].reshape(1, 3))
            counter+=1
    #Clean some round off errors        
    round_off_err=np.where(abs(T)<1e-8)
    T[round_off_err[0],round_off_err[1],round_off_err[2]]=0.

    return T
    
    
#%%    
    
if __name__ == '__main__':
    """
    USAGE : python genschmidtensors material_file material_name uvw hkl output_file_stem
    """
    parser = argparse.ArgumentParser(description='Generate a set of schmid tensors for a given slip direction [uvw] and slip plane (hkl)')


    parser.add_argument('mat_file_loc', type=str)
    parser.add_argument('mat_name', type=str)
    parser.add_argument('uvw', type=str)
    parser.add_argument('hkl', type=str)
    parser.add_argument('out_file', type=str)
    

    args = vars(parser.parse_args(sys.argv[1:]))
    
    
    mat_list = cpl.load(open(args['mat_file_loc'], 'r'))

    # need to find the index of the active material
    # ***PROBABLY WILL CHANGE TO DICT INSTEAD OF LIST
    mat_idx = np.where([mat_list[i].name == args['mat_name'] for i in range(len(mat_list))])[0]
    
    # grab plane data, and useful things hanging off of it
    pd = mat_list[mat_idx[0]].planeData    
    
    uvw=np.zeros([3,1])
    sign=1.
    increment=0
    for ii in args['uvw']:       
        if ii =='-':
            sign=-1.
        else:
            uvw[increment,0]=sign*float(ii)
            sign=1.
            increment+=1
    
    hkl=np.zeros([3,1])
    sign=1.
    increment=0
    for ii in args['hkl']:       
        if ii =='-':
            sign=-1.
        else:
            hkl[increment,0]=sign*float(ii)
            sign=1.
            increment+=1        
    
            
    T=gen_schmid_tensors(pd,uvw,hkl)
    
    
    f=open(args['out_file']+'.txt','w')

    for i in np.arange(T.shape[0]):
        f.write("%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e  \n" % \
            (T[i,0,0],T[i,0,1],T[i,0,2],T[i,1,0],T[i,1,1],T[i,1,2],T[i,2,0],T[i,2,1],T[i,2,2]))
    f.close()
        
