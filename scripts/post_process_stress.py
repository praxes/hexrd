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
import cPickle as cpl

import numpy as np

import argparse

from hexrd              import matrixutil as mutil
from hexrd.xrd          import rotations  as rot

#%% Extract stress data from grains.out data


def post_process_stress(grain_data,c_mat_C,schmid_T_list=None):
    num_grains=grain_data.shape[0]

    stress_S=np.zeros([num_grains,6])
    stress_C=np.zeros([num_grains,6])
    hydrostatic=np.zeros([num_grains,1])
    pressure=np.zeros([num_grains,1])
    von_mises=np.zeros([num_grains,1])

    if schmid_T_list is not None:
        num_slip_systems=schmid_T_list.shape[0]
        RSS=np.zeros([num_grains,num_slip_systems])


    for jj in np.arange(num_grains):    
    
        expMap=np.atleast_2d(grain_data[jj,3:6]).T    
        strainTmp=np.atleast_2d(grain_data[jj,15:21]).T
        
        #Turn exponential map into an orientation matrix
        Rsc=rot.rotMatOfExpMap(expMap)

        strainTenS = np.zeros((3, 3), dtype='float64')
        strainTenS[0, 0] = strainTmp[0]
        strainTenS[1, 1] = strainTmp[1]
        strainTenS[2, 2] = strainTmp[2]
        strainTenS[1, 2] = strainTmp[3]
        strainTenS[0, 2] = strainTmp[4] 
        strainTenS[0, 1] = strainTmp[5]  
        strainTenS[2, 1] = strainTmp[3] 
        strainTenS[2, 0] = strainTmp[4] 
        strainTenS[1, 0] = strainTmp[5] 

                  
        strainTenC=np.dot(np.dot(Rsc.T,strainTenS),Rsc)
        strainVecC = mutil.strainTenToVec(strainTenC)
        
        
        #Calculate stress        
        stressVecC=np.dot(c_mat_C,strainVecC)
        stressTenC = mutil.stressVecToTen(stressVecC)  
        stressTenS = np.dot(np.dot(Rsc,stressTenC),Rsc.T)
        stressVecS = mutil.stressTenToVec(stressTenS)       
        
        #Calculate hydrostatic stress
        hydrostaticStress=(stressVecS[:3].sum()/3)       
        
        
        #Calculate Von Mises Stress
        devStressS=stressTenS-hydrostaticStress*np.identity(3)        
        vonMisesStress=np.sqrt((3/2)*(devStressS**2).sum())        
        
        
        #Project on to slip systems
        if schmid_T_list is not None:
            for ii in np.arange(num_slip_systems):        
                RSS[jj,ii]=np.abs((stressTenC*schmid_T_list[ii,:,:]).sum())
            
            
        stress_S[jj,:]=stressVecS.flatten()
        stress_C[jj,:]=stressVecC.flatten()
        
        hydrostatic[jj,0]=hydrostaticStress
        pressure[jj,0]=-hydrostaticStress
        von_mises[jj,0]=vonMisesStress
    
    stress_data=dict()    
    
    stress_data['stress_S']=stress_S
    stress_data['stress_C']=stress_C
    stress_data['hydrostatic']=hydrostatic
    stress_data['pressure']=pressure
    stress_data['von_mises']=von_mises
    
    if schmid_T_list is not None:
        stress_data['RSS']=RSS
        
    return stress_data

#%% Command Line Access
if __name__ == '__main__':
    """
    USAGE : python post_process_stress grains_file stiffness_mat_file output_file_stem schmid_tensors_file 
    """
    parser = argparse.ArgumentParser(description='Post Process HEXRD Grains File To Extract Stress Tensor and Associated Quantities (Assuming Small Strain, Linear Elasticty)')


    parser.add_argument('grains_file', type=str)
    parser.add_argument('stiffness_mat_file', type=str)
    parser.add_argument('output_file_stem', type=str)
    parser.add_argument('--schmid_tensors_file', type=str, default=None)    

    args = vars(parser.parse_args(sys.argv[1:]))
    
    
    grain_data=np.loadtxt(args['grains_file'])
    
    c_mat=np.loadtxt(args['stiffness_mat_file'])

    #Extract Schmid Tensors from txt file
    if args['schmid_tensors_file'] is not None:
        T_vec = np.atleast_2d(np.loadtxt(args['schmid_tensors_file']))
        num_ten=T_vec.shape[0]
        T=np.zeros([num_ten,3,3])
        for i in np.arange(num_ten):
            T[i,:,:]=T_vec[i,:].reshape([3,3])
    
        stress_data=post_process_stress(grain_data,c_mat,T)
       
    else:
        stress_data=post_process_stress(grain_data,c_mat)


    cpl.dump(stress_data, open( args['output_file_stem']+'.cpl', "wb" ) )
