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
hexrd_script_directory='###' #Needs post_process_stress.py from the scripts directory

c_mat_C_file='###' # text file containing the stiffness matrix in the crystal coordinate system (6x6)

schmid_tensor_directory='###'
schmid_tensor_files=['###',"###']

num_load_steps=###

processed_data_directory='###'
analysis_stem='###'



###############################################################################
#%% #Import Modules
import os, sys

import numpy as np

import copy

from matplotlib import pyplot as plt

import hexrd.fitting.fitpeak as fitpeaks

from scipy.stats import gaussian_kde

sys.path.append(hexrd_script_directory)
import post_process_stress as stress_proc



#%% Function For Extracting Strength
###############################################################################



def extract_strength(stress_data,ss_bnds,completeness_mat=None,threshold=0.8):   
    
    tau_star=np.zeros(len(stress_data))
    w_tau=np.zeros(len(stress_data))    
    for i in np.arange(len(stress_data)):
        
        if completeness_mat is not None:        
            grains_to_use=np.where(completeness_mat[i]>threshold)[0]
            max_rss=np.max(stress_data[i]['RSS'][grains_to_use,ss_bnds[0]:ss_bnds[-1]],1)
        else:
            max_rss=np.max(stress_data[i]['RSS'][:,ss_bnds[0]:ss_bnds[-1]],1)
    
        tau=np.linspace(0., 1.5*np.max(max_rss), num=2000)        

        G = gaussian_kde(max_rss)
        tau_pdf = G.evaluate(tau) 

    
        maxPt=np.argmax(tau_pdf)   
        tmp_pdf=copy.copy(tau_pdf)
        tmp_pdf[:maxPt]=np.max(tau_pdf)
    
        pfit=fitpeaks.fit_pk_parms_1d([np.max(tau_pdf),tau[maxPt]+1e7,45e6],tau,tmp_pdf,pktype='tanh_stepdown')    
    
        tau_star[i]=pfit[1]
        w_tau[i]=pfit[2]
        
    return tau_star, w_tau


#%% Function For Extracting Strength

def plot_strength_curve(tau_star,w_tau,macro_strain=None,plot_color='blue'):  
    if macro_strain is None:
        macro_strain=np.arange(len(tau_star))
    
    strain_fine=np.linspace(macro_strain[0],macro_strain[-1],1000)
    interp_tau_star=np.interp(strain_fine,macro_strain,tau_star)
    interp_w_tau=np.interp(strain_fine,macro_strain,w_tau)    

    
    plt.errorbar(strain_fine,interp_tau_star,yerr=interp_w_tau,color=plot_color, capthick=0)
    plt.plot(macro_strain,tau_star,'s--',markerfacecolor=plot_color,markeredgecolor='k',markeredgewidth=1,color='k')
    plt.plot(strain_fine,interp_tau_star+interp_w_tau,'k--',linewidth=2)
    plt.plot(strain_fine,interp_tau_star-interp_w_tau,'k--',linewidth=2)
    
    plt.grid()
    

#%% Loading Data
###############################################################################

#Load Stiffness Matrix
c_mat_C=np.loadtxt(c_mat_C_file)

#Load Schmid Tensors
num_ss_fams=len(schmid_tensor_files)
num_per_sys=[None]*num_ss_fams
T_vecs=[None]*num_ss_fams
for i in np.arange(num_ss_fams):
    T_vecs[i] = np.atleast_2d(np.loadtxt(os.path.join(schmid_tensor_directory,schmid_tensor_files[i])))
    num_per_sys[i]=T_vecs[i].shape[0]

num_ten=int(np.sum(num_per_sys))
T=np.zeros([num_ten,3,3])
counter=0
for j in np.arange(num_ss_fams):
    for i in np.arange(num_per_sys[j]):  
        T[counter,:,:]=T_vecs[j][i,:].reshape([3,3])
        counter+=1

#Load and Process Stress Data
stress_data=[None]*(num_load_steps)
completeness=[None]*(num_load_steps)

for i in np.arange(num_load_steps):
    print('Processing Load ' + str(i))   
    grain_data=np.loadtxt(os.path.join(processed_data_directory,analysis_stem + '%03d'%(i),'grains.out'))
    completeness[i]=grain_data[:,1]
    stress_data[i]=stress_proc.post_process_stress(grain_data,c_mat_C,T)



#%% Extract Strengths from Different Slip System Families
tau_star=[None]*num_ss_fams
w_tau=[None]*num_ss_fams
for i in np.arange(num_ss_fams):
    tau_star[i], w_tau[i]  = extract_strength(stress_data,[int(np.sum(num_per_sys[:i])),int(np.sum(num_per_sys[:(i+1)]))],completeness,0.8)
    
#%%  Plot Slip System Strength Curves
plt.close('all')
plot_strength_curve(tau_star[2],w_tau[2],macro_strain=None,plot_color='blue')
    
