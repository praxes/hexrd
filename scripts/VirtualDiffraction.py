# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:03:24 2015

@author: pagan2
"""

#%%

import numpy as np
import hexrd.xrd.material as mat
import hexrd.xrd.crystallography as crys
import hexrd.xrd.transforms_CAPI as trans
import multiprocessing as mp


#%%


material=mat.Material()
material.beamEnergy=15
material.sgnum=227
material.latticeParameters=[5.4310,]
material.name='Silicon'

#%%

samplePos=np.array([[0],[0],[0]])
crysPos=np.array([[0],[0],[0]])
rMat_c=np.identity(3)
bMat=material.planeData.latVecOps['B']
wavelength=material.planeData.wavelength

material.planeData.t

#%%
omega0,omega1=trans.oscillAnglesOfHKLs(material.planeData.hkls.T, 0, rMat_c, bMat, wavelength)




#%%

def VirtDiffWorker




