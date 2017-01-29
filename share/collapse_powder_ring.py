import os, copy

import yaml

import numpy as np

from skimage import io

from hexrd.xrd import transforms_CAPI as xfc
from hexrd.xrd import material

from hexrd import constants, imageseries

import instrument

from matplotlib import pyplot as plt
from matplotlib import cm

import hexrd.fitting.fitpeak
import hexrd.fitting.peakfunctions as pkfuncs
import scipy.optimize as optimize

def make_matl(mat_name, sgnum, lparms, hkl_ssq_max=100):
    matl = material.Material(mat_name)
    matl.sgnum = sgnum
    matl.latticeParameters = lparms
    matl.hklMax = hkl_ssq_max

    nhkls = len(matl.planeData.exclusions)
    matl.planeData.set_exclusions(np.zeros(nhkls, dtype=bool))
    return matl

#%%
#instr_cfg_file = open('./dexela2_new.yml', 'r')
#instr_cfg = yaml.load(instr_cfg_file)
#instr = instrument.HEDMInstrument(instr_cfg)
#
#data_path = './'
#img_series_root = 'ff2_00047'
instr_cfg_file = open('./ge_detector_new.yml', 'r')
instr_cfg = yaml.load(instr_cfg_file)
instr = instrument.HEDMInstrument(instr_cfg)

data_path = './'
img_series_root = 'ge_scan_114'
img_series = imageseries.open(
        os.path.join(
                img_series_root + '-fcache-dir', img_series_root + '-fcache.yml'
                )
        , 'frame-cache')

nframes = 240

average_frame = imageseries.stats.average(img_series)

#%%
wlen = constants.keVToAngstrom(instr_cfg['beam']['energy'])

matl = make_matl('LSHR', 225, [3.5905,])

pd = matl.planeData
pd.wavelength = instr_cfg['beam']['energy'] # takes keV
pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)

set_by_tth_max = False
if set_by_tth_max:
    pd.tThMax = np.radians(6.75)
    tth_del = np.radians(0.75)
    tth_avg = np.average(pd.getTTh())
    tth_lo = pd.getTTh()[0] - tth_del
    tth_hi = pd.getTTh()[-1] + tth_del
else:
    tth_lo = np.radians(5.)
    tth_hi = np.radians(7.)
    tth_avg = 0.5*(tth_hi + tth_lo)
    excl = np.logical_or(pd.getTTh() <= tth_lo, pd.getTTh() >= tth_hi) 
    pd.exclusions = excl

panel_id = instr_cfg['detectors'].keys()[0]

d = instr.detectors[panel_id]

pangs, pxys = d.make_powder_rings([tth_avg, ])

#tth, peta = d.pixel_angles
#Y, X = d.pixel_coords
#xy = np.vstack([X.flatten(), Y.flatten()]).T
aps = d.angularPixelSize(pxys[0])

print "min angular pixel sizes: %.4f, %.4f" \
    %(np.degrees(np.min(aps[:, 0])), np.degrees(np.min(aps[:, 1])))

#%% set from looking at GUI
tth_size = np.degrees(np.min(aps[:, 0]))
eta_size = np.degrees(np.min(aps[:, 1]))

tth0 = np.degrees(tth_avg) 
eta0 = 0.

tth_range = np.degrees(tth_hi - tth_lo)
eta_range = 360.

ntth = int(tth_range/tth_size)
neta = int(eta_range/eta_size)

tth_vec = tth_size*(np.arange(ntth) - 0.5*ntth - 1) + tth0
eta_vec = eta_size*(np.arange(neta) - 0.5*neta - 1) + eta0

angpts = np.meshgrid(eta_vec, tth_vec, indexing='ij')
gpts = xfc.anglesToGVec(
    np.vstack([
        np.radians(angpts[1].flatten()), 
        np.radians(angpts[0].flatten()), 
        np.zeros(neta*ntth)
        ]).T, bHat_l=d.bvec)

xypts = xfc.gvecToDetectorXY(
    gpts, 
    d.rmat, np.eye(3), np.eye(3), 
    d.tvec, np.zeros(3), np.zeros(3), 
    beamVec=d.bvec)

img2 = d.interpolate_bilinear(xypts, average_frame).reshape(neta, ntth)
img3 = copy.deepcopy(img2)
borders = np.isnan(img2)
img2[borders] = 0.
img3[borders] = 0.
img3 += np.min(img3) + 1
img3 = np.log(img3)
img3[borders] = np.nan

extent = (
    np.min(angpts[1]), np.max(angpts[1]), 
    np.min(angpts[0]), np.max(angpts[0])
)

fig, ax = plt.subplots(2, 1, sharex=True, sharey=False)
ax[0].imshow(img3.reshape(neta, ntth),
             interpolation='nearest',
             cmap=cm.plasma, vmax=None, 
             extent=extent, 
             origin='lower')
ax[1].plot(angpts[1][0, :], np.sum(img2, axis=0)/img2.size)
ax[0].axis('tight')
ax[0].grid(True)
ax[1].grid(True)
ax[0].set_ylabel(r'$\eta$ [deg]', size=18)
ax[1].set_xlabel(r'$2\theta$ [deg]', size=18)
ax[1].set_ylabel(r'Intensity (arbitrary)', size=18)

plt.show()



#%% Multipeak Kludge

def fit_pk_obj_1d_mpeak(p,x,f0,pktype,num_pks):  
    
    f=np.zeros(len(x))
    p=np.reshape(p,[num_pks,p.shape[0]/num_pks])
    for ii in np.arange(num_pks):
        if pktype == 'gaussian':
            f=f+pkfuncs._gaussian1d_no_bg(p[ii],x)
        elif pktype == 'lorentzian':
            f=f+pkfuncs._lorentzian1d_no_bg(p[ii],x)
        elif pktype == 'pvoigt':
            f=f+pkfuncs._pvoigt1d_no_bg(p[ii],x)
        elif pktype == 'split_pvoigt':
            f=f+pkfuncs._split_pvoigt1d_no_bg(p[ii],x)

    
    resd = f-f0
    return resd



#%%
#plt.close('all')

num_tth=len(pd.getTTh())

x=angpts[1][0, :]
f=np.sum(img2, axis=0)/img2.size
pktype='pvoigt'
num_pks=num_tth

ftol=1e-6
xtol=1e-6  

fitArgs=(x,f,pktype,num_pks)

tth=matl.planeData.getTTh()*180./np.pi


p0=np.zeros([num_tth,4])

for ii in np.arange(num_tth):
    pt=np.argmin(np.abs(x-tth[ii]))
    
    p0[ii,:]=[f[pt],tth[ii],0.1,0.5]                     



p, outflag = optimize.leastsq(fit_pk_obj_1d_mpeak, p0, args=fitArgs,ftol=ftol,xtol=xtol)

p=np.reshape(p,[num_pks,p.shape[0]/num_pks])
f_fit=np.zeros(len(x))

for ii in np.arange(num_pks):
    f_fit=f_fit+pkfuncs._pvoigt1d_no_bg(p[ii],x)


#plt.plot(x,f,'x')
#plt.hold('true')
#plt.plot(x,f_fit)
ax[1].plot(x, f_fit, 'm+', ms=1)

#%%
fit_tths = p[:, 1]
fit_dsps = 0.5*wlen/np.sin(0.5*np.radians(fit_tths))
nrml_strains = fit_dsps/pd.getPlaneSpacings() - 1.

print nrml_strains
print "avg normal strain: %.3e" %np.average(nrml_strains)