# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:09:20 2016

@author: bernier2
"""
import sys
import yaml
import progressbar

import numpy as np
from scipy import sparse

from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from hexrd import matrixutil as mutil

# from hexrd.gridutil import cellIndices, cellConnectivity, cellCentroids, compute_areas
import hexrd.gridutil as gridutil

from hexrd.xrd import rotations as rot
from hexrd.xrd import material
from hexrd.xrd.crystallography import processWavelength
from hexrd.xrd.symmetry import applySym
from hexrd.xrd import transforms as xf
from hexrd.xrd.transforms_CAPI import anglesToGVec, makeRotMatOfExpMap, \
                                      makeDetectorRotMat, makeOscillRotMat, \
                                      gvecToDetectorXY, detectorXYToGvec
#from hexrd.xrd.xrdutil import angularPixelSize, make_reflection_patches, \
#                              simulateGVecs, _project_on_detector_plane, \
#                              _fetch_hkls_from_planedata
import hexrd.xrd.xrdutil as xrdutil

from hexrd.xrd.crystallography import processWavelength
from hexrd import valunits
from hexrd.utils import profiler
import numba
from progressbar import ProgressBar, Percentage, Bar


USE_NUMBA = True

ref_gparams = np.r_[0., 0., 0., 1., 1., 1., 0., 0., 0.]

def make_matl(mat_name, sgnum, lparms, hkl_ssq_max=100):
    matl = material.Material(mat_name)
    matl.sgnum = sgnum
    matl.latticeParameters = lparms
    matl.hklMax = hkl_ssq_max

    nhkls = len(matl.planeData.exclusions)
    matl.planeData.set_exclusions(np.zeros(nhkls, dtype=bool))
    return matl

#==============================================================================
# %% USER OPTIONS
#==============================================================================
n_grains = 2
quats = np.array([[ 0.91836393,  0.90869942],
                  [ 0.33952917,  0.1834835 ],
                  [ 0.17216207,  0.10095837],
                  [ 0.10811041,  0.36111851]])
phis = 2.*np.arccos(quats[0, :])
ns = mutil.unitVector(quats[1:, :])
exp_maps = np.array([phis[i]*ns[:, i] for i in range(n_grains)])
rMat_c = rot.rotMatOfQuat(quats)

# base 1-micron grid for 50 micron cube
cvec = np.arange(-25, 26)
X, Y, Z = np.meshgrid(cvec, cvec, cvec)

crd0 = 1e-3*np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).T
crd1 = crd0 + np.r_[0.100, 0.100, 0]
crds = np.array([crd0, crd1])

# make grain parameters
grain_params = []
for i in range(n_grains):
    for j in range(len(crd0)):
        grain_params.append(
            np.hstack([exp_maps[i, :], crds[i][j, :], xf.vInv_ref.flatten()])
            )

# scan range and period
ome_period = (0, 2*np.pi)
ome_range = [ome_period,]
ome_step = np.radians(1.)
nframes = 0
for i in range(len(ome_range)):
    del_ome = ome_range[i][1]-ome_range[i][0]
    nframes += int((ome_range[i][1]-ome_range[i][0])/ome_step)
    pass
ome_edges = np.arange(nframes+1)*ome_step

#==============================================================================
# %% INSTRUMENT
#==============================================================================
# load config
instr_cfg = yaml.load(open('./retiga.yml', 'r'))

tiltAngles = instr_cfg['detector']['transform']['tilt_angles']
tVec_d = np.array(instr_cfg['detector']['transform']['t_vec_d']).reshape(3, 1)
chi = instr_cfg['oscillation_stage']['chi']
tVec_s = np.array(instr_cfg['oscillation_stage']['t_vec_s']).reshape(3, 1)

rMat_d = makeDetectorRotMat(tiltAngles)
rMat_s = makeOscillRotMat([chi, 0.])

pixel_size = instr_cfg['detector']['pixels']['size']

nrows = instr_cfg['detector']['pixels']['rows']
ncols = instr_cfg['detector']['pixels']['columns']

row_dim = pixel_size[0]*nrows # in mm 
col_dim = pixel_size[1]*ncols # in mm 

x_col_edges = pixel_size[1]*(np.arange(ncols+1) - 0.5*ncols)
y_row_edges = pixel_size[0]*(np.arange(nrows+1) - 0.5*nrows)[::-1]

panel_dims = [(-0.5*ncols*pixel_size[1],
               -0.5*nrows*pixel_size[0]),
              ( 0.5*ncols*pixel_size[1],
                0.5*nrows*pixel_size[0])]

# a bit overkill, but grab max two-theta from all pixel transforms
rx, ry = np.meshgrid(x_col_edges, y_row_edges)
gcrds = detectorXYToGvec(np.vstack([rx.flatten(), ry.flatten()]).T,
                         rMat_d, rMat_s,
                         tVec_d, tVec_s, np.zeros(3))
pixel_tth = gcrds[0][0]

detector_params = np.hstack([tiltAngles, tVec_d.flatten(), chi, tVec_s.flatten()])

distortion = None

#==============================================================================
# %% crystallography data
#==============================================================================
gold = make_matl('gold', 225, [4.0782,], hkl_ssq_max=200)
gold.beamEnergy = valunits.valWUnit("wavelength","ENERGY",52,"keV")
pd = gold.planeData
pd.exclusions = np.zeros(len(pd.exclusions), dtype=bool)
pd.tThMax = np.amax(pixel_tth)

#==============================================================================
# %% DIFFRACTION SIMULATION
#==============================================================================
def simulate_diffractions(grain_params):
    pbar = ProgressBar(widgets=['simulate_diffractions', Percentage(), Bar()],
                       maxval=len(grain_params)).start()

    image_stack = np.zeros((nframes, nrows, ncols), dtype=bool)
    for i in range(len(grain_params)):
        sim_results = xrdutil.simulateGVecs(pd,
                                            detector_params,
                                            grain_params[i],
                                            panel_dims=panel_dims,
                                            pixel_pitch=pixel_size,
                                            ome_range=ome_range,
                                            ome_period=ome_period,
                                            distortion=None)
        valid_ids, valid_hkl, valid_ang, valid_xy, ang_ps = sim_results
        j_pix = gridutil.cellIndices(x_col_edges, valid_xy[:, 0])
        i_pix = gridutil.cellIndices(y_row_edges, valid_xy[:, 1])
        k_frame = gridutil.cellIndices(ome_edges, valid_ang[:, 2])

        # assign intensity
        for j, k in enumerate(k_frame):
            image_stack[k][i_pix[j], j_pix[j]] = True
        pbar.update(i+1)
        pass
    pbar.finish()

    #np.save('gold_cubes.npy', image_stack)
    return image_stack

def get_simulate_diffractions(grain_params):
    filename = 'gold_cubes.npy'
    try:
        image_stack = np.load(filename)
    except Exception:
        image_stack = simulate_diffractions(grain_params)
        np.save(filename, image_stack)

    return image_stack

#==============================================================================
# %% ORIENTATION TESTING
#==============================================================================
def test_orientations(image_stack):
    panel_dims_expanded = [(-10, -10), (10, 10)]
    panel_buffer = 0.05 # mm

    pbar = ProgressBar(widgets=['evaluate diffraction angles', Percentage(), Bar()],
                       maxval=n_grains).start()
    # first evaluate diffraction angles from orientation list (fixed)
    all_angles = []
    for i in range(n_grains):
        gparams = np.hstack([exp_maps[i, :].flatten(), ref_gparams])
        sim_results = xrdutil.simulateGVecs(pd,
                                            detector_params,
                                            gparams,
                                            panel_dims=panel_dims_expanded,
                                            pixel_pitch=pixel_size,
                                            ome_range=ome_range,
                                            ome_period=ome_period,
                                            distortion=None)
        all_angles.append(sim_results[2])
        pbar.update(i+1)
        pass
    pbar.finish()

    # form test grid and make main loop over spatial coordinates.  
    #   ** base 5-micron grid for 250x250x250 micron search space
    cvec_s = 0.001*np.arange(-250, 251)[::5]
    Xs, Ys, Zs = np.meshgrid(cvec_s, cvec_s, cvec_s)
    test_crds = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T

    # biggest projected diameter for 5 micron cube
    max_diameter = np.sqrt(3)*0.005

    row_dilation = np.ceil( 0.5*max_diameter/float(pixel_size[0]) )
    col_dilation = np.ceil( 0.5*max_diameter/float(pixel_size[1]) )

    i_dil, j_dil = np.meshgrid(np.arange(-row_dilation, row_dilation + 1),
                               np.arange(-col_dilation, col_dilation + 1))
    i_dil = np.array([i_dil.flatten()], dtype=int)
    j_dil = np.array([j_dil.flatten()], dtype=int)

    # now the grand loop...
    confidence = np.zeros((n_grains, len(test_crds)))
    n_coords = len(test_crds)
    print('Grand loop over {0} test coords and {1} grains.'.format(len(test_crds), n_grains))
    n_coords = 100
    pbar = ProgressBar(widgets=['grand loop', Percentage(), Bar()],
                       maxval=n_coords).start()
    for icrd in range(n_coords):
        for igrn in range(n_grains):
            det_xy, rMat_ss = xrdutil._project_on_detector_plane(all_angles[igrn],
                                                                 rMat_d, rMat_c[igrn], chi,
                                                                 tVec_d, test_crds[icrd, :], tVec_s, 
                                                                 distortion)

            # find on spatial extent of detector
            xTest = np.logical_and(det_xy[:, 0] >= panel_dims[0][0] + panel_buffer,
                                   det_xy[:, 0] <= panel_dims[1][0] - panel_buffer)
            yTest = np.logical_and(det_xy[:, 1] >= panel_dims[0][1] + panel_buffer,
                                   det_xy[:, 1] <= panel_dims[1][1] - panel_buffer)

            onDetector = np.where(np.logical_and(xTest, yTest))[0]

            # pick who's valid
            row_indices = gridutil.cellIndices(y_row_edges, det_xy[onDetector, 1])
            col_indices = gridutil.cellIndices(x_col_edges, det_xy[onDetector, 0])
            frame_indices = gridutil.cellIndices(ome_edges, all_angles[igrn][onDetector, 2])

            # perform check
            tmp_confidence = np.zeros(len(frame_indices), dtype=bool)
            for iref, indices in enumerate(zip(frame_indices, row_indices, col_indices)):
                i_sup = indices[1] + i_dil
                j_sup = indices[2] + j_dil
                idx_mask = np.where(
                    np.logical_and(np.logical_and(i_sup >= 0, i_sup < nrows),
                                   np.logical_and(j_sup >= 0, j_sup < ncols))
                                   )[0]
                tmp_confidence[iref] = np.any(image_stack[indices[0]][i_sup[idx_mask], j_sup[idx_mask]])
                pass
            confidence[igrn, icrd] = sum(tmp_confidence)/float(len(tmp_confidence))
            pass
        pbar.update(icrd + 1)
        pass
    pbar.finish()


def main():
    image_stack = get_simulate_diffractions(grain_params)
    test_orientations(image_stack)


if __name__=='__main__':
    if len(sys.argv) > 1:
        profiler.instrument_all(sys.argv[1:])

    main()

    if len(sys.argv) > 0:
        profiler.dump_results(sys.argv[1:])
