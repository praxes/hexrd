import sys, os
import argparse

import yaml
import cPickle

import numpy as np

from hexrd import config
from hexrd import matrixutil as mutil

from hexrd.fitgrains import get_instrument_parameters

from hexrd.matrixutil import unitVector

# HARD-CODED DISTORTION!!!
from hexrd.xrd.distortion import GE_41RT

from hexrd.gridutil import cellIndices

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi

from hexrd.xrd import rotations as rot
from hexrd.xrd.symmetry import toFundamentalRegion
from hexrd.xrd.xrdutil import simulateOmeEtaMaps

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
plt.ioff()

def check_indexing_plots(cfg_filename, plot_trials=False, plot_from_grains=False):
    cfg = config.open(cfg_filename)[0]    # use first block, like indexing

    working_dir = cfg.working_dir
    analysis_dir = os.path.join(working_dir, cfg.analysis_name)
    
    #instrument parameters
    icfg = get_instrument_parameters(cfg)
    chi = icfg['oscillation_stage']['chi']

    # load maps that were used
    oem = cPickle.load(
        open(cfg.find_orientations.orientation_maps.file, 'r')
        )
    nmaps = len(oem.dataStore)
    omeEdges = np.degrees(oem.omeEdges); nome = len(omeEdges) - 1
    etaEdges = np.degrees(oem.etaEdges); neta = len(etaEdges) - 1
    delta_ome = abs(omeEdges[1]-omeEdges[0])

    full_ome_range = xf.angularDifference(omeEdges[0], omeEdges[-1]) == 0
    full_eta_range = xf.angularDifference(etaEdges[0], etaEdges[-1]) == 0
    
    # grab plane data and figure out IDs of map HKLS
    pd = oem.planeData
    gvids = [pd.hklDataList[i]['hklID'] for i in np.where(pd.exclusions == False)[0].tolist()]

    # load orientations
    quats = np.atleast_2d(np.loadtxt(os.path.join(working_dir, 'accepted_orientations.dat')))
    if plot_trials:
        compl = np.loadtxt(os.path.join(working_dir, 'completeness.dat'))
        tq = np.atleast_2d(np.loadtxt(os.path.join(working_dir, 'trial_orientations.dat')))
        quats = tq[compl >= 0., :]
        pass
    expMaps = np.tile(2. * np.arccos(quats[:, 0]), (3, 1))*unitVector(quats[:, 1:].T)

    ##########################################
    #      SPECIAL CASE FOR FIT GRAINS       #
    ##########################################
    if plot_from_grains:
        distortion = (GE_41RT, icfg['detector']['distortion']['parameters'])
        #
        grain_table = np.atleast_2d(np.loadtxt(os.path.join(analysis_dir, 'grains.out')))
        ngrains = len(grain_table)
        #
        expMaps = grain_table[:, 3:6]
        tVec_c = grain_table[:, 6:9]
        vInv = grain_table[:, 6:12]
        #
        rMat_d = xf.makeDetectorRotMat(icfg['detector']['transform']['tilt_angles'])
        tVec_d = np.vstack(icfg['detector']['transform']['t_vec_d'])
        #
        chi = icfg['oscillation_stage']['chi']
        tVec_s = np.vstack(icfg['oscillation_stage']['t_vec_s'])
        #
        oes = np.zeros(oem.dataStore.shape)
        for i_grn in range(ngrains):
            spots_table = np.loadtxt(os.path.join(analysis_dir, 'spots_%05d.out' %i_grn))
            idx_m = spots_table[:, 0] >= 0
            for i_map in range(nmaps):
                idx_g = spots_table[:, 1] == gvids[i_map]
                idx = np.logical_and(idx_m, idx_g)
                nrefl = sum(idx)
                
                omes_fit = xf.mapAngle(spots_table[idx, -1], np.radians(cfg.find_orientations.omega.period), units='radians')
                xy_det = spots_table[idx, -3:]
                xy_det[:, 2] = np.zeros(nrefl)
                
                rMat_s_array = xfcapi.makeOscillRotMatArray(chi, omes_fit)
                
                # form in-plane vectors for detector points list in DETECTOR FRAME
                P2_d = xy_det.T
                
                # in LAB FRAME
                P2_l = np.dot(rMat_d, P2_d) + tVec_d # point on detector
                P0_l = np.hstack(
                    [tVec_s + np.dot(rMat_s_array[j], tVec_c[i_grn, :].reshape(3, 1)) for j in range(nrefl)]
                ) # origin of CRYSTAL FRAME

                # diffraction unit vector components in LAB FRAME
                dHat_l = unitVector(P2_l - P0_l)
                P2_l = np.dot(rMat_d, xy_det.T) + tVec_d
                
                # angles for reference frame
                dHat_ref_l = unitVector(P2_l)
    
                # append etas and omes
                etas_fit = np.arctan2(dHat_ref_l[1, :], dHat_ref_l[0, :]).flatten()
           
                # find indices, then truncate or wrap
                i_ome = cellIndices(oem.omeEdges, omes_fit)
                if full_ome_range:
                    i_ome[i_ome < 0] = np.mod(i_ome, nome) + 1
                    i_ome[i_ome >= nome] = np.mod(i_ome, nome)
                else:
                    incl = np.logical_or(i_ome >= 0, i_ome < nome)
                    i_ome = i_ome[incl]
                j_eta = cellIndices(oem.etaEdges, etas_fit)
                if full_eta_range:
                    j_eta[j_eta < 0] = np.mod(j_eta, neta) + 1
                    j_eta[j_eta >= neta] = np.mod(j_eta, neta)
                else:
                    incl = np.logical_or(j_eta >= 0, j_eta < neta)
                    j_eta = j_eta[incl]

                if np.max(i_ome) >= nome or np.min(i_ome) < 0 or np.max(j_eta) >= neta or np.min(j_eta) < 0:
                    import pdb; pdb.set_trace()
                # add to map
                oes[i_map][i_ome, j_eta] = 1
            pass
        pass
    
    # simulate quaternion points
    if not plot_from_grains:
        oes = simulateOmeEtaMaps(omeEdges, etaEdges, pd,
                                 expMaps,
                                 chi=chi,
                                 etaTol=0.01, omeTol=0.01,
                                 etaRanges=None, omeRanges=None,
                                 bVec=xf.bVec_ref, eVec=xf.eta_ref, vInv=xf.vInv_ref)
    
    # tick labling
    omes = np.degrees(oem.omeEdges)
    etas = np.degrees(oem.etaEdges)
    num_ticks = 7
    xmin = np.amin(etas); xmax = np.amax(etas)
    dx = (xmax - xmin) / (num_ticks - 1.); dx1 = (len(etas) - 1) / (num_ticks - 1.)
    xtlab = ["%.0f" % (xmin + i*dx) for i in range(num_ticks)]
    xtloc = np.array([i*dx1 for i in range(num_ticks)]) - 0.5
    ymin = np.amin(omes); ymax = np.amax(omes)
    dy = (ymax - ymin) / (num_ticks - 1.); dy1 = (len(omes) - 1) / (num_ticks - 1.)
    ytlab = ["%.0f" % (ymin + i*dy) for i in range(num_ticks)]
    ytloc = np.array([i*dy1 for i in range(num_ticks)]) - 0.5
    
    # Plot the three kernel density estimates
    n_maps = len(oem.iHKLList)
    
    fig_list =[plt.figure(num=i+1) for i in range(n_maps)]
    ax_list = [fig_list[i].gca() for i in range(n_maps)]
    for i_map in range(n_maps):
        y, x = np.where(oes[i_map] > 0)
        ax_list[i_map].hold(True)
        ax_list[i_map].imshow(oem.dataStore[i_map] > 0.1, cmap=cm.bone)
        ax_list[i_map].set_title(r'Map for $\{%d %d %d\}$' %tuple(pd.hkls[:, i_map]))
        ax_list[i_map].set_xlabel(r'Azimuth channel, $\eta$; $\Delta\eta=%.3f$' %delta_ome)
        ax_list[i_map].set_ylabel(r'Rotation channel, $\omega$; $\Delta\omega=%.3f$' %delta_ome)
        ax_list[i_map].plot(x, y, 'c+')
        ax_list[i_map].xaxis.set_ticks(xtloc)
        ax_list[i_map].xaxis.set_ticklabels(xtlab)
        ax_list[i_map].yaxis.set_ticks(ytloc)
        ax_list[i_map].yaxis.set_ticklabels(ytlab)
        ax_list[i_map].axis('tight')
    plt.show()
    return fig_list, oes
    
if __name__ == '__main__':
    """
    USAGE : python check_indexing <cfg_file> 
    """
    parser = argparse.ArgumentParser(description='Make median dark from cfg file')

    parser.add_argument('cfg', metavar='cfg_filename', type=str, help='a YAML config filename')
    parser.add_argument('-t','--show-trials', help='plot trial orientations', action='store_true', default=False)
    parser.add_argument('-g','--plot-from-grains', help='plot fit orientations', action='store_true', default=False)

    args = vars(parser.parse_args(sys.argv[1:]))
    
    dark = check_indexing_plots(args['cfg'], plot_trials=args['show_trials'], plot_from_grains=args['plot_from_grains'])
   
