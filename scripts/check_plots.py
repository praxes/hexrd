import sys, os

import argparse

import yaml
import cPickle

import numpy as np

from hexrd import config
from hexrd import matrixutil as mutil

from hexrd.fitgrains import get_instrument_parameters

from hexrd.matrixutil import unitVector

from hexrd.xrd import transforms as xf
from hexrd.xrd import rotations as rot
from hexrd.xrd.symmetry import toFundamentalRegion
from hexrd.xrd.xrdutil import simulateOmeEtaMaps

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
plt.ioff()

def check_indexing_plots(cfg_filename, plot_trials=False):
    cfg = config.open(cfg_filename)[0]    # use first block, like indexing
    icfg = get_instrument_parameters(cfg)

    working_dir = cfg.working_dir

    if plot_trials:
        compl = np.loadtxt(os.path.join(working_dir, 'completeness.dat'))
        tq = np.atleast_2d(np.loadtxt(os.path.join(working_dir, 'trial_orientations.dat')))
        quats = tq[compl >= 0., :]
    else:
        quats = np.atleast_2d(np.loadtxt(os.path.join(working_dir, 'accepted_orientations.dat')))
    expMaps = np.tile(2. * np.arccos(quats[:, 0]), (3, 1))*unitVector(quats[:, 1:].T)
    
    # load maps that were used
    oem = cPickle.load(
        open(cfg.find_orientations.orientation_maps.file, 'r')
        )
    omeEdges = np.degrees(oem.omeEdges)
    etaEdges = np.degrees(oem.etaEdges)
    planeData = oem.planeData
    
    delta_ome = abs(omeEdges[1]-omeEdges[0])
    
    chi = icfg['oscillation_stage']['chi']
    
    # simulate quaternion points
    oes = simulateOmeEtaMaps(omeEdges, etaEdges, planeData,
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
        ax_list[i_map].set_title(r'Map for $\{%d %d %d\}$' %tuple(planeData.hkls[:, i_map]))
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
    parser.add_argument('-t','--show-trials', help='plot trial orientations', nargs='?', type=bool, default=False)

    args = vars(parser.parse_args(sys.argv[1:]))
    
    dark = check_indexing_plots(args['cfg'], plot_trials=args['show_trials'])
   