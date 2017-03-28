# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import argparse, os, sys

import cPickle

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

from hexrd import config
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.coreutil import get_instrument_parameters

from sklearn.cluster import dbscan

from scipy import cluster

def adist(ang0, ang1):
    resd = xfcapi.angularDifference(ang0 ,ang1)
    return np.sqrt(sum(resd**2))

def build_overlap_table(cfg, tol_mult=0.5):
    
    icfg = get_instrument_parameters(cfg)
    
    gt = np.loadtxt(
        os.path.join(cfg.analysis_dir, 'grains.out')
    )
    
    ngrains = len(gt)
    
    mat_list = cPickle.load(open(cfg.material.definitions, 'r'))
    mat_names = [mat_list[i].name for i in range(len(mat_list))]
    mat_dict = dict(zip(mat_names, mat_list))
    
    matl = mat_dict[cfg.material.active]
    
    pd = matl.planeData
    pd.exclusions = np.zeros(len(pd.exclusions), dtype=bool)
    pd.tThMax = np.radians(cfg.fit_grains.tth_max)
    pd.tThWidth = np.radians(cfg.fit_grains.tolerance.tth[-1])
    
    # for clustering...
    eps = tol_mult*np.radians(
        min(
            min(cfg.fit_grains.tolerance.eta), 
            2*min(cfg.fit_grains.tolerance.omega)
        )
    )

    # merged two-theta indices
    tth_ranges_merged = pd.getMergedRanges()[0]
    pids = []
    for hklids in tth_ranges_merged:
        pids.append(
            [pd.hklDataList[hklids[i]]['hklID'] for i in range(len(hklids))]
        )
        
    # Make table of unit diffraction vectors
    st = []
    for i in range(ngrains):
        this_st = np.loadtxt(
            os.path.join(cfg.analysis_dir, 'spots_%05d.out' %i)
            )
        #... do all predicted?
        valid_spt = this_st[:, 0] >= 0
        #valid_spt = np.ones(len(this_st), dtype=bool)

        angs = this_st[valid_spt, 7:10]

        dvec = xfcapi.anglesToDVec(
            angs, 
            chi=icfg['oscillation_stage']['chi']
        )

        # [ grainID, reflID, hklID, D_s[0], D_s[1], D_s[2], tth, eta, ome ]
        st.append(
            np.hstack([
                i*np.ones((sum(valid_spt), 1)),
                this_st[valid_spt, :2], 
                dvec, 
                angs,
            ])
        )

    # make overlap table
    # [[range_0], [range_1], ..., [range_n]]
    # range_0 = [grainIDs, reflIDs, hklIDs] that are within tol
    overlap_table = []
    ii = 0
    for pid in pids:
        tmp = []; a = []; b = []; c = []
        for j in range(len(pid)):
            a.append(
                np.vstack(
                    [st[i][st[i][:, 2] == pid[j], 3:6] for i in range(len(st))]
                )
            )
            b.append(
                np.vstack(
                    [st[i][st[i][:, 2] == pid[j], 0:3] for i in range(len(st))]
                )
            )
            c.append(
                np.vstack(
                    [st[i][st[i][:, 2] == pid[j], 6:9] for i in range(len(st))]
                )
            )
            pass
        a = np.vstack(a)
        b = np.vstack(b)
        c = np.vstack(c)    
        if len(a) > 0:
            # run dbscan
            core_samples, labels = dbscan(
                a,
                eps=eps,
                min_samples=2,
                metric='minkowski', p=2,
            )
            
            cl = np.array(labels, dtype=int) # convert to array
            noise_points = cl == -1 # index for marking noise
            cl += 1 # move index to 1-based instead of 0
            cl[noise_points] = -1 # re-mark noise as -1
            
            # extract number of clusters
            if np.any(cl == -1):
                nblobs = len(np.unique(cl)) - 1
            else:
                nblobs = len(np.unique(cl))
            
            for i in range(1, nblobs+1):
                # put in check on omega here
                these_angs = c[np.where(cl == i)[0], :]
                local_cl = cluster.hierarchy.fclusterdata(
                    these_angs[:, 1:],
                    eps,
                    criterion='distance',
                    metric=adist
                    )
                local_nblobs = len(np.unique(local_cl))
                if local_nblobs < len(these_angs):
                    for j in range(1, local_nblobs + 1):
                        npts = sum(local_cl == j)
                        if npts >= 2:
                            cl_idx = np.where(local_cl == j)[0]
                            #import pdb; pdb.set_trace()
                            tmp.append(
                                b[np.where(cl == i)[0][cl_idx], :]
                            )
        print "processing ring set %d" %ii
        ii += 1
        overlap_table.append(tmp)
    return overlap_table
    
def build_discrete_cmap(ngrains):
    
    # define the colormap
    cmap = plt.cm.jet
    
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    
    # define the bins and normalize
    bounds = np.linspace(0, ngrains, ngrains+1)
    norm = BoundaryNorm(bounds, cmap.N)

    return cmap, norm
    
#%%
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Make overlap table from cfg file')
    parser.add_argument(
        'cfg', metavar='cfg_filename', 
        type=str, help='a YAML config filename')
    parser.add_argument(
        '-m', '--multiplier', 
        help='multiplier on angular tolerance', 
        type=float, default=0.5)

    args = vars(parser.parse_args(sys.argv[1:]))
    
    cfg = config.open(args['cfg'])[0]
    print "loaded config file %s" %args['cfg']
    overlap_table = build_overlap_table(cfg)
    np.savez(os.path.join(cfg.analysis_dir, 'overlap_table.npz'), 
             *overlap_table)
#%%
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
# 
#etas = np.radians(np.linspace(0, 359, num=360))
#cx = np.cos(etas)
#cy = np.sin(etas)
#cz = np.zeros_like(etas)
#
#ax.plot(cx, cy, cz, c='b')
#ax.plot(cx, cz, cy, c='g')
#ax.plot(cz, cx, cy, c='r')
#ax.scatter3D(a[:, 0], a[:, 1], a[:, 2], c=b[:, 0], cmap=cmap, norm=norm, marker='o', s=20)
#
#ax.set_xlabel(r'$\mathbf{\mathrm{X}}_s$')
#ax.set_ylabel(r'$\mathbf{\mathrm{Y}}_s$')
#ax.set_zlabel(r'$\mathbf{\mathrm{Z}}_s$')
#
#ax.elev = 124
#ax.azim = -90
#
#ax.axis('equal')
#
##fname = "overlaps_%03d.png"
##for i in range(360):
##    ax.azim += i
##    fig.savefig(
##        fname %i, dpi=200, facecolor='w', edgecolor='w',
##        orientation='landcape')
