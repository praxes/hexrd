# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import cPickle
import yaml
import time
import sys

import numpy as np
import numba

from hexrd import gridutil as gutil
from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd.xrdutil import angularPixelSize
from hexrd.utils import profiler

import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt

# load in config files for convenience
instr_cfg = yaml.load(open('./analysis/ge_detector.yml','r'))
pixel_pitch = instr_cfg['detector']['pixels']['size']

# load material file
mat_list = cPickle.load(open('./include/materials.cpl', 'r'))
pd = mat_list[-1].planeData

# define instrument params
rMat_d = xfcapi.makeDetectorRotMat(instr_cfg['detector']['transform']['tilt_angles'])
tVec_d = np.r_[instr_cfg['detector']['transform']['t_vec_d']]
rMat_s = np.eye(3)
tVec_s = np.zeros(3)
rMat_c = np.eye(3)
tVec_c = np.zeros(3)

try:
    #dfunc = instr_cfg['detector']['distortion']['function_name']
    dparams = instr_cfg['detector']['distortion']['parameters']
    distortion = (xf.dFunc_ref, dparams)
except (KeyError):
    distortion=None

# for defining patches
delta_eta = 0.5
neta = int(360/float(delta_eta))

eta = np.radians(delta_eta*np.linspace(0, neta-1, num=neta))

angs = [np.vstack([i*np.ones(neta), eta, np.zeros(neta)]) for i in pd.getTTh()]

# need xy coords and pixel sizes
gVec_ring_l = xf.anglesToGVec(angs[0].T, xf.bVec_ref, xf.eta_ref)
xydet_ring = xfcapi.gvecToDetectorXY(gVec_ring_l.T,
                                     rMat_d, rMat_s, rMat_c,
                                     tVec_d, tVec_s, tVec_c)

if distortion is not None:
    det_xy = distortion[0](xydet_ring,
                           distortion[1],
                           invert=True)
ang_ps = angularPixelSize(det_xy, pixel_pitch,
                          rMat_d, rMat_s,
                          tVec_d, tVec_s, tVec_c,
                          distortion=distortion)

def compute_areas(xy_eval_vtx, conn):
    areas = np.zeros(len(conn))
    for i in range(len(conn)):
        polygon = [[xy_eval_vtx[conn[i, j], 0],
                    xy_eval_vtx[conn[i, j], 1]] for j in range(4)]
        areas[i] = gutil.computeArea(polygon)
    return areas

@numba.jit
def compute_areas_2(xy_eval_vtx, conn):
    areas = np.empty(len(conn))
    for i in range(len(conn)):
        c0, c1, c2, c3 = conn[i]
        vtx0x, vtx0y = xy_eval_vtx[conn[i,0]]
        vtx1x, vtx1y = xy_eval_vtx[conn[i,1]]
        v0x, v0y = vtx1x-vtx0x, vtx1y-vtx0y
        acc = 0
        for j in range(2, 4):
            vtx_x, vtx_y = xy_eval_vtx[conn[i,j]]
            v1x = vtx_x - vtx0x
            v1y = vtx_y - vtx0y
            acc += v0x*v1y - v1x*v0y

        areas[i] = 0.5 * acc
    return areas

@numba.jit
def compute_areas_3(xy_eval_vtx, conn):
    areas = np.empty(len(conn))
    for i in range(len(conn)):
        c0, c1, c2, c3 = conn[i]
        vtx0x, vtx0y = xy_eval_vtx[conn[i,0]]
        vtx1x, vtx1y = xy_eval_vtx[conn[i,1]]
        vtx2x, vtx2y = xy_eval_vtx[conn[i,2]]
        vtx3x, vtx3y = xy_eval_vtx[conn[i,3]]
        v0x, v0y = vtx1x-vtx0x, vtx1y-vtx0y
        v1x, v1y = vtx2x-vtx0x, vtx2y-vtx0y
        v2x, v2y = vtx3x-vtx0x, vtx3y-vtx0y

        areas[i] = 0.5 * (v0x*v1y - v1x*v0y + v1x*v2y - v2x*v1y)

    return areas


        
def make_reflection_patches(instr_cfg, tth_eta, ang_pixel_size,
                            omega=None,
                            tth_tol=0.2, eta_tol=1.0,
                            rMat_c=np.eye(3), tVec_c=np.c_[0.,0.,0.].T,
                            distortion=distortion,
                            npdiv=1, quiet=False, compute_areas_func=compute_areas):
    """
    prototype function for making angular patches on a detector

    panel_dims are [(xmin, ymin), (xmax, ymax)] in mm

    pixel_pitch is [row_size, column_size] in mm

    DISTORTION HANDING IS STILL A KLUDGE

    patches are:

                 delta tth
   d  ------------- ... -------------
   e  | x | x | x | ... | x | x | x |
   l  ------------- ... -------------
   t                 .
   a                 .
                     .
   e  ------------- ... -------------
   t  | x | x | x | ... | x | x | x |
   a  ------------- ... -------------

    """
    npts = len(tth_eta)

    # detector frame
    rMat_d = xfcapi.makeDetectorRotMat(
        instr_cfg['detector']['transform']['tilt_angles']
        )
    tVec_d = np.r_[instr_cfg['detector']['transform']['t_vec_d']]
    pixel_size = instr_cfg['detector']['pixels']['size']

    frame_nrows = instr_cfg['detector']['pixels']['rows']
    frame_ncols = instr_cfg['detector']['pixels']['columns']

    panel_dims = (-0.5*np.r_[frame_ncols*pixel_size[1],
                             frame_nrows*pixel_size[0]],
                  0.5*np.r_[frame_ncols*pixel_size[1],
                            frame_nrows*pixel_size[0]])
    row_edges = np.arange(frame_nrows + 1)[::-1]*pixel_size[1] + panel_dims[0][1]
    col_edges = np.arange(frame_ncols + 1)*pixel_size[0] + panel_dims[0][0]

    # sample frame
    chi    = instr_cfg['oscillation_stage']['chi']
    tVec_s = np.r_[instr_cfg['oscillation_stage']['t_vec_s']]

    # data to loop
    # ...WOULD IT BE CHEAPER TO CARRY ZEROS OR USE CONDITIONAL?
    if omega is None:
        full_angs = np.hstack([tth_eta, np.zeros((npts, 1))])
    else:
        full_angs = np.hstack([tth_eta, omega.reshape(npts, 1)])
    patches = []
    for angs, pix in zip(full_angs, ang_pixel_size):
        # need to get angular pixel size
        rMat_s = xfcapi.makeOscillRotMat([chi, angs[2]])

        ndiv_tth = npdiv*np.ceil( tth_tol/np.degrees(pix[0]) )
        ndiv_eta = npdiv*np.ceil( eta_tol/np.degrees(pix[1]) )

        tth_del = np.arange(0, ndiv_tth+1)*tth_tol/float(ndiv_tth) - 0.5*tth_tol
        eta_del = np.arange(0, ndiv_eta+1)*eta_tol/float(ndiv_eta) - 0.5*eta_tol

        # store dimensions for convenience
        #   * etas and tths are bin vertices, ome is already centers
        sdims = [ len(eta_del)-1, len(tth_del)-1 ]

        # meshgrid args are (cols, rows), a.k.a (fast, slow)
        m_tth, m_eta = np.meshgrid(tth_del, eta_del)
        npts_patch   = m_tth.size

        # calculate the patch XY coords from the (tth, eta) angles
        # * will CHEAT and ignore the small perturbation the different
        #   omega angle values causes and simply use the central value
        gVec_angs_vtx = np.tile(angs, (npts_patch, 1)) \
                        + np.radians(
                            np.vstack([m_tth.flatten(),
                                       m_eta.flatten(),
                                       np.zeros(npts_patch)
                                       ]).T
                                     )

        # will need this later
        rMat_s = xfcapi.makeOscillRotMat([chi, angs[2]])

        # FOR ANGULAR MESH
        conn = gutil.cellConnectivity( sdims[0], sdims[1], origin='ll')
        gVec_c = xf.anglesToGVec(gVec_angs_vtx,
                                 xf.bVec_ref, xf.eta_ref,
                                 rMat_s=rMat_s,
                                 rMat_c=rMat_c)

        xy_eval_vtx = xfcapi.gvecToDetectorXY(gVec_c.T,
                                              rMat_d, rMat_s, rMat_c,
                                              tVec_d, tVec_s, tVec_c)
        if distortion is not None and len(distortion) == 2:
            xy_eval_vtx = distortion[0](xy_eval_vtx, distortion[1], invert=True)
            pass

        areas = compute_areas_func(xy_eval_vtx, conn)
        
        # EVALUATION POINTS
        #   * for lack of a better option will use centroids
        tth_eta_cen = gutil.cellCentroids( np.atleast_2d(gVec_angs_vtx[:, :2]), conn )
        gVec_angs  = np.hstack([tth_eta_cen,
                                 np.tile(angs[2], (len(tth_eta_cen), 1))])
        gVec_c = xf.anglesToGVec(gVec_angs,
                                 xf.bVec_ref, xf.eta_ref,
                                 rMat_s=rMat_s,
                                 rMat_c=rMat_c)

        xy_eval = xfcapi.gvecToDetectorXY(gVec_c.T,
                                          rMat_d, rMat_s, rMat_c,
                                          tVec_d, tVec_s, tVec_c)
        if distortion is not None and len(distortion) == 2:
            xy_eval = distortion[0](xy_eval, distortion[1], invert=True)
            pass
        row_indices   = gutil.cellIndices(row_edges, xy_eval[:, 1])
        col_indices   = gutil.cellIndices(col_edges, xy_eval[:, 0])

        patches.append(((gVec_angs_vtx[:, 0].reshape(m_tth.shape),
                         gVec_angs_vtx[:, 1].reshape(m_tth.shape)),
                        (xy_eval_vtx[:, 0].reshape(m_tth.shape),
                         xy_eval_vtx[:, 1].reshape(m_tth.shape)),
                        conn,
                        areas.reshape(sdims[0], sdims[1]),
                        (row_indices.reshape(sdims[0], sdims[1]),
                         col_indices.reshape(sdims[0], sdims[1]))
                        )
                    )
        pass
    return patches


"""
#######################
TESTING
"""
# fig = plt.figure()
# ax = fig.add_subplot(111, aspect='equal')
# ax.hold(True)
def run_with_impl(compute_areas_impl):
    start = time.clock()                      # time this
    pts = []
    for i_ring in range(len(pd.getTTh())):
        pts.append(make_reflection_patches(instr_cfg, angs[0].T[:, :2], ang_ps,
                                           omega=None,
                                           tth_tol=0.1, eta_tol=delta_eta,
                                           rMat_c=np.eye(3), tVec_c=np.c_[0.,0.,0.].T,
                                           distortion=distortion,
                                           npdiv=2, quiet=False,
                                           compute_areas_func=compute_areas_impl))
    elapsed = (time.clock() - start)
    print "make_reflection_patches (%s) on %d patches on %d rings: %f" \
        % (compute_areas_impl.__name__, neta, len(pd.getTTh()), elapsed)

def main(args):
    # if there are arguments, try to load them as profile config
    if args:
        profiler.instrument_all(args)

    run_with_impl(compute_areas_3)
    run_with_impl(compute_areas_2)
    run_with_impl(compute_areas_3)
    run_with_impl(compute_areas_2)
    run_with_impl(compute_areas_3)
    run_with_impl(compute_areas_2)
    run_with_impl(compute_areas)

    if args:
        profiler.dump_results(args)
#     for data in pts[i_ring]:
#         patch = data[-1]
#         ax.plot(patch[0], patch[1], 'r.')
# fig.show()
if __name__ == '__main__':
    main(sys.argv[1:])
