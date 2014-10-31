from __future__ import absolute_import

import copy
import logging
import multiprocessing as mp
from multiprocessing.queues import Empty
import os
import sys
import time

import yaml

import numpy as np
from scipy.sparse import coo_matrix

from hexrd.coreutil import initialize_experiment, migrate_detector_config
from hexrd.utils.progressbar import (
    Bar, ETA, Percentage, ProgressBar, ReverseBar
    )

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd.transforms import bVec_ref, eta_ref, vInv_ref
from hexrd.xrd import rotations as rot
from hexrd.xrd.xrdutil import pullSpots


logger = logging.getLogger(__name__)


# grain parameter refinement flags
gFlag = np.array([1, 1, 1,
                  1, 1, 1,
                  1, 1, 1, 1, 1, 1], dtype=bool)
# grain parameter scalings
gScl  = np.array([1., 1., 1.,
                  1., 1., 1.,
                  1., 1., 1., 0.01, 0.01, 0.01])


def read_frames(reader, cfg):
    start = time.time()

    n_frames = reader.getNFrames()
    logger.info("reading %d frames of data", n_frames)
    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=n_frames).start()

    frame_list = []
    for i in range(n_frames):
        frame = reader.read()
        frame[frame <= cfg.fit_grains.threshold] = 0
        frame_list.append(coo_matrix(frame))
        pbar.update(i)
    frame_list = np.array(frame_list)
    omega_start = np.radians(cfg.image_series.omega.start)
    omega_step = np.radians(cfg.image_series.omega.step)
    reader = [frame_list, [omega_start, omega_step]]
    pbar.finish()
    logger.info('read %d frames in %.2f seconds', n_frames, time.time()-start)
    return reader


def fit_grain():
    pass



def fit_grains(cfg, force=False):

    pd, reader, detector = initialize_experiment(cfg)

    # attempt to load the new detector parameter file
    det_p = cfg.detector.parameters
    if not os.path.exists(det_p):
        migrate_detector_config(
            np.loadtxt(cfg.detector.parameters_old),
            cfg.detector.pixels.rows,
            cfg.detector.pixels.columns,
            cfg.detector.pixels.size,
            detID='GE',
            chi=0.,
            tVec_s=np.zeros(3),
            filename=cfg.detector.parameters
            )

    with open(det_p, 'r') as f:
        # only one panel for now
        # TODO: configurize this
        instr_cfg = [instr_cfg for instr_cfg in yaml.load_all(f)][0]
    detector_params = np.hstack([
        instr_cfg['detector']['transform']['tilt_angles'],
        instr_cfg['detector']['transform']['t_vec_d'],
        instr_cfg['oscillation_stage']['chi'],
        instr_cfg['oscillation_stage']['t_vec_s'],
        ])
    # ***FIX***
    # at this point we know we have a GE and hardwire the distortion func;
    # need to pull name from yml file in general case
    distortion = (
        dFuncs.GE_41RT, instr_cfg['detector']['distortion']['parameters']
        )

    tth_max = cfg.fit_grains.tth_max
    if tth_max is True:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() > detector.getTThMax()
    elif tth_max > 0:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() >= np.radians(tth_max)

    # load quaternion file
    quats = np.atleast_2d(
        np.loadtxt(os.path.join(cfg.analysis_dir, 'quats.out'))
        )
    n_quats = len(quats)
    quats = quats.T

    reader = read_frames(reader, cfg)

    logger.info("pulling spots for %d orientations", n_quats)
    pbar = ProgressBar(
        widgets=[Bar('>'), ' ', ETA(), ' ', ReverseBar('<')],
        maxval=n_quats
        ).start()

    # create the job queue
    job_queue = mp.JoinableQueue()
    manager = mp.Manager()
    results = manager.list()

    # load the queue
    phi, n = rot.angleAxisOfRotMat(rot.rotMatOfQuat(quats))
    for i, quat in enumerate(quats.T):
        exp_map = phi[i]*n[:, i]
        grain_params = np.hstack(
            [exp_map.flatten(), 0., 0., 0., 1., 1., 1., 0., 0., 0.]
            )
        job_queue.put((i, grain_params))

    # don't query these in the loop, will spam the logger:
    pkwargs = {
        'distortion': distortion,
        'eta_range': np.radians(cfg.find_orientations.eta.range),
        'ome_period': np.radians(cfg.find_orientations.omega.period),
        'tth_tol': cfg.fit_grains.tolerance.tth,
        'eta_tol': cfg.fit_grains.tolerance.eta,
        'ome_tol': cfg.fit_grains.tolerance.omega,
        'panel_buff': cfg.fit_grains.panel_buffer,
        'npdiv': cfg.fit_grains.npdiv,
        'threshold': cfg.fit_grains.threshold,
        }
    spots_stem = os.path.join(cfg.analysis_dir, 'spots_%05d.out')

    ncpus = cfg.multiprocessing
    logging.info('running pullspots with %d processors')
    for i in range(ncpus):
        w = PullSpotsWorker(
            job_queue, results, reader, pd, detector_params, spots_stem,
            copy.copy(pkwargs)
            )
        w.daemon = True
        w.start()
    while True:
        n_res = len(results)
        pbar.update(n_res)
        if n_res == n_quats:
            break
    job_queue.join()

    pbar.finish()

    # bMat = pd.latVecOps['B']
    # wlen = pd.wavelength

    #print >> grains_file, \
    #  "%d\t%1.7e\t%1.7e\t"                         % (grainID, sum(idx)/float(len(idx)), sum(resd_f2**2)) + \
    #  "%1.7e\t%1.7e\t%1.7e\t"                      % tuple(g_refined[:3]) + \
    #  "%1.7e\t%1.7e\t%1.7e\t"                      % tuple(g_refined[3:6]) + \
    #  "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t" % tuple(g_refined[6:]) + \
    #  "%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e\t%1.7e"   % (eMat[0, 0], eMat[1, 1], eMat[2, 2], eMat[1, 2], eMat[0, 2], eMat[0, 1])
    #elapsed = (time.clock() - start)
    #print "grain %d took %.2f seconds" %(grainID, elapsed)




class PullSpotsWorker(mp.Process):


    def __init__(
        self, jobs, results, reader, plane_data, det_pars, spots_stem, pkwargs
        ):
        super(PullSpotsWorker, self).__init__()
        self._jobs = jobs
        self._results = results
        self._reader = reader
        self._plane_data = plane_data
        self._det_pars = det_pars
        self._spots_stem = spots_stem
        self._eta_tol = pkwargs.pop('eta_tol')
        self._ome_tol = pkwargs.pop('ome_tol')
        self._tth_tol = pkwargs.pop('tth_tol')
        # a dict containing the pullSpots kwargs
        self._pkwargs = pkwargs


    def pull_spots(self, grain_id, grain_params, iteration):
        return pullSpots(
            self._plane_data,
            self._det_pars,
            grain_params,
            self._reader,
            filename=self._spots_stem % grain_id,
            eta_tol=self._eta_tol[iteration],
            ome_tol=self._ome_tol[iteration],
            tth_tol=self._tth_tol[iteration],
            **self._pkwargs
            )


    def fit_grains(self, grain_id, grain_params):
        gtable = np.loadtxt(self._spots_f % grain_id)
        valid_refl_ids = gtable[:, 0] >= 0
        pred_ome = gtable[:, 6]
        if np.sign(ome_delta) < 0:
            idx_ome = np.logical_and(
                pred_ome < d2r*(ome_start + 2*ome_delta),
                pred_ome > d2r*(ome_stop  - 2*ome_delta)
                )
        else:
            idx_ome = np.logical_and(
                pred_ome > d2r*(ome_start + 2*ome_delta),
                pred_ome < d2r*(ome_stop  - 2*ome_delta)
                )

        idx = np.logical_and(idx0, idx_ome)
        hkls = gtable[idx, 1:4].T # must be column vectors
        xyo_det = gtable[idx, -3:] # these are the cartesian centroids + ome
        xyo_det[:, 2] = xf.mapAngle(xyo_det[:, 2], ome_period)
        if sum(idx) <= 12:
            return grain_params
        return fitting.fitGrain(
            xyo_det, hkls, bMat, wlen,
            detector_params,
            grain_params[:3], grain_params[3:6], grain_params[6:],
            beamVec=bVec_ref, etaVec=eta_ref,
            distortion=(dFunc, dParams),
            gFlag=gFlag, gScl=gScl,
            omePeriod=ome_period
            )


    def run(self):
        while True:
            try:
                grain_id, grain_params = self._jobs.get(False)

                for iteration in range(1): # TODO: change to 2
                    self.pull_spots(grain_id, grain_params, iteration)
                    #temp = self.fit_grains()
                    #if temp is grain_params:
                    #    break
                    #grain_params = temp

#                eMat = logm(
#                    np.linalg.inv(mutil.vecMVToSymm(refined_grain_params[6:]))
#                    )
#
#                resd_f2 = fitting.objFuncFitGrain(
#                    grain_params[gFlag], grain_params, gFlag,
#                    detector_params,
#                    xyo_det, hkls, bMat, wlen,
#                    bVec_ref, eta_ref,
#                    dFunc, dParams,
#                    ome_period,
#                    simOnly=False
#                    )


                self._results.append((grain_id, res))
                self._jobs.task_done()
            except Empty:
                break
