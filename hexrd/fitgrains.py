import argparse
import logging
import multiprocessing as mp
from multiprocessing.queues import Empty
import os
import sys
import textwrap
import time

import yaml

import numpy as np
from scipy.sparse import coo_matrix

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except IOError:
    pass

from hexrd.coreutil import initialize_experiment, migrate_detector_config

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd import rotations as rot
from hexrd.xrd.xrdutil import pullSpots


logger = logging.getLogger(__name__)


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

    n_frames = reader.getNFrames()
    logger.info("reading %d frames of data", n_frames)
    if have_progBar:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=n_frames).start()

    frame_list = []
    for i in range(n_frames):
        frame = reader.read()
        frame[frame <= cfg.fit_grains.threshold] = 0
        frame_list.append(coo_matrix(frame))
        if have_progBar:
            pbar.update(i)
    frame_list = np.array(frame_list)
    omega_start = np.radians(cfg.image_series.omega.start)
    omega_step = np.radians(cfg.image_series.omega.step)
    reader = [frame_list, [omega_start, omega_step]]
    if have_progBar:
        pbar.finish()

    logger.info("pulling spots for %d orientations", n_quats)
    if have_progBar:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=n_quats).start()

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
    eta_range = np.radians(cfg.find_orientations.eta.range)
    omega_period = np.radians(cfg.find_orientations.omega.period)
    tth_tol = cfg.fit_grains.tolerance.tth
    eta_tol = cfg.fit_grains.tolerance.eta
    omega_tol = cfg.fit_grains.tolerance.omega
    panel_buffer = cfg.fit_grains.panel_buffer
    npdiv = cfg.fit_grains.npdiv
    threshold = cfg.fit_grains.threshold

    ncpus = cfg.multiprocessing
    logging.info('running pullspots with %d processors')
    for i in range(ncpus):
        w = PullSpotsWorker(
            job_queue,
            results,
            reader, pd, detector_params, distortion,
            eta_range,
            omega_period,
            tth_tol,
            eta_tol,
            omega_tol,
            panel_buffer,
            npdiv,
            threshold,
            os.path.join(cfg.analysis_dir, 'spots_%05d.out')
            )
        w.daemon = True
        w.start()
    while True:
        n_res = len(results)
        if have_progBar:
            pbar.update(n_res)
        if n_res == n_quats:
            break
    job_queue.join()
    if have_progBar:
        pbar.finish()

# bMat = pd.latVecOps['B']
# wlen = pd.wavelength


class PullSpotsWorker(mp.Process):


    def __init__(self,
                 jobs, results,
                 reader, plane_data, detector_params, distortion, eta_range,
                 ome_period, tth_tol, eta_tol, ome_tol, panel_buff, npdiv,
                 pthresh, spots_f
                 ):
        super(PullSpotsWorker, self).__init__()
        self._jobs = jobs
        self._results = results
        self._reader = reader
        self._plane_data = plane_data
        self._detector_params = detector_params
        self._distortion = distortion
        self._eta_range = eta_range
        self._ome_period = ome_period
        self._tth_tol = tth_tol
        self._eta_tol = eta_tol
        self._ome_tol = ome_tol
        self._panel_buff = panel_buff
        self._npdiv = npdiv
        self._pthresh = pthresh
        self._spots_f = spots_f


    def pull_spots(self, grain_id, grain_params, iteration):
        return pullSpots(
            self._plane_data,
            self._detector_params,
            grain_params,
            self._reader,
            distortion=self._distortion,
            eta_range=self._eta_range,
            ome_period=self._ome_period,
            tth_tol=self._tth_tol[iteration],
            eta_tol=self._eta_tol[iteration],
            ome_tol=self._ome_tol[iteration],
            panel_buff=self._panel_buff,
            npdiv=self._npdiv,
            threshold=self._pthresh,
            filename=self._spots_f % grain_id,
            )


    def fit_grains(self):
        pass


    def run(self):
        while True:
            try:
                grain_id, grain_params = self._jobs.get(False)

                for iteration in range(1):
                    res = self.pull_spots(grain_id, grain_params, iteration)
                    # self.fit_grains()


                self._results.append((grain_id, res))
                self._jobs.task_done()
            except Empty:
                break
