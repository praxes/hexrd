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
from scipy.linalg.matfuncs import logm

from hexrd.coreutil import initialize_experiment, migrate_detector_config
from hexrd.matrixutil import vecMVToSymm
from hexrd.utils.progressbar import (
    Bar, ETA, Percentage, ProgressBar, ReverseBar
    )

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd.fitting import fitGrain, objFuncFitGrain
from hexrd.xrd.rotations import angleAxisOfRotMat, rotMatOfQuat
from hexrd.xrd.transforms import bVec_ref, eta_ref, mapAngle, vInv_ref
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


def get_frames(reader, cfg, show_progress=False):
    # TODO: this should be updated to read only the frames requested in cfg
    # either the images start, step, stop, or based on omega start, step, stop
    start = time.time()

    n_frames = reader.getNFrames()
    logger.info("reading %d frames of data", n_frames)
    if show_progress:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=n_frames).start()

    frame_list = []
    for i in range(n_frames):
        frame = reader.read()
        mask = frame > cfg.fit_grains.threshold
        sparse_frame = coo_matrix(
            (frame[mask], mask.nonzero()),
            shape=mask.shape
            )
        frame_list.append(sparse_frame)

        if show_progress:
            pbar.update(i)
    frame_list = np.array(frame_list)
    omega_start = np.radians(cfg.image_series.omega.start)
    omega_step = np.radians(cfg.image_series.omega.step)
    reader = [frame_list, [omega_start, omega_step]]
    if show_progress:
        pbar.finish()
    elapsed = time.time()-start
    logger.info('read %d frames in %g seconds', n_frames, elapsed)
    return reader


def get_instrument_parameters(cfg):
    with open(cfg.detector.parameters, 'r') as f:
        # only one panel for now
        # TODO: configurize this
        return [cfg for cfg in yaml.load_all(f)][0]


def get_detector_parameters(cfg, instr_cfg):
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

    return np.hstack([
        instr_cfg['detector']['transform']['tilt_angles'],
        instr_cfg['detector']['transform']['t_vec_d'],
        instr_cfg['oscillation_stage']['chi'],
        instr_cfg['oscillation_stage']['t_vec_s'],
        ])


def get_distortion_correction(instrument_cfg):
    # ***FIX***
    # at this point we know we have a GE and hardwire the distortion func;
    # need to pull name from yml file in general case
    return (
        dFuncs.GE_41RT,
        instrument_cfg['detector']['distortion']['parameters']
        )


def set_planedata_exclusions(cfg, detector, pd):
    tth_max = cfg.fit_grains.tth_max
    if tth_max is True:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() > detector.getTThMax()
    elif tth_max > 0:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() >= np.radians(tth_max)


def get_job_queue(cfg, max_grains=None):
    job_queue = mp.JoinableQueue()
    # load the queue
    try:
        estimate_f = cfg.fit_grains.estimate
        grain_params_list = np.atleast_2d(np.loadtxt(estimate_f))
        n_quats = len(grain_params_list)
        if max_grains is None:
            max_grains = n_quats
        for grain_params in grain_params_list:
            grain_id = grain_params[0]
            if grain_id == max_grains:
                break
            job_queue.put((grain_id, grain_params[3:15]))
        logger.info(
            'fitting grains using "%s" for the initial estimate',
            estimate_f
            )
    except (ValueError, IOError):
        logger.info('fitting grains using default initial estimate')
        # load quaternion file
        quats = np.atleast_2d(
            np.loadtxt(os.path.join(cfg.working_dir, 'accepted_orientations.dat'))
            )
        n_quats = len(quats)
        quats = quats.T
        phi, n = angleAxisOfRotMat(rotMatOfQuat(quats))
        if max_grains is None:
            max_grains = n_quats
        for i, quat in enumerate(quats.T):
            if i == max_grains:
                break
            exp_map = phi[i]*n[:, i]
            grain_params = np.hstack(
                [exp_map.flatten(), 0., 0., 0., 1., 1., 1., 0., 0., 0.]
                )
            job_queue.put((i, grain_params))
    logger.info("fitting grains for %d of %d orientations", max_grains, n_quats)
    return job_queue


def get_data(cfg, show_progress=False):
    # TODO: this should be refactored somehow to avoid initialize_experiment
    # and avoid using the old reader. Also, the detector is not used here.
    pd, reader, detector = initialize_experiment(cfg)
    reader = get_frames(reader, cfg, show_progress)

    instrument_cfg = get_instrument_parameters(cfg)
    detector_params = get_detector_parameters(cfg, instrument_cfg)
    distortion = get_distortion_correction(instrument_cfg)
    set_planedata_exclusions(cfg, detector, pd)
    pkwargs = {
        'distortion': distortion,
        'omega_start': cfg.image_series.omega.start,
        'omega_step': cfg.image_series.omega.step,
        'omega_stop': cfg.image_series.omega.stop,
        'eta_range': np.radians(cfg.find_orientations.eta.range),
        'omega_period': np.radians(cfg.find_orientations.omega.period),
        'tth_tol': cfg.fit_grains.tolerance.tth,
        'eta_tol': cfg.fit_grains.tolerance.eta,
        'omega_tol': cfg.fit_grains.tolerance.omega,
        'panel_buffer': cfg.fit_grains.panel_buffer,
        'npdiv': cfg.fit_grains.npdiv,
        'threshold': cfg.fit_grains.threshold,
        'spots_stem': os.path.join(cfg.analysis_dir, 'spots_%05d.out'),
        'plane_data': pd,
        'detector_params': detector_params,
        }
    return reader, pkwargs


def fit_grains(cfg, force=False, show_progress=False, max_grains=None):
    # load the data
    reader, pkwargs = get_data(cfg, show_progress)
    job_queue = get_job_queue(cfg, max_grains)
    njobs = job_queue.qsize()

    # log this before starting progress bar
    ncpus = cfg.multiprocessing
    ncpus = ncpus if ncpus < njobs else njobs
    logger.info('running pullspots with %d processors', ncpus)
    if ncpus == 1:
        logger.info('multiprocessing disabled')

    start = time.time()
    pbar = None
    if show_progress:
        pbar = ProgressBar(
            widgets=[Bar('>'), ' ', ETA(), ' ', ReverseBar('<')],
            maxval=njobs
            ).start()

    # finally start processing data
    if ncpus == 1:
        # no multiprocessing
        results = []
        w = FitGrainsWorker(
            job_queue, results, reader, copy.deepcopy(pkwargs),
            progressbar=pbar
            )
        w.run()
    else:
        # multiprocessing
        manager = mp.Manager()
        results = manager.list()
        for i in range(ncpus):
            # lets make a deep copy of the pkwargs, just in case:
            w = FitGrainsWorkerMP(job_queue, results, reader, copy.deepcopy(pkwargs))
            w.daemon = True
            w.start()
    while True:
        n_res = len(results)
        if show_progress:
            pbar.update(n_res)
        if n_res == njobs:
            break
        time.sleep(0.1)
    job_queue.join()

    write_grains_file(cfg, results)

    if show_progress:
        pbar.finish()
    elapsed = time.time() - start
    logger.info('processed %d grains in %g minutes', n_res, elapsed/60)


def write_grains_file(cfg, results):
    # record the results to file
    f = open(os.path.join(cfg.analysis_dir, 'grains.out'), 'w')
    # going to some length to make the header line up with the data
    # while also keeping the width of the lines to a minimum, settled
    # on %14.7g representation.
    header_items = (
        'grain ID', 'completeness', 'sum(resd**2)/nrefl',
        'xi[0]', 'xi[1]', 'xi[2]', 'tVec_c[0]', 'tVec_c[1]', 'tVec_c[2]',
        'vInv_s[0]', 'vInv_s[1]', 'vInv_s[2]', 'vInv_s[4]*sqrt(2)',
        'vInv_s[5]*sqrt(2)', 'vInv_s[6]*sqrt(2)', 'ln(V[0,0])',
        'ln(V[1,1])', 'ln(V[2,2])', 'ln(V[1,2])', 'ln(V[0,2])', 'ln(V[0,1])',
        )
    len_items = []
    for i in header_items[1:]:
        temp = len(i)
        len_items.append(temp if temp > 14 else 14) # for %14.7g
    fmtstr = '#%8s  ' + '  '.join(['%%%ds' % i for i in len_items]) + '\n'
    f.write(fmtstr % header_items)
    for (id, g_refined, compl, eMat, resd) in sorted(results):
        res_items = (
            id, compl, resd, g_refined[0], g_refined[1], g_refined[2],
            g_refined[3], g_refined[4], g_refined[5], g_refined[6],
            g_refined[7], g_refined[8], g_refined[9], g_refined[10],
            g_refined[11], eMat[0, 0], eMat[1, 1], eMat[2, 2], eMat[1, 2],
            eMat[0, 2], eMat[0, 1],
            )
        fmtstr = '%9d  ' + '  '.join(['%%%d.7g' % i for i in len_items]) + '\n'
        f.write(fmtstr % res_items)



class FitGrainsWorker(object):


    def __init__(self, jobs, results, reader, pkwargs, **kwargs):
        self._jobs = jobs
        self._results = results
        self._reader = reader
        # a dict containing the rest of the parameters
        self._p = pkwargs

        # lets make a couple shortcuts:
        self._p['bMat'] = np.ascontiguousarray(
            self._p['plane_data'].latVecOps['B']
            ) # is it still necessary to re-cast?
        self._p['wlen'] = self._p['plane_data'].wavelength
        self._pbar = kwargs.get('progressbar', None)


    def pull_spots(self, grain_id, grain_params, iteration):
        return pullSpots(
            self._p['plane_data'],
            self._p['detector_params'],
            grain_params,
            self._reader,
            distortion=self._p['distortion'],
            eta_range=self._p['eta_range'],
            ome_period=self._p['omega_period'],
            eta_tol=self._p['eta_tol'][iteration],
            ome_tol=self._p['omega_tol'][iteration],
            tth_tol=self._p['tth_tol'][iteration],
            panel_buff=self._p['panel_buffer'],
            npdiv=self._p['npdiv'],
            threshold=self._p['threshold'],
            doClipping=False,
            filename=self._p['spots_stem'] % grain_id,
            )


    def fit_grains(self, grain_id, grain_params):
        ome_start = self._p['omega_start']
        ome_step = self._p['omega_step']
        ome_stop =  self._p['omega_stop']
        gtable = np.loadtxt(self._p['spots_stem'] % grain_id)
        valid_refl_ids = gtable[:, 0] >= 0
        pred_ome = gtable[:, 6]
        if np.sign(ome_step) < 0:
            idx_ome = np.logical_and(
                pred_ome < np.radians(ome_start + 2*ome_step),
                pred_ome > np.radians(ome_stop - 2*ome_step)
                )
        else:
            idx_ome = np.logical_and(
                pred_ome > np.radians(ome_start + 2*ome_step),
                pred_ome < np.radians(ome_stop - 2*ome_step)
                )

        idx = np.logical_and(valid_refl_ids, idx_ome)
        hkls = gtable[idx, 1:4].T # must be column vectors
        self._p['hkls'] = hkls
        xyo_det = gtable[idx, -3:] # these are the cartesian centroids + ome
        xyo_det[:, 2] = mapAngle(xyo_det[:, 2], self._p['omega_period'])
        self._p['xyo_det'] = xyo_det
        if sum(idx) <= 12:
            return grain_params, 0
        grain_params = fitGrain(
            xyo_det, hkls, self._p['bMat'], self._p['wlen'],
            self._p['detector_params'],
            grain_params[:3], grain_params[3:6], grain_params[6:],
            beamVec=bVec_ref, etaVec=eta_ref,
            distortion=self._p['distortion'],
            gFlag=gFlag, gScl=gScl,
            omePeriod=self._p['omega_period']
            )
        completeness = sum(idx)/float(len(idx))
        return grain_params, completeness


    def get_e_mat(self, grain_params):
        # TODO: document what is this?
        return logm(np.linalg.inv(vecMVToSymm(grain_params[6:])))


    def get_residuals(self, grain_params):
        dFunc, dParams = self._p['distortion']
        return objFuncFitGrain(
            grain_params[gFlag], grain_params, gFlag,
            self._p['detector_params'],
            self._p['xyo_det'], self._p['hkls'],
            self._p['bMat'], self._p['wlen'],
            bVec_ref, eta_ref,
            dFunc, dParams,
            self._p['omega_period'],
            simOnly=False
            )


    def loop(self):
        id, grain_params = self._jobs.get(False)

        # skips the first loop if have_estimate is True
        have_estimate = not np.all(grain_params[-9] == [0,0,0,1,1,1,0,0,0])
        iterations = (have_estimate, len(self._p['eta_tol']))
        for iteration in range(*iterations):
            self.pull_spots(id, grain_params, iteration)
            grain_params, compl = self.fit_grains(id, grain_params)
            if compl == 0:
                break

        eMat = self.get_e_mat(grain_params)
        resd = self.get_residuals(grain_params)

        self._results.append((id, grain_params, compl, eMat, sum(resd**2)))
        self._jobs.task_done()


    def run(self):
        n_res = 0
        while True:
            try:
                self.loop()
                n_res += 1
                if self._pbar is not None:
                    self._pbar.update(n_res)
            except Empty:
                break



class FitGrainsWorkerMP(FitGrainsWorker, mp.Process):

    def __init__(self, *args, **kwargs):
        mp.Process.__init__(self)
        FitGrainsWorker.__init__(self, *args, **kwargs)
