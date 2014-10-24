import argparse
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

from hexrd.coreutil import (
    initialize_experiment, iter_cfg_sections, make_eta_ranges, merge_dicts,
    migrate_detector_config,
    )

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd import rotations as rot
from hexrd.xrd.xrdutil import pullSpots as pull_spots


def extract_g_vectors(
    cfg, verbose=False, force=False
    ):
    """ takes a cfg dict, not a file """

    pd, reader, detector = initialize_experiment(cfg)

    #####################
    ## load parameters ##
    #####################

    cwd = cfg.get('working_dir', os.getcwd())
    analysis_root = os.path.join(cwd, cfg['analysis_name'])

    exgcfg = cfg['extract_grains']
    pthresh = exgcfg['threshold']
    tth_max = exgcfg.get('tth_max', True)
    ome_start = cfg['image_series']['ome']['start']
    ome_step = cfg['image_series']['ome']['step']
    try:
        tth_tol = exgcfg['tolerance'].get('tth', None)
        eta_tol = exgcfg['tolerance'].get('eta', None)
        ome_tol = exgcfg['tolerance'].get('ome', None)
    except KeyError:
        tth_tol = eta_tol = ome_tol = None
    if tth_tol is None:
        tth_tol = 0.2
        if verbose:
            print "tth tolerance is %g" % tth_tol
    if eta_tol is None:
        eta_tol = 1
        if verbose:
            print "eta tolerance is %g" % eta_tol
    if ome_tol is None:
        ome_tol = 2*ome_step
        if verbose:
            "ome tolerance is %g" % ome_tol

    try:
        eta_mask = abs(cfg['find_orientations']['eta'].get('mask', 5))
    except (KeyError, AttributeError):
        eta_mask = 5

    try:
        ome_period = cfg['find_orientations']['ome'].get('period', None)
    except (KeyError, AttributeError):
        ome_period = None
    if ome_period is None:
        if ome_step > 0:
            ome_period = [ome_start, ome_start + 360]
        else:
            ome_period = [ome_start, ome_start - 360]
        if verbose:
            print "Omega tolerance: %s" % ome_tol
            print "Omega period: %s" % ome_period
    ome_period = np.radians(ome_period)

    if eta_mask:
        eta_range = make_eta_ranges(eta_mask)
        if verbose:
            print (
                "Masking eta angles within %g degrees of ome rotation axis"
                % eta_mask
                )
    else:
        if verbose:
            print "Using full eta range"

    panel_buffer = exgcfg.get('panel_buffer', 10)
    if isinstance(panel_buffer, int):
        panel_buffer = [panel_buffer, panel_buffer]
    npdiv = exgcfg.get('pixel_subdivisions', 2)

    # determine number of processes to run in parallel
    multiproc = cfg.get('multiprocessing', -1)
    ncpus = mp.cpu_count()
    if multiproc == 'all':
        pass
    elif multiproc == -1:
        ncpus -= 1
    elif int(ncpus) == 'half':
        ncpus /= 2
    elif isinstance(multiproc, int):
        if multiproc < ncpus:
            ncpus = multiproc
    else:
        ncpus -= 1
        if verbose:
            print (
                "Invalid value %s for find_orientations:multiprocessing"
                % multiproc
                )
    ncpus = ncpus if ncpus else 1

    ###########################
    ## Instrument parameters ##
    ###########################

    # attempt to load the new detector parameter file
    det_p = os.path.join(cwd, cfg['detector']['parameters'])
    if not os.path.exists(det_p):
        det_o = os.path.join(cwd, cfg['detector']['parameters_old'])
        nrows = cfg['detector']['pixels']['rows']
        ncols = cfg['detector']['pixels']['columns']
        psize = cfg['detector']['pixels']['size']
        old_par = np.loadtxt(det_o)
        migrate_detector_config(old_par, nrows, ncols, psize,
                                detID='GE', chi=0., tVec_s=np.zeros(3),
                                filename=det_p)

    with open(det_p, 'r') as f:
        # only one panel for now
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

    if tth_max is True:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() > detector.getTThMax()
    elif tth_max > 0:
        pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
        pd.exclusions = pd.getTTh() >= np.radians(tth_max)
    elif tth_max < 0:
        if verbose:
            print "Ignoring invalid tth_max: %g."
            "Must be 'true', 'false', or a non-negative value"

    # load quaternion file
    quats = np.atleast_2d(
        np.loadtxt(os.path.join(analysis_root, 'quats.out'))
        )
    n_quats = len(quats)
    quats = quats.T

    phi, n = rot.angleAxisOfRotMat(rot.rotMatOfQuat(quats))

    cwd = cfg.get('working_dir', os.getcwd())
    analysis_name = cfg['analysis_name']
    spots_f = os.path.join(cwd, analysis_name, 'spots_%05d.out')

    job_queue = mp.JoinableQueue()
    manager = mp.Manager()
    results = manager.list()

    n_frames = reader.getNFrames()
    if verbose:
        print "reading %d frames of data" % n_frames
        if have_progBar:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=n_frames).start()

    frame_list = []
    for i in range(n_frames):
        frame = reader.read()
        frame[frame <= pthresh] = 0
        frame_list.append(coo_matrix(frame))
        if have_progBar and verbose:
            pbar.update(i)
    frame_list = np.array(frame_list)
    reader = [frame_list, [np.radians(ome_start), np.radians(ome_step)]]
    if have_progBar and verbose:
        pbar.finish()

    if verbose:
        print "pulling spots for %d orientations...\n" % n_quats
        if have_progBar:
            widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
            pbar = ProgressBar(widgets=widgets, maxval=n_quats).start()

    for i, quat in enumerate(quats.T):
        exp_map = phi[i]*n[:, i]
        grain_params = np.hstack(
            [exp_map.flatten(), 0., 0., 0., 1., 1., 1., 0., 0., 0.]
            )
        job_queue.put((i, grain_params))

    for i in range(ncpus):
        w = PullSpotsWorker(
            job_queue, results,
            reader, pd, detector_params, distortion, eta_range,
            ome_period, tth_tol[0], eta_tol[0], ome_tol[0], panel_buffer, npdiv,
            pthresh, spots_f
            )
        w.daemon = True
        w.start()
    while True:
        n_res = len(results)
        if have_progBar and verbose:
            pbar.update(n_res)
        if n_res == n_quats:
            break
    job_queue.join()
    if have_progBar and verbose:
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


    def run(self):
        while True:
            try:
                i, grain_params = self._jobs.get(False)
                res = pull_spots(
                    self._plane_data,
                    self._detector_params,
                    grain_params,
                    self._reader,
                    distortion=self._distortion,
                    eta_range=self._eta_range,
                    ome_period=self._ome_period,
                    tth_tol=self._tth_tol,
                    eta_tol=self._eta_tol,
                    ome_tol=self._ome_tol,
                    panel_buff=self._panel_buff,
                    npdiv=self._npdiv,
                    threshold=self._pthresh,
                    filename=self._spots_f % i,
                    )
                self._results.append((i, res))
                self._jobs.task_done()
            except Empty:
                break
