import argparse
import os
import sys
import textwrap
import time

import yaml

import numpy as np

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except IOError:
    pass

from hexrd.coreutil import (
    initialize_experiment, merge_dicts,
    migrate_detector_config, make_eta_ranges
    )

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd import rotations as rot
from hexrd.xrd import xrdutil


def _extract_measured_g_vectors(
    cfg, pd, reader, detector, verbose=False, force=False
    ):
    """ takes a cfg dict, not a file """

    #####################
    ## load parameters ##
    #####################

    cwd = cfg.get('working_dir', os.getcwd())
    analysis_root = os.path.join(cwd, cfg['analysis_name'])

    exgcfg = cfg['extract_measured_g_vectors']
    pthresh = exgcfg['threshold']
    tth_max = exgcfg.get('tth_max', True)
    try:
        tth_tol = exgcfg['tolerance'].get('tth', None)
        eta_tol = exgcfg['tolerance'].get('eta', None)
        ome_tol = exgcfg['tolerance'].get('ome', None)
    except KeyError:
        tth_tol = eta_tol = ome_tol = None
    if tth_tol is None:
        tth_tol = 0.2
    if eta_tol is None:
        eta_tol = 1
    if ome_tol is None:
        ome_tol = 2*cfg['image_series']['ome']['step']

    try:
        eta_mask = abs(cfg['find_orientations']['eta'].get('mask', 5))
    except (KeyError, AttributeError):
        eta_mask = 5

    try:
        ome_period = cfg['find_orientations']['ome'].get('period', None)
    except (KeyError, AttributeError):
        ome_period = None
    if ome_period is None:
        temp = cfg['image_series']['ome']['start']
        if cfg['image_series']['ome']['step'] > 0:
            ome_period = [temp, temp + 360]
        else:
            ome_period = [temp, temp - 360]
        if verbose:
            print "Omega tolerance: %g" % ome_tol
            print "Omega period: %s" % ome_period
    ome_period = np.radians(ome_period)

    if eta_mask:
        eta_range = make_eta_ranges(eta_mask)
        if verbose:
            print (
                "Masking eta angles within %g of ome rotation axis"
                % eta_mask
                )
    else:
        if verbose:
            print "Using full eta range"

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
    quats = np.atleast_2d(np.loadtxt(os.path.join(analysis_root, 'quats.out'))).T

    phi, n = rot.angleAxisOfRotMat(rot.rotMatOfQuat(quats))
    if have_progBar:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=len(quats.T)).start()
    if verbose:
        print "pulling spots for %d orientations...\n" %len(quats.T)

    cwd = cfg.get('working_dir', os.getcwd())
    analysis_name = cfg['analysis_name']
    spots_f = os.path.join(cwd, analysis_name, 'spots_%05d.out')
    for iq, quat in enumerate(quats.T):
        if have_progBar:
            pbar.update(iq)
        exp_map = phi[iq]*n[:, iq]
        grain_params = np.hstack([exp_map.flatten(), 0., 0., 0., 1., 1., 1., 0., 0., 0.])
        sd = xrdutil.pullSpots(
            pd,
            detector_params,
            grain_params,
            reader,
            distortion=distortion,
            eta_range=eta_range,
            ome_period=ome_period,
            tth_tol=tth_tol,
            eta_tol=eta_tol,
            ome_tol=ome_tol,
            panel_buff=[10, 10],
            npdiv=2,
            threshold=pthresh,
            filename=spots_f % iq,
            )
    if have_progBar:
        pbar.finish()


def extract_g_vectors(cfg, verbose=False, force=False):
    if verbose:
        print "Using '%s' configuration file" % cfg

    # need to iterate here
    # for cfg in cfgs
    with open(cfg, 'r') as f:
        cfgs = [cfg for cfg in yaml.load_all(f)]
    cfg = cfgs[0]

    # a goofy call, could be replaced with two more targeted calls
    pd, reader, detector = initialize_experiment(cfg)

    for i, c in enumerate(cfgs):
        cfg = merge_dicts(cfg, c)
        extract_measured_g_vectors(
            cfg, pd, reader, detector,
            verbose=verbose, force=force
        )
