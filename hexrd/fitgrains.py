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

from hexrd.coreutil import (
    initialize_experiment, migrate_detector_to_instrument_config,
    get_instrument_parameters, get_detector_parameters, get_detector_parameters,
    get_distortion_correction, get_saturation_level, set_planedata_exclusions
    )
from hexrd.matrixutil import vecMVToSymm
from hexrd.utils.progressbar import (
    Bar, ETA, Percentage, ProgressBar, ReverseBar
    )

from hexrd.xrd import distortion as dFuncs
from hexrd.xrd.fitting import fitGrain, objFuncFitGrain
from hexrd.xrd.rotations import angleAxisOfRotMat, rotMatOfQuat
from hexrd.xrd.transforms import bVec_ref, eta_ref, mapAngle, vInv_ref, angularDifference
from hexrd.xrd.xrdutil import pullSpots
from .cacheframes import get_frames
from hexrd import USE_NUMBA
if USE_NUMBA:
    import numba

logger = logging.getLogger(__name__)

# grain parameter refinement flags
gFlag = np.array([1, 1, 1,
                  1, 1, 1,
                  1, 1, 1, 1, 1, 1], dtype=bool)
# grain parameter scalings
gScl  = np.array([1., 1., 1.,
                  1., 1., 1.,
                  1., 1., 1., 0.01, 0.01, 0.01])


def get_job_queue(cfg, ids_to_refine=None):
    job_queue = mp.JoinableQueue()
    # load the queue
    try:
        # use an estimate of the grain parameters, if available
        estimate_f = cfg.fit_grains.estimate
        grain_params_list = np.atleast_2d(np.loadtxt(estimate_f))
        n_quats = len(grain_params_list)
        n_jobs = 0
        for grain_params in grain_params_list:
            grain_id = grain_params[0]
            if ids_to_refine is None or grain_id in ids_to_refine:
                job_queue.put((grain_id, grain_params[3:15]))
                n_jobs += 1
        logger.info(
            'fitting grains using "%s" for the initial estimate',
            estimate_f
            )
    except (ValueError, IOError):
        # no estimate available, use orientations and defaults
        logger.info('fitting grains using default initial estimate')
        
        # ...make this an attribute in cfg?
        analysis_id = '%s_%s' %(
            cfg.analysis_name.strip().replace(' ', '-'),
            cfg.material.active.strip().replace(' ', '-'),
            )
        
        # load quaternion file
        quats = np.atleast_2d(
            np.loadtxt(
                os.path.join(
                    cfg.working_dir, 
                    'accepted_orientations_%s.dat' %analysis_id
                    )
                )
            )
        n_quats = len(quats)
        n_jobs = 0
        phi, n = angleAxisOfRotMat(rotMatOfQuat(quats.T))
        for i, (phi, n) in enumerate(zip(phi, n.T)):
            if ids_to_refine is None or i in ids_to_refine:
                exp_map = phi*n
                grain_params = np.hstack(
                    [exp_map, 0., 0., 0., 1., 1., 1., 0., 0., 0.]
                    )
                job_queue.put((i, grain_params))
                n_jobs += 1
    logger.info("fitting grains for %d of %d orientations", n_jobs, n_quats)
    return job_queue, n_jobs


def get_data(cfg, show_progress=False, force=False, clean=False):
    # TODO: this should be refactored somehow to avoid initialize_experiment
    # and avoid using the old reader. Also, the detector is not used here.
    pd, reader, detector = initialize_experiment(cfg)
    if cfg.fit_grains.fit_only:
        reader = None
    else:
        reader = get_frames(reader, cfg, show_progress, force, clean)

    instrument_cfg = get_instrument_parameters(cfg)
    detector_params = get_detector_parameters(instrument_cfg)
    saturation_level = get_saturation_level(instrument_cfg)
    distortion = get_distortion_correction(instrument_cfg)
    set_planedata_exclusions(cfg, detector, pd)
    pkwargs = {
        'detector_params': detector_params,
        'distortion': distortion,
        'eta_range': np.radians(cfg.find_orientations.eta.range),
        'eta_tol': cfg.fit_grains.tolerance.eta,
        'fit_only': cfg.fit_grains.fit_only,
        'ncols': instrument_cfg['detector']['pixels']['columns'],
        'npdiv': cfg.fit_grains.npdiv,
        'nrows': instrument_cfg['detector']['pixels']['rows'],        
        'omega_period': np.radians(cfg.find_orientations.omega.period),
        'omega_start': cfg.image_series.omega.start,
        'omega_step': cfg.image_series.omega.step,
        'omega_stop': cfg.image_series.omega.stop,
        'omega_tol': cfg.fit_grains.tolerance.omega,
        'overlap_table': os.path.join(cfg.analysis_dir, 'overlap_table.npz'),
        'panel_buffer': cfg.fit_grains.panel_buffer,
        'pixel_pitch': instrument_cfg['detector']['pixels']['size'],
        'plane_data': pd,
        'refit_tol': cfg.fit_grains.refit,
        'saturation_level': saturation_level,
        'spots_stem': os.path.join(cfg.analysis_dir, 'spots_%05d.out'),
        'threshold': cfg.fit_grains.threshold,
        'tth_tol': cfg.fit_grains.tolerance.tth,
        }
    return reader, pkwargs


def fit_grains(cfg, force=False, clean=False, show_progress=False, ids_to_refine=None):
    # load the data
    reader, pkwargs = get_data(cfg, show_progress, force, clean)
    job_queue, njobs = get_job_queue(cfg, ids_to_refine)

    # log this before starting progress bar
    ncpus = cfg.multiprocessing
    ncpus = ncpus if ncpus < njobs else njobs
    logger.info(
        'will use %d of %d processors', ncpus, mp.cpu_count()
        )
    if ncpus == 1:
        logger.info('multiprocessing disabled')

    # echo some of the fitting options
    if cfg.fit_grains.fit_only:
        logger.info('\t**fitting only; will not pull spots')
    if cfg.fit_grains.refit is not None:
        msg = 'will perform refit excluding spots > ' + \
              '%.2f pixels and ' %cfg.fit_grains.refit[0] + \
              '%.2f frames from expected values' %cfg.fit_grains.refit[1]
        logger.info(msg)
    
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


def write_grains_file(cfg, results, output_name=None):
    # record the results to file
    if output_name is None:
        f = open(os.path.join(cfg.analysis_dir, 'grains.out'), 'w')
    else:
        f = open(os.path.join(cfg.analysis_dir, output_name), 'w')
    # going to some length to make the header line up with the data
    # while also keeping the width of the lines to a minimum, settled
    # on %19.12g representation.
    header_items = (
        'grain ID', 'completeness', 'chi2',
        'xi[0]', 'xi[1]', 'xi[2]', 'tVec_c[0]', 'tVec_c[1]', 'tVec_c[2]',
        'vInv_s[0]', 'vInv_s[1]', 'vInv_s[2]', 'vInv_s[4]*sqrt(2)',
        'vInv_s[5]*sqrt(2)', 'vInv_s[6]*sqrt(2)', 'ln(V[0,0])',
        'ln(V[1,1])', 'ln(V[2,2])', 'ln(V[1,2])', 'ln(V[0,2])', 'ln(V[0,1])',
        )
    len_items = []
    for i in header_items[1:]:
        temp = len(i)
        len_items.append(temp if temp > 19 else 19) # for %19.12g
    fmtstr = '#%13s  ' + '  '.join(['%%%ds' % i for i in len_items]) + '\n'
    f.write(fmtstr % header_items)
    for (id, g_refined, compl, eMat, resd) in sorted(results):
        res_items = (
            id, compl, resd, g_refined[0], g_refined[1], g_refined[2],
            g_refined[3], g_refined[4], g_refined[5], g_refined[6],
            g_refined[7], g_refined[8], g_refined[9], g_refined[10],
            g_refined[11], eMat[0, 0], eMat[1, 1], eMat[2, 2], eMat[1, 2],
            eMat[0, 2], eMat[0, 1],
            )
        fmtstr = (
            '%14d  ' + '  '.join(['%%%d.12g' % i for i in len_items]) + '\n'
            )
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
        # need to calc panel dims on the fly
        xdim = self._p['pixel_pitch'][1] * self._p['ncols']
        ydim = self._p['pixel_pitch'][0] * self._p['nrows']
        panel_dims = [(-0.5*xdim, -0.5*ydim),
                      ( 0.5*xdim,  0.5*ydim)]
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
            pixel_pitch=self._p['pixel_pitch'],
            panel_dims=panel_dims,
            panel_buff=self._p['panel_buffer'],
            npdiv=self._p['npdiv'],
            threshold=self._p['threshold'],
            doClipping=False,
            filename=self._p['spots_stem'] % grain_id,
            )


    def fit_grains(self, grain_id, grain_params, refit_tol=None):
        """
        Executes lsq fits of grains based on spot files
        
        REFLECTION TABLE
        
        Cols as follows:
            0-6:   ID    PID    H    K    L    sum(int)    max(int)    
            6-9:   pred tth    pred eta    pred ome              
            9-12:  meas tth    meas eta    meas ome              
            12-15: meas X      meas Y      meas ome
        """
        ome_start = self._p['omega_start']
        ome_step = self._p['omega_step']
        ome_stop =  self._p['omega_stop']
        refl_table = np.loadtxt(self._p['spots_stem'] % grain_id)
        valid_refl_ids = refl_table[:, 0] >= 0
        unsat_spots = refl_table[:, 6] < self._p['saturation_level']
        pred_ome = refl_table[:, 9]
        if angularDifference(ome_start, ome_stop, units='degrees') > 0:
            # if here, incomplete have omega range and
            # clip the refelctions very close to the edges to avoid
            # problems with the least squares...
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
            idx = np.logical_and(
                valid_refl_ids,
                np.logical_and(unsat_spots, idx_ome)
                )
        else:
            idx = np.logical_and(valid_refl_ids, unsat_spots)
            pass # end if edge case

        # if an overlap table has been written, load it and use it
        overlaps = np.zeros(len(refl_table), dtype=bool)
        try:
            ot = np.load(self._p['overlap_table'])
            for key in ot.keys():
                for this_table in ot[key]:
                    these_overlaps = np.where(
                        this_table[:, 0] == grain_id)[0]
                    if len(these_overlaps) > 0:
                        mark_these = np.array(this_table[these_overlaps, 1], dtype=int)
                        overlaps[mark_these] = True
            idx = np.logical_and(idx, ~overlaps)
        except IOError, IndexError:
            #print "no overlap table found"
            pass
        
        # completeness from pullspots only; incl saturated and overlaps
        completeness = sum(valid_refl_ids)/float(len(valid_refl_ids))

        # extract data from grain table
        hkls = refl_table[idx, 2:5].T # must be column vectors
        xyo_det = refl_table[idx, -3:] # these are the cartesian centroids + ome

        # set in parameter attribute
        self._p['hkls'] = hkls
        self._p['xyo_det'] = xyo_det
        
        if sum(idx) <= 12: # not enough reflections to fit... exit
            completeness = 0.
        else:
            grain_params = fitGrain(
                xyo_det, hkls, self._p['bMat'], self._p['wlen'],
                self._p['detector_params'],
                grain_params[:3], grain_params[3:6], grain_params[6:],
                beamVec=bVec_ref, etaVec=eta_ref,
                distortion=self._p['distortion'],
                gFlag=gFlag, gScl=gScl,
                omePeriod=self._p['omega_period']
                )
            if refit_tol is not None:
                xpix_tol = refit_tol[0]*self._p['pixel_pitch'][1]
                ypix_tol = refit_tol[0]*self._p['pixel_pitch'][0]
                fome_tol = refit_tol[1]*self._p['omega_step']
                
                xyo_det_fit = objFuncFitGrain(
                    grain_params[gFlag], grain_params, gFlag,
                    self._p['detector_params'],
                    xyo_det, hkls, self._p['bMat'], self._p['wlen'],
                    bVec_ref, eta_ref,
                    self._p['distortion'][0], self._p['distortion'][1],
                    self._p['omega_period'], simOnly=True
                    )

                # define difference vectors for spot fits
                x_diff = abs(xyo_det[:, 0] - xyo_det_fit[:, 0])
                y_diff = abs(xyo_det[:, 1] - xyo_det_fit[:, 1])
                ome_diff = np.degrees(
                    angularDifference(xyo_det[:, 2], xyo_det_fit[:, 2])
                    )

                # filter out reflections with centroids more than 
                # a pixel and delta omega away from predicted value
                idx_1 = np.logical_and(
                    x_diff <= xpix_tol,
                    np.logical_and(y_diff <= ypix_tol,
                                   ome_diff <= fome_tol)
                                   )
                idx_new = np.zeros_like(idx, dtype=bool)
                idx_new[np.where(idx == 1)[0][idx_1]] = True

                if sum(idx_new) > 12 and (sum(idx_new) > 0.5*sum(idx)):
                    # have enough reflections left
                    #   ** the check that we have more than half of what
                    #      we started with is a hueristic
                    hkls = refl_table[idx_new, 2:5].T
                    xyo_det = refl_table[idx_new, -3:]

                    # set in parameter attribute
                    self._p['hkls'] = hkls
                    self._p['xyo_det'] = xyo_det

                    # do fit
                    grain_params = fitGrain(
                        xyo_det, hkls,
                        self._p['bMat'], self._p['wlen'],
                        self._p['detector_params'],
                        grain_params[:3], grain_params[3:6], grain_params[6:],
                        beamVec=bVec_ref, etaVec=eta_ref,
                        distortion=self._p['distortion'],
                        gFlag=gFlag, gScl=gScl,
                        omePeriod=self._p['omega_period']
                        )
                    pass # end check on num of refit refls
                pass # end refit loop
            pass # end on num of refls
        return grain_params, completeness
        

    def get_e_mat(self, grain_params):
        """
        strain tensor calculation
        """
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
            simOnly=False, return_value_flag=2)

    def loop(self):
        id, grain_params = self._jobs.get(False)
        iterations = (0, len(self._p['eta_tol']))
        for iteration in range(*iterations):
            # pull spots if asked to, otherwise just fit
            if not self._p['fit_only']:
                self.pull_spots(id, grain_params, iteration)
            # FITTING HERE
            grain_params, compl = self.fit_grains(id, grain_params,
                                                  refit_tol=self._p['refit_tol'])
            if compl == 0:
                break
            pass
        
        # final pull spots if enabled
        if not self._p['fit_only']:
            self.pull_spots(id, grain_params, -1)

        eMat = self.get_e_mat(grain_params)
        resd = self.get_residuals(grain_params)

        self._results.append((id, grain_params, compl, eMat, resd))
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
