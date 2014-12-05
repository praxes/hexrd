from __future__ import absolute_import

import copy
import logging
import os
import sys
import time

import yaml

import numpy as np
from scipy.sparse import coo_matrix

from hexrd.coreutil import initialize_experiment
from hexrd.utils.progressbar import (
    Bar, ETA, Percentage, ProgressBar, ReverseBar
    )

from hexrd import USE_NUMBA
if USE_NUMBA:
    import numba

logger = logging.getLogger(__name__)

if USE_NUMBA:
    @numba.njit
    def extract_ijv(in_array, threshold, out_i, out_j, out_v):
        n = 0
        w, h = in_array.shape

        for i in range(w):
            for j in range(h):
                v = in_array[i,j]
                if v > threshold:
                    out_v[n] = v
                    out_i[n] = i
                    out_j[n] = j
                    n += 1

        return n

    class CooMatrixBuilder(object):
        def __init__(self):
            self.v_buff = np.empty((2048*2048,), dtype=np.int16)
            self.i_buff = np.empty((2048*2048,), dtype=np.int16)
            self.j_buff = np.empty((2048*2048,), dtype=np.int16)

        def build_matrix(self, frame, threshold):
            count = extract_ijv(frame, threshold,
                                self.i_buff, self.j_buff, self.v_buff)
            return coo_matrix((self.v_buff[0:count].copy(),
                               (self.i_buff[0:count].copy(),
                                self.j_buff[0:count].copy())),
                              shape=frame.shape)

else: # not USE_NUMBA
    class CooMatrixBuilder(object):
        def build_matrix(self, frame, threshold):
            mask = frame > threshold
            return coo_matrix((frame[mask], mask.nonzero()),
                              shape=frame.shape)

def load_frames(reader, cfg, show_progress=False):
    # TODO: this should be updated to read only the frames requested in cfg
    # either the images start, step, stop, or based on omega start, step, stop
    start = time.time()

    n_frames = reader.getNFrames()
    logger.info("reading %d frames of data, storing values > %.1f", n_frames, cfg.fit_grains.threshold)
    if show_progress:
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=n_frames).start()

    frame_list = []
    coo_builder = CooMatrixBuilder()
    for i in range(n_frames):
        frame = reader.read()
        frame_list.append(coo_builder.build_matrix(frame, cfg.fit_grains.threshold))

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

def cache_frames(reader, cfg, show_progress=False, overwrite=True):
    cache_file = os.path.join(cfg.analysis_dir, 'frame_cache.npz')
    # load the data
    reader = load_frames(reader, cfg, show_progress)
    # save all the data to a .npz file
    arrs = {}
    arrs['omega'] = np.array(reader[1])
    arrs['shape'] = np.array(reader[0][0].shape)
    for i, coo in enumerate(reader[0]):
        arrs['%d_data' % i] = coo.data
        arrs['%d_row' % i] = coo.row
        arrs['%d_col' % i] = coo.col
    start = time.time()
    np.savez_compressed(cache_file, **arrs)
    elapsed = time.time()-start
    logger.info('wrote %d frames to cache in %g seconds', len(reader[0]), elapsed)
    return reader

def get_frames(reader, cfg, show_progress=False, force=False, clean=False):
    cache_file = os.path.join(cfg.analysis_dir, 'frame_cache.npz')
    if not os.path.exists(cache_file) or clean:
        if clean:
            msg = 'no frame cache file %s found, generating cache' % cache_file
        else:
            msg = 'clean specified, generating cache'
        logger.info(msg)
        return cache_frames(reader, cfg, show_progress)

    start = time.time()

    n_frames = reader.getNFrames()
    logger.info("reading %d frames from cache", n_frames)

    with np.load(cache_file) as npz:
        omega = npz['omega'].tolist()
        shape = npz['shape'].tolist()
        frame_list = []
        for i in range(n_frames):
            frame_list.append(coo_matrix((npz['%d_data' % i],
                                          (npz['%d_row' % i],
                                           npz['%d_col' % i])),
                                         shape=shape))
    omega_start = np.radians(cfg.image_series.omega.start)
    omega_step = np.radians(cfg.image_series.omega.step)
    if omega_start != omega[0] or omega_step != omega[1]:
        logger.warning('inconsistent omegas %s in cache file vs config %s'
                       % (omega, [omega_start, omega_step]))
    reader = [frame_list, [omega_start, omega_step]]
    elapsed = time.time()-start
    logger.info('read %d cached frames in %g seconds', n_frames, elapsed)
    return reader
