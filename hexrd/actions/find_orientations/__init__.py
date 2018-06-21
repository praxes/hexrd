"""find_orientations command"""
from __future__ import print_function, division, absolute_import

import timeit

from .utils import get_eta_ome, generate_orientation_fibers

def find_orientations(cfg, hkls=None, clean=False, profile=False):
    print('ready to run find_orientations')
    # %%
    # =============================================================================
    # SEARCH SPACE GENERATION
    # =============================================================================

    hedm = cfg.instrument.hedm

    ncpus = cfg.multiprocessing

    # for indexing
    fiber_ndiv = cfg.find_orientations.seed_search.fiber_ndiv
    fiber_seeds = cfg.find_orientations.seed_search.hkl_seeds
    on_map_threshold = cfg.find_orientations.threshold

    print("INFO:\tgenerating search quaternion list using %d processes" % ncpus)
    start = timeit.default_timer()

    eta_ome = get_eta_ome(cfg)
    qfib = generate_orientation_fibers(
        eta_ome, hedm.chi, on_map_threshold,
        fiber_seeds, fiber_ndiv,
        ncpus=ncpus
    )
    print("INFO:\t\t...took %f seconds" % (timeit.default_timer() - start))
    print("INFO: will test %d quaternions using %d processes"
          % (qfib.shape[1], ncpus))
