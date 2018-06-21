"""Functions used in find_orientations"""
from __future__ import print_function, division, absolute_import

import os
import time
import logging
import multiprocessing as mp

import numpy as np
from scipy import ndimage
import timeit

from hexrd import instrument
from hexrd import matrixutil as mutil
from hexrd.xrd.xrdutil import EtaOmeMaps
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import rotations as rot

print (__name__)
logger = logging.getLogger(__name__)

# ==================== Hardwired options
# maps options
clobber_maps = False
show_maps = False

# ==================== Functions

def analysis_id(cfg):
    return '%s_%s' % (
        cfg.analysis_name.strip().replace(' ', '-'),
        cfg.material.active.strip().replace(' ', '-')
    )


def get_eta_ome(cfg):
    """Return eta-omega maps"""
    # Use existing ones if available
    maps_fname = analysis_id(cfg) + "_maps.npz"
    if os.path.exists(maps_fname) and not clobber_maps:
        print("INFO: loading existing eta_ome maps")
        eta_ome = EtaOmeMaps(maps_fname)
        return eta_ome

    print("INFO: building eta_ome maps")
    start = timeit.default_timer()

    # make eta_ome maps
    imsd = cfg.image_series
    instr = cfg.instrument.hedm
    plane_data = cfg.material.plane_data
    active_hkls = cfg.find_orientations.orientation_maps.active_hkls
    build_map_threshold = cfg.find_orientations.orientation_maps.threshold
    ome_period = np.radians(cfg.find_orientations.omega.period)

    eta_ome = instrument.GenerateEtaOmeMaps(
        imsd, instr, plane_data,
        active_hkls=active_hkls,
        threshold=build_map_threshold,
        ome_period=cfg.find_orientations.omega.period
    )

    print("INFO:  ...took %f seconds" % (timeit.default_timer() - start))

    # save them
    eta_ome.save(maps_fname)

    return eta_ome

# ============================== Fibers

def generate_orientation_fibers(
        eta_ome, chi, threshold, seed_hkl_ids, fiber_ndiv,
        filt_stdev=0.8, ncpus=1):
    """
    From ome-eta maps and hklid spec, generate list of
    quaternions from fibers
    """
    # seed_hkl_ids must be consistent with this...
    pd_hkl_ids = eta_ome.iHKLList[seed_hkl_ids]

    # grab angular grid infor from maps
    del_ome = eta_ome.omegas[1] - eta_ome.omegas[0]
    del_eta = eta_ome.etas[1] - eta_ome.etas[0]

    # labeling mask
    structureNDI_label = ndimage.generate_binary_structure(2, 1)

    # crystallography data from the pd object
    pd = eta_ome.planeData
    hkls = pd.hkls
    tTh = pd.getTTh()
    bMat = pd.latVecOps['B']
    csym = pd.getLaueGroup()

    params = dict(
        bMat=bMat,
        chi=chi,
        csym=csym,
        fiber_ndiv=fiber_ndiv)

    # =========================================================================
    # Labeling of spots from seed hkls
    # =========================================================================

    qfib = []
    input_p = []
    numSpots = []
    coms = []
    for i in seed_hkl_ids:
        # First apply filter
        this_map_f = -ndimage.filters.gaussian_laplace(
            eta_ome.dataStore[i], filt_stdev)

        labels_t, numSpots_t = ndimage.label(
            this_map_f > threshold,
            structureNDI_label
            )
        coms_t = np.atleast_2d(
            ndimage.center_of_mass(
                this_map_f,
                labels=labels_t,
                index=np.arange(1, np.amax(labels_t)+1)
                )
            )
        numSpots.append(numSpots_t)
        coms.append(coms_t)
        pass

    for i in range(len(pd_hkl_ids)):
        for ispot in range(numSpots[i]):
            if not np.isnan(coms[i][ispot][0]):
                ome_c = eta_ome.omeEdges[0] + (0.5 + coms[i][ispot][0])*del_ome
                eta_c = eta_ome.etaEdges[0] + (0.5 + coms[i][ispot][1])*del_eta
                input_p.append(
                    np.hstack(
                        [hkls[:, pd_hkl_ids[i]],
                         tTh[pd_hkl_ids[i]], eta_c, ome_c]
                    )
                )
                pass
            pass
        pass

    # do the mapping
    start = time.time()
    qfib = None
    if ncpus > 1:
        # multiple process version
        # QUESTION: Need a chunksize?
        pool = mp.Pool(ncpus, discretefiber_init, (params, ))
        qfib = pool.map(discretefiber_reduced, input_p)  # chunksize=chunksize)
        pool.close()
    else:
        # single process version.
        global paramMP
        discretefiber_init(params)  # sets paramMP
        qfib = map(discretefiber_reduced, input_p)
        paramMP = None  # clear paramMP
    elapsed = (time.time() - start)
    logger.info("fiber generation took %.3f seconds", elapsed)
    return np.hstack(qfib)


def discretefiber_init(params):
    global paramMP
    paramMP = params


def discretefiber_reduced(params_in):
    """
    input parameters are [hkl_id, com_ome, com_eta]
    """
    bMat       = paramMP['bMat']
    chi        = paramMP['chi']
    csym       = paramMP['csym']
    fiber_ndiv = paramMP['fiber_ndiv']

    hkl = params_in[:3].reshape(3, 1)

    gVec_s = xfcapi.anglesToGVec(
        np.atleast_2d(params_in[3:]),
        chi=chi,
        ).T

    tmp = mutil.uniqueVectors(
        rot.discreteFiber(
            hkl,
            gVec_s,
            B=bMat,
            ndiv=fiber_ndiv,
            invert=False,
            csym=csym
            )[0]
        )
    return tmp
