from __future__ import print_function

#
# mockup for usage of paintGridThis from find-orientation
#
# requires pickle files for inputs and expected outputs
# The pickle files are not in the repo due to their size. In order to
# generate them, the best way is running find-orientatins instrumenting
# PaintGrid such as the code that spawns "paintGridThis" pickles out the
# arguments. Something like:
#
#    # do the mapping
#    import cPickle as pickle
#
#    with open('paintgrid-inputs.pickled', 'wb') as f:
#        pickle.dump(paramMP, f)
#        pickle.dump(quats, f)
#
#    start = time.time()
#    retval = None
#    if multiProcMode:
#        pool = multiprocessing.Pool(nCPUs, paintgrid_init, (paramMP, ))
#        retval = pool.map(paintGridThis, quats.T, chunksize=chunksize)
#    else:
#        retval = map(paintGridThis, quats.T)
#    elapsed = (time.time() - start)
#    logger.info("paintGrid took %.3f seconds", elapsed)
#
#    with open('paintgrid-outputs.pickled', 'wb') as f:
#        pickle.dump(retval, f)
#
# This file contains various alternative implementations with increasing
# performance (in fact, it shows the process of optimization).
# Every optimization is checked against the original results.

import argparse
import cPickle as pickle
from hexrd.cli.main import profile_instrument_all, profile_dump_results
import hexrd.xrd.indexer as indexer
import numpy as np
num = np # original funcion uses num instead of np

from hexrd.xrd import transforms      as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import rotations
import numba


def _meshgrid2d(x, y):
    """
    A special-cased implementation of np.meshgrid, for just
    two arguments. Found to be about 3x faster on some simple
    test arguments.
    """
    x, y = (num.asarray(x), num.asarray(y))
    shape = (len(y), len(x))
    dt = num.result_type(x, y)
    r1, r2 = (num.empty(shape, dt), num.empty(shape, dt))
    r1[...] = x[num.newaxis, :]
    r2[...] = y[:, num.newaxis]
    return (r1, r2)

def paintGridThis_original(quat):
    """
    """
    # Unpack common parameters into locals
    symHKLs = paramMP['symHKLs']
    symHKLs_ix = paramMP['symHKLs_ix']
    wavelength = paramMP['wavelength']
    hklList = paramMP['hklList']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeTol = paramMP['omeTol']
    omePeriod = paramMP['omePeriod']
    omeIndices = paramMP['omeIndices']
    omeEdges = paramMP['omeEdges']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaTol = paramMP['etaTol']
    etaIndices = paramMP['etaIndices']
    etaEdges = paramMP['etaEdges']
    etaOmeMaps = paramMP['etaOmeMaps']
    bMat = paramMP['bMat']
    threshold = paramMP['threshold']

    # need this for proper index generation

    delOmeSign = num.sign(omeEdges[1] - omeEdges[0])

    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])

    dpix_ome = round(omeTol / del_ome)
    dpix_eta = round(etaTol / del_eta)

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    nHKLs = len(symHKLs_ix) - 1

    rMat = rotations.rotMatOfQuat(quat)

    nPredRefl = 0
    nMeasRefl = 0
    reflInfoList = []
    dummySpotInfo = num.nan * num.ones(3)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(
        symHKLs, 0., rMat, bMat, wavelength)
    # Interleave the two produced oang solutions to simplify later processing
    oangs = num.empty((len(symHKLs)*2, 3), dtype=oangs_pair[0].dtype)
    oangs[0::2] = oangs_pair[0]
    oangs[1::2] = oangs_pair[1]

    # Map all of the angles at once
    oangs[:, 1] = xf.mapAngle(oangs[:, 1])
    oangs[:, 2] = xf.mapAngle(oangs[:, 2], omePeriod)

    # Create a mask of the good ones
    oangMask = num.logical_and(
                ~num.isnan(oangs[:, 0]),
            num.logical_and(
                xf.validateAngleRanges(oangs[:, 1], etaMin, etaMax),
                xf.validateAngleRanges(oangs[:, 2], omeMin, omeMax)))

    hklCounterP = 0 # running count of excpected (predicted) HKLs
    hklCounterM = 0 # running count of "hit" HKLs
    for iHKL in range(nHKLs):
        start, stop = symHKLs_ix[iHKL:iHKL+2]
        start, stop = (2*start, 2*stop)

        angList = oangs[start:stop]
        angMask = oangMask[start:stop]

        allAngs_m = angList[angMask, :]
        if len(allAngs_m) > 0:
            # not output # # duplicate HKLs
            # not output # allHKLs_m = num.vstack(
            #     [these_hkls, these_hkls]
            #     )[angMask, :]

            culledTTh  = allAngs_m[:, 0]
            culledEta  = allAngs_m[:, 1]
            culledOme  = allAngs_m[:, 2]
            # not output # culledHKLs = allHKLs_m.T

            nThisPredRefl = len(culledTTh)
            hklCounterP += nThisPredRefl
            for iTTh in range(nThisPredRefl):
                culledEtaIdx = num.where(etaEdges - culledEta[iTTh] > 0)[0]
                if len(culledEtaIdx) > 0:
                    culledEtaIdx = culledEtaIdx[0] - 1
                    if culledEtaIdx < 0:
                        culledEtaIdx = None
                else:
                    culledEtaIdx = None
                culledOmeIdx = num.where(omeEdges - culledOme[iTTh] > 0)[0]
                if len(culledOmeIdx) > 0:
                    if delOmeSign > 0:
                        culledOmeIdx = culledOmeIdx[0] - 1
                    else:
                        culledOmeIdx = culledOmeIdx[-1]
                    if culledOmeIdx < 0:
                        culledOmeIdx = None
                else:
                    culledOmeIdx = None

                if culledEtaIdx is not None and culledOmeIdx is not None:
                    if dpix_ome > 0 or dpix_eta > 0:
                        i_dil, j_dil = _meshgrid2d(
                            num.arange(-dpix_ome, dpix_ome + 1),
                            num.arange(-dpix_eta, dpix_eta + 1)
                            )
                        i_sup = omeIndices[culledOmeIdx] + num.array(
                            [i_dil.flatten()], dtype=int
                            )
                        j_sup = etaIndices[culledEtaIdx] + num.array(
                            [j_dil.flatten()], dtype=int
                            )
                        # catch shit that falls off detector...
                        # ...maybe make this fancy enough to wrap at 2pi?
                        i_max, j_max = etaOmeMaps[iHKL].shape
                        idx_mask = num.logical_and(
                            num.logical_and(i_sup >= 0, i_sup < i_max),
                            num.logical_and(j_sup >= 0, j_sup < j_max)
                            )
                        pixelVal = etaOmeMaps[iHKL][
                            i_sup[idx_mask], j_sup[idx_mask]
                            ]
                    else:
                        pixelVal = etaOmeMaps[iHKL][
                            omeIndices[culledOmeIdx], etaIndices[culledEtaIdx]
                            ]
                    isHit = num.any(pixelVal >= threshold[iHKL])
                    if isHit:
                        hklCounterM += 1
                # if debug:
                #     print "hkl %d -->\t%d\t%d\t%d\t" % (
                #         iHKL, culledHKLs[0, iTTh], culledHKLs[1, iTTh],
                #         culledHKLs[2, iTTh]
                #         ) \
                #         + "isHit=%d\tpixel value: %g\n" % (isHit, pixelVal) \
                #         + "eta index: %d,%d\tetaP: %g\n" % (
                #         culledEtaIdx, etaIndices[culledEtaIdx],
                #         r2d*culledEta[iTTh]
                #         ) \
                #         + "ome index: %d,%d\tomeP: %g\n" % (
                #         culledOmeIdx, omeIndices[culledOmeIdx],
                #         r2d*culledOme[iTTh]
                #         )
                #
                # close conditional on valid reflections
            # close loop on signed reflections
        # close loop on HKL
    if hklCounterP == 0:
        retval = 0.
    else:
        retval = float(hklCounterM) / float(hklCounterP)
    return retval

def _cull(edges, val):
    """search val in the sets of ranges defined by 'edges'
    'edges' is a sorted array containing the points defining extremes of angles.
            that defines a partition.

    returns the index of the last value in 'edges' that is less than 'val'. If
            'val' is out of the range defined by 'edges', returns None.
    """
    idx = num.where(edges - val > 0)[0]
    if len(idx) > 0:
        idx = idx[0] - 1
        if idx < 0:
            idx = None
    else:
        idx = None

    return idx

def paintGridThis_refactor_1(quat):
    """
    """
    # Unpack common parameters into locals
    symHKLs = paramMP['symHKLs']
    symHKLs_ix = paramMP['symHKLs_ix']
    wavelength = paramMP['wavelength']
    hklList = paramMP['hklList']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeTol = paramMP['omeTol']
    omePeriod = paramMP['omePeriod']
    omeIndices = paramMP['omeIndices']
    omeEdges = paramMP['omeEdges']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaTol = paramMP['etaTol']
    etaIndices = paramMP['etaIndices']
    etaEdges = paramMP['etaEdges']
    etaOmeMaps = paramMP['etaOmeMaps']
    bMat = paramMP['bMat']
    threshold = paramMP['threshold']

    # need this for proper index generation

    delOmeSign = num.sign(omeEdges[1] - omeEdges[0])

    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])

    dpix_ome = round(omeTol / del_ome)
    dpix_eta = round(etaTol / del_eta)


    if dpix_ome > 0 or dpix_eta > 0:
        i_dil, j_dil = _meshgrid2d(
            num.arange(-dpix_ome, dpix_ome + 1),
            num.arange(-dpix_eta, dpix_eta + 1)
        )
        i_dil = num.array([i_dil.flatten()], dtype=int)
        j_dil = num.array([j_dil.flatten()], dtype=int)

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    nHKLs = len(symHKLs_ix) - 1

    rMat = rotations.rotMatOfQuat(quat)

    nPredRefl = 0
    nMeasRefl = 0
    reflInfoList = []
    dummySpotInfo = num.nan * num.ones(3)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(symHKLs, 0., rMat, bMat, wavelength)
    # Interleave the two produced oang solutions to simplify later processing
    oangs = num.empty((len(symHKLs)*2, 3), dtype=oangs_pair[0].dtype)
    oangs[0::2] = oangs_pair[0]
    oangs[1::2] = oangs_pair[1]

    # Map all of the angles at once
    oangs[:, 1] = xf.mapAngle(oangs[:, 1])
    oangs[:, 2] = xf.mapAngle(oangs[:, 2], omePeriod)

    # Create a mask of the good ones
    oangMask = num.logical_and(
                ~num.isnan(oangs[:, 0]),
            num.logical_and(
                xf.validateAngleRanges(oangs[:, 1], etaMin, etaMax),
                xf.validateAngleRanges(oangs[:, 2], omeMin, omeMax)))

    hklCounterP = 0 # running count of expected (predicted) HKLs
    hklCounterM = 0 # running count of "hit" HKLs
    for iHKL in range(nHKLs):
        start, stop = 2*symHKLs_ix[iHKL], 2*symHKLs_ix[iHKL+1]

        angList = oangs[start:stop]
        angMask = oangMask[start:stop]

        allAngs_m = angList[angMask, :]
        if len(allAngs_m) > 0:
            # not output # # duplicate HKLs
            # not output # allHKLs_m = num.vstack(
            #     [these_hkls, these_hkls]
            #     )[angMask, :]

            culledTTh  = allAngs_m[:, 0]
            culledEta  = allAngs_m[:, 1]
            culledOme  = allAngs_m[:, 2]
            # not output # culledHKLs = allHKLs_m.T

            nThisPredRefl = len(culledTTh)
            hklCounterP += nThisPredRefl
            for iTTh in range(nThisPredRefl):
                # find the first
                culledEtaIdx = _cull(etaEdges, culledEta[iTTh])
                culledOmeIdx = _cull(omeEdges, culledOme[iTTh])

                if culledEtaIdx is not None and culledOmeIdx is not None:
                    # got a result
                    if dpix_ome > 0 or dpix_eta > 0:
                        i_sup = omeIndices[culledOmeIdx] + i_dil
                        j_sup = etaIndices[culledEtaIdx] + j_dil
                        # catch shit that falls off detector...
                        # ...maybe make this fancy enough to wrap at 2pi?
                        i_max, j_max = etaOmeMaps[iHKL].shape
                        idx_mask = num.logical_and(
                            num.logical_and(i_sup >= 0, i_sup < i_max),
                            num.logical_and(j_sup >= 0, j_sup < j_max)
                            )
                        pixelVal = etaOmeMaps[iHKL][
                            i_sup[idx_mask], j_sup[idx_mask]
                            ]
                    else:
                        pixelVal = etaOmeMaps[iHKL][
                            omeIndices[culledOmeIdx], etaIndices[culledEtaIdx]
                            ]
                    isHit = num.any(pixelVal >= threshold[iHKL])
                    if isHit:
                        hklCounterM += 1
                # if debug:
                #     print "hkl %d -->\t%d\t%d\t%d\t" % (
                #         iHKL, culledHKLs[0, iTTh], culledHKLs[1, iTTh],
                #         culledHKLs[2, iTTh]
                #         ) \
                #         + "isHit=%d\tpixel value: %g\n" % (isHit, pixelVal) \
                #         + "eta index: %d,%d\tetaP: %g\n" % (
                #         culledEtaIdx, etaIndices[culledEtaIdx],
                #         r2d*culledEta[iTTh]
                #         ) \
                #         + "ome index: %d,%d\tomeP: %g\n" % (
                #         culledOmeIdx, omeIndices[culledOmeIdx],
                #         r2d*culledOme[iTTh]
                #         )
                #
                # close conditional on valid reflections
            # close loop on signed reflections
        # close loop on HKL
    if hklCounterP == 0:
        retval = 0.
    else:
        retval = float(hklCounterM) / float(hklCounterP)
    return retval



#
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def _mk_dil(dpix_ome, dpix_eta):
    if dpix_ome > 0 or dpix_eta > 0:
        range_ome = num.arange(-dpix_ome, dpix_ome + 1)
        range_eta = num.arange(-dpix_eta, dpix_eta + 1)
        i_dil, j_dil = _meshgrid2d(range_ome, range_eta)
        i_dil = num.atleast_2d(i_dil.flatten())
        j_dil = num.atleast_2d(j_dil.flatten())
        return i_dil, j_dil
    else:
        return None, None

def _count_hits(angList, hkl_idx, etaOmeMaps, etaEdges, omeEdges,
                etaIndices, omeIndices, dpix_eta, dpix_ome, threshold):
    # pre-generate the indexing for dpix_eta dpix_ome handling
    i_dil, j_dil = _mk_dil(dpix_ome, dpix_eta)

    predicted = len(angList)
    hits = 0
    for curr_ang in range(predicted):
        culledEtaIdx = _cull(etaEdges, angList[curr_ang, 1])

        if culledEtaIdx is None:
            continue

        culledOmeIdx = _cull(omeEdges, angList[curr_ang, 2])
        if culledOmeIdx is None:
            continue

        iHKL = hkl_idx[curr_ang]
        # got a result
        if dpix_ome > 0 or dpix_eta > 0:
            # check if there is any hit in the grid defined by
            # i_dil and j_dil (i-dilatation, j-dilatation)
            i_sup = omeIndices[culledOmeIdx] + i_dil
            j_sup = etaIndices[culledEtaIdx] + j_dil
            # catch shit that falls off detector...
            # ...maybe make this fancy enough to wrap at 2pi?
            i_max, j_max = etaOmeMaps[iHKL].shape
            idx_mask = num.logical_and(
                num.logical_and(i_sup >= 0, i_sup < i_max),
                num.logical_and(j_sup >= 0, j_sup < j_max)
            )
            i_sup = i_sup[idx_mask]
            j_sup = j_sup[idx_mask]
            pixelVal = etaOmeMaps[iHKL][i_sup, j_sup]
        else:
            pixelVal = etaOmeMaps[iHKL][
                omeIndices[culledOmeIdx],
                etaIndices[culledEtaIdx]
            ]
        isHit = num.any(pixelVal >= threshold[iHKL])
        if isHit:
            hits += 1

    return hits


def paintGridThis_refactor_2(quat):
    """
    """
    # Unpack common parameters into a Params object
    #params = Params(**paramMP)

    symHKLs = paramMP['symHKLs'] # the HKLs
    symHKLs_ix = paramMP['symHKLs_ix'] # index partitioning of symHKLs
    bMat = paramMP['bMat']
    wavelength = paramMP['wavelength']
    #hklList = paramMP['hklList'] ### *UNUSED* ###
    omeEdges = paramMP['omeEdges']
    omeTol = paramMP['omeTol'] # used once
    omePeriod = paramMP['omePeriod']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeIndices = paramMP['omeIndices']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaEdges = paramMP['etaEdges']
    etaTol = paramMP['etaTol'] # used once
    etaIndices = paramMP['etaIndices']
    etaOmeMaps = paramMP['etaOmeMaps']
    threshold = paramMP['threshold']

    # dpix_ome and dpix_eta are the number of pixels for the tolerance in
    # ome/eta. Maybe we should compute this per run instead of per-quaternion
    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])
    dpix_ome = int(round(omeTol / del_ome))
    dpix_eta = int(round(etaTol / del_eta))

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    # get the equivalent rotation of the quaternion in matrix form (as expected
    # by oscillAnglesOfHKLs

    rMat = rotations.rotMatOfQuat(quat)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(symHKLs, 0., rMat, bMat, wavelength)

    # Interleave the two produced oang solutions to simplify later processing
    oangs = num.empty((len(symHKLs)*2, 3), dtype=oangs_pair[0].dtype)
    oangs[0::2] = oangs_pair[0]
    oangs[1::2] = oangs_pair[1]

    # Map all of the angles at once
    oangs[:, 1] = xf.mapAngle(oangs[:, 1])
    oangs[:, 2] = xf.mapAngle(oangs[:, 2], omePeriod)

    # Create a mask of the good ones
    oangMask = num.logical_and(~num.isnan(oangs[:, 0]),
                               num.logical_and(
                                   xf.validateAngleRanges(oangs[:, 1],
                                                          etaMin, etaMax),
                                   xf.validateAngleRanges(oangs[:, 2],
                                                          omeMin, omeMax)))

    symHKLs_ix = symHKLs_ix*2
    hkl_idx = num.empty((symHKLs_ix[-1],), dtype=int)
    start = symHKLs_ix[0]
    idx=0
    for end in symHKLs_ix[1:]:
        hkl_idx[start:end] = idx
        start = end
        idx+=1

    angList = oangs[oangMask,:]
    hkl_idx = hkl_idx[oangMask]

    if len(angList) > 0:
        hits = _count_hits(angList, hkl_idx, etaOmeMaps, etaEdges, omeEdges,
                           etaIndices, omeIndices, dpix_eta, dpix_ome,
                           threshold)
        retval = float(hits) / float(len(angList))
    else:
        retval = 0

    return retval



def _count_hits_2(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                  etaIndices, omeIndices, dpix_eta, dpix_ome, threshold):
    # pre-generate the indexing for dpix_eta dpix_ome handling
    i_dil, j_dil = _mk_dil(dpix_ome, dpix_eta)

    predicted = len(hkl_idx)
    hits = 0

    for curr_ang in range(predicted):
        culledEtaIdx = eta_idx[curr_ang]
        culledOmeIdx = ome_idx[curr_ang]
        iHKL = hkl_idx[curr_ang]
        # got a result
        if dpix_ome > 0 or dpix_eta > 0:
            # check if there is any hit in the grid defined by
            # i_dil and j_dil (i-dilatation, j-dilatation)
            i_sup = omeIndices[culledOmeIdx] + i_dil
            j_sup = etaIndices[culledEtaIdx] + j_dil
            # catch shit that falls off detector...
            # ...maybe make this fancy enough to wrap at 2pi?
            i_max, j_max = etaOmeMaps[iHKL].shape
            idx_mask = num.logical_and(
                num.logical_and(i_sup >= 0, i_sup < i_max),
                num.logical_and(j_sup >= 0, j_sup < j_max)
            )
            i_sup = i_sup[idx_mask]
            j_sup = j_sup[idx_mask]
            pixelVal = etaOmeMaps[iHKL][i_sup, j_sup]
        else:
            pixelVal = etaOmeMaps[iHKL][
                omeIndices[culledOmeIdx],
                etaIndices[culledEtaIdx]
            ]
        isHit = num.any(pixelVal >= threshold[iHKL])
        if isHit:
            hits += 1

    return hits


def _filter_angs(angs_0, angs_1, symHKLs_ix, etaEdges, etaMin, etaMax,
                 omeEdges, omeMin, omeMax, omePeriod):
    """
    bakes data in a way that invalid (nan or out-of-bound) is discarded.
    returns:
      - hkl_idx, array of associated hkl indices
      - eta_idx, array of associated eta indices of predicted
      - ome_idx, array of associated ome indices of predicted
    """
    # Interleave the two produced oang solutions to simplify later processing
    oangs = num.empty((len(angs_0)*2, 3), dtype=angs_0.dtype)
    oangs[0::2] = angs_0
    oangs[1::2] = angs_1

    # Map all of the angles at once
    oangs[:, 1] = xf.mapAngle(oangs[:, 1])
    oangs[:, 2] = xf.mapAngle(oangs[:, 2], omePeriod)

    symHKLs_ix = symHKLs_ix*2
    hkl_idx = num.empty((symHKLs_ix[-1],), dtype=int)
    start = symHKLs_ix[0]
    idx=0
    for end in symHKLs_ix[1:]:
        hkl_idx[start:end] = idx
        start = end
        idx+=1

    # using "right" side to make sure we always get an index *past* the value
    # if it happens to be equal. That is... we search the index of the first
    # value that is "greater than" rather than "greater or equal"
    culled_eta_indices = num.searchsorted(etaEdges, oangs[:, 1], side='right')
    culled_ome_indices = num.searchsorted(omeEdges, oangs[:, 2], side='right')
    valid_eta = xf.validateAngleRanges(oangs[:, 1], etaMin, etaMax)
    valid_ome = xf.validateAngleRanges(oangs[:, 2], omeMin, omeMax)
    # Create a mask of the good ones
    valid = ~num.isnan(oangs[:, 0]) # tth not NaN
    valid = num.logical_and(valid, valid_eta)
    valid = num.logical_and(valid, valid_ome)
    valid = num.logical_and(valid, culled_eta_indices > 0)
    valid = num.logical_and(valid, culled_eta_indices < len(etaEdges))
    valid = num.logical_and(valid, culled_ome_indices > 0)
    valid = num.logical_and(valid, culled_ome_indices < len(omeEdges))

    hkl_idx = hkl_idx[valid]
    eta_idx = culled_eta_indices[valid] - 1
    ome_idx = culled_ome_indices[valid] - 1

    return hkl_idx, eta_idx, ome_idx


def paintGridThis_refactor_3(quat):
    """
    """
    # Unpack common parameters into a Params object
    #params = Params(**paramMP)

    symHKLs = paramMP['symHKLs'] # the HKLs
    symHKLs_ix = paramMP['symHKLs_ix'] # index partitioning of symHKLs
    bMat = paramMP['bMat']
    wavelength = paramMP['wavelength']
    #hklList = paramMP['hklList'] ### *UNUSED* ###
    omeEdges = paramMP['omeEdges']
    omeTol = paramMP['omeTol'] # used once
    omePeriod = paramMP['omePeriod']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeIndices = paramMP['omeIndices']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaEdges = paramMP['etaEdges']
    etaTol = paramMP['etaTol'] # used once
    etaIndices = paramMP['etaIndices']
    etaOmeMaps = paramMP['etaOmeMaps']
    threshold = paramMP['threshold']

    # dpix_ome and dpix_eta are the number of pixels for the tolerance in
    # ome/eta. Maybe we should compute this per run instead of per-quaternion
    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])
    dpix_ome = int(round(omeTol / del_ome))
    dpix_eta = int(round(etaTol / del_eta))

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    # get the equivalent rotation of the quaternion in matrix form (as expected
    # by oscillAnglesOfHKLs

    rMat = rotations.rotMatOfQuat(quat)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(symHKLs, 0., rMat, bMat, wavelength)
    hkl_idx, eta_idx, ome_idx = _filter_angs(oangs_pair[0],
                                             oangs_pair[1],
                                             symHKLs_ix,
                                             etaEdges, etaMin, etaMax,
                                             omeEdges, omeMin, omeMax, omePeriod)

    if len(hkl_idx) > 0:
        hits = _count_hits_2(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                             etaIndices, omeIndices, dpix_eta, dpix_ome,
                             threshold)
        retval = float(hits) / float(len(hkl_idx))
    else:
        retval = 0

    return retval


################################################################################
# In the following version, try to use a numba function to count hits

@numba.jit
def check_dilated(eta, ome, dpix_eta, dpix_ome, etaOmeMap, threshold):
    i_max, j_max = etaOmeMap.shape
    for i in range(max(ome - dpix_ome, 0), min(ome + dpix_ome + 1, i_max)):
        for j in range(max(eta - dpix_eta, 0), min(eta + dpix_eta + 1, j_max)):
            if etaOmeMap[i,j] > threshold:
                return 1
    return 0

def _count_hits_3(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                        etaIndices, omeIndices, dpix_eta, dpix_ome, threshold):
    """
    for every eta, ome, hkl check if there is a sample that surpasses the
    threshold in the eta ome map.
    """
    predicted = len(hkl_idx)
    hits = 0

    for curr_ang in range(predicted):
        culledEtaIdx = eta_idx[curr_ang]
        culledOmeIdx = ome_idx[curr_ang]
        iHKL = hkl_idx[curr_ang]
        # got a result
        if dpix_ome > 0 or dpix_eta > 0:
            eta = etaIndices[culledEtaIdx]
            ome = omeIndices[culledOmeIdx]
            isHit = check_dilated(eta, ome, dpix_eta, dpix_ome,
                                  etaOmeMaps[iHKL], threshold[iHKL])
            
        else:
            pixelVal = etaOmeMaps[iHKL][
                omeIndices[culledOmeIdx],
                etaIndices[culledEtaIdx]
            ]
            isHit = pixelVal >= threshold[iHKL]
        if isHit:
            hits += 1

    return hits


def paintGridThis_refactor_4(quat):
    """
    """
    # Unpack common parameters into a Params object
    #params = Params(**paramMP)

    symHKLs = paramMP['symHKLs'] # the HKLs
    symHKLs_ix = paramMP['symHKLs_ix'] # index partitioning of symHKLs
    bMat = paramMP['bMat']
    wavelength = paramMP['wavelength']
    #hklList = paramMP['hklList'] ### *UNUSED* ###
    omeEdges = paramMP['omeEdges']
    omeTol = paramMP['omeTol'] # used once
    omePeriod = paramMP['omePeriod']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeIndices = paramMP['omeIndices']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaEdges = paramMP['etaEdges']
    etaTol = paramMP['etaTol'] # used once
    etaIndices = paramMP['etaIndices']
    etaOmeMaps = paramMP['etaOmeMaps']
    threshold = paramMP['threshold']

    # dpix_ome and dpix_eta are the number of pixels for the tolerance in
    # ome/eta. Maybe we should compute this per run instead of per-quaternion
    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])
    dpix_ome = int(round(omeTol / del_ome))
    dpix_eta = int(round(etaTol / del_eta))

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    # get the equivalent rotation of the quaternion in matrix form (as expected
    # by oscillAnglesOfHKLs

    rMat = rotations.rotMatOfQuat(quat)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(symHKLs, 0., rMat, bMat, wavelength)
    hkl_idx, eta_idx, ome_idx = _filter_angs(oangs_pair[0],
                                             oangs_pair[1],
                                             symHKLs_ix,
                                             etaEdges, etaMin, etaMax,
                                             omeEdges, omeMin, omeMax, omePeriod)

    if len(hkl_idx) > 0:
        hits = _count_hits_3(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                             etaIndices, omeIndices, dpix_eta, dpix_ome,
                             threshold)
        retval = float(hits) / float(len(hkl_idx))
    else:
        retval = 0

    return retval

################################################################################
# jit also the loop

@numba.njit
def _count_hits_4(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                  etaIndices, omeIndices, dpix_eta, dpix_ome, threshold):
    """
    for every eta, ome, hkl check if there is a sample that surpasses the
    threshold in the eta ome map.
    """
    predicted = len(hkl_idx)
    hits = 0

    for curr_ang in range(predicted):
        culledEtaIdx = eta_idx[curr_ang]
        culledOmeIdx = ome_idx[curr_ang]
        iHKL = hkl_idx[curr_ang]
        # got a result
        eta = etaIndices[culledEtaIdx]
        ome = omeIndices[culledOmeIdx]
        isHit = check_dilated(eta, ome, dpix_eta, dpix_ome,
                              etaOmeMaps[iHKL], threshold[iHKL])
            
        if isHit:
            hits += 1

    return hits


def paintGridThis_refactor_5(quat):
    """
    """
    # Unpack common parameters into a Params object
    #params = Params(**paramMP)

    symHKLs = paramMP['symHKLs'] # the HKLs
    symHKLs_ix = paramMP['symHKLs_ix'] # index partitioning of symHKLs
    bMat = paramMP['bMat']
    wavelength = paramMP['wavelength']
    #hklList = paramMP['hklList'] ### *UNUSED* ###
    omeEdges = paramMP['omeEdges']
    omeTol = paramMP['omeTol'] # used once
    omePeriod = paramMP['omePeriod']
    omeMin = paramMP['omeMin']
    omeMax = paramMP['omeMax']
    omeIndices = paramMP['omeIndices']
    etaMin = paramMP['etaMin']
    etaMax = paramMP['etaMax']
    etaEdges = paramMP['etaEdges']
    etaTol = paramMP['etaTol'] # used once
    etaIndices = paramMP['etaIndices']
    etaOmeMaps = paramMP['etaOmeMaps']
    threshold = paramMP['threshold']

    # dpix_ome and dpix_eta are the number of pixels for the tolerance in
    # ome/eta. Maybe we should compute this per run instead of per-quaternion
    del_ome = abs(omeEdges[1] - omeEdges[0])
    del_eta = abs(etaEdges[1] - etaEdges[0])
    dpix_ome = int(round(omeTol / del_ome))
    dpix_eta = int(round(etaTol / del_eta))

    debug = False
    if debug:
        print( "using ome, eta dilitations of (%d, %d) pixels" \
              % (dpix_ome, dpix_eta))

    # get the equivalent rotation of the quaternion in matrix form (as expected
    # by oscillAnglesOfHKLs

    rMat = rotations.rotMatOfQuat(quat)

    # Compute the oscillation angles of all the symHKLs at once
    oangs_pair = xfcapi.oscillAnglesOfHKLs(symHKLs, 0., rMat, bMat, wavelength)
    hkl_idx, eta_idx, ome_idx = _filter_angs(oangs_pair[0],
                                             oangs_pair[1],
                                             symHKLs_ix,
                                             etaEdges, etaMin, etaMax,
                                             omeEdges, omeMin, omeMax, omePeriod)

    if len(hkl_idx > 0):
        hits = _count_hits_4(eta_idx, ome_idx, hkl_idx, etaOmeMaps,
                             etaIndices, omeIndices, dpix_eta, dpix_ome,
                             threshold)
        retval = float(hits) / float(len(hkl_idx))
    else:
        retval = 0

    return retval


################################################################################

def checked_run(function, params, expected, name=None):
    if name is None:
        name = function.__name__

    print("running", name)

    results = function(*params)

    results = np.array(results)
    if np.allclose(expected, results):
        print("Results MATCH")
    else:
        print("Results DO NOT MATCH!")
        print("Max error is :", np.abs((expected-results)/expected))


def load_inputs(input_file):
    with open(input_file, "rb") as f:
        paramMP = pickle.load(f)
        indexer.paintgrid_init(paramMP)
        quats = pickle.load(f)

    return paramMP, quats


def load_outputs(output_file):
    with open(output_file, "rb") as f:
        expected = pickle.load(f)
        expected = np.array(expected)

    return expected


def main():
    global paramMP

    p = argparse.ArgumentParser()
    p.add_argument('-p', '--inst-profile',
                   action='append')
    p.add_argument('-i', '--dump-info',
                   action='store_true')
    p.add_argument('-c', '--count',
                   type=int)
    p.add_argument('experiment',
                   help='experiment to mock. ' \
                        'experiment.inputs.pickled and '\
                        'experiment.outputs.pickled must exist')

    args = p.parse_args()

    if args.inst_profile is not None:
        profile_instrument_all(args.inst_profile)

    input_file = args.experiment + '.inputs.pickled'
    output_file = args.experiment + '.outputs.pickled'

    try:
        paramMP, quats = load_inputs(input_file)
    except Exception:
        print('Error loading inputs file: {0}'.format(input_file))
        return 65

    try:
        expected = load_outputs(output_file)
    except Exception:
        print('Error loading outputs file: {0}'.format(output_file))
        return 65


    # limit the number of arguments to make iterating faster
    if args.count is not None:
        COUNT = args.count
        quats = quats[:,:COUNT]
        expected = expected[:COUNT]

    if args.dump_info:
        print(' {0} experiment '.format(args.experiment).center(72, '='))
        for key, val in sorted(paramMP.items()):
            print('{0}:\n{1}\n'.format(key, val))
        print('Experiment for {0} quats.'.format(len(expected)))

    checked_run(hexrd_version, (quats,), expected)
    checked_run(original, (quats,), expected)
    checked_run(refactor_1, (quats,), expected)
    checked_run(refactor_2, (quats,), expected)
    checked_run(refactor_3, (quats,), expected)
    checked_run(refactor_4, (quats,), expected)
    checked_run(refactor_5, (quats,), expected)

    if args.inst_profile is not None:
        profile_dump_results(args.inst_profile)


def hexrd_version(quats):
    """The current paintGridThis in hexrd"""
    return map(indexer.paintGridThis, quats.T)

def original(quats):
    """paintGridThis as was found in the original"""
    return map(paintGridThis_original, quats.T)

def refactor_1(quats):
    """dilation setup extracted out of loop"""
    return map(paintGridThis_refactor_1, quats.T)

def refactor_2(quats):
    """factored out _count_hits to a separate function"""
    return map(paintGridThis_refactor_2, quats.T)

def refactor_3(quats):
    """factored out _filter_angles to a separate function"""
    return map(paintGridThis_refactor_3, quats.T)

def refactor_4(quats):
    """alternative implementation for checking a hit with dilation, numbafied"""
    return map(paintGridThis_refactor_4, quats.T)

def refactor_5(quats):
    """numba jit also the enclosing loop for checking all hits"""
    return map(paintGridThis_refactor_5, quats.T)


if __name__ == '__main__':
    exit(main())
