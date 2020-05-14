



#%%

import time
import os
import logging
import numpy as np
import copy

import numba
import argparse
import contextlib
import multiprocessing
import tempfile
import shutil

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi
from hexrd.xrd import xrdutil

from hexrd.xrd import rotations as rot
from hexrd import valunits

from hexrd.xrd.transforms_CAPI import anglesToGVec, \
                                      makeRotMatOfExpMap, makeDetectorRotMat, makeOscillRotMat, \
                                      gvecToDetectorXY, detectorXYToGvec

import yaml
import cPickle as cpl

import scipy.ndimage as img
try:
    import imageio as imgio
except(ImportError):
    from skimage import io as imgio
import matplotlib.pyplot as plt

# ==============================================================================
# %% SOME SCAFFOLDING
# ==============================================================================

class ProcessController(object):
    """This is a 'controller' that provides the necessary hooks to
    track the results of the process as well as to provide clues of
    the progress of the process"""

    def __init__(self, result_handler=None, progress_observer=None, ncpus = 1,
                 chunk_size = 100):
        self.rh = result_handler
        self.po = progress_observer
        self.ncpus = ncpus
        self.chunk_size = chunk_size
        self.limits = {}
        self.timing = []


    # progress handling --------------------------------------------------------

    def start(self, name, count):
        self.po.start(name, count)
        t = time.time()
        self.timing.append((name, count, t))


    def finish(self, name):
        t = time.time()
        self.po.finish()
        entry = self.timing.pop()
        assert name==entry[0]
        total = t - entry[2]
        logging.info("%s took %8.3fs (%8.6fs per item).", entry[0], total, total/entry[1])


    def update(self, value):
        self.po.update(value)

    # result handler -----------------------------------------------------------

    def handle_result(self, key, value):
        logging.debug("handle_result (%(key)s)", locals())
        self.rh.handle_result(key, value)

    # value limitting ----------------------------------------------------------
    def set_limit(self, key, limit_function):
        if key in self.limits:
            logging.warn("Overwritting limit funtion for '%(key)s'", locals())

        self.limits[key] = limit_function

    def limit(self, key, value):
        try:
            value = self.limits[key](value)
        except KeyError:
            pass
        except Exception:
            logging.warn("Could not apply limit to '%(key)s'", locals())

        return value

    # configuration  -----------------------------------------------------------
    def get_process_count(self):
        return self.ncpus

    def get_chunk_size(self):
        return self.chunk_size


def null_progress_observer():
    class NullProgressObserver(object):
        def start(self, name, count):
            pass

        def update(self, value):
            pass

        def finish(self):
            pass

    return NullProgressObserver()


def progressbar_progress_observer():
    from progressbar import ProgressBar, Percentage, Bar

    class ProgressBarProgressObserver(object):
        def start(self, name, count):
            self.pbar = ProgressBar(widgets=[name, Percentage(), Bar()],
                                    maxval=count)
            self.pbar.start()

        def update(self, value):
            self.pbar.update(value)

        def finish(self):
            self.pbar.finish()

    return ProgressBarProgressObserver()


def forgetful_result_handler():
    class ForgetfulResultHandler(object):
        def handle_result(self, key, value):
            pass # do nothing

    return ForgetfulResultHandler()


def saving_result_handler(filename):
    """returns a result handler that saves the resulting arrays into a file
    with name filename"""
    class SavingResultHandler(object):
        def __init__(self, file_name):
            self.filename = file_name
            self.arrays = {}

        def handle_result(self, key, value):
            self.arrays[key] = value

        def __del__(self):
            logging.debug("Writing arrays in %(filename)s", self.__dict__)
            try:
                np.savez_compressed(open(self.filename, "wb"), **self.arrays)
            except IOError:
                logging.error("Failed to write %(filename)s", self.__dict__)

    return SavingResultHandler(filename)


def checking_result_handler(filename):
    """returns a return handler that checks the results against a
    reference file.

    The Check will consider a FAIL either a result not present in the
    reference file (saved as a numpy savez or savez_compressed) or a
    result that differs. It will consider a PARTIAL PASS if the
    reference file has a shorter result, but the existing results
    match. A FULL PASS will happen when all existing results match

    """
    class CheckingResultHandler(object):
        def __init__(self, reference_file):
            """Checks the result against those save in 'reference_file'"""
            logging.info("Loading reference results from '%s'", reference_file)
            self.reference_results = np.load(open(reference_file, 'rb'))

        def handle_result(self, key, value):
            if key in ['experiment', 'image_stack']:
                return #ignore these

            try:
                reference = self.reference_results[key]
            except KeyError as e:
                logging.warning("%(key)s: %(e)s", locals())
                reference = None

            if reference is None:
                msg = "'{0}': No reference result."
                logging.warn(msg.format(key))

            try:
                if key=="confidence":
                    reference = reference.T
                    value = value.T

                check_len = min(len(reference), len(value))
                test_passed = np.allclose(value[:check_len], reference[:check_len])

                if not test_passed:
                    msg = "'{0}': FAIL"
                    logging.warn(msg.format(key))
                    lvl = logging.WARN
                elif len(value) > check_len:
                    msg = "'{0}': PARTIAL PASS"
                    lvl = logging.WARN
                else:
                    msg = "'{0}': FULL PASS"
                    lvl = logging.INFO
                logging.log(lvl, msg.format(key))
            except Exception as e:
                msg = "%(key)s: Failure trying to check the results.\n%(e)s"
                logging.error(msg, locals())

    return CheckingResultHandler(filename)


# ==============================================================================
# %% OPTIMIZED BITS
# ==============================================================================

# Some basic 3d algebra ========================================================
@numba.njit
def _v3_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


@numba.njit
def _m33_v3_multiply(m, v, dst):
    v0 = v[0]; v1 = v[1]; v2 = v[2]
    dst[0] = m[0, 0]*v0 + m[0, 1]*v1 + m[0, 2]*v2
    dst[1] = m[1, 0]*v0 + m[1, 1]*v1 + m[1, 2]*v2
    dst[2] = m[2, 0]*v0 + m[2, 1]*v1 + m[2, 2]*v2

    return dst


@numba.njit
def _v3_normalized(src, dst):
    v0 = src[0]
    v1 = src[1]
    v2 = src[2]
    sqr_norm = v0*v0 + v1*v1 + v2*v2
    inv_norm = 1.0 if sqr_norm == 0.0 else 1./np.sqrt(sqr_norm)

    dst[0] = v0 * inv_norm
    dst[1] = v1 * inv_norm
    dst[2] = v2 * inv_norm

    return dst


@numba.njit
def _make_binary_rot_mat(src, dst):
    v0 = src[0]; v1 = src[1]; v2 = src[2]

    dst[0,0] = 2.0*v0*v0 - 1.0
    dst[0,1] = 2.0*v0*v1
    dst[0,2] = 2.0*v0*v2
    dst[1,0] = 2.0*v1*v0
    dst[1,1] = 2.0*v1*v1 - 1.0
    dst[1,2] = 2.0*v1*v2
    dst[2,0] = 2.0*v2*v0
    dst[2,1] = 2.0*v2*v1
    dst[2,2] = 2.0*v2*v2 - 1.0

    return dst


# code transcribed in numba from transforms module =============================

# This is equivalent to the transform module anglesToGVec, but written in
# numba. This should end in a module to share with other scripts
@numba.njit
def _anglesToGVec(angs, rMat_ss, rMat_c):
    """From a set of angles return them in crystal space"""
    result = np.empty_like(angs)
    for i in range(len(angs)):
        cx = np.cos(0.5*angs[i, 0])
        sx = np.sin(0.5*angs[i, 0])
        cy = np.cos(angs[i,1])
        sy = np.sin(angs[i,1])
        g0 = cx*cy
        g1 = cx*sy
        g2 = sx

        # with g being [cx*xy, cx*sy, sx]
        # result = dot(rMat_c, dot(rMat_ss[i], g))
        t0_0 = rMat_ss[ i, 0, 0]*g0 + rMat_ss[ i, 1, 0]*g1 + rMat_ss[ i, 2, 0]*g2
        t0_1 = rMat_ss[ i, 0, 1]*g0 + rMat_ss[ i, 1, 1]*g1 + rMat_ss[ i, 2, 1]*g2
        t0_2 = rMat_ss[ i, 0, 2]*g0 + rMat_ss[ i, 1, 2]*g1 + rMat_ss[ i, 2, 2]*g2

        result[i, 0] = rMat_c[0, 0]*t0_0 + rMat_c[ 1, 0]*t0_1 + rMat_c[ 2, 0]*t0_2
        result[i, 1] = rMat_c[0, 1]*t0_0 + rMat_c[ 1, 1]*t0_1 + rMat_c[ 2, 1]*t0_2
        result[i, 2] = rMat_c[0, 2]*t0_0 + rMat_c[ 1, 2]*t0_1 + rMat_c[ 2, 2]*t0_2

    return result


# This is equivalent to the transform's module gvecToDetectorXYArray, but written in
# numba.
# As of now, it is not a good replacement as efficient allocation of the temporary
# arrays is not competitive with the stack allocation using in the C version of the
# code (WiP)

# tC varies per coord
# gvec_cs, rSm varies per grain
#
# gvec_cs
beam = xf.bVec_ref[:, 0]
Z_l = xf.Zl[:,0]
@numba.jit()
def _gvec_to_detector_array(vG_sn, rD, rSn, rC, tD, tS, tC):
    """ beamVec is the beam vector: (0, 0, -1) in this case """
    ztol = xrdutil.epsf
    p3_l = np.empty((3,))
    tmp_vec = np.empty((3,))
    vG_l = np.empty((3,))
    tD_l = np.empty((3,))
    norm_vG_s = np.empty((3,))
    norm_beam = np.empty((3,))
    tZ_l = np.empty((3,))
    brMat = np.empty((3,3))
    result = np.empty((len(rSn), 2))

    _v3_normalized(beam, norm_beam)
    _m33_v3_multiply(rD, Z_l, tZ_l)

    for i in xrange(len(rSn)):
        _m33_v3_multiply(rSn[i], tC, p3_l)
        p3_l += tS
        p3_minus_p1_l = tD - p3_l

        num = _v3_dot(tZ_l, p3_minus_p1_l)
        _v3_normalized(vG_sn[i], norm_vG_s)

        _m33_v3_multiply(rC, norm_vG_s, tmp_vec)
        _m33_v3_multiply(rSn[i], tmp_vec, vG_l)

        bDot = -_v3_dot(norm_beam, vG_l)

        if bDot < ztol or bDot > 1.0 - ztol:
            result[i, 0] = np.nan
            result[i, 1] = np.nan
            continue

        _make_binary_rot_mat(vG_l, brMat)
        _m33_v3_multiply(brMat, norm_beam, tD_l)
        denom = _v3_dot(tZ_l, tD_l)

        if denom < ztol:
            result[i, 0] = np.nan
            result[i, 1] = np.nan
            continue

        u = num/denom
        tmp_res = u*tD_l - p3_minus_p1_l
        result[i,0] = _v3_dot(tmp_res, rD[:,0])
        result[i,1] = _v3_dot(tmp_res, rD[:,1])

    return result


@numba.njit
def _quant_and_clip_confidence(coords, angles, image, base, inv_deltas, clip_vals,bshw):
    """quantize and clip the parametric coordinates in coords + angles

    coords - (..., 2) array: input 2d parametric coordinates
    angles - (...) array: additional dimension for coordinates
    base   - (3,) array: base value for quantization (for each dimension)
    inv_deltas - (3,) array: inverse of the quantum size (for each dimension)
    clip_vals - (2,) array: clip size (only applied to coords dimensions)
    bshw - (1,) half width of the beam stop in mm

    clipping is performed on ranges [0, clip_vals[0]] for x and
    [0, clip_vals[1]] for y

    returns an array with the quantized coordinates, with coordinates
    falling outside the clip zone filtered out.

    """
    count = len(coords)

    in_sensor = 0
    matches = 0
    for i in range(count):
        xf = coords[i, 0]
        yf = coords[i, 1]

        xf = np.floor((xf - base[0]) * inv_deltas[0])
        if not xf >= 0.0:
            continue
        if not xf < clip_vals[0]:
            continue

        if not np.abs(yf)>bshw:
            continue

        yf = np.floor((yf - base[1]) * inv_deltas[1])



        if not yf >= 0.0:
            continue
        if not yf < clip_vals[1]:
            continue

        zf = np.floor((angles[i] - base[2]) * inv_deltas[2])

        in_sensor += 1

        x, y, z = int(xf), int(yf), int(zf)

        #x_byte = x // 8
        #x_off = 7 - (x % 8)
        #if image[z, y, x_byte] (1<<x_off):
        if image[z, y, x]:
            matches += 1

    return 0 if in_sensor == 0 else float(matches)/float(in_sensor)


# ==============================================================================
# %% ORIENTATION TESTING
# ==============================================================================
def test_orientations(image_stack, experiment, test_crds, controller,multiprocessing_start_method):
    """grand loop precomputing the grown image stack

    image-stack -- is the image stack to be tested against.

    experiment  -- A bunch of experiment related parameters.

    controller  -- An external object implementing the hooks to notify progress
                   as well as figuring out what to do with results.
    """

    # extract some information needed =========================================
    # number of grains, number of coords (maybe limited by call), projection
    # function to use, chunk size to use if multiprocessing and the number
    # of cpus.
    n_grains = experiment.n_grains
    chunk_size = controller.get_chunk_size()
    ncpus = controller.get_process_count()


    # generate angles =========================================================
    # all_angles will be a list containing arrays for the different angles to
    # use, one entry per grain.
    #
    # Note that the angle generation is driven by the exp_maps in the experiment
    all_angles = evaluate_diffraction_angles(experiment, controller)

    # generate coords =========================================================
    # The grid of coords to use to test
    #test_crds = generate_test_grid(-0.25, 0.25, 101)
    n_coords = controller.limit('coords', len(test_crds))
#
#    # first, perform image dilation ===========================================
#    # perform image dilation (using scikit_image dilation)
#    subprocess = 'dilate image_stack'
#    dilation_shape = np.ones((2*experiment.row_dilation + 1,
#                              2*experiment.col_dilation + 1),
#                             dtype=np.uint8)
#    image_stack_dilated = np.empty_like(image_stack)
#    dilated = np.empty((image_stack.shape[-2], image_stack.shape[-1]<<3),
#                       dtype=np.bool)
#    n_images = len(image_stack)
#    controller.start(subprocess, n_images)
#    for i_image in range(n_images):
#        to_dilate = np.unpackbits(image_stack[i_image], axis=-1)
#        ski_dilation(to_dilate, dilation_shape,
#                     out=dilated)
#        image_stack_dilated[i_image] = np.packbits(dilated, axis=-1)
#        controller.update(i_image+1)
#    controller.finish(subprocess)

    # precompute per-grain stuff ==============================================
    # gVec_cs and rmat_ss can be precomputed, do so.
    subprocess = 'precompute gVec_cs'
    controller.start(subprocess, len(all_angles))
    precomp = []
    for i, angs in enumerate(all_angles):
        rmat_ss = xfcapi.makeOscillRotMatArray(experiment.chi, angs[:,2])
        gvec_cs = _anglesToGVec(angs, rmat_ss, experiment.rMat_c[i])
        precomp.append((gvec_cs, rmat_ss))
    controller.finish(subprocess)

    # grand loop ==============================================================
    # The near field simulation 'grand loop'. Where the bulk of computing is
    # performed. We are looking for a confidence matrix that has a n_grains
    chunks = xrange(0, n_coords, chunk_size)
    subprocess = 'grand_loop'
    controller.start(subprocess, n_coords)
    finished = 0
    ncpus = min(ncpus, len(chunks))

    logging.info('Checking confidence for %d coords, %d grains.',
                 n_coords, n_grains)
    confidence = np.empty((n_grains, n_coords))
    if ncpus > 1:
        global _multiprocessing_start_method
        _multiprocessing_start_method=multiprocessing_start_method
        logging.info('Running multiprocess %d processes (%s)',
                     ncpus, _multiprocessing_start_method)
        with grand_loop_pool(ncpus=ncpus, state=(chunk_size,
                                                 image_stack,
                                                 all_angles, precomp, test_crds,
                                                 experiment)) as pool:
            for rslice, rvalues in pool.imap_unordered(multiproc_inner_loop,
                                                       chunks):
                count = rvalues.shape[1]
                confidence[:, rslice] = rvalues
                finished += count
                controller.update(finished)
    else:
        logging.info('Running in a single process')
        for chunk_start in chunks:
            chunk_stop = min(n_coords, chunk_start+chunk_size)
            rslice, rvalues = _grand_loop_inner(image_stack, all_angles,
                                                precomp, test_crds, experiment,
                                                start=chunk_start,
                                                stop=chunk_stop)
            count = rvalues.shape[1]
            confidence[:, rslice] = rvalues
            finished += count
            controller.update(finished)

    controller.finish(subprocess)
    controller.handle_result("confidence", confidence)

    del _multiprocessing_start_method

    pool.close()

    return confidence


def evaluate_diffraction_angles(experiment, controller=None):
    """Uses simulateGVecs to generate the angles used per each grain.
    returns a list containg one array per grain.

    experiment -- a bag of experiment values, including the grains specs and other
                  required parameters.
    """
    # extract required data from experiment
    exp_maps = experiment.exp_maps
    plane_data = experiment.plane_data
    detector_params = experiment.detector_params
    pixel_size = experiment.pixel_size
    ome_range = experiment.ome_range
    ome_period = experiment.ome_period

    panel_dims_expanded = [(-10, -10), (10, 10)]
    subprocess='evaluate diffraction angles'
    pbar = controller.start(subprocess,
                            len(exp_maps))
    all_angles = []
    ref_gparams = np.array([0., 0., 0., 1., 1., 1., 0., 0., 0.])
    for i, exp_map in enumerate(exp_maps):
        gparams = np.hstack([exp_map, ref_gparams])
        sim_results = xrdutil.simulateGVecs(plane_data,
                                            detector_params,
                                            gparams,
                                            panel_dims=panel_dims_expanded,
                                            pixel_pitch=pixel_size,
                                            ome_range=ome_range,
                                            ome_period=ome_period,
                                            distortion=None)
        all_angles.append(sim_results[2])
        controller.update(i+1)
        pass
    controller.finish(subprocess)

    return all_angles


def _grand_loop_inner(image_stack, angles, precomp,
                      coords, experiment, start=0, stop=None):
    """Actual simulation code for a chunk of data. It will be used both,
    in single processor and multiprocessor cases. Chunking is performed
    on the coords.

    image_stack -- the image stack from the sensors
    angles -- the angles (grains) to test
    coords -- all the coords to test
    precomp -- (gvec_cs, rmat_ss) precomputed for each grain
    experiment -- bag with experiment parameters
    start -- chunk start offset
    stop -- chunk end offset
    """

    t = time.time()
    n_coords = len(coords)
    n_angles = len(angles)

    # experiment geometric layout parameters
    rD = experiment.rMat_d
    rCn = experiment.rMat_c
    tD = experiment.tVec_d[:,0]
    tS = experiment.tVec_s[:,0]

    # experiment panel related configuration
    base = experiment.base
    inv_deltas = experiment.inv_deltas
    clip_vals = experiment.clip_vals
    distortion = experiment.distortion
    bshw=experiment.bsw/2.

    _to_detector = xfcapi.gvecToDetectorXYArray
    #_to_detector = _gvec_to_detector_array
    stop = min(stop, n_coords) if stop is not None else n_coords

    distortion_fn = None
    if distortion is not None and len(distortion > 0):
        distortion_fn, distortion_args = distortion

    acc_detector = 0.0
    acc_distortion = 0.0
    acc_quant_clip = 0.0
    confidence = np.zeros((n_angles, stop-start))
    grains = 0
    crds = 0

    if distortion_fn is None:
        for igrn in xrange(n_angles):
            angs = angles[igrn]; rC = rCn[igrn]
            gvec_cs, rMat_ss = precomp[igrn]
            grains += 1
            for icrd in xrange(start, stop):
                t0 = time.time()
                det_xy = _to_detector(gvec_cs, rD, rMat_ss, rC, tD, tS, coords[icrd])
                t1 = time.time()
                c = _quant_and_clip_confidence(det_xy, angs[:,2], image_stack,
                                               base, inv_deltas, clip_vals,bshw)
                t2 = time.time()
                acc_detector += t1 - t0
                acc_quant_clip += t2 - t1
                crds += 1
                confidence[igrn, icrd - start] = c
    else:
        for igrn in xrange(n_angles):
            angs = angles[igrn]; rC = rCn[igrn]
            gvec_cs, rMat_ss = precomp[igrn]
            grains += 1
            for icrd in xrange(start, stop):
                t0 = time.time()
                tmp_xys = _to_detector(gvec_cs, rD, rMat_ss, rC, tD, tS, coords[icrd])
                t1 = time.time()
                det_xy = distortion_fn(tmp_xys, distortion_args, invert=True)
                t2 = time.time()
                c = _quant_and_clip_confidence(det_xy, angs[:,2], image_stack,
                                               base, inv_deltas, clip_vals,bshw)
                t3 = time.time()
                acc_detector += t1 - t0
                acc_distortion += t2 - t1
                acc_quant_clip += t3 - t2
                crds += 1
                confidence[igrn, icrd - start] = c

    t = time.time() - t
    return slice(start, stop), confidence


def multiproc_inner_loop(chunk):
    """function to use in multiprocessing that computes the simulation over the
    task's alloted chunk of data"""

    chunk_size = _mp_state[0]
    n_coords = len(_mp_state[4])
    chunk_stop = min(n_coords, chunk+chunk_size)
    return _grand_loop_inner(*_mp_state[1:], start=chunk, stop=chunk_stop)


def worker_init(id_state, id_exp):
    """process initialization function. This function is only used when the
    child processes are spawned (instead of forked). When using the fork model
    of multiprocessing the data is just inherited in process memory."""
    import joblib

    global _mp_state
    state = joblib.load(id_state)
    experiment = joblib.load(id_exp)
    _mp_state = state + (experiment,)

@contextlib.contextmanager
def grand_loop_pool(ncpus, state):
    """function that handles the initialization of multiprocessing. It handles
    properly the use of spawned vs forked multiprocessing. The multiprocessing
    can be either 'fork' or 'spawn', with 'spawn' being required in non-fork
    platforms (like Windows) and 'fork' being preferred on fork platforms due
    to its efficiency.
    """
    # state = ( chunk_size,
    #           image_stack,
    #           angles,
    #           precomp,
    #           coords,
    #           experiment )
    global _multiprocessing_start_method
    if _multiprocessing_start_method == 'fork':
        # Use FORK multiprocessing.

        # All read-only data can be inherited in the process. So we "pass" it as
        # a global that the child process will be able to see. At the end of the
        # processing the global is removed.
        global _mp_state
        _mp_state = state
        pool = multiprocessing.Pool(ncpus)
        yield pool
        del (_mp_state)
    else:
        # Use SPAWN multiprocessing.

        # As we can not inherit process data, all the required data is
        # serialized into a temporary directory using joblib. The
        # multiprocessing pool will have the "worker_init" as initialization
        # function that takes the key for the serialized data, which will be
        # used to load the parameter memory into the spawn process (also using
        # joblib). In theory, joblib uses memmap for arrays if they are not
        # compressed, so no compression is used for the bigger arrays.
        import joblib
        tmp_dir = tempfile.mkdtemp(suffix='-nf-grand-loop')
        try:
            # dumb dumping doesn't seem to work very well.. do something ad-hoc
            logging.info('Using "%s" as temporary directory.', tmp_dir)

            id_exp = joblib.dump(state[-1],
                                 os.path.join(tmp_dir,
                                              'grand-loop-experiment.gz'),
                                 compress=True)
            id_state = joblib.dump(state[:-1],
                                   os.path.join(tmp_dir, 'grand-loop-data'))
            pool = multiprocessing.Pool(ncpus, worker_init,
                                        (id_state[0], id_exp[0]))
            yield pool
        finally:
            logging.info('Deleting "%s".', tmp_dir)
            shutil.rmtree(tmp_dir)





#%% Loading Utilities


def gen_trial_exp_data(grain_out_file,det_file,mat_file, x_ray_energy, mat_name, max_tth, comp_thresh, chi2_thresh, misorientation_bnd, \
                       misorientation_spacing,ome_range_deg, nframes, beam_stop_width):

    print('Loading Grain Data.....')
    #gen_grain_data
    ff_data=np.loadtxt(grain_out_file)

    #ff_data=np.atleast_2d(ff_data[2,:])

    exp_maps=ff_data[:,3:6]
    t_vec_ds=ff_data[:,6:9]


    #
    completeness=ff_data[:,1]

    chi2=ff_data[:,2]

    n_grains=exp_maps.shape[0]

    rMat_c = rot.rotMatOfExpMap(exp_maps.T)




    cut=np.where(np.logical_and(completeness>comp_thresh,chi2<chi2_thresh))[0]
    exp_maps=exp_maps[cut,:]
    t_vec_ds=t_vec_ds[cut,:]
    chi2=chi2[cut]


    # Add Misorientation
    mis_amt=misorientation_bnd*np.pi/180.
    spacing=misorientation_spacing*np.pi/180.

    mis_steps = int(misorientation_bnd/misorientation_spacing)

    ori_pts = np.arange(-mis_amt, (mis_amt+(spacing*0.999)),spacing)
    num_ori_grid_pts=ori_pts.shape[0]**3
    num_oris=exp_maps.shape[0]


    XsO, YsO, ZsO = np.meshgrid(ori_pts, ori_pts, ori_pts)

    grid0 = np.vstack([XsO.flatten(), YsO.flatten(), ZsO.flatten()]).T


    exp_maps_expanded=np.zeros([num_ori_grid_pts*num_oris,3])
    t_vec_ds_expanded=np.zeros([num_ori_grid_pts*num_oris,3])


    for ii in np.arange(num_oris):
        pts_to_use=np.arange(num_ori_grid_pts)+ii*num_ori_grid_pts
        exp_maps_expanded[pts_to_use,:]=grid0+np.r_[exp_maps[ii,:] ]
        t_vec_ds_expanded[pts_to_use,:]=np.r_[t_vec_ds[ii,:] ]


    exp_maps=exp_maps_expanded
    t_vec_ds=t_vec_ds_expanded

    n_grains=exp_maps.shape[0]

    rMat_c = rot.rotMatOfExpMap(exp_maps.T)


    print('Loading Instrument Data.....')
    instr_cfg = yaml.load(open(det_file, 'r'))

    tiltAngles = instr_cfg['detector']['transform']['tilt_angles']
    tVec_d = np.array(instr_cfg['detector']['transform']['t_vec_d']).reshape(3, 1)
    #tVec_d[0] = -0.05
    chi = instr_cfg['oscillation_stage']['chi']
    tVec_s = np.array(instr_cfg['oscillation_stage']['t_vec_s']).reshape(3, 1)

    rMat_d = makeDetectorRotMat(tiltAngles)
    rMat_s = makeOscillRotMat([chi, 0.])

    pixel_size = instr_cfg['detector']['pixels']['size']

    nrows = instr_cfg['detector']['pixels']['rows']
    ncols = instr_cfg['detector']['pixels']['columns']

#    row_dim = pixel_size[0]*nrows # in mm
#    col_dim = pixel_size[1]*ncols # in mm

    x_col_edges = pixel_size[1]*(np.arange(ncols+1) - 0.5*ncols)
    y_row_edges = pixel_size[0]*(np.arange(nrows+1) - 0.5*nrows)[::-1]

    panel_dims = [(-0.5*ncols*pixel_size[1],
                   -0.5*nrows*pixel_size[0]),
                  ( 0.5*ncols*pixel_size[1],
                    0.5*nrows*pixel_size[0])]

    # a bit overkill, but grab max two-theta from all pixel transforms
    rx, ry = np.meshgrid(x_col_edges, y_row_edges)
    gcrds = detectorXYToGvec(np.vstack([rx.flatten(), ry.flatten()]).T,
                             rMat_d, rMat_s,
                             tVec_d, tVec_s, np.zeros(3))
    pixel_tth = gcrds[0][0]

    detector_params = np.hstack([tiltAngles, tVec_d.flatten(), chi, tVec_s.flatten()])


    ome_period_deg=(ome_range_deg[0][0], (ome_range_deg[0][0]+360.)) #degrees
    ome_step_deg=(ome_range_deg[0][1]-ome_range_deg[0][0])/nframes #degrees


    ome_period = (ome_period_deg[0]*np.pi/180.,ome_period_deg[1]*np.pi/180.)
    ome_range = [(ome_range_deg[0][0]*np.pi/180.,ome_range_deg[0][1]*np.pi/180.)]
    ome_step = ome_step_deg*np.pi/180.



    ome_edges = np.arange(nframes+1)*ome_step+ome_range[0][0]#fixed 2/26/17


    base = np.array([x_col_edges[0],
                     y_row_edges[0],
                     ome_edges[0]])
    deltas = np.array([x_col_edges[1] - x_col_edges[0],
                       y_row_edges[1] - y_row_edges[0],
                       ome_edges[1] - ome_edges[0]])
    inv_deltas = 1.0/deltas
    clip_vals = np.array([ncols, nrows])

    print('Loading Material Data.....')
    #Load Material Data
    materials=cpl.load(open( mat_file, "rb" ))


    check=np.zeros(len(materials))
    for ii in np.arange(len(materials)):
        #print materials[ii].name
        check[ii]=materials[ii].name==mat_name

    mat_used=materials[np.where(check)[0][0]]

    #niti_mart.beamEnergy = valunits.valWUnit("wavelength","ENERGY",61.332,"keV")
    mat_used.beamEnergy = valunits.valWUnit("wavelength","ENERGY",x_ray_energy,"keV")
    mat_used.planeData.exclusions = np.zeros(len(mat_used.planeData.exclusions), dtype=bool)


    if max_tth>0.:
         mat_used.planeData.tThMax = np.amax(np.radians(max_tth))
    else:
        mat_used.planeData.tThMax = np.amax(pixel_tth)

    pd=mat_used.planeData


    print('Final Assembly.....')
    experiment = argparse.Namespace()
    # grains related information
    experiment.n_grains = n_grains # this can be derived from other values...
    experiment.rMat_c = rMat_c # n_grains rotation matrices (one per grain)
    experiment.exp_maps = exp_maps # n_grains exp_maps -angle * rotation axis- (one per grain)

    experiment.plane_data = pd
    experiment.detector_params = detector_params
    experiment.pixel_size = pixel_size
    experiment.ome_range = ome_range
    experiment.ome_period = ome_period
    experiment.x_col_edges = x_col_edges
    experiment.y_row_edges = y_row_edges
    experiment.ome_edges = ome_edges
    experiment.ncols = ncols
    experiment.nrows = nrows
    experiment.nframes = nframes# used only in simulate...
    experiment.rMat_d = rMat_d
    experiment.tVec_d = np.atleast_2d(detector_params[3:6]).T
    experiment.chi = detector_params[6] # note this is used to compute S... why is it needed?
    experiment.tVec_s = np.atleast_2d(detector_params[7:]).T
    experiment.rMat_c = rMat_c
    experiment.distortion = None
    experiment.panel_dims = panel_dims # used only in simulate...
    experiment.base = base
    experiment.inv_deltas = inv_deltas
    experiment.clip_vals = clip_vals
    experiment.bsw = beam_stop_width

    if mis_steps ==0:
        nf_to_ff_id_map = cut
    else:
        nf_to_ff_id_map=np.tile(cut,27*mis_steps)

    return experiment, nf_to_ff_id_map

#%%


def gen_nf_test_grid_tomo(x_dim_pnts, z_dim_pnts, v_bnds, voxel_spacing):

    if v_bnds[0]==v_bnds[1]:
        Xs,Ys,Zs=np.meshgrid(np.arange(x_dim_pnts),v_bnds[0],np.arange(z_dim_pnts))
    else:
        Xs,Ys,Zs=np.meshgrid(np.arange(x_dim_pnts),np.arange(v_bnds[0]+voxel_spacing/2.,v_bnds[1],voxel_spacing),np.arange(z_dim_pnts))
        #note numpy shaping of arrays is goofy, returns(length(y),length(x),length(z))


    Zs=(Zs-(z_dim_pnts/2))*voxel_spacing
    Xs=(Xs-(x_dim_pnts/2))*voxel_spacing


    test_crds = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T
    n_crds = len(test_crds)

    return test_crds, n_crds, Xs, Ys, Zs


#%%

def gen_nf_dark(data_folder,img_nums,num_for_dark,nrows,ncols,dark_type='median',stem='nf_',num_digits=5,ext='.tif'):

    dark_stack=np.zeros([num_for_dark,nrows,ncols])

    print('Loading data for dark generation...')
    for ii in np.arange(num_for_dark):
        print('Image #: ' + str(ii))
        dark_stack[ii,:,:]=imgio.imread(data_folder+'%s'%(stem)+str(img_nums[ii]).zfill(num_digits)+ext)
        #image_stack[ii,:,:]=np.flipud(tmp_img>threshold)

    if dark_type=='median':
        print('making median...')
        dark=np.median(dark_stack,axis=0)
    elif dark_type=='min':
        print('making min...')
        dark=np.min(dark_stack,axis=0)

    return dark


#%%
def gen_nf_image_stack(data_folder,img_nums,dark,num_erosions,num_dilations,ome_dilation_iter,threshold,nrows,ncols,stem='nf_',num_digits=5,ext='.tif'):


    image_stack=np.zeros([img_nums.shape[0],nrows,ncols],dtype=bool)

    print('Loading and Cleaning Images...')
    for ii in np.arange(img_nums.shape[0]):
        print('Image #: ' + str(ii))
        tmp_img=imgio.imread(data_folder+'%s'%(stem)+str(img_nums[ii]).zfill(num_digits)+ext)-dark
        #image procesing
        image_stack[ii,:,:]=img.morphology.binary_erosion(tmp_img>threshold,iterations=num_erosions)
        image_stack[ii,:,:]=img.morphology.binary_dilation(image_stack[ii,:,:],iterations=num_dilations)

    #%A final dilation that includes omega
    print('Final Dilation Including Omega....')
    image_stack=img.morphology.binary_dilation(image_stack,iterations=ome_dilation_iter)

    return image_stack


#%%
def scan_detector_parm(image_stack, experiment,test_crds,controller,parm_to_opt,parm_vector,slice_shape):
    #0-distance
    #1-x center
    #2-xtilt
    #3-ytilt
    #4-ztilt

    multiprocessing_start_method = 'fork' if hasattr(os, 'fork') else 'spawn'

    #current detector parameters, note the value for the actively optimized parameters will be ignored
    distance=experiment.detector_params[5]#mm
    x_cen=experiment.detector_params[3]#mm
    xtilt=experiment.detector_params[0]
    ytilt=experiment.detector_params[1]
    ztilt=experiment.detector_params[2]

    num_parm_pts=len(parm_vector)

    trial_data=np.zeros([num_parm_pts,slice_shape[0],slice_shape[1]])

    tmp_td=copy.copy(experiment.tVec_d)
    for jj in np.arange(num_parm_pts):
        print('cycle %d of %d'%(jj+1,num_parm_pts))


        if parm_to_opt==0:
            tmp_td[2]=parm_vector[jj]
        else:
            tmp_td[2]=distance

        if parm_to_opt==1:
            tmp_td[0]=parm_vector[jj]
        else:
            tmp_td[0]=x_cen

        if  parm_to_opt==2:
            rMat_d_tmp=makeDetectorRotMat([parm_vector[jj],ytilt,ztilt])
        elif parm_to_opt==3:
            rMat_d_tmp=makeDetectorRotMat([xtilt,parm_vector[jj],ztilt])
        elif parm_to_opt==4:
            rMat_d_tmp=makeDetectorRotMat([xtilt,ytilt,parm_vector[jj]])
        else:
            rMat_d_tmp=makeDetectorRotMat([xtilt,ytilt,ztilt])

        experiment.rMat_d = rMat_d_tmp
        experiment.tVec_d = tmp_td



        conf=test_orientations(image_stack, experiment, test_crds,
                      controller,multiprocessing_start_method)


        trial_data[jj]=np.max(conf,axis=0).reshape(slice_shape)

    return trial_data

#%%

def extract_max_grain_map(confidence,grid_shape,binary_recon_bin=None):
    if binary_recon_bin == None:
        binary_recon_bin=np.ones([grid_shape[1],grid_shape[2]])


    conf_squeeze=np.max(confidence,axis=0).reshape(grid_shape)
    grains=np.argmax(confidence,axis=0).reshape(grid_shape)
    out_bounds=np.where(binary_recon_bin==0)
    conf_squeeze[:,out_bounds[0],out_bounds[1]] =-0.001

    return conf_squeeze,grains
#%%

def process_raw_confidence(raw_confidence,vol_shape,tomo_mask=None,id_remap=None):

    print('Compiling Confidence Map...')
    confidence_map=np.max(raw_confidence,axis=0).reshape(vol_shape)
    grain_map=np.argmax(raw_confidence,axis=0).reshape(vol_shape)


    if tomo_mask is not None:
        print('Applying tomography mask...')
        out_bounds=np.where(tomo_mask==0)
        confidence_map[:,out_bounds[0],out_bounds[1]] =-0.001
        grain_map[:,out_bounds[0],out_bounds[1]] =-1


    if id_remap is not None:
        max_grain_no=np.max(grain_map)
        grain_map_copy=copy.copy(grain_map)
        print('Remapping grain ids to ff...')
        for ii in np.arange(max_grain_no):
            this_grain=np.where(grain_map==ii)
            grain_map_copy[this_grain]=id_remap[ii]
        grain_map=grain_map_copy

    return grain_map, confidence_map

#%%

def save_raw_confidence(save_dir,save_stem,raw_confidence,id_remap=None):
    print('Saving raw confidence, might take a while...')
    if id_remap is not None:
        np.savez(save_dir+save_stem+'_raw_confidence.npz',raw_confidence=raw_confidence,id_remap=id_remap)
    else:
        np.savez(save_dir+save_stem+'_raw_confidence.npz',raw_confidence=raw_confidence)
#%%

def save_nf_data(save_dir,save_stem,grain_map,confidence_map,Xs,Ys,Zs,ori_list,id_remap=None):
    print('Saving grain map data...')
    if id_remap is not None:
        np.savez(save_dir+save_stem+'_grain_map_data.npz',grain_map=grain_map,confidence_map=confidence_map,Xs=Xs,Ys=Ys,Zs=Zs,ori_list=ori_list,id_remap=id_remap)
    else:
        np.savez(save_dir+save_stem+'_grain_map_data.npz',grain_map=grain_map,confidence_map=confidence_map,Xs=Xs,Ys=Ys,Zs=Zs,ori_list=ori_list)

#%%

def plot_ori_map(grain_map, confidence_map, exp_maps, layer_no,id_remap=None):

    grains_plot=np.squeeze(grain_map[layer_no,:,:])
    conf_plot=np.squeeze(confidence_map[layer_no,:,:])
    n_grains=len(exp_maps)

    rgb_image=np.zeros([grains_plot.shape[0],grains_plot.shape[1],4], dtype='float32')
    rgb_image[:,:,3]=1.

    for ii in np.arange(n_grains):
        if id_remap is not None:
            this_grain=np.where(np.squeeze(grains_plot)==id_remap[ii])
        else:
            this_grain=np.where(np.squeeze(grains_plot)==ii)
        if np.sum(this_grain[0])>0:

            ori=exp_maps[ii,:]

            #cubic mapping
            rgb_image[this_grain[0],this_grain[1],0]=(ori[0]+(np.pi/4.))/(np.pi/2.)
            rgb_image[this_grain[0],this_grain[1],1]=(ori[1]+(np.pi/4.))/(np.pi/2.)
            rgb_image[this_grain[0],this_grain[1],2]=(ori[2]+(np.pi/4.))/(np.pi/2.)



    plt.imshow(rgb_image,interpolation='none')
    plt.hold(True)
    plt.imshow(conf_plot,vmin=0.0,vmax=1.,interpolation='none',cmap=plt.cm.gray,alpha=0.5)
