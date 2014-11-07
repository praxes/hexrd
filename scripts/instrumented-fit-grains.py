#
#
# Runner of pulls_spots_mp instrumenting for profile
#

import hexrd.fitgrains as target
from hexrd import config

import sys, os, time
import functools
from contextlib import contextmanager
import numpy as np
from hexrd import USE_NUMBA


def add_nvtx_instrumentation(nvtx, refinement=1):
    from scipy import sparse
    from scipy import ndimage
    from scipy import optimize
    from hexrd.xrd import xrdutil
    from hexrd.xrd import fitting
    from hexrd.xrd import xrdutil
    from hexrd.xrd import distortion
    from hexrd import gridutil
    from hexrd.xrd import transforms_CAPI as xfcapi
    from hexrd.xrd import transforms as xf
    from hexrd.xrd import detector
    import numpy
    if USE_NUMBA:
        from numba import dispatcher

    _locals = locals()
    def PROFILE(func_path, color):
        exec '{0} = nvtx.profiled("{0}", color={1})({0})'.format(func_path, color) in globals(), _locals


    if refinement < 0:
        return

    # This 3 show a bit of structure:
    # config.open is initialization time that is outside of multiprocessing
    # target.get_frames is the time taken by reading frame data *and* converting to sparse
    # target.FitGrainsWorker.loop is the "per grain" process, which should be handled by multiproc
    if USE_NUMBA:
        PROFILE('dispatcher.Overloaded.compile', nvtx.colors.black)
    PROFILE('target.config.open', nvtx.colors.red)
    PROFILE('target.FitGrainsWorker.loop', nvtx.colors.blue)


    if refinement < 1:
        return

    # This shows where numba is compiling the @jit decorated functions
    PROFILE('target.get_frames', nvtx.colors.green)
    PROFILE('xrdutil.pullSpots', nvtx.colors.blue)
    PROFILE('fitting.fitGrain', nvtx.colors.yellow)
    PROFILE('target.extract_ijv', nvtx.colors.magenta)
    PROFILE('numpy.nonzero', nvtx.colors.red)

    if refinement < 2:
        return

    PROFILE('optimize.leastsq', nvtx.colors.blue)
    PROFILE('detector.ReadGE.read', nvtx.colors.cyan)

    # some key functions
    PROFILE('numpy.meshgrid', nvtx.colors.red)
    PROFILE('numpy.arange', nvtx.colors.red)
    PROFILE('xrdutil.simulateGVecs', nvtx.colors.yellow)
    PROFILE('xrdutil._coo_build_window', nvtx.colors.yellow)
    PROFILE('fitting.objFuncFitGrain', nvtx.colors.magenta)
    PROFILE('sparse.coo_matrix', nvtx.colors.yellow)
    PROFILE('sparse.coo.coo_matrix.todense', nvtx.colors.magenta)
    PROFILE('gridutil.cellIndices', nvtx.colors.red)
    PROFILE('gridutil.computeArea', nvtx.colors.red)
    PROFILE('gridutil.sutherlandHodgman', nvtx.colors.red)
    PROFILE('gridutil.cellCentroids', nvtx.colors.red)
    PROFILE('gridutil.cellConnectivity', nvtx.colors.red)
    PROFILE('ndimage.label', nvtx.colors.white)
    PROFILE('distortion._ge_41rt_distortion', nvtx.colors.yellow)
    PROFILE('distortion._ge_41rt_inverse_distortion', nvtx.colors.yellow)
    PROFILE('xfcapi.makeDetectorRotMat', nvtx.colors.blue)
    PROFILE('xfcapi.makeRotMatOfExpMap', nvtx.colors.blue)
    PROFILE('xfcapi.makeOscillRotMat', nvtx.colors.blue)
    PROFILE('xfcapi.gvecToDetectorXY', nvtx.colors.blue)
    PROFILE('xf.anglesToGVec', nvtx.colors.green)
    PROFILE('xf.angularDifference', nvtx.colors.green)


def print_nvtx_profile(nvtx):
    print " STATS ".center(72, '=')
    fmt = "{2:>14}, {1:>8}, {0:<40}"
    print fmt.format("FUNCTION", "CALLS", "TIME")
    fmt = "{2:>14F}, {1:>8}, {0:<40}"
    sorted_by_time = sorted(nvtx.getstats().iteritems(), key=lambda tup: tup[1][1])
    for key, val in sorted_by_time:
        print fmt.format(key, *val)

@contextmanager
def profiling(profile=False, use_nvtx=False):
    if use_nvtx:
        try:
            from numbapro import nvtx
            add_nvtx_instrumentation(nvtx)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print "Could not import nvtx, skipping nvtx profile"
            use_nvtx = False

    if profile:
        import cProfile as profile
        pr = profile.Profile()
        pr.enable()

    yield

    if profile:
        pr.disable()
        pr.print_stats(sort='cumulative')

    if use_nvtx:
        print_nvtx_profile(nvtx)


def usage():
    print "USAGE: {0} [-p] [-n] [-c <max_grains>] <experiment_file>".format(sys.argv[0])

if __name__ == '__main__':
    import getopt
    max_grains = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'pnc:')
    except getopt.GeoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    prof_dict = {}
    for o, val in opts:
        if o == "-p":
            prof_dict["profile"] = True
        elif o == "-n":
            prof_dict["use_nvtx"] = True
        elif o=="-c":
            try:
                max_grains = int(val)
            except ValueError:
                print "invalid grain count"
                usage()
                sys.exit(2)
        else:
            assert False, "unhandled option '{0}'".format(o)

    if len(args) < 1:
        usage()
        sys.exit(2)
    cfg_filename = args[0]

    with profiling(**prof_dict):
        start = time.time()
        print "Using cfg file '%s'" % (cfg_filename)
        config = config.open(cfg_filename)[0]
        config._cfg['multiproc'] = 1 # force sequential run

        target.fit_grains(config, force=True, max_grains=max_grains)

        elapsed = time.time() - start
        print "\nTotal processing time %.2f seconds" % elapsed
