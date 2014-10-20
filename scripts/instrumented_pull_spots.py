#
#
# Runner of pulls_spots_mp instrumenting for profile
#

import pull_spots_mp as target

import sys, os, time
import functools
from contextlib import contextmanager
import numpy as np

def add_nvtx_instrumentation(nvtx):
    from scipy import sparse
    from scipy import ndimage
    from hexrd.xrd import xrdutil
    from hexrd.xrd import fitting
    from hexrd.xrd import xrdutil
    from hexrd.xrd import distortion
    from hexrd import gridutil
    from numba import dispatcher

    _locals = locals()
    def PROFILE(func_path, color):
        exec '{0} = nvtx.profiled("{0}", color={1})({0})'.format(func_path, color) in globals(), _locals


    # This 3 show a bit of structure:
    # target.Config.__init__ is initialization time that is outside of multiprocessing
    # target.read_frames is the time taken by reading frame data *and* converting to sparse
    # target.process_grain is the "per grain" process, which should be handled by multiproc
    PROFILE('target.Config.__init__', nvtx.colors.red)
    PROFILE('target.read_frames', nvtx.colors.green)
    PROFILE('target.process_grain', nvtx.colors.blue)

    # This shows where numba is compiling the @jit decorated functions
    PROFILE('dispatcher.Overloaded.compile', nvtx.colors.black)

    # some key functions
    PROFILE('xrdutil.pullSpots', nvtx.colors.blue)
    PROFILE('fitting.fitGrain', nvtx.colors.yellow)
    PROFILE('fitting.objFuncFitGrain', nvtx.colors.magenta)
    PROFILE('sparse.coo_matrix', nvtx.colors.yellow)
    PROFILE('sparse.coo.coo_matrix.todense', nvtx.colors.magenta)
    PROFILE('gridutil.computeArea', nvtx.colors.green)
    PROFILE('gridutil.sutherlandHodgman', nvtx.colors.blue)
    PROFILE('xrdutil.simulateGVecs', nvtx.colors.yellow)
    PROFILE('gridutil.cellCentroids', nvtx.colors.cyan)
    PROFILE('gridutil.cellConnectivity', nvtx.colors.red)
    PROFILE('ndimage.label', nvtx.colors.white)
    PROFILE('distortion._ge_41rt_distortion', nvtx.colors.yellow)
    PROFILE('distortion._ge_41rt_inverse_distortion', nvtx.colors.yellow)


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
    

if __name__ == '__main__':
    cfg_filename = sys.argv[1]
    quats_filename = sys.argv[2] if len(sys.argv) >= 3 else None

    with profiling(use_nvtx=True):
        start = time.time()
        print "Using cfg file '%s'" % (cfg_filename)
        config = target.Config(cfg_filename)
        config.multiproc = False # force sequential run
        quats_filename = quats_filename if quats_filename is not None else config.analysis_name+'-quats.out'
    
        quats = np.loadtxt(os.path.join(config.working_dir, quats_filename))

        target.process_all_grains(config, quats)

        elapsed = time.time() - start
        print "\nTotal processing time %.2f seconds" % elapsed
