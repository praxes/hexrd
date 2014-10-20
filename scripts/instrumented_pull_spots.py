#
#
# Runner of pulls_spots_mp instrumenting for profile
#

import pull_spots_mp as target

import sys, os, time
import cProfile as profile
import functools
from contextlib import contextmanager


import numpy as np

@contextmanager
def profiling(pr_file=None, use_nvtx=False):
    pr = profile.Profile()
    pr.enable()
    yield
    pr.disable()
    pr.print_stats(sort='cumulative')
    

if __name__ == '__main__':
    cfg_filename = sys.argv[1]
    quats_filename = sys.argv[2] if len(sys.argv) >= 3 else None

    with profiling():
        start = time.time()
        print "Using cfg file '%s'" % (cfg_filename)
        config = target.Config(cfg_filename)
        config.multiproc = False # force sequential run
        quats_filename = quats_filename if quats_filename is not None else config.analysis_name+'-quats.out'
    
        quats = np.loadtxt(os.path.join(config.working_dir, quats_filename))

        target.process_all_grains(config, quats)

        elapsed = time.time() - start
        print "\nTotal processing time %.2f seconds" % elapsed
