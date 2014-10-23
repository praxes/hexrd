import os, subprocess, sys, time

from ConfigParser import SafeConfigParser

import numpy as np
from scipy.sparse import dok_matrix

from hexrd             import coreutil
from hexrd.xrd         import material
from hexrd.xrd.xrdbase import multiprocessing

# from pull_spots_base import d2r, pull_spots_block

have_progBar = False
try:
    from progressbar import ProgressBar, Bar, ETA, ReverseBar
    have_progBar = True
except:
    pass

if __name__ == '__main__':
    cfg_filename = sys.argv[1]
    gp_fileroot  = sys.argv[2]
    
    print "Using cfg file '%s'" % (cfg_filename)

    # pd, reader, detector = coreutil.initialize_experiment(cfg_filename)
    
    parser = SafeConfigParser()
    parser.read(cfg_filename)

    working_dir   = parser.get('base', 'working_dir')

    n_jobs        = parser.get('paint_grid', 'ncpus')
    if n_jobs.strip() == '':
        n_jobs = multiprocessing.cpu_count()
    elif int(n_jobs) == -1:
        n_jobs = multiprocessing.cpu_count() - 1
    elif int(n_jobs) == -2:
        n_jobs = multiprocessing.cpu_count() / 2
    else:
        n_jobs = int(n_jobs)

    print "going to run %d jobs" % (n_jobs)
    
    for job in range(n_jobs):
        print "starting job %d" %job
        subprocess.Popen("python pull_spots_base_2.py %s %d %s" % (cfg_filename, job, gp_fileroot), shell=True)
        
