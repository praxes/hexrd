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

    
    print "Using cfg file '%s'" % (cfg_filename)

    # pd, reader, detector = coreutil.initialize_experiment(cfg_filename)
    
    parser = SafeConfigParser()
    parser.read(cfg_filename)

    working_dir   = parser.get('base', 'working_dir')
    analysis_name = parser.get('base', 'analysis_name')

    if len(sys.argv) < 3:
        quats_filename = analysis_name+'-quats.out'
    else:
        quats_filename = sys.argv[2]
    
    
    # # material class
    # material_name = parser.get('material', 'material_name')
    # matl = material.loadMaterialList(os.path.join(working_dir, material_name+'.ini'))[0]
    # 
    # # planeData and reader
    # pd = matl.planeData
    # pd.exclusions = np.zeros_like(pd.exclusions, dtype=bool)
    # pd.exclusions = pd.getTTh() >= d2r*15

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
    
    """
    LOAD QUATS AND CHUNK THEM
    """
    quats = np.loadtxt(os.path.join(working_dir, quats_filename))

    min_block = len(quats) / n_jobs
    last_block = min_block + np.remainder(len(quats), n_jobs)

    print "%d jobs have %d orientations, and the last has %d" % (n_jobs-1, min_block, last_block)
    
    start_indices = np.arange(n_jobs)*min_block
    num_qblock    = np.hstack([np.ones(n_jobs-1)*min_block, last_block])
    
    for i in range(n_jobs):
        fbname = os.path.join(working_dir, analysis_name+"_block_%03d-quats.out" %i)
        fid = open(fbname, 'w')
        ii = start_indices[i]
        for j in range(int(num_qblock[i])):
            print >> fid, "%1.12e\t%1.12e\t%1.12e\t%1.12e" % tuple(quats[ii+j, :])
            pass
        fid.close()
        pass

    """
    ####### READER --> FRAME LIST
    """
    ome_start = parser.getfloat('reader', 'ome_start')     # in DEGREES
    ome_delta = parser.getfloat('reader', 'ome_delta')     # in DEGREES
    threshold = parser.getfloat('pull_spots', 'threshold')
    
    # print "reading frames into memory"
    # if have_progBar:
    #     widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    #     pbar = ProgressBar(widgets=widgets, maxval=reader.getNFrames()).start()
    #     pass
    # frame_list = []
    # for i in range(reader.getNFrames()):
    #     if have_progBar:
    #         pbar.update(i)
    #     frame = reader.read()
    #     frame[frame <= threshold] = 0
    #     frame_list.append(dok_matrix(frame))
    # frame_list = np.array(frame_list)
    # reader = [frame_list, [ome_start*d2r, ome_delta*d2r]]
    # if have_progBar:
    #     pbar.finish()
    #     pass

    # pull_spots_block(cfg_filename, 0, pd, reader, detector)
    # thread_list = []
    # for job in range(n_jobs):
    #     thread_list.append(threading.Thread(target=pull_spots_block, args=(cfg_filename, job, pd, reader, detector)))
    #     thread_list[0].start()
    for job in range(n_jobs):
        print "starting job %d" %job
        subprocess.Popen("python pull_spots_base.py %s %d" % (cfg_filename, job), shell=True)
        
