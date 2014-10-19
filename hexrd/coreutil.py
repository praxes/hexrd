import sys, os, time

import shelve, cPickle
import numpy as np

from scipy.optimize import leastsq
from scipy.linalg import solve

from ConfigParser import SafeConfigParser

from hexrd.xrd import experiment      as expt
from hexrd.xrd import transforms      as xf
from hexrd.xrd import transforms_CAPI as xfcapi

from hexrd.xrd.detector import ReadGE


# #################################################
# PARAMETERS
d2r = np.pi/180.
r2d = 180./np.pi

bVec_ref = xf.bVec_ref
vInv_ref = xf.vInv_ref
# #################################################

# #################################################
# LEGACY FUCTIONS FOR WORKING WITH OLD HEXRD PAR
def tVec_d_from_old_detector_params(old_par, det_origin):
    """
    calculate tVec_d from old [xc, yc, D] parameter spec

    as loaded the params have 2 columns: [val, bool]
    """
    rMat_d = xf.makeDetectorRotMat(old_par[3:6, 0])
    #
    P2_d   = np.c_[ old_par[0, 0] - det_origin[0],
                    old_par[1, 0] - det_origin[1],
                    0.].T
    P2_l   = np.c_[ 0.,
                    0.,
                   -old_par[2, 0]].T
    return P2_l - np.dot(rMat_d, P2_d)

def old_detector_params_from_new(new_par, det_origin):
    """
    calculate beam position in old parameter spec

    *) may need to consider distortion...
    """
    rMat_d = xf.makeDetectorRotMat(new_par[:3])
    tVec_d = new_par[3:6].reshape(3, 1)
    #
    A = np.eye(3); A[:, :2] = rMat_d[:, :2]
    #
    return solve(A, -tVec_d) + np.vstack([det_origin[0], det_origin[1], 0])

def make_old_detector_parfile(results, det_origin=(204.8, 204.8), filename=None):
    tiltAngles = np.array(results['tiltAngles'])
    rMat_d     = xf.makeDetectorRotMat(tiltAngles)
    tVec_d     = results['tVec_d'] - results['tVec_s']
    #
    beamXYD = old_detector_params_from_new(np.hstack([tiltAngles.flatten(), tVec_d.flatten()]), det_origin)
    #
    det_plist = np.zeros(12)
    det_plist[:3]  = beamXYD.flatten()
    det_plist[3:6] = tiltAngles
    det_plist[6:]  = results['dParams']
    if filename is not None:
        if isinstance(filename, file):
            fid = filename
        elif isinstance(filename, str) or isinstance(filename, unicode):
            fid = open(filename, 'w')
        print >> fid, "# DETECTOR PARAMETERS (from new geometry model fit)"
        print >> fid, "# \n# <class 'hexrd.xrd.detector.DetectorGeomGE'>\n#"
        for i in range(len(det_plist)):
            print >> fid, "%1.8e\t%d" % (det_plist[i], 0)
        fid.close()
    return det_plist

def convert_detector_params(old_par, det_origin, chi=0., tVec_s=np.zeros(3), filename=None):
    """
    takes old ge detector parfile from hexrd and converts to the new 10
    parameters
    """
    tVec_d = tVec_d_from_old_detector_params(old_par, det_origin)
    detector_params = np.hstack([old_par[3:6, 0], tVec_d.flatten(), chi, tVec_s.flatten()])
    if filename is not None:
        if isinstance(filename, file):
            fid = filename
        elif isinstance(filename, str) or isinstance(filename, unicode):
            fid = open(filename, 'w')
        print >> fid, "# DETECTOR PARAMETERS"
        print >> fid, \
            "%1.7e\t# tiltAngles[0]\n" % (detector_params[0]) + \
            "%1.7e\t# tiltAngles[1]\n" % (detector_params[1]) + \
            "%1.7e\t# tiltAngles[2]\n" % (detector_params[2]) + \
            "%1.7e\t# tVec_d[0]    \n" % (detector_params[3]) + \
            "%1.7e\t# tVec_d[1]    \n" % (detector_params[4]) + \
            "%1.7e\t# tVec_d[2]    \n" % (detector_params[5]) + \
            "%1.7e\t# chi          \n" % (detector_params[6]) + \
            "%1.7e\t# tVec_s[0]    \n" % (detector_params[7]) + \
            "%1.7e\t# tVec_s[1]    \n" % (detector_params[8]) + \
            "%1.7e\t# tVec_s[2]    \n" % (detector_params[9])
        fid.close()
    return detector_params

def make_grain_params(quat, tVec_c=np.zeros(3), vInv=vInv_ref, filename=None):
    phi      = 2*np.arccos(quat[0])
    n        = xf.unitVector(quat[1:, :])
    expMap_c = phi*n
    if filename is not None:
        if isinstance(filename, file):
            fid = filename
        elif isinstance(filename, str) or isinstance(filename, unicode):
            fid = open(filename, 'w')
        print >> fid, \
            "%1.7e\t# expMap_c[0]\n" % (expMap_c[0]) + \
            "%1.7e\t# expMap_c[1]\n" % (expMap_c[1]) + \
            "%1.7e\t# expMap_c[2]\n" % (expMap_c[2]) + \
            "%1.7e\t# tVec_c[0]  \n" % (tVec_c[0])   + \
            "%1.7e\t# tVec_c[1]  \n" % (tVec_c[1])   + \
            "%1.7e\t# tVec_c[2]  \n" % (tVec_c[2])   + \
            "%1.7e\t# vInv_s[0]  \n" % (vInv[0])     + \
            "%1.7e\t# vInv_s[1]  \n" % (vInv[1])     + \
            "%1.7e\t# vInv_s[2]  \n" % (vInv[2])     + \
            "%1.7e\t# vInv_s[3]  \n" % (vInv[3])     + \
            "%1.7e\t# vInv_s[4]  \n" % (vInv[4])     + \
            "%1.7e\t# vInv_s[5]  \n" % (vInv[5])
        fid.close()
    return np.hstack([expMap_c.flatten(), tVec_c.flatten(), vInv.flatten()])
# #################################################

# #################################################
#
def initialize_experiment(cfg_file):
    """
    """
    parser = SafeConfigParser()
    parser.read(cfg_file)

    hexrd_root = parser.get('base', 'hexrd_root')

    # make experiment
    print __file__
    ws = expt.Experiment()

    working_dir = parser.get('base', 'working_dir')

    materials_fname = parser.get('material', 'materials_fname')
    material_name   = parser.get('material', 'material_name')
    detector_fname  = parser.get('detector', 'parfile_name')

    # MATERIALS
    ws.loadMaterialList(os.path.join(working_dir, materials_fname))
    ws.activeMaterial = material_name
    print "setting active material to '%s'" % (material_name)

    pd = ws.activeMaterial.planeData

    image_dir = parser.get('reader', 'image_dir')

    # number of files ASSUMING SEQUENTIAL SCAN NUMBERS
    file_start  = parser.getint('reader', 'file_start')
    file_stop   = parser.getint('reader', 'file_stop')
    file_suffix = parser.get('reader', 'file_suffix')
    nfiles      = file_stop - file_start + 1
    zpad_str    = '%0' + parser.get('reader', 'file_zpad') + 'd' # int
    fileInfo = []
    for i in [file_start + i for i in range(nfiles)]:
        if file_suffix == '':
            image_filename = parser.get('reader', 'file_root') + '_' + zpad_str % (i)
        else:
            image_filename = parser.get('reader', 'file_root') + '_' + zpad_str % (i) + '.' + file_suffix
        fileInfo.append( ( os.path.join(image_dir, image_filename), parser.getint('reader', 'nempty') ) )
    ome_start = parser.getfloat('reader', 'ome_start') * d2r
    ome_delta = parser.getfloat('reader', 'ome_delta') * d2r
    darkName  = parser.get('reader', 'dark')
    if darkName.strip() == '':
        dark         = None
        subtractDark = False
    else:
        # dark = os.path.join(image_dir, darkName)
        dark = darkName
        subtractDark = True
    doFlip  = parser.getboolean('reader', 'doFlip')
    flipArg = parser.get('reader', 'flipArg')

    # make frame reader
    reader   = ReadGE(fileInfo, ome_start, ome_delta,
                      subtractDark=subtractDark, dark=dark,
                      doFlip=doFlip, flipArg=flipArg)

    # DETECTOR
    ws.loadDetector(os.path.join(working_dir, detector_fname))

    return pd, reader, ws.detector
