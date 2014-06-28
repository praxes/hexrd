import sys, os, time

import shelve, cPickle
import numpy as np

from scipy.optimize import leastsq

from ConfigParser import SafeConfigParser

from hexrd.xrd import transforms      as xf
from hexrd.xrd import transforms_CAPI as xfcapi

# #################################################
# PARAMETERS
d2r = np.pi/180.
r2d = 180./np.pi

bVec_ref = xf.bVec_ref
vInv_ref = xf.vInv_ref
# #################################################

# #################################################
# LEGACY FUCTIONS FOR WORKING WITH OLD HEXRD PAR
def objFun_tVec_d(tvd_xy, rMat_d, beamXYD, det_origin, bHat_l):
    """
    """
    xformed_xy = beamXYD[:2] - det_origin
    tVec_d = np.hstack([tvd_xy, -beamXYD[2]]).T
    n_d    = rMat_d[:, 2]
    bVec_l = (np.dot(n_d, tVec_d) / np.dot(n_d, bHat_l)) * bHat_l
    bVec_d = np.hstack([xformed_xy, 0.]).T
    return np.dot(rMat_d, bVec_d).flatten() + tVec_d.flatten() - bVec_l.flatten()

def tVec_d_from_old_detector_params(old_par, det_origin):
    """
    as loaded the params have 2 columns: val bool
    """
    beamXYD = old_par[:3, 0]
    rMat_d  = xf.makeDetectorRotMat(old_par[3:6, 0])
    args=(rMat_d, beamXYD, det_origin, bVec_ref)
    tvd_xy = leastsq(objFun_tVec_d, -beamXYD[:2], args=args)[0]
    return np.hstack([tvd_xy, -beamXYD[2]]).reshape(3, 1)

def beamXYD_from_tVec_d(rMat_d, tVec_d, bVec_ref, det_origin):
    # calculate beam position
    n_d = np.dot(rMat_d, np.c_[0., 0., 1.].T)
    u = np.dot(n_d.flatten(), tVec_d) / np.dot(n_d.flatten(), bVec_ref)
    det_origin - tVec_d[:2].flatten()
    return np.hstack([det_origin - tVec_d[:2].flatten(), u.flatten()])

def make_old_detector_parfile(results, filename=None):
    rMat_d = xf.makeDetectorRotMat(results['tiltAngles'])
    tVec_d = results['tVec_d'] - results['tVec_s']
    beamXYD = beamXYD_from_tVec_d(rMat_d, tVec_d, bVec_ref, det_origi)
    det_plist = np.zeros(12)
    det_plist[:3]  = beamXYD.flatten()
    det_plist[3:6] = results['tiltAngles']
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
    ws = expt.Experiment(cfgFile=os.path.join(hexrd_root, "hexrd/data/materials.cfg"),
                         matFile=os.path.join(hexrd_root, "hexrd/data/all_materials.cfg"))

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
