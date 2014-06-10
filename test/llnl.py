#!/usr/bin/env python
from numbapro import vectorize, float64, jit, guvectorize, autojit
import numpy as np
from hexrd.xrd import nb_distortion as dFuncs

epsf = np.finfo(float).eps      # ~2.2e-16

#array declarations
Xl = np.array([[1., 0., 0.]], order='C').T     # X in the lab frame

Zl = np.array([[0., 0., 1.]]).T     # Z in the lab frame

bVec_ref = -Zl 
eta_ref = Xl


# distortion for warping detector coords
dFunc_ref   = dFuncs.GE_41RT
dParams_ref = [0., 0., 0., 2., 2., 2]


def detectorXYToGvec(xy_det,
                     rMat_d, rMat_s,
                     tVec_d, tVec_s, tVec_c,
                     distortion=(dFunc_ref, dParams_ref),
                     beamVec=bVec_ref, etaVec=eta_ref):
    """
    Takes a list cartesian (x, y) pairs in the detector coordinates and calculates
    the associated reciprocal lattice (G) vectors and (bragg angle, azimuth) pairs
    with respect to the specified beam and azimth (eta) reference directions

    Required Arguments:
    xy_det -- (n, 2) ndarray or list-like input of n detector (x, y) points
    rMat_d -- (3, 3) ndarray, the COB taking DETECTOR FRAME components to LAB FRAME
    rMat_s -- (3, 3) ndarray, the COB taking SAMPLE FRAME components to LAB FRAME
    tVec_d -- (3, 1) ndarray, the translation vector connecting LAB to DETECTOR
    tVec_s -- (3, 1) ndarray, the translation vector connecting LAB to SAMPLE
    tVec_c -- (3, 1) ndarray, the translation vector connecting SAMPLE to CRYSTAL

    Optional Keyword Arguments:
    beamVec -- (3, 1) mdarray containing the incident beam direction components in the LAB FRAME
    etaVec  -- (3, 1) mdarray containing the reference azimuth direction components in the LAB FRAME

    Outputs:
    (n, 2) ndarray containing the (tTh, eta) pairs associated with each (x, y)
    (3, n) ndarray containing the associated G vector directions in the LAB FRAME
    associated with gVecs
    """
    npts   = len(xy_det)                       # number of input (x, y) pairs
    bHat_l = unitVector(beamVec.reshape(3, 1)) # make sure beam direction is a unit vector
    eHat_l = unitVector(etaVec.reshape(3, 1))  # make sure eta=0 direction is a unit vector

    xy_det = distortion[0](xy_det, distortion[1])

    # form in-plane vectors for detector points list in DETECTOR FRAME
    P2_d = np.hstack([np.atleast_2d(xy_det), np.zeros((npts, 1))]).T

    # in LAB FRAME
    P2_l = np.dot(rMat_d, P2_d) + tVec_d
    P0_l = tVec_s + np.dot(rMat_s, tVec_c)   # origin of CRYSTAL FRAME

    # diffraction unit vector components in LAB FRAME
    dHat_l = unitVector(P2_l - P0_l)

    # ###############################################################
    # generate output

    # DEBUGGING
    assert abs(np.dot(bHat_l.T, eHat_l)) < 1. - sqrt_epsf, "eta ref and beam cannot be parallel!"

    rMat_e = makeEtaFrameRotMat(bHat_l, eHat_l)
    dHat_e = np.dot(rMat_e.T, dHat_l)

    tTh = np.arccos(np.dot(bHat_l.T, dHat_l)).flatten()
    eta = np.arctan2(dHat_e[1, :], dHat_e[0, :]).flatten()

    # get G-vectors by rotating d by 90-theta about b x d (numpy 'cross' works on row vectors)
    n_g = unitVector(np.cross(bHat_l.T, dHat_l.T).T)

    gVec_l = rotate_vecs_about_axis(0.5 * (np.pi - tTh), n_g, dHat_l)

    return (tTh, eta), gVec_l

def nb_detectorXYToGvec(xy_det,
                     rMat_d, rMat_s,
                     tVec_d, tVec_s, tVec_c,
                     distortion=(dFunc_ref, dParams_ref),
                     beamVec=bVec_ref, etaVec=eta_ref):
    """
    Takes a list cartesian (x, y) pairs in the detector coordinates and calculates
    the associated reciprocal lattice (G) vectors and (bragg angle, azimuth) pairs
    with respect to the specified beam and azimth (eta) reference directions

    Required Arguments:
    xy_det -- (n, 2) ndarray or list-like input of n detector (x, y) points
    rMat_d -- (3, 3) ndarray, the COB taking DETECTOR FRAME components to LAB FRAME
    rMat_s -- (3, 3) ndarray, the COB taking SAMPLE FRAME components to LAB FRAME
    tVec_d -- (3, 1) ndarray, the translation vector connecting LAB to DETECTOR
    tVec_s -- (3, 1) ndarray, the translation vector connecting LAB to SAMPLE
    tVec_c -- (3, 1) ndarray, the translation vector connecting SAMPLE to CRYSTAL

    Optional Keyword Arguments:
    beamVec -- (3, 1) mdarray containing the incident beam direction components in the LAB FRAME
    etaVec  -- (3, 1) mdarray containing the reference azimuth direction components in the LAB FRAME

    Outputs:
    (n, 2) ndarray containing the (tTh, eta) pairs associated with each (x, y)
    (3, n) ndarray containing the associated G vector directions in the LAB FRAME
    associated with gVecs
    """
    npts   = len(xy_det)                       # number of input (x, y) pairs

    # make sure beam direction is a unit vector# make sure beam direction is a unit vector
    bHat_l = normalize(beamVec.reshape(3, 1), create_unit_vector(beamVec.reshape(3, 1))) 
    # make sure eta=0 direction is a unit vector
    eHat_l = normalize(etaVec.reshape(3, 1), create_unit_vector(etaVec.reshape(3,1)))  
    xy_det = distortion[0](xy_det, distortion[1])

    # form in-plane vectors for detector points list in DETECTOR FRAME
    P2_d = np.hstack([np.atleast_2d(xy_det), np.zeros((npts, 1))]).T

    # in LAB FRAME
    P2_l = np.dot(rMat_d, P2_d) + tVec_d
    P0_l = tVec_s + np.dot(rMat_s, tVec_c)   # origin of CRYSTAL FRAME

    # diffraction unit vector components in LAB FRAME
    dHat_l = unitVector(P2_l - P0_l)

    # ###############################################################
    # generate output

    # DEBUGGING
    assert abs(np.dot(bHat_l.T, eHat_l)) < 1. - sqrt_epsf, "eta ref and beam cannot be parallel!"

    rMat_e = makeEtaFrameRotMat(bHat_l, eHat_l)
    dHat_e = np.dot(rMat_e.T, dHat_l)

    tTh = np.arccos(np.dot(bHat_l.T, dHat_l)).flatten()
    eta = np.arctan2(dHat_e[1, :], dHat_e[0, :]).flatten()

    # get G-vectors by rotating d by 90-theta about b x d (numpy 'cross' works on row vectors)
    n_g = unitVector(np.cross(bHat_l.T, dHat_l.T).T)

    gVec_l = rotate_vecs_about_axis(0.5 * (np.pi - tTh), n_g, dHat_l)

    return (tTh, eta), gVec_l




def create_unit_vector(a):
    assert a.ndim in [1, 2], "incorrect arg shape; must be 1-d or 2-d, yours is %d-d" % (a.ndim)
    m = a.shape[0]
    n = 1
    nrm = np.tile(np.sqrt(np.sum(np.asarray(a)**2)), (m, n))
    # prevent divide by zero
    zchk = (nrm <= epsf)
    nrm[zchk] = 1.
    return nrm


#@autojit
@vectorize(['f8(f8, f8)'])
def normalize(a, nrm):
    """
    normalize array of column vectors (hstacked, axis = 0)
    """
    return a / nrm


def unitVector(a):
    """
    normalize array of column vectors (hstacked, axis = 0)
    """
    #assert a.ndim in [1, 2], "incorrect arg shape; must be 1-d or 2-d, yours is %d-d" % (a.ndim)

    m = a.shape[0]
    n = 1
    
    #why asarray
    nrm = np.tile(np.sqrt(sum(np.asarray(a)**2)), (m, n))
    
    #nrm = np.tile(np.sqrt(sum(a**2)), (m, n))
    #print type(nrm)
    # prevent divide by zero
    zchk = (nrm <= epsf)
    nrm[zchk] = 1.

    nrma = a / nrm

    return nrma


if __name__ == '__main__':

   # nrm = np.tile(np.sqrt(sum(np.asarray(Zl**2)), (3, 1))
    z = unitVector(Zl)
    print z
    a = normalize(Zl.reshape(3,1), create_unit_vector(Zl.reshape(3, 1)))
    print a
    #create array and then pass off to numba

