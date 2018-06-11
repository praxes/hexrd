"""Beam parameters"""
import numpy as np

from hexrd import constants

class Beam(object):

    def __init__(self, energy, vector):
        self._energy = energy
        self._vector = vector

    @property
    def energy(self):
        return self._energy

    @property
    def vector(self):
        return self._vector

    @property
    def wavelength(self):
        return constants.keVToAngstrom(self.energy)


def calc_beam_vec(azim, pola):
    """
    Calculate unit beam propagation vector from
    spherical coordinate spec in DEGREES

    ...MAY CHANGE; THIS IS ALSO LOCATED IN XRDUTIL!
    """
    tht = np.radians(azim)
    phi = np.radians(pola)
    bv = np.r_[
        np.sin(phi)*np.cos(tht),
        np.cos(phi),
        np.sin(phi)*np.sin(tht)]
    return -bv


def calc_angles_from_beam_vec(bvec):
    """
    Return the azimuth and polar angle from a beam
    vector
    """
    bvec = np.atleast_2d(bvec).reshape(3, 1)
    nvec = mutil.unitVector(-bvec)
    azim = float(
        np.degrees(np.arctan2(nvec[2], nvec[0]))
    )
    pola = float(np.degrees(np.arccos(nvec[1])))
    return azim, pola
