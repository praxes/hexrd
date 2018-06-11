"""Beam parameters"""
import numpy as np

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
