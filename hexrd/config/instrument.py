import os

import numpy as np

from hexrd import constants
from hexrd import instrument
from hexrd.xrd.distortion import GE_41RT as dfunc  # !!! UGH, FIXME !!!

from .config import Config


class Instrument(Config):

    def __init__(self, icfg):
        self._hedm = instrument.HEDMInstrument(icfg)

    # Note: instrument is instantiated with a yaml dictionary; use self
    #       to instantiate classes based on this one
    @property
    def hedm(self):
        return self._hedm
    @hedm.setter
    def hedm(self, yml):
        with open(yml, 'r') as f:
            icfg = yaml.safe_load(f)
        self._hedm = instrument.HEDMInstrument(icfg)

    @property
    def detector_dict(self):
        """returns dictionary of detectors"""
        return self.hedm.detectors
