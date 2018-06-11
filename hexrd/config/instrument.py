import os

from hexrd import instrument
from .config import Config

from hexrd import constants

class Instrument(Config):

    @property
    def hedm(self):
        return instrument.HEDMInstrument(self._cfg)

    def beam(self):
        bcfg = Beam(self._cfg)
        return(instrument.beam(bcfg.energy, bcfg.vector))

class Beam(Config):

    beam_energy_DFLT = 65.351
    beam_vec_DFLT = constants.beam_vec

    BASEKEY = 'beam'

    def get(self, key, **kwargs):
        """get item with given key"""
        return self._cfg.get(':'.join([self.BASEKEY, key]), **kwargs)

    @property
    def energy(self):
        return self.get('energy', default=self.beam_energy_DFLT)

    @property
    def vector(self):
        d = self.get('vector', default=None)
        if d is None:
            return self.beam_vec_DFLT

        az = d['azimuth']
        pa = d['polar_angle']
        return instrument.beam.calc_beam_vec(az, pa)
