import os

from hexrd.instrument import HEDMInstrument
from .config import Config


class InstrumentConfig(Config):

    @property
    def hedm(self):
        return HEDMInstrument(self._cfg)
