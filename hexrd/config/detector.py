import os

from .config import Config


class DetectorConfig(Config):


    @property
    def columns(self):
        return self._cfg.get('detector:pixels:columns')


    @property
    def parameters_old(self):
        return self._cfg.get(
            'detector:parameters_old', default=None, path_exists=True
            )


    @property
    def parameters(self):
        return self._cfg.get('detector:parameters', path_exists=True)


    @property
    def pixel_size(self):
        temp = self._cfg.get('detector:pixels:size')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def rows(self):
        return self._cfg.get('detector:pixels:rows')
