import glob
import os

from .config import Config
from hexrd import imageseries


class ImageSeriesConfig(Config):

    def _open(self):
        self._imser = imageseries.open(self.filename, self.format, **self.args)

    def _meta(self):
        pass # to be done later


    @property
    def filename(self):
        temp = self._cfg.get('image_series:filename')
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        return temp

    @property
    def format(self):
        return self._cfg.get('image_series:format')

    @property
    def args(self):
        return self._cfg.get('image_series:args')

    @property
    def omega(self):
        return OmegaConfig(self._cfg)


class OmegaConfig(Config):

    @property
    def step(self):
        return self._cfg.get('image_series:omega:step')

    @property
    def start(self):
        return self._cfg.get('image_series:omega:start')

    @property
    def stop(self):
        return self._cfg.get('image_series:omega:stop')
