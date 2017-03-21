import glob
import os

import numpy as np

from .config import Config
from hexrd import imageseries


class ImageSeriesConfig(Config):

    def __init__(self, cfg):
        super(ImageSeriesConfig, self).__init__(cfg)
        self._imser = None
        self._omseries = None

    @property
    def imageseries(self):
        """return the imageseries without checking for omega metadata"""
        if self._imser is None:
            self._imser = imageseries.open(self.filename, self.format, **self.args)
        plist = self.process.process_list
        if plist:
            self._imser = imageseries.process.ProcessedImageSeries(self._imser, plist)
        return self._imser

    @property
    def omegaseries(self):
        """return the imageseries and ensure it has omega metadata"""
        if self._omseries is None:
            self._omseries = imageseries.omega.OmegaImageSeries(self.imageseries)
        return self._omseries

    # ========== yaml inputs

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
        return self._cfg.get('image_series:args', default={})

    # ========== Other Configs

    @property
    def omega(self):
        return OmegaConfig(self._cfg)

    @property
    def process(self):
        return ProcessConfig(self._cfg)

    @property
    def stop(self):
        return self._cfg.get('image_series:omega:stop', default=None)


class ProcessConfig(Config):

    @property
    def process_list(self):
        plist = []
        dark = self.dark
        if self.dark is not None:
            plist.append(('dark', dark))
        flip = self.flip
        if self.flip is not None:
            plist.append(('flip', flip))

        return plist

    @property
    def flip(self):
        return self._cfg.get('image_series:process:flip', default=None)

    @property
    def dark(self):
        # fixed bug that returned np.load(None)
        fname = self._cfg.get('image_series:process:dark', default=None)
        if fname is not None:
            return np.load(fname)


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
