import glob
import os

from .config import Config



class ImageSeriesConfig(Config):


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
