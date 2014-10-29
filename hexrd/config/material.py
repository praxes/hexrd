import os

from .config import Config


class MaterialConfig(Config):


    @property
    def definitions(self):
        return self._cfg.get('material:definitions', path_exists=True)


    @property
    def active(self):
        return self._cfg.get('material:active')
