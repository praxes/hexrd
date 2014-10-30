import os

from .config import Config


class MaterialConfig(Config):


    @property
    def definitions(self):
        temp = self._cfg.get('material:definitions')
        if os.path.exists(temp):
            return temp
        raise IOError(
            '"material:definitions": "%s" does not exist'
            )


    @property
    def active(self):
        return self._cfg.get('material:active')
