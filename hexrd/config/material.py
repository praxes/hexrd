import os

from .config import Config


class MaterialConfig(Config):


    @property
    def definitions(self):
        try:
            temp = self._cfg['definitions']
            if not os.path.isfile(temp):
                raise IOError(
                    'material:definitions specified at "%s" do not exist' % temp
                    )
            return temp
        except KeyError:
            raise RuntimeError(
                'material:definitions must be specified in the config file'
                )


    @property
    def active(self):
        try:
            return self._cfg['active']
        except KeyError:
            raise RuntimeError(
                'material:active must be specified in the config file'
                )
