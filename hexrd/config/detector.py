import os

from .config import Config


class DetectorConfig(Config):


    @property
    def columns(self):
        return self._cfg.get('detector:pixels:columns')


    @property
    def parameters_old(self):
        temp = self._cfg.get('detector:parameters_old', default=None)
        if temp is None:
            return
        if not os.path.isfile(temp):
            raise IOError(
                'detector:parameters_old specified at "%s" do not exist' % temp
                )
        return temp


    @property
    def parameters(self):
        temp = self._cfg.get('detector:parameters')
        if not os.path.isfile(temp) and self.parameters_old is None:
            raise IOError(
                'detector:parameters specified at "%s" do not exist'
                % temp
                )
        return temp


    @property
    def pixel_size(self):
        temp = self._cfg.get('detector:pixels:size')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def rows(self):
        return self._cfg.get('detector:pixels:rows')
