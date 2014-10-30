import os

from .config import Config



class PixelsConfig(Config):


    @property
    def columns(self):
        return self._cfg.get('detector:pixels:columns')


    @property
    def size(self):
        temp = self._cfg.get('detector:pixels:size')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def rows(self):
        return self._cfg.get('detector:pixels:rows')



class DetectorConfig(Config):


    @property
    def parameters_old(self):
        temp = self._cfg.get('detector:parameters_old', default=None)
        if temp is None:
            return temp
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        if os.path.exists(temp):
            return temp
        raise IOError(
            '"detector:parameters_old": "%s" does not exist' % temp
            )


    @property
    def parameters(self):
        temp = self._cfg.get('detector:parameters')
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        if os.path.exists(temp):
            return temp
        if self.parameters_old is not None:
            # converting old parameter file
            return temp
        raise IOError(
            '"detector:parameters": "%s" does not exist' % temp
            )



    @property
    def pixels(self):
        return PixelsConfig(self._cfg)
