import os

from .config import Config



class PixelsConfig(Config):


    @property
    def columns(self):
        return self._cfg.get('instrument:detector:pixels:columns')


    @property
    def size(self):
        temp = self._cfg.get('instrument:detector:pixels:size')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def rows(self):
        return self._cfg.get('instrument:detector:pixels:rows')



class DetectorConfig(Config):


    @property
    def parameters_old(self):
        key = 'instrument:detector:parameters_old'
        temp = self._cfg.get(key, default=None)
        if temp is None:
            return temp
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        if os.path.exists(temp):
            return temp
        raise IOError(
            '"%s": "%s" does not exist' % (key, temp)
            )


    @property
    def pixels(self):
        return PixelsConfig(self._cfg)



class InstrumentConfig(Config):


    @property
    def detector(self):
        return DetectorConfig(self._cfg)


    @property
    def parameters(self):
        key = 'instrument:parameters'
        temp = self._cfg.get(key)
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        if os.path.exists(temp):
            return temp
        if self.detector.parameters_old is not None:
            # converting old parameter file
            return temp
        raise IOError(
            '"%s": "%s" does not exist' % (key, temp)
            )
