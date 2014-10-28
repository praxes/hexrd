import os

from .config import Config


class DetectorConfig(Config):


    @property
    def columns(self):
        try:
            return self._get_nested_val('pixels', 'columns')
        except KeyError:
            raise RuntimeError(
                'detector:pixels:colums must be specified in the config file'
                )


    @property
    def parameters_old(self):
        temp = self._cfg.get('parameters_old', None)
        if temp is None:
            return
        if not os.path.isfile(temp):
            raise IOError(
                'detector:parameters_old specified at "%s" do not exist' % temp
                )
        return temp


    @property
    def parameters(self):
        try:
            temp = self._cfg['parameters']
            if not os.path.isfile(temp) and self.parameters_old is None:
                raise IOError(
                    'detector:parameters specified at "%s" do not exist'
                    % temp
                    )
            return temp
        except KeyError:
            raise RuntimeError(
                'detector:parameters must be specified in the config file.'
                )


    @property
    def pixel_size(self):
        try:
            temp = self._get_nested_val('pixels', 'size')
            if isinstance(temp, (int, float)):
                temp = [temp, temp]
            return temp
        except KeyError:
            raise RuntimeError(
                'detector:pixels:size must be specified in the config file'
                )


    @property
    def rows(self):
        try:
            return self._get_nested_val('pixels', 'rows')
        except KeyError:
            raise RuntimeError(
                'detector:pixels:rows must be specified in the config file'
                )
