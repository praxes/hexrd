import os

from .config import Config



class ToleranceConfig(Config):


    @property
    def eta(self):
        temp = self._cfg.get('fit_grains:tolerance:eta')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def omega(self):
        temp = self._cfg.get('fit_grains:tolerance:omega')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def tth(self):
        temp = self._cfg.get('fit_grains:tolerance:tth')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp



class FitGrainsConfig(Config):


    @property
    def do_fit(self):
        return self._cfg.get('fit_grains:do_fit', True)


    @property
    def npdiv(self):
        return self._cfg.get('fit_grains:npdiv', 2)


    @property
    def panel_buffer(self):
        temp = self._cfg.get('fit_grains:panel_buffer')
        if isinstance(temp, (int, float)):
            temp = [temp, temp]
        return temp


    @property
    def threshold(self):
        return self._cfg.get('fit_grains:threshold')


    @property
    def tolerance(self):
        return ToleranceConfig(self._cfg)


    @property
    def tth_max(self):
        key = 'fit_grains:tth_max'
        temp = self._cfg.get(key, True)
        if temp in (True, False):
            return temp
        if isinstance(temp, (int, float)):
            if temp > 0:
                return temp
        raise RuntimeError(
            '"%s" must be > 0, true, or false, got "%s"' % (key, temp)
            )
