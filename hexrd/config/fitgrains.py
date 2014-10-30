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
    def ome(self):
        temp = self._cfg.get('fit_grains:tolerance:ome')
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
    def parameters_old(self):
        pass


    @property
    def parameters(self):
        pass


    @property
    def tolerance(self):
        return ToleranceConfig(self._cfg)
