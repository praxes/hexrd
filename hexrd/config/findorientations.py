import logging
import os

from .config import Config


class FindOrientationsConfig(Config):


    @property
    def orientation_maps(self):
        return OrientationMapsConfig(self._cfg)


class OrientationMapsConfig(Config):


    @property
    def active_hkls(self):
        temp = self._cfg.get(
            'find_orientations:orientation_maps:active_hkls', default='all'
            )
        return [temp] if isinstance(temp, int) else temp


    @property
    def bin_frames(self):
        return self._cfg.get(
            'find_orientations:orientation_maps:bin_frames', default=1
            )


    @property
    def file(self):
        return self._cfg.get('find_orientations:orientation_maps:file')


    @property
    def threshold(self):
        return self._cfg.get('find_orientations:orientation_maps:threshold')
