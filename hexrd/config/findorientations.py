import logging
import os

from .config import Config


class MaterialConfig(Config):


    @property
    def active_hkls(self):
        temp = self._get_nested_val(
            'orientation_maps', 'active_hkls', default=None
            )
        if temp is None:
            temp = 'all'
            logging.info(
                'find_orientations:orientation_maps:active_hkls not specified.'
                ' defaulting to all hkls defined for active material'
                )
        return temp


    @property
    def bin_frames(self):
        temp = self._get_nested_val(
            'orientation_maps', 'bin_frames', default=None
            )
        if temp is None:
            temp = 1
            logging.info(
                'find_orientations:orientation_maps:bin_frames not specified.'
                ' defaulting to 1'
                )
        return temp


    @property
    def orientation_maps_file(self):
        try:
            return self._get_nested_val('orientation_maps', 'file')
        except KeyError:
            raise RuntimeError(
                'find_orientations:orientation_maps:file must be specified in'
                ' the config file'
                )


    @property
    def threshold(self):
        try:
            return self._get_nested_val('orientation_maps', 'threshold')
        except KeyError:
            raise RuntimeError(
                'find_orientations:orientation_maps:threshold not specified'
                )
