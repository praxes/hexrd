import glob
import os

from .config import Config


class ImageSeriesConfig(Config):


    @property
    def dark(self):
        temp = self._cfg.get('dark')
        if temp and not os.path.isfile(temp):
            raise IOError(
                'image_series:dark file "%s" not found' % temp
                )
        return temp


    @property
    def file_stem(self):
        try:
            return self._get_nested_val('file', 'stem')
        except KeyError:
            raise RuntimeError(
                'image_series:file:stem must be specified in the config file'
                )


    @property
    def file_ids(self):
        try:
            temp = self._get_nested_val('file', 'ids')
            return temp if isinstance(temp, list) else [temp]
        except KeyError:
            raise RuntimeError(
                'image_series:file:ids must be specified in the config file'
                )


    @property
    def files(self):
        res = []
        for id in self.file_ids:
            id = self.file_stem % id
            res.extend(glob.glob(id))
        return res


    @property
    def flip(self):
        temp = self._cfg.get('flip')
        if temp is None:
            return
        temp = temp.lower()
        if temp not in ['h', 'v', 'hv', 'vh', 'cw', 'ccw']:
            raise RuntimeError(
                'image_series:flip setting "%s" is not valid' % temp
                )
        return temp


    @property
    def im_start(self):
        return self._get_nested_val('images', 'start', default=0)


    @property
    def im_step(self):
        return self._get_nested_val('images', 'step', default=1)


    @property
    def im_stop(self):
        return self._get_nested_val('images', 'stop', default=None)


    @property
    def ome_start(self):
        try:
            return self._get_nested_val('ome', 'start')
        except KeyError:
            raise RuntimeError(
                'image_series:ome:start must be specified not found'
                )


    @property
    def ome_step(self):
        try:
            return self._get_nested_val('ome', 'step')
        except KeyError:
            raise RuntimeError(
                'image_series:ome:start must be specified not found'
                )
