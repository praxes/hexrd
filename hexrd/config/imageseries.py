import glob
import os

from .config import Config


class ImageSeriesConfig(Config):


    @property
    def dark(self):
        return self._cfg.get(
            'image_series:dark', default=None, path_exists=True
            )


    @property
    def file_stem(self):
        return self._cfg.get('image_series:file:stem')


    @property
    def file_ids(self):
        temp = self._cfg.get('image_series:file:ids')
        return temp if isinstance(temp, list) else [temp]


    @property
    def files(self):
        res = []
        for id in self.file_ids:
            id = self.file_stem % id
            res.extend(glob.glob(id))
        return res


    @property
    def flip(self):
        temp = self._cfg.get('image_series:flip', default=None)
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
        return self._cfg.get('image_series:images:start', default=0)


    @property
    def im_step(self):
        return self._cfg.get('image_series:images:step', default=1)


    @property
    def im_stop(self):
        return self._cfg.get('image_series:images:stop', default=None)


    @property
    def ome_start(self):
        return self._cfg.get('image_series:ome:start')


    @property
    def ome_step(self):
        return self._cfg.get('image_series:ome:step')
