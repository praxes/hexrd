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
        pass


    @property
    def flip(self):
        pass


    @property
    def im_start(self):
        pass


    @property
    def im_step(self):
        pass


    @property
    def im_stop(self):
        pass


    @property
    def ome_start(self):
        pass


    @property
    def ome_step(self):
        pass


    @property
    def path(self):
        try:
            temp = self._cfg['path']
            if not os.path.isdir(temp):
                raise IOError(
                    'image_series:path specified as "%s" does not exist' % temp
                    )
            return temp
        except KeyError:
            raise RuntimeError(
                'image_series:path must be specified in the config file'
                )
