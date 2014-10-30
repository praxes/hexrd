import glob
import os

from .config import Config



class FileConfig(Config):


    @property
    def stem(self):
        return self._cfg.get('image_series:file:stem')


    @property
    def ids(self):
        temp = self._cfg.get('image_series:file:ids')
        return temp if isinstance(temp, list) else [temp]



class ImagesConfig(Config):


    @property
    def start(self):
        return self._cfg.get('image_series:images:start', default=0)


    @property
    def step(self):
        return self._cfg.get('image_series:images:step', default=1)


    @property
    def stop(self):
        return self._cfg.get('image_series:images:stop', default=None)



class OmegaConfig(Config):


    @property
    def start(self):
        return self._cfg.get('image_series:omega:start')


    @property
    def step(self):
        return self._cfg.get('image_series:omega:step')



class ImageSeriesConfig(Config):


    @property
    def dark(self):
        temp = self._cfg.get(
            'image_series:dark', default=None
            )
        if temp is None or os.path.exists(temp):
            return temp
        raise IOError(
            '"image_series:dark": "%s" does not exist' % temp
            )


    @property
    def file(self):
        return FileConfig(self._cfg)


    @property
    def files(self):
        res = []
        for id in self._cfg.image_series.file.ids:
            id = self._cfg.image_series.file.stem % id
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
    def images(self):
        return ImagesConfig(self._cfg)


    @property
    def omega(self):
        return OmegaConfig(self._cfg)
