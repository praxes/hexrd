import glob
import os

from .config import Config



class FileConfig(Config):


    @property
    def stem(self):
        temp = self._cfg.get('image_series:file:stem')
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        return temp


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


    @property
    def stop(self):
        return self._cfg.get('image_series:omega:stop', default=None)




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
        stem = self._cfg.image_series.file.stem
        res = []
        missing = []
        for id in self._cfg.image_series.file.ids:
            try:
                id = stem % id
            except TypeError:
                # string interpolation failed, join stem and id:
                id = stem + id
            temp = glob.glob(id)
            if temp:
                res.extend(temp)
            else:
                missing.append(id)
        if missing:
            raise IOError(
                'Image files not found: %s' % (', '.join(missing))
                )
        return sorted(res)


    @property
    def flip(self):
        temp = self._cfg.get('image_series:flip', default=None)
        if temp is None:
            return
        temp = temp.lower()
        if temp not in ['h', 'v', 'hv', 'vh', 'cw90', 'ccw90', 't']:
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


    @property
    def n_frames(self):
        return (self.omega.stop - self.omega.start)/self.omega.step
