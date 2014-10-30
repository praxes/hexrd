import logging
import os

from .config import Config



class FindOrientationsConfig(Config):


    @property
    def extract_measured_g_vectors(self):
        return self._cfg.get(
            'find_orientations:extract_measured_g_vectors',
            False
            )


    @property
    def omega(self):
        return OmegaConfig(self._cfg)


    @property
    def orientation_maps(self):
        return OrientationMapsConfig(self._cfg)


    @property
    def seed_search(self):
        return SeedSearchConfig(self._cfg)


    @property
    def threshold(self):
        return self._cfg.get('find_orientations:threshold', 1)


    @property
    def use_quaternion_grid(self):
        key = 'find_orientations:use_quaternion_grid'
        temp = self._cfg.get(key, None)
        if temp is None or os.path.isfile(temp):
            return temp
        raise IOError(
            '"%s": "%s" does not exist' % (key, temp)
            )


class OmegaConfig(Config):


    @property
    def period(self):
        temp = self._cfg.get('find_orientations:omega:period', None)
        if temp is None:
            temp = self._cfg.image_series.omega.start
            temp = [temp, temp+360]
        return temp


    @property
    def tolerance(self):
        return self._cfg.get(
            'find_orientations:omega:tolerance',
            2 * self._cfg.image_series.omega.step
            )



class SeedSearchConfig(Config):


    @property
    def hkl_seeds(self):
        key = 'find_orientations:seed_search:hkl_seeds'
        try:
            temp = self._cfg.get(key)
            if isinstance(temp, int):
                temp = [temp]
            return temp
        except KeyError:
            if self._cfg.find_orientations.use_quaternion_grid is None:
                raise RuntimeError(
                    '"%s" must be defined for seeded search' % key
                    )


    @property
    def fiber_step(self):
        return self._cfg.get(
            'find_orientations:seed_search:fiber_step',
            self._cfg.find_orientations.omega.tolerance
            )



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
