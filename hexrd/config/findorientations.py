import os

import numpy as np

from .config import Config


class FindOrientationsConfig(Config):

    @property
    def clustering(self):
        return ClusteringConfig(self._cfg)

    @property
    def eta(self):
        return EtaConfig(self._cfg)

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
        if temp is None:
            return temp
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        if os.path.isfile(temp):
            return temp
        raise IOError(
            '"%s": "%s" does not exist' % (key, temp)
            )


class ClusteringConfig(Config):

    @property
    def algorithm(self):
        key = 'find_orientations:clustering:algorithm'
        choices = ['dbscan', 'ort-dbscan', 'sph-dbscan', 'fclusterdata']
        temp = self._cfg.get(key, 'dbscan').lower()
        if temp in choices:
            return temp
        raise RuntimeError(
            '"%s": "%s" not recognized, must be one of %s'
            % (key, temp, choices)
            )

    @property
    def completeness(self):
        key = 'find_orientations:clustering:completeness'
        temp = self._cfg.get(key, None)
        if temp is not None:
            return temp
        raise RuntimeError(
            '"%s" must be specified' % key
            )

    @property
    def radius(self):
        key = 'find_orientations:clustering:radius'
        temp = self._cfg.get(key, None)
        if temp is not None:
            return temp
        raise RuntimeError(
            '"%s" must be specified' % key
            )


class OmegaConfig(Config):

    @property
    def period(self):
        key = 'find_orientations:omega:period'
        ome_start = self._cfg.image_series.omega.start
        range = 360 if self._cfg.image_series.omega.step > 0 else -360
        try:
            temp = self._cfg.get(key, [ome_start, ome_start + range])
        except(TypeError):
            temp = self._cfg.get(key, None)

        try:
            range = np.abs(temp[1] - temp[0])
        except(TypeError):
            raise RuntimeError(
                "without imageseries spec, must specify period"
            )
        if range != 360:
            raise RuntimeError(
                '"%s": range must be 360 degrees, range of %s is %g'
                % (key, temp, range)
                )
        return temp

    @property
    def tolerance(self):
        return self._cfg.get(
            'find_orientations:omega:tolerance',
            2 * self._cfg.image_series.omega.step
            )


class EtaConfig(Config):

    @property
    def tolerance(self):
        return self._cfg.get(
            'find_orientations:eta:tolerance',
            2 * self._cfg.image_series.omega.step
            )

    @property
    def mask(self):
        return self._cfg.get('find_orientations:eta:mask', 5)

    @property
    def range(self):
        mask = self.mask
        if mask is None:
            return mask
        return [[-90 + mask, 90 - mask],
                [90 + mask, 270 - mask]]


class SeedSearchConfig(Config):

    @property
    def hkl_seeds(self):
        key = 'find_orientations:seed_search:hkl_seeds'
        try:
            temp = self._cfg.get(key)
            if isinstance(temp, int):
                temp = [temp, ]
            return temp
        except:
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

    @property
    def fiber_ndiv(self):
        return int(360.0 / self.fiber_step)


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
        temp = self._cfg.get('find_orientations:orientation_maps:file')
        if not os.path.isabs(temp):
            temp = os.path.join(self._cfg.working_dir, temp)
        return temp

    @property
    def threshold(self):
        return self._cfg.get('find_orientations:orientation_maps:threshold')
