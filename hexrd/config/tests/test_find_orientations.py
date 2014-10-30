from .common import TestConfig, test_data


reference_data = \
"""
analysis_name: analysis
image_series:
  omega:
    start: -180
    step: 0.25
---
find_orientations:
  orientation_maps:
    active_hkls: 1
    bin_frames: 2
    file: %(nonexistent_file)s
    threshold: 100
  use_quaternion_grid: %(nonexistent_file)s
  threshold: 5
  extract_measured_g_vectors: true
  omega:
    period: [0, 360]
    tolerance: 1.0
---
find_orientations:
  orientation_maps:
    active_hkls: [1, 2]
    file: %(existing_file)s
  use_quaternion_grid: %(existing_file)s
  seed_search:
    hkl_seeds: 1
    fiber_step: 2.0
  omega:
    tolerance: 3.0
---
find_orientations:
  seed_search:
    hkl_seeds: [1, 2]
""" % test_data



class TestFindOrientationsConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_threshold(self):
        self.assertFalse(
            self.cfgs[0].find_orientations.extract_measured_g_vectors
            )
        self.assertTrue(
            self.cfgs[1].find_orientations.extract_measured_g_vectors
            )
        self.assertFalse(
            self.cfgs[2].find_orientations.extract_measured_g_vectors
            )


    def test_threshold(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.threshold,
            1
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.threshold,
            5
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.threshold,
            1
            )


    def test_use_quaternion_grid(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.use_quaternion_grid,
            None
            )
        self.assertRaises(
            IOError,
            getattr, self.cfgs[1].find_orientations, 'use_quaternion_grid'
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.use_quaternion_grid,
            test_data['existing_file']
            )



class TestOmegaConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_period(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.omega.period,
            [-180, 180]
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.omega.period,
            [0, 360]
            )


    def test_tolerance(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.omega.tolerance,
            0.5
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.omega.tolerance,
            1.0
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.omega.tolerance,
            3.0
            )



class TestSeedSearchConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_hkl_seeds(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].find_orientations.seed_search, 'hkl_seeds'
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.seed_search.hkl_seeds,
            [1]
            )
        self.assertEqual(
            self.cfgs[3].find_orientations.seed_search.hkl_seeds,
            [1, 2]
            )


    def test_fiber_step(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.seed_search.fiber_step,
            0.5
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.seed_search.fiber_step,
            1.0
            )

        self.assertEqual(
            self.cfgs[2].find_orientations.seed_search.fiber_step,
            2.0
            )



class TestOrientationMapsConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_active_hkls(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.orientation_maps.active_hkls,
            'all'
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.orientation_maps.active_hkls,
            [1]
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.orientation_maps.active_hkls,
            [1, 2]
            )


    def test_bin_frames(self):
        self.assertEqual(
            self.cfgs[0].find_orientations.orientation_maps.bin_frames,
            1
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.orientation_maps.bin_frames,
            2
            )


    def test_file(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].find_orientations.orientation_maps, 'file'
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.orientation_maps.file,
            test_data['nonexistent_file']
            )
        self.assertEqual(
            self.cfgs[2].find_orientations.orientation_maps.file,
            test_data['existing_file']
            )


    def test_threshold(self):
        self.assertRaises(
            RuntimeError,
            getattr,
            self.cfgs[0].find_orientations.orientation_maps,
            'threshold'
            )
        self.assertEqual(
            self.cfgs[1].find_orientations.orientation_maps.threshold,
            100
            )
