from .common import TestConfig, test_data


reference_data = \
"""
analysis_name: analysis
---
find_orientations:
  orientation_maps:
    active_hkls: 1
    bin_frames: 2
    file: %(nonexistent_file)s
    threshold: 100
---
find_orientations:
  orientation_maps:
    active_hkls: [1, 2]
    file: %(existing_file)s
""" % test_data


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
