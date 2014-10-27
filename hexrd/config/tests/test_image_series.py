import os

from .common import TestConfig, test_data


reference_data = \
"""
image_series:
#  path: %(existing_path)s
#  dark: # not specified to test default is None
---
image_series:
  path: %(nonexistent_path)s
  dark: %(existing_file)s
  file:
---
image_series:
  path: %(existing_path)s
  dark: %(nonexistent_file)s
  file:
    stem: %(file_stem)s
""" % test_data

print test_data


class TestImageSeriesConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_dark(self):
        self.assertEqual(self.cfgs[0].image_series.dark, None)
        self.assertEqual(
            self.cfgs[1].image_series.dark,
            test_data['existing_file']
            )
        self.assertRaises(
            IOError,
            getattr, self.cfgs[2].image_series, 'dark'
            )


    def test_file_stem(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'file_stem'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].image_series, 'file_stem'
            )
        self.assertEqual(
            self.cfgs[2].image_series.file_stem,
            test_data['file_stem']
            )


    def test_path(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'path'
            )
        self.assertRaises(
            IOError,
            getattr, self.cfgs[1].image_series, 'path'
            )
        self.assertEqual(
            self.cfgs[2].image_series.path,
            test_data['existing_path']
            )
