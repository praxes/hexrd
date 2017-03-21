import os
import tempfile

from .common import TestConfig, test_data


reference_data = \
"""
analysis_name: analysis
working_dir: %(tempdir)s
---
image_series:
  filename: %(nonexistent_file)s
  format: hdf5
  args:
     path: %(nonexistent_path)s
---
image_series:
  filename: %(nonexistent_file)s
  format: frame-cache
  args:
""" % test_data


class TestImageSeriesConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_filename(self):

        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'filename'
            )

        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'format'
            )

        self.assertEqual(
            self.cfgs[1].image_series.filename,
            os.path.join(test_data['tempdir'], test_data['nonexistent_file'])
            )

        self.assertEqual(
            self.cfgs[1].image_series.format, 'hdf5'
            )

        a = self.cfgs[1].image_series.args
        self.assertEqual(
            a['path'], test_data['nonexistent_path']
            )

        a = self.cfgs[2].image_series.args
        self.assertEqual(a, None)
