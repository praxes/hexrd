import os

from .common import TestConfig, test_data


reference_data = \
"""
analysis_name: foo
---
detector:
---
detector:
  parameters: %(nonexistent_file)s
  parameters_old: %(nonexistent_file)s
  pixels:
---
detector:
  parameters: %(existing_file)s
  pixels:
    size: 1
    rows: 1024
    columns: 2048
---
detector:
  parameters_old: %(existing_file)s
  parameters: %(nonexistent_file)s
  pixels:
    size: [1, 2]
""" % test_data


class TestDetectorConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_parameters(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].detector, 'parameters'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].detector, 'parameters'
            )
        self.assertRaises(
            IOError,
            getattr, self.cfgs[2].detector, 'parameters'
            )
        self.assertEqual(
            self.cfgs[3].detector.parameters,
            test_data['existing_file']
            )
        # next test should succeed, converting from old parameters
        self.assertEqual(
            self.cfgs[4].detector.parameters,
            test_data['nonexistent_file']
            )


    def test_parameters_old(self):
        self.assertEqual(self.cfgs[0].detector.parameters_old, None)
        self.assertEqual(self.cfgs[1].detector.parameters_old, None)
        self.assertRaises(
            IOError,
            getattr, self.cfgs[2].detector, 'parameters_old'
            )
        self.assertEqual(
            self.cfgs[4].detector.parameters_old,
            test_data['existing_file']
            )



class TestDetectorPixelsConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_columns(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].detector.pixels, 'columns'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].detector.pixels, 'columns'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].detector.pixels, 'columns'
            )
        self.assertEqual(self.cfgs[3].detector.pixels.columns, 2048)



    def test_size(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].detector.pixels, 'size'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].detector.pixels, 'size'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].detector.pixels, 'size'
            )
        self.assertEqual(self.cfgs[3].detector.pixels.size, [1, 1])
        self.assertEqual(self.cfgs[4].detector.pixels.size, [1, 2])


    def test_rows(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].detector.pixels, 'rows'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].detector.pixels, 'rows'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].detector.pixels, 'rows'
            )
        self.assertEqual(self.cfgs[3].detector.pixels.rows, 1024)
