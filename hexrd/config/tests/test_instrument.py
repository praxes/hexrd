import os

from .common import TestConfig, test_data


reference_data = \
"""
analysis_name: foo
working_dir: %(tempdir)s
---
instrument:
---
instrument:
  parameters: %(nonexistent_file)s
  detector:
    parameters_old: %(nonexistent_file)s
    pixels:
---
instrument:
  parameters: %(existing_file)s
  detector:
    pixels:
      size: 1
      rows: 1024
      columns: 2048
---
instrument:
  parameters: %(nonexistent_file)s
  detector:
    parameters_old: %(existing_file)s
    pixels:
      size: [1, 2]
""" % test_data


class TestInstrumentConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_parameters(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].instrument, 'parameters'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].instrument, 'parameters'
            )
        self.assertRaises(
            IOError,
            getattr, self.cfgs[2].instrument, 'parameters'
            )
        self.assertEqual(
            self.cfgs[3].instrument.parameters,
            test_data['existing_file']
            )
        # next test should succeed, converting from old parameters
        self.assertEqual(
            self.cfgs[4].instrument.parameters,
            os.path.join(test_data['tempdir'], test_data['nonexistent_file'])
            )



class TestDetectorConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_parameters_old(self):
        self.assertEqual(self.cfgs[0].instrument.detector.parameters_old, None)
        self.assertEqual(self.cfgs[1].instrument.detector.parameters_old, None)
        self.assertRaises(
            IOError,
            getattr, self.cfgs[2].instrument.detector, 'parameters_old'
            )
        self.assertEqual(
            self.cfgs[4].instrument.detector.parameters_old,
            os.path.join(test_data['tempdir'], test_data['existing_file'])
            )



class TestDetectorPixelsConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_columns(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].instrument.detector.pixels, 'columns'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].instrument.detector.pixels, 'columns'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].instrument.detector.pixels, 'columns'
            )
        self.assertEqual(self.cfgs[3].instrument.detector.pixels.columns, 2048)



    def test_size(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].instrument.detector.pixels, 'size'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].instrument.detector.pixels, 'size'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].instrument.detector.pixels, 'size'
            )
        self.assertEqual(self.cfgs[3].instrument.detector.pixels.size, [1, 1])
        self.assertEqual(self.cfgs[4].instrument.detector.pixels.size, [1, 2])


    def test_rows(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].instrument.detector.pixels, 'rows'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].instrument.detector.pixels, 'rows'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].instrument.detector.pixels, 'rows'
            )
        self.assertEqual(self.cfgs[3].instrument.detector.pixels.rows, 1024)
