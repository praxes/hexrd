import os
import tempfile

from .common import TestConfig, test_data


reference_data = \
"""
image_series:
#  dark: # not specified to test default is None
---
image_series:
  dark: %(existing_file)s
  file:
  flip: V
  images:
---
image_series:
  dark: %(nonexistent_file)s
  file:
    stem: %(file_stem)s
    ids: [1]
  flip: triple_lindy
  images:
    start: 1
    step: 2
    stop: -1
  ome:
    start: 0
    step: 0.25
---
image_series:
  file:
    ids: 2
---
image_series:
  file:
    ids: [1,2]
---
image_series:
  file:
    stem: %(tempdir)s%(pathsep)s%%s.dat
    ids: ["*001*"]
""" % test_data


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


    def test_file_ids(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'file_ids'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].image_series, 'file_ids'
            )
        self.assertEqual(
            self.cfgs[2].image_series.file_ids,
            [1]
            )
        self.assertEqual(
            self.cfgs[3].image_series.file_ids,
            [2]
            )
        self.assertEqual(
            self.cfgs[4].image_series.file_ids,
            [1, 2]
            )


    def test_files(self):
        files = []
        for i in ['00011.dat', '00012.dat', '00021.dat']:
            with tempfile.NamedTemporaryFile(delete=False, suffix=i) as f:
                files.append(f.name)
                f.file.write('foo')
        try:
            self.assertEqual(
                sorted(self.cfgs[5].image_series.files),
                sorted(files[:2])
                )
        finally:
            for f in files:
                os.remove(f)


    def test_flip(self):
        self.assertEqual(self.cfgs[0].image_series.flip, None)
        self.assertEqual(self.cfgs[1].image_series.flip, 'v')
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[2].image_series, 'flip'
            )


    def test_im_start(self):
        self.assertEqual(self.cfgs[0].image_series.im_start, 0)
        self.assertEqual(self.cfgs[1].image_series.im_start, 0)
        self.assertEqual(self.cfgs[2].image_series.im_start, 1)


    def test_im_step(self):
        self.assertEqual(self.cfgs[0].image_series.im_step, 1)
        self.assertEqual(self.cfgs[1].image_series.im_step, 1)
        self.assertEqual(self.cfgs[2].image_series.im_step, 2)


    def test_im_stop(self):
        self.assertEqual(self.cfgs[0].image_series.im_stop, None)
        self.assertEqual(self.cfgs[1].image_series.im_stop, None)
        self.assertEqual(self.cfgs[2].image_series.im_stop, -1)


    def test_ome_start(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'ome_start'
            )
        self.assertEqual(self.cfgs[2].image_series.ome_start, 0)



    def test_ome_step(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series, 'ome_start'
            )
        self.assertEqual(self.cfgs[2].image_series.ome_step, 0.25)
