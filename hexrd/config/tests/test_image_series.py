import os
import tempfile

from .common import TestConfig, test_data


reference_data = \
"""
working_dir: %(tempdir)s
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
  omega:
    start: 0
    step: 0.25
    stop: 360
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
---
image_series:
  file:
    stem: %(tempdir)s%(pathsep)s%%05d.dat
    ids: [0, 1]
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
            self.assertRaises(
                IOError,
                getattr, self.cfgs[6].image_series, 'files'
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



class TestFileConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_stem(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series.file, 'stem'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].image_series.file, 'stem'
            )
        self.assertEqual(
            self.cfgs[2].image_series.file.stem,
            os.path.join(test_data['tempdir'], test_data['file_stem'])
            )


    def test_ids(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series.file, 'ids'
            )
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[1].image_series.file, 'ids'
            )
        self.assertEqual(
            self.cfgs[2].image_series.file.ids,
            [1]
            )
        self.assertEqual(
            self.cfgs[3].image_series.file.ids,
            [2]
            )
        self.assertEqual(
            self.cfgs[4].image_series.file.ids,
            [1, 2]
            )



class TestImagesConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_start(self):
        self.assertEqual(self.cfgs[0].image_series.images.start, 0)
        self.assertEqual(self.cfgs[1].image_series.images.start, 0)
        self.assertEqual(self.cfgs[2].image_series.images.start, 1)


    def test_step(self):
        self.assertEqual(self.cfgs[0].image_series.images.step, 1)
        self.assertEqual(self.cfgs[1].image_series.images.step, 1)
        self.assertEqual(self.cfgs[2].image_series.images.step, 2)


    def test_stop(self):
        self.assertEqual(self.cfgs[0].image_series.images.stop, None)
        self.assertEqual(self.cfgs[1].image_series.images.stop, None)
        self.assertEqual(self.cfgs[2].image_series.images.stop, -1)



class TestOmegaConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_start(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series.omega, 'start'
            )
        self.assertEqual(self.cfgs[2].image_series.omega.start, 0)



    def test_step(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series.omega, 'step'
            )
        self.assertEqual(self.cfgs[2].image_series.omega.step, 0.25)



    def test_stop(self):
        self.assertRaises(
            RuntimeError,
            getattr, self.cfgs[0].image_series.omega, 'stop'
            )
        self.assertEqual(self.cfgs[2].image_series.omega.stop, 360)
