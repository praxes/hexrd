#! /usr/bin/env python
#
import sys
import os
import argparse
import unittest
import tempfile

import numpy as np

from hexrd import imageseries
from hexrd.imageseries import save, process, stats, ImageSeries

# ========== Test Data

_NFXY = (3, 7, 5)

class TestImageSeriesProperties(unittest.TestCase):
    def setUp(self):
        self._a = make_array()
        self._is_a = make_array_ims()

    def test_prop_nframes(self):
        self.assertEqual(self._a.shape[0], len(self._is_a))

    def test_prop_shape(self):
        self.assertEqual(self._a.shape[1:], self._is_a.shape)

    def test_prop_dtype(self):
        self.assertEqual(self._a.dtype, self._is_a.dtype)

class TestImageSeriesFmts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        os.rmdir(cls.tmpdir)

    # ==================== Tests

    def test_fmth5(self):
        """save/load HDF5 format"""
        h5file = os.path.join(self.tmpdir, 'test_ims.h5')
        h5path = 'array-data'
        fmt = 'hdf5'

        is_a = make_array_ims()
        save.write(is_a, h5file, fmt, path=h5path)
        is_h = imageseries.open(h5file, fmt, path=h5path)
        diff = compare(is_a, is_h)
        self.assertAlmostEqual(diff, 0., "h5 reconstruction failed")
        self.assertTrue(compare_meta(is_a, is_h))
        del is_h
        os.remove(h5file)

    def test_fmth5_nparray(self):
        """HDF5 format with numpy array metadata"""
        h5file = os.path.join(self.tmpdir, 'test_ims.h5')
        h5path = 'imagedata'
        fmt = 'hdf5'
        key = 'np-array'
        npa = np.array([0,2.0,1.3])

        is_a = make_array_ims()
        is_a.metadata[key] = npa
        save.write(is_a, h5file, fmt, path=h5path)
        is_h = imageseries.open(h5file, fmt, path=h5path)
        meta = is_h.metadata
        diff = np.linalg.norm(meta[key] - npa)
        self.assertAlmostEqual(diff, 0., "h5 numpy array metadata failed")

        del is_h
        os.remove(h5file)

    def test_fmtfc(self):
        """save/load frame-cache format"""
        fcfile = os.path.join(self.tmpdir,  'frame-cache.yml')
        fmt = 'frame-cache'
        thresh = 0.5

        is_a = make_array_ims()
        save.write(is_a, fcfile, fmt,
            threshold=thresh, cache_file='frame-cache.npz')
        is_fc = imageseries.open(fcfile, fmt)
        diff = compare(is_a, is_fc)
        self.assertAlmostEqual(diff, 0., "frame-cache reconstruction failed")
        self.assertTrue(compare_meta(is_a, is_fc))

        del is_fc
        os.remove(fcfile)

    def test_fmtfc_nparray(self):
        """frame-cache format with numpy array metadata"""
        fcfile = os.path.join(self.tmpdir,  'frame-cache.yml')
        fmt = 'frame-cache'
        thresh = 0.5
        key = 'np-array'
        npa = np.array([0,2.0,1.3])

        is_a = make_array_ims()
        is_a.metadata[key] = npa
        save.write(is_a, fcfile, fmt,
            threshold=thresh, cache_file='frame-cache.npz')
        is_fc = imageseries.open(fcfile, fmt)
        meta = is_fc.metadata
        diff = np.linalg.norm(meta[key] - npa)
        self.assertAlmostEqual(diff, 0.,
                               "frame-cache numpy array metadata failed")


        del is_fc
        os.remove(fcfile)

class TestImageSeriesProcess(unittest.TestCase):

    def _runfliptest(self, a, flip, aflip):
        is_a = imageseries.open(None, 'array', data=a)
        ops = [('flip', flip)]
        is_p = process.ProcessedImageSeries(is_a, ops)
        is_aflip = imageseries.open(None, 'array', data=aflip)
        diff = compare(is_aflip, is_p)
        msg = "flipped [%s] image series failed" % flip
        self.assertAlmostEqual(diff, 0., msg=msg)

    def test_process(self):
        """Processed image series"""
        is_a = make_array_ims()
        is_p = process.ProcessedImageSeries(is_a, [])
        diff = compare(is_a, is_p)
        msg = "processed image series failed to reproduce original"
        self.assertAlmostEqual(diff, 0., msg)

    def test_process_flip_t(self):
        """Processed image series: flip transpose"""
        flip = 't'
        a = make_array()
        aflip = np.transpose(a, (0, 2, 1))
        self._runfliptest(a, flip, aflip)

    def test_process_flip_v(self):
        """Processed image series: flip vertical"""
        flip = 'v'
        a = make_array()
        aflip = a[:, :, ::-1]
        self._runfliptest(a, flip, aflip)

    def test_process_flip_h(self):
        """Processed image series: flip horizontal"""
        flip = 'h'
        a = make_array()
        aflip = a[:, ::-1, :]
        self._runfliptest(a, flip, aflip)

    def test_process_flip_vh(self):
        """Processed image series: flip horizontal"""
        flip = 'vh'
        a = make_array()
        aflip = a[:, ::-1, ::-1]
        self._runfliptest(a, flip, aflip)

    def test_process_flip_r90(self):
        """Processed image series: flip horizontal"""
        flip = 'ccw90'
        a = make_array()
        aflip = np.transpose(a, (0, 2, 1))[:, :, ::-1]
        self._runfliptest(a, flip, aflip)

    def test_process_flip_r270(self):
        """Processed image series: flip horizontal"""
        flip = 'cw90'
        a = make_array()
        aflip = np.transpose(a, (0, 2, 1))[:, ::-1, :]
        self._runfliptest(a, flip, aflip)

    def test_process_dark(self):
        """Processed image series: dark image"""
        a = make_array()
        dark = np.ones_like(a[0])
        is_a = imageseries.open(None, 'array', data=a)
        apos = np.where(a >= 1, a-1, 0)
        is_a1 = imageseries.open(None, 'array', data=apos)
        ops = [('dark', dark)]
        is_p = process.ProcessedImageSeries(is_a, ops)
        diff = compare(is_a1, is_p)
        self.assertAlmostEqual(diff, 0., msg="dark image failed")

class TestImageSeriesStats(unittest.TestCase):

    def test_stats_median(self):
        """Processed imageseries: median"""
        a = make_array()
        is_a = imageseries.open(None, 'array', data=a)
        ismed = stats.median(is_a)
        amed = np.median(a, axis=0)
        err = np.linalg.norm(amed - ismed)
        self.assertAlmostEqual(err, 0., msg="median image failed")

    def test_stats_max(self):
        """Processed imageseries: median"""
        a = make_array()
        is_a = imageseries.open(None, 'array', data=a)
        ismax = stats.max(is_a)
        amax = np.max(a, axis=0)
        err = np.linalg.norm(amax - ismax)
        self.assertAlmostEqual(err, 0., msg="max image failed")

# ==================== Utility functions

def make_array():
    a = np.zeros(_NFXY)
    ind = np.array([0,1,2])
    a[ind, 1,2] = 1 + ind
    return a

def make_meta():
    return {'testing': '1,2,3'}

def make_array_ims():
    is_a = imageseries.open(None, 'array', data=make_array(),
                            meta=make_meta())
    return is_a

def compare(ims1, ims2):
    """compare two imageseries"""
    if len(ims1) != len(ims2):
        raise ValueError("lengths do not match")

    if ims1.dtype is not ims2.dtype:
        raise ValueError("types do not match")

    maxdiff = 0.0
    for i in range(len(ims1)):
        f1 = ims1[i]
        f2 = ims2[i]
        fdiff = np.linalg.norm(f1 - f2)
        maxdiff = np.maximum(maxdiff, fdiff)

    return maxdiff

def compare_meta(ims1, ims2):
    # check metadata (simple immutable cases only for now)

    m1 = set(ims1.metadata.items())
    m2 = set(ims2.metadata.items())
    return m1.issubset(m2) and m2.issubset(m1)

# ================================================== Execution
#
def set_options():
    """Set options for command line"""
    parser = argparse.ArgumentParser(description='test hexrd.imageseries')

    return parser

def execute(args):
    """Main execution"""
    p = set_options()
    unittest.main()

    return


if __name__ == '__main__':
    execute(sys.argv[1:])
