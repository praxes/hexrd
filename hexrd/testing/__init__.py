import shutil
import sys
import tempfile
import unittest as ut

import numpy as np


class TestCase(ut.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.mkdtemp(prefix='hexrd-test_')


    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tempdir)


    def mktemp(self, suffix='', prefix='', dir=None):
        if dir is None:
            dir = self.tempdir
        return tempfile.mktemp(suffix, prefix, dir=self.tempdir)


    def assertArrayEqual(self, a1, a2, msg=None, delta=None):
        """
        Make sure a1 and a2 have the same shape and contents to within the
        given precision.
        """
        if delta is None:
            delta = 1e-5
        if msg is None:
            msg = ''
        else:
            msg = ' (%s)' % msg
        a1 = np.asanyarray(a1)
        a2 = np.asanyarray(a2)

        if a1.shape != a2.shape:
            raise self.failureException(
                "Shape mismatch (%s vs %s)%s" % (a1.shape, a2.shape, msg)
                )
        if not np.all(np.abs(a1 - a2) < delta):
            raise self.failureException(
                "Arrays differ by more than %g%s" % (delta, msg)
                )
