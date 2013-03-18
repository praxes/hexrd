#! /usr/bin/env python
#
import sys
import argparse
import unittest

import numpy as np

from hexrd.quadrature import q1db

class TestQuadrature(unittest.TestCase):
    """Verify quadrature rules"""
    def setUp(self):
        self.tol = 1.0e-14
        self.monomial = lambda d: np.poly1d([1] + d*[0])
        self.integral = lambda x, w, f: \
          np.sum([w[i]*f(x[i]) for i in range(len(x))])

    def test_1d(self):
        """1D Quadrature Rules"""
        for npts in [1,2,3,4,5,8]:
            self._check(npts, 2*npts -1)

    def _check(self, npts, dmax):
        """Check for all monomials up to degree dmax"""
        print 'Checking %d-point rule:' % npts
        xi, wt = q1db.qLoc(npts)
        self.assertEqual(len(xi), len(wt))

        for d in range(dmax + 1):
            integral = self.integral(xi, wt, self.monomial(d))
            trueval = (float(1)/float(d + 1))
            print '    degree: %d, value: %f, true = %f' % (d, integral, trueval)
            self.assertTrue(np.abs(integral - trueval) < self.tol,
                            'degree: %d, value: %f, true = %f' % (d, integral, trueval))

        print '... %d-point rule:  OK\n' % npts

        return
#
# ================================================== Execution
#
def set_options():
    """Set options for command line"""
    parser = argparse.ArgumentParser(description='test MDE')

    return parser

def execute(args):
    """Main execution"""
    p = set_options()
    unittest.main()

    return

if __name__ == '__main__':
    execute(sys.argv[1:])
