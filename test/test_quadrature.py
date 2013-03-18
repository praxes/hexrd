#! /usr/bin/env python
#
import sys
import argparse
import unittest

import numpy as np

from hexrd.quadrature import q1db
from hexrd.quadrature import q2db
from hexrd.quadrature import q3db

class TestQuadrature(unittest.TestCase):
    """Verify quadrature rules"""
    def setUp(self):
        self.tol = 1.0e-14
        self.monomial = lambda d: np.poly1d([1] + d*[0])
        self.integral = lambda x, w, f: \
          np.sum([w[i]*f(x[i]) for i in range(len(x))])

    def test_1d(self):
        """1D Quadrature Rules"""
        for npts in q1db.available:
            self.__check1d(npts, 2*npts -1)

    def test_2d(self):
        """2D Quadrature Rules"""
        print "\nChecking 2D rules"
        rules1d = q1db.available
        for n1 in rules1d:
            for n2 in rules1d:
                self.__checknd([n1, n2])

    def test_3d(self):
        """3D Quadrature Rules"""
        print "\nChecking 3D rules"
        rules1d = q1db.available
        for n1 in rules1d:
            for n2 in rules1d:
                for n3 in rules1d:
                    self.__checknd([n1, n2, n3])

    def __check1d(self, npts, dmax):
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

    def __checknd(self, npts):
        """Check product rules"""
        dim = len(npts)
        btmpl = '%d' + (dim - 1) * 'b%d'
        xtmpl = '%dx%d'
        if dim == 2:
            n0 = npts[0]
            n1 = npts[1]
            if n0 == n1:
                rule = xtmpl % (dim,n0) if (n0 > 3) else n0
            else:
                rule = btmpl % (n0, n1)
        elif dim == 3:
            n0 = npts[0]
            n1 = npts[1]
            n2 = npts[2]
            if n0 == n1 and n0 == n2:
                rule = xtmpl % (dim,n0) if (n0 > 3) else n0
            else:
                rule = btmpl % (n0, n1, n2)

        qclass = q2db if dim == 2 else q3db
        xi, wt = qclass.qLoc(rule)

        fi = [self.monomial(2*n - 1) for n in npts]
        f = lambda(x): np.prod([fi[i](x[i]) for i in range(dim)])
        integral = self.integral(xi, wt, f)
        tmp = [1.0/(2*float(n)) for n in npts]
        trueval = np.prod(tmp)

        self.assertTrue(np.abs(integral - trueval) < self.tol,
                        '%s  value: %f, true = %f' % (repr(rule), integral, trueval))
        print '   nD rule: %s (%s) OK' % (repr(rule), repr(npts))

        return
#
# ================================================== Execution
#
def set_options():
    """Set options for command line"""
    parser = argparse.ArgumentParser(description='test hexrd.quadrature')

    return parser

def execute(args):
    """Main execution"""
    p = set_options()
    unittest.main()

    return

if __name__ == '__main__':
    execute(sys.argv[1:])
