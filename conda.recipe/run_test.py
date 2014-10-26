import unittest

suite = unittest.TestLoader().discover('hexrd')
unittest.TextTestRunner(verbosity=self.verbosity+1).run(suite)
