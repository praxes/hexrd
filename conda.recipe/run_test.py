import unittest

suite = unittest.TestLoader().discover('hexrd')
unittest.TextTestRunner().run(suite)
