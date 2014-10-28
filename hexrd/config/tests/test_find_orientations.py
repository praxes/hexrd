from .common import TestConfig, test_data


reference_data = \
"""
""" % test_data


class TestMaterialConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


    def test_definitions(self):
        self.assertEqual(
            )
        self.assertRaises(
            )


    def test_active(self):
        self.assertRaises(
            )
