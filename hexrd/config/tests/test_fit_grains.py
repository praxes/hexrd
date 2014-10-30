import os

from .common import TestConfig, test_data


reference_data = \
"""
""" % test_data


class TestFitGrainsConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data



class TestToleranceConfig(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data
