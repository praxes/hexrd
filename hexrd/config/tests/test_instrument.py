import os

from .common import TestConfig, test_data
from ..instrument import Instrument, Beam

reference_data = \
"""
beam: {}
---
beam:
  energy: 2.0
  vector: {azimuth: 0.0, polar_angle: 0.0}
---
instrument: instrument.yaml
""" % test_data


class TestInstrument(TestConfig):


    @classmethod
    def get_reference_data(cls):
        return reference_data


class TestBeam(TestConfig):

    @classmethod
    def get_reference_data(cls):
        return reference_data

    def test_beam_energy_dflt(self):
        bcfg = Beam(self.cfgs[0])
        energy = bcfg.energy
        self.assertEqual(energy, Beam.beam_energy_DFLT, "Incorrect default beam energy")

    def test_beam_energy(self):
        bcfg = Beam(self.cfgs[1])
        energy = bcfg.energy
        self.assertEqual(energy, 2.0, "Incorrect beam energy")

    def test_beam_vector_dflt(self):
        bcfg = Beam(self.cfgs[0])
        bvecdflt = Beam.beam_vec_DFLT
        bvec = bcfg.vector

        self.assertEqual(bvec[0], bvecdflt[0], "Incorrect default beam vector")
        self.assertEqual(bvec[1], bvecdflt[1], "Incorrect default beam vector")
        self.assertEqual(bvec[2], bvecdflt[2], "Incorrect default beam vector")

    def test_beam_vector_dflt(self):
        bcfg = Beam(self.cfgs[0])
        bvecdflt = Beam.beam_vec_DFLT
        bvec = bcfg.vector

        self.assertEqual(bvec[0], bvecdflt[0], "Incorrect default beam vector")
        self.assertEqual(bvec[1], bvecdflt[1], "Incorrect default beam vector")
        self.assertEqual(bvec[2], bvecdflt[2], "Incorrect default beam vector")

    def test_beam_vector(self):
        bcfg = Beam(self.cfgs[1])
        bvec = bcfg.vector

        self.assertEqual(bvec[0], 0.0, "Incorrect default beam vector")
        self.assertEqual(bvec[1], -1.0, "Incorrect default beam vector")
        self.assertEqual(bvec[2], 0.0, "Incorrect default beam vector")
