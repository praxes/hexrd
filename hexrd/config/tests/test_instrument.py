import os

import hexrd.instrument
from .common import TestConfig, test_data
from ..instrument import Instrument, Beam, OscillationStage

reference_data = \
"""
beam: {}
---
beam:
  energy: 2.0
  vector: {azimuth: 0.0, polar_angle: 0.0}
---
oscillation_stage:
  chi: 0.05
  t_vec_s: [1., 2., 3.]
---
instrument: instrument.yaml
""" % test_data


class TestInstrument(TestConfig):

    @classmethod
    def get_reference_data(cls):
        return reference_data

    def test_beam(self):
        icfg = Instrument(self.cfgs[1])
        b = icfg.beam
        self.assertTrue(isinstance(b, hexrd.instrument.beam.Beam), "Failed to produce a Beam instance")

    def test_oscillation_stage(self):
        icfg = Instrument(self.cfgs[2])
        ostage = icfg.oscillation_stage
        self.assertTrue(isinstance(ostage, hexrd.instrument.oscillation_stage.OscillationStage),
                        "Failed to produce an OscillationStage instance")


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

    def test_beam_vector(self):
        bcfg = Beam(self.cfgs[1])
        bvec = bcfg.vector

        self.assertEqual(bvec[0], 0.0, "Incorrect default beam vector")
        self.assertEqual(bvec[1], -1.0, "Incorrect default beam vector")
        self.assertEqual(bvec[2], 0.0, "Incorrect default beam vector")


class TestOscillationStage(TestConfig):

    @classmethod
    def get_reference_data(cls):
        return reference_data

    def test_chi_dflt(self):
        oscfg = OscillationStage(self.cfgs[0])
        self.assertEqual(oscfg.chi, OscillationStage.chi_DFLT, "Incorrect default chi for oscillation stage")

    def test_chi(self):
        oscfg = OscillationStage(self.cfgs[2])
        self.assertEqual(oscfg.chi, 0.05, "Incorrect default chi for oscillation stage")

    def test_tvec_dflt(self):
        oscfg = OscillationStage(self.cfgs[0])
        tvec_dflt = OscillationStage.tvec_DFLT
        tvec = oscfg.tvec

        self.assertEqual(tvec[0], tvec_dflt[0], "Incorrect default translation vector")
        self.assertEqual(tvec[1], tvec_dflt[1], "Incorrect default translation vector")
        self.assertEqual(tvec[2], tvec_dflt[2], "Incorrect default translation vector")

    def test_tvec(self):
        oscfg = OscillationStage(self.cfgs[2])
        tvec = oscfg.tvec

        self.assertEqual(tvec[0], 1., "Incorrect translation vector")
        self.assertEqual(tvec[1], 2., "Incorrect translation vector")
        self.assertEqual(tvec[2], 3., "Incorrect translation vector")
