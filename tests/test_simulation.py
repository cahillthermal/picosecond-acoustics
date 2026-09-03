"""
Unit tests for picosecond_acoustics package using standard unittest.
"""

import unittest
import numpy as np
from pathlib import Path

from picosecond_acoustics import (
    Material,
    FilmLayer,
    Structure,
    LightBeam,
    SimulationConfig,
    PicosecondAcousticsSimulation,
    load_simulation_from_legacy_inp
)


class TestPicosecondAcoustics(unittest.TestCase):

    def test_material_load(self):
        al = Material.from_inp_file("Al.inp")
        self.assertEqual(al.density, 2.7)
        self.assertEqual(al.sound_velocity, 64.2)
        self.assertEqual(al.n_pump, 2.8)
        self.assertEqual(al.kappa_pump, 8.4)

    def test_simulation_run(self):
        al = Material.from_inp_file("Al.inp")
        si = Material.from_inp_file("Si.inp")

        structure = Structure()
        structure.add_layer(FilmLayer(material=al, thickness_mean=760.0, roughness=0.0))
        structure.add_layer(FilmLayer(material=si, thickness_mean=7000.0, roughness=0.0))

        pump = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization=1)
        probe = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization=2)
        config = SimulationConfig(tau_data=0.2, n_tau_divide=1, t_stop=100.0, t_save=50.0, n_samples=1)

        sim = PicosecondAcousticsSimulation(
            structure=structure,
            pump_beam=pump,
            probe_beam=probe,
            sim_config=config
        )

        res = sim.run()

        self.assertEqual(len(res.time_ps), 501)  # 100 / 0.2 + 1
        self.assertEqual(len(res.reflectivity_change), 501)
        self.assertGreater(res.r_pump_magnitude, 0)
        self.assertGreater(res.r_probe_magnitude, 0)

    def test_legacy_inp_load(self):
        if Path("r_t_6.inp").exists():
            sim = load_simulation_from_legacy_inp("r_t_6.inp")
            res = sim.run()
            self.assertGreater(len(res.time_ps), 0)


if __name__ == "__main__":
    unittest.main()
