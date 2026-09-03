"""
High-level simulation driver for picosecond acoustics.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any
import numpy as np

from .config import LightBeam, SimulationConfig
from .layer import Structure
from .optics import (
    compute_beam_matrices,
    compute_stack_transfer_matrix,
    compute_initial_strain,
    compute_sensitivity_function
)
from .acoustics import AcousticSolver


@dataclass
class SimulationResult:
    """
    Simulation output container.
    """
    time_ps: np.ndarray
    reflectivity_change: np.ndarray
    z_coords: np.ndarray
    sensitivity_function: np.ndarray
    strain_snapshots: Dict[float, np.ndarray] = field(default_factory=dict)
    r_pump_magnitude: float = 0.0
    r_probe_magnitude: float = 0.0
    sampled_layer_thicknesses: List[List[float]] = field(default_factory=list)

    def save_txt(self, filepath: str) -> None:
        """Save time vs reflectivity change to text file."""
        data = np.column_stack((self.time_ps, self.reflectivity_change))
        np.savetxt(filepath, data, fmt="%.6e   %.6e", header="Time_ps   dR/R")


class PicosecondAcousticsSimulation:
    """
    Main simulator class for picosecond acoustic response in multilayer structures.
    """

    def __init__(
        self,
        structure: Structure,
        pump_beam: LightBeam,
        probe_beam: LightBeam,
        sim_config: SimulationConfig
    ):
        self.structure = structure
        self.pump = pump_beam
        self.probe = probe_beam
        self.config = sim_config

    def run(self, seed: Optional[int] = 42) -> SimulationResult:
        """
        Run the simulation with optional roughness Monte Carlo sampling.
        """
        if seed is not None:
            np.random.seed(seed)

        num_data_steps = self.config.num_data_steps
        time_ps = np.arange(num_data_steps + 1) * self.config.tau_data
        total_r_t = np.zeros(num_data_steps + 1)

        all_sampled_thicknesses = []
        last_z_coords = None
        last_f_sens = None
        r_pump_mag = 0.0
        r_probe_mag = 0.0
        strain_snapshots = {}

        for sample_idx in range(self.config.n_samples):
            # Sample layer thicknesses
            if self.config.n_samples == 1:
                rnd_vals = [0.0] * self.structure.num_layers
            else:
                rnd_vals = np.random.normal(size=self.structure.num_layers).tolist()

            thicknesses = self.structure.sample_structure(rnd_vals)
            all_sampled_thicknesses.append(thicknesses)

            # 1. Setup acoustic discretization & GSA solver
            solver = AcousticSolver(
                structure=self.structure,
                layer_thicknesses=thicknesses,
                tau_sim=self.config.tau_sim
            )
            last_z_coords = solver.z_coords

            # 2. Pump optical matrix calculation
            m_pump_bins, m_pump_end, _ = compute_beam_matrices(
                beam=self.pump,
                structure=self.structure,
                bin_sizes=solver.bin_sizes,
                is_pump=True
            )
            m_pump_tot, r_pump, max_pump_bin, _ = compute_stack_transfer_matrix(
                m_bins=m_pump_bins,
                m_end=m_pump_end,
                n_bins_per_film=solver.n_bins_per_film
            )
            r_pump_mag = abs(r_pump)**2

            # 3. Initial strain calculation
            eta_left, eta_right = compute_initial_strain(
                structure=self.structure,
                n_bins_per_film=solver.n_bins_per_film,
                bin_sizes=solver.bin_sizes,
                r_pump=r_pump,
                m_pump_bins=m_pump_bins,
                max_pump_bin=max_pump_bin
            )

            # 4. Probe optical matrix & sensitivity function calculation
            m_probe_bins, m_probe_end, dm_probe_bins = compute_beam_matrices(
                beam=self.probe,
                structure=self.structure,
                bin_sizes=solver.bin_sizes,
                is_pump=False
            )
            m_probe_tot, r_probe, max_probe_bin, _ = compute_stack_transfer_matrix(
                m_bins=m_probe_bins,
                m_end=m_probe_end,
                n_bins_per_film=solver.n_bins_per_film
            )
            r_probe_mag = abs(r_probe)**2

            f_sens = compute_sensitivity_function(
                structure=self.structure,
                n_bins_per_film=solver.n_bins_per_film,
                m_probe_bins=m_probe_bins,
                dm_probe_bins=dm_probe_bins,
                m_probe_end=m_probe_end,
                m_probe_total=m_probe_tot,
                r_probe=r_probe,
                max_probe_bin=max_probe_bin
            )
            last_f_sens = f_sens

            # 5. Time stepping loop
            sample_r_t = np.zeros(num_data_steps + 1)
            next_save_time = 0.0

            for i_data in range(num_data_steps + 1):
                t_curr = i_data * self.config.tau_data

                # Save strain profile snapshot if requested
                if sample_idx == 0 and abs(t_curr - next_save_time) < 0.5 * self.config.tau_sim:
                    strain_snapshots[t_curr] = eta_left + eta_right
                    next_save_time += self.config.t_save

                # Reflectivity change at current strain profile
                # Sum over probe light penetration depth bins
                sample_r_t[i_data] = np.dot(
                    eta_left[:max_probe_bin] + eta_right[:max_probe_bin],
                    f_sens[:max_probe_bin]
                )

                # Step strain forward by n_tau_divide sub-steps
                for _ in range(self.config.n_tau_divide):
                    eta_left, eta_right = solver.step(eta_left, eta_right)

            total_r_t += sample_r_t

        avg_r_t = total_r_t / float(self.config.n_samples)

        return SimulationResult(
            time_ps=time_ps,
            reflectivity_change=avg_r_t,
            z_coords=last_z_coords,
            sensitivity_function=last_f_sens,
            strain_snapshots=strain_snapshots,
            r_pump_magnitude=r_pump_mag,
            r_probe_magnitude=r_probe_mag,
            sampled_layer_thicknesses=all_sampled_thicknesses
        )
