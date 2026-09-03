"""
Configuration dataclasses for pump/probe optical beams and simulation parameters.
"""

from dataclasses import dataclass
from typing import Literal, Union
import numpy as np


@dataclass
class LightBeam:
    """
    Optical beam configuration (Pump or Probe).

    Attributes:
        wavelength: Light wavelength in vacuum (same length units as layer thickness, e.g. Angstroms or nm).
        angle_deg: Angle of incidence in vacuum in degrees.
        polarization: Polarization state ('s' / 'sigma' or 'p' / 'pi', or numeric 1/2).
    """
    wavelength: float
    angle_deg: float = 0.0
    polarization: Union[str, int] = 's'

    @property
    def wavevector(self) -> float:
        """Free space wavevector k = 2 * pi / wavelength."""
        return 2.0 * np.pi / self.wavelength

    @property
    def angle_rad(self) -> float:
        """Angle of incidence in radians."""
        return np.radians(self.angle_deg)

    @property
    def polarization_code(self) -> int:
        """
        1 for s-polarization (sigma), 2 for p-polarization (pi).
        """
        if isinstance(self.polarization, int):
            return self.polarization
        pol = str(self.polarization).lower()
        if pol in ('s', 'sigma', '1'):
            return 1
        elif pol in ('p', 'pi', '2'):
            return 2
        else:
            raise ValueError(f"Unknown polarization: {self.polarization}")


@dataclass
class SimulationConfig:
    """
    Time stepping and output options for simulation.

    Attributes:
        tau_data: Time step in recorded data acquisition.
        n_tau_divide: Number of simulation sub-steps per data time step.
        t_stop: Total simulation end time.
        t_save: Time interval between saving spatial strain distribution snapshots.
        n_samples: Number of thickness fluctuation configurations to average.
    """
    tau_data: float
    n_tau_divide: int = 1
    t_stop: float = 700.0
    t_save: float = 100.0
    n_samples: int = 1

    @property
    def tau_sim(self) -> float:
        """Simulation numerical time step tau = tau_data / n_tau_divide."""
        return self.tau_data / float(self.n_tau_divide)

    @property
    def num_data_steps(self) -> int:
        """Total number of data time steps."""
        return int(np.round(self.t_stop / self.tau_data))
