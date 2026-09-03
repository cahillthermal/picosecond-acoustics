"""
Material class and loader for picosecond acoustics simulation.
"""

from dataclasses import dataclass
from typing import Union, Dict, Any, Optional
from pathlib import Path


@dataclass
class Material:
    """
    Material properties for acoustic and optical simulation.

    Attributes:
        name: Name of the material (e.g., 'Al', 'Si', 'SiO2').
        density: Density in g/cm^3 or equivalent units (e.g. 2.7 for Al).
        sound_velocity: Sound velocity in nm/ps or km/s (e.g. 64.2 for Al).
        n_pump: Real part of pump refractive index.
        kappa_pump: Imaginary part of pump refractive index (extinction coefficient).
        n_probe: Real part of probe refractive index.
        kappa_probe: Imaginary part of probe refractive index.
        d_n_d_strain: Strain derivative of real refractive index dn/d_eta.
        d_kappa_d_strain: Strain derivative of imaginary refractive index d_kappa/d_eta.
        thermal_expansion: Thermal expansion coefficient divided by specific heat (alpha).
        diffusion_length: Thermal/carrier diffusion length in bin units or nm.
        attenuation: Acoustic attenuation parameter divided by omega^2.
    """
    name: str
    density: float
    sound_velocity: float
    n_pump: float
    kappa_pump: float
    n_probe: float
    kappa_probe: float
    d_n_d_strain: float
    d_kappa_d_strain: float
    thermal_expansion: float
    diffusion_length: float
    attenuation: float

    @property
    def epsilon_pump(self) -> complex:
        """Dielectric constant for pump light: (n + i*kappa)^2."""
        return (self.n_pump + 1j * self.kappa_pump) ** 2

    @property
    def epsilon_probe(self) -> complex:
        """Dielectric constant for probe light: (n + i*kappa)^2."""
        return (self.n_probe + 1j * self.kappa_probe) ** 2

    @property
    def d_epsilon_d_strain(self) -> complex:
        """Strain derivative of probe dielectric constant: 2*(n + i*kappa)*(dn/d_eta + i*d_kappa/d_eta)."""
        probe_index = self.n_probe + 1j * self.kappa_probe
        d_index_d_strain = self.d_n_d_strain + 1j * self.d_kappa_d_strain
        return 2.0 * probe_index * d_index_d_strain

    @property
    def acoustic_impedance(self) -> float:
        """Acoustic impedance Z = density * sound_velocity."""
        return self.density * self.sound_velocity

    @classmethod
    def from_inp_file(cls, filepath: Union[str, Path]) -> "Material":
        """
        Parse material parameters from a legacy MATLAB .inp file.

        Format:
        Line 1: density
        Line 2: sound velocity
        Line 3: pump n, kappa
        Line 4: probe n, kappa
        Line 5: probe d_n/d_strain, d_kappa/d_strain
        Line 6: expansion coefficient
        Line 7: diffusion length
        Line 8: attenuation coefficient
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Material input file not found: {filepath}")

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = [line.strip() for line in f if line.strip()]

        def parse_floats(line_str: str):
            # Strip potential comments in parentheses or text
            clean_str = line_str.split('(')[0].strip()
            parts = clean_str.split()
            return [float(p) for p in parts if p.replace('.', '', 1).replace('-', '', 1).replace('e', '', 1).replace('E', '', 1).isdigit() or p.replace('-', '', 1).isdigit()]

        # Helper to extract numbers strictly
        import re
        def extract_numbers(line_str: str):
            return [float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line_str)]

        nums = [extract_numbers(line) for line in lines]

        density = nums[0][0]
        velocity = nums[1][0]
        n_pump, kappa_pump = nums[2][0], nums[2][1]
        n_probe, kappa_probe = nums[3][0], nums[3][1]
        d_n, d_kappa = nums[4][0], nums[4][1]
        alpha = nums[5][0]
        diffusion = nums[6][0]
        attenuation = nums[7][0]

        return cls(
            name=path.stem,
            density=density,
            sound_velocity=velocity,
            n_pump=n_pump,
            kappa_pump=kappa_pump,
            n_probe=n_probe,
            kappa_probe=kappa_probe,
            d_n_d_strain=d_n,
            d_kappa_d_strain=d_kappa,
            thermal_expansion=alpha,
            diffusion_length=diffusion,
            attenuation=attenuation
        )
