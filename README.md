# Picosecond Acoustics Simulation in Python

A fast, vectorized Python simulation package for picosecond ultrasonics / picosecond acoustics in thin-film multilayer structures.

This project simulates:
1. **Optical Pump Absorption**: Transfer matrix method for pump light, energy deposition, and initial strain profile generation including electron/heat diffusion.
2. **Probe Sensitivity Function**: Optical perturbation matrices and depth-dependent sensitivity function $f_{\text{sens}}(z)$ for probe light reflectivity changes.
3. **Acoustic Wave Propagation**: Acoustic boundary conditions, transmission/reflection coefficients, and Generalized Stepping Algorithm (GSA) time-stepping.
4. **Thickness Fluctuations & Roughness**: Monte Carlo sampling across thickness variations.

---

## Features & Optimizations

- **High Speed**: Replaced MATLAB loops with vectorized NumPy array slicing and matrix multiplication (`@`). Simulations run in under 0.2 seconds.
- **Modern Object-Oriented Architecture**: Clean dataclasses for materials, layers, structures, optical beams, and simulation parameters.
- **Legacy Compatibility**: Built-in parser [`load_simulation_from_legacy_inp()`](picosecond_acoustics/io_legacy.py:27) to directly run existing legacy MATLAB `.inp` simulation input files.
- **Interactive Plotting**: Command-line flag `--plot` (`-p`) to display simulation plots before exiting, and `--save-plot` to export plot images.
- **Extensible & Maintainable**: Easily modify layer stacks, material properties, or numerical parameters directly in Python scripts.

---

## Installation

Ensure Python 3.8+ and dependencies are installed:

```bash
pip install numpy matplotlib
```

---

## Quick Start

### 1. Running a Legacy Simulation File via CLI

Run a simulation using an existing `.inp` file (e.g., [`r_t_6.inp`](r_t_6.inp:1)):

```bash
# Run with defaults (loads r_t_6.inp, saves r_t_py.dat, saves simulation_result.png, and displays plot):
python run_simulation.py

# Custom input file and output file:
python run_simulation.py r_t_6.inp --output r_t_py.dat --save-plot plot.png

# Run headless without displaying plot window:
python run_simulation.py r_t_6.inp --no-plot
```

#### CLI Command Options for [`run_simulation.py`](run_simulation.py:1):
- `inp_file`: Path to legacy input file (default: `r_t_6.inp`).
- `--output <filename>`, `-o <filename>`: Path to save numerical simulation result data text file (default: `r_t_py.dat`).
- `--save-plot <filename>`: Path to save plot image file (default: `simulation_result.png`).
- `--no-save-plot`: Disable saving the plot image to file.
- `--no-plot`: Disable displaying the interactive plot window.

---

## Input & Parameter File Specifications

To run a simulation via legacy `.inp` files, you prepare a main input file (e.g., [`r_t_6.inp`](r_t_6.inp:1)) along with individual material parameter files for each film in your stack (e.g., [`Al.inp`](Al.inp:1), [`Si.inp`](Si.inp:1), [`pmma.inp`](pmma.inp:1)).

### 1. Main Input File (`r_t_6.inp`)

The main input file configures time stepping, optics, multilayer film structure, and roughness Monte Carlo parameters line by line:

1. **`Time_step`** *(float)*: Data acquisition time step $\tau_{\text{data}}$ in picoseconds (ps).
2. **`Ratio of data to simulation`** *(float/int)*: Number of simulation sub-steps per data time step (`n_tau_divide`, usually `1.0`).
3. **`Strain save interval`** *(float)*: Time interval (ps) at which to snapshot spatial strain distributions vs. depth ($\eta(z)$). Recorded in `strain_snapshots` (or concatenated in legacy `r_t_eta.dat`). *Note: Avoid setting this interval too small to prevent simulation slowdowns.*
4. **`Time to stop`** *(float)*: Total simulation run time $t_{\text{stop}}$ in picoseconds (ps).
5. **`Pump Beam Parameters`** *(float float int)*: Vacuum wavelength ($\text{\AA}$), incidence angle ($^\circ$), and polarization (`1` = perpendicular to film / s-pol, `2` = parallel / p-pol).
6. **`Probe Beam Parameters`** *(float float int)*: Vacuum wavelength ($\text{\AA}$), incidence angle ($^\circ$), and polarization (`1` = perpendicular / s-pol, `2` = parallel / p-pol).
7. **`Number of films`** *(int)*: Total number of layers in the stack. *Recommendation: Include a nominal layer thickness (e.g., $100\,\text{\AA}$) of substrate material as the final layer to ensure accurate substrate boundary reflection coefficients.*
8. **`Configurations to average over`** *(int)*: Number of random roughness configurations to average ($N_{\text{samples}}$). Set to `1` for nominal thickness; enter e.g. `3` to run Monte Carlo sampling over random thickness variations and average $\Delta R / R$.
9. **`Film Layers`** *(1 line per film)*: `Thickness` ($\text{\AA}$), `Roughness amplitude` ($\text{\AA}$), and `Material file` name (e.g., `760.0  0.0  Al.inp`).

### 2. Material Parameter Files (`.inp`, e.g. `Al.inp`, `Si.inp`, `pmma.inp`)

Each film layer references a material parameter file formatted as follows:

1. **`Density`** *(float)*: Density in $\text{g/cm}^3$ ($\text{g}\cdot 10^{-24}/\text{\AA}^3$).
2. **`Sound Velocity`** *(float)*: Acoustic velocity in $\text{\AA/ps}$ ($\text{km/s}$).
3. **`Pump Refractive Index`** *(float float)*: Real ($n$) and imaginary ($\kappa$) refractive index at the pump wavelength.
4. **`Probe Refractive Index`** *(float float)*: Real ($n$) and imaginary ($\kappa$) refractive index at the probe wavelength.
5. **`Photoelastic Derivatives`** *(float float)*: Derivatives of $n$ and $\kappa$ with respect to strain ($\partial n / \partial \eta$, $\partial \kappa / \partial \eta$), used as fitting parameters for acoustic echo pulse shapes.
6. **`Thermal Expansion`** *(float)*: Linear thermal expansion coefficient ($\times 10^{-6}\,\text{K}^{-1}$).
7. **`Electron Diffusion Length`** *(float)*: Electronic energy diffusion length in Angstroms ($\text{\AA}$). Expands initial strain pulse beyond pure optical absorption depth before energy transfers to the lattice.
8. **`Attenuation Coefficient`** *(float)*: Acoustic attenuation rate in $\text{ps}^{-1}$.

### 3. Custom Python Script Example

You can define materials, layer stacks, optical beams, and run simulations directly in Python:

```python
from picosecond_acoustics import (
    Material,
    FilmLayer,
    Structure,
    LightBeam,
    SimulationConfig,
    PicosecondAcousticsSimulation,
)

# Load or define materials
al = Material.from_inp_file("Al.inp")
si = Material.from_inp_file("Si.inp")

# Define multilayer structure
structure = Structure()
structure.add_layer(FilmLayer(material=al, thickness_mean=760.0, roughness=0.0))
structure.add_layer(FilmLayer(material=si, thickness_mean=7000.0, roughness=0.0))

# Configure optics
pump = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization="s")
probe = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization="p")

# Configure simulation time-stepping
config = SimulationConfig(
    tau_data=0.2,       # Data step in ps
    n_tau_divide=1,     # Sub-steps per data step
    t_stop=700.0,       # Simulation stop time in ps
    t_save=100.0,       # Snapshot save interval
    n_samples=1         # Roughness samples
)

# Instantiate and run
sim = PicosecondAcousticsSimulation(
    structure=structure,
    pump_beam=pump,
    probe_beam=probe,
    sim_config=config
)

result = sim.run()

# Save results
result.save_txt("output.dat")
print(f"Calculated {len(result.time_ps)} time steps.")
```

---

## Running Unit Tests

Run the included test suite:

```bash
python -m unittest discover tests
```

---

## Repository Structure

- [`picosecond_acoustics/`](picosecond_acoustics/__init__.py:1) - Main Python package module
  - [`material.py`](picosecond_acoustics/material.py:1) - [`Material`](picosecond_acoustics/material.py:10) dataclass and `.inp` loader
  - [`layer.py`](picosecond_acoustics/layer.py:1) - [`FilmLayer`](picosecond_acoustics/layer.py:10) and [`Structure`](picosecond_acoustics/layer.py:35) representation
  - [`config.py`](picosecond_acoustics/config.py:1) - [`LightBeam`](picosecond_acoustics/config.py:10) and [`SimulationConfig`](picosecond_acoustics/config.py:50) configuration
  - [`optics.py`](picosecond_acoustics/optics.py:1) - Transfer matrix method, absorption, diffusion, sensitivity function
  - [`acoustics.py`](picosecond_acoustics/acoustics.py:1) - Acoustic discretization, boundary coefficients, and GSA solver
  - [`simulation.py`](picosecond_acoustics/simulation.py:1) - Simulation driver and [`SimulationResult`](picosecond_acoustics/simulation.py:18) container
  - [`io_legacy.py`](picosecond_acoustics/io_legacy.py:1) - Legacy MATLAB `.inp` parser
- [`matlab_legacy_code/`](matlab_legacy_code/) - Original legacy MATLAB source scripts (`.m`) preserved for reference
- [`run_simulation.py`](run_simulation.py:1) - CLI script to execute simulations
- [`tests/`](tests/test_simulation.py:1) - Unit tests

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [`LICENSE`](LICENSE:1) file for details.
