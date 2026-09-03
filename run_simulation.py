#!/usr/bin/env python3
"""
CLI entry point to run picosecond acoustics simulations.
Usage:
    python run_simulation.py [path_to_inp_file] [--plot] [--output output.dat] [--save-plot plot.png]
"""

import sys
import argparse
import time
from pathlib import Path

from picosecond_acoustics import (
    load_simulation_from_legacy_inp,
    PicosecondAcousticsSimulation,
    Structure,
    FilmLayer,
    Material,
    LightBeam,
    SimulationConfig,
    SimulationResult
)


def display_plot(result: SimulationResult, title: str = "Picosecond Acoustics Simulation", save_path: str = None) -> None:
    """Display interactive plot of reflectivity change vs time."""
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 5))
        plt.plot(result.time_ps, result.reflectivity_change, 'b-', linewidth=1.5, label="dR/R (Reflectivity Change)")
        plt.xlabel("Time (ps)", fontsize=11)
        plt.ylabel("Reflectivity Change dR/R", fontsize=11)
        plt.title(title, fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.legend(loc="upper right")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"Plot saved to: {save_path}")

        try:
            plt.show()
        except Exception:
            pass
    except ImportError:
        print("Warning: matplotlib is not installed. Unable to display or save plot.")
    except Exception as e:
        print(f"Display error: {e}")


def run_demo_custom_structure() -> SimulationResult:
    """Run a sample simulation defined purely in Python code."""
    print("--- Running Demo Custom Structure Simulation ---")
    al = Material.from_inp_file("Al.inp")
    si = Material.from_inp_file("Si.inp")

    structure = Structure()
    structure.add_layer(FilmLayer(material=al, thickness_mean=760.0, roughness=0.0))
    structure.add_layer(FilmLayer(material=si, thickness_mean=7000.0, roughness=0.0))

    pump = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization='s')
    probe = LightBeam(wavelength=7850.0, angle_deg=0.0, polarization='p')
    config = SimulationConfig(tau_data=0.2, n_tau_divide=1, t_stop=700.0, t_save=100.0, n_samples=1)

    sim = PicosecondAcousticsSimulation(
        structure=structure,
        pump_beam=pump,
        probe_beam=probe,
        sim_config=config
    )

    t0 = time.time()
    result = sim.run()
    t1 = time.time()

    print(f"Simulation completed in {t1 - t0:.4f} seconds!")
    print(f"Time steps calculated: {len(result.time_ps)}")
    print(f"Pump Reflection Magnitude R_pump: {result.r_pump_magnitude:.6e}")
    print(f"Probe Reflection Magnitude R_probe: {result.r_probe_magnitude:.6e}")

    return result


def main():
    parser = argparse.ArgumentParser(description="Picosecond Acoustics Simulator")
    parser.add_argument("inp_file", nargs="?", default="r_t_6.inp", help="Path to legacy .inp file")
    parser.add_argument("--output", "-o", default="r_t_py.dat", help="Path to save output data")
    parser.add_argument("--plot", "-p", action="store_true", help="Display plot of simulation reflectivity change before exiting")
    parser.add_argument("--save-plot", default=None, help="Optional image filepath to save plot (e.g. plot.png)")

    args = parser.parse_args()

    inp_path = Path(args.inp_file)
    if inp_path.exists():
        print(f"Loading simulation from: {inp_path}")
        sim = load_simulation_from_legacy_inp(inp_path)
        
        t0 = time.time()
        result = sim.run()
        t1 = time.time()

        print(f"Simulation finished in {t1 - t0:.4f} seconds.")
        print(f"Saving output data to: {args.output}")
        result.save_txt(args.output)

        if args.plot or args.save_plot:
            display_plot(result, title=f"Picosecond Acoustics Simulation ({inp_path.name})", save_path=args.save_plot)
    else:
        print(f"Input file '{inp_path}' not found. Running custom python demo...")
        result = run_demo_custom_structure()
        result.save_txt(args.output)

        if args.plot or args.save_plot:
            display_plot(result, title="Picosecond Acoustics Simulation (Demo Structure)", save_path=args.save_plot)


if __name__ == "__main__":
    main()
