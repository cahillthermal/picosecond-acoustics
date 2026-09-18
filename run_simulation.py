#!/usr/bin/env python3
"""
CLI entry point to run picosecond acoustics simulations.
Usage:
    python run_simulation.py [path_to_inp_file] [--output output.dat] [--save-plot plot.png] [--no-plot] [--no-save-plot]
"""

import sys
import argparse
import time
from pathlib import Path

from picosecond_acoustics import (
    load_simulation_from_legacy_inp,
    SimulationResult
)


def display_plot(result: SimulationResult, title: str = "Picosecond Acoustics Simulation", save_path: str = None, show: bool = True) -> None:
    """Display and/or save plot of reflectivity change vs time."""
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

        if show:
            try:
                plt.show()
            except Exception:
                pass
        plt.close()
    except ImportError:
        print("Warning: matplotlib is not installed. Unable to display or save plot.")
    except Exception as e:
        print(f"Display error: {e}")


def main():
    parser = argparse.ArgumentParser(description="Picosecond Acoustics Simulator")
    parser.add_argument("inp_file", nargs="?", default="r_t_6.inp", help="Path to legacy .inp file (default: r_t_6.inp)")
    parser.add_argument("--output", "-o", default="r_t_py.dat", help="Path to save output data (default: r_t_py.dat)")
    parser.add_argument("--save-plot", default="simulation_result.png", help="Path to save plot image (default: simulation_result.png)")
    parser.add_argument("--no-save-plot", action="store_true", help="Do not save plot image to file")
    parser.add_argument("--no-plot", action="store_true", help="Do not display interactive plot window")
    parser.add_argument("--plot", "-p", dest="show_plot", action="store_true", default=True, help="Display plot of reflectivity change (default: True)")

    args = parser.parse_args()

    inp_path = Path(args.inp_file)
    if not inp_path.exists() and not inp_path.is_absolute():
        script_dir_path = Path(__file__).resolve().parent / args.inp_file
        if script_dir_path.exists():
            inp_path = script_dir_path

    if not inp_path.exists():
        print(f"Error: Input file '{args.inp_file}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading simulation from: {inp_path}")
    sim = load_simulation_from_legacy_inp(inp_path)
    
    t0 = time.time()
    result = sim.run()
    t1 = time.time()

    print(f"Simulation finished in {t1 - t0:.4f} seconds.")
    print(f"Saving output data to: {args.output}")
    result.save_txt(args.output)

    save_path = None if args.no_save_plot else args.save_plot
    show = False if args.no_plot else args.show_plot

    if show or save_path:
        display_plot(result, title=f"Picosecond Acoustics Simulation ({inp_path.name})", save_path=save_path, show=show)


if __name__ == "__main__":
    main()
