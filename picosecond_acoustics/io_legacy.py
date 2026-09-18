"""
Legacy loader for MATLAB simulation input files (.inp format).
"""

import re
from pathlib import Path
from typing import Tuple, Union, List

from .config import LightBeam, SimulationConfig
from .layer import FilmLayer, Structure
from .material import Material
from .simulation import PicosecondAcousticsSimulation


def parse_numbers_and_text(line: str) -> Tuple[List[float], str]:
    """Extract all float numbers from line and any non-numeric trailing text."""
    clean = line.strip()
    # Find trailing word filename if present (e.g., Al.inp)
    match_file = re.search(r"([A-Za-z0-9_\-]+\.inp)$", clean, re.IGNORECASE)
    filename = match_file.group(1) if match_file else ""

    nums = [float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", clean)]
    return nums, filename


def load_simulation_from_legacy_inp(
    inp_filepath: Union[str, Path]
) -> PicosecondAcousticsSimulation:
    """
    Load a complete PicosecondAcousticsSimulation instance from a legacy .inp file.
    """
    inp_path = Path(inp_filepath)
    base_dir = inp_path.parent

    with open(inp_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]

    tau_data = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[0])[0])
    ratio = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[1])[0])
    n_tau_divide = max(1, int(round(ratio)))
    t_save = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[2])[0])
    t_stop = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[3])[0])

    pump_vals = [float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[4])]
    pump_beam = LightBeam(
        wavelength=pump_vals[0],
        angle_deg=pump_vals[1],
        polarization=int(pump_vals[2])
    )

    probe_vals = [float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[5])]
    probe_beam = LightBeam(
        wavelength=probe_vals[0],
        angle_deg=probe_vals[1],
        polarization=int(probe_vals[2])
    )

    num_films_spec = int(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[6])[0])
    n_samples = int(re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", lines[7])[0])

    # Parse film lines starting at line index 8
    film_lines = lines[8:]
    parsed_film_entries = []
    for line in film_lines:
        nums, mat_filename = parse_numbers_and_text(line)
        if not mat_filename:
            for token in line.split():
                if token.lower().endswith(".inp"):
                    mat_filename = token
                    break
        if nums and mat_filename:
            parsed_film_entries.append((nums, mat_filename))

    # Fallback to range(8, 8 + num_films_spec) if no pattern matched
    if not parsed_film_entries:
        for line_idx in range(8, min(len(lines), 8 + num_films_spec)):
            line = lines[line_idx]
            nums, mat_filename = parse_numbers_and_text(line)
            if not mat_filename:
                for token in line.split():
                    if token.lower().endswith(".inp"):
                        mat_filename = token
                        break
            if nums:
                parsed_film_entries.append((nums, mat_filename))

    structure = Structure()
    for nums, mat_filename in parsed_film_entries:
        thickness = nums[0]
        roughness = nums[1] if len(nums) > 1 else 0.0

        mat_path = base_dir / mat_filename
        if not mat_path.exists():
            # Try current directory
            mat_path = Path(mat_filename)

        material = Material.from_inp_file(mat_path)
        film_layer = FilmLayer(
            material=material,
            thickness_mean=thickness,
            roughness=roughness
        )
        structure.add_layer(film_layer)

    sim_config = SimulationConfig(
        tau_data=tau_data,
        n_tau_divide=n_tau_divide,
        t_stop=t_stop,
        t_save=t_save,
        n_samples=n_samples
    )

    return PicosecondAcousticsSimulation(
        structure=structure,
        pump_beam=pump_beam,
        probe_beam=probe_beam,
        sim_config=sim_config
    )
