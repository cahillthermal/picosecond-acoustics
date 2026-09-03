"""
Optical calculations: 2x2 transfer matrix method, pump absorption,
initial thermal strain generation with diffusion, and probe sensitivity functions.
"""

from typing import Tuple, List, Dict
import numpy as np

from .config import LightBeam
from .layer import Structure
from .material import Material


def compute_layer_bin_matrices(
    beam: LightBeam,
    materials: List[Material],
    bin_sizes: List[float]
) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Compute 2x2 bin transfer matrix for each film layer and the end transfer matrix.

    Returns:
        m_bins: List of 2x2 complex numpy arrays (one per film layer).
        m_end: 2x2 complex numpy array for the exit/substrate interface.
    """
    k_z_beam = []
    c = np.cos(beam.angle_rad)
    s = np.sin(beam.angle_rad)
    i_pol = beam.polarization_code

    m_bins = []
    m_end = np.identity(2, dtype=complex)

    for film_idx, (mat, bin_sz) in enumerate(zip(materials, bin_sizes)):
        eps = mat.epsilon_pump if beam.wavelength == beam.wavelength else mat.epsilon_probe
        # We pass epsilon explicitly below, but for generic beam:
        # Use mat.epsilon_pump if beam is pump or mat.epsilon_probe if probe.
        # Handled in caller functions.
        pass

    raise NotImplementedError("Use compute_beam_matrices")


def compute_beam_matrices(
    beam: LightBeam,
    structure: Structure,
    bin_sizes: List[float],
    is_pump: bool = True
) -> Tuple[List[np.ndarray], np.ndarray, List[np.ndarray]]:
    """
    Compute bin transfer matrices for pump or probe light for each film.
    Also returns derivative matrices dM_bin per film for probe light.

    Returns:
        m_bins: List of 2x2 complex arrays (one per film layer).
        m_end: 2x2 complex array for exit boundary.
        dm_bins: List of 2x2 complex arrays (dM/d_eta derivative per film layer).
    """
    c = np.cos(beam.angle_rad)
    s = np.sin(beam.angle_rad)
    i_pol = beam.polarization_code
    k0 = beam.wavevector

    m_bins = []
    dm_bins = []
    m_end = np.identity(2, dtype=complex)

    num_films = structure.num_layers
    strain_pert = 0.01  # Small strain value for numerical derivative

    for film_idx, (layer, bin_sz) in enumerate(zip(structure.layers, bin_sizes)):
        mat = layer.material
        eps = mat.epsilon_pump if is_pump else mat.epsilon_probe

        c1 = np.sqrt(eps - s**2 + 0j)
        c2 = np.sqrt(eps + 0j)
        c3 = eps * c
        kz = k0 * c1

        if i_pol == 1:  # s-polarization
            t_in = 2.0 * c / (c + c1)
            t_out = 2.0 * c1 / (c + c1)
            r_in = (c - c1) / (c + c1)
            r_out = -r_in
        else:  # p-polarization
            t_in = 2.0 * c2 * c / (eps * c + c1)
            t_out = 2.0 * c2 * c1 / (eps * c + c1)
            r_in = (c3 - c1) / (c3 + c1)
            r_out = -r_in

        c4 = np.exp(1j * kz * bin_sz)
        c5 = c4 * c4

        t = t_in * t_out * c4 / (1.0 - r_out**2 * c5)
        r = r_in * (1.0 - c5) / (1.0 - r_out**2 * c5)

        m_bin = np.array([
            [t - r**2 / t, r / t],
            [-r / t, 1.0 / t]
        ], dtype=complex)
        m_bins.append(m_bin)

        if film_idx == num_films - 1:
            m_end = np.array([
                [t_in - r_out * r_in / t_out, r_out / t_out],
                [-r_in / t_out, 1.0 / t_out]
            ], dtype=complex)

        # Calculate perturbed matrix dM for probe light
        if not is_pump:
            eps_prime = eps + mat.d_epsilon_d_strain * strain_pert
            bin_sz_prime = bin_sz * (1.0 + strain_pert)

            c1_p = np.sqrt(eps_prime - s**2 + 0j)
            c2_p = np.sqrt(eps_prime + 0j)
            c3_p = eps_prime * c
            kz_p = k0 * c1_p

            if i_pol == 1:
                t_in_p = 2.0 * c / (c + c1_p)
                t_out_p = 2.0 * c1_p / (c + c1_p)
                r_in_p = (c - c1_p) / (c + c1_p)
                r_out_p = -r_in_p
            else:
                t_in_p = 2.0 * c2_p * c / (eps_prime * c + c1_p)
                t_out_p = 2.0 * c2_p * c1_p / (eps_prime * c + c1_p)
                r_in_p = (c3_p - c1_p) / (c3_p + c1_p)
                r_out_p = -r_in_p

            c4_p = np.exp(1j * kz_p * bin_sz_prime)
            c5_p = c4_p * c4_p

            t_p = t_in_p * t_out_p * c4_p / (1.0 - r_out_p**2 * c5_p)
            r_p = r_in_p * (1.0 - c5_p) / (1.0 - r_out_p**2 * c5_p)

            m_bin_prime = np.array([
                [t_p - r_p**2 / t_p, r_p / t_p],
                [-r_p / t_p, 1.0 / t_p]
            ], dtype=complex)

            dm_bin = (m_bin_prime - m_bin) / strain_pert
            dm_bins.append(dm_bin)

    return m_bins, m_end, dm_bins


def compute_stack_transfer_matrix(
    m_bins: List[np.ndarray],
    m_end: np.ndarray,
    n_bins_per_film: List[int]
) -> Tuple[np.ndarray, complex, int, List[np.ndarray]]:
    """
    Multiply bin matrices sequentially from surface to end.
    Includes numerical trap for highly absorbing opaque films.

    Returns:
        m_total: 2x2 total transfer matrix.
        r_amplitude: Amplitude reflection coefficient r = -M[1,0] / M[1,1].
        max_bin: Index of max bin reached before opacity cutoff.
        m_cumulative: Cumulative product matrices for each bin.
    """
    m_curr = np.identity(2, dtype=complex)
    m_cumulative = [m_curr.copy()]
    max_bin = 0
    total_bins = sum(n_bins_per_film)
    trap_triggered = False

    bin_counter = 0
    for film_idx, n_bins in enumerate(n_bins_per_film):
        m_bin = m_bins[film_idx]
        for _ in range(n_bins):
            m_curr = m_bin @ m_curr
            bin_counter += 1
            m_cumulative.append(m_curr.copy())

            # Absorption trap check: |M[0,0] - M[1,0]*M[0,1]/M[1,1]| <= 1.0e-3
            val = abs(m_curr[0, 0] - m_curr[1, 0] * m_curr[0, 1] / m_curr[1, 1])
            if val <= 1.0e-3:
                max_bin = bin_counter
                trap_triggered = True
                break
        if trap_triggered:
            break

    if not trap_triggered:
        max_bin = total_bins
        m_total = m_end @ m_curr
    else:
        m_total = m_curr

    r_amplitude = -m_total[1, 0] / m_total[1, 1]
    return m_total, r_amplitude, max_bin, m_cumulative


def compute_initial_strain(
    structure: Structure,
    n_bins_per_film: List[int],
    bin_sizes: List[float],
    r_pump: complex,
    m_pump_bins: List[np.ndarray],
    max_pump_bin: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate initial strain distribution (left- and right-propagating waves)
    due to pump light absorption and thermal/carrier diffusion.

    Returns:
        eta_left: Initial strain array going left.
        eta_right: Initial strain array going right.
    """
    total_bins = sum(n_bins_per_film)
    eta_raw = np.zeros(total_bins)

    a1_left = 1.0 + 0j
    a2_left = r_pump

    bin_idx = 0
    for film_idx, (layer, n_bins) in enumerate(zip(structure.layers, n_bins_per_film)):
        mat = layer.material
        bin_sz = bin_sizes[film_idx]
        m_bin = m_pump_bins[film_idx]

        for _ in range(n_bins):
            if bin_idx < max_pump_bin:
                a1_right = m_bin[0, 0] * a1_left + m_bin[0, 1] * a2_left
                a2_right = m_bin[1, 0] * a1_left + m_bin[1, 1] * a2_left

                energy = (abs(a1_left)**2 - abs(a1_right)**2 +
                          abs(a2_right)**2 - abs(a2_left)**2)
                eta_raw[bin_idx] = -0.5 * mat.thermal_expansion * energy / bin_sz

                a1_left = a1_right
                a2_left = a2_right
            else:
                eta_raw[bin_idx] = 0.0

            bin_idx += 1

    # Apply spatial diffusion within each film layer
    eta_diffused = np.zeros(total_bins)
    x_max = 40.0

    film_start_bin = 0
    for film_idx, (layer, n_bins) in enumerate(zip(structure.layers, n_bins_per_film)):
        mat = layer.material
        bin_sz = bin_sizes[film_idx]
        diff_len = mat.diffusion_length

        x0 = x_max if diff_len == 0.0 else bin_sz / diff_len

        for local_bin in range(n_bins):
            curr_bin = film_start_bin + local_bin
            if eta_raw[curr_bin] == 0.0:
                continue

            # Construct spread profile across this film
            spread = np.zeros(n_bins)
            for p_bin in range(n_bins):
                p_dist = p_bin + 1
                b_dist = local_bin + 1

                if p_dist <= b_dist:
                    if b_dist * x0 < x_max:
                        spread[p_bin] = np.cosh(p_dist * x0) / np.cosh(b_dist * x0)
                    else:
                        x1 = min(x_max, (b_dist - p_dist) * x0)
                        spread[p_bin] = np.exp(-x1)
                else:
                    rem_b = n_bins - local_bin
                    rem_p = n_bins - p_bin
                    if rem_b * x0 < x_max:
                        spread[p_bin] = np.cosh(rem_p * x0) / np.cosh(rem_b * x0)
                    else:
                        x1 = min(x_max, (p_bin - local_bin) * x0)
                        spread[p_bin] = np.exp(-x1)

            # Normalize spread function across film
            spread_sum = np.sum(spread)
            if spread_sum > 0:
                spread /= spread_sum

            # Accumulate diffused strain
            eta_diffused[film_start_bin:film_start_bin + n_bins] += eta_raw[curr_bin] * spread

        film_start_bin += n_bins

    eta_left = eta_diffused.copy()
    eta_right = eta_diffused.copy()
    return eta_left, eta_right


def compute_sensitivity_function(
    structure: Structure,
    n_bins_per_film: List[int],
    m_probe_bins: List[np.ndarray],
    dm_probe_bins: List[np.ndarray],
    m_probe_end: np.ndarray,
    m_probe_total: np.ndarray,
    r_probe: complex,
    max_probe_bin: int
) -> np.ndarray:
    """
    Compute sensitivity function f_sens(z) for probe reflectivity changes.

    Returns:
        f_sens: 1D float array of size total_bins.
    """
    total_bins = sum(n_bins_per_film)
    f_sens = np.zeros(total_bins)

    # Build backward chain M_A
    m_a = [np.identity(2, dtype=complex) for _ in range(total_bins + 2)]
    if max_probe_bin == total_bins:
        m_a[total_bins] = m_probe_end.copy()
    else:
        m_a[total_bins] = np.identity(2, dtype=complex)

    # Propagate M_A backwards
    bin_counter = total_bins - 1
    for film_idx in range(structure.num_layers - 1, -1, -1):
        n_bins = n_bins_per_film[film_idx]
        m_bin = m_probe_bins[film_idx]
        for _ in range(n_bins):
            if bin_counter < max_probe_bin:
                m_a[bin_counter] = m_a[bin_counter + 1] @ m_bin
            else:
                m_a[bin_counter] = m_a[bin_counter + 1].copy()
            bin_counter -= 1

    # Build forward chain M_B
    m_b = [np.identity(2, dtype=complex) for _ in range(total_bins + 1)]
    bin_counter = 0
    for film_idx in range(structure.num_layers):
        n_bins = n_bins_per_film[film_idx]
        m_bin = m_probe_bins[film_idx]
        for _ in range(n_bins):
            if bin_counter < max_probe_bin:
                m_b[bin_counter + 1] = m_bin @ m_b[bin_counter]
            else:
                m_b[bin_counter + 1] = m_b[bin_counter].copy()
            bin_counter += 1

    # Compute sensitivity f_sens per bin
    m_22_sq = m_probe_total[1, 1]**2
    bin_counter = 0
    for film_idx in range(structure.num_layers):
        n_bins = n_bins_per_film[film_idx]
        dm_bin = dm_probe_bins[film_idx]
        for _ in range(n_bins):
            if bin_counter < max_probe_bin:
                m_work = dm_bin @ m_b[bin_counter]
                dm_total = m_a[bin_counter + 1] @ m_work

                dr_probe = (
                    m_probe_total[1, 0] * dm_total[1, 1] -
                    m_probe_total[1, 1] * dm_total[1, 0]
                ) / m_22_sq

                f_sens[bin_counter] = 2.0 * np.real(r_probe * np.conj(dr_probe))
            else:
                f_sens[bin_counter] = 0.0

            bin_counter += 1

    return f_sens
