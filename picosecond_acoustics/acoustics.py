"""
Acoustic discretization, reflection/transmission coefficients,
Generalized Stepping Algorithm (GSA) matrix generation, and fast vectorized strain time-stepping.
"""

from typing import Tuple, List
import numpy as np

from .layer import Structure


class AcousticSolver:
    """
    Solves acoustic wave propagation in multilayer structures using GSA.
    """

    def __init__(
        self,
        structure: Structure,
        layer_thicknesses: List[float],
        tau_sim: float
    ):
        self.structure = structure
        self.layer_thicknesses = layer_thicknesses
        self.tau_sim = tau_sim

        self.n_bins_per_film = []
        self.bin_sizes = []
        self._discretize()
        self._compute_boundary_coefficients()
        self._compute_gsa_coefficients()

    def _discretize(self) -> None:
        """Discretize film layers into spatial bins matching acoustic time step."""
        total_bins = 0
        bin_z_list = [0.0]
        densities = []
        velocities = []

        for film_idx, (layer, thick) in enumerate(zip(self.structure.layers, self.layer_thicknesses)):
            mat = layer.material
            vel = mat.sound_velocity
            bin_sz = self.tau_sim * vel

            # Number of bins in film (at least 1)
            n_bins = int(np.floor(thick / bin_sz + 0.5))
            if n_bins == 0:
                raise ValueError(
                    f"Film {film_idx} ({mat.name}) thickness {thick} is too small for bin size {bin_sz}"
                )

            adjusted_bin_sz = thick / n_bins
            self.n_bins_per_film.append(n_bins)
            self.bin_sizes.append(adjusted_bin_sz)

            for b in range(n_bins):
                total_bins += 1
                bin_z_list.append(bin_z_list[-1] + adjusted_bin_sz)
                densities.append(mat.density)
                velocities.append(vel)

        self.total_bins = total_bins
        self.z_coords = np.array(bin_z_list)
        self.bin_z_centers = 0.5 * (self.z_coords[:-1] + self.z_coords[1:])
        self.densities = np.array(densities)
        self.velocities = np.array(velocities)
        self.acoustic_impedances = self.densities * self.velocities

    def _compute_boundary_coefficients(self) -> None:
        """
        Compute acoustic reflection and transmission coefficients at boundaries.
        Left surface (interface 0): Free boundary (refl_left = -1, tran_left = 0).
        Internal interfaces 1..N-1: Impedance mismatch.
        Right substrate (interface N): Matched or free.
        """
        N = self.total_bins
        Z = self.acoustic_impedances
        v = self.velocities

        self.refl_left = np.zeros(N)
        self.refl_right = np.zeros(N)
        self.tran_left = np.zeros(N)
        self.tran_right = np.zeros(N)

        # Interface 0 (Surface z=0, free boundary for bin 0 going left)
        self.refl_left[0] = -1.0
        self.tran_left[0] = 0.0

        # Internal interfaces
        # Interface i is between bin i-1 and bin i
        for i in range(1, N):
            z_prev, z_curr = Z[i - 1], Z[i]
            v_prev, v_curr = v[i - 1], v[i]

            # Waves traveling left from bin i to i-1
            self.refl_left[i] = (z_prev - z_curr) / (z_prev + z_curr)
            self.tran_left[i] = (2.0 * z_curr / (z_prev + z_curr)) * (v_curr / v_prev)

            # Waves traveling right from bin i-1 to i
            self.refl_right[i - 1] = (z_curr - z_prev) / (z_curr + z_prev)
            self.tran_right[i - 1] = (2.0 * z_prev / (z_curr + z_prev)) * (v_prev / v_curr)

        # Substrate boundary at right end of last bin (N-1)
        # Continuation into semi-infinite substrate or rigid boundary:
        # Assuming last bin matched impedance boundary:
        self.refl_right[N - 1] = 0.0
        self.tran_right[N - 1] = 1.0

    def _compute_gsa_coefficients(self) -> None:
        """Compute 6-component GSA coefficient arrays g_left and g_right for each bin."""
        N = self.total_bins
        g0 = np.zeros(N)
        g1 = np.zeros(N)
        g2 = np.zeros(N)

        bin_idx = 0
        for film_idx, (layer, n_bins) in enumerate(zip(self.structure.layers, self.n_bins_per_film)):
            mat = layer.material
            x = mat.attenuation / self.tau_sim
            if x > 0.25:
                raise ValueError(
                    f"Acoustic attenuation in film {film_idx} ({mat.name}) is too large: {x} > 0.25"
                )
            y = mat.sound_velocity * self.tau_sim / self.bin_sizes[film_idx]

            for _ in range(n_bins):
                g0[bin_idx] = x + 0.5 * (y - 2.0) * (y - 1.0)
                g1[bin_idx] = 1.0 - 2.0 * x - (y - 1.0)**2
                g2[bin_idx] = x + 0.5 * y * (y - 1.0)
                bin_idx += 1

        self.g_left = np.zeros((6, N))
        self.g_right = np.zeros((6, N))

        rl = self.refl_left
        rr = self.refl_right
        tl = self.tran_left
        tr = self.tran_right

        for n in range(N):
            # g_left components for bin n
            self.g_left[0, n] = g0[n] + g2[n] * rl[n] * rr[n]
            self.g_left[1, n] = g1[n + 1] * tl[n + 1] if n + 1 < N else 0.0
            self.g_left[2, n] = g2[n + 2] * tl[n + 2] * tl[n + 1] if n + 2 < N else 0.0
            self.g_left[3, n] = g2[n - 1] * tr[n - 1] * rr[n] if n - 1 >= 0 else 0.0
            self.g_left[4, n] = g1[n] * rr[n]
            self.g_left[5, n] = g2[n + 1] * rr[n + 1] * tl[n + 1] if n + 1 < N else 0.0

            # g_right components for bin n
            self.g_right[0, n] = g0[n] + g2[n] * rr[n] * rl[n]
            self.g_right[1, n] = g1[n - 1] * tr[n - 1] if n - 1 >= 0 else 0.0
            self.g_right[2, n] = g2[n - 2] * tr[n - 2] * tr[n - 1] if n - 2 >= 0 else 0.0
            self.g_right[3, n] = g2[n + 1] * tl[n + 1] * rl[n] if n + 1 < N else 0.0
            self.g_right[4, n] = g1[n] * rl[n]
            self.g_right[5, n] = g2[n - 1] * rl[n - 1] * tr[n - 1] if n - 1 >= 0 else 0.0

    def step(self, eta_left: np.ndarray, eta_right: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Advance strain state by one acoustic time step using vectorized GSA stencil.
        """
        N = self.total_bins
        # Pad strain arrays with 2 zeros on left and 2 on right for vectorized stencil offsets
        eL = np.pad(eta_left, (2, 2), mode='constant')
        eR = np.pad(eta_right, (2, 2), mode='constant')

        # eL slice indices relative to padded array:
        # eL[i]: eL[2:N+2]
        # eL[i+1]: eL[3:N+3]
        # eL[i+2]: eL[4:N+4]
        # eR[i-2]: eR[0:N]
        # eR[i-1]: eR[1:N+1]
        # eR[i]: eR[2:N+2]
        # eR[i+1]: eR[3:N+3]

        eta_left_new = (
            eL[2:N+2] * self.g_left[0] +
            eL[3:N+3] * self.g_left[1] +
            eL[4:N+4] * self.g_left[2] +
            eR[1:N+1] * self.g_left[3] +
            eR[2:N+2] * self.g_left[4] +
            eR[3:N+3] * self.g_left[5]
        )

        eta_right_new = (
            eR[2:N+2] * self.g_right[0] +
            eR[1:N+1] * self.g_right[1] +
            eR[0:N]   * self.g_right[2] +
            eL[3:N+3] * self.g_right[3] +
            eL[2:N+2] * self.g_right[4] +
            eL[1:N+1] * self.g_right[5]
        )

        return eta_left_new, eta_right_new
