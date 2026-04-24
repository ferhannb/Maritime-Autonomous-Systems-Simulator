"""
bem_solver.py
-------------
High-level interface wrapping the Capytaine BEM solver.

If Capytaine is not installed, importing the module does not raise an error;
an ImportError is raised only when solve() is called — so workflows that
run without BEM are not affected.

Inputs
------
  - ShipGeometry (box hull or STL mesh)

Outputs
-------
  - BEMResults dataclass
      frequencies  : list[float]   [rad/s]
      added_mass   : ndarray[N,6,6]
      damping      : ndarray[N,6,6]
      exc_re       : ndarray[N,6]
      exc_im       : ndarray[N,6]

Matrices from BEM can be written directly to seakeeping_coeffs.yaml format.

References
----------
  Babarit & Delhommeau (2015), EWTEC — Capytaine code
  Ogilvie (1964) — added mass / damping theory
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .geometry_loader import ShipGeometry


@dataclass
class BEMResults:
    """Output of Capytaine or empirical estimation."""
    frequencies: List[float]           # [rad/s]
    added_mass:  np.ndarray            # shape (N, 6, 6)  [kg] or [kg·m²]
    damping:     np.ndarray            # shape (N, 6, 6)  [N·s/m] or [N·m·s/rad]
    exc_re:      np.ndarray            # shape (N, 6) — wave excitation force real part
    exc_im:      np.ndarray            # shape (N, 6) — wave excitation force imaginary part

    source: str = "unknown"            # "capytaine" or "empirical"

    def __post_init__(self):
        n = len(self.frequencies)
        assert self.added_mass.shape == (n, 6, 6), "added_mass shape mismatch"
        assert self.damping.shape == (n, 6, 6), "damping shape mismatch"
        assert self.exc_re.shape == (n, 6), "exc_re shape mismatch"
        assert self.exc_im.shape == (n, 6), "exc_im shape mismatch"

    def as_dict_list(self) -> list:
        """Returns a list compatible with seakeeping_coeffs.yaml format."""
        result = []
        for i, omega in enumerate(self.frequencies):
            result.append({
                "omega": float(omega),
                "added_mass": self.added_mass[i].tolist(),
                "damping":    self.damping[i].tolist(),
                "excitation_re": self.exc_re[i].tolist(),
                "excitation_im": self.exc_im[i].tolist(),
            })
        return result


class BEMSolver:
    """
    BEM analysis via Capytaine, or empirical estimation when Capytaine is unavailable.

    Usage
    -----
    solver = BEMSolver(geo, omega_min=0.3, omega_max=3.0, n_omega=10)
    results = solver.solve()            # falls back to empirical if Capytaine missing
    results = solver.solve_empirical()  # always empirical
    """

    def __init__(
        self,
        geo: ShipGeometry,
        omega_min: float = 0.3,     # [rad/s]
        omega_max: float = 3.14,    # [rad/s]
        n_omega: int = 10,
        wave_heading: float = 0.0,  # [deg]  0 = head seas
        fluid_density: float = 1000.0,
        gravity: float = 9.81,
        mesh_resolution: float = 0.10,   # [m] — target element size
    ):
        self.geo = geo
        self.omegas = np.linspace(omega_min, omega_max, n_omega).tolist()
        self.wave_heading = wave_heading
        self.rho = fluid_density
        self.g = gravity
        self.mesh_res = mesh_resolution

    # ------------------------------------------------------------------
    # Main solve method — full BEM if Capytaine available, else empirical
    # ------------------------------------------------------------------
    def solve(self, force_empirical: bool = False) -> BEMResults:
        """
        Attempts a BEM solution; returns empirical results if unavailable.
        """
        if force_empirical:
            return self.solve_empirical()

        try:
            import capytaine  # type: ignore  # noqa: F401
        except ImportError:
            print("[BEMSolver] capytaine not installed → using empirical estimation.")
            return self.solve_empirical()

        return self._solve_capytaine()

    # ------------------------------------------------------------------
    # Capytaine BEM
    # ------------------------------------------------------------------
    def _solve_capytaine(self) -> BEMResults:
        """Full panel BEM solution using Capytaine."""
        import capytaine as cpt  # type: ignore

        geo = self.geo

        # ---- Build mesh ----
        if geo.mesh_path is not None:
            mesh = cpt.load_mesh(str(geo.mesh_path), name="hull")
        else:
            # Box hull — Capytaine OffshorePlatform = rectangular prism
            mesh = cpt.mesh_parallelepiped(
                size=(geo.length, geo.beam, geo.draft),
                center=(0.0, 0.0, -geo.draft / 2),
                name="box_hull",
                resolution=(
                    max(2, round(geo.length / self.mesh_res)),
                    max(2, round(geo.beam / self.mesh_res)),
                    max(2, round(geo.draft / self.mesh_res)),
                ),
            )
            mesh.keep_immersed_part()

        body = cpt.FloatingBody(mesh, name="usv_hull")
        body.add_all_rigid_body_dofs()

        # ---- Define problems ----
        problems = []
        for omega in self.omegas:
            for dof in body.dofs:
                problems.append(
                    cpt.RadiationProblem(body=body, radiating_dof=dof, omega=omega,
                                        rho=self.rho, g=self.g)
                )
            # Excitation force
            problems.append(
                cpt.DiffractionProblem(
                    body=body, omega=omega,
                    wave_direction=math.radians(self.wave_heading),
                    rho=self.rho, g=self.g,
                )
            )

        solver = cpt.BEMSolver()
        results = solver.solve_all(problems, keep_details=False)
        ds = cpt.assemble_dataset(results)

        # ---- Convert matrices to numpy ----
        n = len(self.omegas)
        am  = np.zeros((n, 6, 6))
        dam = np.zeros((n, 6, 6))
        fre = np.zeros((n, 6))
        fim = np.zeros((n, 6))

        dof_order = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

        for i_o, omega in enumerate(self.omegas):
            for i_r, dof_r in enumerate(dof_order):
                for i_c, dof_c in enumerate(dof_order):
                    try:
                        val = float(ds["added_mass"].sel(
                            omega=omega, radiating_dof=dof_r, influenced_dof=dof_c,
                            method="nearest"
                        ).values)
                        am[i_o, i_r, i_c] = val
                    except Exception:
                        pass
                    try:
                        val = float(ds["radiation_damping"].sel(
                            omega=omega, radiating_dof=dof_r, influenced_dof=dof_c,
                            method="nearest"
                        ).values)
                        dam[i_o, i_r, i_c] = val
                    except Exception:
                        pass

            for i_d, dof in enumerate(dof_order):
                try:
                    val = complex(ds["Froude_Krylov_force"].sel(
                        omega=omega, influenced_dof=dof, method="nearest"
                    ).values)
                    fre[i_o, i_d] = val.real
                    fim[i_o, i_d] = val.imag
                except Exception:
                    pass

        return BEMResults(
            frequencies=self.omegas,
            added_mass=am, damping=dam,
            exc_re=fre, exc_im=fim,
            source="capytaine",
        )

    # ------------------------------------------------------------------
    # Empirical estimation — without Capytaine
    # ------------------------------------------------------------------
    def solve_empirical(self) -> BEMResults:
        """
        Simple closed-form models for frequency-dependent A(ω), B(ω), F_exc(ω).

        Models
        ------
        A_ii(ω)   : decreases from Lewis value at ω→0 limit to zero at ω→∞
                    A_ii(ω) = A_ii^0 · [1 - tanh((ω-ω_c)/Δω)] / 2
        B_ii(ω)   : proportional to Haskind relation from potential theory
        F_exc(ω)  : simple Froude-Krylov (hydrostatic + inertia)
        """
        from .viscous_drag_estimator import ViscousDragEstimator

        geo = self.geo
        rho, g = self.rho, self.g

        # Zero-frequency added mass estimate (ω→0 limit)
        est = ViscousDragEstimator(geo, fluid_density=rho, gravity=g).estimate()
        a0 = [est.a11, est.a22, est.a33, est.a44, est.a55, est.a66]
        mass_list = [geo.mass, geo.mass, geo.mass, geo.ixx, geo.iyy, geo.izz]

        # Characteristic frequency (in the vicinity of heave natural frequency)
        C33 = rho * g * geo.waterplane_area
        omega_c = math.sqrt(C33 / (geo.mass + a0[2]))

        n = len(self.omegas)
        am  = np.zeros((n, 6, 6))
        dam = np.zeros((n, 6, 6))
        fre = np.zeros((n, 6))
        fim = np.zeros((n, 6))

        dof_names = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

        for i_o, omega in enumerate(self.omegas):
            for k in range(6):
                # A_kk(ω): decreases from a0→0 as ω increases (transition near ω_c)
                width = 0.5 * omega_c
                blend = 0.5 * (1.0 - math.tanh((omega - omega_c) / width))
                am[i_o, k, k] = a0[k] * blend

                # B_kk(ω): peak in resonance region
                # Simple Gaussian: B_peak = ζ_cr · 2·(m+A)·ω_c
                m_eff = mass_list[k] + a0[k]
                b_peak = 0.12 * 2.0 * m_eff * omega_c   # ~ζ=0.12 shift
                dam[i_o, k, k] = b_peak * math.exp(-((omega - omega_c) ** 2) / (2 * width**2))

            # F_exc: simple Froude-Krylov (wave + inertia, small body)
            # Heave excitation: F3 = ρg·A_wp·ζ_a·e^{ikx} — for unit amplitude
            k_wave = omega ** 2 / g
            # Heave (dof=2): |F_3| = ρ·g·A_wp  (unit wave amplitude)
            F3 = rho * g * geo.waterplane_area
            # Phase: for head seas at x=0 → purely real component
            fre[i_o, 2] = F3
            fim[i_o, 2] = 0.0

            # Surge excitation (Froude-Krylov): F1 ≈ ρ·g·k·V_disp
            F1 = rho * g * k_wave * geo.displacement_volume
            fre[i_o, 0] = 0.0
            fim[i_o, 0] = -F1    # 90° phase shift (inertia component)

            # Pitch excitation: M5 ≈ -F3 · L/12 (fore-aft phase difference)
            fre[i_o, 4] = -F3 * geo.length * k_wave / 12.0
            fim[i_o, 4] = 0.0

        return BEMResults(
            frequencies=self.omegas,
            added_mass=am, damping=dam,
            exc_re=fre, exc_im=fim,
            source="empirical",
        )
