"""
viscous_drag_estimator.py
--------------------------
Computes viscous resistance coefficients from ship geometry using empirical methods.

Computed quantities:
  - cd[3]              : Morison quadratic drag coefficients [C_dx, C_dy, C_dz]
  - linear_total[3]    : Linear damping coefficients [b_x, b_y, b_z]  [N·s/m]
  - added_mass_approx  : Simple empirical added mass estimate (used when BEM is unavailable)

Methods:
  - C_d  : Hoerner (1965) form drag table + C_B correction
  - b_i  : From ITTC-1957 frictional resistance or critical damping ratio
  - A_ii : Motora (1959) / Lewis section method approximation
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

from .geometry_loader import ShipGeometry


@dataclass
class DragCoefficients:
    """Output of ViscousDragEstimator."""
    # Morison quadratic drag coefficients [-]
    cd_x: float   # Surge — bow/end cross-section
    cd_y: float   # Sway  — side surface
    cd_z: float   # Heave — bottom surface

    # Linear damping coefficients [N·s/m]
    b_x: float    # Surge
    b_y: float    # Sway
    b_z: float    # Heave

    # Empirical added mass estimate [kg]
    a11: float    # Surge
    a22: float    # Sway
    a33: float    # Heave
    a44: float    # Roll  [kg·m²]
    a55: float    # Pitch [kg·m²]
    a66: float    # Yaw   [kg·m²]

    # Source tracking
    method_cd: str = ""
    method_b: str = ""

    @property
    def cd(self) -> List[float]:
        return [self.cd_x, self.cd_y, self.cd_z]

    @property
    def linear_total(self) -> List[float]:
        return [self.b_x, self.b_y, self.b_z]

    def summary(self) -> str:
        lines = [
            "=== DragCoefficients ===",
            f"  Cd  [x,y,z]           : [{self.cd_x:.3f}, {self.cd_y:.3f}, {self.cd_z:.3f}]",
            f"  b   [x,y,z] N·s/m     : [{self.b_x:.1f}, {self.b_y:.1f}, {self.b_z:.1f}]",
            f"  A11 (surge AM)  kg     : {self.a11:.1f}",
            f"  A22 (sway  AM)  kg     : {self.a22:.1f}",
            f"  A33 (heave AM)  kg     : {self.a33:.1f}",
            f"  A44 (roll  AM)  kg·m²  : {self.a44:.2f}",
            f"  A55 (pitch AM)  kg·m²  : {self.a55:.2f}",
            f"  A66 (yaw   AM)  kg·m²  : {self.a66:.2f}",
            f"  Cd method    : {self.method_cd}",
            f"  b  method    : {self.method_b}",
        ]
        return "\n".join(lines)


class ViscousDragEstimator:
    """
    Computes empirical viscous resistance coefficients from ship geometry.

    Usage
    -----
    est = ViscousDragEstimator(geo, zeta_heave=0.15, zeta_roll=0.10)
    coeffs = est.estimate()
    print(coeffs.summary())
    """

    # Kinematic viscosity (fresh water at 15°C)  [m²/s]
    NU_WATER = 1.14e-6

    def __init__(
        self,
        geo: ShipGeometry,
        design_speed: float = 1.5,      # m/s — surge design speed
        zeta_surge: float = 0.05,        # Surge critical damping ratio (very light)
        zeta_sway: float = 0.15,         # Sway
        zeta_heave: float = 0.15,        # Heave
        zeta_roll: float = 0.10,         # Roll (usually kept low)
        fluid_density: float = 1000.0,
        gravity: float = 9.81,
    ):
        self.geo = geo
        self.V_des = design_speed
        self.zeta = {"surge": zeta_surge, "sway": zeta_sway,
                     "heave": zeta_heave, "roll": zeta_roll}
        self.rho = fluid_density
        self.g = gravity

    # ------------------------------------------------------------------
    # Main method
    # ------------------------------------------------------------------
    def estimate(self) -> DragCoefficients:
        cd_x, cd_y, cd_z = self._cd_hoerner()
        a11, a22, a33, a44, a55, a66 = self._added_mass_lewis()
        b_x = self._linear_damping_surge(a11)
        b_y, b_z = self._linear_damping_transverse(a22, a33)

        return DragCoefficients(
            cd_x=cd_x, cd_y=cd_y, cd_z=cd_z,
            b_x=b_x, b_y=b_y, b_z=b_z,
            a11=a11, a22=a22, a33=a33, a44=a44, a55=a55, a66=a66,
            method_cd="Hoerner (1965) + C_B correction",
            method_b="Critical damping ratio + Ogilvie approximation",
        )

    # ------------------------------------------------------------------
    # C_d — Hoerner form drag (Hoerner 1965, Table 3-1)
    # ------------------------------------------------------------------
    def _cd_hoerner(self):
        geo = self.geo
        cb = geo.block_coefficient
        cb = max(0.3, min(0.95, cb))   # safe range

        # Surge: bow section — larger C_B generally implies larger Cd
        # Source: Hoerner Fig.3-24, Molland & Turnock Table 3.2
        cd_x = 0.10 + 1.2 * (cb - 0.50) ** 2   # C_B=0.5 → ~0.10, C_B=0.9 → ~0.29
        cd_x = max(0.07, min(1.5, cd_x))

        # Sway: broad side surface — less sensitive to block coefficient
        # Rectangular-section baseline: ~1.2, correction: ±0.3
        cd_y = 1.10 + 0.5 * (cb - 0.60)
        cd_y = max(0.6, min(2.0, cd_y))

        # Heave: flat bottom plate (normal-to-flow) — typically highest Cd
        # Aspect ratio correction: A_wp = L × B
        ar = geo.length / geo.beam  # length-to-beam ratio
        cd_z = 1.10 + 0.8 / ar     # broad flat plate: ~1.5–2.0
        cd_z = max(0.8, min(2.5, cd_z))

        return cd_x, cd_y, cd_z

    # ------------------------------------------------------------------
    # Added mass — Lewis section method approximation
    # ------------------------------------------------------------------
    def _added_mass_lewis(self):
        geo = self.geo
        rho = self.rho

        # Surge: very small — slender body theory (Munk 1924)
        # A_11 ≈ (C_B - 0.2)·m   (larger for C_B > 0.5)
        a11 = max(0.02, geo.block_coefficient - 0.45) * geo.mass

        # Sway + Heave: Lewis form (semi-cylinder approximation)
        # Per unit length: a'(ω→∞) = ρ · π/2 · (B/2)²  [kg/m]
        a22 = rho * math.pi / 2 * (geo.beam / 2) ** 2 * geo.length * 0.80
        a33 = rho * math.pi / 4 * (geo.beam / 2) ** 2 * geo.length * 0.70

        # Roll: proportional to moment of inertia (Kato 1958)
        a44 = 0.30 * geo.ixx

        # Pitch and Yaw: Munk estimate
        a55 = (a11 - a22) / 2 * geo.length ** 2 / 12  # Munk moment ≈ around zero
        a55 = max(0.10 * geo.iyy, abs(a55))
        a66 = 0.25 * geo.izz

        return a11, a22, a33, a44, a55, a66

    # ------------------------------------------------------------------
    # b_x — Surge linear damping (ITTC-1957 + critical damping)
    # ------------------------------------------------------------------
    def _linear_damping_surge(self, a11: float) -> float:
        geo = self.geo

        # ITTC-1957 frictional resistance
        Re = self.V_des * geo.length / self.NU_WATER
        if Re > 0:
            Cf = 0.075 / (math.log10(Re) - 2.0) ** 2
        else:
            Cf = 0.004  # default

        # Form factor (Holtrop-Mennen) — depends on C_B
        k1 = 0.93 + 0.487 * (geo.block_coefficient / (geo.beam / geo.length)) ** 1.07
        k1 = max(0.1, min(0.5, k1 - 0.93))   # correction term

        # Frictional resistance [N]
        Rf = 0.5 * self.rho * self.V_des ** 2 * geo.wetted_surface_area * Cf * (1 + k1)

        # Linear coefficient [N·s/m] = R_f / V
        b_friction = Rf / self.V_des if self.V_des > 0 else 0.0

        # Cross-check against critical damping (take the larger of the two)
        m_total = geo.mass + a11
        # No natural frequency for surge; assume τ=30s
        tau_surge = 30.0
        b_decay = m_total / tau_surge

        return max(b_friction, b_decay)

    # ------------------------------------------------------------------
    # b_y, b_z — Sway and heave linear damping (critical damping ratio)
    # ------------------------------------------------------------------
    def _linear_damping_transverse(self, a22: float, a33: float):
        geo = self.geo

        # Heave natural frequency
        C33 = self.rho * self.g * geo.waterplane_area
        m_heave = geo.mass + a33
        omega_n_heave = math.sqrt(C33 / m_heave)
        b_z = 2.0 * self.zeta["heave"] * omega_n_heave * m_heave

        # No restoring stiffness for sway → use decay time τ_sway
        # τ_sway ≈ B / (2 · V_des) — characteristic sweep time
        tau_sway = max(5.0, geo.beam / (2.0 * self.V_des + 1e-6))
        m_sway = geo.mass + a22
        b_y = m_sway / tau_sway

        return b_y, b_z
