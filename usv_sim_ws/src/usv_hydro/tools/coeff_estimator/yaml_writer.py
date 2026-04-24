"""
yaml_writer.py
--------------
Writes computed coefficients to the YAML format expected by the usv_hydro plugin.

Generated files:
  1. hydro_profile_<name>.yaml  — cd, linear_total, cell counts, metadata
  2. seakeeping_coeffs_<name>.yaml — BEM A(ω), B(ω), F_exc(ω)

These files can be copied directly to the usv_hydro config folder, or printed
as a profile snippet ready to be appended to hydro_profiles.yaml.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

import yaml

from .geometry_loader import ShipGeometry
from .viscous_drag_estimator import DragCoefficients
from .bem_solver import BEMResults


# ------------------------------------------------------------------
# YAML normalization: Python float / numpy float -> native Python float
# ------------------------------------------------------------------
def _clean(obj):
    """Convert all numpy/np.float64 types to Python float."""
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean(v) for v in obj]
    try:
        return float(obj)
    except (TypeError, ValueError):
        return obj


def _r3(v: float, n: int = 3) -> float:
    """Round and return a Python float."""
    return round(float(v), n)


class YamlWriter:
    """
    Writes computation results to YAML files.

    Usage
    -----
    writer = YamlWriter(output_dir="output/", profile_name="my_usv")
    writer.write_profile(geo, drag)
    writer.write_seakeeping(bem_results)
    """

    def __init__(
        self,
        output_dir: str | Path = ".",
        profile_name: str = "generated",
    ):
        self.out_dir = Path(output_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.name = profile_name

    # ------------------------------------------------------------------
    # 1. Profile YAML — cd, b, cells
    # ------------------------------------------------------------------
    def write_profile(
        self,
        geo: ShipGeometry,
        drag: DragCoefficients,
        seakeeping_file: Optional[str] = None,
    ) -> Path:
        """
        Writes hydro_profile_<name>.yaml.
        Also prints a profile snippet ready to be appended to
        the existing hydro_profiles.yaml.
        """
        profile_data = {
            # ---- Metadata ----
            "_generated_by": "usv_hydro coeff_estimator",
            "_geometry": {
                "L": _r3(geo.length),
                "B": _r3(geo.beam),
                "D": _r3(geo.depth),
                "T": _r3(geo.draft, 4),
                "mass_kg": _r3(geo.mass, 1),
                "C_B": _r3(geo.block_coefficient, 4),
                "drag_method": drag.method_cd,
                "damping_method": drag.method_b,
            },

            # ---- Cell counts ----
            "cells": {
                "x": geo.nx,
                "y": geo.ny,
                "z": geo.nz,
            },

            # ---- Morison quadratic drag ----
            "cd": [_r3(drag.cd_x), _r3(drag.cd_y), _r3(drag.cd_z)],

            # ---- Linear damping ----
            "drag": {
                "linear_total": [
                    _r3(drag.b_x, 2),
                    _r3(drag.b_y, 2),
                    _r3(drag.b_z, 2),
                ],
            },

            # ---- Estimated added mass (informational) ----
            "_added_mass_estimate": {
                "A11_surge_kg":  _r3(drag.a11, 2),
                "A22_sway_kg":   _r3(drag.a22, 2),
                "A33_heave_kg":  _r3(drag.a33, 2),
                "A44_roll_kgm2": _r3(drag.a44, 3),
                "A55_pitch_kgm2":_r3(drag.a55, 3),
                "A66_yaw_kgm2":  _r3(drag.a66, 3),
            },
        }

        # Append seakeeping file reference if provided
        if seakeeping_file:
            profile_data["seakeeping"] = {
                "use_linear_seakeeping": True,
                "coeffs_file": str(seakeeping_file),
            }

        out_path = self.out_dir / f"hydro_profile_{self.name}.yaml"
        with open(out_path, "w") as f:
            yaml.dump(_clean(profile_data), f,
                      default_flow_style=False, sort_keys=False, allow_unicode=True)

        # ---- Print ready-to-use profile snippet ----
        self._print_profile_snippet(geo, drag)

        print(f"[YamlWriter] Profile file written: {out_path}")
        return out_path

    # ------------------------------------------------------------------
    # 2. Seakeeping YAML — A(ω), B(ω), F_exc(ω)
    # ------------------------------------------------------------------
    def write_seakeeping(self, results: BEMResults) -> Path:
        """
        Writes seakeeping_coeffs_<name>.yaml.
        Format is fully compatible with usv_hydro seakeeping_coeffs_example.yaml.
        """
        data = {
            "schema_version": 1,
            "_source": results.source,
            "frequencies": _clean(results.as_dict_list()),
        }

        out_path = self.out_dir / f"seakeeping_coeffs_{self.name}.yaml"
        with open(out_path, "w") as f:
            yaml.dump(data, f,
                      default_flow_style=False, sort_keys=False, allow_unicode=True)

        print(f"[YamlWriter] Seakeeping file written: {out_path}")
        return out_path

    # ------------------------------------------------------------------
    # 3. Snippet to be added to hydro_profiles.yaml
    # ------------------------------------------------------------------
    def _print_profile_snippet(
        self, geo: ShipGeometry, drag: DragCoefficients
    ):
        """
        Prints a ready-to-paste profile block for hydro_profiles.yaml.
        """
        snippet = f"""
# -----------------------------------------------------------------
# Profile: {self.name}
# L={geo.length:.2f}m  B={geo.beam:.2f}m  D={geo.depth:.2f}m  m={geo.mass:.0f}kg
# -----------------------------------------------------------------
{self.name}:
  cells: [{geo.nx}, {geo.ny}, {geo.nz}]
  cd: [{drag.cd_x:.3f}, {drag.cd_y:.3f}, {drag.cd_z:.3f}]
  drag:
    linear_total: [{drag.b_x:.1f}, {drag.b_y:.1f}, {drag.b_z:.1f}]
"""
        print("\n[hydro_profiles.yaml snippet]")
        print(snippet)
