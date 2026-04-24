#!/usr/bin/env python3
"""
estimate_coeffs.py
------------------
usv_hydro coefficient estimation tool — command-line interface.

Usage examples
--------------
# From box parameters (usv_minimal.sdf values)
python3 estimate_coeffs.py box \
    --L 2.2 --B 1.1 --D 0.5 --mass 350 \
    --ixx 42.58 --iyy 148.46 --izz 176.46 \
    --name usv_minimal --out output/

# Automatic read from SDF file
python3 estimate_coeffs.py sdf \
    --sdf ../usv_description/sdf/usv_minimal.sdf \
    --name usv_minimal --out output/

# STL mesh + BEM (Capytaine required)
python3 estimate_coeffs.py stl \
    --stl hull.stl --mass 350 \
    --bem --omega-min 0.3 --omega-max 3.14 --n-omega 12 \
    --name my_hull --out output/

Output files
------------
  output/hydro_profile_<name>.yaml
  output/seakeeping_coeffs_<name>.yaml  (if --bem flag is set)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add module directory to sys.path (when run directly)
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))

from coeff_estimator.geometry_loader import GeometryLoader
from coeff_estimator.viscous_drag_estimator import ViscousDragEstimator
from coeff_estimator.bem_solver import BEMSolver
from coeff_estimator.yaml_writer import YamlWriter


# ======================================================================
# Sub-commands
# ======================================================================

def run_box(args: argparse.Namespace):
    geo = GeometryLoader.from_box(
        L=args.L, B=args.B, D=args.D, mass=args.mass,
        ixx=args.ixx, iyy=args.iyy, izz=args.izz,
        fluid_density=args.density,
    )
    _process(geo, args)


def run_sdf(args: argparse.Namespace):
    geo = GeometryLoader.from_sdf(args.sdf, fluid_density=args.density)
    _process(geo, args)


def run_stl(args: argparse.Namespace):
    geo = GeometryLoader.from_stl(
        stl_path=args.stl, mass=args.mass,
        ixx=args.ixx, iyy=args.iyy, izz=args.izz,
        fluid_density=args.density,
    )
    _process(geo, args)


def _process(geo, args: argparse.Namespace):
    """Compute coefficients from geometry and write output files."""

    print(geo.summary())
    print()

    # ---- Viscous drag estimation ----
    estimator = ViscousDragEstimator(
        geo,
        design_speed=args.speed,
        zeta_surge=args.zeta_surge,
        zeta_sway=args.zeta_sway,
        zeta_heave=args.zeta_heave,
        fluid_density=args.density,
        gravity=args.gravity,
    )
    drag = estimator.estimate()
    print(drag.summary())
    print()

    writer = YamlWriter(output_dir=args.out, profile_name=args.name)

    # ---- BEM (optional) ----
    bem_results = None
    seakeeping_file = None

    if args.bem:
        solver = BEMSolver(
            geo,
            omega_min=args.omega_min,
            omega_max=args.omega_max,
            n_omega=args.n_omega,
            fluid_density=args.density,
            gravity=args.gravity,
        )
        bem_results = solver.solve(force_empirical=args.empirical_bem)
        out_sk = writer.write_seakeeping(bem_results)
        seakeeping_file = out_sk.name

    writer.write_profile(geo, drag, seakeeping_file=seakeeping_file)

    print("\nDone. Next steps:")
    print(f"  Paste the profile snippet above into hydro_profiles.yaml.")
    if seakeeping_file:
        print(f"  Copy seakeeping_coeffs_{args.name}.yaml to the usv_hydro/config/ folder.")
        print(f"  Update the seakeeping.coeffs_file path under the profile in hydro_profiles.yaml.")


# ======================================================================
# Parser
# ======================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="estimate_coeffs.py",
        description="Estimate usv_hydro coefficients from USV hull geometry",
    )

    sub = p.add_subparsers(dest="command", required=True)

    # ---- Common arguments ----
    def _add_common(parser):
        parser.add_argument("--name",    default="generated", help="Profile name")
        parser.add_argument("--out",     default="output/",   help="Output directory")
        parser.add_argument("--density", type=float, default=1000.0, help="Fluid density [kg/m³]")
        parser.add_argument("--gravity", type=float, default=9.81,   help="Gravitational acceleration [m/s²]")
        parser.add_argument("--speed",   type=float, default=1.5,    help="Design speed [m/s]")
        parser.add_argument("--zeta-surge", type=float, default=0.05, dest="zeta_surge")
        parser.add_argument("--zeta-sway",  type=float, default=0.15, dest="zeta_sway")
        parser.add_argument("--zeta-heave", type=float, default=0.15, dest="zeta_heave")
        parser.add_argument("--bem",     action="store_true",  help="Run BEM seakeeping analysis")
        parser.add_argument("--empirical-bem", action="store_true", dest="empirical_bem",
                            help="Use empirical BEM instead of Capytaine")
        parser.add_argument("--omega-min", type=float, default=0.3,  dest="omega_min", help="[rad/s]")
        parser.add_argument("--omega-max", type=float, default=3.14, dest="omega_max", help="[rad/s]")
        parser.add_argument("--n-omega",   type=int,   default=10,   dest="n_omega")

    def _add_inertia(parser):
        parser.add_argument("--ixx", type=float, default=0.0, help="Roll moment of inertia [kg·m²]")
        parser.add_argument("--iyy", type=float, default=0.0, help="Pitch moment of inertia [kg·m²]")
        parser.add_argument("--izz", type=float, default=0.0, help="Yaw moment of inertia [kg·m²]")

    # ---- box command ----
    box_p = sub.add_parser("box", help="Estimate from box parameters")
    box_p.add_argument("--L",    type=float, required=True, help="Length [m]")
    box_p.add_argument("--B",    type=float, required=True, help="Beam [m]")
    box_p.add_argument("--D",    type=float, required=True, help="Depth [m]")
    box_p.add_argument("--mass", type=float, required=True, help="Mass [kg]")
    _add_inertia(box_p)
    _add_common(box_p)
    box_p.set_defaults(func=run_box)

    # ---- sdf command ----
    sdf_p = sub.add_parser("sdf", help="Estimate from SDF file")
    sdf_p.add_argument("--sdf",  required=True, help="SDF file path")
    _add_common(sdf_p)
    sdf_p.set_defaults(func=run_sdf)

    # ---- stl command ----
    stl_p = sub.add_parser("stl", help="Estimate from STL mesh file")
    stl_p.add_argument("--stl",  required=True, help="STL file path")
    stl_p.add_argument("--mass", type=float, required=True, help="Mass [kg]")
    _add_inertia(stl_p)
    _add_common(stl_p)
    stl_p.set_defaults(func=run_stl)

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
