"""
geometry_loader.py
------------------
Loads and normalises ship geometry from various sources.

Supported sources:
  - Box parameters (L, B, D, mass) — from SDF/YAML or as arguments
  - STL file (via trimesh)
  - SDF <box> tag (via xml.etree)

Output: ShipGeometry dataclass — used by all modules.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class ShipGeometry:
    """Normalised geometric properties of a ship hull."""

    # ---- Dimensions ----
    length: float       # L  [m]  — ship length (x direction)
    beam: float         # B  [m]  — beam (y direction)
    depth: float        # D  [m]  — total depth (z direction)

    # ---- Mass properties ----
    mass: float         # m  [kg]
    ixx: float = 0.0    # Roll moment of inertia  [kg·m²]
    iyy: float = 0.0    # Pitch moment of inertia [kg·m²]
    izz: float = 0.0    # Yaw moment of inertia   [kg·m²]

    # ---- Waterline properties ----
    draft: float = 0.0  # T  [m] — static draft; 0 → computed automatically

    # ---- Grid suggestion ----
    nx: int = 0         # 0 → computed automatically
    ny: int = 0
    nz: int = 0

    # ---- Mesh (optional, for BEM) ----
    mesh_path: Optional[Path] = field(default=None, repr=False)

    # ---- Derived properties ----
    fluid_density: float = 1000.0

    def __post_init__(self):
        # Automatic draft calculation: equilibrium condition m = ρ·L·B·T
        if self.draft == 0.0:
            self.draft = self.mass / (self.fluid_density * self.length * self.beam)

        # Grid suggestion: at least 4 cells per metre, maximum 20
        if self.nx == 0:
            self.nx = max(4, min(20, round(self.length * 6)))
        if self.ny == 0:
            self.ny = max(4, min(16, round(self.beam * 6)))
        if self.nz == 0:
            self.nz = 4  # Vertical: 4 is generally sufficient

        # Box formula for inertia if not provided
        if self.ixx == 0.0:
            self.ixx = (1.0 / 12.0) * self.mass * (self.beam**2 + self.depth**2)
        if self.iyy == 0.0:
            self.iyy = (1.0 / 12.0) * self.mass * (self.length**2 + self.depth**2)
        if self.izz == 0.0:
            self.izz = (1.0 / 12.0) * self.mass * (self.length**2 + self.beam**2)

    @property
    def waterplane_area(self) -> float:
        """Static waterplane area A_wp = L × B  [m²]"""
        return self.length * self.beam

    @property
    def displacement_volume(self) -> float:
        """Displacement volume V = m / ρ  [m³]"""
        return self.mass / self.fluid_density

    @property
    def wetted_surface_area(self) -> float:
        """Approximate wetted surface area (Denny-Mumford empirical formula)  [m²]"""
        return 1.025 * math.sqrt(self.displacement_volume * self.length)

    @property
    def block_coefficient(self) -> float:
        """Block coefficient C_B = V / (L·B·T)  [-]"""
        return self.displacement_volume / (self.length * self.beam * self.draft)

    @property
    def midship_area(self) -> float:
        """Midship section area (trapezoidal approximation) A_M = B × T  [m²]"""
        return self.beam * self.draft

    @property
    def lateral_plane_area(self) -> float:
        """Lateral plane area A_L = L × T  [m²]"""
        return self.length * self.draft

    def summary(self) -> str:
        lines = [
            "=== ShipGeometry ===",
            f"  L × B × D  : {self.length:.3f} × {self.beam:.3f} × {self.depth:.3f} m",
            f"  Mass       : {self.mass:.1f} kg",
            f"  Draft T    : {self.draft:.4f} m",
            f"  C_B        : {self.block_coefficient:.3f}",
            f"  A_wp       : {self.waterplane_area:.4f} m²",
            f"  V_disp     : {self.displacement_volume:.5f} m³",
            f"  S_wet      : {self.wetted_surface_area:.4f} m²",
            f"  Grid       : {self.nx} × {self.ny} × {self.nz} = {self.nx*self.ny*self.nz} cells",
            f"  I_xx/yy/zz : {self.ixx:.2f} / {self.iyy:.2f} / {self.izz:.2f} kg·m²",
        ]
        return "\n".join(lines)


class GeometryLoader:
    """
    Creates a ShipGeometry object from various sources.

    Usage examples
    --------------
    # 1. From direct parameters
    geo = GeometryLoader.from_box(L=2.2, B=1.1, D=0.5, mass=350)

    # 2. From an SDF file
    geo = GeometryLoader.from_sdf("usv_minimal.sdf")

    # 3. From an STL mesh + parameters
    geo = GeometryLoader.from_stl("hull.stl", mass=350)
    """

    @staticmethod
    def from_box(
        L: float,
        B: float,
        D: float,
        mass: float,
        ixx: float = 0.0,
        iyy: float = 0.0,
        izz: float = 0.0,
        fluid_density: float = 1000.0,
    ) -> ShipGeometry:
        """Create a ShipGeometry from box parameters."""
        return ShipGeometry(
            length=L,
            beam=B,
            depth=D,
            mass=mass,
            ixx=ixx,
            iyy=iyy,
            izz=izz,
            fluid_density=fluid_density,
        )

    @staticmethod
    def from_sdf(sdf_path: str | Path, fluid_density: float = 1000.0) -> ShipGeometry:
        """
        Load geometry from the first <box> + <inertial> tags in an SDF file.
        Supports both Gazebo 1.x and 1.10 formats.
        """
        tree = ET.parse(str(sdf_path))
        root = tree.getroot()

        # ---- Dimensions ----
        size_el = root.find(".//*[box]/box/size")
        if size_el is None:
            raise ValueError(f"SDF file does not contain <box><size> tag: {sdf_path}")
        dims = [float(v) for v in size_el.text.strip().split()]
        if len(dims) != 3:
            raise ValueError(f"<size> must contain 3 values, found: {dims}")
        L, B, D = dims

        # ---- Mass and inertia ----
        def _float(tag: str) -> float:
            el = root.find(f".//*[{tag}]/{tag}")
            if el is None:
                el = root.find(f".//{tag}")
            return float(el.text.strip()) if el is not None else 0.0

        mass = _float("mass")
        ixx  = _float("ixx")
        iyy  = _float("iyy")
        izz  = _float("izz")

        if mass <= 0.0:
            raise ValueError("No valid <mass> value found in SDF file.")

        return ShipGeometry(
            length=L, beam=B, depth=D,
            mass=mass, ixx=ixx, iyy=iyy, izz=izz,
            fluid_density=fluid_density,
        )

    @staticmethod
    def from_stl(
        stl_path: str | Path,
        mass: float,
        ixx: float = 0.0,
        iyy: float = 0.0,
        izz: float = 0.0,
        fluid_density: float = 1000.0,
    ) -> ShipGeometry:
        """
        Extract AABB dimensions from an STL mesh file and create a ShipGeometry.
        Requires the trimesh package: pip install trimesh
        """
        try:
            import trimesh  # type: ignore
        except ImportError:
            raise ImportError(
                "trimesh is required to load STL files: pip install trimesh"
            )

        mesh = trimesh.load(str(stl_path), force="mesh")
        bounds = mesh.bounding_box.extents   # [dx, dy, dz]
        L, B, D = float(bounds[0]), float(bounds[1]), float(bounds[2])

        geo = ShipGeometry(
            length=L, beam=B, depth=D,
            mass=mass, ixx=ixx, iyy=iyy, izz=izz,
            fluid_density=fluid_density,
            mesh_path=Path(stl_path),
        )
        return geo
