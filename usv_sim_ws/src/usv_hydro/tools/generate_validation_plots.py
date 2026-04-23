#!/usr/bin/env python3
import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@dataclass
class HydroParams:
    hull_length: float = 2.2
    hull_width: float = 1.1
    hull_height: float = 0.5
    cells_x: int = 14
    cells_y: int = 8
    cells_z: int = 4
    fluid_density: float = 1000.0
    gravity: float = 9.81
    mass: float = 350.0
    water_level: float = 0.0
    cd: tuple[float, float, float] = (0.9, 1.3, 2.0)
    linear_drag: tuple[float, float, float] = (18.0, 55.0, 80.0)
    current_velocity: tuple[float, float, float] = (0.0, 0.0, 0.0)


def make_box_grid(params: HydroParams):
    count = params.cells_x * params.cells_y * params.cells_z
    offsets = np.zeros((count, 3), dtype=np.float64)

    min_x = -0.5 * params.hull_length
    min_y = -0.5 * params.hull_width
    min_z = -0.5 * params.hull_height

    idx = 0
    for ix in range(params.cells_x):
        for iy in range(params.cells_y):
            for iz in range(params.cells_z):
                x = min_x + (ix + 0.5) * params.hull_length / params.cells_x
                y = min_y + (iy + 0.5) * params.hull_width / params.cells_y
                z = min_z + (iz + 0.5) * params.hull_height / params.cells_z
                offsets[idx, :] = np.array([x, y, z], dtype=np.float64)
                idx += 1

    cell_volume = (
        params.hull_length * params.hull_width * params.hull_height / count
    )
    cell_height = params.hull_height / params.cells_z
    drag_area = np.array(
        [
            (params.hull_width * params.hull_height) / count,
            (params.hull_length * params.hull_height) / count,
            (params.hull_length * params.hull_width) / count,
        ],
        dtype=np.float64,
    )
    return offsets, cell_volume, cell_height, drag_area


def rotz(yaw_rad: float):
    c = math.cos(yaw_rad)
    s = math.sin(yaw_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


def aggregate_forces(
    offsets_local: np.ndarray,
    cell_volume: float,
    cell_height: float,
    drag_area_cell: np.ndarray,
    params: HydroParams,
    pos_world: np.ndarray,
    rot_world: np.ndarray,
    linear_vel_world: np.ndarray,
    angular_vel_world: np.ndarray,
    current_velocity_world: np.ndarray,
):
    total_buoyancy = np.zeros(3, dtype=np.float64)
    total_drag = np.zeros(3, dtype=np.float64)
    total_force = np.zeros(3, dtype=np.float64)
    total_moment = np.zeros(3, dtype=np.float64)

    cd = np.array(params.cd, dtype=np.float64)
    linear_drag = np.array(params.linear_drag, dtype=np.float64)

    for offset_local in offsets_local:
        offset_world = rot_world @ offset_local
        cell_pos_world = pos_world + offset_world

        depth = params.water_level - cell_pos_world[2]
        if depth <= 0.0:
            continue

        submergence = max(0.0, min(1.0, depth / cell_height))
        if submergence <= 0.0:
            continue

        displaced_volume = cell_volume * submergence
        buoyancy_world = np.array(
            [0.0, 0.0, params.fluid_density * params.gravity * displaced_volume],
            dtype=np.float64,
        )

        cell_vel_world = linear_vel_world + np.cross(angular_vel_world, offset_world)
        rel_vel_world = cell_vel_world - current_velocity_world
        rel_vel_body = rot_world.T @ rel_vel_world

        drag_body = np.zeros(3, dtype=np.float64)
        for i in range(3):
            v = rel_vel_body[i]
            drag_body[i] = (
                -linear_drag[i] * v
                -0.5
                * params.fluid_density
                * cd[i]
                * drag_area_cell[i]
                * abs(v)
                * v
            )
        drag_world = rot_world @ (drag_body * submergence)
        force = buoyancy_world + drag_world

        total_buoyancy += buoyancy_world
        total_drag += drag_world
        total_force += force
        total_moment += np.cross(offset_world, force)

    return {
        "buoyancy": total_buoyancy,
        "drag": total_drag,
        "force": total_force,
        "moment": total_moment,
    }


def solve_discretized_equilibrium_z(
    offsets_local, cell_volume, cell_height, drag_area_cell, params
):
    low = -0.5
    high = 0.5
    rot = rotz(0.0)
    vel0 = np.zeros(3, dtype=np.float64)

    def net_fz(z):
        agg = aggregate_forces(
            offsets_local,
            cell_volume,
            cell_height,
            drag_area_cell,
            params,
            np.array([0.0, 0.0, z], dtype=np.float64),
            rot,
            vel0,
            vel0,
            np.array(params.current_velocity, dtype=np.float64),
        )
        return agg["force"][2] - params.mass * params.gravity

    if not (net_fz(low) > 0.0 and net_fz(high) < 0.0):
        raise RuntimeError("Could not bracket vertical equilibrium root.")

    for _ in range(80):
        mid = 0.5 * (low + high)
        if net_fz(mid) > 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def run_hydrostatic_sim(
    offsets_local, cell_volume, cell_height, drag_area_cell, params, z0=0.40, dt=0.01, steps=4000
):
    times = np.arange(steps, dtype=np.float64) * dt
    z = z0
    vz = 0.0
    z_hist = np.zeros(steps, dtype=np.float64)
    vz_hist = np.zeros(steps, dtype=np.float64)

    rot = rotz(0.0)
    omega = np.zeros(3, dtype=np.float64)
    current = np.array(params.current_velocity, dtype=np.float64)

    for i in range(steps):
        agg = aggregate_forces(
            offsets_local,
            cell_volume,
            cell_height,
            drag_area_cell,
            params,
            np.array([0.0, 0.0, z], dtype=np.float64),
            rot,
            np.array([0.0, 0.0, vz], dtype=np.float64),
            omega,
            current,
        )
        net_fz = agg["force"][2] - params.mass * params.gravity
        az = net_fz / params.mass
        vz += az * dt
        z += vz * dt
        z_hist[i] = z
        vz_hist[i] = vz

    return times, z_hist, vz_hist


def run_surge_decay(
    offsets_local,
    cell_volume,
    cell_height,
    drag_area_cell,
    params,
    z_eq,
    vx0=2.0,
    dt=0.01,
    steps=2000,
):
    times = np.arange(steps, dtype=np.float64) * dt
    vx_hist = np.zeros(steps, dtype=np.float64)
    fx_hist = np.zeros(steps, dtype=np.float64)

    vx = vx0
    rot = rotz(0.0)
    omega = np.zeros(3, dtype=np.float64)
    current = np.array(params.current_velocity, dtype=np.float64)

    for i in range(steps):
        agg = aggregate_forces(
            offsets_local,
            cell_volume,
            cell_height,
            drag_area_cell,
            params,
            np.array([0.0, 0.0, z_eq], dtype=np.float64),
            rot,
            np.array([vx, 0.0, 0.0], dtype=np.float64),
            omega,
            current,
        )
        ax = agg["force"][0] / params.mass
        vx += ax * dt
        vx_hist[i] = vx
        fx_hist[i] = agg["force"][0]

    return times, vx_hist, fx_hist


def sweep_current_cancellation(
    offsets_local, cell_volume, cell_height, drag_area_cell, params, z_eq
):
    sweep_params = HydroParams(**params.__dict__)
    sweep_params.current_velocity = (1.2, -0.4, 0.0)

    u_values = np.linspace(-1.0, 3.0, 120)
    drag_x = np.zeros_like(u_values)

    rot = rotz(0.0)
    omega = np.zeros(3, dtype=np.float64)
    current = np.array(sweep_params.current_velocity, dtype=np.float64)

    for i, u in enumerate(u_values):
        agg = aggregate_forces(
            offsets_local,
            cell_volume,
            cell_height,
            drag_area_cell,
            sweep_params,
            np.array([0.0, 0.0, z_eq], dtype=np.float64),
            rot,
            np.array([u, current[1], 0.0], dtype=np.float64),
            omega,
            current,
        )
        drag_x[i] = agg["drag"][0]

    return u_values, drag_x, sweep_params.current_velocity[0]


def sweep_symmetry_moments(
    offsets_local, cell_volume, cell_height, drag_area_cell, params, z_eq
):
    yaw_deg = np.linspace(0.0, 360.0, 181)
    mx = np.zeros_like(yaw_deg)
    my = np.zeros_like(yaw_deg)

    vel0 = np.zeros(3, dtype=np.float64)
    current = np.array(params.current_velocity, dtype=np.float64)

    for i, deg in enumerate(yaw_deg):
        yaw = math.radians(deg)
        rot = rotz(yaw)
        agg = aggregate_forces(
            offsets_local,
            cell_volume,
            cell_height,
            drag_area_cell,
            params,
            np.array([0.0, 0.0, z_eq], dtype=np.float64),
            rot,
            vel0,
            vel0,
            current,
        )
        mx[i] = agg["moment"][0]
        my[i] = agg["moment"][1]

    return yaw_deg, mx, my


def save_csv(path: Path, header: list[str], rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Generate hydrodynamics validation plots.")
    parser.add_argument(
        "--out-dir",
        default="usv_sim_ws/src/usv_hydro/validation_output",
        help="Directory for generated PNG/CSV/JSON outputs.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    params = HydroParams()
    offsets_local, cell_volume, cell_height, drag_area_cell = make_box_grid(params)

    z_eq_discrete = solve_discretized_equilibrium_z(
        offsets_local, cell_volume, cell_height, drag_area_cell, params
    )
    z_eq_analytical = (
        0.5 * params.hull_height
        - params.mass / (params.fluid_density * params.hull_length * params.hull_width)
    )

    t_h, z_h, vz_h = run_hydrostatic_sim(
        offsets_local, cell_volume, cell_height, drag_area_cell, params
    )
    t_s, vx_s, fx_s = run_surge_decay(
        offsets_local,
        cell_volume,
        cell_height,
        drag_area_cell,
        params,
        z_eq_discrete,
    )
    u_vals, drag_x, u_current = sweep_current_cancellation(
        offsets_local, cell_volume, cell_height, drag_area_cell, params, z_eq_discrete
    )
    yaw_deg, mx, my = sweep_symmetry_moments(
        offsets_local, cell_volume, cell_height, drag_area_cell, params, z_eq_discrete
    )

    plt.style.use("seaborn-v0_8-whitegrid")

    # Hydrostatic figure
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(t_h, z_h, label="z(t) [m]", color="#005f73", linewidth=2.0)
    ax1.axhline(z_eq_discrete, color="#0a9396", linestyle="--", label="discrete equilibrium z")
    ax1.axhline(z_eq_analytical, color="#ee9b00", linestyle=":", label="analytical equilibrium z")
    ax1.set_title("Hydrostatic Equilibrium Response")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("z [m]")
    ax1.legend(loc="best")
    fig1.tight_layout()
    fig1.savefig(out_dir / "hydrostatic_equilibrium.png", dpi=180)
    plt.close(fig1)

    # Surge decay figure
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    ax2.plot(t_s, vx_s, color="#9b2226", linewidth=2.0, label="surge velocity vx")
    ax2.set_title("Surge Decay Without Thrust")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("vx [m/s]")
    ax2.legend(loc="best")
    fig2.tight_layout()
    fig2.savefig(out_dir / "surge_decay.png", dpi=180)
    plt.close(fig2)

    # Current cancellation figure
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    ax3.plot(u_vals, drag_x, color="#bb3e03", linewidth=2.0, label="drag_x")
    ax3.axvline(u_current, color="#005f73", linestyle="--", label="current-matched u")
    ax3.axhline(0.0, color="black", linewidth=1.0)
    ax3.set_title("Current Cancellation Check")
    ax3.set_xlabel("vehicle surge velocity u [m/s]")
    ax3.set_ylabel("aggregated drag_x [N]")
    ax3.legend(loc="best")
    fig3.tight_layout()
    fig3.savefig(out_dir / "current_cancellation.png", dpi=180)
    plt.close(fig3)

    # Symmetry moments figure
    fig4, ax4 = plt.subplots(figsize=(10, 5))
    ax4.plot(yaw_deg, mx, color="#0a9396", linewidth=2.0, label="Mx (roll moment)")
    ax4.plot(yaw_deg, my, color="#ca6702", linewidth=2.0, label="My (pitch moment)")
    ax4.axhline(0.0, color="black", linewidth=1.0)
    ax4.set_title("Symmetry Check: Roll/Pitch Moments vs Yaw")
    ax4.set_xlabel("yaw [deg]")
    ax4.set_ylabel("moment [N*m]")
    ax4.legend(loc="best")
    fig4.tight_layout()
    fig4.savefig(out_dir / "symmetry_moments.png", dpi=180)
    plt.close(fig4)

    # CSV outputs
    save_csv(
        out_dir / "hydrostatic_equilibrium.csv",
        ["time_s", "z_m", "vz_mps"],
        zip(t_h.tolist(), z_h.tolist(), vz_h.tolist()),
    )
    save_csv(
        out_dir / "surge_decay.csv",
        ["time_s", "vx_mps", "fx_N"],
        zip(t_s.tolist(), vx_s.tolist(), fx_s.tolist()),
    )
    save_csv(
        out_dir / "current_cancellation.csv",
        ["u_mps", "drag_x_N"],
        zip(u_vals.tolist(), drag_x.tolist()),
    )
    save_csv(
        out_dir / "symmetry_moments.csv",
        ["yaw_deg", "mx_Nm", "my_Nm"],
        zip(yaw_deg.tolist(), mx.tolist(), my.tolist()),
    )

    summary = {
        "discrete_equilibrium_z_m": z_eq_discrete,
        "analytical_equilibrium_z_m": z_eq_analytical,
        "hydrostatic_final_z_m": float(z_h[-1]),
        "hydrostatic_final_vz_mps": float(vz_h[-1]),
        "surge_final_vx_mps": float(vx_s[-1]),
        "current_match_u_mps": float(u_current),
        "drag_at_current_match_N": float(np.interp(u_current, u_vals, drag_x)),
        "max_abs_mx_Nm": float(np.max(np.abs(mx))),
        "max_abs_my_Nm": float(np.max(np.abs(my))),
    }
    with (out_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"Generated validation plots and data in: {out_dir}")


if __name__ == "__main__":
    main()
