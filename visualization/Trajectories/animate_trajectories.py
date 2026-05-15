#!/usr/bin/env python3
"""
Animate particle trajectories stored in a trajectories.dat file.

This script loads a whitespace-separated trajectories file (numpy loadtxt) where
each row corresponds to a frame and columns are organized as: [t, x1, y1, x2, y2, ...].
It displays an animated scatter plot with trailing paths and the confinement
boundary. The confinement type is inferred from the data file path and used to
plot both the boundary and an effective boundary.

Inputs:
- trajectories file path (default: the last example path used in development)
- optional frames-per-second (fps) for the animation

Outputs:
- An interactive matplotlib window showing the animation.

Usage:
python animate_trajectories.py [PATH_TO_TRAJECTORIES] [--fps FPS]
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from confinement_setup import confinement_functions_int, get_conf_type_key


def parse_args(default_path: Path):
    parser = argparse.ArgumentParser(
        description="Animate particle trajectories from a trajectories.dat file"
    )
    parser.add_argument(
        "file",
        nargs="?",
        default=str(default_path),
        help="Path to trajectories file (trajectories.dat)",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=50.0,
        help="Frames per second for the animation (default: 50)",
    )
    return parser.parse_args()


def load_data(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    return data


def prepare_boundaries(conf_type_key, y_eff_func, y_bound_func):
    # x values for plotting boundaries
    x_single = np.linspace(0, 1.0, 500)
    n_periods = 3
    x_all = np.concatenate([x_single + i * 1.0 for i in range(n_periods)])

    y_bound_upper = np.array([y_bound_func(x % 1.0) for x in x_all])
    y_bound_lower = -y_bound_upper

    y_bound_upper_eff = np.array([y_eff_func(x % 1.0, 0.1) for x in x_all])
    y_bound_lower_eff = -y_bound_upper_eff

    return x_all, y_bound_upper, y_bound_lower, y_bound_upper_eff, y_bound_lower_eff


def main():
    # Default path (keeps the last development example as the default)
    test_data_dir = Path("test_dirs_trajectories")
    default_path = (
        test_data_dir
        / "splitter_LJ_pot_2026_05_12_23h_19min_34sec_R_0.050_LJMIN_0.100_EPS_2.00_F_40.00_setn_1"
        / "trajectories.dat"
    )

    args = parse_args(default_path)
    data_file_path = Path(args.file)

    if not data_file_path.exists():
        print(f"Error: trajectories file not found: {data_file_path}")
        sys.exit(1)

    data = load_data(data_file_path)

    # infer confinement type and get boundary functions
    conf_type_key = get_conf_type_key(data_file_path)
    print(f"Detected confinement type: {conf_type_key}")

    y_eff_func, y_bound_func = confinement_functions_int(conf_type_key)
    if y_eff_func is None:
        print(f"Error: No effective confinement function for {conf_type_key}")
        sys.exit(1)

    # prepare plotting boundaries
    x_all, y_bound_upper, y_bound_lower, y_bound_upper_eff, y_bound_lower_eff = (
        prepare_boundaries(conf_type_key, y_eff_func, y_bound_func)
    )

    fig, ax = plt.subplots()
    ax.set_xlim(-0.1, 2.3)
    ax.set_ylim(-1.2, 1.2)

    ax.plot(x_all, y_bound_upper, "-", label="Boundary", color="black")
    ax.plot(x_all, y_bound_lower, "-", label="", color="black")

    ax.plot(x_all, y_bound_upper_eff, "--", label="Eff. Boundary", color="blue")
    ax.plot(x_all, y_bound_lower_eff, "--", label="", color="blue")

    # infer number of particles from data columns if possible
    if data.ndim != 2 or data.shape[1] < 3 or (data.shape[1] - 1) % 2 != 0:
        print("Error: Unexpected data shape in trajectories file")
        sys.exit(1)

    n_particles = (data.shape[1] - 1) // 2

    # visualization parameters
    trail_length = 60
    particle_rad = 120

    # Scatter for all particles
    scat = ax.scatter([], [], s=particle_rad)

    # Trails: line and buffer for each particle
    lines = [ax.plot([], [], lw=1)[0] for _ in range(n_particles)]
    trail_x = [[] for _ in range(n_particles)]
    trail_y = [[] for _ in range(n_particles)]

    def init():
        scat.set_offsets(np.empty((0, 2)))
        for line in lines:
            line.set_data([], [])
        return [scat, *lines]

    def update(frame):
        # extract particle coordinates
        coords = data[frame, 1:].reshape(n_particles, 2)

        # update scatter
        scat.set_offsets(coords)
        # update trails
        for i in range(n_particles):
            x, y = coords[i]

            trail_x[i].append(x)
            trail_y[i].append(y)

            if len(trail_x[i]) > trail_length:
                trail_x[i].pop(0)
                trail_y[i].pop(0)

            lines[i].set_data(trail_x[i], trail_y[i])

        return [scat, *lines]

    # compute interval from fps
    interval_ms = max(1.0, 1000.0 / float(args.fps))

    ani = FuncAnimation(fig, update, frames=len(data), init_func=init, interval=interval_ms)

    plt.show()


if __name__ == "__main__":
    main()
