#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


from confinement_setup import confinement_functions_int, get_conf_type_key

data_file_path = Path(
    "./sept_hard_spheres_2026_03_25_21h_38min_10sec_R_0.050_F_40.00_setn_1/trajectories.dat"
)  # noqa
data_file_path = Path(
    "./cos_hard_spheres_2026_04_20_23h_19min_40sec_R_0.050_F_15.00_setn_20/trajectories.dat"
)  # noqa
data_file_path = Path(
    "./splitter_LJ_pot_2026_04_20_23h_35min_53sec_R_0.050_LJMIN_0.100_EPS_2.00_F_40.00_setn_20/trajectories.dat"
)  # noqa
data = np.loadtxt(data_file_path)
# data = np.loadtxt("../../runs/test_dir/trajectories.dat")


conf_type_key = get_conf_type_key(data_file_path)
print(f"Detected confinement type: {conf_type_key}")

# Get the functions for the detected confinement
y_eff_func, y_bound_func = confinement_functions_int(conf_type_key)

if y_eff_func is None:
    print(f"Error: No effective confinement function for {conf_type_key}")
    sys.exit(1)
# x values for plotting boundaries
x_single = np.linspace(0, 1.0, 500)
n_periods = 3
x_all = np.concatenate([x_single + i * 1.0 for i in range(n_periods)])

y_bound_upper = np.array([y_bound_func(x % 1.0) for x in x_all])
y_bound_lower = -y_bound_upper

y_bound_upper_eff = np.array([y_eff_func(x % 1.0, 0.1) for x in x_all])
y_bound_lower_eff = -y_bound_upper_eff

fig, ax = plt.subplots()

ax.set_xlim(-0.1, 1.8)
ax.set_ylim(-1.2, 1.2)

ax.plot(x_all, y_bound_upper, "-", label="Boundary", color="black")
ax.plot(x_all, y_bound_lower, "-", label="", color="black")

ax.plot(x_all, y_bound_upper_eff, "--", label="Eff. Boundary", color="blue")
ax.plot(x_all, y_bound_lower_eff, "--", label="", color="blue")

# number of particles and trail length
n_particles = 20
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


ani = FuncAnimation(fig, update, frames=len(data), init_func=init, interval=100)

plt.show()
