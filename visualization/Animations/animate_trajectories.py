#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt("./test_dir_N_1_F_10pt0/positions.dat")

fig, ax = plt.subplots()

# Channel (fix)
# ax.plot(channel_x, channel_y, 'k-')
ax.set_xlim(-0.01, 1)
ax.set_ylim(-1, 1)

# Teilchen
scat = ax.scatter([], [])

trail_x = []
trail_y = []
line, = ax.plot([], [], lw=1)


def init():
    scat.set_offsets([])
    return scat,


def ring_buffer(buffer, new_value, max_length):
    buffer.append(new_value)
    if len(buffer) > max_length:
        buffer.pop(0)
    return buffer


def update(frame):
    x, y = data[frame]
    trail_x.append(x)
    if len(trail_x) > 20:
        trail_x.pop(0)
    trail_y.append(y)
    if len(trail_y) > 20:
        trail_y.pop(0)
    # trail_x = ring_buffer(trail_x, x, 50)
    # trail_y = ring_buffer(trail_y, y, 50)

    scat.set_offsets([[x, y]])
    line.set_data(trail_x, trail_y)

    return scat, line

ani = FuncAnimation(fig, update, frames=len(data), interval=150)
plt.show()
