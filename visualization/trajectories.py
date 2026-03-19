# -*- coding: iso-8859-1 -*-
from __future__ import division
import ctypes
import numpy as np
import sys

sys.path.insert(0, "/ltmp/PyX-0.12.1/")
from pyx import *
from pyx.graph import axis
from pyx.deco import barrow, earrow
import math

text.set(mode="latex")

text.preamble(r"\usepackage{amsmath}")

c = canvas.canvas()

colors = [
    color.rgb.blue,
    color.rgb.blue,
    color.gray.black,
    color.gray.black,
    color.rgb.blue,
    color.rgb.red,
]

line = graph.style.line(
    lineattrs=[
        attr.changelist(colors),
        attr.changelist(
            [
                style.linestyle.solid,
                style.linestyle.solid,
            ]
        ),
        attr.changelist(
            [
                style.linewidth.normal,
                style.linewidth.normal,
            ]
        ),
    ]
)

colorstraj = [
    color.cmyk.RubineRed,
    color.cmyk.Cerulean,
    color.cmyk.Plum,
    color.cmyk.Green,
]

traj = graph.style.line(
    lineattrs=[attr.changelist(colorstraj), attr.changelist([style.linestyle.solid])]
)
xticks = [
    graph.axis.tick.tick(0, label="$0$"),
    graph.axis.tick.tick(0.1, label="$R$"),
    graph.axis.tick.tick(1, label="$L$"),
    graph.axis.tick.tick(2, label="$2\,L$"),
]

yticks = [
    graph.axis.tick.tick(1, label="$B$"),
    graph.axis.tick.tick(0, label="$0$"),
    graph.axis.tick.tick(-1, label="$-B$"),
]


# p = graph.axis.painter.regular(basepathattrs=[deco.barrow.normal],)
p = graph.axis.painter.regular(
    basepathattrs=[deco.earrow.normal],
)


g = c.insert(
    graph.graphxy(
        width=7.5,
        x=graph.axis.lin(
            min=-0.001, max=3.9, painter=p, parter=None, manualticks=xticks, title="$x$"
        ),
        y=graph.axis.lin(
            min=-1.2, max=1.2, painter=p, parter=None, manualticks=yticks, title="$y$"
        ),
        x2=None,
        y2=None,
        key=graph.key.key(pos="tl", dist=0.1),
    )
)

# Load shared library
lib = ctypes.CDLL("./libconf_splitter.so")

# Declare argument and return types
lib.CONF_yuef_wrapper.argtypes = [ctypes.c_double, ctypes.c_double]
lib.CONF_yuef_wrapper.restype = ctypes.c_double

x_single_period = np.linspace(0, 1.0, 500)
n_periods = 3

# Build repeated x with shifts
x_all_periods = np.concatenate([x_single_period + i * 1.0 for i in range(n_periods)])

# Compute y only once, then tile
y_single = np.array([lib.CONF_yuef_wrapper(xi, 0.7) for xi in x_single_period])

y_upper = np.tile(y_single, n_periods)
y_lower = -y_upper

g.plot(
    [
        graph.data.points(
            list(zip(x_all_periods, y_upper)), x=1, y=2, title="Upper Boundary"
        ),
        graph.data.points(
            list(zip(x_all_periods, y_lower)), x=1, y=2, title="Lower Boundary"
        ),
    ],
    styles=[line],
)

data_file_path = "./test_dir_N_1_F_-1pt0/trajectories.dat"
data_file_path = "../runs/test_dir/trajectories.dat"

g.plot(
    [
        graph.data.file(data_file_path, x=2, y=3, title=None),
    ],
    styles=[line],
)
file_name = "trajectories"
c.writePDFfile(file_name)
