# -*- coding: iso-8859-1 -*-
from __future__ import division
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


data_file_path = "./test_dir_N_1_F_-1pt0/trajectories.dat"
data_file_path = "../runs/test_dir/trajectories.dat"

g.plot(
    [
        graph.data.file(
            data_file_path, x=2, y=3, title=None
        ),
    ],
    styles=[line],
)
file_name = "trajectories"
c.writePDFfile(file_name)
