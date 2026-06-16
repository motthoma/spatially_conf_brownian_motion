# -*- coding: iso-8859-1 -*-
from __future__ import division
import os
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
                style.linestyle.dashed,
                style.linestyle.dashed,
                style.linestyle.solid,
                style.linestyle.solid,
                style.linestyle.solid,
                style.linestyle.solid,
            ]
        ),
        attr.changelist(
            [
                style.linewidth.normal,
                style.linewidth.normal,
                style.linewidth.Thick,
                style.linewidth.Thick,
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
            min=-0.001, max=2.9, painter=p, parter=None, manualticks=xticks, title="$x$"
        ),
        y=graph.axis.lin(
            min=-1.2, max=1.2, painter=p, parter=None, manualticks=yticks, title="$y$"
        ),
        x2=None,
        y2=None,
        key=graph.key.key(pos="tl", dist=0.1),
    )
)


g.plot(
    [
        graph.data.file("confinementc_R_0.05.dat", x=1, y=3, title=None),
        graph.data.file("confinementc_R_0.05.dat", x=1, y=4, title=None),
        graph.data.file("confinementc.dat", x=1, y=5, title=None),
        graph.data.file("confinementc.dat", x=1, y=6, title=None),
        graph.data.file("brownian2d_nonint_modified.dat", x=1, y=2, title=None),
        graph.data.file("brownian2d_nonint_modified.dat", x=3, y=4, title=None),
    ],
    styles=[line],
)


x1, y1 = g.pos(0.6, -0.38)
x2, y2 = g.pos(1.52, 0.118)
x3, y3 = g.pos(1.56, 0.02)
x4, y4 = g.pos(1.005, -0.09)
x5, y5 = g.pos(1.005, 0.09)
x6, y6 = g.pos(1.05, -0.03)
x6, y6 = g.pos(1.05, -0.045)
x7, y7 = g.pos(0.46, 0.19)
x8, y8 = g.pos(0.66, -0.02)
x9, y9 = g.pos(0.50, 0.09)
x10, y10 = g.pos(0.56, 0.05)
# x9, y9 = g.pos(2.45, -0.4)


r = 0.1
c.fill(path.circle(x1, y1, r), [color.rgb.blue, deco.filled([color.rgb.blue])])

c.fill(path.circle(x2, y2, r), [color.rgb.blue, deco.filled([color.rgb.blue])])

c.fill(path.circle(x3, y3, r), [color.rgb.red, deco.filled([color.rgb.red])])

c.fill(path.circle(x9, y9, r), [color.rgb.blue, deco.filled([color.rgb.blue])])

c.fill(path.circle(x10, y10, r), [color.rgb.red, deco.filled([color.rgb.red])])

g.stroke(path.line(x4, y4, x5, y5), [barrow.small, earrow.small])
g.text(x6, y6, r"\small{$2b$}")

c.stroke(path.line(x7, y7, x8, y8), [style.linewidth.Thick])
c.stroke(path.line(x7, y8, x8, y7), [style.linewidth.Thick])

file_name = "confinementc"
c.writePDFfile("../" + file_name)
c.writePDFfile(file_name)
