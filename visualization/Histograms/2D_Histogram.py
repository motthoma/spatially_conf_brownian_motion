#-*- coding: iso-8859-1 -*-
from __future__ import division
import os
import sys
sys.path.insert(0, "/ltmp/PyX-0.12.1/")
from pyx import *
from pyx.graph import axis
import math

text.set(mode="latex")

text.preamble(r"\usepackage{amsmath}")

c=canvas.canvas()

xmin = 0
xmax = 0.99999
ymin = -0.9999
ymax = 0.999999

w=6
h=4

g1 = c.insert(graph.graphxy(width=w, height=h, ypos=0,
                            x=graph.axis.lin(min=xmin, max=xmax, title=r"$x$"),
                            y=graph.axis.lin(min=ymin, max=ymax, title=r"$y$"),
                            key=graph.key.key(pos="tr",dist=0.1, hdist=3, vdist=0.3),))



normfac = 1.0**2
normfac = 1.0

g1.plot([
        graph.data.file("./test_dir_N_1_F_-10pt0/meanpos_Histogram2d_F_-10.000.dat", x="($1+$2)/2", y="($3+$4)/2", color="$5/(11000*1.0)", title=None),
	], 
	[graph.style.density(gradient=color.rgbgradient.Rainbow, keygraph=None)])

 
x, y = g1.pos(0.4,0.7)
# g1.text(x,y, r"$N=5$")
file_name = "2D_Histogram_f_-10pt0"
c.writePDFfile(file_name) 
