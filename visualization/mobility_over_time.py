#!/usr/bin/env python3
"""
Plot mobility over time from Brownian dynamics simulation output.
"""

from pathlib import Path
import sys

# ---------------------------------------------------------------------------
# PyX setup (pinned legacy version)
# ---------------------------------------------------------------------------

sys.path.insert(0, "/ltmp/PyX-0.12.1/")

from pyx import color, graph, text
from pyx.graph import axis


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

FORCE = 10.0
FORCE_STR = f"{FORCE:.3f}".replace(".", "pt")

DATA_DIR = Path("test_dir")
DATA_FILE = DATA_DIR / f"muovert_F_{FORCE_STR}.dat"

OUTPUT_FILE = "mobility_over_time"
X_MAX = 0.6
Y_MAX = 1.5 


# ---------------------------------------------------------------------------
# LaTeX configuration
# ---------------------------------------------------------------------------

# NOTE: PyX 0.12 uses "mode"; newer versions use "engine"
text.set(mode="latex")
text.preamble(r"\usepackage{amsmath}")


# ---------------------------------------------------------------------------
# Plot construction
# ---------------------------------------------------------------------------

g = graph.graphxy(
    width=7,
    height=4.5,
    x=axis.lin(min=0, max=X_MAX, title=r"$t$"),
    y=axis.lin(min=0, max=Y_MAX),
    key=graph.key.key(
        pos="tl",
        dist=0.15,
        vdist=0.5,
        symbolwidth=0.1,
    ),
)

symbol_attrs = [
    color.gradient.BlueRed,
]

g.plot(
    [
        graph.data.file(
            str(DATA_FILE),
            x="abs($1)",
            y=4,
            title=r"$\mu(t)$",
        ),
    ],
    [
        graph.style.symbol(
            graph.style.symbol.changecross,
            size=0.1,
            symbolattrs=symbol_attrs,
        )
    ],
)

# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

x_text, y_text = g.pos(0.75 * X_MAX, 0.8*Y_MAX)
g.text(x_text, y_text, rf"$f_0 = {FORCE}$")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

g.writePDFfile(OUTPUT_FILE)

