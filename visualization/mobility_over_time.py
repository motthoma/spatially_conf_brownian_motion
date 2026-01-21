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

FORCE_POS = 10.0
FORCE_STR_POS = f"{FORCE_POS:.3f}".replace(".", "pt")
DATA_DIR_POS = Path("test_dir_N_1_F_10pt0")
DATA_FILE_POS = DATA_DIR_POS / f"muovert_F_{FORCE_STR_POS}.dat"

FORCE_NEG = -10.0
FORCE_STR_NEG = f"{FORCE_NEG:.3f}".replace(".", "pt")
DATA_DIR_NEG = Path("test_dir_N_1_F_-10pt0")
DATA_FILE_NEG = DATA_DIR_NEG / f"muovert_F_{FORCE_STR_NEG}.dat"
print('DATA_FILE_NEG:', DATA_FILE_NEG)

OUTPUT_FILE = "mobility_over_time"

X_MAX = 1.2
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
            str(DATA_FILE_POS),
            x="abs($1)",
            y=4,
            title="$\mu_{-}(t)$, " + f"$f_+ = {FORCE_POS}$"
        ),
        graph.data.file(
            str(DATA_FILE_NEG),
            x="abs($1)",
            y=4,
            title="$\mu_{+}(t)$, " + f"$f_-={FORCE_NEG}$",
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

# x_text, y_text = g.pos(0.7 * X_MAX, 0.8*Y_MAX)
# g.text(x_text, y_text, rf"$f_+ = {FORCE_POS}$")
# x_text, y_text = g.pos(0.7 * X_MAX, 0.7*Y_MAX)
# g.text(x_text, y_text, rf"$f_-={FORCE_NEG}$")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

g.writePDFfile(OUTPUT_FILE)

