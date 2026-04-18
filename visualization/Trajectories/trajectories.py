#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
"""
trajectories.py
----------------
Script to visualize Brownian motion trajectories with confinement boundaries
using the PyX plotting library.
"""

from __future__ import division
import sys
from pathlib import Path
import numpy as np
import re

from confinement_setup import confinement_functions_int, \
                              get_conf_type_key


# PyX setup for LaTeX-quality plots
from pyx import canvas, color, graph, attr, style, deco, text, path

# Set LaTeX mode and include amsmath for mathematical typesetting
text.set(mode="latex")
text.preamble(r"\usepackage{amsmath}")


def setup_plotting_styles():
    """Defines colors and line styles for plot boundaries and trajectories."""
    boundary_colors = [
        color.gray.black,
        color.gray.black,
        color.rgb.red,
        color.rgb.red,
        color.rgb.blue,
    ]

    boundary_styles = [
        style.linestyle.solid,
        style.linestyle.solid,
        style.linestyle.dashed,
        style.linestyle.dashed,
        style.linestyle.solid,
    ]

    boundary_widths = [style.linewidth.normal] * 5

    line_style = graph.style.line(
        lineattrs=[
            attr.changelist(boundary_colors),
            attr.changelist(boundary_styles),
            attr.changelist(boundary_widths),
        ]
    )

    trajectory_colors = [
        color.cmyk.RubineRed,
        color.cmyk.Cerulean,
        color.cmyk.Plum,
        color.cmyk.Green,
    ]

    trajectory_style = graph.style.line(
        lineattrs=[
            attr.changelist(trajectory_colors),
            attr.changelist([style.linestyle.solid]),
        ]
    )

    return line_style, trajectory_style


def main():
    # Data file path (default or from command-line argument)
    if len(sys.argv) > 1:
        data_file_path = Path(sys.argv[1])
    else:
        data_file_path = Path("../../runs/test_dir/trajectories.dat")

    if not data_file_path.exists():
        print(f"Error: Data file {data_file_path} not found.")
        sys.exit(1)

    conf_type_key = get_conf_type_key(data_file_path)
    print(f"Detected confinement type: {conf_type_key}")

    #  Configuration and Setup
    line_style, _ = setup_plotting_styles()
    # lib = get_shared_library()

    # Get the functions for the detected confinement
    y_eff_func, y_bound_func = confinement_functions_int(conf_type_key)

    if y_eff_func is None:
        print(f"Error: No effective confinement function for {conf_type_key}")
        sys.exit(1)

    #  Coordinate System Axes
    xticks = [
        graph.axis.tick.tick(0, label=r"$0$"),
        graph.axis.tick.tick(0.1, label=r"$R$"),
        graph.axis.tick.tick(1, label=r"$L$"),
        graph.axis.tick.tick(2, label=r"$2\,L$"),
    ]

    yticks = [
        graph.axis.tick.tick(1 + 0.1, label=r"$B$"),
        graph.axis.tick.tick(0, label=r"$0$"),
        graph.axis.tick.tick(-1 - 0.1, label=r"$-B$"),
    ]

    painter = graph.axis.painter.regular(basepathattrs=[deco.earrow.normal])

    # 3. Compute Boundaries
    x_single = np.linspace(0, 1.0, 500)
    n_periods = 3
    x_all = np.concatenate([x_single + i * 1.0 for i in range(n_periods)])

    # Effective boundaries (using y=0.7 as in original script)
    y_eff_single = np.array([y_eff_func(xi, 0.1) for xi in x_single])
    y_upper_eff = np.tile(y_eff_single, n_periods)
    y_lower_eff = -y_upper_eff

    # Physical boundaries
    if y_bound_func:
        print("y_bound_func is available, computing physical boundaries...")
        y_bound_single = np.array([y_bound_func(xi) for xi in x_single])
        y_upper = np.tile(y_bound_single, n_periods)
        y_lower = -y_upper
    else:
        # Fallback if no boundary function is available
        y_upper = np.zeros_like(x_all)
        y_lower = np.zeros_like(x_all)

    # 4. Create Plot
    c = canvas.canvas()
    g = c.insert(
        graph.graphxy(
            width=7.5,
            x=graph.axis.lin(
                min=-0.001,
                max=3.0,
                painter=painter,
                parter=None,
                manualticks=xticks,
                title="$x$",
            ),
            y=graph.axis.lin(
                min=-1.3,
                max=1.3,
                painter=painter,
                parter=None,
                manualticks=yticks,
                title="$y$",
            ),
            x2=None,
            y2=None,
            key=graph.key.key(pos="tr", dist=0.1, hdist=-0.1),
        )
    )

    # 5. Plot Data and Boundaries
    plot_items = []

    if y_bound_func:
        plot_items.extend(
            [
                graph.data.points(
                    list(zip(x_all, y_upper)),
                    x=1,
                    y=2,
                    title=None,  # r"\small{$y_\mathrm{conf}(x)$}",
                ),
                graph.data.points(list(zip(x_all, y_lower)), x=1, y=2, title=None),
            ]
        )

    plot_items.extend(
        [
            graph.data.points(
                list(zip(x_all, y_upper_eff)),
                x=1,
                y=2,
                title=None,  # r"\small{$y_\mathrm{conf, eff}(x)$}",
            ),
            graph.data.points(list(zip(x_all, y_lower_eff)), x=1, y=2, title=None),
        ]
    )

    match = re.search(r"F_(-?\d+\.\d{2})", str(data_file_path))
    force = None
    if match:
        force = match.group(1)
    print("force:", force)
    force_string = f"\\small{{$f={force}$}}"
    force = None
    x_col = 2
    y_col = 3
    plot_items.extend(
        [
            graph.data.file(
                str(data_file_path),
                x=x_col,
                y=y_col,
                title=force_string if force else None,
            ),
        ]
    )

    g.plot(plot_items, styles=[line_style])

    # 6. Save Output
    print("Saving plot to PDF...")
    output_filename = f"trajectories_{conf_type_key}"
    c.writePDFfile(output_filename)


if __name__ == "__main__":
    main()
