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
import argparse

from confinement_setup import confinement_functions_int, get_conf_type_key


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


def parse_args(default_path: Path):
    parser = argparse.ArgumentParser(description="Visualize Brownian motion trajectories.")
    parser.add_argument("file", nargs="?", default=str(default_path), help="Path to the trajectories data file (default: %(default)s)")
    return parser.parse_args()


def detect_conf_and_get_funcs(data_file_path: Path):
    conf_type_key = get_conf_type_key(data_file_path)
    print(f"Detected confinement type: {conf_type_key}")
    y_eff_func, y_bound_func = confinement_functions_int(conf_type_key)
    if y_eff_func is None:
        print(f"Error: No effective confinement function for {conf_type_key}")
        sys.exit(1)
    return conf_type_key, y_eff_func, y_bound_func


def compute_boundaries(y_eff_func, y_bound_func, n_periods=3, x_res=500):
    """Compute x and y arrays for plotting physical and effective boundaries."""
    x_single = np.linspace(0, 1.0, x_res)
    x_all = np.concatenate([x_single + i * 1.0 for i in range(n_periods)])

    if y_bound_func:
        y_bound_upper = np.array([y_bound_func(x % 1.0) for x in x_all])
    else:
        y_bound_upper = np.zeros_like(x_all)
    y_bound_lower = -y_bound_upper

    y_bound_upper_eff = np.array([y_eff_func(x % 1.0, 0.1) for x in x_all])
    y_bound_lower_eff = -y_bound_upper_eff

    return x_all, y_bound_upper, y_bound_lower, y_bound_upper_eff, y_bound_lower_eff


def create_graph_and_canvas(painter, xticks, yticks):
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
    return c, g


def build_plot_items(data_file_path: Path, x_all, y_upper, y_lower, y_upper_eff, y_lower_eff, y_bound_func):
    plot_items = []

    if y_bound_func:
        plot_items.extend(
            [
                graph.data.points(list(zip(x_all, y_upper)), x=1, y=2, title=None),
                graph.data.points(list(zip(x_all, y_lower)), x=1, y=2, title=None),
            ]
        )

    plot_items.extend(
        [
            graph.data.points(list(zip(x_all, y_upper_eff)), x=1, y=2, title=None),
            graph.data.points(list(zip(x_all, y_lower_eff)), x=1, y=2, title=None),
        ]
    )

    match = re.search(r"F_(-?\d+\.\d{2})", str(data_file_path))
    force = match.group(1) if match else None
    print("force:", force)
    force_string = f"\\small{{${{f={force}}}$}}" if force else None

    x_col = 2
    y_col = 3
    plot_items.append(
        graph.data.file(str(data_file_path), x=x_col, y=y_col, title=force_string if force else None)
    )

    return plot_items


def save_plot(c, g, plot_items, line_style, conf_type_key: str):
    g.plot(plot_items, styles=[line_style])
    print("Saving plot to PDF...")
    output_filename = f"trajectories_{conf_type_key}"
    c.writePDFfile(output_filename)


def main():
    # Data file path (default or from command-line argument)
    default_path = Path("./test_dirs_trajectories/cos_hard_spheres_2026_04_20_23h_19min_40sec_R_0.050_F_15.00_setn_20/trajectories.dat")
    args = parse_args(default_path)
    data_file_path = Path(args.file)

    if not data_file_path.exists():
        print(f"Error: Data file {data_file_path} not found.")
        sys.exit(1)

    conf_type_key, y_eff_func, y_bound_func = detect_conf_and_get_funcs(data_file_path)

    # Configuration and Setup
    line_style, _ = setup_plotting_styles()

    #  Coordinate System Axes
    xticks = [
        graph.axis.tick.tick(0, label=r"$0$"),
        graph.axis.tick.tick(0.1, label=r"$R$"),
        graph.axis.tick.tick(1, label=r"$L$"),
        graph.axis.tick.tick(2, label=r"$2\\,L$"),
    ]

    yticks = [
        graph.axis.tick.tick(1 + 0.1, label=r"$B$"),
        graph.axis.tick.tick(0, label=r"$0$"),
        graph.axis.tick.tick(-1 - 0.1, label=r"$-B$"),
    ]

    painter = graph.axis.painter.regular(basepathattrs=[deco.earrow.normal])

    # Compute boundaries
    x_all, y_upper, y_lower, y_upper_eff, y_lower_eff = compute_boundaries(y_eff_func, y_bound_func)

    # Create plot canvas/graph
    c, g = create_graph_and_canvas(painter, xticks, yticks)

    # Build plot items and save
    plot_items = build_plot_items(data_file_path, x_all, y_upper, y_lower, y_upper_eff, y_lower_eff, y_bound_func)

    save_plot(c, g, plot_items, line_style, conf_type_key)


if __name__ == "__main__":
    main()
