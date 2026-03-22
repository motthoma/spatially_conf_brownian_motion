#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
"""
trajectories.py
----------------
Script to visualize Brownian motion trajectories with confinement boundaries
using the PyX plotting library.
"""

from __future__ import division
import ctypes
import sys
from pathlib import Path
import numpy as np


# PyX setup for LaTeX-quality plots
from pyx import canvas, color, graph, attr, style, deco, text

# Set LaTeX mode and include amsmath for mathematical typesetting
text.set(mode="latex")
text.preamble(r"\usepackage{amsmath}")


def setup_plotting_styles():
    """Defines colors and line styles for plot boundaries and trajectories."""
    boundary_colors = [
        color.rgb.blue,
        color.gray.black,
        color.gray.black,
        color.rgb.red,
        color.rgb.red,
    ]

    boundary_styles = [
        style.linestyle.solid,
        style.linestyle.solid,
        style.linestyle.solid,
        style.linestyle.dashed,
        style.linestyle.dashed,
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


# def get_shared_library():
#     """Loads the shared library."""
#     try:
#         # We use one shared library libconf_splitter.so which contains all implementations
#         lib = ctypes.CDLL("./libconf_splitter.so")
#         return lib
#     except OSError as e:
#         print(f"Error loading libconf_splitter.so: {e}")
#         sys.exit(1)


def confinement_functions_int(key="splitter"):
    """
    Configures and returns the requested confinement functions from the shared library.
    Returns: (yuef_func, yu_func)
    """
    if key == "splitter":
        # Setup y_effective wrapper: double CONF_yuef_wrapper(double x, double y)
        lib = ctypes.CDLL("./libconf_splitter.so")
        lib.CONF_yuef_wrapper.argtypes = [ctypes.c_double, ctypes.c_double]
        lib.CONF_yuef_wrapper.restype = ctypes.c_double

        # Setup y_boundary wrapper: double CONF_yu_splitter(double x)
        lib.CONF_yu_splitter.argtypes = [ctypes.c_double]
        lib.CONF_yu_splitter.restype = ctypes.c_double

        return lib.CONF_yuef_wrapper, lib.CONF_yu_splitter

    elif key == "cos":
        # Setup y_effective cosine wrapper: double CONF_yuef_cos(double x, double y)
        lib = ctypes.CDLL("./libconf_cos.so")
        lib.CONF_yuef_cos.argtypes = [ctypes.c_double, ctypes.c_double]
        lib.CONF_yuef_cos.restype = ctypes.c_double

        # Setup y_effective_boundary wrapper: double CONF_yu_eff_cos(double x)
        if hasattr(lib, "CONF_yu_eff_cos"):
            print(
                "CONF_yu_eff_cos found in libconf_cos.so, setting up effective boundary function..."
            )
            lib.CONF_yu_eff_cos.argtypes = [ctypes.c_double]
            lib.CONF_yu_eff_cos.restype = ctypes.c_double
            # Wrapper to ignore second argument (y) used in main's plot logic
            y_eff_func_for_plot = lambda x, y: lib.CONF_yu_eff_cos(x)
        else:
            # Fallback if no specific effective boundary function is defined for cosine
            # We use a lambda that calls the original check function with a dummy y
            # (though this is not very useful for plotting)
            y_eff_func_for_plot = lambda x, y: lib.CONF_yuef_cos(x, y)

        # Setup y_boundary wrapper: double CONF_yu_cos(double x)
        if hasattr(lib, "CONF_yu_cos"):
            print(
                "CONF_yu_cos found in libconf_cos.so, setting up boundary function..."
            )
            lib.CONF_yu_cos.argtypes = [ctypes.c_double]
            lib.CONF_yu_cos.restype = ctypes.c_double
            y_bound_func = lib.CONF_yu_cos
        else:
            y_bound_func = None

        return y_eff_func_for_plot, y_bound_func

    elif key == "sept":
        # Setup y_effective septated channel wrapper
        lib = ctypes.CDLL("./libconf_sept.so")
        lib.CONF_yuef_sept.argtypes = [ctypes.c_double, ctypes.c_double]
        lib.CONF_yuef_sept.restype = ctypes.c_double

        if hasattr(lib, "CONF_yu_sept"):
            lib.CONF_yu_sept.argtypes = [ctypes.c_double]
            lib.CONF_yu_sept.restype = ctypes.c_double
            return lib.CONF_yuef_sept, lib.CONF_yu_sept
        else:
            return lib.CONF_yuef_sept, None

    else:
        raise ValueError(f"Unknown confinement key: {key}")


def get_conf_type_key(data_file_path):
    """Determines the confinement type key based on existing files or parameters."""
    # 1. Try to find hint in filenames in the data directory
    # Check for conf_*.c files which seem to indicate the type
    for file in data_file_path.parent.glob("conf_*.c"):
        if "splitter" in file.name:
            return "splitter"
        elif "cos" in file.name:
            return "cos"
        elif "sept" in file.name:
            return "sept"

    # 2. Fallback: check parameters_confinement.dat
    params_file = data_file_path.parent / "parameters_confinement.dat"
    if params_file.exists():
        with open(params_file, "r") as f:
            content = f.read().lower()
            if "splitter" in content:
                return "splitter"
            elif "cosine" in content or "cos-shape" in content:
                return "cos"
            elif "septated" in content or "sept" in content:
                return "sept"

    # Default fallback
    return "splitter"


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
        graph.axis.tick.tick(1, label=r"$B$"),
        graph.axis.tick.tick(0, label=r"$0$"),
        graph.axis.tick.tick(-1, label=r"$-B$"),
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
                max=3.9,
                painter=painter,
                parter=None,
                manualticks=xticks,
                title="$x$",
            ),
            y=graph.axis.lin(
                min=-1.2,
                max=1.2,
                painter=painter,
                parter=None,
                manualticks=yticks,
                title="$y$",
            ),
            x2=None,
            y2=None,
            key=graph.key.key(pos="tr", dist=0.1),
        )
    )

    # 5. Plot Data and Boundaries
    plot_items = [graph.data.file(str(data_file_path), x=2, y=3, title=None)]

    if y_bound_func:
        plot_items.extend(
            [
                graph.data.points(
                    list(zip(x_all, y_upper)), x=1, y=2, title=r"$y_\mathrm{u}(x)$"
                ),
                graph.data.points(list(zip(x_all, y_lower)), x=1, y=2, title=None),
            ]
        )

    plot_items.extend(
        [
            graph.data.points(
                list(zip(x_all, y_upper_eff)), x=1, y=2, title=r"$y_\mathrm{ueff}(x)$"
            ),
            graph.data.points(list(zip(x_all, y_lower_eff)), x=1, y=2, title=None),
        ]
    )

    g.plot(plot_items, styles=[line_style])

    # 6. Save Output
    print("Saving plot to PDF...")
    output_filename = f"trajectories_{conf_type_key}"
    c.writePDFfile(output_filename)


if __name__ == "__main__":
    main()
