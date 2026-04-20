import ctypes

def setup_splitter_functions():
    """Loads the splitter confinement functions from the shared library."""
    lib = ctypes.CDLL("./libconf_splitter.so")

    # Setup y_effective wrapper: double CONF_yuef_wrapper(double x, double y)
    lib.CONF_yuef_wrapper.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.CONF_yuef_wrapper.restype = ctypes.c_double

    # Setup y_boundary wrapper: double CONF_yu_splitter(double x)
    lib.CONF_yu_splitter.argtypes = [ctypes.c_double]
    lib.CONF_yu_splitter.restype = ctypes.c_double

    return lib.CONF_yuef_wrapper, lib.CONF_yu_splitter


def setup_cosine_functions():
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
        print("CONF_yu_cos found in libconf_cos.so, setting up boundary function...")  # noqa
        lib.CONF_yu_cos.argtypes = [ctypes.c_double]
        lib.CONF_yu_cos.restype = ctypes.c_double
        y_bound_func = lib.CONF_yu_cos
    else:
        y_bound_func = None

    return y_eff_func_for_plot, y_bound_func


def setup_septated_functions():
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


def confinement_functions_int(key="splitter"):
    """
    Configures and returns the requested confinement functions from the shared library.
    Returns: (yuef_func, yu_func)
    """
    if key == "splitter":
        CONF_yuef_wrapper, CONF_yu_splitter = setup_splitter_functions()
        return CONF_yuef_wrapper, CONF_yu_splitter

    elif key == "cos":
        y_eff_func_for_plot, y_bound_func = setup_cosine_functions()
        return y_eff_func_for_plot, y_bound_func

    elif key == "sept":
        bottleneck_half_width = 0.1
        CONF_yuef_sept, CONF_yu_sept = setup_septated_functions()

        return (
            CONF_yuef_sept,
            lambda x: CONF_yu_sept(x) if x % 1 > 1e-10 else bottleneck_half_width,
        )

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

