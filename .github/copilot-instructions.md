# Copilot Instructions for Spatially Confined Brownian Motion

This document provides GitHub Copilot agents with essential context for working on the spatially confined Brownian motion simulator. Consult referenced documentation for detailed guidelines.

## Project Overview

**Spatially Confined Brownian Motion** is a scientific computing project that simulates overdamped Brownian particles in 2D periodic channels using an Euler-Maruyama stochastic integration scheme.

- **Simulation Core**: C (authoritative implementation of physical model)
- **Supporting Tools**: Python (visualization, post-processing, analysis, profiling)
- **Physics**: Overdamped Brownian dynamics with configurable interactions and confinement geometries
- **Build System**: Modular makefiles with compile-time selection of interaction and geometry models

## Critical Repository Rules

### Code Placement
- **C code** (`*.c`, `*.h`, `*.conf`): Only in `src/`
- **Python code**: `python_tools/` (build helpers, profiling) or `visualization/` (plots, animations)
- **Documentation**: `docs/` or `doxygen/`
- **DO NOT commit**: Contents of `runs/` directory (simulation output)

### Architecture Constraints

1. **Simulation Logic is Immutable**: Core physics algorithms live in C and must never be ported to Python. Python exists only to support (visualize, analyze, build) the C simulation.

2. **Modular Physics**: Interaction models (hard-sphere, Lennard-Jones) and confinement geometries (cosine, septagon, splitter) are selected at **compile time** via Makefile and linked modules. Do not hardcode physics choices or move selection logic to runtime.

3. **Numerical Reproducibility**: Preserve exact numerical behavior. Higher-order integration schemes are unnecessary for this system (no inertia, no time-dependent forces). Do not optimize at the expense of reproducibility.

4. **Minimum-Image Convention**: All distance calculations in periodic boundaries must use the minimum-image convention to ensure correct periodic boundary behavior.

## File-Specific Guidelines

Detailed coding standards and workflows are maintained in:

- **`instructions/c.instructions.md`** — C code standards, memory management, error handling, numerics
- **`instructions/python-viz.instructions.md`** — Python coding style, visualization patterns, post-processing patterns
- **`docs/architecture.md`** — Module structure, build system, initialization flow
- **`docs/numerics.md`** — Stochastic integration, boundary conditions, physical model
- **`docs/data-format.md`** — Output file formats, directory structure, HDF5 schema
- **`SKILLS.md`** — Reusable workflows for common tasks (adding physics models, profiling, validation)

## Python Integration Philosophy

Python scripts must support but never duplicate C logic:

- **Visualization**: Generate plots and animations from simulation output (never resimulate in Python)
- **Post-processing**: Extract, transform, and summarize C output for analysis
- **Build Helpers**: Generate makefiles, compile configurations, parameter sweeps
- **Profiling**: Instrument and measure C performance; profile output processing only
- **Analysis**: Statistical analysis of output; hypothesis testing; visualization of results

Use `pathlib.Path`, type annotations, and modern Python idioms. Maximum line length is 80 characters (with `# noqa` for unavoidable longer lines).

## Compiler Requirements

All C code must compile cleanly without warnings:
```
gcc -Wall -Wextra -Wpedantic -Werror
```

Do not introduce compiler-specific extensions unless already present in codebase.

## Performance and Optimization

- Prefer correctness and reproducibility over micro-optimizations
- Measure before optimizing; do not introduce premature optimizations
- Future parallelization: prefer OpenMP and CUDA over MPI
- Do not add new MPI-based features unless explicitly requested

## Future Considerations

Current parallelization uses MPI. Future development should transition to:
1. OpenMP (preferred for shared-memory parallelization)
2. CUDA (for GPU acceleration)

Do not implement new MPI features unless explicitly requested.

## Testing Strategy

Currently, there is no automated test suite. When implementing changes:
- Do not claim tests exist unless they are actually present in the repository
- Document testing approach in `SKILLS.md`
- Manual validation: compare output against known reference simulations
- Numerical regression testing: verify reproducibility across runs

## Key Contacts & Resources

- **Physics Model Details**: See `docs/numerics.md`
- **Architecture & Modules**: See `docs/architecture.md`
- **Output Data Formats**: See `docs/data-format.md`
- **Common Workflows**: See `SKILLS.md`
- **Build System**: See main `README.md` and `docs/architecture.md`

---

**Version**: 1.0  
**Last Updated**: June 2026
