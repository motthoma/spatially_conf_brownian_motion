# Python Coding Standards & Guidelines (Visualization & Post-Processing)

This document specifies standards for Python code in `python_tools/` and `visualization/` directories of the Spatially Confined Brownian Motion simulator.

## Golden Rule: Never Replicate Simulation Logic

**Python must never duplicate or reimplement physical simulation logic.** Python's role is to:
- Visualize C simulation output
- Post-process results for analysis
- Generate build configurations and sweep parameters
- Profile C performance
- Analyze statistical properties of output

All fundamental physics computations belong in C. Do not port stochastic integration, interaction models, or confinement geometry to Python.

## Code Style & Formatting

### Line Length

- **Maximum line length**: 80 characters
- **Unavoidable exceptions** (long strings, URLs, file paths): append `# noqa` to the line

```python
# Good: splits naturally at 80 chars
trajectory_file = os.path.join(
    output_dir, f"run_{config_id}_{timestamp}.h5"
)

# Long string: use # noqa
long_label = "Particle Distribution for N=1024, T=300K, Force=40.0"  # noqa
```

### Naming Conventions

- **Functions & variables**: `snake_case_functions()`, `snake_case_variables`
- **Classes**: `PascalCaseClasses`
- **Constants**: `UPPER_CASE_CONSTANTS`
- **Private methods/vars**: prefix with `_underscore`

```python
def compute_mean_displacement(trajectory: np.ndarray) -> float:
    """Compute mean squared displacement from trajectory array."""
    return np.mean(np.sum(trajectory**2, axis=1))

class TrajectoryAnalyzer:
    """Analyze particle trajectories from simulation output."""
    
    _DEFAULT_WINDOW_SIZE = 100  # Private constant
```

### Indentation & Formatting

- Use 4 spaces for indentation (never tabs)
- One import per line (except `from x import a, b` for closely related items)
- Alphabetize imports within groups: stdlib, third-party, local

```python
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from trajectory_analyzer import TrajectoryAnalyzer
from config_parser import parse_simulation_config
```

### Comments

- Use `#` comments only for explaining non-obvious logic
- Write docstrings instead of inline comments for function behavior
- Use type hints and clear naming instead of comments for variable purpose

```python
def apply_minimum_image_convention(
    positions: np.ndarray, 
    channel_width: float
) -> np.ndarray:
    """Apply minimum-image convention to x-coordinates.
    
    Assumes periodic boundary conditions in x-direction only.
    See docs/numerics.md for details.
    
    Parameters
    ----------
    positions : np.ndarray
        Shape (n_particles, 2) array of particle positions
    channel_width : float
        Width of periodic domain in x-direction
    
    Returns
    -------
    np.ndarray
        Adjusted positions with periodic boundaries applied
    """
    adjusted = positions.copy()
    # Adjust x-coordinates to be in [-L/2, L/2]
    adjusted[:, 0] = np.where(
        adjusted[:, 0] > channel_width / 2.0,
        adjusted[:, 0] - channel_width,
        adjusted[:, 0]
    )
    adjusted[:, 0] = np.where(
        adjusted[:, 0] < -channel_width / 2.0,
        adjusted[:, 0] + channel_width,
        adjusted[:, 0]
    )
    return adjusted
```

## Type Annotations

**Use type annotations everywhere.** They catch errors early and make code self-documenting.

### Function Signatures

```python
def read_trajectory_file(
    filename: str,
    start_frame: int = 0,
    end_frame: int | None = None,
) -> np.ndarray:
    """Read particle trajectories from HDF5 file.
    
    Parameters
    ----------
    filename : str
        Path to HDF5 file containing trajectory data
    start_frame : int, optional
        First frame to read (default: 0)
    end_frame : int | None, optional
        Last frame to read (default: None, read all)
    
    Returns
    -------
    np.ndarray
        Shape (n_frames, n_particles, 2) trajectory array
    """
    import h5py
    
    with h5py.File(filename, 'r') as f:
        if end_frame is None:
            end_frame = f['trajectory'].shape[0]
        return f['trajectory'][start_frame:end_frame]
```

### Container Type Hints

```python
from typing import Optional, List, Dict, Tuple

def parse_config_directory(
    config_dir: str,
) -> Dict[str, Dict[str, float | int | str]]:
    """Parse all .conf files in directory.
    
    Returns dictionary mapping config filename to parameter dict.
    """
    configs: Dict[str, Dict[str, float | int | str]] = {}
    return configs

def compute_statistics(
    data: np.ndarray,
    percentiles: List[float] = [25, 50, 75],
) -> Tuple[float, float, np.ndarray]:
    """Compute mean, std, and given percentiles of data.
    
    Returns (mean, std, percentile_values).
    """
    pass
```

## Docstrings

Every function, class, and module shall have a concise docstring. Use NumPy docstring format.

### Module Docstring

```python
"""Trajectory analysis and visualization utilities.

This module provides functions for reading, analyzing, and visualizing
particle trajectories from Brownian dynamics simulations. All trajectory
data must come from the C simulation in src/; no physics simulation logic
is implemented here.

See Also
--------
docs/data-format.md : Description of trajectory HDF5 file format
docs/architecture.md : Overview of simulation pipeline
"""
```

### Function Docstring

```python
def compute_mean_squared_displacement(
    trajectory: np.ndarray,
    frame_skip: int = 1,
) -> np.ndarray:
    """Compute mean squared displacement for each time lag.
    
    The MSD is computed as <(r(t+Δt) - r(t))^2> where the average is
    over all particles and all starting times.
    
    Parameters
    ----------
    trajectory : np.ndarray
        Shape (n_frames, n_particles, 2) array of positions
    frame_skip : int, optional
        Sample every frame_skip frames to reduce computation (default: 1)
    
    Returns
    -------
    np.ndarray
        Shape (max_lag,) array of MSD for each time lag
    """
```

### Class Docstring

```python
class SimulationRun:
    """Container for a single simulation run and its metadata.
    
    Attributes
    ----------
    config : dict
        Simulation parameters from .conf file
    trajectory : np.ndarray
        Shape (n_frames, n_particles, 2) trajectory array
    metadata : dict
        Additional metadata (wall-clock time, random seed, etc.)
    """
    
    def __init__(self, config: dict, trajectory: np.ndarray):
        """Initialize SimulationRun."""
        pass
```

## File and Path Handling

**Prefer `pathlib.Path` over `os.path`.**

```python
from pathlib import Path

# Prefer this:
config_dir = Path("runs") / "exp_001"
output_file = config_dir / "trajectory.h5"
backup_file = output_file.with_suffix(".h5.bak")

# Not this:
import os
config_dir = os.path.join("runs", "exp_001")
output_file = os.path.join(config_dir, "trajectory.h5")
backup_file = os.path.splitext(output_file)[0] + ".h5.bak"

# Safe directory operations:
output_dir = Path("results")
output_dir.mkdir(parents=True, exist_ok=True)

# Safe file reading:
config_file = Path("simulation.conf")
if config_file.exists() and config_file.is_file():
    content = config_file.read_text()
```

## Executable Scripts

Every executable Python script must have a `main()` function and proper entry point.

```python
#!/usr/bin/env python3
"""Generate particle visualization from trajectory file.

Usage:
    python visualize_trajectory.py <trajectory_file> <output.png>
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def visualize_trajectory(
    trajectory_file: str,
    output_file: str,
    frame_idx: int = 0,
) -> None:
    """Load trajectory and create visualization.
    
    Parameters
    ----------
    trajectory_file : str
        Path to HDF5 file with trajectory data
    output_file : str
        Path for output PNG image
    frame_idx : int, optional
        Which frame to visualize (default: 0)
    """
    import h5py
    
    with h5py.File(trajectory_file, 'r') as f:
        frame = f['trajectory'][frame_idx]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(frame[:, 0], frame[:, 1], s=50, alpha=0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Particle Configuration (Frame {frame_idx})")
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    """Entry point for command-line execution."""
    if len(sys.argv) != 3:
        print("Usage: python visualize_trajectory.py <input.h5> <output.png>")
        return 1
    
    trajectory_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        visualize_trajectory(trajectory_file, output_file)
        print(f"Successfully wrote visualization to {output_file}")
        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
```

## Visualization Patterns

### Separating Data Processing from Plotting

Keep data processing logic separate from plotting code. This enables testing and reuse.

```python
# data_processing.py
def compute_radial_distribution(
    trajectory: np.ndarray,
    bins: int = 50,
    max_radius: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute radial distribution function from trajectory.
    
    Returns (bin_centers, histogram).
    """
    pass


# visualization.py
import data_processing
import matplotlib.pyplot as plt

def plot_rdf(trajectory_file: str, output_file: str) -> None:
    """Plot radial distribution function."""
    trajectory = load_trajectory(trajectory_file)
    bin_centers, hist = data_processing.compute_radial_distribution(trajectory)
    
    fig, ax = plt.subplots()
    ax.plot(bin_centers, hist)
    ax.set_xlabel("Distance")
    ax.set_ylabel("RDF")
    fig.savefig(output_file)
    plt.close(fig)
```

### Common Visualization Tasks

**Particle positions (scatter plot):**
```python
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(positions[:, 0], positions[:, 1], s=100, alpha=0.6, c='blue')
ax.set_xlim(-CHANNEL_WIDTH/2, CHANNEL_WIDTH/2)
ax.set_ylim(0, CHANNEL_HEIGHT)
ax.set_aspect('equal')
ax.set_xlabel("x (periodic)")
ax.set_ylabel("y (confined)")
ax.set_title("Particle Configuration")
fig.savefig("particles.png", dpi=150)
```

**Trajectory over time (line plot):**
```python
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# x displacement vs time
axes[0].plot(trajectory[:, 0], label='x(t)', linewidth=0.5)
axes[0].set_ylabel("x Position (periodic)")
axes[0].legend()

# y displacement vs time
axes[1].plot(trajectory[:, 1], label='y(t)', linewidth=0.5)
axes[1].set_ylabel("y Position (confined)")
axes[1].set_xlabel("Time (frames)")
axes[1].legend()

fig.tight_layout()
fig.savefig("trajectory.png", dpi=150)
```

**Animation (for time-dependent data):**
```python
import matplotlib.animation as animation

def animate_trajectory(trajectory: np.ndarray, output_file: str) -> None:
    """Create animation of particle motion."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(-CHANNEL_WIDTH/2, CHANNEL_WIDTH/2)
    ax.set_ylim(0, CHANNEL_HEIGHT)
    ax.set_aspect('equal')
    
    scatter = ax.scatter([], [], s=100, alpha=0.6)
    title = ax.set_title("")
    
    def update(frame):
        positions = trajectory[frame]
        scatter.set_offsets(positions)
        title.set_text(f"Frame {frame}")
        return scatter, title
    
    anim = animation.FuncAnimation(
        fig, update, frames=len(trajectory),
        interval=50, blit=True
    )
    anim.save(output_file, writer='ffmpeg', dpi=100)
    plt.close(fig)
```

## Post-Processing Patterns

### Reading Simulation Output

```python
import h5py
from pathlib import Path

def load_simulation_run(run_dir: str) -> dict:
    """Load complete simulation output from run directory.
    
    Returns dict with trajectory, parameters, and metadata.
    """
    run_path = Path(run_dir)
    
    # Load trajectory
    traj_file = run_path / "trajectory.h5"
    with h5py.File(traj_file, 'r') as f:
        trajectory = f['trajectory'][:]
    
    # Load configuration
    config_file = run_path / "sim_params.conf"
    config = parse_config_file(config_file)
    
    # Load metadata
    metadata_file = run_path / "metadata.txt"
    metadata = parse_metadata_file(metadata_file)
    
    return {
        'trajectory': trajectory,
        'config': config,
        'metadata': metadata,
    }
```

### Statistical Analysis

```python
def analyze_displacement_distribution(
    trajectories: list[np.ndarray],
    dt: float = 0.01,
) -> dict:
    """Analyze displacement statistics across multiple runs.
    
    Returns dict with mean, std, and other statistics.
    """
    displacements = []
    for traj in trajectories:
        # Compute displacement for each frame
        disp = np.diff(traj, axis=0)
        # Distance in each direction
        distances = np.sqrt(np.sum(disp**2, axis=1))
        displacements.extend(distances)
    
    return {
        'mean': np.mean(displacements),
        'std': np.std(displacements),
        'median': np.median(displacements),
    }
```

## Dependencies & Imports

### Required Libraries

The project uses (see `requirements.txt`):
- **numpy**: Numerical array operations
- **matplotlib**: Plotting and visualization
- **scipy**: Statistical functions and scientific computing
- **h5py**: HDF5 file I/O for trajectory data

### Using NumPy

Prefer NumPy operations over explicit loops:

```python
# Good: vectorized operations
distances = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
mean_value = np.mean(data)
valid_particles = data[data[:, 2] > 0.5]

# Avoid: explicit loops
distances = []
for i in range(len(x1)):
    d = np.sqrt((x1[i] - x2[i])**2 + (y1[i] - y2[i])**2)
    distances.append(d)
```

### Lazy Imports for Optional Features

```python
def save_animation(trajectory: np.ndarray, output_file: str) -> None:
    """Save trajectory animation (requires ffmpeg)."""
    import matplotlib.animation as animation  # Import only if used
    
    # ... animation code ...
```

## Avoiding Global Mutable State

Don't use module-level mutable objects that persist across function calls.

```python
# BAD: Global mutable state
_cache = {}  # This persists across calls!

def compute_something(key: str) -> float:
    if key not in _cache:
        _cache[key] = expensive_computation(key)
    return _cache[key]


# GOOD: Pass cache as parameter
def compute_something(
    key: str,
    cache: dict | None = None,
) -> float:
    if cache is None:
        cache = {}
    if key not in cache:
        cache[key] = expensive_computation(key)
    return cache[key]


# Or use a class for stateful computations:
class ComputationCache:
    def __init__(self):
        self._cache = {}
    
    def compute(self, key: str) -> float:
        if key not in self._cache:
            self._cache[key] = expensive_computation(key)
        return self._cache[key]
```

## Profiling & Performance

When optimizing Python code:

1. **Profile first**: Use `cProfile` to identify bottlenecks
2. **Measure performance**: Use `timeit` for timing critical sections
3. **Optimize judiciously**: Only optimize code identified as a bottleneck

```python
import cProfile
import pstats
from pathlib import Path

def profile_analysis_pipeline(trajectory_file: str) -> None:
    """Profile the analysis pipeline."""
    profiler = cProfile.Profile()
    
    profiler.enable()
    # Run analysis code
    trajectory = load_trajectory(trajectory_file)
    stats = compute_statistics(trajectory)
    profiler.disable()
    
    # Print results
    ps = pstats.Stats(profiler)
    ps.sort_stats('cumulative')
    ps.print_stats(20)  # Print top 20 functions
```

## Testing & Validation

Since the C simulation is authoritative:

1. **Validate output**: Check that trajectory data satisfies physical constraints
2. **Cross-check**: Compare Python post-processing results with manual calculations
3. **Edge cases**: Test with single particles, high/low forces, different geometries

```python
def validate_trajectory(trajectory: np.ndarray) -> bool:
    """Validate trajectory for physical plausibility.
    
    Returns True if trajectory passes basic sanity checks.
    """
    # Check shapes
    if len(trajectory.shape) != 3:
        return False
    
    # Check y-coordinates within channel
    if np.any(trajectory[:, :, 1] < 0) or np.any(trajectory[:, :, 1] > CHANNEL_HEIGHT):
        return False
    
    # Check for NaN/inf values
    if np.any(np.isnan(trajectory)) or np.any(np.isinf(trajectory)):
        return False
    
    return True
```

## Common Pitfalls to Avoid

1. **Replicating C physics logic in Python** — Never! All physics in C only.
2. **Global mutable state** — Pass state as parameters or use classes.
3. **Type confusion** — Use type hints to catch errors early.
4. **Hardcoded paths** — Use `pathlib.Path` and configuration.
5. **Missing error handling** — Validate input files and catch I/O errors.
6. **Inefficient loops** — Use NumPy vectorization when possible.

## References

- **Data Format**: See `docs/data-format.md` for trajectory HDF5 schema
- **Architecture**: See `docs/architecture.md` for pipeline overview
- **Numerics**: See `docs/numerics.md` for physical model details
- **Workflows**: See `SKILLS.md` for common post-processing tasks

---

**Version**: 1.0  
**Last Updated**: June 2026
