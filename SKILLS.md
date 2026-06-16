# Reusable Workflows & Common Tasks (SKILLS)

This document provides step-by-step instructions for Copilot agents to perform common tasks in the Spatially Confined Brownian Motion simulator.

## Skill 1: Adding a New Interaction Model

**Goal**: Implement a new particle-particle interaction potential (e.g., soft-core, dipole-dipole).

**Prerequisites**: Understand the interaction module interface (see `docs/architecture.md`).

### Step 1: Design the Physical Model

1. **Define the potential energy**:
   - Write U(r) in mathematical form
   - Specify parameters (energy scales, length scales)
   - Document asymptotic behavior (r → 0, r → ∞)

2. **Derive the force**:
   - F(r) = -dU/dr (magnitude, radial)
   - Vector form: **F** = -∇U(r)
   - Check limits and singularities

3. **Computational efficiency**:
   - Determine cutoff radius (if any)
   - Estimate computational cost (O(N²) pairs)
   - Consider numerical stability

### Step 2: Create Module Files

```bash
cd src/

# Create header file
touch int_mymodel.h

# Create implementation file
touch int_mymodel.c
```

### Step 3: Implement Header File (int_mymodel.h)

```c
/* int_mymodel.h - My new interaction model */
#ifndef INT_MYMODEL_H
#define INT_MYMODEL_H

/* Model parameters (populate from config) */
typedef struct {
    double param1;  /* Energy scale */
    double param2;  /* Length scale */
    double cutoff;  /* Interaction cutoff */
} MyModelParams;

/* Initialize model with parameters */
int init_mymodel(const MyModelParams *params);

/* Compute force between two particles at distance r */
double compute_interaction_force_magnitude(double r);

/* Compute energy between two particles at distance r */
double compute_interaction_energy(double r);

/* Check if interaction is valid (no singularities, etc.) */
int check_interaction_validity(double r);

#endif
```

### Step 4: Implement Source File (int_mymodel.c)

```c
#include <math.h>
#include <stdio.h>
#include "int_mymodel.h"
#include "code_handling.h"

static MyModelParams g_params = {0};
static int g_initialized = 0;

int init_mymodel(const MyModelParams *params) {
    if (!params) {
        fprintf(stderr, "ERROR: init_mymodel() called with NULL params\n");
        return ERR_INVALID_PARAMS;
    }
    
    /* Validate parameters */
    if (params->param1 <= 0 || params->param2 <= 0) {
        fprintf(stderr, "ERROR: MyModel parameters must be positive\n");
        return ERR_INVALID_PARAMS;
    }
    
    g_params = *params;
    g_initialized = 1;
    return ERR_SUCCESS;
}

/* Example: soft-core interaction
   U(r) = epsilon * (sigma / r)^12 for r > r_cutoff
   U(r) = 0 for r <= r_cutoff
*/
double compute_interaction_energy(double r) {
    if (!g_initialized) {
        fprintf(stderr, "ERROR: MyModel not initialized\n");
        return 0.0;
    }
    
    if (r > g_params.cutoff) {
        return 0.0;  /* Beyond cutoff */
    }
    
    /* Your model: e.g., power-law repulsion */
    double ratio = g_params.param2 / r;
    return g_params.param1 * pow(ratio, 12.0);
}

double compute_interaction_force_magnitude(double r) {
    if (!g_initialized) {
        fprintf(stderr, "ERROR: MyModel not initialized\n");
        return 0.0;
    }
    
    if (r > g_params.cutoff) {
        return 0.0;
    }
    
    /* F = -dU/dr (magnitude, repulsive is positive) */
    double ratio = g_params.param2 / r;
    return 12.0 * g_params.param1 * pow(ratio, 12.0) / r;
}

int check_interaction_validity(double r) {
    /* Check for NaN, inf, singularities */
    if (!isfinite(r) || r <= 0) {
        return ERR_INVALID_PARAMS;
    }
    return ERR_SUCCESS;
}
```

### Step 5: Update Makefile

Edit `src/makefile` to include the new module:

```makefile
# Add to INTERACTION options
INTERACTION = mymodel  # or hardspheres or lennardjones

# Add object file to OBJECTS
OBJECTS = main_brownconf.o ... int_$(INTERACTION).o

# Recompile
make clean
make
```

### Step 6: Update Build Helper

Edit `python_tools/dtool_create_header_makefile.py`:

```python
INTERACTION_MODELS = [
    "hardspheres",
    "lennardjones",
    "mymodel",  # Add here
]
```

### Step 7: Test Against Reference Simulation

1. **Create test configuration** (`test_mymodel.conf`):
   ```conf
   N_PARTICLES 10
   INTERACTION_MODEL mymodel
   PARAM1 1.0
   PARAM2 0.1
   N_PRODUCTION_STEPS 1000
   RNG_SEED 12345
   ```

2. **Run simulation**:
   ```bash
   cd src/
   make INTERACTION=mymodel
   ./main_brownconf ../test_mymodel.conf
   ```

3. **Validate numerically**:
   - Run twice with same seed: verify bit-identical output
   - Compare single-particle behavior (should diffuse normally without interactions)
   - Check energy conservation (if applicable)

4. **Verify forces**:
   - Numerically differentiate U to verify F = -dU/dr
   - Check force symmetry: F_12 = -F_21

5. **Benchmark**:
   - Time per step and verify performance is acceptable
   - Compare against reference interaction models

### Step 8: Document in Comments

Update `docs/architecture.md` and interaction module header to document:
- Physical interpretation of parameters
- Cutoff radius and approximations
- Computational complexity
- Known limitations

---

## Skill 2: Adding a New Confinement Geometry

**Goal**: Implement a new channel geometry (e.g., triangular, exponential width variation).

### Step 1: Define Geometry

1. **Mathematical description**:
   - Channel width as function W(x)
   - Periodicity and domain
   - Wall positions (boundaries)

2. **Example - Exponential Channel**:
   ```
   W(x) = W_0 * exp(A * sin(2*pi*x/lambda))
   where:
     W_0 = mean width
     A = amplitude parameter
     lambda = wavelength
   ```

### Step 2: Create Module Files

```bash
cd src/
touch conf_exponential.h conf_exponential.c
```

### Step 3: Implement Header File

```c
/* conf_exponential.h - Exponential width variation channel */
#ifndef CONF_EXPONENTIAL_H
#define CONF_EXPONENTIAL_H

typedef struct {
    double width_0;      /* Mean channel width */
    double amplitude;    /* Amplitude of variation */
    double wavelength;   /* Wavelength of periodicity */
    double period;       /* Domain period */
} ExponentialChannelParams;

int init_exponential_channel(const ExponentialChannelParams *params);
double get_channel_width(double x);
int is_position_valid(const double *pos);
void compute_wall_force(const double *pos, double *force);
void apply_periodic_bc(double *pos);

#endif
```

### Step 4: Implement Source File

Focus on:
1. Validating parameters (width > 0, all finite)
2. Efficient width calculation (memoize if needed)
3. Correct wall force calculation
4. Proper periodic boundary handling

### Step 5: Update Build System

- Modify `src/makefile` to include new geometry
- Update `python_tools/dtool_create_header_makefile.py`

### Step 6: Validate Geometry

1. **Visual inspection**: Plot W(x) over one period
2. **Particle placement**: Verify particles stay within walls
3. **Reference comparison**: Compare against known simulations
4. **Periodicity**: Verify no spurious effects at x=0, x=period

---

## Skill 3: Profiling Performance

**Goal**: Identify bottlenecks and measure scaling performance of the simulator.

### Step 1: Build with Profiling

```bash
cd src/
make clean
gcc -Wall -Wextra -Wpedantic -Werror -O2 -g -pg \
    -o main_brownconf_profile *.c
```

### Step 2: Run with Profiling

```bash
cd profiling/
../src/main_brownconf_profile params_profiling.conf

# Generates gmon.out
```

### Step 3: Analyze Profile

```bash
# View flat profile
gprof ../src/main_brownconf_profile gmon.out

# Top routines by CPU time
gprof ../src/main_brownconf_profile gmon.out | head -30

# Save to file
gprof ../src/main_brownconf_profile gmon.out > profile_report.txt
```

### Step 4: Identify Bottlenecks

Look for:
- Functions consuming >50% time (focus areas)
- Functions called very frequently (optimization targets)
- Unexpected overhead (synchronization, I/O)

**Common bottlenecks**:
- Force calculation loop (expected)
- Random number generation
- Collision detection
- File I/O (if too frequent)

### Step 5: Profile Scaling

```bash
# Test with varying particle counts
for N in 10 50 100 500 1000; do
    echo "N=$N"
    sed -i "s/N_PARTICLES.*/N_PARTICLES $N/" params_profiling.conf
    ./main_brownconf_profile params_profiling.conf
    gprof main_brownconf_profile gmon.out | head -5
done
```

**Expected scaling**:
- O(N) per step: linear increase in CPU time with N
- If super-linear (O(N²)): collision detection may be unoptimized

### Step 6: Python Performance Profiling

For post-processing scripts:

```python
import cProfile
import pstats

def profile_analysis():
    profiler = cProfile.Profile()
    profiler.enable()
    
    # Your analysis code here
    trajectory = load_trajectory("trajectory.h5")
    stats = compute_statistics(trajectory)
    
    profiler.disable()
    ps = pstats.Stats(profiler)
    ps.sort_stats('cumulative')
    ps.print_stats(20)  # Top 20 functions
```

---

## Skill 4: Validating Simulation Output

**Goal**: Verify that simulation results are physically reasonable and numerically correct.

### Step 1: Reproducibility Test

```bash
# Run same simulation twice
./main_brownconf test.conf
mv trajectory.h5 trajectory_run1.h5

./main_brownconf test.conf
mv trajectory.h5 trajectory_run2.h5

# Compare trajectories
python3 << 'PYTHON'
import h5py
import numpy as np

with h5py.File("trajectory_run1.h5", "r") as f1, \
     h5py.File("trajectory_run2.h5", "r") as f2:
    traj1 = f1["trajectory"][:]
    traj2 = f2["trajectory"][:]
    
    diff = np.abs(traj1 - traj2)
    print(f"Max difference: {np.max(diff):.2e}")
    print(f"Mean difference: {np.mean(diff):.2e}")
    
    # Should be < 1e-15 (floating-point precision)
    if np.max(diff) < 1e-12:
        print("✓ REPRODUCIBILITY OK")
    else:
        print("✗ REPRODUCIBILITY FAILED")
PYTHON
```

### Step 2: Sanity Checks

```python
def validate_trajectory(trajectory, channel_width, channel_height):
    """Perform sanity checks on trajectory."""
    
    # Check shape
    assert len(trajectory.shape) == 3, "Shape mismatch"
    n_frames, n_particles, dims = trajectory.shape
    assert dims == 2, "Expected 2D positions"
    
    # Check y-coordinates within channel
    y_coords = trajectory[:, :, 1]
    assert np.all(y_coords >= 0), "Particles below channel"
    assert np.all(y_coords <= channel_height), "Particles above channel"
    
    # Check for NaN/Inf
    assert np.all(np.isfinite(trajectory)), "NaN or Inf values found"
    
    # Check particle density
    density = n_particles / (channel_width * channel_height)
    print(f"✓ Particle density: {density:.2f}")
    
    # Check mean displacement
    mean_disp = np.mean(np.sum((trajectory[-1] - trajectory[0])**2, axis=1))
    print(f"✓ Mean displacement: {mean_disp:.4f}")
    
    return True
```

### Step 3: Physical Constraints

1. **Particle Conservation**: Same particles throughout run
2. **Energy Bounds**: Total energy reasonable (if conserved)
3. **Density**: Particle density constant
4. **Confinement**: No particles outside channel bounds

### Step 4: Statistical Validation

```python
def compute_diffusion_coefficient(trajectory, dt, skip=1):
    """Compute diffusion coefficient from trajectory."""
    # Mean squared displacement: <x²> = 2*D*t
    
    msd = []
    for lag in range(1, len(trajectory), skip):
        disp = trajectory[lag:] - trajectory[:-lag]
        msd.append(np.mean(np.sum(disp**2, axis=1)))
    
    # Fit linear region
    t_vals = np.arange(len(msd)) * skip * dt
    coeffs = np.polyfit(t_vals, msd, 1)
    D = coeffs[0] / 2.0  # slope / 2 = D
    
    print(f"✓ Diffusion coefficient: {D:.6f}")
    return D
```

### Step 5: Convergence Testing

```bash
# Run with different time steps
for dt in 0.001 0.0005 0.00025; do
    echo "Testing dt=$dt"
    sed -i "s/TIME_STEP.*/TIME_STEP $dt/" test.conf
    ./main_brownconf test.conf
    
    # Extract statistics
    python3 compute_stats.py trajectory.h5 >> convergence.txt
done

# Verify statistics converge as dt → 0
```

---

## Skill 5: Creating Visualization Scripts

**Goal**: Generate high-quality plots and animations of simulation output.

### Step 1: Simple Particle Snapshot

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_particle_snapshot(trajectory_file, frame_idx, output_file):
    """Plot particle configuration at given frame."""
    
    with h5py.File(trajectory_file, 'r') as f:
        frame = f['trajectory'][frame_idx]
        channel_width = f['trajectory'].attrs['channel_width']
        channel_height = f['trajectory'].attrs['channel_height']
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot particles
    ax.scatter(frame[:, 0], frame[:, 1], s=100, alpha=0.6, 
               color='blue', edgecolors='black')
    
    # Plot channel boundaries
    ax.axhline(y=0, color='red', linestyle='--', linewidth=2)
    ax.axhline(y=channel_height, color='red', linestyle='--', linewidth=2)
    
    ax.set_xlim(-channel_width/2, channel_width/2)
    ax.set_ylim(-0.2, channel_height + 0.2)
    ax.set_aspect('equal')
    ax.set_xlabel('x (periodic)', fontsize=12)
    ax.set_ylabel('y (confined)', fontsize=12)
    ax.set_title(f'Particle Configuration (Frame {frame_idx})', fontsize=14)
    
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_file}")
```

### Step 2: Trajectory Time Series

```python
def plot_particle_trajectory(trajectory_file, particle_idx, output_file):
    """Plot single particle trajectory over time."""
    
    with h5py.File(trajectory_file, 'r') as f:
        traj = f['trajectory'][:, particle_idx, :]
        dt = f['trajectory'].attrs['time_step']
    
    time = np.arange(len(traj)) * dt
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # x vs time
    axes[0].plot(time, traj[:, 0], linewidth=0.5, color='blue')
    axes[0].set_ylabel('x position', fontsize=11)
    axes[0].set_title(f'Particle {particle_idx} Trajectory')
    axes[0].grid(True, alpha=0.3)
    
    # y vs time
    axes[1].plot(time, traj[:, 1], linewidth=0.5, color='green')
    axes[1].set_ylabel('y position', fontsize=11)
    axes[1].set_xlabel('Time', fontsize=11)
    axes[1].grid(True, alpha=0.3)
    
    fig.tight_layout()
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_file}")
```

### Step 3: Animation

```python
import matplotlib.animation as animation

def animate_trajectory(trajectory_file, output_file, skip=10):
    """Create animation of particle motion."""
    
    with h5py.File(trajectory_file, 'r') as f:
        traj = f['trajectory'][::skip]
        channel_width = f['trajectory'].attrs['channel_width']
        channel_height = f['trajectory'].attrs['channel_height']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter([], [], s=100, alpha=0.6, color='blue')
    title = ax.set_title("")
    
    # Set limits
    ax.set_xlim(-channel_width/2, channel_width/2)
    ax.set_ylim(-0.2, channel_height + 0.2)
    ax.set_aspect('equal')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    ax.axhline(y=channel_height, color='red', linestyle='--', alpha=0.5)
    
    def update(frame):
        positions = traj[frame]
        scatter.set_offsets(positions)
        title.set_text(f"Frame {frame}")
        return scatter, title
    
    anim = animation.FuncAnimation(
        fig, update, frames=len(traj),
        interval=50, blit=True, repeat=True
    )
    
    anim.save(output_file, writer='ffmpeg', dpi=100)
    plt.close(fig)
    print(f"Saved: {output_file}")
```

### Step 4: Statistical Plots

```python
def plot_mean_squared_displacement(trajectory_file, output_file):
    """Plot mean squared displacement."""
    
    with h5py.File(trajectory_file, 'r') as f:
        traj = f['trajectory'][:]
        dt = f['trajectory'].attrs['time_step']
    
    # Compute MSD for each lag
    max_lag = len(traj) // 2
    msd = []
    for lag in range(1, max_lag):
        disp = traj[lag:] - traj[:-lag]
        msd_lag = np.mean(np.sum(disp**2, axis=1))
        msd.append(msd_lag)
    
    time_lags = np.arange(1, max_lag) * dt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(time_lags, msd, 'o-', linewidth=2, markersize=4)
    ax.set_xlabel('Time lag (log)', fontsize=12)
    ax.set_ylabel('MSD (log)', fontsize=12)
    ax.set_title('Mean Squared Displacement', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_file}")
```

---

## Skill 6: Setting up a Parameter Sweep

**Goal**: Run multiple simulations with varying parameters to explore parameter space.

### Step 1: Create Parameter Template

```conf
# sweep_template.conf - Base configuration for sweep

N_PARTICLES 100
PARTICLE_RADIUS 0.050
CHANNEL_WIDTH 2.0
CHANNEL_LENGTH 10.0

TEMPERATURE 300.0
INTERACTION_MODEL lennardjones
LENNARD_JONES_EPSILON 2.0
LENNARD_JONES_SIGMA_MIN 0.1

TIME_STEP 0.001
N_PRODUCTION_STEPS 50000
OUTPUT_FILENAME trajectory.h5
WRITE_INTERVAL 50

RNG_SEED 12345

# These will be swept:
# EXTERNAL_FORCE_X: 10, 20, 40
# PARTICLE_RADIUS: 0.025, 0.050, 0.100
```

### Step 2: Write Sweep Script

```python
#!/usr/bin/env python3
"""Run parameter sweep of Brownian dynamics simulations."""

import subprocess
import tempfile
from pathlib import Path
from itertools import product

def run_sweep():
    """Execute parameter sweep."""
    
    # Parameter ranges
    forces = [10.0, 20.0, 40.0]
    radii = [0.025, 0.050, 0.100]
    seeds = [12345, 54321, 99999]
    
    base_config = Path("sweep_template.conf").read_text()
    
    for force, radius, seed in product(forces, radii, seeds):
        # Generate unique directory name
        run_dir = Path("runs") / (
            f"sweep_F{force:.1f}_R{radius:.3f}_seed{seed}"
        )
        run_dir.mkdir(parents=True, exist_ok=True)
        
        # Create config file
        config_text = base_config + f"\n"
        config_text += f"EXTERNAL_FORCE_X {force}\n"
        config_text += f"PARTICLE_RADIUS {radius}\n"
        config_text += f"RNG_SEED {seed}\n"
        
        config_file = run_dir / "sim_params.conf"
        config_file.write_text(config_text)
        
        # Run simulation
        print(f"Running: F={force}, R={radius}, seed={seed}")
        result = subprocess.run(
            ["./src/main_brownconf", str(config_file)],
            cwd=run_dir, capture_output=True
        )
        
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr.decode()}")
        else:
            print(f"  OK: {run_dir}")

if __name__ == "__main__":
    run_sweep()
```

### Step 3: Post-Process Sweep Results

```python
#!/usr/bin/env python3
"""Analyze sweep results."""

import h5py
import numpy as np
from pathlib import Path

def analyze_sweep():
    """Extract statistics from all sweep runs."""
    
    results = []
    
    for run_dir in sorted(Path("runs").glob("sweep_*")):
        traj_file = run_dir / "trajectory.h5"
        if not traj_file.exists():
            continue
        
        with h5py.File(traj_file, 'r') as f:
            traj = f['trajectory'][:]
            force = f['metadata']['external_force_x']
            radius = f['metadata']['particle_radius']
        
        # Compute statistics
        mean_disp = np.mean(np.sum((traj[-1] - traj[0])**2, axis=1)**0.5)
        
        results.append({
            'force': force,
            'radius': radius,
            'mean_displacement': mean_disp,
        })
    
    # Print table
    print("Force  Radius  Mean Displacement")
    for r in sorted(results, key=lambda x: (x['force'], x['radius'])):
        print(f"{r['force']:.1f}    {r['radius']:.3f}    {r['mean_displacement']:.4f}")

if __name__ == "__main__":
    analyze_sweep()
```

---

## Skill 7: Debugging & Troubleshooting

**Goal**: Diagnose and fix common issues.

### Issue: Simulation crashes with segmentation fault

1. Compile with debug symbols: `gcc -g ...`
2. Run under debugger: `gdb ./main_brownconf`
3. Set breakpoint: `break main`
4. Run: `run config_file`
5. Get backtrace on crash: `bt`

### Issue: Particles "fly away" or diverge

1. Check overlap handling mode (reject vs. resample)
2. Reduce time step Δt
3. Verify wall force is repulsive
4. Check force directions and signs

### Issue: Output file corrupt or empty

1. Verify output path is writable
2. Check disk space available
3. Ensure HDF5 library is installed
4. Try writing to `/tmp/` to isolate the issue

### Issue: Non-reproducible results

1. Verify RNG seed is set and not randomized
2. Check for uninitialized memory (use `calloc`, not `malloc`)
3. Disable OpenMP/threading if enabled
4. Run on same hardware (floating-point can be platform-dependent)

---

## Skill 8: Future Testing Strategy

**Current State**: No automated test suite.

**Recommended Future Approach**:

### 1. Unit Tests (C)

```c
/* test_conf_utils.c - Test confinement utilities */
#include <assert.h>
#include "conf_*.h"

void test_minimum_image_convention() {
    double result = apply_minimum_image(2.5, 2.0);
    assert(abs(result - (-0.5)) < 1e-12);
}

void test_wall_validity() {
    double pos[2] = {1.0, 0.5};
    assert(is_position_valid(pos) == 1);
    
    double pos_invalid[2] = {1.0, -0.1};
    assert(is_position_valid(pos_invalid) == 0);
}

int main() {
    test_minimum_image_convention();
    test_wall_validity();
    printf("All tests passed!\n");
    return 0;
}
```

Compile: `gcc test_*.c src/*.c -o test_suite`  
Run: `./test_suite`

### 2. Regression Tests (Python)

```python
# regression_test.py - Compare against reference outputs

import h5py
import numpy as np

def test_single_particle_diffusion():
    """Single particle should exhibit free diffusion."""
    # Run with N=1
    # Verify: <x²> ∝ t
    pass

def test_lennard_jones_force():
    """Two particles should repel."""
    # Run with N=2, far apart
    # Verify: distance decreases, then equilibrates
    pass

if __name__ == "__main__":
    test_single_particle_diffusion()
    test_lennard_jones_force()
    print("All regression tests passed!")
```

### 3. Continuous Integration

Use GitHub Actions to automatically run tests on every commit:

```yaml
name: Regression Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Compile tests
        run: make test
      - name: Run C unit tests
        run: ./test_suite
      - name: Run Python regression tests
        run: python3 regression_test.py
```

---

## References

- **Architecture Details**: See `docs/architecture.md`
- **Coding Standards**: See `instructions/c.instructions.md` and `instructions/python-viz.instructions.md`
- **Numerical Methods**: See `docs/numerics.md`
- **Data Formats**: See `docs/data-format.md`

---

**Version**: 1.0  
**Last Updated**: June 2026
