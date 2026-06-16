# Data Format & Output Structure

This document specifies the formats of simulation output files, directory structure, and data schemas used in the Spatially Confined Brownian Motion simulator.

## Directory Structure

### Run Output Directory

Each simulation creates a run directory with the following structure:

```
runs/
  experiment_name_YYYY_MM_DD_HHh_MMmin_SSsec_PARAM1_VAL1_PARAM2_VAL2/
  ├── trajectory.h5              (HDF5 file with particle trajectories)
  ├── sim_params.conf            (copy of configuration file used)
  ├── metadata.txt               (human-readable metadata & statistics)
  └── info.json                  (machine-readable metadata)
```

**Example directory name**:
```
runs/exp_hardspheres_2026_06_16_21h_50min_31sec_N_100_R_0.050_F_40.00_setn_3/
```

Format: `<experiment_name>_<timestamp>_<param1>_<val1>_<param2>_<val2>_setn_<set_number>`

### Archive Structure (After Analysis)

Post-processing may create subdirectories:

```
runs/
  experiment_run_YYYY_MM_DD_HHh_MMmin_SSsec_..../
  ├── trajectory.h5
  ├── sim_params.conf
  ├── metadata.txt
  ├── analysis/
  │   ├── mean_squared_displacement.txt
  │   ├── radial_distribution.txt
  │   ├── velocity_distribution.txt
  │   └── particle_flux.txt
  └── plots/
      ├── particles_frame_0.png
      ├── particles_frame_1000.png
      ├── trajectory_particle_0.png
      ├── msd_all_particles.png
      └── animation_trajectory.mp4
```

## Trajectory Output File Format (HDF5)

### File Structure

The primary output is an **HDF5 file** (`trajectory.h5`) containing particle trajectories and metadata.

**Why HDF5?**
- Hierarchical: organizes related datasets
- Efficient: fast I/O for large trajectories
- Portable: platform-independent binary format
- Self-describing: metadata embedded in file
- Python integration: h5py library simplifies access

### HDF5 Schema

```
trajectory.h5
├── trajectory/
│   ├── dataset shape: (n_frames, n_particles, 2)
│   ├── dtype: float64
│   ├── attributes:
│   │   ├── "time_step" (float64): Δt used in simulation
│   │   ├── "n_particles" (int32): number of particles
│   │   ├── "n_frames" (int32): number of frames written
│   │   ├── "channel_width" (float64): W(x=0)
│   │   ├── "channel_length" (float64): period in x
│   │   └── "confinement_type" (string): "cosine" | "splitter" | "sept" | ...
│   └── chunked storage (for efficient access)
│
├── forces/
│   ├── dataset shape: (n_frames, n_particles, 2)
│   ├── dtype: float64
│   ├── attributes: (same as trajectory)
│   └── [optional: only if force tracking enabled]
│
├── metadata/
│   ├── "rng_seed" (int32): random seed used
│   ├── "simulation_time" (float64): elapsed simulation time
│   ├── "wall_clock_time" (float64): elapsed wall-clock time
│   ├── "temperature" (float64): system temperature
│   ├── "n_equilibration_steps" (int32): pre-equilibration steps
│   ├── "n_production_steps" (int32): data collection steps
│   ├── "interaction_model" (string): "lennardjones" | "hardspheres" | ...
│   ├── "rejection_mode" (int32): 0 = resample, 1 = reject
│   └── [other scalar metadata]
│
└── parameters/
    ├── "external_force_x" (float64): applied force in x
    ├── "external_force_y" (float64): applied force in y
    ├── "particle_radius" (float64): σ
    ├── "lennard_jones_epsilon" (float64): ε (if applicable)
    ├── "lennard_jones_sigma_min" (float64): σ_min (if applicable)
    └── [other force field parameters]
```

### Data Access in Python

```python
import h5py
import numpy as np

# Open file
with h5py.File("trajectory.h5", "r") as f:
    # Read trajectory
    traj = f["trajectory"][:]  # Load entire array
    
    # Access attributes
    dt = f["trajectory"].attrs["time_step"]
    n_frames = f["trajectory"].attrs["n_frames"]
    
    # Selective reading (memory efficient)
    frame_0 = f["trajectory"][0, :, :]  # First frame
    particle_0_traj = f["trajectory"][:, 0, :]  # Trajectory of particle 0
    
    # Read metadata
    seed = f["metadata"]["rng_seed"][()]  # Scalar dataset
    
    print(f"Loaded {n_frames} frames of {len(particle_0_traj)} particles")
```

### Storage Details

**Chunking strategy** (for efficient I/O):
- Chunk size: (1, n_particles, 2) — one frame per chunk
- Enables fast frame-by-frame access
- Enables compression (optional GZIP)

**Compression** (optional, configured at write time):
```python
# In output routine (C or Python wrapper)
create_dataset('trajectory', data=traj_data, 
               chunks=(1, n_particles, 2),
               compression='gzip', compression_opts=4)
```

## Configuration File Format (sim_params.conf)

### Text-Based Configuration

Human-readable configuration file using KEY VALUE format.

**Location**: `src/sim_params.conf` (main), or `runs/<run>/sim_params.conf` (copy with output)

### Parameter Specifications

```conf
# ==============================================================================
# SIMULATION PARAMETERS FOR BROWNIAN DYNAMICS
# ==============================================================================
# Each parameter is on a separate line: KEY VALUE
# Comments begin with '#' and extend to end of line
# Blank lines are ignored

# PARTICLE PARAMETERS
N_PARTICLES 100                    # Number of particles
PARTICLE_RADIUS 0.050              # Particle radius σ (in simulation units)

# CHANNEL GEOMETRY (confinement-specific parameters)
CHANNEL_WIDTH 2.0                  # Mean channel width W₀ (for cosine/sept)
CHANNEL_LENGTH 10.0                # Channel period in x-direction
COSINE_AMPLITUDE 0.2               # A in W(x) = W₀ + A·cos(2πx/λ)
COSINE_WAVELENGTH 5.0              # λ in cosine width variation

# TEMPERATURE & DIFFUSION
TEMPERATURE 300.0                  # System temperature (K)
                                    # Sets thermal noise magnitude: √(2Ddt)

# INTERACTION PARAMETERS (model-specific)
INTERACTION_MODEL lennardjones     # "hardspheres" or "lennardjones"
LENNARD_JONES_EPSILON 2.0          # ε (well depth)
LENNARD_JONES_SIGMA_MIN 0.1        # σ_min (characteristic length)

# EXTERNAL FORCES
EXTERNAL_FORCE_X 40.0              # Applied force in x-direction
EXTERNAL_FORCE_Y 0.0               # Applied force in y-direction

# TIME-STEPPING PARAMETERS
TIME_STEP 0.001                    # Δt (integration step size)
N_EQUILIBRATION_STEPS 10000        # Pre-equilibration steps (not recorded)
N_PRODUCTION_STEPS 100000          # Production/data collection steps

# OUTPUT CONFIGURATION
OUTPUT_FILENAME trajectory.h5       # Output HDF5 file path
WRITE_INTERVAL 100                 # Write every N steps (skip steps for efficiency)
WRITE_FORCES 0                      # 0 = no, 1 = write force arrays to HDF5

# COLLISION HANDLING
OVERLAP_MODE 1                      # 1 = reject-and-restore, 0 = resample
OVERLAP_TOLERANCE 1e-12            # Tolerance for overlap detection

# RANDOM NUMBER GENERATION
RNG_SEED 12345                     # Seed for reproducibility
                                    # 0 = use system time (non-reproducible)

# INITIAL CONDITIONS
INITIAL_CONDITION uniform           # "uniform", "random", "boltzmann", "file"
INITIAL_SEED_FILE initial.conf     # If INITIAL_CONDITION = "file"

# PROFILING & DIAGNOSTICS
VERBOSE 1                           # 1 = print progress, 0 = silent
PROFILE_INTERVAL 1000              # Print timing every N steps (if VERBOSE)
```

### Configuration Parsing Rules

- **Whitespace**: Leading/trailing whitespace trimmed
- **Comments**: `#` and everything after ignored
- **Case sensitivity**: Keys are case-insensitive; values may be case-sensitive
- **Type conversion**: 
  - Numbers parsed as int or double as appropriate
  - Strings remain as text
  - Boolean: 0/1 or yes/no
- **Default values**: Missing parameters use safe defaults (see code)
- **Validation**: Invalid ranges trigger errors with helpful messages

### Parameter Dependencies

Some parameters are only used with specific configurations:

```
INTERACTION_MODEL = "lennardjones"
  ├─→ Requires: LENNARD_JONES_EPSILON, LENNARD_JONES_SIGMA_MIN
  ├─→ Ignores: (none specific to hard-sphere)

INTERACTION_MODEL = "hardspheres"
  ├─→ Requires: PARTICLE_RADIUS
  └─→ Ignores: LENNARD_JONES_* parameters

CONFINEMENT_TYPE = "cosine"
  ├─→ Requires: COSINE_AMPLITUDE, COSINE_WAVELENGTH
  └─→ Ignores: SPLITTER_* parameters
```

## Metadata Files

### metadata.txt (Human-Readable)

Text summary of simulation parameters and results.

```
================================================================================
SIMULATION METADATA
================================================================================

Simulation Parameters:
  Number of particles:           100
  Particle radius (σ):           0.0500
  Channel geometry:              cosine
  Channel width (mean):          2.0000
  Channel period (x):            10.0000
  Temperature:                   300.0 K
  
Physics Parameters:
  Interaction model:             lennard-jones
  LJ epsilon (ε):                2.0000
  LJ sigma_min (σ_min):          0.1000
  External force (x):            40.0000
  External force (y):            0.0000

Time-Stepping:
  Time step (Δt):                0.0010
  Equilibration steps:           10000
  Production steps:              100000
  Total simulation time:         100.0000 (in simulation units)
  Total frames written:          1000 (at interval 100)

Output:
  Trajectory file:               trajectory.h5
  RNG seed:                      12345
  Overlap handling:              reject-and-restore
  Collision rejections:          1234 (1.23% of attempts)

Performance:
  Wall-clock time:               234.56 seconds
  Mean time per step:            0.00235 ms
  Particle-steps/sec:            42.7 million

Random Number Generation:
  Generator type:                Mersenne Twister (MT19937)
  Seed:                          12345
  Reproducibility:               YES (same seed → same trajectory)

Quality Assurance:
  Particle density:              50.0 (constant)
  Mean displacement:             12.34 (in simulation units)
  Diffusion coefficient:         0.123
  Mean interaction energy:       -0.456
  
Data Files Generated:
  trajectory.h5                  (1.2 MB)
  sim_params.conf                (0.5 KB)
  metadata.txt                   (this file)
  info.json                      (0.3 KB)
```

### info.json (Machine-Readable)

JSON format for programmatic parsing.

```json
{
  "experiment": {
    "name": "exp_splitter_LJ",
    "timestamp": "2026-06-16T21:50:31Z",
    "git_commit": "a1b2c3d",
    "git_branch": "main",
    "version": "1.0"
  },
  "particles": {
    "n_particles": 100,
    "radius": 0.05
  },
  "geometry": {
    "type": "splitter",
    "width": 2.0,
    "length": 10.0,
    "parameters": {
      "barrier_width": 0.5,
      "splitter_x": 5.0
    }
  },
  "physics": {
    "interaction_model": "lennardjones",
    "interaction_params": {
      "epsilon": 2.0,
      "sigma_min": 0.1
    },
    "external_force": [40.0, 0.0],
    "temperature": 300.0
  },
  "integration": {
    "time_step": 0.001,
    "n_equilibration": 10000,
    "n_production": 100000,
    "total_time": 100.0
  },
  "output": {
    "trajectory_file": "trajectory.h5",
    "write_interval": 100,
    "n_frames": 1000,
    "write_forces": false
  },
  "rng": {
    "seed": 12345,
    "generator": "MT19937"
  },
  "performance": {
    "wall_clock_time_sec": 234.56,
    "particle_steps_per_sec": 42700000,
    "collision_rejections": 1234,
    "rejection_rate": 0.0123
  }
}
```

## Post-Processing Output Formats

### Analysis Results (Text)

Standard format for statistical analysis results:

```
# Mean Squared Displacement Analysis
# Generated: 2026-06-16T21:50:31Z
# Trajectory file: trajectory.h5
# Time step: 0.001
# 
time_lag  msd_mean  msd_std   n_samples
0.1       0.00123   0.00045   1000
0.2       0.00456   0.00089   990
0.5       0.01234   0.00234   950
1.0       0.02456   0.00456   900
2.0       0.04567   0.00890   800
5.0       0.10234   0.01234   500
10.0      0.20456   0.02345   100
```

### Visualization Output

- **PNG images**: `plots/*.png` (static visualizations)
- **MP4 videos**: `animation_*.mp4` (trajectory animations)
- **PDF reports**: `report_*.pdf` (multi-panel analysis summary)

## Conventions & Units

### Dimensionless Units

All simulation quantities are **dimensionless** (in simulation units):

- **Length**: [L₀] where L₀ is a characteristic particle separation
- **Time**: [T₀] where T₀ is a characteristic thermal diffusion time
- **Force**: [F₀] = [E₀/L₀] where E₀ is thermal energy
- **Temperature**: Implicit in diffusion coefficient D

**Conversion to physical units**:

If L₀ = 1 μm and T₀ = 1 ms, then:
```
Simulation distance 1.0 = 1 μm
Simulation time 1.0 = 1 ms
Simulation velocity 1.0 = 1 μm/ms = 1 mm/s
```

### Coordinate System

- **x-axis**: Periodic direction (transport), wraps at x = channel_length
- **y-axis**: Confined direction, 0 ≤ y ≤ W(x)

## Backup & Archival

### Commit Strategy

**Never commit**:
- `runs/` directory (too large, simulation-specific)
- `*.h5` files (binary data)
- Analysis output (regenerable)

**Always commit**:
- Configuration templates (`.conf.template`)
- Python analysis scripts
- Documentation

### Backup Procedure

```bash
# Archive interesting runs
tar czf runs/backup_exp_001_2026_06_16.tar.gz \
  runs/exp_*2026_06_16*

# Store off-site
# (add to external backup system)
```

## References

- **HDF5 Format**: https://www.hdfgroup.org/
- **h5py Documentation**: https://docs.h5py.org/
- **JSON Format**: https://www.json.org/

---

**Version**: 1.0  
**Last Updated**: June 2026
