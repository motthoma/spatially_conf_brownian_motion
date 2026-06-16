# Architecture & Module Structure

This document describes the high-level architecture, module organization, build system, and initialization flow of the Spatially Confined Brownian Motion simulator.

## High-Level Overview

```
┌─────────────────────────────────────────────────────────────┐
│                  Simulation Pipeline                        │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Configuration (sim_params.conf)                            │
│        │                                                     │
│        ▼                                                     │
│  ┌─────────────────────────────────────────┐               │
│  │  C Simulation (main_brownconf)           │               │
│  │  ├─ Load config parameters               │               │
│  │  ├─ Initialize particles & channel       │               │
│  │  ├─ Time-stepping loop                   │               │
│  │  │   ├─ Euler-Maruyama integration       │               │
│  │  │   ├─ Interaction force calculation    │               │
│  │  │   ├─ Confinement checking             │               │
│  │  │   └─ Write trajectory output          │               │
│  │  └─ Cleanup & finalization               │               │
│  └─────────────────────────────────────────┘               │
│        │                                                     │
│        ▼                                                     │
│  Output: trajectory.h5, metadata.txt                        │
│                                                              │
│  Python Post-Processing & Visualization                    │
│  ├─ Read trajectory data                                    │
│  ├─ Compute statistics                                      │
│  ├─ Generate plots & animations                            │
│  └─ Create analysis reports                                │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

## C Module Structure (src/)

The C codebase is organized into **core modules** (always present) and **pluggable modules** (selected at compile time).

### Core Modules

These modules are always compiled and linked:

#### `main_brownconf.c / main_brownconf.h`
**Main simulation driver.**

- Parses command-line arguments and configuration file
- Initializes simulation state
- Implements time-stepping loop
- Handles I/O for configuration and trajectory output
- Coordinates physics calculations via callbacks to pluggable modules

Key functions:
```c
int main(int argc, char *argv[]);  // Entry point
int read_configuration(const char *config_file, SimParams *params);
int initialize_simulation(const SimParams *params, SimulationState *state);
int perform_time_step(SimulationState *state, const SimParams *params);
int write_trajectory_frame(FILE *output, const SimulationState *state);
```

**Numerical scheme**:
- Implements Euler-Maruyama integration
- Updates positions according to: `x_{n+1} = x_n + D∇U*dt + √(2Ddt)*N(0,1)`
- Handles rejection/resampling of invalid moves (configurable)

#### `random_numb_gen.c / random_numb_gen.h`
**Pseudorandom number generation for stochastic dynamics.**

- Provides standard normal and uniform random numbers
- Uses MT19937 (Mersenne Twister) for reproducibility
- Supports seeding for deterministic runs

Key functions:
```c
int init_rng(unsigned long seed);
double randn_std_normal(void);  // N(0,1)
double randu_uniform(void);     // U(0,1)
```

**Reproducibility**: Same seed produces identical random sequences across runs and platforms (when possible).

#### `print_routines.c / print_routines.h`
**Output and logging utilities.**

- Writes trajectory frames to output file (or stdout)
- Formats configuration summary
- Provides error logging and diagnostic output

Key functions:
```c
int write_frame_header(FILE *f, int frame_index, double time);
int write_particle_positions(FILE *f, const Particle *particles, int n);
int print_configuration_summary(const SimParams *params);
```

#### `array_utils.c / array_utils.h`
**Utility functions for array operations.**

- Memory-efficient array manipulation
- Distance calculations respecting periodic boundaries
- Utility functions for array setup and cleanup

Key functions:
```c
double compute_distance(const double *pos1, const double *pos2);
double *allocate_positions(int n_particles);
void free_positions(double *positions);
```

#### `code_handling.c / code_handling.h`
**Error codes and handling utilities.**

- Centralized error code definitions
- Error message lookup
- Error propagation and logging

Error codes:
```c
#define ERR_SUCCESS           0
#define ERR_ALLOC_FAILED      1
#define ERR_INVALID_PARAMS    2
#define ERR_IO_FAILURE        3
#define ERR_OVERLAP_COLLISION 4
#define ERR_WALL_VIOLATION    5
```

#### `equilibration_manager.c / equilibration_manager.h`
**Equilibration phase management.**

- Runs particles to equilibrium before production run
- Monitors energy/force evolution for convergence
- Configurable via `sim_params.conf`

Key functions:
```c
int equilibrate_system(SimulationState *state, const SimParams *params);
int is_equilibrated(const SimulationState *state, double tol);
```

### Pluggable Physics Modules

These modules implement specific physics and are selected at compile time via Makefile.

#### Confinement Geometry Modules

Each module implements the same interface but provides a different 2D channel geometry.

**Available geometries:**

1. **Cosine Channel** (`conf_cos.c / conf_cos.h`)
   - Channel width varies smoothly as W(x) = W₀ + A·cos(2πx/λ)
   - Smooth periodic geometry
   - Parameters: W₀ (mean width), A (amplitude), λ (wavelength)

2. **Septagon Channel** (`conf_sept.c / conf_sept.h`)
   - Septagonal (7-sided) geometry repeated periodically
   - Sharp corners and edges
   - Parameters: polygon size, orientation

3. **Splitter Channel** (`conf_splitter.c / conf_splitter.h`)
   - Two parallel channels with barrier in center region
   - Models particle splitter dynamics
   - Parameters: barrier width, splitter position

**Confinement Module Interface:**

Every confinement module must implement:

```c
/* Initialize confinement geometry with parameters from config */
int init_confinement(const ConfinementParams *params);

/* Check if position is valid (within channel, not in wall) */
int is_position_valid(const double *pos);

/* Apply wall force (repulsive) at position */
void compute_wall_force(const double *pos, double *force);

/* Get channel width at x-coordinate */
double get_channel_width(double x);

/* Apply periodic boundary conditions in x-direction */
void apply_periodic_bc(double *pos);
```

#### Interaction Modules

Each module implements particle-particle interactions.

**Available interactions:**

1. **Hard-Sphere Repulsion** (`int_hardspheres.c / int_hardspheres.h`)
   - Particles repel at contact (radius σ)
   - Zero interaction energy if not overlapping
   - Overlaps handled by collision detection
   - Parameter: particle radius σ

2. **Lennard-Jones Potential** (`int_lennardjones.c / int_lennardjones.h`)
   - LJ(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
   - Soft repulsion + weak attraction
   - Smoother dynamics than hard-sphere
   - Parameters: well depth ε, minimum distance σ_min

**Interaction Module Interface:**

Every interaction module must implement:

```c
/* Compute interaction force between two particles */
void compute_interaction_force(const double *pos1, const double *pos2,
                               double *force_on_1);

/* Compute interaction energy between two particles */
double compute_interaction_energy(const double *pos1, const double *pos2);

/* Check for overlaps (hard-sphere) or other validity checks */
int check_interaction_validity(const double *pos1, const double *pos2);
```

## Build System

### Makefile Configuration

The Makefile (in `src/`) selects which confinement and interaction modules to link:

```makefile
# Select geometry
GEOMETRY = splitter  # Options: splitter, cos, sept

# Select interaction
INTERACTION = lennardjones  # Options: lennardjones, hardspheres

# Compile with selections
OBJECTS = main_brownconf.o random_numb_gen.o print_routines.o \
          array_utils.o code_handling.o equilibration_manager.o \
          conf_$(GEOMETRY).o int_$(INTERACTION).o

main_brownconf: $(OBJECTS)
	gcc $(CFLAGS) -o main_brownconf $(OBJECTS) $(LDFLAGS)
```

### Python Makefile Generator

**Tool**: `python_tools/dtool_create_header_makefile.py`

This interactive script guides the user through:
1. Selecting C compiler (gcc, clang)
2. Selecting confinement geometry
3. Selecting interaction model
4. Setting optimization flags

Output: Custom Makefile in `src/` ready for compilation.

## Configuration System

### Configuration File Format (sim_params.conf)

Text-based configuration file specifying all simulation parameters:

```conf
# Simulation parameters for Brownian dynamics

# Particle parameters
N_PARTICLES 100
PARTICLE_RADIUS 0.05

# System parameters
CHANNEL_WIDTH 2.0
CHANNEL_LENGTH 10.0
TEMPERATURE 300.0

# Interaction parameters (geometry & interaction specific)
LENNARD_JONES_EPSILON 2.0
LENNARD_JONES_SIGMA_MIN 0.1

# Time-stepping
DT 0.001
N_STEPS 100000
N_EQUILIBRATION_STEPS 10000

# Output
OUTPUT_FILE trajectory.h5
WRITE_INTERVAL 100

# Physical forces
EXTERNAL_FORCE_X 40.0
EXTERNAL_FORCE_Y 0.0

# Collision handling
REJECT_OVERLAP 1       # 1: reject moves, 0: resample random numbers
OVERLAP_TOLERANCE 1e-12

# Random number seed
RNG_SEED 12345
```

### Configuration Parsing

- Configuration file is human-readable and editable
- Parsed line-by-line in `main_brownconf.c`
- Unknown parameters are logged as warnings
- Invalid values trigger errors with meaningful messages

## Initialization Flow

```
main()
  │
  ├─→ parse command-line args
  │   └─→ get config filename (or use default)
  │
  ├─→ read_configuration()
  │   └─→ parse sim_params.conf
  │   └─→ validate parameter ranges
  │   └─→ populate SimParams struct
  │
  ├─→ init_rng(RNG_SEED)
  │   └─→ seed Mersenne Twister
  │
  ├─→ init_confinement()
  │   └─→ load geometry module (compile-time selected)
  │   └─→ validate geometry parameters
  │
  ├─→ initialize_simulation()
  │   ├─→ allocate particle arrays
  │   ├─→ randomize initial positions (uniform or Boltzmann)
  │   └─→ initialize velocity (from Boltzmann distribution)
  │
  ├─→ equilibrate_system() [if N_EQUILIBRATION_STEPS > 0]
  │   ├─→ run equilibration loop
  │   ├─→ monitor energy/force evolution
  │   └─→ check for convergence
  │
  ├─→ TIME-STEPPING LOOP:
  │   │
  │   ├─→ for step = 0 to N_STEPS:
  │   │   │
  │   │   ├─→ for particle i:
  │   │   │   ├─→ generate random normal N(0,1) for each dimension
  │   │   │   ├─→ compute wall forces (confinement module)
  │   │   │   ├─→ compute interaction forces (interaction module)
  │   │   │   ├─→ compute external forces (from config)
  │   │   │   ├─→ Euler-Maruyama: x_new = x + F*dt + sqrt(2Ddt)*dW
  │   │   │   └─→ handle boundary conditions (periodic x, walls y)
  │   │   │
  │   │   ├─→ check for collisions (module-dependent)
  │   │   │   ├─→ if reject mode: restore previous positions
  │   │   │   └─→ if resample mode: regenerate random numbers
  │   │   │
  │   │   └─→ if step % WRITE_INTERVAL == 0:
  │   │       └─→ write_trajectory_frame()
  │   │           └─→ append positions to output file
  │   │
  │   └─→ cleanup & finalize
  │       ├─→ free particle arrays
  │       ├─→ close output files
  │       └─→ print summary statistics
  │
  └─→ exit(SUCCESS)
```

## Data Flow

### Simulation Run

```
sim_params.conf ─→ read_configuration() ─→ SimParams struct
                    
SimParams + RNG ─→ initialize_simulation() ─→ SimulationState
                    (particle positions, velocities, forces)

SimulationState ──→ [time-stepping loop] ──→ Updated state
                        per-step updates       (next time slice)

Updated state ────→ write_trajectory_frame() ──→ trajectory.h5
(every N steps)         (using h5py wrapper)      (HDF5 format)
```

### Post-Processing (Python)

```
trajectory.h5 ─→ load_trajectory() ─→ numpy array
                                       (n_frames, n_particles, 2)

numpy array ───→ compute_statistics() ─→ mean disp, RDF, etc.
                 analyze_trajectory()
                 

statistics ────→ plot_trajectory() ─→ PNG / PDF visualization
               → animate_trajectory() ─→ MP4 video animation
               → generate_report() ────→ HTML report
```

## Memory Layout

### Particle State Structure

```c
typedef struct {
    double x, y;           /* Current position */
    double vx, vy;         /* Velocity (used for debugging) */
    int id;                /* Particle ID for tracking */
} Particle;
```

**Storage**: Particles stored in contiguous arrays for cache efficiency:
- `double *positions` — shape (n, 2), x and y coords interleaved
- `double *velocities` — shape (n, 2), or not stored if velocities not tracked
- `double *forces` — shape (n, 2), forces used for next time step

### Thread Safety & Parallelization

**Current state**: Single-threaded C with MPI for distributed runs.

**Future parallelization**:
1. OpenMP: parallelize per-particle loops with `#pragma omp parallel for`
2. CUDA: offload Euler-Maruyama integration and force calculations to GPU

**Constraints**:
- Random number sequences must remain reproducible after parallelization
- Thread-safe RNG (e.g., per-thread seeds or thread-safe streams)

## Extension Points

### Adding a New Confinement Geometry

1. Create `conf_newgeometry.{c,h}`
2. Implement required interface functions (see above)
3. Update Makefile with geometry selection
4. Update `dtool_create_header_makefile.py` with new option
5. Test against reference simulations

**Detailed workflow**: See `SKILLS.md` → "Adding a New Confinement Geometry"

### Adding a New Interaction Model

1. Create `int_newinteraction.{c,h}`
2. Implement required interface functions
3. Update Makefile with interaction selection
4. Update `dtool_create_header_makefile.py` with new option
5. Validate against reference simulations

**Detailed workflow**: See `SKILLS.md` → "Adding a New Interaction Model"

## Performance Considerations

### Cache Efficiency
- Particles stored in contiguous arrays for spatial locality
- Force/position arrays accessed sequentially in time-stepping loop
- Avoid unnecessary copies between temporary buffers

### Numerical Stability
- Preserve double-precision throughout calculations
- Respect minimum-image convention for all distances
- Handle particle overlaps consistently (reject or resample)

### Reproducibility
- Identical random seed → identical trajectory (within floating-point precision)
- No platform-specific operations in core simulation
- Time-stepping order matters (force calculation before position update)

## References

- **Detailed Instructions**: See `instructions/c.instructions.md`
- **Numerical Methods**: See `docs/numerics.md`
- **Output Formats**: See `docs/data-format.md`
- **Common Workflows**: See `SKILLS.md`

---

**Version**: 1.0  
**Last Updated**: June 2026
