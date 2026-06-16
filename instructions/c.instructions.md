# C Coding Standards & Guidelines

This document specifies standards for C code in `src/` directory of the Spatially Confined Brownian Motion simulator.

## Compilation Requirements

All C code **must** compile cleanly without warnings using:

```bash
gcc -Wall -Wextra -Wpedantic -Werror
```

- Do not suppress warnings with pragmas unless absolutely necessary and documented
- Do not introduce compiler-specific extensions (e.g., GCC attributes) unless already present in the codebase
- Avoid platform-specific code unless necessary for cross-platform support

## Code Style

### Naming Conventions

- **Functions**: `snake_case_for_functions()`
- **Variables**: `snake_case_local_vars`, `global_with_prefix_g_`
- **Constants & Macros**: `UPPER_CASE_CONSTANTS`
- **Types**: `PascalCaseTypes` for struct/typedef names
- **Files**: `lower_case_with_underscores.{c,h}`

### Formatting

- Use 4 spaces for indentation (never tabs)
- Line length: preferably ≤80 characters (maximum ≤100 for unavoidable cases)
- One declaration per line (except loop counters)
- Braces: opening brace on same line for functions; opening brace on new line for control structures

```c
int compute_interaction(double x, double y) {
    if (x > y) {
        return calculate_value(x);
    }
    return calculate_value(y);
}
```

### Comments

- Comment **why**, not **what** (the code shows what it does)
- Use block comments (`/* */`) for multi-line explanations
- Use line comments (`//`) only for single-line clarifications
- Document non-obvious numerical constants and unit assumptions

```c
/* Minimum-image convention: adjust distance to account for periodic
   boundary conditions along x-axis. See docs/numerics.md for details. */
double dx_adjusted = dx > channel_width / 2.0 
    ? dx - channel_width 
    : dx;

// Pre-allocate 10% extra for thermal fluctuations
double *velocity = malloc(n_particles * 1.1 * sizeof(double));
```

## Memory Management

### Allocation Ownership

Every dynamic allocation must have a **clear, documented owner**:

1. **Ownership is explicit**: Document in function signatures and comments who allocates and who frees
2. **Check failures**: Always check return value of `malloc`, `calloc`, `realloc`
3. **Free in reverse order**: Free in the opposite order of acquisition
4. **No leaks**: Every allocation path must have a corresponding free path

```c
/* Allocate particle state structure. Caller is responsible for freeing
   via free_particle_state(). Returns NULL on failure. */
ParticleState *init_particle_state(int n_particles, int *error_code) {
    ParticleState *ps = malloc(sizeof(ParticleState));
    if (!ps) {
        *error_code = ERR_ALLOC_FAILED;
        return NULL;
    }
    
    ps->positions = malloc(n_particles * 2 * sizeof(double));
    if (!ps->positions) {
        free(ps);  // Free in reverse order
        *error_code = ERR_ALLOC_FAILED;
        return NULL;
    }
    
    ps->n_particles = n_particles;
    return ps;
}

void free_particle_state(ParticleState *ps) {
    if (!ps) return;
    free(ps->positions);  // Free components first
    free(ps);             // Then free container
}
```

### Stack vs Heap

- Prefer stack allocation for fixed-size, small temporary objects
- Use heap only for data with unbounded or large size
- Avoid unnecessary copies; pass pointers to large structures
- Document why each heap allocation is necessary

```c
/* Small temporary: use stack */
double buffer[3];  // Fixed size, small

/* Large unbounded data: use heap */
double *trajectories = malloc(n_steps * n_particles * 2 * sizeof(double));

/* Return by pointer to avoid copy */
void compute_distances(const Particle *p1, const Particle *p2, 
                       double *result) {
    result[0] = p1->x - p2->x;
    result[1] = p1->y - p2->y;
}
```

## Error Handling

### Return Codes

Functions that can fail must return explicit error codes:

```c
/* Return codes for error handling */
#define ERR_SUCCESS           0
#define ERR_ALLOC_FAILED      1
#define ERR_INVALID_PARAMS    2
#define ERR_IO_FAILURE        3
#define ERR_OVERLAP_COLLISION 4
#define ERR_WALL_VIOLATION    5
```

- **Return 0 on success**, non-zero on error
- **Do not use exit() or abort()** in library code—propagate errors to caller
- **Output pointer parameters** can return additional data; use return value for error status

```c
int update_particle_positions(Particle *particles, int n, double dt,
                               const SimParams *params) {
    if (!particles || !params) {
        return ERR_INVALID_PARAMS;
    }
    
    for (int i = 0; i < n; i++) {
        int err = validate_position(&particles[i], params);
        if (err != ERR_SUCCESS) {
            return err;  // Propagate error
        }
    }
    
    return ERR_SUCCESS;
}
```

### Error Logging

- Log meaningful error messages with context (e.g., which particle, which iteration)
- Use `fprintf(stderr, ...)` for errors in library functions
- Include function name and likely cause

```c
if (fread(buffer, sizeof(double), expected_count, file) != expected_count) {
    fprintf(stderr, 
            "ERROR: read_trajectory_frame() failed to read %d values from %s\n",
            expected_count, filename);
    return ERR_IO_FAILURE;
}
```

## Numerics & Physics

### Reproducibility is Paramount

- Preserve exact numerical behavior across runs and platforms (when possible)
- Do not modify stochastic integration algorithms without explicit request
- Document all unit assumptions in code and comments
- Avoid hidden unit conversions; be explicit about scaling factors

### Precision

- Use `double` by default for floating-point calculations
- Use `float` only when explicitly justified (storage, not computation)
- Document precision choices, especially for distance thresholds

```c
/* Use double for all physical calculations to preserve reproducibility.
   Tolerance thresholds must be consistent with double precision. */
double OVERLAP_TOLERANCE = 1e-12;  /* Must match double precision limits */

/* Never mix float and double in same calculation; always promote to double */
double distance = sqrt((double)dx * dx + (double)dy * dy);
```

### Stochastic Integration

The codebase uses **Euler-Maruyama scheme** for overdamped Brownian dynamics:

```
x_{n+1} = x_n + D * ∇U(x_n) * dt + √(2 * D * dt) * N(0, 1)
```

Where:
- `D` is the diffusion coefficient
- `∇U` is the force field
- `dt` is the time step
- `N(0, 1)` is a standard normal random number

**Do not change this scheme without explicit request.** Higher-order schemes (Milstein, etc.) are unnecessary for this system (no inertia, no time-dependent forces).

### Periodic Boundary Conditions & Minimum-Image Convention

All distance calculations must respect **minimum-image convention**:

```c
double apply_minimum_image(double dx, double L) {
    /* Adjust distance to account for periodic boundaries.
       Ensures we use the shortest path across the periodic boundary. */
    if (dx > L / 2.0) {
        return dx - L;
    }
    if (dx < -L / 2.0) {
        return dx + L;
    }
    return dx;
}
```

The channel is periodic in **x-direction only**; y-direction has hard walls.

## Module Organization

### File Structure

Each physics module (confinement geometry, interaction model) has:
- **Header file** (`.h`): Public interface, struct definitions, function declarations
- **Implementation file** (`.c`): Function implementations, static helpers

### Header Files

- Include guard using `#ifndef FILENAME_H` convention
- Declare only public API and essential types
- Use `static inline` for simple utility functions
- Document all public functions with purpose, parameters, return value

```c
/* conf_splitter.h - Splitter channel confinement geometry */
#ifndef CONF_SPLITTER_H
#define CONF_SPLITTER_H

#include <stdint.h>

/* Initialize splitter channel geometry with given parameters */
int init_splitter_channel(double width, double length, 
                          SplitterChannel *channel);

/* Check if particle position is valid (within channel bounds) */
int is_position_valid_splitter(const double *pos, 
                               const SplitterChannel *channel);

#endif
```

### Implementation

- Keep functions focused and single-purpose
- Use static functions for module-internal helpers
- Minimize global state; pass configuration via parameters
- Document tricky algorithms with references to papers or docs/

```c
static double compute_interaction_energy_hardsphere(const Particle *p1,
                                                     const Particle *p2,
                                                     double sigma) {
    /* Hard-sphere repulsion: zero if overlapping, otherwise infinite.
       In practice, we handle overlaps via collision detection. */
    double dx = apply_minimum_image(p1->x - p2->x, CHANNEL_WIDTH);
    double dy = p1->y - p2->y;
    double r_sq = dx * dx + dy * dy;
    return r_sq < sigma * sigma ? 1e308 : 0.0;  // Represent infinity
}
```

## Interfacing with Confinement & Interaction Modules

The build system selects one confinement and one interaction module at compile time:

- **Confinement modules**: `conf_splitter.{c,h}`, `conf_cos.{c,h}`, `conf_sept.{c,h}`
- **Interaction modules**: `int_hardspheres.{c,h}`, `int_lennardjones.{c,h}`

Each module implements a standard interface (see `docs/architecture.md`). When adding a new module:

1. Implement all required functions matching the interface
2. Add selection to `makefile` and `dtool_create_header_makefile.py`
3. Test against reference simulations to verify numerical correctness

## Common Patterns

### Safe Iteration

```c
for (int i = 0; i < n_particles; i++) {
    if (particles[i].id < 0) {
        fprintf(stderr, "ERROR: Invalid particle ID at index %d\n", i);
        return ERR_INVALID_PARAMS;
    }
    /* Process particle */
}
```

### Safe File I/O

```c
FILE *f = fopen(filename, "rb");
if (!f) {
    fprintf(stderr, "ERROR: Cannot open file %s for reading\n", filename);
    return ERR_IO_FAILURE;
}

size_t n_read = fread(buffer, sizeof(double), count, f);
if (n_read != count) {
    fprintf(stderr, "ERROR: Expected %zu doubles, read %zu\n", count, n_read);
    fclose(f);
    return ERR_IO_FAILURE;
}

if (fclose(f) != 0) {
    fprintf(stderr, "ERROR: Failed to close file %s\n", filename);
    return ERR_IO_FAILURE;
}
```

### Configuration Parsing

Use the `sim_params.conf` format for all configuration. Parse via provided utility functions; do not hardcode parameters:

```c
int read_simulation_parameters(const char *config_file, 
                               SimParams *params) {
    if (!config_file || !params) {
        return ERR_INVALID_PARAMS;
    }
    /* Parse config_file line-by-line, populating params struct */
    return ERR_SUCCESS;
}
```

## Testing & Validation

### Manual Validation

Since there is no automated test suite, validate manually by:

1. **Numerical reproducibility**: Run same simulation twice, verify bit-for-bit identical output
2. **Sanity checks**: Verify energy conservation, particle densities, trajectory bounds
3. **Reference comparison**: Compare output against known correct simulations
4. **Edge cases**: Test zero-force, single-particle, high-density scenarios

### Debugging

- Add detailed tracing output under conditional compilation (`#ifdef DEBUG`)
- Use `fprintf(stderr, ...)` for diagnostics
- Never leave debug output enabled in committed code

```c
#ifdef DEBUG
fprintf(stderr, "DEBUG: Particle %d at (%.6f, %.6f), force=(%.6e, %.6e)\n",
        i, particles[i].x, particles[i].y, forces[i][0], forces[i][1]);
#endif
```

## References

- **Physics Model**: See `docs/numerics.md` for stochastic integration and boundary conditions
- **Architecture**: See `docs/architecture.md` for module structure and initialization flow
- **Build System**: See main `README.md` and `Makefile` for compilation options
- **Data Formats**: See `docs/data-format.md` for output file specifications

---

**Version**: 1.0  
**Last Updated**: June 2026
