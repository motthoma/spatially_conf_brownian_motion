# Numerical Methods & Physical Model

This document describes the physical model, numerical integration scheme, boundary conditions, and key numerical assumptions underlying the Spatially Confined Brownian Motion simulator.

## Physical Model Overview

### Overdamped Brownian Dynamics

The simulator models **overdamped Brownian motion** where particles move under the influence of:
1. **Thermal random forces** (from molecular collisions)
2. **Interaction forces** (particle-particle and wall confinement)
3. **External forces** (e.g., electric or gravitational fields)
4. **Viscous drag** (implicit through overdamped approximation)

**Key assumption**: **No inertia** — particles instantly respond to forces without acceleration. The velocity is proportional to force:

```
v = F / ζ

where ζ is the friction coefficient (related to viscosity)
```

This **overdamped limit** is valid when:
- Particle size is small (nm scale)
- System is viscous (e.g., water or polymer solutions)
- Relaxation time is microseconds (much smaller than simulation timescale)

### Generalized Langevin Equation

The position evolution follows the **generalized Langevin equation**:

```
m dv/dt = -ζ v + F_interaction + F_external + F_thermal

Overdamped limit (m → 0):
0 = -ζ v + F_interaction + F_external + F_thermal

Solving for v:
v = (F_interaction + F_external + F_thermal) / ζ

dx/dt = (F_interaction + F_external + F_thermal) / ζ

Or equivalently:
dx/dt = D ∇U + F_external/ζ + √(2D) dW/dt
```

Where:
- **D** is the diffusion coefficient (related to ζ by Einstein relation: D = k_B T / ζ)
- **U** is the interaction potential energy
- **∇U** is the negative interaction force
- **dW/dt** is Gaussian white noise

## Time Integration: Euler-Maruyama Scheme

### Stochastic Differential Equation

The position evolution is governed by:

```
dx = [D ∇U(x) + a(x)] dt + √(2D) dW_t
```

Where:
- First term: drift (deterministic forces)
- Second term: diffusion (random thermal motion)
- `a(x)` includes external forces and confinement repulsion
- `dW_t` is a Wiener process increment (standard Brownian increment)

### Discretization: Euler-Maruyama Method

**Discrete time-stepping**:

```
x_{n+1} = x_n + D ∇U(x_n) Δt + a(x_n) Δt + √(2 D Δt) ε_n

where:
  - ε_n ~ N(0, 1) independent standard normal random number
  - Δt is time step
  - ∇U(x_n) is force evaluated at current position
```

**Why Euler-Maruyama?**

1. **Simplicity**: Easy to implement and verify
2. **Accuracy**: Weak order 1.0 (sufficient for our purposes)
3. **Stability**: Unconditionally stable for our system (overdamped, bounded forces)
4. **Reproducibility**: Deterministic given RNG seed

**Why not higher-order schemes?**

- Milstein scheme (order 1.5) requires derivatives of force (∇²U), not needed here
- Runge-Kutta schemes add complexity without benefit for this system
- For overdamped dynamics without time-dependent forces, Euler-Maruyama convergence is adequate

### Time Step Selection

The time step must be small enough to:
1. **Resolve dynamics**: Δt << characteristic timescale of system
2. **Preserve stability**: Ensure forces don't cause instability
3. **Maintain accuracy**: Error budget acceptable for statistics

**Typical choices**:
- Δt ~ 0.001 (in simulation units)
- Must be re-validated when changing system parameters

**Validation**:
- Run ensemble of simulations with Δt and Δt/2
- Check that statistics (mean displacement, diffusion coefficient) converge
- If not converged, reduce Δt

## Boundary Conditions

### Periodic Boundary Conditions (x-direction)

The channel is periodic along the transport direction (x-axis):

```
If x > L_x:   x_new = x - L_x  (particle wraps to left side)
If x < 0:     x_new = x + L_x  (particle wraps to right side)

Domain: x ∈ [0, L_x) with wraparound
```

**Why periodic?**

- Models an infinite channel by replication
- Eliminates boundary effects at inlet/outlet
- Simplifies statistics (no need for entrance/exit effects)

**Minimum-Image Convention**

When computing distances between particles, use the **shortest distance** accounting for periodicity:

```
dx_raw = x1 - x2
dx_imaged = dx_raw if |dx_raw| ≤ L_x/2
dx_imaged = dx_raw - L_x if dx_raw > L_x/2
dx_imaged = dx_raw + L_x if dx_raw < -L_x/2

distance = sqrt(dx_imaged² + dy²)
```

**Critical**: This must be applied in **all distance calculations** for forces and overlap checking.

### Hard Wall Boundaries (y-direction)

The channel has hard walls at y = 0 and y = W(x), where W(x) is the local channel width:

```
Confinement: 0 ≤ y ≤ W(x)

Wall repulsion:
- If y < 0 or y > W(x): apply strong repulsive force
- Particles cannot penetrate walls (hard sphere approximation)
```

**Wall force** (implementation-dependent):

Option 1: Hard sphere model
- Particle radius σ
- Wall contact occurs at y = σ or y = W(x) - σ
- Repulsive force (infinite stiffness) applied on violation

Option 2: Soft wall model (Lennard-Jones like)
- Gradual repulsion as particle approaches wall
- More realistic but computationally more expensive

## Particle Interactions

### Interaction Potential Energy

Particles interact via a pairwise potential U(r):

```
E_interaction = Σ_{i<j} U(r_ij)

where r_ij = distance between particles i and j
```

### Hard-Sphere Model

**Potential**:
```
U_HS(r) = ∞ if r < 2σ  (particles cannot overlap)
U_HS(r) = 0 if r ≥ 2σ  (no interaction beyond contact)

where σ is particle radius
```

**Force**:
```
F_HS(r) = -dU/dr = 0 everywhere (except discontinuity at 2σ)
```

**Implementation**:
- No attractive force
- Overlaps detected and rejected (configuration space exploration)
- Particles modeled as hard impenetrable spheres

**Collision handling** (configurable):

1. **Reject-on-overlap mode**:
   - If new position causes overlap, keep old position
   - Particle "bounces back"
   - Efficiency: moderate; rejection rate increases with density

2. **Resample-on-overlap mode**:
   - If new position causes overlap, resample random numbers
   - Regenerate Brownian noise until valid position found
   - Efficiency: better at high density; can be slow if very crowded

### Lennard-Jones Model

**Potential**:
```
U_LJ(r) = 4ε [ (σ/r)¹² - (σ/r)⁶ ]

where:
  - ε = well depth (energy scale)
  - σ = characteristic length scale
```

**Force**:
```
F_LJ(r) = -dU/dr = 24ε/r [ 2(σ/r)¹² - (σ/r)⁶ ]

Repulsive at small r (positive F, pushes apart)
Attractive at large r (negative F, pulls together)
Zero at r = 2^(1/6) σ ≈ 1.122 σ
```

**Characteristics**:
- Soft repulsion: particles can slightly overlap
- Long-range attraction: weak binding at distance
- Smooth force landscape: numerically stable
- Physical relevance: approximates van der Waals interactions

**Computational considerations**:
- Cutoff radius needed for efficiency (typically r_cut ~ 2.5σ)
- Pairs with r > r_cut: U ≈ 0, force ≈ 0
- Reduces computational complexity from O(N²) to O(N)

## External Forces

The simulator supports static external forces (e.g., electric field):

```
F_external = (F_x, F_y)  [constant across all particles]

Drift term in Euler-Maruyama:
v_drift = F_external / ζ
```

**Common usage**:
- Apply force F_x along transport direction (x) to induce flow
- Measure response as mean drift velocity or particle flux
- Study driven systems far from equilibrium

## Confinement Geometries

### Available Geometries

Each geometry is a 2D channel with varying width W(x).

#### 1. Rectangular Channel
- Constant width: W(x) = W
- Simplest case for validation
- Parameters: W (width), L (length)

#### 2. Cosine Channel
- Sinusoidal width variation: W(x) = W₀ + A cos(2πx/λ)
- Smooth, periodic geometry
- Parameters: W₀ (mean width), A (amplitude), λ (wavelength)
- Applications: smooth constrictions, barriers

#### 3. Septagon Channel
- Septagonal (7-sided) polygon repeated periodically
- Sharp corners and edges (challenging for particle dynamics)
- Parameters: side length, orientation
- Applications: study corner effects, numerical stability

#### 4. Splitter Channel
- Two parallel channels with barrier in middle region
- Models a flow splitter or partition
- Parameters: barrier width, splitter position along x
- Applications: study particle routing, bifurcation dynamics

### Distance Calculations with Geometry

Wall distance determines wall forces:

```c
/* For cosine channel */
double wall_distance_top = W(x) - y;
double wall_distance_bot = y - 0;

/* Apply repulsive force if close to wall */
double f_wall_top = repulsive_force(wall_distance_top);
double f_wall_bot = repulsive_force(wall_distance_bot);
```

## Numerical Reproducibility

### Key Principles

1. **Same seed → Same trajectory** (bit-for-bit reproducible)
2. **Deterministic RNG**: Use Mersenne Twister with fixed seed
3. **Fixed arithmetic**: Use double precision throughout
4. **No adaptive algorithms**: Fixed time step, no adaptive refinement
5. **Consistent force evaluation**: Forces computed in consistent order

### Achieving Reproducibility

**In C code**:
- Initialize RNG with fixed seed before each run
- Use sequential RNG calls (no random-access to stream)
- Evaluate forces in consistent particle order (0 to N-1)
- No thread parallelization (would require careful RNG management)

**Verification**:
- Run same simulation twice with same seed
- Compute difference in trajectories
- Should be < 1e-15 (floating-point precision limit)

**Pitfalls to avoid**:
- Reading RNG from file (can have platform-dependent formatting)
- Reordering particle force calculations (breaks reproducibility)
- Mixing single/double precision (promotes/demotes differently)
- Uninitialized memory (use calloc, not malloc)

## Numerical Stability & Error Bounds

### Error Analysis

**Local truncation error** (per step):
```
ε_local ~ O(Δt^(3/2))  [Euler-Maruyama local error]
```

**Global error** (accumulated over time):
```
ε_global ~ O(Δt^1/2)   [Weak convergence order]
```

**For Δt = 0.001 and T_sim = 10000 steps (time = 10 units)**:
- Relative error ~ √(0.001) ~ 0.03 (3%)

### Stability Condition

For overdamped dynamics with bounded forces, Euler-Maruyama is **unconditionally stable**:
- No CFL condition required
- No upper limit on Δt for stability (accuracy is limiting factor)

### Force Clipping

If forces become very large (e.g., near wall overlap), consider:

1. **Soft wall model**: Use smooth repulsion instead of hard walls
2. **Smaller time step**: Reduce Δt if forces near limit
3. **Pre-equilibration**: Allow system to relax before data collection

## Physical Units

### Characteristic Scales

Define a unit system:

```
Length scale:   L₀ = 1 μm  (characteristic particle separation)
Time scale:     T₀ = 1 ms  (characteristic diffusion time)
Energy scale:   E₀ = k_B T (thermal energy)
Force scale:    F₀ = E₀/L₀
```

**Parameters expressed in these units**:

```
Dimensionless time step:  Δt* = Δt / T₀
Diffusion coefficient:    D* = D / (L₀²/T₀)
External force:           F_ext* = F_ext / F₀
Particle radius:          σ* = σ / L₀
```

### Assumptions

- All simulation units are dimensionless (relative to characteristic scales)
- Document your choice of characteristic scales in output metadata
- When interpreting results, remember to multiply by L₀, T₀, etc.

## Validation & Benchmarking

### Analytical Benchmarks

1. **Single particle in uniform channel** (no interactions):
   - Mean displacement scales as √t (random walk)
   - Diffusion coefficient D = <x²> / (2t)

2. **Harmonic potential**:
   - Equilibrium distribution: Boltzmann exp(-U/k_B T)
   - Diffusion coefficient measurable

3. **Periodic boundary conditions**:
   - Verify no spurious correlations at boundary
   - Check force on opposite sides of domain

### Numerical Tests

1. **Convergence**: Run with Δt and Δt/2, verify statistics converge
2. **Energy conservation**: Check total mechanical energy (if defined)
3. **Density constancy**: Verify particle density remains constant
4. **Reproducibility**: Run twice with same seed, compare trajectories

### Reference Simulations

Store reference trajectories for regression testing:
- Document hardware and RNG seed used
- Include bit-exact trajectory (or bit-different tolerance)
- Use for validation after code changes

## References

### Physics Background
- Gardiner, C. W. (1997). *Handbook of Stochastic Methods*. Springer.
- van Kampen, N. G. (2007). *Stochastic Processes in Physics and Chemistry*. Elsevier.

### Numerical Methods
- Kloeden, P. E., & Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations*. Springer.
- Mannella, R. (2002). Integration of stochastic differential equations on a computer. *International Journal of Modern Physics C*, 13(09), 1177-1194.

### Simulation Techniques
- Allen, M. P., & Tildesley, D. J. (2017). *Computer Simulation of Liquids*. Oxford University Press.
- Frenkel, D., & Smit, B. (2002). *Understanding Molecular Simulation*. Academic Press.

---

**Version**: 1.0  
**Last Updated**: June 2026
