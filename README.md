# Spatially Confined Brownian Motion 

Simulation code for **interacting** and **non-interacting** Brownian particles in 2D periodic channels.  
The particles undergo overdamped Brownian motion simulated via a stochastic Euler scheme.  
External static forces can be applied, but long-range hydrodynamic interactions are **not** included.  
Periodic boundary conditions are used, and particle density is kept constant in the interacting case.

<p align="center">
  <img src="./logo.png" alt="logo_conf_brownian_motion" width="400">
</p>


Exemplary system: Two hardspheres move within a confinement with a saw-tooth profile
as implemented in the conf_splitter module

---

## 📂 Project Structure

The repository consists of three main directories:

1. **`src/`**
   - Contains the main simulation code, Makefile, and latest results.
   - Each simulation run creates a new directory containing the code and data.

2. **`doxygen/`**
   - Contains the [Doxygen](https://www.doxygen.nl/) configuration with a Doxyfile that can be used to generate documentation.
   - Run `doxygen` to generate:
     - HTML documentation
     - PDF (via LaTeX)

3. **`profiling/`**
   - Provides an environment for profiling runtime performance with **gprof**.
   - See the README inside for details.

3. **`visualization/`**
   - Contains python scripts for data visualization. Plotting is done with `PyX`, or `Matplotlib` or `Seaborn` (TBD).

---

## ⚙️ Compilation & Execution

The code is modular to support different models:

- Channel shapes: prefixed with `conf_`
- Particle interactions:
  - Hard-core interactions (`int_...`)
  - Lennard-Jones interactions (`int_...`)

If **MPI** is available (e.g., on the *Albeniz* cluster), the code can be parallelized.

### 🔨 Building
A Makefile is used for compilation.  
To simplify setup, use:

```bash
python dtool_create_makefile.py


