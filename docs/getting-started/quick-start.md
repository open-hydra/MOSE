# Quick Start: Prandtl–Meyer Expansion

This guide walks you through running the **Prandtl–Meyer expansion** test case — a canonical 2-D supersonic flow verification that demonstrates corner-based isentropic expansion with analytic validation.

!!! note "Prerequisites"
    - MOSE built and available at `bin/MOSE` (see [Installation](installation.md))
    - Python 3 with NumPy (for post-run verification)

<figure>
  {% include "getting-started/mach.svg" %}
  Contours of Mach number from the Prandtl–Meyer expansion test case.
</figure>

---

## Physics Overview

The Prandtl–Meyer expansion describes the isentropic turning of a supersonic flow over a smooth convex corner. As the flow expands, the Mach number increases and pressure decreases while entropy remains constant.

**Canonical parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| Upstream Mach | $M_1 = 2.0$ | Supersonic inlet condition |
| Expansion angle | $\theta \approx 35°$ | Corner deflection |
| Specific heat ratio | $\gamma = 1.4$ | Diatomic gas (air) |
| Upstream T, p | $T_1 = 1$ K, $p_1 = 1$ Pa | Normalized reference state |

**Expected result:** The downstream Mach $M_2$ computed by MOSE matches the analytic Prandtl–Meyer prediction to within **0.5%** relative error.

The analytic downstream Mach is found by solving:

$$\nu(M_2) = \nu(M_1) + \theta$$

where the Prandtl–Meyer function is:

$$\nu(M) = \sqrt{\frac{\gamma + 1}{\gamma - 1}} \arctan\sqrt{\frac{\gamma - 1}{\gamma + 1}(M^2 - 1)} - \arctan\sqrt{M^2 - 1}$$

---

## Case Directory Structure

The case is located at `test/2D/euler/prandtl-meyer/`:

```
test/2D/euler/prandtl-meyer/
├── input.ini           # Solver configuration
├── MOSE.sh             # Run script
├── check.py            # Verification script
├── INPUT/
│   ├── ic.tec          # Initial condition (mesh + field data)
│   └── bc.txt          # Boundary condition table
└── OUTPUT/             # (Created at runtime)
```

Shared thermodynamic data (phase table, transport properties) lives in `test/common/Air/` and is symlinked at runtime.

!!! info "Where do input files come from?"
    The input files (`ic.tec`, `bc.txt`, `phase.txt`, `thermo.dat`) were produced by **ATLAS**, the pre-processor of [Hydra](../overview.md#hydra-cfd-suite). ATLAS is distributed separately and is not part of the MOSE package. For this bundled test case all necessary files are already provided. If ATLAS is not available, files can be prepared manually following the formats described in the [User Guide](../user/using.md).

### Input Files Description

- **`input.ini`** — Solver configuration specifying physics (Euler equations), numerics (time scheme, Riemann solver, reconstruction), and I/O preferences
- **`ic.tec`** — Initial condition file in Tecplot format containing both the computational mesh **and** the flow field (density, pressure, velocity). MOSE reads the grid coordinates directly from this file — no separate mesh file is needed
- **`bc.txt`** — Boundary condition table defining edge or face-wise boundary treatment (e.g., inlet, exit, slip wall) referenced by the solver
- **`phase.txt`, `thermo.dat`** — Pre-compiled species properties (located in `test/common/Air/`). These files are handled via the **FLINT** library.

---

## Solver Configuration

The `input.ini` file specifies the physics and numerics:

```ini
[MOSE-Physics]
equations = euler            ; Inviscid compressible Euler

[MOSE-Numerics]
cfl = 1.0                    ; Courant–Friedrichs–Lewy number
time-scheme = RK3            ; 3rd-order Runge–Kutta temporal scheme
time-accurate = false        ; Steady-state integration (no physical time)
irs-beta = 0.5               ; Implicit residual smoothing (steady runs only)
riemann-solver = HLLE        ; HLLE numerical flux
space-reconstruction = MUSCL ; 2nd-order spatial reconstruction
flux-limiter = vanleer       ; van Leer slope limiter
integration-variables = cons ; Conserved variables (ρ, ρu, ρv, E)

[MOSE-IO]
shell-diter = 100            ; Diagnostic output every 100 iterations
sol-format = tecplot ascii   ; Tecplot ASCII output format

[MOSE-Probes]
probe1 = exit                ; Probe at the exit plane

[exit]
index-position = 1 50 33 1   ; Probe location (block, i, j, k)
variables = M                ; Mach number output
diter = 100                  ; Probe output frequency

```

Key choices:  
- Steady-state integration (`time-accurate = false`) — physical time is not needed; we iterate to steady state  
- HLLE solver — robust, entropy-satisfying choice for expansion waves  
- van Leer limiter — smooth and non-oscillatory  
- Probe output — captures Mach number at the exit plane for verification

---

## Running the Simulation

**Basic Execution**

Navigate to the case directory and run:

```bash
cd test/2D/euler/prandtl-meyer
./MOSE.sh solve
```

The `MOSE.sh` script:

1. Symlinks shared thermodynamic data (`phase.txt`, `thermo.dat`) into `INPUT/`
2. Creates `OUTPUT/` and `bin/` directories if needed
3. Launches the solver executable

Typical runtime: 5–30 seconds (depending on mesh size and hardware).

**Parallel Execution**

To run with N OpenMP threads:

```bash
./MOSE.sh solve -p 4    # Run with 4 threads
```

---

## Understanding the Output

### Diagnostics on screen

Before the iteration starts, MOSE prints a summary of the case set-up. This includes the number of blocks, cells, boundary faces, types of boundary conditions, physical model details (equations, gas model), and numerical scheme choices. This diagnostic is useful for verifying that the input files were read correctly and that the solver is configured as expected.

!!! tip "Loading checks"
    If any issue occurs during loading (e.g., missing files, format errors), MOSE will print an error message and exit, so the successful loading message is a good sign that the case is set up correctly.

```
 =========================================================================================
 Loading
 =========================================================================================
   Thermodynamics                 OK
   Transport                      OK
   Chemistry                      OK
   Input file                     OK
   Initial conditions             OK
   Boundary conditions            OK
 =========================================================================================

 =========================================================================================
 Set-up
 =========================================================================================
 Domain
   Blocks                         1
   Cells                          2500
   Boundary faces                 5200
   Mesh type                      2-D, plane

 Boundary conditions
   Inflow/outflow                 50
   Symmetry                       50
   Extrapolation                  100

 Input/Output
   Initial conditions file        INPUT/ic.tec
   Solution format                tecplot ascii
   Probes number                  1          
 =========================================================================================

 =========================================================================================
 Physical model
 =========================================================================================
 Gas model
   Equations                      Euler
   Equation of state              Ideal
   Thermodynamics                 Thermally perfect gas
   Species                        1
 =========================================================================================

 =========================================================================================
 Numerical scheme
 =========================================================================================
 Space
   Reconstruction                 MUSCL
   Flux limiter                   Van Leer

 Time
   Scheme                         Third-order Runge-Kutta
   Integration variables          Conservative
   Implicit residual smoothing    Beta set to 0.500000E+00

 Fluxes
   Riemann solver                 HLLE
 =========================================================================================
```

During execution, MOSE prints iteration count and residuals. In this case, the density residual decreases smoothly, indicating stable convergence to the steady-state solution:

```
 400 boundary of block 1, face 3, gas entering at iteration 1, faces: 33, largest |w|/a:  9.05293E-08
 400 boundary of block 1, face 3, gas entering at iteration 15, faces with |w|/a above 10 times its last line: 2, largest |w|/a:  1.07239E-06
 400 boundary of block 1, face 3, gas entering at iteration 55, faces with |w|/a above 10 times its last line: 3, largest |w|/a:  1.17654E-05
 400 boundary of block 1, face 3, gas entering at iteration 61, faces with |w|/a above 10 times its last line: 4, largest |w|/a:  1.39872E-04
 400 boundary of block 1, face 3, gas entering at iteration 68, faces with |w|/a above 10 times its last line: 1, largest |w|/a:  1.41811E-03
 400 boundary of block 1, face 3, gas entering at iteration 78, faces with |w|/a above 10 times its last line: 1, largest |w|/a:  1.52527E-02
MOSE | Iter =      100 | Global iter =      100 | Density residual = 0.293637E-03
 400 boundary of block 1, face 3, gas entering at iteration 101, faces with |w|/a above 10 times its last line: 1, largest |w|/a:  1.56579E-01
MOSE | Iter =      200 | Global iter =      200 | Density residual = 0.181007E-04
MOSE | Iter =      300 | Global iter =      300 | Density residual = 0.183076E-06
MOSE | Iter =      400 | Global iter =      400 | Density residual = 0.517294E-08
MOSE | Iter =      500 | Global iter =      500 | Density residual = 0.485840E-09
MOSE | Iter =      600 | Global iter =      600 | Density residual = 0.383678E-10
MOSE | Iter =      700 | Global iter =      700 | Density residual = 0.293705E-11
MOSE | Iter =      767 | Global iter =      767 | Density residual = 0.627384E-12

 =========================================================================================
 Timing
 =========================================================================================
   Iterations                     767
   Solver                          8.91734E+00 s
   Solver per iteration            1.16263E-02 s
   Elapsed                         8.94817E+00 s
 =========================================================================================
 Ghost cells of walls and 400 faces, fills with the extrapolation cut: 0, with the 302 temperature at T_w/2: 0
 400 faces with gas entering, flux evaluations: 49162, largest |w|/a:  2.31614E-01
```

The lines that start with `400 boundary`, and the line `400 faces with gas entering` after the timing, report the gas
that enters the domain through the lower boundary of the case, face 3, of type `400` (extrapolation): the expansion fan
reaches the last part of that boundary, and behind the fan the flow draws gas in through it, up to $|w|/a$ of about
$0.2$. A `400` face lets that gas in with the state of its boundary cell, and the lines change nothing in the run (see
[Boundary Conditions](../user/boundary-conditions.md#euler-boundaries)). `check.py` compares the Mach number at the
probe at the exit with the analytic value.

### Output Files

| File | Content |
|------|---------|
| `OUTPUT/field.tec` | Full 2-D field data (ρ, u, v, p, T, M, etc.) at convergence |
| `OUTPUT/exit.txt` | Probe data (Mach vs. probe location) from the exit boundary |
| `OUTPUT/residual-history.dat` | Convergence history |
| `OUTPUT/diagnostic.tec` | Some diagnostic quantities |

---

## Post-Run Verification

After the run completes, verify that the CFD result matches the analytic prediction:

```bash
python3 check.py
```

This script:

1. Reads `OUTPUT/exit.txt` (probe data from the exit plane)
2. Computes the expected Mach using the Prandtl–Meyer function with $M_1 = 2.0$ and the expansion angle
3. Compares the CFD Mach (`M_CFD`) against the analytic Mach (`M_analytical`)
4. Reports **success** (exit code 0) if relative error < 0.5%, **failure** (exit code 1) otherwise

You can also visualize the full field solution in `OUTPUT/field.tec` using Tecplot or Paraview to confirm that the flow expands smoothly and that the Mach number increases as expected.

---

## Next Steps

Now that you understand the basic workflow:

- **Run other 2-D cases:** See `test/2D/euler/` for additional verification problems.
- **Study the code:** Refer to [Code Structure](../development/structure.md) and [Contributing](../development/contributing.md) for developer documentation.

---

## References

- **Prandtl–Meyer expansion:** [NACA Report 1135](https://ntrs.nasa.gov/citations/19930090976)
- **Riemann solvers:** Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics* (3rd ed.)
- **MUSCL reconstruction:** van Leer, B. "Towards the Ultimate Conservative Difference Scheme. V. A Second-Order Sequel to Godunov's Method" (1979)

---