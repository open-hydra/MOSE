# Running MOSE

This page describes the complete workflow for running a MOSE simulation: preparing the case directory, defining the needed files, launching the solver, and inspecting the output.

All configuration is driven by:  
- `input.ini` — the main INI file that defines the simulation parameters, including the numerical settings and physical models.  
- `INPUT/` — a directory containing the supporting data files for the simulation, such as thermodynamic properties, initial conditions, and boundary conditions.

!!! tip
    All files within `INPUT/` can be built using the MOSE pre-processing tools provided by ATLAS. ATLAS is a very powerful tool, but it is not required to run MOSE. All files can be created manually or with custom scripts as long as they follow the expected format.

For the detailed format of thermodynamic, initial and boundary conditions see the dedicated reference pages:

- [Gas mixture data (FLINT library)](https://github.com/MarcoGrossi92/FLINT)
- [Initial Conditions](initial-conditions.md)
- [Boundary Conditions](boundary-conditions.md)

---

## Simulation Workflow

A typical MOSE run follows the following phases (using ATLAS):

1. **Setup** — define the mesh, thermodynamic data, initial and boundary conditions
2. **Preprocess** — run ATLAS pre-processing tools to generate the necessary input files in `INPUT/`
3. **Solve** — advance the solution in time (explicit time integration)
4. **Postprocess** — write solution files, wall data, probes, and diagnostics

All configuration is driven by the INI file (`input.ini`).

---

## Case Directory Structure

Every MOSE case follows a standard layout:

```
my_case/
├── input.ini                 ← case configuration
├── MOSE.sh                   ← convenience run script (optional)
├── INPUT/
│   ├── phase.txt             ← species and mixture definition
│   ├── thermo.dat            ← thermodynamic property tables
│   ├── ...                   ← other mixture files (chemistry, transport, etc.)
│   ├── ic.*                  ← initial condition grid + flow field
│   └── bc.txt                ← face-by-face boundary condition data
└── OUTPUT/                   ← created at runtime
    ├── field.*               ← solution fields (when sol-overwrite = true)
    ├── field1.*              ← solution snapshots (when sol-overwrite = false)
    ├── field2.*
    ├── ...
    ├── wall.*                ← wall quantities (viscous cases only)
    ├── residual-history.dat  ← residual norm vs. iteration
    └── <name>.txt            ← time-history at probe locations
```

!!! warning
    MOSE reads the grid geometry directly from the initial condition file (the `x`, `y`, `z` node coordinates). There is no separate mesh file.

---

## Running the Solver

Every simulation can be launched with the following steps:

```bash
export OMP_NUM_THREADS=4
ulimit -s unlimited
./bin/MOSE
```

!!! warning
    MOSE may require a large stack size. Always run `ulimit -s unlimited` (or set `export KMP_STACKSIZE=100M` for Intel compilers) before launching the solver.

!!! tip
    A basic script like `MOSE.sh` is recommended to ensure the simulation is launched smoothly.

---

## Logging and Diagnostics

MOSE prints diagnostic information to the console during loading and iteration. This includes:  
- Loading status of input files (thermodynamics, transport, chemistry, initial and boundary conditions)
- Summary of the case set-up (number of blocks, cells, boundary faces, types of boundary conditions, physical model details, numerical scheme choices)
- Iteration count and residuals during the run

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

MOSE | Iter =      100 | Global iter =      100 | Density residual = 0.125027E-03
MOSE | Iter =      200 | Global iter =      200 | Density residual = 0.231340E-06
MOSE | Iter =      300 | Global iter =      300 | Density residual = 0.221740E-08
MOSE | Iter =      384 | Global iter =      384 | Density residual = 0.998459E-10

   Time of operation was   7.9302149999999988E-002 min
```

## Output

After a run completes, the `OUTPUT/` directory contains:

| File | Description |
|------|-------------|
| `field.*` | Solution fields ($\rho$, $u$, $v$, $w$, $p$, $T$, …) |
| `residual-history.dat` | Iteration-by-iteration residual norms |
| `wall.tec` | Wall quantities (skin friction, heat flux) — viscous cases only |
| `<probe>.txt` | Time-history data at probe locations called <probe>|

### Input/Output Formats

Controlled by `sol-format` in `[MOSE-IO]`:

| Value | Format |
|-------|--------|
| `tecplot ascii` | Tecplot ASCII (`.tec`) — default |
| `tecplot binary` | Tecplot binary (`.szplt`) — requires TecIO |
| `vtk ascii` | VTK ASCII (`.vts`) |
| `vtk raw` | VTK binary (`.vts`) |

### Output Variables

By default, MOSE writes primitive variables ($\rho$, $u$, $v$, $w$, $p$, $T$). Additional variable groups can be enabled in `[MOSE-IO]` via `sol-variables`:

- `thermo` — specific heat ratio $\gamma$, gas constant $R$
- `transport` — laminar viscosity $\mu_l$, thermal conductivity $k_l$, turbulent viscosity $\mu_t$

### Probes

Point-measurement time histories are configured with probe sections:

```ini
[MOSE-Probes]
probe1 = exit

[exit]
variables = p T u
dtime = 1e-4
```

Each probe writes a text file `OUTPUT/probe_exit.txt` with columns for time and the requested variables.

---
