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
│   ├── bc.txt                ← face-by-face boundary condition data
│   └── decomposition.map     ← block-to-rank assignment (optional, MPI runs)
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

!!! note "Mesh type"
    MOSE takes the type of the mesh from the nodes of every block. The mesh is 3-D if a block has more than one cell
    in $k$, 2-D if not and a block has more than one cell in $j$, 1-D otherwise. A 2-D mesh is an axisymmetric wedge
    when, in every block, the $k = 1$ node plane is the $k = 0$ plane turned about the $x$ axis by the same angle ($x$
    and $r = \sqrt{y^2+z^2}$ kept, nodes on the axis included); it is plane when no block is so turned, and a mesh
    with blocks of both kinds stops the run. In 3-D, type `200` on faces 5 and 6 closes a sector: the last $k$ node
    plane of the block is the first one turned about the $x$ axis by the same angle in every such block, or translated
    (a periodicity with no rotation); otherwise the run stops. Every node is checked against the average move of the
    nodes per cell in $k$, the largest move of a node from the first $k$ node plane to the last over the cells between
    them (the thickness of a 2-D block): a block is turned when its nodes match the turn to within 1 % of that move
    and better than a translation, translated when they match a translation to within 1 % of it and better than a
    turn, and the angles of the blocks agree to 1 %. So the round-off of the coordinates in the file has to stay below
    1 % of a cell in $k$, and a plane slab far from the $x$ axis, which a small turn nearly fits, stays plane. In a
    wedge whose two planes lie symmetric about $z = 0$ ($y$ the same in both), four significant digits are enough at
    any angle, three are not in general. A 2-D block nearer a turn than a translation, or nearer
    a translation than a turn, but not to 1 %, stops the run: write its coordinates with more digits, or make its $k$
    planes parallel for a plane mesh. The set-up summary names the type, as `Mesh type`: 1-D, 2-D plane, 2-D
    axisymmetric wedge with its angle in degrees, 3-D, or 3-D with the angle of its sectors. The periodic images of
    the wall in the wall distance of the turbulence models follow the same rule: the translations of the periodic
    connections (`201`), and of the slabs closed in $k$ by type `200`, are the same when they differ by at most 1 % of
    the thinnest cell across the periodic boundaries, and translations that do not form a lattice to that stop the
    run: the two sides of each periodic boundary must match to that, with enough digits. A channel periodic in $x$
    and turned by 30 degrees about $z$ gives translations that differ by 1.1e-4 of its cells with five digits, by
    8.8e-4 with four, and by 0.15 with two, which stops the run.

---

## Running the Solver

Every simulation can be launched with the following steps:

```bash
export OMP_NUM_THREADS=4
ulimit -s unlimited
./bin/MOSE
```

An MPI build is launched through the usual job launcher, with each rank running
its own OpenMP threads:

```bash
export OMP_NUM_THREADS=10
export OMP_PLACES=cores I_MPI_PIN_DOMAIN=omp
ulimit -s unlimited
mpirun -n 8 ./bin/MOSE
```

!!! warning
    MOSE may require a large stack size. Always run `ulimit -s unlimited` (or set `export KMP_STACKSIZE=100M` for Intel compilers) before launching the solver.

!!! tip
    A basic script like `MOSE.sh` is recommended to ensure the simulation is launched smoothly.

!!! tip "How many ranks, how many threads?"
    The split matters — the same cores can be more than 2× slower in a bad
    layout. As a rule, put at least one MPI rank on every socket and keep each
    rank's threads inside one socket. See
    [Parallel Execution](parallel.md) for the full rule and the measured
    scaling.

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

 =========================================================================================
 Timing
 =========================================================================================
   Iterations                     384
   Solver                          4.75813E+00 s
   Solver per iteration            1.23910E-02 s
   Elapsed                         4.79402E+00 s
 =========================================================================================
```

### Stopping criteria

A run ends after `iter-threshold` iterations, at `time-threshold` (time-accurate
runs) or, in steady runs only, when the density residual falls below
`res-threshold`. The density residual is the L2 norm, over all the cells, of
the change of density in one iteration: it is not normalized, so a suitable
`res-threshold` depends on the density level and on the number of cells. In
time-accurate runs it is reported but does not end the run.

### Timing

Every run closes with the block above. `Solver` is the time spent inside the
iteration loop and is the figure to quote when timing the code; `Elapsed` adds
what happens between iterations, which is dominated by solution output. On more
than one rank the block also reports load imbalance and the fraction of time
blocked on halo exchanges and on collectives.

For a running report rather than a single summary, set `timer-diter` in
`[MOSE-IO]` to the number of iterations between reports:

```ini
[MOSE-IO]
timer-diter = 100
```

```
 MOSE Timing | Level 1 | Iter 100 | 100 iters | wall/iter  9.2450E-01 s | rank min  9.1120E-01 avg  9.1934E-01 | imbalance     0.6 % | exchange wait     3.1 % | collective wait     0.9 %
 MOSE Ranks  | compute/iter max  8.8600E-01 s (rank 3) | min  8.6010E-01 s (rank 7) | mean  8.7412E-01 s | spread     1.4 %
```

Times are the **maximum** over the ranks — the critical path, which is what
sets the time to solution. `imbalance` is how far the slowest rank is above the
mean; a large `compute/iter` spread with a small `imbalance` means the ranks are
unevenly loaded and the fast ones are absorbing it in the exchange wait.
Reporting is off by default (`timer-diter = 0`) and costs a handful of clock
reads per iteration when on.

## Output

After a run completes, the `OUTPUT/` directory contains:

| File | Description |
|------|-------------|
| `field.*` | Solution fields ($\rho$, $u$, $v$, $w$, $p$, $T$, …) |
| `residual-history.dat` | Iteration-by-iteration residual norms |
| `wall.tec` | Wall quantities (skin friction, heat flux) — viscous cases only |
| `diagnostic.*` | Per-cell residuals of $\rho$, $\rho u$, $\rho v$, $\rho w$, $\rho E$, of the soot variables and of the RANS variables, local time step `dt`, shock-sensor flag `beta` — steady runs only; with a single species it needs `sol-variables` to name `thermo` (the default) or `transport`, and the run stops at its first solution output otherwise |
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

By default, MOSE writes primitive variables ($\rho$, $u$, $v$, $w$, $p$, $T$) and, under a RANS model, the density-weighted turbulence variables it advances: `rhoNuTilde` $= \rho\tilde\nu$ (Spalart-Allmaras), `rhoK` and `rhoOmega` ($k$-$\omega$ models), the names of the initial condition ([Initial Conditions](initial-conditions.md#required-variables)). Additional variable groups can be enabled in `[MOSE-IO]` via `sol-variables`:

- `thermo` — specific heat ratio $\gamma$, gas constant $R$
- `transport` — laminar viscosity $\mu_l$ (`mil`), thermal conductivity $k_l$ (`kl`) and, under a RANS model, eddy viscosity $\mu_t$ (`mut`)

### Probes

Point-measurement time histories are configured with probe sections:

```ini
[MOSE-Probes]
probe1 = exit

[exit]
variables = p T u
dtime = 1e-4
```

Each probe writes a text file named after the probe — here `OUTPUT/exit.txt` — with columns for time (or the iteration, with `diter`) and the requested variables.

Accepted variables: `rho(n)` (density of species `n`), `u`, `v`, `w`, `p`, `T`, `M` (Mach number) and, with a RANS model, `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde` (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model). An unknown name stops the run at set-up.

A new run replaces the probe files; a restart appends to them, and creates any that are missing (for example a probe added before the restart).

---
