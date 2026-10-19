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

!!! note "Exit code"
    A run that MOSE stops on an error — at loading, at set-up, or because the solution diverges («Integration failed at
    b, i, j, k», then «NaN, infinity or p<0 detected») — ends with a non-zero exit code, so that scripts and job
    schedulers see the failure: 1 where MOSE ends with `stop 1` or `error stop 1`; otherwise that of `ERROR STOP` (128
    with ifx, 1 with gfortran) or, when a rank stops every rank of an MPI run with `MPI_ABORT`, the one the MPI launcher
    gives for `MPI_ABORT` with error code 1. With MPI, an error that every rank finds (for example a decomposition map
    that does not fit, an `input.ini` that cannot be read, a line of `bc.txt` that cannot be read, a time table that
    cannot be read or whose times decrease, a mesh whose blocks do not make one type) stops every rank together, and so
    do a value of `input.ini` read again out of its rule, an initial condition that cannot be read, also when one rank
    alone finds it, and an error that rank 0 alone finds at set-up (in the thermodynamic, transport or chemistry data,
    in the rules of `input.ini`, in the boundary conditions), which rank 0 writes and then passes to every rank: each rank writes its
    message, the ranks finalize MPI and end with exit code 1. A divergence is caught on any value that is not finite
    (NaN, infinity) or on a negative pressure, also in RELEASE builds with fast math; a value of the initial condition
    that is not finite stops the run at set-up, in every build ([Initial Conditions](initial-conditions.md#required-variables)).

    Known limit: errors raised inside the ORION library, which reads and writes the solution files, still end with
    code 0, because ORION stops with `STOP`. For example, a binary Tecplot file (`ic-format = tecplot binary`) in a
    build without TecIO ends with «You can not read in binary format without compiling against TecIO» and exit code
    0. MOSE cannot refuse that input beforehand: whether TecIO is compiled in is known only inside ORION (its `TECIO`
    macro). Check the end of the log, not only the exit code, when a run reads or writes binary files.

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
runs) or, in steady runs only, when the density residual has fallen by the
factor `res-threshold` from the first iteration of the run (of each grid level
of a multigrid sequence). The density residual is the L2 norm, over all the
cells, of the change of density in one iteration; `residual-history.dat` and
the screen report it as it is, while the convergence test takes it relative to
the first iteration, so that the density level and the number of cells do not
set the threshold. An input written for the absolute test of earlier versions
keeps its stop with `res-threshold` divided by the density residual of the
first iteration, R₁: Rₙ ≤ (ε/R₁) R₁ is Rₙ ≤ ε. R₁ is the first row of
`residual-history.dat` when `res-diter` is 1 (ten digits) and the
`MOSE_RES_REFERENCE` of the vtk solutions (all its bits); this holds on one
grid level, since a sequence of grids measures each level from its own first
residual. A residual that rises before it falls is still measured
from the first iteration: the run goes on until the residual is
`res-threshold` times that of the first iteration, not times its peak, so it
ends later than a test relative to the peak would, never earlier. A restart
from a solution keeps the residual of the first iteration of the
uninterrupted run, with all its bits in every format (see
[Restart](initial-conditions.md#restart)). In
time-accurate runs the residual is reported but does not end the run.

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
| `diagnostic.*` | Per-cell residuals of $\rho$, $\rho u$, $\rho v$, $\rho w$, $\rho E$, of the soot variables and of the RANS variables (those of the last stage of the iteration written; the wall output, written just before, leaves them alone), local time step `dt`, shock-sensor flag `beta` — steady runs only; with a single species it needs `sol-variables` to name `thermo` (the default) or `transport`, and the run stops at its first solution output otherwise |
| `<probe>.txt` | Time-history data at probe locations called <probe>|

### Input/Output Formats

Controlled by `sol-format` in `[MOSE-IO]`:

| Value | Format |
|-------|--------|
| `tecplot ascii` | Tecplot ASCII (`.tec`) — default |
| `tecplot binary` | Tecplot binary (`.szplt`) — requires TecIO |
| `vtk ascii` | VTK ASCII (`.vts`) |
| `vtk raw` | VTK binary (`.vts`): the data of the field, wall and diagnostic outputs in raw binary, appended to each file |

With the vtk formats each output is a multi-block file, `OUTPUT/field.vtm` for instance, that lists its blocks relative to its own folder: `vtk/fieldBlock1.vts`, ... In every format a restart can reproduce the uninterrupted run exactly (see [Restart](initial-conditions.md#restart)): the binary formats store the bits of every value, the text formats write 17 significant digits. A field solution also holds the counters of a restart: field data in vtk, data set auxiliary data in Tecplot, which a Tecplot ASCII file writes as `DATASETAUXDATA` records after its zones; a script that reads every number after a zone header must stop at them.

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

A probe records one cell: the cell that `index-position` gives (block, i, j, k), or, when `index-position` is left out, the cell whose centre, the mean of its eight nodes, is the nearest to `position`.

Accepted variables: `rho(n)` (density of species `n`), `u`, `v`, `w`, `p`, `T`, `M` (Mach number) and, with a RANS model, `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde` (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model). An unknown name stops the run at set-up.

A new run replaces the probe files; a restart appends to them, and creates any that are missing (for example a probe added before the restart). Time-based rows fall on the multiples of `dtime`, from the first. A restart that goes on in time writes its next time-based row at the first multiple of `dtime` its time has not reached; a restart that goes on from the counters of a solution, in any format, first drops the rows written after that solution (see [Restart](initial-conditions.md#restart)).

---
