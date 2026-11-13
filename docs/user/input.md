# Input File

MOSE is configured through a single **INI-format** file called `input.ini`, located in the case root directory.

The file is organized into parameter blocks — each block controls a specific aspect of the simulation (flow model, numerics, output, turbulence, chemistry, etc.). Parameters not specified take their default values.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Example:

```ini
[MOSE-Parameters]
res-threshold = 1e-8
iter-threshold = 500000

[MOSE-Physics]
equations = navier-stokes
turbulence = SST

[MOSE-Numerics]
cfl = 0.9
time-scheme = RK3
space-reconstruction = MUSCL
```

!!! note "Comments"
    A line that starts with `!`, `;` or `#` is a comment. After a value only `;` starts a comment
    (`cfl = 0.9 ; CFL number`): a `#` or a `!` after a value is part of the value. A number followed by one
    does not read and stops the run; a name followed by one stops the run when its parameter takes one of a
    list of names, and is otherwise taken with the comment in it (a probe named so finds no section of that name,
    and the run stops).

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Incorrect names are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

!!! warning "Keys that are not read"
    MOSE warns at the start for every line `key = value` of `input.ini` whose key the reader does not give in its section, with its line: such a key is an error in the input, and the run goes on with the default of its parameter in place of the value written. A parameter of the registry written with no value is named here, and stops the run (below). Comment lines (lines that start with `!`, `;` or `#`) are skipped, also when they hold `=`.

!!! warning "An input.ini that cannot be read stops the run at the start"
    If `input.ini` cannot be opened (missing, or without read permission), holds no line that is not blank, or the INI reader returns an error while it loads the file (a file with no section, a line it cannot parse, which the reader names on its error output), the run stops at the start with the reason and exit code 1; with MPI every rank stops. Without the stop the parameters would keep their defaults in silence. For a reading during the run, see the note «Errors of a reading during the run» below.

!!! warning "Values outside their rules stop the run at the start"
    A value outside the rule of its parameter (column Allowed of the [Parameter Registry](registry.md): a bound such as `> 0`, an interval such as `in (0, 1)`, where `(` and `)` leave the end out and `[` and `]` take it in, each held by every value of a list of numbers; a list of names; `logical`) stops the run with a message at the start. MOSE reads `input.ini` again during the run, every `ini-diter` iterations (10000 by default), to take the values changed meanwhile; what an error of such a reading does is in the note «Errors of a reading during the run» below. In an MPI run rank 0 reads `input.ini`, at the start and at every reading, and passes its text to every rank, so that every rank takes the same values. With one working folder per rank, only the `input.ini` of rank 0 is read. The parameters that `input.ini` leaves out keep their values, the `stationary-face-<n>` keys of the rotating frame too; what a reading changes in the rotating frame, and what it does not take, is in [Rotating frame](#rotating-frame).

!!! warning "Values that do not read stop the run at the start"
    A value that does not read as the type of its parameter in the [Parameter Registry](registry.md) stops the run at the start with a message that gives the file, the section, the key and the value as written (a value of `levels` of `[MOSE-Multigrid]` that is not an integer too; for a reading during the run, see the note below): a real must be a finite number (not a word, `NaN`, `0,5`, or `1e999`, which gfortran reads as an infinity), an integer has no decimal point and no exponent (not `2.0` or `1e3`), a logical starts with `T` or `F` after an optional period (not `yes` or `1`; a word that starts with one of them reads as true or false), a list holds no more values than its parameter, and a key has a value. A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity) stops the run at the start, with a message that names the key.

!!! warning "Errors of a reading during the run"
    When MOSE reads `input.ini` again during the run, an error of the reading (a file that cannot be read, a value that does not read as the type of its parameter, is out of its rule or is not finite, a value of a key fixed for the run other than the start read, a value of `levels` that is not an integer or is not the one the start read or more `probe<n>` keys than at the start, a change of `rotational-frame` or of the stationary faces of the rotating frame, a value of `axis` or `origin` that is not finite or does not read as numbers, a `stationary-face-<n>` value that does not read, an axis of zero length, residual smoothing switched on in a run started without it, values that a check of the set-up between keys refuses) does what `ini-error` of `[MOSE-IO]` says. With `warn`, the default, the reading is not taken: MOSE writes a warning with the iteration and the error, which names the key and its value as written, and the run goes on with every value it had before the reading, those of the keys of the same reading that were right too. With `stop` the run stops with the message and exit code 1, as at the start; with MPI every rank stops. The value of `ini-error` in force when the reading starts decides. At the start every error stops the run, whatever `ini-error` says.

!!! warning "What a reading during the run changes"
    The column During the run of the [Parameter Registry](registry.md) gives what a reading does with each key. A key *taken*, such as `cfl`, `vnn` or the intervals of the outputs, takes the value written. A key *fixed* for the run, such as `equations`, `turbulence`, `time-scheme`, `time-accurate`, `space-reconstruction` or `sol-format`, is one from which the set-up makes a choice that a reading cannot make again: the run keeps the value that the start read, and a reading that writes another one is an error of the reading, also when the value names what the run does (`shock-detector = Tramel` in a run that started without the key, whose sensor comes from `riemann-solver = HLLC+Tramel`); the message gives the value of the start, or says that the start left the key out. The grid levels, the probes and the stationary faces of the rotating frame are those of the start too: another `levels`, more `probe<n>` keys, or a `probe<n>` or a `stationary-face-<n>` that names another probe or face, is an error of the reading. A key of the *start only* (`newrun`, `ic-format`, `decomposition-map`) serves the start of the run, and a reading leaves it as it was, whatever value of its rule the file writes (a value that does not read, or is out of its rule, is an error of the reading, as for every key): `newrun = false` written during the run, to prepare a restart, does nothing to the run. A key that the file leaves out keeps its value. The checks between keys that the set-up makes are made again on the values read: `riemann-options-Mco = 0` with a low-Mach Riemann solver (`AUSM+M`, `HLLC+Chen`, `LMRoe`), `schmidt` = 0 with more than one species in a viscous run without `INPUT/diffusion.dat`, or, with the multirate integrator, a `cfl` above its bound or `cfl-rise-threshold` above 0 (section `[MOSE-Multirate]` below), is an error of the reading as it is at the start. Residual smoothing is on with `irs = true` or `irs-beta` above 0, at the start and at every reading ([Models and methods](#models-and-methods)): a reading that switches it on in a run started without it is an error of the reading, since its work arrays exist only when it was on at the start.

---

## Sections

| Section | Description | Reference |
|---------|-------------|-----------|
| `[MOSE-Parameters]` | Simulation control: convergence thresholds, restart | [→](registry.md#mose-parameters) |
| `[MOSE-Physics]` | Equations, turbulence, chemistry, soot, rotating frame | [→](registry.md#mose-physics) |
| `[MOSE-IO]` | Output formats, frequencies, and variables | [→](registry.md#mose-io) |
| `[MOSE-Probes]` | Point probes for time-history recording | [→](registry.md#mose-probes) |
| `[MOSE-Numerics]` | Spatial/temporal discretization, Riemann solver, CFL | [→](registry.md#mose-numerics) |
| `[MOSE-Multigrid]` | Multigrid acceleration settings | [→](registry.md#mose-multigrid) |
| `[MOSE-Multirate]` | Multirate (local time stepping) integration of time-accurate runs | [→](registry.md#mose-multirate) |
| `[MOSE-Chemistry]` | ODE solver and tolerances for finite-rate chemistry | [→](registry.md#mose-chemistry) |
| `[MOSE-Turbulence]` | Turbulent Prandtl/Schmidt numbers | [→](registry.md#mose-turbulence) |
| `[MOSE-Rotating-Frame]` | Rotating reference frame parameters | [→](registry.md#mose-rotating-frame) |

`levels` in `[MOSE-Multigrid]` is the number of grid levels, 1 (a single grid) by default; it is not
in the registry. MOSE reads it at set-up as an integer: a value that is not an integer, or no value,
stops the run, and a value below 1 counts as 1. With `levels = n` the number of cells of every block
must be a multiple of 2^(n−1) in i, and in j and k unless it is 1, or the run stops at set-up.
`level1-iter`, `level2-iter`, ... give the iterations of each level. They are read again with `input.ini`
during the run: lowering one below the iterations already done ends that level at the iteration after that reading.

---

## Multirate time stepping

In a time-accurate run, `[MOSE-Multirate] enabled = .true.` lets each cell advance with a sub-step close to its own
stability limit instead of the global minimum. The method is in [Time Integration](../theory/multirate.md#multirate-time-stepping);
its checks are in [Multirate Verification](../vv/multirate.md).

```ini
[MOSE-Numerics]
time-scheme = RK2
time-accurate = true
integration-variables = cons

[MOSE-Multirate]
enabled = .true.
max-level = 10
safety = 2.0
levels = recompute
recompute-every = 1
```

The defaults suit production runs. Set `cfl` and `vnn` within the limits [below](#limits-on-cfl-and-vnn).

### Keys

One iteration is one macro step $H$, and the residual printed is the change of the state over it.

| Key | Default | Meaning |
|-----|---------|---------|
| `max-level` | 10 | Deepest level allowed. The finest cells advance with $H/2^L$, $L \le$ `max-level`; fewer levels are used when the local steps span fewer octaves. |
| `safety` | 2.0 | Safety factor $\sigma \ge 1$: every cell advances with a sub-step of at most $\Delta t^\text{cap}/\sigma$, $\Delta t^\text{cap}$ being its local CFL/VNN step. With `max-level = 0`, $\sigma = 1$ gives the single-rate step. |
| `macro-step-target` | 0 | Caps $H$, in seconds; 0 = no cap. |
| `levels` | `recompute` | `recompute` rebuilds the level map for the $H$ taken: use it for production. `fixed` keeps the map of the first macro step, so that the error of a run measures the order of the scheme alone: use it for order studies only. |
| `recompute-every` | 1 | Macro steps between two rebuilds of the map (`levels = recompute`). |
| `fixed-macro-step` | 0 | A constant $H$, as a multiple of the smallest local step of the first macro step, for verification runs; 0 = adaptive. |
| `nonlocal-cfl` | `redo` | Check of condition (4) of Gnedin et al.: `off`, `count` or `redo` ([below](#the-check-nonlocal-cfl)). |
| `interface-margin` | `.false.` | Keeps the coarse step at a face between levels within $\Delta t^\text{cap}/\kappa$ of the finer cell. Always on with `redo`. |
| `max-attempts` | 21 | With `redo`: the most attempts of one macro step. |
| `force-redo` | `.false.` | Test device: rejects the first attempt of every macro step; the results must not change. |
| `conservative-correction` | `.true.` | `.false.` removes the flux correction between levels, as a positive control of the conservation tests. |
| `trace` | `.false.` | Prints the stage sequence, level map and halo lists of the first macro step (lines `MR-TRACE`, `MR-LEVEL`, `MR-HALO`, `MR-VISC`). |

`fixed-macro-step` is for verification runs. With `levels = fixed` and an $H$ above the step of the map, the finest
cells start beyond their limit, and the run warns once.

### The check `nonlocal-cfl`

- `off`: no check.
- `count`: the violations are counted and printed at the end of a run that has some. Nothing else changes.
- `redo`: a macro step that violates the condition is redone from its start, with finer levels (at most `max-level`) or
  with $H$ halved when a cell needs more levels than that.

With `redo` a macro step in which a cell or a face takes a state that is not physical (not finite, a total density that is
not positive, a negative pressure, a temperature outside the tables of the gas) is redone with $H$ halved on the same levels. A
face that is not physical at the start of the macro step stops the run, as in the single-rate step.

The run stops, with the rejections by cause, when $H$ falls below $2^{-10}$ of its first attempt or when `max-attempts`
attempts have been rejected. With `count`, `off`, `fixed-macro-step` or `levels = fixed`, `max-attempts` gives a warning and does nothing.

At the start of every macro step a local step that is zero, negative, infinite or not a number stops the run, with the
number of such cells and the first one, with its local step and primitive variables.

### Limits on `cfl` and `vnn`

With `max-level` > 0, `nonlocal-cfl = redo` and `levels = recompute` and no `fixed-macro-step`, the run stops at set-up
when `cfl` exceeds 1 or, in a Navier–Stokes run, `vnn` exceeds 0.5 (the stability limits of
[Time step](#time-step-cfl-and-vnn)). Otherwise the set-up checks neither value.

### What the run prints

Every `shell-diter` iterations the run prints:

- the macro step $H$, the depth $L$ of its map and the cells per level;
- the wall time per level, and on several processes the part spent in the exchanges (these are the times of process 0);
- $r_\text{post}$ and $r_\text{nb}$, the largest ratios of the step each cell took to the step admitted on the state it
  reached, and to that of a face neighbour; a line with the number of cells above 1 for each.

When a macro step is rejected, more lines say why. At the end of the run the largest ratios are printed, with the number of
macro steps that had a cell above 1. A run in which none of this happened prints no such lines. A run on several processes
ends with the minimum, mean and maximum over the processes of the time of the stages and of their exchanges, per level
(see [Timing](using.md#timing)). These lines only report: they change neither the steps nor the results.

### Output times, probes and restarts

- The macro step lands exactly on the output instants $k\,$`sol-dtime` and on `time-threshold`: when one falls within
  $H$, $H$ is cut to end there, and the time-based solution is written at that macro step. A cut macro step counts as an iteration.
- Time-based probe rows come at the end of the macro steps. A probe whose `dtime` is below $H$ writes one row per macro step.
- A restart (`newrun = false`) from a `vtk raw` or `tecplot ascii` solution goes on bit for bit as the uninterrupted run
  with `levels = recompute`. With `recompute-every` = $k$ > 1 it does so when it falls on a multiple of $k$. `levels = fixed` and
  `fixed-macro-step` cannot restart, and stop the run at set-up.
- The multirate integrator writes no file of its own.

### Reading `input.ini` during the run

The multirate keys are fixed for the run (column *During the run* of the [Parameter Registry](registry.md)): another
value written during the run, `enabled` included, is an error of the reading, which `ini-error` decides. `time-accurate` is
fixed too. The checks of `cfl`, `vnn` and `cfl-rise-threshold` are made again at every reading: a `cfl` or `vnn` above its
limit, or `cfl-rise-threshold` > 0, is an error of the reading. A `cfl` within the limit is taken, and the interface
margin follows it from the next macro step.

### Scope

The integrator runs time-accurate Heun (`time-scheme = RK2`) with `integration-variables = cons`, for the Euler and
Navier–Stokes equations, laminar or with:

- `SA`, `SA-R` or `SST`, with `k-energy-coupling` and either `sst-production`, and `sst-dilatation-dissipation`;
- with `max-level = 0` also `SA-RC` and `SST-RC`.

It runs on several blocks (the 2-D axisymmetric wedge included), with an SRM grain (`502`) in mode 2, and on several MPI
processes with OpenMP threads, with the same results as on one process. With `SA`, `SA-R` or `SST` the local step
includes a bound of the stiffness of the turbulence sources.

Outside its scope the run stops at set-up or at its first macro step, with a message `[MOSE-Multirate] <item>` and a
non-zero exit code:

- other time schemes, `prim` or `prec`, steady runs, `cfl-rise-threshold` > 0, multigrid sequences;
- implicit residual smoothing, a shock detector (also implied by `HLLC+Tramel`, `HLLC+Chen`, `HLLE++`, `AUSM+M`);
- chemistry, soot, the rotating frame;
- `Wilcox2006`, `SAcomp`, and with `max-level` > 0 `SA-RC` and `SST-RC`;
- chimera (`102`), coupled walls (`103`), an SRM grain in mode 1, gas–surface interaction (`503`–`505`), rotational
  periodicity (`200`) in 3-D, and boundary types without a multirate boundary flux;
- a `cfl` or `vnn` above the limits [above](#limits-on-cfl-and-vnn);
- a restart with `levels = fixed` or a `fixed-macro-step`.

---

## Time step: `cfl` and `vnn`

Each cell takes the largest explicit step that keeps it stable, from one rule
([theory](../theory/time-integration.md#stability-conditions)):

- `cfl` is the Courant number of the **sum over the directions** of the convective
  spectral radii, $|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$, measured on the
  width of the cell across its faces;
- `vnn` is the diffusion number of the **largest diffusivity** of the equations, also summed
  over the directions: in air the energy one, $\gamma\nu/\text{Pr} \approx 1.94\,\nu$;
  with Spalart–Allmaras also $(\nu + \tilde\nu)/\sigma$;
- the convective and the diffusive rates **add up**: a cell where both matter gets a step
  smaller than either limit alone;
- in time-accurate RANS runs the stiffness of the turbulence sources adds up as well, with no
  parameter: near walls, where the destruction of $k$, $\omega$ or $\tilde\nu$ is fast, the step can
  be smaller than the convective and diffusive limits
  ([theory](../theory/turbulence-numerics.md#time-accurate-runs-explicit-sources)). With SST or Wilcox 2006,
  a time-accurate start far from the wall value of $\omega$ (100 times or more below it, in wall cells
  where diffusion sets the step) needs the CFL ramp `cfl-rise-threshold`: without it the first step
  can make $\omega$ negative.

Stability limits of the linearised schemes (MUSCL with upwind fluxes, central diffusion);
the additive rule stays stable for any `cfl` and `vnn` up to these values:

| `time-scheme` | `cfl` ≤ | `vnn` ≤ |
|---|---|---|
| `RK2` (Heun) | 1 | 0.5 |
| `RK3` | 1.17 | 0.62 |
| `euler` | 1 with `first-order` reconstruction (unstable with `MUSCL`) | 0.5 |

These are linear limits. With the limiter, positivity (a maximum principle; the scheme is not TVD in more than one
dimension) is guaranteed
only up to `cfl = 0.5` and `vnn = 0.5`, with Heun and with RK3
([theory](../theory/time-integration.md#local-time-step)): the transport keeps the turbulence variables, and with
`integration-variables = cons` (the default) the species, non-negative on any mesh, with every limiter but `vanalbada` and
without implicit residual smoothing (`irs`, which mixes the updates of neighbouring cells), when the Courant number of the
flow leaving each cell is at most 1/2, which `cfl = 0.5` gives in smooth flow; the bound is lower in the first cell next
to a wall, and positivity is not guaranteed next to faces whose
ghost cells extrapolate the cells (types `408` and `502` in its mode 1). The `vnn` limits are sharp: a run
a little beyond them blows up. Suggested values:

- `RK3`: `cfl = 0.8` (above the guaranteed 0.5) and `vnn = 0.5`;
- `RK2` (Heun): `cfl = 0.5` with moving shocks, up to `cfl = 0.7` where the flow is smooth, and
  `vnn = 0.4`; never `cfl` ≥ 0.95 or `vnn` ≥ 0.5.

The defaults are `cfl = 0.5` and `vnn = 0.3`; in a viscous time-accurate run raise `vnn`
within the limits of the table before `cfl`. With the multirate integrator (`[MOSE-Multirate] max-level` > 0 and
`nonlocal-cfl = redo`, the default, with `levels = recompute` and no `fixed-macro-step`) the set-up and every reading of `input.ini` reject a `cfl` above 1 and, in a
Navier–Stokes run, a `vnn` above 0.5, the linear limits of Heun in the table
([Multirate time stepping](#multirate-time-stepping)); the values suggested above stay suggestions.

---

## Models and methods

The [Theory Guide](../theory/index.md) names the methods and their parameters; these are the keys that select them
(values and defaults in the [Parameter Registry](registry.md)).

| Method (Theory Guide) | Key | Section |
|---|---|---|
| Time scheme: forward Euler, Heun, SSP-RK3 | `time-scheme = euler`, `RK2`, `RK3` (required) | `[MOSE-Numerics]` |
| Conservative, primitive or preconditioned update | `integration-variables = cons`, `prim`, `prec` | `[MOSE-Numerics]` |
| Preconditioning: $U_\text{ref}$, target Mach $M_t$, floor $\varepsilon$, factor $f$, $U_{\min}$, smoothing passes $n$ | `preconditioning-Uref`, `-Mach`, `-eps-min`, `-Ur-factor`, `-Ur-min`, `-Ur-smooth` | `[MOSE-Numerics]` |
| CFL, VNN; CFL ramp over $N$ iterations | `cfl`, `vnn`; `cfl-rise-threshold` | `[MOSE-Numerics]` |
| Reconstruction (MUSCL, first order) and limiter | `space-reconstruction`, `flux-limiter` | `[MOSE-Numerics]` |
| Riemann solver; low-Mach cutoff $M_{co}$ | `riemann-solver`; `riemann-options-Mco` | `[MOSE-Numerics]` |
| Shock sensor (Tramel, Chen) | `shock-detector = Tramel`, `Chen` | `[MOSE-Numerics]` |
| Face metric of the viscous fluxes (stencil, mean) | `diffusive-metric = stencil`, `mean` | `[MOSE-Numerics]` |
| Residual smoothing, coefficient $\varepsilon$; flow or all residuals | `irs`, `irs-beta`; `irs-variables = flow`, `all` | `[MOSE-Numerics]` |
| Grid sequencing: levels MGL, iterations per level | `levels`; `level<n>-iter` | `[MOSE-Multigrid]` |
| Multirate integrator | see [Multirate time stepping](#multirate-time-stepping) | `[MOSE-Multirate]` |
| SA, SA-R, SA-RC, SST, SST-RC, Wilcox 2006 | `turbulence = SA`, `SA-R`, `SA-RC`, `SST`, `SST-RC`, `Wilcox2006` | `[MOSE-Physics]` |
| SSG/LRR: SSGLRR-RSM-w2012, -w2012-SD, -w2019, -w2019-SD | `turbulence = SSGLRR`, `SSGLRR-SD`, `SSGLRR-w2019`, `SSGLRR-w2019-SD` | `[MOSE-Physics]` |
| SST production (incompressible, compressible); energy coupling | `sst-production`; `k-energy-coupling` | `[MOSE-Turbulence]` |
| Wall condition of $\omega$ (practical, asymptotic); point-implicit factor | `omega-wall-bc`; `point-implicit` | `[MOSE-Turbulence]` |
| Turbulent Prandtl and Schmidt numbers $Pr_t$, $Sc_t$ | `Prt`, `Sct` | `[MOSE-Turbulence]` |
| Frozen, finite-rate or equilibrium chemistry | `chemistry = frozen`, `finite-rate`, `equilibrium` | `[MOSE-Physics]` |
| Chemistry integrator | `ode-solver`, `ode-max-steps` | `[MOSE-Chemistry]` |
| Rotating frame: $\boldsymbol\Omega$, axis, origin $\boldsymbol{x}_o$ | `rotational-frame = rigid-body`; `omega`, `axis`, `origin` | `[MOSE-Physics]`, `[MOSE-Rotating-Frame]` |

Notes on these keys:

- `time-scheme` is required: a run without it stops at set-up with `[ERROR] Required parameter not set: time-scheme`.
  The Default column of the registry lists `euler` for it, a placeholder that a required key never uses.
- `integration-variables = prec` is for steady runs: a time-accurate run with it stops at set-up. `HLLC-PC` needs
  `prec`, and stops at set-up without it. `prec` with another flux keeps a dissipation scaled with the speed of sound:
  with `HLLC` the low-Mach cylinder at Mach 0.01 diverges. The modes of $U_r$ take, in this order, `preconditioning-Uref`
  > 0, then `preconditioning-Mach` > 0, then neither (Weiss–Smith); take `preconditioning-eps-min` about a tenth of the
  Mach number of the flow ([Time Integration](../theory/time-integration.md#reference-velocity-of-the-preconditioning)).
- `space-reconstruction = MUSCL` without `flux-limiter` takes van Leer, with a warning; a `flux-limiter` given with
  `first-order` is ignored, with a warning.
- Without `shock-detector`, `HLLC+Tramel` and `HLLE++` take the Tramel sensor, `HLLC+Chen` and `AUSM+M` that of Chen,
  and the other Riemann solvers none.
- Residual smoothing is on with `irs = true` or `irs-beta` > 0, and MOSE refuses it at set-up in a time-accurate run.
  Its arrays are allocated at the start only: a reading of `input.ini` that switches it on in a run started without it
  is an error of the reading. On a run started with it, `irs-beta` and `irs-variables` can change during the run, and
  it goes off when the values in force give neither `irs = true` nor `irs-beta` > 0: `irs-beta = 0` switches it off
  unless `irs = true` is in force, and `irs = false` leaves it on while `irs-beta` > 0.
- In time-accurate runs `point-implicit` is ignored; when `input.ini` sets `point-implicit = .true.`, a warning at
  start-up recalls it. Time-accurate runs with `SA-RC` or `SST-RC` print a warning: the correction is for steady flows.
- `ode-solver`: the default `H-radau5` is the Radau IIA method of order 5; `cvode` needs OSLO built with SUNDIALS, and
  `H-dopri5` is explicit, for non-stiff mechanisms only. `ode-max-steps` bounds the integrators of Hairer and of FATODE;
  `cvode` and `dodesol` do not take it.

### Rotating frame

A wall at rest in the laboratory (a casing or a hub that does not turn) is listed in `[MOSE-Rotating-Frame]` as `stationary-face-<n> = <block> <face>`, numbered from 1; the walls `301` and `302` of that block face take it. The other records of that face may be of the types that never read whether their face is stationary, whose boundary condition is then the same as without the key: block connections, periodic and chimera faces (`101`, `102`, `200`, `201`), extrapolation (`400`), the ambient boundary (`407`), the mapped boundary (`410`) and unused faces (`0`); so a casing face that is a wall in part and a block connection in part can be stationary. A stationary face that covers a type that would take the velocity of the wall at rest in part only (the inviscid flux of a symmetry `300`, of the inlets and outlets `401`–`406`, `408` and `420`, and of the walls `103`, `104` and `502`–`505`), that has no wall `301` or `302`, or that names a face with no record in `bc.txt` stops the run at set-up.

The axis and angular velocity magnitude/direction must be specified at initialization. When `input.ini` is read again during the run (every `ini-diter` iterations), the `omega`, `axis` and `origin` read again replace those of the set-up in every term at once (the log of the run prints them): the Coriolis and centrifugal sources, the velocity of the stationary walls, of the inlets and of the far field, and the power of the torque file. The frame is taken at constant $\boldsymbol{\Omega}$ between two readings: the force $-\rho\,\dot{\boldsymbol{\Omega}} \times \boldsymbol{r}$ of a frame that speeds up is not modelled. A change of `rotational-frame` or of the stationary faces, which changes the equations or the boundary conditions, is an error of the reading: a `stationary-face-<n>` key that names another face, or one that the set-up does not have. A key that the file read again leaves out keeps its face, as every parameter left out keeps its value. The keys are counted from 1 up to the first one missing, at the start and at every reading: a key after a gap of the numbering is not read. A value of `axis` or `origin` read again that is not finite (NaN or infinity) or does not read as numbers, a `stationary-face-<n>` value that does not read, or an axis of zero length, is an error of the reading too. With `ini-error = warn`, the default, a reading with an error is not taken and the frame goes on as it was; with `ini-error = stop` the run stops (note «Errors of a reading during the run» above).

---

## Blocks without chemistry: `exclude-blocks`

`exclude-blocks` in `[MOSE-Chemistry]` lists the blocks that take no chemistry. With
`chemistry = finite-rate` or `chemistry = equilibrium` the chemistry step skips them, on every grid
level of a multigrid sequence, so their composition changes by transport only; the set-up summary
lists them as inactive blocks. The value is `none`, the default, or the numbers of the blocks
separated by blanks, for example `exclude-blocks = 2 5`, the blocks numbered from 1 in the order of
the mesh. A word that is not a block number, or a number beyond the blocks of the mesh, stops the
run at set-up.

---

## Parameter Registry

The **[Parameter Registry](registry.md)** lists the parameters with their defaults and allowed values. It is generated automatically from the source code at every release — do not edit it manually. Three kinds of keys are read at set-up before the registry is built, and the page shows them only in part: `levels` of `[MOSE-Multigrid]`, which it does not list (see [Sections](#sections)); `probe1`, `probe2`, ... of `[MOSE-Probes]`, which it shows as `probe1` with the section of one probe ([Probes](using.md#probes)); `stationary-face-1`, `stationary-face-2`, ... of `[MOSE-Rotating-Frame]`, which it does not list ([Rotating frame](#rotating-frame)).
