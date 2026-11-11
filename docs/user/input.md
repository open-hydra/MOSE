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
    A value outside the rule of its parameter (column Allowed of the [Parameter Registry](registry.md): a bound such as `> 0`, an interval such as `in (0, 1)`, where `(` and `)` leave the end out and `[` and `]` take it in, each held by every value of a list of numbers; a list of names; `logical`) stops the run with a message at the start. MOSE reads `input.ini` again during the run, every `ini-diter` iterations (10000 by default), to take the values changed meanwhile; what an error of such a reading does is in the note «Errors of a reading during the run» below. In an MPI run rank 0 reads `input.ini`, at the start and at every reading, and passes its text to every rank, so that every rank takes the same values. With one working folder per rank, only the `input.ini` of rank 0 is read. The parameters that `input.ini` leaves out keep their values, the `stationary-face-<n>` keys of the rotating frame too; what a reading changes in the rotating frame, and what it does not take, is in [Rotating Frame](../theory/rotating-frame.md#implementation-notes).

!!! warning "Values that do not read stop the run at the start"
    A value that does not read as the type of its parameter in the [Parameter Registry](registry.md) stops the run at the start with a message that gives the file, the section, the key and the value as written (a value of `levels` of `[MOSE-Multigrid]` that is not an integer too; for a reading during the run, see the note below): a real must be a finite number (not a word, `NaN`, `0,5`, or `1e999`, which gfortran reads as an infinity), an integer has no decimal point and no exponent (not `2.0` or `1e3`), a logical starts with `T` or `F` after an optional period (not `yes` or `1`; a word that starts with one of them reads as true or false), a list holds no more values than its parameter, and a key has a value. A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity) stops the run at the start, with a message that names the key.

!!! warning "Errors of a reading during the run"
    When MOSE reads `input.ini` again during the run, an error of the reading (a file that cannot be read, a value that does not read as the type of its parameter, is out of its rule or is not finite, a value of a key fixed for the run other than the start read, a value of `levels` that is not an integer or is not the one the start read or more `probe<n>` keys than at the start, a change of `rotational-frame` or of the stationary faces of the rotating frame, a value of `axis` or `origin` that is not finite or does not read as numbers, a `stationary-face-<n>` value that does not read, an axis of zero length, residual smoothing switched on in a run started without it, values that a check of the set-up between keys refuses) does what `ini-error` of `[MOSE-IO]` says. With `warn`, the default, the reading is not taken: MOSE writes a warning with the iteration and the error, which names the key and its value as written, and the run goes on with every value it had before the reading, those of the keys of the same reading that were right too. With `stop` the run stops with the message and exit code 1, as at the start; with MPI every rank stops. The value of `ini-error` in force when the reading starts decides. At the start every error stops the run, whatever `ini-error` says.

!!! warning "What a reading during the run changes"
    The column During the run of the [Parameter Registry](registry.md) gives what a reading does with each key. A key *taken*, such as `cfl`, `vnn` or the intervals of the outputs, takes the value written. A key *fixed* for the run, such as `equations`, `turbulence`, `time-scheme`, `time-accurate`, `space-reconstruction` or `sol-format`, is one from which the set-up makes a choice that a reading cannot make again: the run keeps the value that the start read, and a reading that writes another one is an error of the reading, also when the value names what the run does (`shock-detector = Tramel` in a run that started without the key, whose sensor comes from `riemann-solver = HLLC+Tramel`); the message gives the value of the start, or says that the start left the key out. The grid levels, the probes and the stationary faces of the rotating frame are those of the start too: another `levels`, more `probe<n>` keys, or a `probe<n>` or a `stationary-face-<n>` that names another probe or face, is an error of the reading. A key of the *start only* (`newrun`, `ic-format`, `decomposition-map`) serves the start of the run, and a reading leaves it as it was, whatever value of its rule the file writes (a value that does not read, or is out of its rule, is an error of the reading, as for every key): `newrun = false` written during the run, to prepare a restart, does nothing to the run. A key that the file leaves out keeps its value. The checks between keys that the set-up makes are made again on the values read: `riemann-options-Mco = 0` with a low-Mach Riemann solver (`AUSM+M`, `HLLC+Chen`, `LMRoe`), `schmidt` = 0 with more than one species in a viscous run without `INPUT/diffusion.dat`, or, with the multirate integrator, a `cfl` above its bound or `cfl-rise-threshold` above 0 (section `[MOSE-Multirate]` below), is an error of the reading as it is at the start. Residual smoothing is on with `irs = true` or `irs-beta` above 0, at the start and at every reading, as [Time Integration](../theory/time-integration.md) says: a reading that switches it on in a run started without it is an error of the reading, since its work arrays exist only when it was on at the start.

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
stability limit instead of the global minimum (method and scope: [Time Integration](../theory/time-integration.md#multirate-time-stepping)):

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

`max-level` is the deepest level allowed: the map of the adaptive step is built for the largest depth up to `max-level`
at which the cells with the largest local step stay at level 0, so that $H$ never exceeds their step divided by $\sigma$
and the formula of the levels leaves no coarse level empty (nesting, the interface margin and the equal levels across
connections only raise levels, and can raise every cell of the slowest octave: the line `Multirate: cells per level`
then starts with 0); without `macro-step-target` the map has $L = \min(\texttt{max-level}, \lfloor K \rfloor)$ levels
below level 0, with $K = \log_2(\max_c \Delta t^{\text{cap}}_c/\min_c \Delta t^{\text{cap}}_c)$
the octaves spanned by the local steps (the finest cells advance with $H/2^L$; on a uniform mesh $L = 0$ and
$H = \Delta t_{min}/\sigma$); the default, 10, is a choice, not a bound of the scheme
([Levels and macro step](../theory/time-integration.md#levels-and-macro-step)).  `safety` is the factor $\sigma$ (every cell
advances with $h \le \Delta t^{\text{cap}}/\sigma$, at least 1; the default 2 keeps every cell, at the start of the macro
step, within the admissible step of each face neighbour without the interface margin, which with `nonlocal-cfl = redo`
gives it for every $\sigma$; with `max-level` = 0, $\sigma = 1$ gives the single-rate step), `levels = recompute` with `recompute-every = 1` recomputes the level
map from the local steps at every macro step.  `macro-step-target` caps $H$ (in seconds; 0 = no cap).  With a target at
or below $\max_c \Delta t^{\text{cap}}_c/\sigma$ the target already keeps $H$ within the step of the slowest cells: the
map is built for the shorter of the target and $2^{\texttt{max-level}}\,\Delta t_{min}/\sigma$, so its depth $L$ is the
one the target gives, at most `max-level`.  With a target above it the map is the one without a target.
The two map policies promise different things.  `levels = recompute` (for production runs) rebuilds the map for the $H$
actually taken, so every cell runs at the coarsest level within its own limit while the state changes.  `levels = fixed`
(for studies of the temporal order) keeps the map of the adaptive step of the first macro step, so that with a smaller
$H$ every cell's step $H/2^{\ell}$ scales with $H$ and the error of a run measures the order of the scheme alone; the
map does not follow the state, so it is not meant for production runs.
`fixed-macro-step` imposes a constant $H$, as a multiple of the minimum local step of the first macro step, for
verification runs.  With `levels = recompute` the level map is built for that $H$ (each cell at the lowest level within
its limit, at most `max-level`), and the run reports the largest ratio $h/(\Delta t^{\text{cap}}/\sigma)$, which exceeds
1 only where `max-level` levels are not enough.  With `levels = fixed` the map is that of the adaptive step: it keeps
every cell within its limit only up to its own step ($2^L\,\Delta t_{min}/\sigma$, or `macro-step-target`), so with an
$H$ above it the finest cells start the first macro step beyond their limit, and every later one until their admissible
step grows by the ratio of $H$ to the step of the map, which the run reports with a warning, once, and with the same
ratio; with `levels = fixed` the certificate of `nonlocal-cfl` only counts.  Every `shell-diter`
iterations the run prints the macro step $H$, the depth $L$ of its map,
the cells per level, the wall time per level spent in the stages since the start of the run and, when process 0 exchanges
data with other processes, the part of it spent waiting for the messages of the exchanges of the stages and copying them
from the receive buffers (these two times are those of process 0), and $r_{post}$, the largest ratio
between the step $h$ each cell took in the last macro step and the step $\Delta t^{\text{cap}}$ admitted on the state it
reached, and $r_{nb}$, the same ratio with the $\Delta t^{\text{cap}}$ of a face neighbour of the cell (condition (4) of
Gnedin, Semenov and Kravtsov, see [Non-local CFL condition](../theory/time-integration.md#non-local-cfl-condition));
when a cell is above 1 or above $1/\text{cfl}$ (a Courant number above 1 where the convective limit rules) a line with
the number of such cells for each.  The lines of the certificate of `nonlocal-cfl` come in a macro step whose last
attempt violated condition (4) or that took more than one attempt, with a line of the attempts rejected for a state
that is not physical when there were some, and the line of the cells whose level the term of the
interface margin raises before nesting, when there are some.  $r_{post}$ and $r_{nb}$ are partial measures: only the state at the end of the macro step, and for $r_{nb}$ only the face
neighbours in the block.  At the end of the run it prints, with 15 digits, the largest $h/(\Delta t^{\text{cap}}/\sigma)$ at the start
of the macro steps when a macro step started with a cell above its limit, and the largest $r_{post}$ and $r_{nb}$ over the run when
a macro step had a cell above 1, each with the number of such macro steps; a run in which none of this happened prints none of
these lines.  A run on several processes ends with the minimum, mean and maximum over the processes, per level, of the wall
time of the stages without the time of their exchanges and of the time of the exchanges: a stage of a level exchanges with
the other processes what it reads from them, so on a level where the processes are unevenly loaded the processes with less
work wait longer in those exchanges (see also [Timing](using.md#timing)); a process with neither cells nor exchanges at a level counts zero
there, in the minimum and in the mean.  These lines only report: they change neither the steps nor the results.
`nonlocal-cfl` checks condition (4) of Gnedin, Semenov and Kravtsov at the start of every sub-step of every cell
([Non-local CFL condition](../theory/time-integration.md#non-local-cfl-condition)): `off`, `count` (the
violations are counted, and printed at the end of a run that has some; nothing else changes) or `redo` (the default: a macro
step that violates it is
redone from its start with finer levels (deeper than the depth $L$ of its map when a cell needs it, at most `max-level`),
or a halved $H$ when a cell would need more than `max-level`; a macro step in which a cell or a face takes a state that
is not physical (a value that is not finite, a total density that is not positive, a negative pressure, a temperature
outside the tables of the gas) is redone from its start with $H$ halved on the same levels, and a face that is not
physical at the first stage of a level, at the start of the macro step, stops the run as in the single-rate step; with
`fixed-macro-step` or `levels = fixed` it only counts, and there, as with a map of one level, the stops of the single-rate
step stay (a cell with a value that is not finite or a negative pressure after its update, a face that is not physical)).  `interface-margin = .true.` keeps the coarse step of a face between levels within
$\Delta t^{\text{cap}}/\kappa$ of the finer cell, $\kappa = \max(\text{cfl}\,(1+\theta),\ 1.1)$ with $\theta$ of the limiter
(1 for `vanleer`, `mc`, `superbee`; 1/2 for `minmod`; $\sqrt2/2$ for `vanalbada`; 0 with `first-order`); it is always on with
`redo`.  Below $\text{cfl}\,(1+\theta) = 1.1$ the margin does not depend on cfl, and no key turns the floor off
([Non-local CFL condition](../theory/time-integration.md#non-local-cfl-condition)).  With
`max-level` > 0 and `nonlocal-cfl = redo`, with `levels = recompute` and no `fixed-macro-step` (where the certificate
redoes the attempts), the run stops at set-up when cfl exceeds 1 or, in a Navier–Stokes run, vnn
exceeds 0.5: within a macro step the certificate accepts a sub-step up to the whole step that cfl and vnn admit on the
state where it starts, and on the checkerboard mode, with every limiter, Heun is linearly stable only within these values
(the table of [Time step](#time-step-cfl-and-vnn)); at the start of the macro step the cells run within $1/\sigma$ of
their step, where a larger value would not show.  With `off` and `count` nothing bounds the sub-steps after the start of
the macro step, and the set-up checks neither value.
`force-redo = .true.` is a test device
of `redo`: the first attempt of every macro step is rejected even when it satisfies (4) and is redone with the same $H$
and level map, so the results must not change.  With `redo` a macro step stops the run, with its rejections by cause,
when its $H$ falls below $2^{-10}$ of its first attempt or when `max-attempts` of its attempts (21 by default; the
attempt rejected by `force-redo` is not counted) have been rejected; written with `count`, `off`, `fixed-macro-step` or
`levels = fixed`, `max-attempts` gives a warning and does nothing.  Whatever the keys, at the start of every macro step
a local step $\Delta t^{\text{cap}}$ that is zero, negative, infinite or not a number (from a state that is not physical,
such as a pressure and a density of opposite signs, a zero pressure at rest in an inviscid run or a density and a
pressure both negative in a viscous run, or from a negative eddy viscosity that makes every diffusivity negative) stops
the run, with the number of such cells and the first one with its local step and its primitive variables.
`conservative-correction = .false.` removes the flux correction between levels and is meant only as a positive control
of the conservation tests.  One iteration is one macro step $H$; the residual printed is the change of the state over
the macro step.  The macro step lands exactly on the output instants $t_k = k\,$`sol-dtime` ($k$ integer) and on
`time-threshold`: when the next of them falls within $H$, $H$ is cut to end there (or stretched by at most
$10^{-12} H$ when the step would end just short of it), the time is set to that instant, and the time-based solution
is written at that macro step, also when it is the last of the run; a cut macro step counts as an iteration.  A macro step whose $H$ a redo of
`nonlocal-cfl = redo` halves does not land: a later macro step reaches the instant.  With the multirate disabled the single-rate
step is unchanged (it writes the solution at the first step past the time of the last output, or of the start, plus
`sol-dtime`).  Time-based probe rows (`[MOSE-Probes]`) come
at the end of the macro steps, which do not land on the probe instants: a probe whose `dtime` is below the macro step
writes one row per macro step, at the time of its end.  A restart (`newrun = false`) from a `vtk raw` or a `tecplot ascii`
solution goes on as the uninterrupted run, bit for bit, with `levels = recompute`: the level map depends only on the state, and
with `recompute-every` = $k > 1$ it is rebuilt at the iterations $1, 1 + k, \dots$ of the run counter, which the
restart carries on, so the restart is bit for bit when it falls on a multiple of $k$ (it rebuilds the map at its first
macro step).  `levels = fixed` and `fixed-macro-step` cannot restart and stop the run at set-up.  The largest ratios
printed at the end of the run cover the macro steps of that run only.  The
multirate keys are fixed for the run (column During the run of the [Parameter Registry](registry.md)): a reading of
`input.ini` during the run that writes another value of one, `enabled` included, is an error of the reading, which
`ini-error` decides (with `warn` the reading is not taken, with `stop` the run stops).  Of the keys of other sections
that the multirate integrator checks at set-up or at its first macro step, `time-accurate` is fixed for every run, and
the checks of `cfl`, `vnn` and `cfl-rise-threshold`, which a reading takes, are made again at every reading, with the
multirate keys of the set-up: with `max-level` > 0 and the redo of `nonlocal-cfl = redo` a cfl above 1 or a vnn above
0.5 in a
Navier–Stokes run, and `cfl-rise-threshold` > 0, is an error of the reading too; a cfl within the bound is taken, and
$\kappa$ of the interface margin follows it from the next macro step, with the $\theta$ of the set-up.  The configurations that the multirate
integrator does not support stop
the run at set-up or at its first macro step with a message
`[MOSE-Multirate] <item>` and a non-zero exit code: among them `Wilcox2006`, whose sources the stages do not evaluate,
and, with `max-level` > 0, `SA-RC` and `SST-RC`, whose rotation and curvature terms are a field of the whole domain that
the cell-wise sources of a level do not compute (`SAcomp` the set-up of MOSE refuses for both integrators); `SA`,
`SA-R`, `SST`, `k-energy-coupling` and both forms of `sst-production` run.  An SRM grain (`502`) runs in mode 2, an
inflow, and stops the run in mode 1, a wall that blows the products in.  The boundary types `409`, `501` and `506`,
which MOSE does not apply, stop every run at
the reading of `bc.txt` ([Boundary Conditions](boundary-conditions.md)).  With `SA`, `SA-R` or `SST` the local step includes a bound of the stiffness of the turbulence
sources (see the time step below), and so do the levels.  The multirate integrator runs on several MPI processes, each with
OpenMP threads, with the same results as on one process.  `trace = .true.` prints on the standard output the sequence of
the stages of the first macro step (lines `MR-TRACE`), its level map (lines `MR-LEVEL`), its halo lists (lines
`MR-HALO`) and, in viscous runs, at every stage the check of the viscous variables that the fluxes read (lines
`MR-VISC`); with several processes, those of the blocks of process 0.  The multirate integrator writes no file of its
own.

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
  ([theory](../theory/turbulence.md#time-accurate-runs-explicit-sources)). With SST or Wilcox 2006,
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

The **[Parameter Registry](registry.md)** lists the parameters with their defaults and allowed values. It is generated automatically from the source code at every release — do not edit it manually. Three kinds of keys are read at set-up before the registry is built, and the page shows them only in part: `levels` of `[MOSE-Multigrid]`, which it does not list (see [Sections](#sections)); `probe1`, `probe2`, ... of `[MOSE-Probes]`, which it shows as `probe1` with the section of one probe ([Probes](using.md#probes)); `stationary-face-1`, `stationary-face-2`, ... of `[MOSE-Rotating-Frame]`, which it does not list ([Wall Boundaries](../theory/rotating-frame.md#wall-boundaries)).
