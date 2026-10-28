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
    When MOSE reads `input.ini` again during the run, an error of the reading (a file that cannot be read, a value that does not read as the type of its parameter, is out of its rule or is not finite, a value of a key fixed for the run other than the start read, a value of `levels` that is not an integer or is not the one the start read or more `probe<n>` keys than at the start, a change of `rotational-frame` or of the stationary faces of the rotating frame, a value of `axis` or `origin` that is not finite or does not read as numbers, a `stationary-face-<n>` value that does not read, an axis of zero length) does what `ini-error` of `[MOSE-IO]` says. With `warn`, the default, the reading is not taken: MOSE writes a warning with the iteration and the error, which names the key and its value as written, and the run goes on with every value it had before the reading, those of the keys of the same reading that were right too. With `stop` the run stops with the message and exit code 1, as at the start; with MPI every rank stops. The value of `ini-error` in force when the reading starts decides. At the start every error stops the run, whatever `ini-error` says.

!!! warning "What a reading during the run changes"
    The column During the run of the [Parameter Registry](registry.md) gives what a reading does with each key. A key *taken*, such as `cfl`, `vnn` or the intervals of the outputs, takes the value written. A key *fixed* for the run, such as `equations`, `turbulence`, `time-scheme`, `time-accurate`, `space-reconstruction` or `sol-format`, is one from which the set-up makes a choice that a reading cannot make again: the run keeps the value that the start read, and a reading that writes another one is an error of the reading, also when the value names what the run does (`shock-detector = Tramel` in a run that started without the key, whose sensor comes from `riemann-solver = HLLC+Tramel`); the message gives the value of the start, or says that the start left the key out. The grid levels, the probes and the stationary faces of the rotating frame are those of the start too: another `levels`, more `probe<n>` keys, or a `probe<n>` or a `stationary-face-<n>` that names another probe or face, is an error of the reading. A key of the *start only* (`newrun`, `ic-format`, `decomposition-map`) serves the start of the run, and a reading leaves it as it was, whatever value of its rule the file writes (a value that does not read, or is out of its rule, is an error of the reading, as for every key): `newrun = false` written during the run, to prepare a restart, does nothing to the run. A key that the file leaves out keeps its value. The checks between keys that the set-up makes are not made again. `irs = true`, or `irs-beta` above 0, written during a run started without residual smoothing is ignored with a warning, and the rest of the reading is taken: its work arrays exist only when it was on at the start.

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
([theory](../theory/time-integration.md#local-time-step)), and the `vnn` limits are sharp: a run
a little beyond them blows up. Suggested values:

- `RK3`: `cfl = 0.8` (above the guaranteed 0.5) and `vnn = 0.5`;
- `RK2` (Heun): `cfl = 0.5` with moving shocks, up to `cfl = 0.7` where the flow is smooth, and
  `vnn = 0.4`; never `cfl` ≥ 0.95 or `vnn` ≥ 0.5.

The defaults are `cfl = 0.5` and `vnn = 0.3`; in a viscous time-accurate run raise `vnn`
within the limits of the table before `cfl`.

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
