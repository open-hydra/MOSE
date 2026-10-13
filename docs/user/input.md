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

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Incorrect names are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

!!! warning "Keys that are not read"
    MOSE warns at the start for every line `key = value` of `input.ini` whose key the reader does not give in its section, with its line: such a key is an error in the input, and the run goes on with the default of its parameter in place of the value written. A parameter of the registry written with no value is named here, and stops the run (below). Comment lines (lines that start with `!`, `;` or `#`) are skipped, also when they hold `=`.

!!! warning "An input.ini that cannot be read stops the run"
    If `input.ini` cannot be opened (missing, or without read permission), holds no line that is not blank, or the INI reader returns an error while it loads the file (a file with no section, a line it cannot parse, which the reader names on its error output), the run stops with the reason and exit code 1: at the start, and when MOSE reads `input.ini` again during the run. With MPI every rank stops. Without the stop the parameters would keep their defaults, or the values of the previous reading, in silence.

!!! warning "Values outside their rules stop the run"
    A value outside the rule of its parameter (column Allowed of the [Parameter Registry](registry.md): a bound such as `> 0`, an interval such as `in (0, 1)`, where `(` and `)` leave the end out and `[` and `]` take it in, each held by every value of a list of numbers; a list of names; `logical`) stops the run with a message: at the start, and when MOSE reads `input.ini` again during the run, every `ini-diter` iterations (10000 by default), to take the values changed meanwhile. The parameters that `input.ini` leaves out keep their values. Of the rotating frame, `omega`, `axis` and `origin` read again enter every term of the frame together (the log of the run gives them); a change of `rotational-frame` or of the `stationary-face-<n>` keys, or a value of `axis` or `origin` that is not finite (NaN or infinity), stops the run ([Rotating Frame](../theory/rotating-frame.md#implementation-notes)).

!!! warning "Values that do not read stop the run"
    A parameter of the [Parameter Registry](registry.md) whose value in `input.ini` does not read as its type stops the run with a message that gives the file, the section, the key and the value as written: a number that is not one (a word, `NaN`, `Infinity`, a decimal comma such as `0,5`), an integer written with a decimal point or an exponent (`2.0`, `1e3`: one compiler refuses it, another cuts the fraction), a logical that does not start with `T` or `F` after an optional period (`yes`, `1`; a word that starts with one of them reads as true or false), a list with more values than its parameter holds, a key with no value. So does a real read as an infinity (gfortran reads a number too large for a real, such as `1e999`, as one), with a rule or without: a bound such as `> 0` let it through. This holds at the start and when MOSE reads `input.ini` again during the run; `levels` of `[MOSE-Multigrid]`, which MOSE reads at the start only, stops the run as well. Before, such a value left its parameter at the default, or at the value of the previous reading, in silence. A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity) stops the run at the start, with a message that names the key.

!!! warning "What a reading during the run changes"
    When MOSE reads `input.ini` again during the run, only the parameters that the solver reads at every step take the new values, such as `cfl`, `vnn` and the intervals of the outputs. The choices that the set-up makes once from the input are not made again, nor are its checks between keys: a new value of `equations`, `turbulence`, `chemistry`, `space-reconstruction`, `flux-limiter`, `riemann-solver`, `shock-detector` or `time-scheme` is checked against its rule and leaves the run as it was. Do not turn residual smoothing on during a run (`irs = true`): its work arrays exist only when it was on at the start, and the run stops with a segmentation fault.

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
`level1-iter`, `level2-iter`, ... give the iterations of each level.

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

The rule used to take the smallest of the direction limits, the edge length of the cell
instead of its width, and the smaller of the convective and the viscous limit, with the
viscous one built on $\nu$. Each of these let the step exceed the stability limit of the
scheme: on square 2-D cells by up to a factor 2, on skewed cells more, in a cell set by heat
conduction by the factor 1.94. **The same `cfl` and `vnn` therefore give a smaller step than
before**, most of all in viscous cells: on a turbulent nozzle case (Spalart–Allmaras, RK3,
`cfl = 0.8`, `vnn = 0.3`) the global time-accurate step is 0.59 times the old one.

Stability limits of the linearised schemes (MUSCL with upwind fluxes, central diffusion);
the additive rule stays stable for any `cfl` and `vnn` up to these values:

| `time-scheme` | `cfl` ≤ | `vnn` ≤ |
|---|---|---|
| `RK2` (Heun) | 1 | 0.5 |
| `RK3` | 1.17 | 0.62 |
| `euler` | 1 with `first-order` reconstruction (unstable with `MUSCL`) | 0.5 |

These are linear limits. With the limiter, positivity and the TVD property are guaranteed
only up to `cfl = 0.5` and `vnn = 0.5`, with Heun and with RK3
([theory](../theory/time-integration.md#local-time-step)), and the `vnn` limits are sharp: a run
a little beyond them blows up. Suggested values:

- `RK3`: `cfl = 0.8` (above the guaranteed 0.5) and `vnn = 0.5`;
- `RK2` (Heun): `cfl = 0.5` with moving shocks, up to `cfl = 0.7` where the flow is smooth, and
  `vnn = 0.4`; never `cfl` ≥ 0.95 or `vnn` ≥ 0.5.

To recover part of the old step of a viscous time-accurate run, raise `vnn` rather than
`cfl`: in the nozzle case above, `cfl = 0.8` with `vnn = 0.5` gives 0.84 of the old global
step with RK3, and `vnn = 0.52`, beyond the positivity limit, hardly more (0.84): the step of
a time-accurate RANS run also holds the bound of the stiffness of the turbulence sources,
which `vnn` does not scale. With Heun, `cfl = 0.5` and `vnn = 0.4` give 0.67 of it. The
defaults stay `cfl = 0.5` and `vnn = 0.3`.

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
