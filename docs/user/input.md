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

!!! warning "A comment line that holds `=` ends the reading of its section"
    The reader of `input.ini` stops reading a section at a line it cannot parse, such as a comment line (a line that starts with `!`, `;` or `#`) that holds `=`: the keys after it in that section are not read and keep their defaults. MOSE warns at the start for every such key, with its line. A key lost this way is an error in the input: the run goes on with the default in place of the value written. Keep `=` out of comment lines.

!!! warning "Values outside their rules stop the run"
    A value outside the rule of its parameter (column Allowed of the [Parameter Registry](registry.md)) stops the run with a message: at the start, and when MOSE reads `input.ini` again during the run, every `ini-diter` iterations (10000 by default), to take the values changed meanwhile. The parameters that `input.ini` leaves out keep their values.

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

To recover the old step of a viscous time-accurate run, raise `vnn` rather than `cfl`: in the
nozzle case above, `cfl = 0.8` with `vnn = 0.5` gives 0.96 of the old global step with RK3
(the whole old step needs `vnn = 0.52`, beyond the positivity limit). With Heun, `cfl = 0.5`
and `vnn = 0.4` give 0.76 of it. The defaults stay `cfl = 0.5` and `vnn = 0.3`.

---

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.
