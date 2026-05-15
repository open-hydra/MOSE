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
    Incorrect names or values are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

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

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.
