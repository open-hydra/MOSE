# Common switches

Find what you want to change, copy the key, and follow the link if you need the details. Every key goes in `input.ini`; omitted keys use the default shown. For all keys and exact allowed values, see the [Parameter Registry](../user/registry.md).

!!! tip "How to read the tables"
    **Mid-run** says what happens when MOSE re-reads `input.ini` during a run: **yes** = the new value is used, **no** = chosen at startup and cannot change, **start** = used only at startup. Logicals accept `true` or `false`.

## Required keys

These five have no usable default. A case needs all of them.

| Key | Section | Values |
|---|---|---|
| `equations` | `[MOSE-Physics]` | `euler`, `navier-stokes` |
| `time-scheme` | `[MOSE-Numerics]` | `euler`, `RK2`, `RK3` |
| `cfl` | `[MOSE-Numerics]` | number > 0 |
| `time-accurate` | `[MOSE-Numerics]` | `true`, `false` |
| `space-reconstruction` | `[MOSE-Numerics]` | `MUSCL`, `first-order` |

## Steady or time-accurate

| I want to… | Set | Mid-run |
|---|---|---|
| Converge to a steady state | `time-accurate = false` | no |
| Follow a transient | `time-accurate = true` | no |
| Stop after N iterations | `iter-threshold = N` | yes |
| Stop at a physical time | `time-threshold = t` (seconds) | yes |
| Stop when converged (steady only) | `res-threshold = 1e-6` (drop relative to iteration 1; default `1e-10`) | yes |

## Physics

| I want to… | Set | Mid-run | More |
|---|---|---|---|
| Solve viscous flow | `equations = navier-stokes` | no | |
| Add a turbulence model | `turbulence = SA` (or `SST`, `SSGLRR`, …; default `none`) in `[MOSE-Physics]`; tune in `[MOSE-Turbulence]` | no | [Registry](../user/registry.md#mose-turbulence) |
| Add finite-rate chemistry | `chemistry = finite-rate` (default `frozen`); needs `chemistry-info.txt` and `chemistry-*.dat` in `INPUT/` | no | [Case files](../user/using.md#case-directory-structure) |
| Use equilibrium chemistry | `chemistry = equilibrium`; needs `composition.txt` in `INPUT/` | no | [Case files](../user/using.md#case-directory-structure) |
| Choose the chemistry ODE solver | `ode-solver = H-radau5` (default) in `[MOSE-Chemistry]` | no | [Registry](../user/registry.md#mose-chemistry) |
| Run in a rotating frame | `rotational-frame = rigid-body`, then `omega`, `axis`, `origin` in `[MOSE-Rotating-Frame]` | no (`omega`, `axis`, `origin`: yes) | [Registry](../user/registry.md#mose-rotating-frame) |

## Numerics

| I want to… | Set | Mid-run | More |
|---|---|---|---|
| Change the time step | `cfl` (convective) and `vnn` (diffusive, default `0.3`) | yes | [Time-step limits](../user/input.md#time-step-cfl-and-vnn) |
| Pick the Riemann solver | `riemann-solver = HLLC` (default; many others, see registry) | no | [Registry](../user/registry.md#mose-numerics) |
| Go first-order | `space-reconstruction = first-order` | no | |
| Change the limiter (MUSCL) | `flux-limiter = minmod` (also `vanalbada`, `superbee`, `vanleer`, `mc`; van Leer if omitted) | no | |
| Speed up steady convergence | `irs = true` or `irs-beta > 0` (steady only; cannot be switched on mid-run) | yes, if on at start | [Models and methods](../user/input.md#models-and-methods) |
| Precondition for low Mach (steady) | `integration-variables = prec` with `riemann-solver = HLLC-PC` | no | [Registry](../user/registry.md#mose-numerics) |

## Multigrid and multirate

| I want to… | Set | Mid-run | More |
|---|---|---|---|
| Start on coarse grids (grid sequencing) | `[MOSE-Multigrid]` `levels = 3`, plus `level<n>-iter` per level; one `bc<n>.txt` per level; cell counts divisible by 2^(levels−1) | `levels`: no; `level<n>-iter`: yes | [Multigrid & FAS](../user/multigrid.md) |
| Add a V or W cycle (steady, no chemistry) | `cycle = v` or `w` in `[MOSE-Multigrid]` | yes, but not switched on mid-run | [Multigrid & FAS](../user/multigrid.md#grid-sequencing-or-fas) |
| Use local time stepping (time-accurate) | `[MOSE-Multirate]` `enabled = true`, with `time-scheme = RK2`, `time-accurate = true`, `cfl ≤ 1` | no | [Multirate](../user/multirate.md) |

## Output

| I want to… | Set | Mid-run |
|---|---|---|
| Choose the solution format | `sol-format = tecplot ascii` (default), `tecplot binary`, `vtk ascii`, `vtk raw` | no |
| Write the solution every N iterations | `sol-diter = N` | yes |
| Write the solution every Δt (time-accurate) | `sol-dtime = Δt` | yes |
| Keep numbered snapshots instead of overwriting | `sol-overwrite = false` | yes |
| Choose what is written | `sol-variables`, `wall-variables` | no |
| Record a point over time | `[MOSE-Probes]` `probe1 = name`, then a `[name]` section | no |
| Change residual-history frequency | `res-diter = N` (default `1`) | yes |
| Print a timing report | `timer-diter = N` (default `0` = off) | yes |

By default MOSE writes the solution only when the run ends. For long jobs, set `sol-diter` or `sol-dtime` and `sol-overwrite = false` before launching: see [Restart a simulation](restart.md#keep-intermediate-checkpoints).

## Start, restart and re-reading

| I want to… | Set | Mid-run |
|---|---|---|
| Start from `INPUT/ic.*` (default) | `newrun = true` | start |
| Continue a saved run | `newrun = false` | start |
| Read the initial condition in another format | `ic-format = vtk raw` | start |
| Change how often `input.ini` is re-read | `ini-diter = N` (default `10000`) | yes |
| Stop instead of ignoring a bad re-read | `ini-error = stop` (default `warn`) | yes |

A value outside its allowed range, or the wrong type, stops MOSE at startup. See [Input file essentials](input-file.md) and [Validation & runtime changes](../user/input-validation.md).
