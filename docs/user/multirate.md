# Multirate time stepping

Multirate time stepping lets cells advance with different local sub-steps instead of forcing every cell to use the smallest stable time step in the domain. In MOSE, it is available for **time-accurate runs** and uses a recursive Heun (SSP-RK2) method.

For the numerical method and its analysis, see [Multirate time stepping theory](../theory/multirate.md) and the [multirate verification suite](../vv/multirate/index.md).

## Minimal configuration

```ini
[MOSE-Numerics]
time-scheme = RK2
time-accurate = true
integration-variables = cons

[MOSE-Multirate]
enabled = true
max-level = 10
safety = 2.0
levels = recompute
recompute-every = 1
```

`enabled` is fixed at startup. Multirate is not supported for a steady run; use the single-rate integrator for steady simulations.

## Main settings

| Parameter | Default | Purpose |
|---|---:|---|
| `enabled` | `.false.` | Enable multirate stepping. |
| `max-level` | `10` | Maximum depth of the local step hierarchy. The finest cells can take steps of `H/2^L`. |
| `safety` | `2.0` | Safety factor applied to each cell's local stability limit. Must be at least one. |
| `macro-step-target` | `0` | Optional upper target for the macro step `H`, in seconds; zero leaves it uncapped by this setting. |
| `levels` | `recompute` | `recompute` rebuilds the local level map; `fixed` retains the first map for verification studies. |
| `recompute-every` | `1` | Number of macro steps between map rebuilds when `levels = recompute`. |
| `fixed-macro-step` | `0` | Fixed macro step as a multiple of the smallest local step from the first macro step; intended for verification. Zero is adaptive. |
| `nonlocal-cfl` | `redo` | Check mode: `off`, `count`, or `redo`. |
| `interface-margin` | `.false.` | Adds a margin to the coarse step at faces between levels; always enabled with `nonlocal-cfl = redo`. |
| `max-attempts` | `21` | Maximum number of rejected attempts of one macro step when `redo` is active. |
| `force-redo` | `.false.` | Test control that deliberately rejects the first attempt of every macro step. |
| `conservative-correction` | `.true.` | Conservative flux correction across faces between levels. Keep enabled for production. |
| `trace` | `.false.` | Print the stage sequence, level map, and halo lists of the first macro step. |

## The `nonlocal-cfl` check

| Mode | Behaviour |
|---|---|
| `off` | Do not check the non-local CFL condition. |
| `count` | Count and report violations; do not alter the step. |
| `redo` | Reject a macro step that violates the condition and retry it from the beginning with finer levels or a smaller macro step. |

With `redo`, MOSE also retries a macro step if a cell or face reaches a non-physical state during that step. A face that is already non-physical at the beginning of the macro step stops the run.

The run stops if the macro step becomes smaller than `2^-10` of the first attempt or if `max-attempts` rejected attempts are reached. The end-of-run report gives rejection counts by cause.

## CFL and VNN limits

When the production combination is active (`max-level > 0`, `nonlocal-cfl = redo`, `levels = recompute`, and no fixed macro step), setup and runtime input checks reject:

- `cfl > 1`;
- `vnn > 0.5` for Navier–Stokes runs.

These are the Heun linear-stability limits. See [Time step: CFL and VNN](input.md#time-step-cfl-and-vnn) for the broader single-rate stability discussion.

## Output and restart behaviour

- A macro step is cut so that it lands on scheduled `sol-dtime` output instants and `time-threshold`.
- Time-based probe rows are written at the end of macro steps. If a probe's `dtime` is smaller than the macro step, it writes one row per macro step.
- The multirate integrator does not write a separate file of its own; normal solution, wall, diagnostic, and probe outputs apply.
- `levels = fixed` and a non-zero `fixed-macro-step` are verification settings and cannot be used for restart. The documented bitwise-restart case uses a `vtk raw` or `tecplot ascii` solution with `levels = recompute`; for `recompute-every = k`, the restart must fall on a multiple of `k`. Do not assume bitwise identity for other combinations.

For how to configure checkpoint output and resume a job, see [Restart & Recovery](restart.md). For the exact limitations, see [multirate verification](../vv/multirate/index.md).
