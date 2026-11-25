# Run diagnostics and troubleshooting

This page explains what to look for in the MOSE log, how a run stops, and where its output is written. For the normal launch procedure, see [Running MOSE](using.md).

## Read the startup summary first

During startup, MOSE reports the files it loads and prints a setup summary. Check that the summary matches the case you intended to run: mesh dimensions and block count, boundary-condition types, physical models, and numerical settings.

A missing or malformed required input normally produces an `[ERROR]` message and stops setup. Start with the **first error message**, not the last line printed by the launcher. The relevant file reference is usually named in the message.

| Message points to… | Check… |
|---|---|
| `input.ini` | Section and parameter spelling, value type, and allowed range. See [Input validation](input-validation.md). |
| `INPUT/ic.*` | File selected by `ic-format`, required variables, mesh coordinates, and finite values. See [Initial conditions](initial-conditions.md). |
| `INPUT/bc.txt` | Block, face, boundary type, and type-specific data line. See [Boundary conditions](boundaries/index.md). |
| Thermodynamic, transport, or chemistry data | The files required by the selected model and the gas tables used by the case. |

## Exit status and parallel runs

A successful run exits normally. An error detected during input loading or setup usually gives a non-zero exit status. Integration failures and invalid solution states also stop the run.

With MPI, many setup and input errors are coordinated so that all ranks stop together. Some errors detected by an individual rank, and some failures that use `MPI_ABORT`, may produce launcher-specific exit codes or incomplete output from other ranks. The message written by MOSE is more useful than assuming a particular numeric exit code across compilers and MPI implementations.

**Practical check:** inspect the MOSE log for the first `[ERROR]` or `[STOP]` message, then check the named block, cell, file, or parameter. Treat the exit code as a success/failure signal, not as a portable description of the cause.

## Finite-rate chemistry integration failure

When the chemistry integrator reports a failed cell integration, MOSE stops the current step before the flux calculation can use that cell. All threads and MPI ranks stop together; MOSE does not write a solution for that failed step.

The diagnostic message identifies the integrator and return code, block and grid level, cell, chemistry step, temperature, pressure, density, and (with MPI) rank. It reports one failing cell deterministically: the lowest block, then the lowest `k`, `j`, and `i` indices.

Common things to inspect are the integrator's return-code meaning, the mechanism and thermodynamic tables, the chemistry tolerances, and `ode-max-steps`. Do not assume that increasing `ode-max-steps` alone will fix a mechanism or state that is outside the model's usable range.

## End-of-run diagnostic counters

MOSE prints additional lines only when the corresponding condition occurred:

| Log line | Meaning |
|---|---|
| `Turbulence variables below zero after a stage (set to the floor): <n>` | A turbulence stage produced a negative working variable before the realizability floor was applied. |
| `Ghost cells of walls and 400 faces, fills with the extrapolation cut: <n>, with the 302 temperature at T_w/2: <m>` | The ghost-cell safeguards for walls or extrapolation (`400`) faces were used. Either count may be zero. |

No such line means that the corresponding safeguard did not need to act during that run. These counts are diagnostics, not necessarily fatal errors. See [Turbulence numerical treatment](../theory/turbulence-numerics.md) and [ghost cells at walls and extrapolation faces](../theory/numerics-ghost-cells.md).

## Stopping criteria

A run stops when it reaches one of these limits:

- `iter-threshold`: maximum iterations;
- `time-threshold`: maximum physical time, for time-accurate runs;
- `res-threshold`: steady runs only, when the density residual has fallen below this factor times the density residual of the first iteration of that grid level.

The density residual is an L2 norm over the cells. The value printed to the screen and written to `residual-history.dat` is the residual itself; the convergence check compares it with the first-iteration reference. A residual that rises before falling is still compared with the first iteration, not with its later peak. Time-accurate runs report residuals but do not stop on `res-threshold`.

A compatible restart preserves the original residual reference. Older solutions without that stored value use the first iteration of the restarted run as a new reference; see [Restart a simulation](restart.md).

## Timing reports

Set `timer-diter` in `[MOSE-IO]` to print a timing report every specified number of iterations:

```ini
[MOSE-IO]
timer-diter = 100
```

The main values are:

| Timing field | How to read it |
|---|---|
| `wall/iter` | Wall-clock time per iteration; useful for comparing configurations. |
| `imbalance` | How much slower the slowest rank is than the mean rank time. |
| `exchange wait` | Time waiting for MPI halo/connection exchanges. |
| `collective wait` | Time spent in MPI reductions and broadcasts. |
| `compute/iter` spread | Difference between the fastest and slowest rank's compute time. |

The reported critical-path times are the maximum over ranks. A high compute spread with a small overall imbalance can mean that faster ranks are waiting for communication with a slower rank. Timing output is disabled by default (`timer-diter = 0`).

## Output files

The files written depend on the selected output format and output settings.

| Output | Contents |
|---|---|
| `OUTPUT/field.*` | Final solution, or the overwritten solution when `sol-overwrite = true`. |
| `OUTPUT/field<N>.*` | Numbered solution snapshots when `sol-overwrite = false`. |
| `OUTPUT/residual-history.dat` | Residual history. |
| `OUTPUT/wall.*` | Wall quantities when wall output is enabled (normally viscous cases). |
| `OUTPUT/diagnostic.*` | Per-cell diagnostics for steady runs. |
| `OUTPUT/<probe>.txt` | Time history at a configured probe. |
| `OUTPUT/field-level<N>.*` | Solution on a coarse multigrid level. |
| `OUTPUT/field-prolongated.*` | Coarse-level solution prolongated to the finest grid. |

The extension and any VTK companion files depend on `sol-format`. See [Running MOSE](using.md#output) for the format list and [Input Parameters](registry.md) for output intervals.

## Probe quick start

A probe is configured by naming it in `[MOSE-Probes]` and adding a section with that name:

```ini
[MOSE-Probes]
probe1 = exit

[exit]
variables = p T u
dtime = 1e-4
```

This example writes `OUTPUT/exit.txt`. A probe records one cell, selected either by `index-position` (block and cell indices) or by the nearest cell centre to `position`. The allowed variables and restart behaviour are described in [Running MOSE — probes](using.md#probes).

## Next steps

- [Input validation](input-validation.md) for startup errors and runtime reads of `input.ini`.
- [Boundary conditions](boundaries/index.md) for errors that name a face or BC type.
- [Restart a simulation](restart.md) for continuing after a job stops or is interrupted.
