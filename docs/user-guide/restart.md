# Restart a simulation

A restart continues an existing run from a solution file that MOSE has already written. It is useful when a simulation is interrupted, reaches a queue time limit, or needs to continue for more iterations.

**You do not need a restart to run a case for the first time.** A fresh run uses `newrun = true`, which is the default.

## Restart or start fresh?

| What you want | Setting | Starting point |
|---|---|---|
| Start a new simulation | `newrun = true` | `INPUT/ic.*`; time and iteration start at zero |
| Continue a previous simulation | `newrun = false` | A saved solution in `OUTPUT/`; saved run state is restored when available |
| Re-use an old field but begin a new simulation | `newrun = true` | Copy the old solution to `INPUT/ic.*`; counters and time start from zero |

## How to resume

1. **Keep the case directory and its `OUTPUT/` files.** A restart reads a solution already written by MOSE.
2. **Check which solution MOSE will read.** MOSE first looks for `OUTPUT/field.*` in the format selected by `sol-format`. If that file exists, it is used. Only when it is absent does MOSE look for numbered snapshots (`field1.*`, `field2.*`, and so on).
3. **Set the startup option** in `input.ini`:

   ```ini
   [MOSE-Parameters]
   newrun = false
   ```

4. **Launch MOSE again** from the same case directory, using the same launch command as before.

The `newrun` setting is read at startup. Changing it while a solver process is running has no effect on that process.

!!! warning "Check for an older `field.*` file"
    `OUTPUT/field.*` always takes priority over numbered snapshots. If the current job was interrupted but `field.*` is left over from an older completed run, MOSE will read that older file instead of a newer `field<N>.*`. Move the stale solution aside before restarting from numbered snapshots (for VTK, keep the referenced block files with the `.vtm` file).

    If `field.*` is absent, MOSE checks numbered snapshots in sequence starting at `field1.*` and stops at the first missing number. Keep the sequence contiguous; a gap means later-numbered files are not considered.

## Keep intermediate checkpoints

To write restartable snapshots before a run completes, set both an output interval and `sol-overwrite = false` in `[MOSE-IO]`.

For a steady run:

```ini
[MOSE-IO]
sol-diter = 500
sol-overwrite = false
```

For a time-accurate run:

```ini
[MOSE-IO]
sol-dtime = 0.01
sol-overwrite = false
```

With `sol-overwrite = false`, MOSE writes intermediate solutions as numbered `field<N>.*` files. A normally completed run also writes `field.*` as its final solution. With the default `sol-overwrite = true`, intermediate solution writes are not kept as restart files; `field.*` is written at normal completion.

The defaults `sol-diter = 1000000000` and `sol-dtime = 1e30` do not provide practical periodic checkpoints for most jobs.

## What is restored?

A restart restores the flow field and its progress marker: physical time for a time-accurate run, or the global iteration count for a steady run.

Solutions written by current MOSE versions also store run counters. When those counters are present and the restart is of the same kind (steady or time-accurate), MOSE can continue the output numbering, output schedule, multigrid position, and steady residual reference. This is why a split run can continue like one uninterrupted run in the cases covered by the restart tests.

For a compatible restart, MOSE also removes rows written after the checkpoint from `residual-history.dat`, probe files, and (when used) `torque.dat`, then appends the rows from the continued run.

## Cases that need attention

- **No time in the solution:** a VTK solution without `TIME` or a Tecplot ASCII solution without `SOLUTIONTIME` cannot be restarted. Use it as an initial condition instead.
- **Steady solution used to start a time-accurate run:** this starts a new transient at time zero; it does not continue a transient.
- **Solution from an older MOSE version:** older solutions do not contain the current run counters. MOSE can restore the field and time/iteration, but output numbering and other counters start as for a new run. Existing history files are not trimmed automatically in this case.
- **Not every configuration is bitwise restart-identical:** see the detailed [restart limitations](../user/restart.md#limits), especially for multigrid, probe schedules, and some preconditioning options.

For reliable recovery from a long-running job, configure periodic output before the job starts. See [Save checkpoints](run-a-case.md#5-save-checkpoints-for-long-jobs).
