# Run a case

This page covers the normal path from a prepared case directory to a running simulation.

## 1. Check the case files

Run MOSE from the case directory that contains `input.ini`. The case needs:

| Path | Purpose |
|---|---|
| `input.ini` | Simulation settings |
| `INPUT/` | Initial condition, boundary conditions, and any data required by the selected physical models |
| `OUTPUT/` | Created by MOSE for solution files, residual history, probes, and other output |

The initial-condition file also contains the mesh coordinates. There is no separate mesh file.

The required files inside `INPUT/` depend on the models selected in `input.ini`. For example, finite-rate chemistry needs its mechanism data. See [Initial conditions](initial-conditions.md) and the [full running guide](../user/using.md#case-directory-structure) for the file list.

## 2. Launch MOSE

For an OpenMP run, set the number of threads and launch the executable:

```bash
export OMP_NUM_THREADS=4
ulimit -s unlimited
/path/to/MOSE/bin/MOSE
```

Run this command from the case directory, replacing the executable path with the location of your MOSE build.

MOSE may need a large stack. On Intel compiler builds, `KMP_STACKSIZE=100M` can be used instead of increasing the shell stack limit.

For MPI + OpenMP, launch one process per MPI rank:

```bash
export OMP_NUM_THREADS=10
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp   # Intel MPI
ulimit -s unlimited
mpirun -n 8 /path/to/MOSE/bin/MOSE
```

The process-binding option depends on the MPI implementation and the cluster scheduler. See [Parallel execution](parallel.md) before moving to a production allocation.

## 3. Read the run log

During startup, MOSE reports which input files it loads and prints a summary of the mesh, boundary conditions, physical models, and numerical settings.

During the run, the log reports iteration count and residuals for steady simulations, or time and time step for time-accurate simulations. A loading or setup error stops the run; do not treat the presence of an output directory alone as proof that setup succeeded.

A non-zero process exit status indicates that MOSE or the MPI launcher reported an error. The detailed meanings of individual errors and exit codes are collected in the [full running guide](../user/using.md#logging-and-diagnostics).

## 4. Find the output

Common files in `OUTPUT/` include:

| File | What it contains |
|---|---|
| `field.*` | Final solution written at normal completion, or the solution file selected for restart when present |
| `field<N>.*` | Numbered intermediate solution snapshots when `sol-overwrite = false` |
| `residual-history.dat` | Residual history |
| `wall.*` | Wall quantities for viscous cases |
| Probe files | Time histories at configured probe locations |

The file extension depends on `sol-format` in `[MOSE-IO]`.

## 5. Save checkpoints for long jobs

A restart can only continue from a solution that has already been written. **To keep intermediate checkpoints, set `sol-overwrite = false`.** With the default `sol-overwrite = true`, intermediate solution writes are not kept as restart files; MOSE writes `field.*` when the run completes normally.

For a steady run, for example:

```ini
[MOSE-IO]
sol-diter = 500
sol-overwrite = false
```

For a time-accurate run, you can schedule output by physical time instead:

```ini
[MOSE-IO]
sol-dtime = 0.01
sol-overwrite = false
```

Choose an interval that balances the amount of work you are willing to lose against output cost and disk usage. Intermediate snapshots are written as numbered `field<N>.*` files. A normally completed run also writes the final `field.*` solution.

The default values (`sol-diter = 1000000000` and `sol-dtime = 1e30`) do not provide practical periodic checkpoints for most jobs. See [Restart a simulation](restart.md) before relying on restart for fault recovery.

## 6. Know when the run stops

A run stops when it reaches `iter-threshold` or, for a time-accurate run, `time-threshold`. A steady run can also stop when the density residual falls below `res-threshold` times the density residual of its first iteration. Residual convergence is not a stopping criterion for time-accurate runs.

For a complete description of timing reports, diagnostics, and output cadence, see the [full running guide](../user/using.md).
