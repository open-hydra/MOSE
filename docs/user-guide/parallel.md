# Parallel execution

MOSE supports OpenMP, MPI, and hybrid MPI + OpenMP execution.

- **OpenMP** uses multiple threads within one process.
- **MPI** distributes mesh blocks across processes (ranks).
- **Hybrid execution** uses MPI across ranks and OpenMP threads within each rank.

The same number of CPU cores can give different performance depending on where ranks and threads are placed.

## A practical starting point

On a NUMA machine, start with at least one MPI rank per socket and keep each rank's OpenMP threads within a single socket. If possible, choose a rank count that divides the number of mesh blocks.

For example, an 8-rank run with 10 threads per rank uses 80 CPU cores:

```bash
export OMP_NUM_THREADS=10
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp   # Intel MPI binding
mpirun -n 8 /path/to/MOSE/bin/MOSE
```

This is an example layout, not a universal optimum. Check the socket/core layout and the binding options supported by your MPI runtime and scheduler.

## Before scaling up

1. Confirm the number of sockets and physical cores per socket (`lscpu` or `numactl --hardware`).
2. Bind MPI ranks and OpenMP threads explicitly.
3. Verify the placement reported by the MPI runtime.
4. Benchmark a representative case using MOSE's `timer-diter` timing report.

If performance or memory use is poor, check the number of blocks per rank, load balance, and output frequency before adding more ranks.

For MPI launch variants, block-to-rank maps, affinity checks, measured scaling results, and known limits, see [Parallel execution details](../user/parallel.md).
