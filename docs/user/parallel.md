# Parallel Execution and Scalability

MOSE supports **hybrid MPI + OpenMP parallel execution**. MPI distributes the computational domain across processes, while OpenMP parallelizes the work performed within each MPI process.

This allows the available CPU cores to be distributed between MPI processes and OpenMP threads in several ways. For example, on a 40-core allocation, MOSE can be executed as:

```text
1 MPI rank  × 40 OpenMP threads
4 MPI ranks × 10 OpenMP threads
8 MPI ranks × 5 OpenMP threads
40 MPI ranks × 1 OpenMP thread
```

These configurations use the same number of CPU cores, but they do **not** necessarily provide the same performance. The choice of MPI/OpenMP decomposition affects NUMA locality, memory bandwidth, domain decomposition, communication overhead, and memory consumption.

> **Recommended configuration**
>
> On a NUMA system, use **at least one MPI rank per socket** and keep the OpenMP threads belonging to each rank within a single socket whenever possible.

The remainder of this page explains how to select a parallel configuration, launch MOSE, verify process and thread placement, and evaluate scaling across multiple cores and nodes.

## 1. Selecting a Parallel Configuration

Three properties of the machine and simulation determine a suitable MPI/OpenMP configuration:

| Symbol | Description                               | How to obtain it                       |
| ------ | ----------------------------------------- | -------------------------------------- |
| `S`    | Number of sockets (NUMA domains) per node | `numactl --hardware`                   |
| `C`    | Number of physical cores per socket       | `lscpu`                                |
| `B`    | Number of mesh blocks                     | MOSE setup output or decomposition log |

For a given core budget `N`, the preferred configuration follows these rules:

1. **Cover all sockets.** Use at least one MPI rank per socket.
2. **Keep OpenMP teams within a socket.** Choose a thread count that divides the number of cores per socket.
3. **Avoid unnecessary MPI ranks.** Prefer a rank count that divides the number of mesh blocks when possible.

For example, on a node with **4 sockets × 20 cores**, the recommended configurations include:

| CPU cores | MPI ranks | OpenMP threads/rank | Configuration |
| --------: | --------: | ------------------: | ------------- |
|        20 |         4 |                   5 | `4 × 5`       |
|        40 |         4 |                  10 | `4 × 10`      |
|        80 |         8 |                  10 | `8 × 10`      |
|       160 |        16 |                  10 | `16 × 10`     |
|       240 |        24 |                  10 | `24 × 10`     |

The important property is not the MPI rank count by itself, but **socket coverage and thread locality**.

### Configurations to Avoid

Avoid configurations in which:

* a single MPI rank spans multiple sockets;
* an OpenMP team crosses a socket boundary;
* an unnecessarily large number of MPI ranks is used for a relatively small number of mesh blocks.

For example, on a 4-socket node, `1 × 40` places one MPI rank across all four sockets. A configuration such as `4 × 10` keeps each MPI rank within one socket and generally provides better locality.

## 2. Running MOSE in Parallel

### Single-node execution

A typical full-node configuration on the reference machine uses:

```bash
export OMP_NUM_THREADS=10
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp

mpirun -n 8 -ppn 8 ./bin/MOSE
```

This launches:

```text
8 MPI ranks
×
10 OpenMP threads/rank
=
80 CPU cores
```

The MPI runtime should bind each rank to a distinct group of 10 cores.

For Intel MPI, `I_MPI_PIN_DOMAIN=omp` associates each MPI rank with an OpenMP-sized CPU domain. With OpenMPI, an equivalent configuration is:

```bash
mpirun \
    --map-by socket:PE=$OMP_NUM_THREADS \
    --bind-to core \
    ./bin/MOSE
```

CPU affinity is important: without explicit process and thread binding, different MPI ranks or OpenMP threads may be placed on the same cores, invalidating the intended parallel configuration.

### Multi-node execution

The same strategy extends naturally across nodes.

For example, with 2 nodes containing 4 sockets × 20 cores each:

```text
Node 1: 8 MPI ranks × 10 threads = 80 cores
Node 2: 8 MPI ranks × 10 threads = 80 cores

Total: 16 MPI ranks × 10 threads = 160 cores
```

The important constraint is to preserve the same **per-node socket coverage** as the single-node configuration.

Thus, increasing from 80 to 160 cores should preferably add a complete, identically configured node rather than changing the MPI/OpenMP layout on the original node.

## 3. Verifying CPU Placement

After configuring a parallel run, verify that MPI ranks and OpenMP threads are actually placed as intended.

For Intel MPI, enable placement diagnostics:

```bash
export I_MPI_DEBUG=4
```

The resulting pinning information should show that each MPI rank receives a contiguous set of cores belonging to a single socket.

This check is especially important when moving to a new cluster or scheduler configuration. An incorrect CPU binding may not produce an error, but can significantly reduce parallel performance.

## 4. Measuring Parallel Scaling

MOSE should be benchmarked using its internal timing instrumentation rather than the total wall-clock time of the batch job.

Enable:

```text
timer-diter
```

in the `[MOSE-IO]` configuration and use the reported:

```text
Solver per iteration
```

or the corresponding `wall/iter` timing.

For meaningful scaling measurements:

* run enough iterations to eliminate initialization effects;
* measure the final timing window;
* disable periodic solution output;
* disable residual computation when measuring solver performance;
* repeat each configuration and use the minimum measured time.

These precautions separate the computational scaling of MOSE from setup, I/O, and other external costs.

## 5. Strong Scaling

Strong scaling measures how the execution time changes when the **same problem size** is solved using an increasing number of CPU cores.

The reference benchmark uses a 6.8-million-cell mesh distributed over 24 blocks. The reference system contains 4 sockets × 20 physical cores per node.

The measured parallel efficiency is:

| Cores | Nodes | MPI × OpenMP | Parallel efficiency |
| ----: | ----: | -----------: | ------------------: |
|     4 |     1 |      `4 × 1` |               99.6% |
|     8 |     1 |      `4 × 2` |               99.0% |
|    16 |     1 |      `4 × 4` |               97.7% |
|    20 |     1 |      `4 × 5` |               97.7% |
|    40 |     1 |     `4 × 10` |               94.9% |
|    80 |     1 |     `4 × 20` |                 88% |
|   160 |     2 |     `8 × 20` |                 89% |
|   240 |     3 |    `12 × 20` |                 87% |

<figure>
  {% include "user/images/scaling-strong.svg" %}
</figure>

The results demonstrate that MOSE maintains high parallel efficiency as the problem is distributed across multiple cores and nodes, provided that the MPI/OpenMP configuration preserves socket locality.

In this benchmark, efficiency remains above approximately **95% up to half a node** and around **87–89% from one to three full nodes**.

The particularly important observation is that crossing a node boundary does not introduce a significant additional scaling penalty when each node uses the same well-balanced configuration.

## 6. MPI, OpenMP, and Hybrid Execution

MOSE can be executed using pure OpenMP, pure MPI, or a hybrid MPI/OpenMP configuration.

For the reference benchmark:

| Cores | OpenMP |   MPI | Hybrid |
| ----: | -----: | ----: | -----: |
|     8 |  96.3% | 96.9% |  98.5% |
|    16 |  94.0% | 94.5% |  96.1% |
|    20 |  92.2% | 89.1% |  95.2% |
|    40 |  83.4% | 83.7% |  92.0% |
|    80 |  64.0% | 73.1% |  83.7% |

<figure>
  {% include "user/images/scaling-modes.svg" %}
</figure>

Up to approximately one socket, the three execution models perform similarly. Beyond one socket, the hybrid configuration becomes preferable because it combines MPI-based domain decomposition with OpenMP parallelism inside each socket.

At full-node scale, the hybrid configuration therefore provides the best balance between:

* NUMA locality;
* MPI communication;
* domain-decomposition overhead;
* per-rank memory consumption;
* OpenMP parallelism.

## 7. Scaling Across Nodes

When the problem size is fixed, using additional nodes can increase available memory bandwidth and reduce the number of cores sharing each node's memory subsystem.

However, adding nodes is not automatically beneficial. The per-node MPI/OpenMP configuration must continue to cover the available sockets.

For example, on the reference machine:

```text
80 cores:
    1 node × 8 ranks × 10 threads

160 cores:
    2 nodes × 8 ranks × 10 threads/node

240 cores:
    3 nodes × 8 ranks × 10 threads/node
```

This preserves the same local configuration while increasing the total number of nodes.

For workloads where wall-clock time is the priority, distributing the computation across additional nodes can be advantageous. If allocation cost is important, however, the faster execution may come at the expense of substantially higher node-hours.

## 8. Practical Recommendations

For most production runs, use the following procedure:

1. Determine the number of sockets and physical cores per socket.
2. Choose at least one MPI rank per socket.
3. Choose an OpenMP thread count that fits within one socket.
4. Prefer an MPI rank count that divides the number of mesh blocks.
5. Bind MPI ranks and OpenMP threads explicitly.
6. Verify the resulting CPU placement.
7. Benchmark a representative case before committing to a large production run.
8. When scaling to multiple nodes, preserve the same per-node configuration.
9. Prefer additional OpenMP threads over very large MPI rank counts when memory consumption becomes limiting.

### Recommended default

For a node with **4 sockets × 20 cores**, a good starting point is:

```text
MPI ranks/node     = 8
OpenMP threads     = 10
cores/node         = 80
```

For multiple nodes, replicate this configuration on every node.

## 9. Known Scalability Limits

Parallel scalability is ultimately limited by several factors.

**Memory consumption:** MOSE allocates some domain-wide structures on every MPI rank. Very large MPI rank counts can therefore exhaust memory before all available CPU cores are used. Increasing OpenMP parallelism is preferable when this becomes a limitation.

**Mesh decomposition:** Increasing the number of MPI ranks can increase the number of ghost cells. If the rank count does not divide the number of mesh blocks, additional decomposition overhead can be introduced.

**Serial output:** Solution output is gathered through rank 0 and can become a significant bottleneck for large parallel jobs. Output frequency should therefore be kept low for large-scale runs.

**Problem size:** The measurements presented here are strong-scaling results for a specific 6.8-million-cell case. They should not be interpreted as universal scaling limits for every MOSE application. Larger meshes, different physics, and reactive cases may exhibit different scaling behaviour.

## 10. Quick Reference

| Goal                          | Recommended approach                             |
| ----------------------------- | ------------------------------------------------ |
| Use all sockets               | ≥ 1 MPI rank/socket                              |
| Keep NUMA locality            | Keep each OpenMP team within one socket          |
| Reduce decomposition overhead | Prefer rank counts dividing the mesh block count |
| Scale to multiple nodes       | Replicate the same per-node layout               |
| Avoid CPU oversubscription    | Explicitly bind MPI ranks and OpenMP threads     |
| Reduce memory pressure        | Prefer more OpenMP threads and fewer MPI ranks   |
| Measure solver scaling        | Use MOSE's internal `wall/iter` timing           |
| Validate placement            | Use MPI pinning diagnostics                      |
| Large production runs         | Benchmark the intended MPI/OpenMP layout first   |

**In short:** MOSE scales effectively across cores and nodes when parallelism is organized around the machine's NUMA topology. The recommended strategy is to use MPI to distribute work across sockets and OpenMP to exploit the cores within each socket, while maintaining the same layout on every node.
