# User Guide

This guide is organized around the jobs most MOSE users need to do: prepare a case, run it, inspect the output, and continue from a saved solution.

## Choose a task

| I want to… | Start here |
|---|---|
| Prepare the mesh and starting flow field | [Initial conditions](initial-conditions.md) |
| Set up boundaries | [Boundary conditions](boundary-conditions.md), or [by family](../user/boundaries/index.md) |
| Understand or edit `input.ini` | [Input file essentials](input-file.md) |
| Turn on a model or numerical option | [Common switches](switches.md) |
| Run a prepared case | [Run a case](run-a-case.md) |
| Use MPI and OpenMP | [Parallel execution](parallel.md) |
| Find out why a run stopped or looks wrong | [Diagnostics & troubleshooting](../user/diagnostics.md) |
| Continue a simulation from a saved solution | [Restart a simulation](restart.md) |
| Use multigrid or multirate stepping | [Multigrid & FAS](../user/multigrid.md), [Multirate](../user/multirate.md) |

New to MOSE? Start with [Getting Started](../getting-started/index.md), then follow the [Quick Start](../getting-started/quick-start.md).

## The case directory

Run MOSE from the directory that contains `input.ini`. A typical case looks like this:

```text
my-case/
├── input.ini
├── INPUT/
│   ├── ic.tec
│   ├── bc.txt
│   └── thermodynamic and model data
└── OUTPUT/                 created by MOSE
```

The initial-condition file contains both the mesh coordinates and the cell-centred flow field. MOSE does not read a separate mesh file.

ATLAS can generate the input files, but it is not required. You can also create them manually or with your own preprocessing scripts, provided they follow the documented formats.

## Detailed reference

The task pages keep the common workflow short. The original reference pages remain available for complete parameter tables, output diagnostics, file layouts, and boundary-condition records:

- [Full running guide](../user/using.md)
- [Full input-file guide](../user/input.md)
- [Parameter registry](../user/registry.md)
- [Initial-condition file format](../user/initial-conditions.md)
- [Restart details and limitations](../user/restart.md)
- [Parallel execution details and benchmarks](../user/parallel.md)
- [Boundary-condition file reference](../user/boundary-conditions.md)
