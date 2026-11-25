# Initial conditions

MOSE starts a fresh run from a file in `INPUT/`. The file contains both the mesh node coordinates and the cell-centred flow variables for every block; no separate mesh file is read.

ATLAS commonly generates this file, but ATLAS is optional. You can prepare the file yourself or use another preprocessing tool if it writes a supported format.

## File name and format

The `ic-format` setting in `[MOSE-IO]` selects the reader. The expected file name is:

| `ic-format` | File |
|---|---|
| `tecplot ascii` | `INPUT/ic.tec` |
| `tecplot binary` | `INPUT/ic.szplt` |
| `vtk ascii` or `vtk raw` | `INPUT/ic.vtm` and the `.vts` block files it references |

For VTK, keep the `.vtm` file and its referenced block files together. The references in the `.vtm` are resolved relative to that file.

## Required flow variables

Every initial-condition file must provide:

| Variable | Meaning |
|---|---|
| `x, y, z` | Mesh node coordinates |
| `rho(i)` | Partial density of each species |
| `u, v, w` | Velocity components |
| `p` | Pressure |

The selected turbulence model adds variables:

| Model | Additional variables |
|---|---|
| Spalart–Allmaras | `rhoNuTilde` |
| k–omega or SST | `rhoK`, `rhoOmega` |
| SSG/LRR Reynolds-stress | `ru'u'`, `rv'v'`, `rw'w'`, `ru'v'`, `ru'w'`, `rv'w'`, `rhoOmega` |

The turbulence variables are density-weighted quantities. For example, `rhoNuTilde` is not the eddy viscosity `mut`.

MOSE reads variables by position. For Tecplot ASCII files, it also checks the turbulence-variable names. All coordinates and flow values that MOSE reads must be finite; a NaN or infinity stops setup.

## Starting from a previous solution

There are two different operations:

- **Restart:** set `newrun = false`; MOSE resumes from a saved solution in `OUTPUT/` and continues its time/iteration state. See [Restart a simulation](restart.md).
- **New run from an old field:** copy the solution into `INPUT/ic.*`, set `ic-format` to match it, and keep `newrun = true`. MOSE uses the field but starts time, iteration, and run counters from zero.

For the complete file layout, legacy variable names, and format-specific details, see the [initial-condition reference](../user/initial-conditions.md).
