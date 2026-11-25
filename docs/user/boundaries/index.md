# Boundary conditions

MOSE reads boundary conditions from `INPUT/bc.txt`. Each record assigns a boundary type to one boundary cell on one face of one block. The type determines whether a data line follows the record.

ATLAS can generate the file. Manual preparation is supported, but the block order, face indices, record layout, and any type-specific data must match the mesh and the MOSE format.

## Choose a boundary-condition family

| Purpose | Type codes | Guide |
|---|---|---|
| Join block faces or connect periodic/overset regions | `101`, `102`, `103`, `201` | [Connections](connections.md) |
| Symmetry, axisymmetry, extrapolation | `200`, `300`, `400` | [Euler boundaries](euler.md) |
| Prescribed inflow, back-pressure outlet, ambient, assigned state, choked nozzle | `401`–`410`, `420` | [Inlets and outlets](inlets-outlets.md) |
| Viscous wall conditions | `301`, `302` | [Viscous walls](viscous-walls.md) |
| Lumped plenum joining outlet and inlet faces | `501` | [Manifold](manifold.md) |
| Solid propellant, melting, pyrolysis, surface reactions | `502`–`506` | [Gas–surface interaction](gas-surface.md) |

## File names for multigrid

A multigrid case needs one boundary-condition file per grid level:

| File | Level |
|---|---|
| `INPUT/bc.txt` | 1, finest grid |
| `INPUT/bc2.txt` | 2 |
| `INPUT/bc3.txt` | 3 |
| `INPUT/bc4.txt` | 4 |

The coarse-grid files list the boundary cells on their own grid. See [Multigrid](../multigrid.md).

## Face numbering

For a block with `ni × nj × nk` cells:

| Face | Location | Indices that vary |
|:---:|---|---|
| 1 | `i = 1` (imin) | `j`, `k` |
| 2 | `i = ni` (imax) | `j`, `k` |
| 3 | `j = 1` (jmin) | `i`, `k` |
| 4 | `j = nj` (jmax) | `i`, `k` |
| 5 | `k = 1` (kmin) | `i`, `j` |
| 6 | `k = nk` (kmax) | `i`, `j` |

In 2-D, faces 5 and 6 are degenerate. Use the null type (`0`) or symmetry (`300`), except when the case is an axisymmetric wedge or a rotationally periodic sector using type `200`.

## Record format

The first line of a record identifies the boundary cell and its type:

```text
<block>  <i>  <j>  <k>  <face>  <type>
<data line, when required by the type>
```

Block and cell indices are one-based. A second line is present only for types that require data. Blank lines between records are skipped.

MOSE checks that required values are present, numeric, and finite. Some optional fields have defaults; an empty optional value can retain its default, but an empty required field is an error. The exact fields and defaults are listed in the [complete BC record reference](../boundary-conditions.md).

## Complete reference

The original [complete boundary-condition reference](../boundary-conditions.md) remains the authoritative source for every record layout, field description, special-case check, and file example. The pages in this section are shorter guides to help choose the right family and understand the main behaviour before editing an individual record.
