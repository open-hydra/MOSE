# Block and periodic connections

Connection types join mesh faces or transfer a state between blocks. They are not ordinary external boundaries.

| Type | Meaning | Data line |
|---:|---|---|
| `0` | Null / unused face | None |
| `101` | Block connection | Target block, cell, face, and index-direction mapping |
| `201` | Periodic connection by translation | Same record structure as `101` |
| `102` | Chimera / overset connection | Donor counts followed by donor block, cell, and volume-fraction records |
| `103` | Coupled wall | Target block, cell, face, and index-direction mapping |

For `101` and `201`, the two faces must point to each other: each record identifies the cell and face on the opposite side. The index-direction cosine matrix maps source indices to target indices.

A block connected across a face must have at least two cells across that face at every multigrid level. MOSE uses those cells for the two ghost-cell layers; a block only one cell thick across a connected face is not supported.

The face flux is computed once and applied with opposite signs to the two adjoining cells, including when the blocks are assigned to different MPI ranks. The result does not depend on the order of the records in `bc.txt`.

For the exact data-line order, donor record syntax, and validation messages, see the [connection records in the complete reference](../boundary-conditions.md#connection). For block decomposition and MPI assignment, see [Parallel execution](../parallel.md#block-to-rank-assignment).
