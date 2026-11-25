# Euler boundaries

These boundary types do not require a second data line.

| Type | Name | Main behaviour |
|---:|---|---|
| `200` | Axisymmetry / rotational periodicity | Joins the two `k` faces of a block. With one cell in `k`, it defines an axisymmetric wedge; with multiple cells in `k`, it closes a rotationally periodic sector. |
| `300` | Symmetry / slip wall | Zero normal flow through the face. |
| `400` | Extrapolation | Uses the state in the boundary cell to compute the face flux, including when that state moves into the domain. |

## Choosing `400` or a prescribed inlet

Type `400` is suitable for an extrapolating boundary or for continuing a uniform inflow whose state is already in the boundary cell. It is **not** a prescribed inflow: gas that enters uses the boundary-cell state. When the incoming state is known, use an inlet type (`401`–`405`, `408`, or `410`); use ambient type `407` when the surrounding flow can drive gas into the domain.

MOSE reports when gas enters through a `400` face. The detailed log reports the affected block and face and the largest normal velocity relative to sound speed; a final count is printed when such inflow occurred.

## Type `200` and mesh geometry

Type `200` is a geometric connection, not a generic symmetry condition. On faces 5 and 6 of a single-cell-`k` block it creates a thin axisymmetric wedge. With multiple `k` cells it creates a rotationally periodic sector. The paired `k` faces must belong to the same block and both carry type `200`.

A sector or wedge that reaches the axis uses symmetry (`300`) on the axis face. A block closed by `200` must not be split across `k` by the ATLAS mesh decomposer; use its `split-directions` setting to keep the block whole in `k`.

## Directions with one cell

A 1-D or 2-D mesh has one-cell directions. For such a direction, use `0` on both faces when no flux is intended, or use a matched pair of `300` faces where symmetry must keep the transverse velocity at zero. With inclined symmetry faces or an extrapolation (`400`) face, MOSE integrates the transverse velocity; thin cells can then restrict stability.

For the detailed geometry checks, degenerate-direction limits, and type-specific log behaviour, see [Euler boundaries in the complete reference](../boundary-conditions.md#euler-boundaries) and the [boundary-condition theory](../../theory/boundary-conditions.md).
