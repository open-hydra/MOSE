# Mesh Geometry and Block Interfaces

This page details the geometry behind [Spatial Discretization](numerics.md): the volume of a cell, the checks of the cells at set-up,
and the flux of a face between two blocks.

## Cell volume

A cell is the image of the unit cube by the trilinear map
$\mathbf{x}(\xi,\eta,\zeta)$ of its eight nodes. Each face is the bilinear
(doubly ruled) surface through its four nodes, the same for the two cells that
share it, and its area vector, the one the fluxes take, is half the cross
product of its diagonals [8, pp. 15–16]. The volume of the cell is the integral
of the Jacobian of the map,

$$
V = \int_0^1\!\!\int_0^1\!\!\int_0^1 \det J\,\mathrm{d}\xi\,\mathrm{d}\eta\,\mathrm{d}\zeta
  = \frac{1}{6}\sum_{d=1}^{3}\bar{\mathbf{e}}_d\cdot\left(\mathbf{S}_d^{-}+\mathbf{S}_d^{+}\right),
$$

with $\bar{\mathbf{e}}_d$ the mean of the four edges of the cell along index $d$
and $\mathbf{S}_d^{\mp}$ the area vectors of the two faces across index $d$,
both along $+d$: an identity of polynomials in the 24 coordinates, which can be
checked term by term ($\det J$ has degree at most two in each of $\xi$, $\eta$,
$\zeta$). It is also
$\frac13\sum_f\bar{\mathbf{x}}_f\cdot\mathbf{S}_f$ over the six faces, outward,
with $\bar{\mathbf{x}}_f$ the mean of the nodes of the face: the divergence
theorem for $\nabla\cdot\mathbf{x} = 3$, with
$\int_f\mathbf{x}\cdot\mathrm{d}\mathbf{S} = \bar{\mathbf{x}}_f\cdot\mathbf{S}_f$ on a
bilinear face.

So the discrete divergence of a field that grows as $\mathbf{x}$
is exact on every cell and, as two neighbours share their face, the volumes of
the cells sum to the volume enclosed by the boundary faces, with no overlap and
no gap [8, p. 16]. On a cell with plane faces (a 2-D mesh, the wedge, a 3-D
mesh of plane faces) it is the volume of the polyhedron. A sum of tetrahedra,
each face cut along a diagonal, gives that volume only where the faces are
plane: the two cells that share a twisted face cut it along different
diagonals, and each volume misses or counts twice
$\mathbf{T}\cdot\mathbf{S}_f/12$, with
$\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$ the
twist of the face.

The volume has the sign of the block (negative on a left-handed cell) and is
stored in absolute value. A cell whose volume has not the sign of the
determinant at its centre is turned over somewhere and stops the run at setup
(below). The volume is computed from differences of the nodes, so its rounding
is that of the size of the cell, not that of its coordinates; the compiler must
keep the order of these differences.

## Checks of the cells at setup

A cell whose determinant is zero to $10^{-12}$ of the product of its edge lengths
stops the run at setup. So does a cell folded at
one of its corners, or turned against its block: a cell is the image of the unit cube by the trilinear map of its eight
nodes, and at a corner the Jacobian of that map is
$J = \mathbf{e}_i\cdot(\mathbf{e}_j\times\mathbf{e}_k)$, with
$\mathbf{e}_i, \mathbf{e}_j, \mathbf{e}_k$ the three edges from that corner along
$i$, $j$ and $k$.

At every corner it must have the sign of the block (the sign of
the determinant of the first cell, from which the face normals of the block take
their orientation) and exceed $10^{-12}\,|\mathbf{e}_i||\mathbf{e}_j||\mathbf{e}_k|$;
a corner with an edge of zero length (nodes that coincide on an axis) is not
checked. On a plane 2-D mesh $J = \Delta z\, j(\xi, \eta)$, with $j$ the Jacobian
of the bilinear map of the section, $\mathbf{x} = \mathbf{a} + \mathbf{b}\xi +
\mathbf{c}\eta + \mathbf{d}\xi\eta$: $j = (\mathbf{b} + \mathbf{d}\eta)\times
(\mathbf{c} + \mathbf{d}\xi)$ is affine in $(\xi, \eta)$ because
$\mathbf{d}\times\mathbf{d} = 0$, so its sign in the cell is decided at the
corners.

On the 2-D wedge $J = R\,\sin\delta\; j_m(\xi, \eta)$, with $R > 0$ the
bilinear radius when no node is on the axis and $j_m$ the Jacobian of the
meridional section, affine as above: the corners decide again. In 3-D the corners
are a necessary condition only; a cell whose volume (above) has not the
sign of the determinant at its centre is turned over somewhere, and stops the
run at setup as well ($\det J$ has the sign of the block in every point of a
cell that is not turned over, and so has its integral).

These checks run on every
multigrid level too, whose nodes are every other node of the level above. Behind a face of zero
area (an axis) the ghost cells are the reflection of the first cell through the
axis, whose determinant has the opposite sign by construction: their sign is
not checked (a wall on such a face stops, see *Gradient Computation for
Diffusive Fluxes*). The first layer of ghost cells behind a chimera face is
checked too (see *Chimera faces* in that section).

---

## Faces Between Blocks

A face shared by two blocks (a block connection or a periodic connection)
has one flux, computed once and given to both cells with opposite signs, as an
interior face. Of the two boundary records of the face, the one of the cell on
the low-index side (whose face is 2, 4 or 6 while the other is 1, 3 or 5)
computes the inviscid and the viscous flux with the shock sensor of its cell,
as the interior face between cells $n$ and $n+1$ uses that of cell $n$, and with
the preconditioning reference velocity $U_r$ of both cells (the ghost cell
carries that of the other block).

When the two faces have the same parity, as
where blocks meet with their index directions permuted, the record with the
smaller (block, face, $k$, $j$, $i$) computes it, so that the result does not
depend on the order in which the boundary faces are listed. The face therefore conserves
mass, momentum and energy exactly whatever the Riemann solver, the shock sensor
and the preconditioning, on any number of MPI ranks; between ranks the other
cell's share travels in one more message per Runge–Kutta stage.

The inviscid flux of a face between blocks is, bit for bit, that of the same
face in the unsplit mesh when three conditions hold: the index directions of
the two blocks agree across the face; the two blocks store the nodes of the
common face with the same bits (the ghost cell takes its geometry from the
nodes of the other block, so nodes that differ in the last digits give a
slightly different face); and $U_r$ is not smoothed.

This is an
identity of the fluxes, not of the solution: in 2-D and 3-D a cell next to the
cut adds the flux of the face between blocks after those of its interior faces,
so its sum is taken in another order than in the unsplit mesh and the two
solutions differ at round-off. Only in 1-D, where a cell sums two face fluxes,
is the solution of an inviscid run on the split mesh the same bit for bit (the
viscous flux of a face between blocks has its own stencil, below).

With smoothing passes of $U_r$ the max-dilation of $U_r$ in a
cell next to a face between blocks also sees the cell across the face, through
the ghost cell, with its value at the start of the stage. In the unsplit mesh
the smoothing sees the new value of the neighbour, so with the smoothing on a
split mesh does not give the fluxes of the unsplit one bit for bit.

With the stencil metric the viscous flux of a face between blocks is
that of an interior face: one face gradient, from the cell, the ghost cell and
their tangential neighbours (those of the ghost cell are cells of the other
block), with the stencil face metric. With
the mean metric it uses its own stencil (two gradients, each with the tangential
differences of one of the two cells, and the mean of their fluxes) and the
inverse-mean face metric.

The ghost cells and the neighbours along the face of those of the faces
between blocks are filled at every Runge–Kutta stage, before the
fluxes, from the state of that stage, in an order such that every value the
fluxes read is that of the stage, with any number of threads and ranks, under
five conditions:

1. the source block of every face between blocks has at least two cells across
   that face;
2. every block with a face whose ghost values are extrapolated
   ($3u_1 - 3u_2 + u_3$) has at least three cells across it;
3. every cell face of the boundary has its own boundary record;
4. there is no chimera face nor face with an assigned state;
5. the tangential index steps of every block and periodic connection are unit steps
   along one index.

The first is checked at setup: a connection whose source block has one cell
across its face stops the run. The others are not checked: where one does not
hold, a ghost value can be that of the previous stage, or depend on the order
of the records or of the threads.

---

## References

[8] M. Vinokur, "An analysis of finite-difference and finite-volume formulations of conservation laws," NASA CR-177416, 1986: pp. 15–16, Eq. (60a) (the area vector of a face of a hexahedron as half the cross product of its diagonals; the shape of each face "precisely defined and consistently used by the two neighboring cells" so that the volumes sum to the total volume; the doubly ruled face); p. 18 (the value of a cell taken as the value at some average point of the cell to compute a flux; the centroid and the vectorial average of the vertices); p. 26, Eq. (101) ("the two difference approximations must be the same"); p. 23, Eq. (90)–(91) (the transport terms on an auxiliary cell, with means of four cells along the face); p. 40, Eq. (153) (the difference across a wall). Journal version: *J. Comput. Phys.*, 81, 1989.
