# Face Gradients and Metrics

This page details how the gradients of the viscous fluxes are computed on a face: the stencil, the face metric, the metric of the
ghost cells, and the treatment of chimera faces. The summary is in [Spatial Discretization](numerics.md#gradient-computation-for-diffusive-fluxes).

## Gradient Computation for Diffusive Fluxes

Velocity, temperature, and species gradients needed for the viscous flux
are computed with a 10-point stencil at each cell interface:

- 2 points in the face-normal direction (the two cells sharing the face)
- 8 points in the two tangential directions (four per direction)

The gradient in computational space $(\xi, \eta, \zeta)$ is mapped to
Cartesian coordinates $(x, y, z)$ with a face metric tensor $M_f$. With
$D\phi = (D_a\phi, D_b\phi, D_c\phi)$ the difference across the face and the
two tangential differences (the mean of the central differences of the two
cells), the face gradient is $\nabla\phi_f = M_f^{T} D\phi$. Each cell
stores its metric $M = A^{-T}$, where the rows of $A$ are the cell edge vectors
$\partial\mathbf{x}/\partial\xi$, $\partial\mathbf{x}/\partial\eta$,
$\partial\mathbf{x}/\partial\zeta$ (each the mean of the four edges of the
cell along that index); the edge lengths also give the physical spacings
$\Delta l$ of the MUSCL reconstruction.

The ghost cells behind a boundary face
take their metric from nodes extrapolated linearly along the index lines,
$\mathbf{N}_{-g} = 2\mathbf{N}_{-g+1} - \mathbf{N}_{-g+2}$ (faces between
blocks: the cell across; symmetry: the mirror image; the wedge: the rotated
cell): both layers have the edge across the face of the first cell, so their
size across the face never vanishes nor turns over however stretched the mesh,
and the first is centred at $2\mathbf{x}_f - \mathbf{x}_1$ ($\mathbf{x}_f$ the
mean of the nodes of the face, $\mathbf{x}_1$ that of the cell).

Their edges
along the face are extrapolated as well, and turn a layer over where they change
sharply across the first cell, as where the index lines across the face spread
apart away from it: on a plane 2-D mesh where they grow 3 times or more the
first layer, from 5/3 times the second. Behind a chimera face the first layer
enters the viscous flux (*Chimera faces*, below). Behind the other faces no flux
reads the metric of either layer: walls take the
metric of the cell, or the stencil of positions, the other faces carry no
viscous flux that needs it, and the wedge rotates the cell.

The solver reads
only the length across the face, in the MUSCL reconstruction, and that is the
length of the first cell whatever the edges along the face: their sign is not
checked. The cells and the first layer of ghost cells behind chimera faces are
checked at setup ([Checks of the cells at set-up](numerics-geometry.md#checks-of-the-cells-at-setup)). The face metric is one of:

| Face metric | | |
|---|---|---|
| stencil (default) | $M_f = \left(J_f^{T}\right)^{-1}$, $J_f^{T} = \left[D_a\mathbf{x},\ D_b\mathbf{x},\ D_c\mathbf{x}\right]$ | the same differences of the cell centres as of the field: exact for a linear field on any mesh (at the symmetry and rotationally periodic faces for the linear fields that they allow), except at the chimera faces, and at the symmetry faces with one cell across them, the solid propellant grain and the ablating gas–surface walls, which take the metric of the cell |
| inverse-mean (not a choice of its own) | $M_f = \left[\tfrac12\left(A_L + A_R\right)\right]^{-T}$ | the face-normal edge is the centre-to-centre vector: exact for a linear field on any stretching normal to the face; taken by the chimera faces always, and with the mean metric by the faces between blocks and the rotationally periodic faces |
| mean | $M_f = \tfrac12\left(M_L + M_R\right)$ | faster; biased on stretched meshes |

With the stencil metric, the columns of $J_f^{T}$ are the differences of the cell
centres (the mean of the eight nodes) that $D_a\phi$, $D_b\phi$ and $D_c\phi$
take of the field. For a linear field $\phi = \phi_0 + \mathbf{g}\cdot\mathbf{x}$
every difference is $D\phi = J_f\,\mathbf{g}$, so the face gradient is
$\mathbf{g}$ on skewed, stretched and non-smooth meshes and across folds at the
faces between blocks: the difference approximations of the field and of the
coordinates must be the same [8]. Vinokur states the rule for the gradient of a
potential (p. 26) and writes the face gradient of the viscous terms on the same
values as here, Eq. (90)–(91), with metrics from averaged areas and volumes,
Eq. (93)–(94), as in the cross-derivative terms of CFL3D [12]; the stencil metric
applies the rule to that stencil.

Its component along the centre-to-centre
vector is the difference of the two cell values over their distance, as in the
face-tangent and node-averaging schemes [9, 10], with the derivatives along the
stencil's tangential directions taken between means of four cells, at the mean
of their centres (near, not in, the plane of the face); on a 2-D mesh of
parallelograms it is the diamond scheme with least-squares weights [11], which
are exactly 1/4 there. Where the tangential differences of the stencil lie in
the plane of the face, its split of the face area vector is the over-relaxed
one of [13]. Exactness for linear fields does not certify accuracy on curved,
highly stretched cells [14]; that case is not measured here.

One effect of the
off-face position of the tangential derivatives is measured (manufactured-solution tests):
on a wedge whose cells grow along $r$ by a ratio $q$, the axial viscous force of
Homann's flow converges at first order on the first rows next to the axis, as
the derivative of $u_r$ along $x$ at the
faces $j$ is that of the mean of the two rows of cells, off the face by about
$\Delta r\,(q-1)/4$. The inverse-mean metric
takes the tangential columns
from the mean edges of the two cells instead.

On smooth meshes they approach
the centre differences as the cells shrink, but where the mesh is not smooth
(cells of alternating size, nodes moved at random, folds at the faces between
blocks) the face gradient of a linear field stays wrong by several per cent or
more at any resolution. The column across the face is kept as $\tfrac12(a_L + a_R)$, which
is $\mathbf{x}_R - \mathbf{x}_L$, and so are the columns of the index
directions with one cell (1-D runs, the $k$ direction of 2-D meshes and
wedges).

On 1-D runs and plane 2-D meshes the field does not vary along them,
and 1-D runs give the bits of the inverse-mean metric. On a wedge the $k$ neighbours are
the cell turned by $\mp\delta$, so the velocity differs across them, and that
column, the mean edge of the cells, is the difference of the turned positions
divided by $\cos^2(\delta/2)$: the error of order $\delta^2$ of the azimuthal
derivatives (see [Axisymmetric wedge](numerics-wedge.md#axisymmetric-wedge)). The ghost cells in the stencil take the
position at which their value is right for a linear field:

- faces between blocks (direct and periodic connections): the centre of the cell across, plus the
  translation of a periodic pair;
- walls (adiabatic or heat-flux, isothermal, gas–surface interaction, and the walls coupled to another
  solver) and extrapolation faces: the position of the class
  of each quantity, $2\mathbf{x}_f - \mathbf{x}_1$ for the quantities with a
  wall value and $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ for the
  extrapolated ones (see [Ghost Cells at Walls and Extrapolation Faces](numerics-ghost-cells.md)); a face
  with such a ghost cell in its stencil takes one metric per group of
  quantities (species and scalars, temperature, velocity, turbulence
  variables);
- faces whose ghost value is $3u_1 - 3u_2 + u_3$ (the velocity inlet and the other boundaries not
  listed here): $3\mathbf{x}_1 - 3\mathbf{x}_2 + \mathbf{x}_3$
  ($2\mathbf{x}_1 - \mathbf{x}_2$, or $\mathbf{x}_1$, in a block with two
  cells, or one, across the face; there the value reads the ghost cell of the
  opposite face, so value and position agree only with three cells or more);
  chimera faces take the same position, an approximation, since their
  values come from donor cells; the tangential differences of their faces take
  as neighbours of a ghost cell along the face the ghost cells of the
  neighbouring faces, and beyond the edges of the block face the row of its
  ghost cells extrapolated linearly (in the k direction of a 2-D mesh, the
  ghost cell turned by the wedge angle); the MUSCL reconstruction of the first
  interior face takes the distance $l_0 = \Delta l_1$ of the ghost cell instead,
  so on cells stretched across these faces the slope behind the face is scaled
  by a factor of order one, which the limiter bounds;
- faces with an assigned state: the file gives the state of the cell across the face,
  taken as the ghost cell, the first cell mirrored through the face; the
  quantities of the file sit at its centre, $2\mathbf{x}_f - \mathbf{x}_1$, and
  the turbulence variables, which the file does not give (the values of the
  first cell per unit mass), at $\mathbf{x}_1$; a face with such a ghost cell in
  its stencil takes one metric per group of quantities, as next to walls. The
  soot, copied in the same way, takes the position of the species: no flux takes
  its gradient. Where the cells of the solver that wrote the file have another
  width across the face, its values sit at half that width from the face, and
  the derivatives across the first row carry the error of that difference;
- zero-gradient faces (unused faces, and the inlets and outlets but the velocity inlet): the cell itself;
- symmetry: the mirror image of the cell across the face; on a face of
  zero area (an axis at $r = 0$) the cell itself, since the ghost value there is
  a copy of the cell, and the viscous flux of that face is zero;
- the rotationally periodic $k$ faces of a wedge or of a 3-D sector: the source of the ghost
  cell, the cell at the other end of the block (on a wedge the cell itself),
  rotated by the angle of the wedge or sector; on a translated sector (no angle)
  the cell at the other end translated by the period, as across a periodic pair: plus the
  vector from the centre of its outer face to the centre of this face. The
  mesh check at setup has found every node of the last $k$ plane of the block
  to be the node of the first plane translated by one vector, within 1 % of
  the move of the nodes per cell in $k$ and better than turned (the round-off
  of the coordinates, in a mesh written with enough digits), so that vector is
  the period of the block to the same tolerance.

Behind a symmetry face, a periodic pair and a rotational periodic face the
gradient is exact for the linear fields that the condition allows:
$\mathbf{g}\cdot\mathbf{n} = 0$ for a scalar and $G = RGR$, $R\mathbf{u}_0 =
\mathbf{u}_0$ for the velocity at a symmetry plane ($R = I - 2\mathbf{n}\mathbf{n}^{T}$),
$\mathbf{g}\cdot\mathbf{T} = 0$ across a periodic pair of translation $\mathbf{T}$
and across the $k$ faces of a translated sector of period $\mathbf{T}$,
$\mathbf{u}(R_\delta\mathbf{x}) = R_\delta\mathbf{u}(\mathbf{x})$ across a rotational
periodic face.

For a field that breaks the condition the value of the ghost cell
(reflected, copied or rotated) is not the value of the field at any point, so no
position makes the stencil exact: the error is proportional to the part that the
condition forbids ($\mathbf{g}\cdot\mathbf{n}$, $G - RGR$, $\mathbf{g}\cdot\mathbf{T}$),
it stays on the faces next to the boundary, and it vanishes as the solution
satisfies the condition. Behind the other faces between blocks and the
zero-gradient faces the gradient is exact for every linear field; behind walls
and extrapolation faces for the linear fields of the next section; behind mapped faces
for the linear fields whose values in the file are those at the centre of the
ghost cell.

Behind the faces whose ghost value is $3u_1 - 3u_2 + u_3$ it is
exact for every linear velocity, but for the mass fractions, the temperature and
the specific turbulence and scalar variables only where the quantities that the
rule extrapolates, $\rho_s$, $p$ and $(\rho\phi)$, are linear too (for $Y_s$ and
$\phi$ where the density is uniform, for $T$ where the composition is uniform
as well): their quotients $Y_s = \rho_s/\rho$, $T = p/(\rho R)$ and
$\phi = (\rho\phi)/\rho$ are not linear where the density varies.

There the
derivative across the first row has a relative error of second order in the
cell size on uniform cells and of first order on cells that grow by a fixed
ratio across the face.

$M_f$ is computed once at setup and stored on every face of the blocks of each
rank, 27 values per cell (the mesh does not move); the cell centres come from
the nodes, which every rank holds, so no message is needed. Where a face or a
cell is degenerate, the differences of its stencil can stop spanning the space
and the metric is not defined there: the run stops at setup, a declared limit of the stencil metric
(the mean metric does not need it, but a viscous run with it stops on a
3-D mesh at a boundary face of zero area, below).

A wall on faces of zero area,
as on
the axis of a 3-D mesh around it (a polar line), stops too: its gradient is not
defined there, and the axis takes the symmetry condition. With symmetry on
the axis the polar line runs, and the face gradients next to it are exact for
linear fields. With the mean metric a viscous run on a 3-D mesh with a
boundary face of zero area, as on a polar line with the axis exactly at $r = 0$,
stops at setup: next to the axis
the viscous force of the first row of cells grows as the mesh is refined (a
linear field, whose force is zero).

With the axis nodes exactly at $r = 0$ the polar line runs with the stencil metric, but
the solution next to the axis is of first order. On smooth fields the residual
of the two rows of cells next to the axis does not converge as the mesh is
refined: its viscous part grows as $1/h$ and its convective part does not
decrease. With the mean metric the viscous residual of the cells next
to a symmetry plane that bounds a sector of a polar mesh does not
converge either: on smooth fields that the plane allows it stays of order one as
the mesh is refined, while with the stencil metric those cells converge as the others.

In one dimension $M = 1/h$, so the mean metric divides the difference across the face by
the mean of $1/h$ instead of the centre-to-centre distance. With a cell-size
ratio $r = h_R/h_L$ it overstates the face-normal gradient, and so the viscous,
conductive and diffusive fluxes, by

$$
\frac{(1+r)^2}{4r}:\qquad
+0.2\,\% \;(r=1.1),\quad +0.8\,\% \;(r=1.2),\quad +4.2\,\% \;(r=1.5),\quad +8.9\,\% \;(r=1.8)
$$

This bias does not decrease with grid refinement at a fixed growth ratio, and it
is present in the converged steady solution. The mean metric is adequate where adjacent
cells grow by about 1.2 or less; boundary-layer meshes often grow faster near
the wall.

With the stencil metric, a face between blocks takes one face gradient with its stencil
metric, as an interior face, and a symmetry face the gradient between the cell
and its mirror image. With the mean metric, block-connection faces
use the inverse-mean metric, as chimera faces do always: across a large
cell-size jump at a block interface, the mean makes the normal gradient so steep
that the explicit diffusion becomes unstable.

At a wall face (adiabatic or heat-flux, isothermal, and gas–surface interaction walls when they do not
ablate) the gradient takes the difference between the cell and the wall value across the face.
With the stencil metric it also takes the differences of the wall values along the
wall, between the neighbouring wall faces of the same block face and type
(central, or one-sided where the wall ends or goes on in another block face,
whose faces it does not take: exact there for a linear wall value, first order
for a curved one), with the metric of the same differences of the positions:
the vector from the face centre (the mean of its four nodes) to the cell centre
across the face, the differences of the face centres along it.

For a linear
field whose wall values are its values at the face centres the wall gradient
is exact: the temperature of an isothermal wall that
varies along the wall, the velocity of a no-slip wall on cells of any shape.
The wall temperatures of two neighbouring records are differenced as they are
given: a jump between them enters the derivative along the wall at the faces
next to it, of the order of the jump over the spacing of the faces, and on
cells skewed by θ their heat flux through $\tan\theta$.

The difference across the face is that of a wall of a finite-volume grid
[8, Eq. (153)]; the differences along the wall and the metric of the same
positions are MOSE's own. On a 2-D mesh the wall stencil has the structure of
the boundary diamond of [11] (a value at the midpoint of the boundary edge and
values along the edge), there at the vertices of the edge and here at the
neighbouring face centres.

A boundary face that keeps only the normal distance,
as in [13] (Eq. 3.57–3.59) and in OpenFOAM (OpenFOAM-dev, commit dfa6c7e,
`fvPatch.C`), has the error of the mean metric given below.
With the mean metric, and for the quantities whose wall value is not
known along the wall (the temperature of adiabatic or heat-flux walls and of the gas–surface
interaction models), the derivatives along the wall are zero and the metric is that of
the wall cell: on a cell skewed by θ the normal gradient is then wrong by
$(\partial\phi_w/\partial s)\tan\theta$, which refinement does not reduce.

The turbulence variables take no difference of their wall values along the
wall, with every metric: their wall values are uniform along a smooth wall
(zero for $\tilde\nu$ of Spalart–Allmaras, for $k$ and for the Reynolds
stresses), but not those of $\omega$ of the $k$–$\omega$ models and of
SSG/LRR, $C\mu_l/d^2$ with $d$ the distance of the cell centre from
the wall, which follow the height of the first cell, nor that of $\tilde\nu$
on a rough wall, taken from the value of the cell. Their derivative along the
wall is left out, with the error above on skewed cells.

A no-slip wall at rest in the laboratory, in a rotating frame, moves in the frame
with $\mathbf{u}_w = -\boldsymbol{\Omega}\times(\mathbf{x}-\mathbf{x}_o)$
($\mathbf{x}_o$ the origin of the frame), a velocity linear in $\mathbf{x}$
with the constant gradient $-[\boldsymbol{\Omega}\times]$. With every metric
the wall gradient takes the differences of $\mathbf{u}_w$ along the position
differences of its metric, the columns of $M^{-1}$: those of the wall stencil
above with the stencil metric, the edges of the wall cell with
the mean metric.

The part of $\mathbf{u}_w$ in the wall gradient is then exact with
the stencil metric on cells of any shape, and with the mean metric on
parallelepiped cells, where the metric of the wall cell is exact; the rest is
that of a wall at rest. Without these differences the wall stress would be
wrong by the order of $\mu\,|\boldsymbol{\Omega}|$ with every metric and on
any mesh, orthogonal cells included, because the stress takes the derivatives
of the velocity along the wall: a gas at rest in the absolute frame, a rigid
rotation in the rotating frame, would feel a shear on the wall.

The part of
$\mathbf{u}_w$ is exact because the values sit where the metric places them:
the value of a cell at the centre of the cell (the mean of its eight nodes),
the wall velocity at the centre of the wall face (the mean of its four nodes).
A change of either position has to be made in the metric too.

**Chimera faces.** The first layer of ghost cells behind a chimera face is the
ghost cell of the face, whose metric enters its viscous flux: a first layer that
turns over (index lines that spread apart sharply across the first cell, above)
or whose determinant is not that of the first cell in sign and to $10^{-12}$ in
size stops the run at setup; the remedy is the
mesh, with index lines that change smoothly across the first two node layers at
that face.

The ghost values of a chimera face come from the donor cells and weights that the overset
pre-processor computed for its own ghost cells; the face metric of their
viscous flux takes these linear ones, so where the pre-processor builds its ghost cells
otherwise (on cells stretched across the face) the normal gradient is taken over
a distance that differs by as much (no test with chimera faces).
The viscous flux of a chimera face is the mean of two fluxes: one with the
central differences along the face of the row of the cell, one with those of
the row of the ghost cells, filled from the donors.

With a stress linear in the
gradient (laminar, Spalart–Allmaras) the mean is the flux of the mean gradient,
whose differences along the face are the mean of those of the two rows, as on
an interior face (Vinokur 1986, NASA CR-177416, eqs. 90–91, p. 23). Beyond the
edges of the block face the row of the ghost cells is extrapolated linearly, in
the variables whose differences the stencil takes: the difference there is the
one-sided one, exact for a linear field (a copy of the ghost cell would halve
it, and the face would take 3/4 of the derivative of a linear field along the
edge).

In the k direction of a 2-D mesh the neighbours of a ghost cell are the
ghost cell turned by the wedge angle, as those of a cell are; across a block
one cell wide in another direction they are the ghost cell itself.

!!! note "Tangential derivatives"
    With the mean metric the tangential derivatives on a face are
    central differences over two cells, divided by the cell's own tangential
    edge. With stretching in a tangential direction they carry the same
    $(1+r)^2/(4r)$ factor, which the mean metric does not correct. In a thin boundary
    layer this enters only through the streamwise variation of the shear on
    faces normal to the flow, so its effect is much smaller than that of the
    face-normal gradient. The stencil metric divides them by the difference of the
    centres they span, and has no such factor.

---

## References

[8] M. Vinokur, "An analysis of finite-difference and finite-volume formulations of conservation laws," NASA CR-177416, 1986: pp. 15–16, Eq. (60a) (the area vector of a face of a hexahedron as half the cross product of its diagonals; the shape of each face "precisely defined and consistently used by the two neighboring cells" so that the volumes sum to the total volume; the doubly ruled face); p. 18 (the value of a cell taken as the value at some average point of the cell to compute a flux; the centroid and the vectorial average of the vertices); p. 26, Eq. (101) ("the two difference approximations must be the same"); p. 23, Eq. (90)–(91) (the transport terms on an auxiliary cell, with means of four cells along the face); p. 40, Eq. (153) (the difference across a wall). Journal version: *J. Comput. Phys.*, 81, 1989.

[9] B. Diskin, J. L. Thomas, E. J. Nielsen, H. Nishikawa, J. A. White, "Comparison of Node-Centered and Cell-Centered Unstructured Finite-Volume Discretizations: Viscous Fluxes," *AIAA J.*, 48(7), 1326–1338, 2010, doi:10.2514/1.44940, Sec. III, III.A and V.A.

[10] H. Nishikawa, "Beyond Interface Gradient: A General Principle for Constructing Diffusion Schemes," AIAA Paper 2010-5093, 2010, App. C.

[11] Y. Coudière, J.-P. Vila, P. Villedieu, "Convergence rate of a finite volume scheme for a two dimensional convection-diffusion problem," *ESAIM: M2AN*, 33(3), 493–516, 1999, Eq. (9) and Lemma 6.2.

[12] CFL3D, `source/cfl3d/libs/gfluxv1.F` (cross-derivative terms of the viscous fluxes), https://github.com/nasa/CFL3D, commit 776d516; S. L. Krist, R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual (Version 5.0)*, NASA/TM-1998-208444, 1998, App. C.2, p. 249 (the volume at a cell interface as the mean of its neighbours, for the thin-layer terms; the manual, p. 27, leaves out the cross-derivative terms, which gfluxv1.F adds).

[13] H. Jasak, *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*, PhD thesis, Imperial College, University of London, 1996, Eq. (3.1), (3.12) and (3.32).

[14] D. J. Mavriplis, "Revisiting the Least-squares Procedure for Gradient Reconstruction on Unstructured Meshes," NASA/CR-2003-212683, 2003, p. 14.
