# Spatial Discretization

MOSE uses a cell-centred finite volume method (FVM) on structured
multi-block grids.  This page describes the discretization framework,
the MUSCL reconstruction that provides second-order accuracy, the
slope limiters that ensure monotonicity, and the shock-detection
sensor that adaptively reduces the order near discontinuities.

---

## Finite Volume Framework

The integral form of the conservation law over a cell $\Omega_i$ with
boundary $\partial\Omega_i$ reads

$$
\frac{\mathrm{d}}{\mathrm{d}t}
\int_{\Omega_i}\!\mathbf{U}\,\mathrm{d}V
+ \oint_{\partial\Omega_i}\!\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
= \oint_{\partial\Omega_i}\!\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}}\,\mathrm{d}A
+ \int_{\Omega_i}\!\mathbf{S}\,\mathrm{d}V
$$

Approximating the integrals with midpoint quadrature yields the
semi-discrete update

$$
\frac{\mathrm{d}\mathbf{U}_i}{\mathrm{d}t}
= -\frac{1}{V_i}\sum_{f}\bigl(\mathbf{F}^{c}_f
  - \mathbf{F}^{v}_f\bigr)\,A_f
+ \mathbf{S}_i
$$

where $V_i$ is the cell volume and the sum runs over the six faces $f$
of each hexahedral cell. The numerical algorithm consists of the following steps:  
- Convective fluxes: left and right states at each face are
  reconstructed, then an approximate (or exact) Riemann solver returns
  the numerical convective flux $\mathbf{F}^{c}_f$ based on these states and the face normal vector (see [Riemann Solvers](riemann-solvers.md)).  
- Diffusive fluxes: velocity, temperature, and species gradients are
  computed at each face with a 10-point stencil and mapped to Cartesian
  coordinates using the face metric tensor.

Faces of type `200` join the two $k$ faces of one block, turned about the $x$
axis or translated (see
[Rotational Periodicity and Axisymmetric Wedge](#rotational-periodicity-and-axisymmetric-wedge-type-200)).

### Cell volume

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
bilinear face. So the discrete divergence of a field that grows as $\mathbf{x}$
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
keep the order of these differences (test `VolatileOrder`).

### Checks of the cells at setup

A cell whose determinant is zero to $10^{-12}$ of the product of its edge lengths
stops the run at setup, with its block and cell. So does a cell folded at
one of its corners, or turned against its block, with its block, cell and corner
node: a cell is the image of the unit cube by the trilinear map of its eight
nodes, and at a corner the Jacobian of that map is
$J = \mathbf{e}_i\cdot(\mathbf{e}_j\times\mathbf{e}_k)$, with
$\mathbf{e}_i, \mathbf{e}_j, \mathbf{e}_k$ the three edges from that corner along
$i$, $j$ and $k$. At every corner it must have the sign of the block (the sign of
the determinant of the first cell, from which the face normals of the block take
their orientation) and exceed $10^{-12}\,|\mathbf{e}_i||\mathbf{e}_j||\mathbf{e}_k|$;
a corner with an edge of zero length (nodes that coincide on an axis) is not
checked. On a plane 2-D mesh $J = \Delta z\, j(\xi, \eta)$, with $j$ the Jacobian
of the bilinear map of the section, $\mathbf{x} = \mathbf{a} + \mathbf{b}\xi +
\mathbf{c}\eta + \mathbf{d}\xi\eta$: $j = (\mathbf{b} + \mathbf{d}\eta)\times
(\mathbf{c} + \mathbf{d}\xi)$ is affine in $(\xi, \eta)$ because
$\mathbf{d}\times\mathbf{d} = 0$, so its sign in the cell is decided at the
corners. On the 2-D wedge $J = R\,\sin\delta\; j_m(\xi, \eta)$, with $R > 0$ the
bilinear radius when no node is on the axis and $j_m$ the Jacobian of the
meridional section, affine as above: the corners decide again. In 3-D the corners
are a necessary condition only; a cell whose volume (above) has not the
sign of the determinant at its centre is turned over somewhere, and stops the
run at setup as well ($\det J$ has the sign of the block in every point of a
cell that is not turned over, and so has its integral). These checks run on every
multigrid level too, whose nodes are every other node of the level above: the stop of
the corners names the level, the other stops give the block and the cell of that
level. Behind a face of zero
area (an axis) the ghost cells are the reflection of the first cell through the
axis, whose determinant has the opposite sign by construction: their sign is
not checked (a wall on such a face stops, see *Gradient Computation for
Diffusive Fluxes*). The first layer of ghost cells behind a chimera face is
checked too (see *Chimera faces* in that section).

---

## MUSCL Reconstruction

First- and second-order spatial accuracy is available. The second-order accuracy is achieved through piecewise-linear reconstruction (MUSCL — Monotone Upstream-centred Schemes for Conservation Laws).

At each cell interface the algorithm uses a four-point stencil
$\{i{-}1,\;i,\;i{+}1,\;i{+}2\}$ and proceeds as follows:

1. Compute slopes from consecutive cell values and physical spacings
   $\Delta l_0,\,\Delta l_1,\,\Delta l_2$:

$$
s_0 = \frac{P_i - P_{i-1}}{\Delta l_0},\qquad
s_1 = \frac{P_{i+1} - P_i}{\Delta l_1},\qquad
s_2 = \frac{P_{i+2} - P_{i+1}}{\Delta l_2}
$$

2. Limit the slopes to prevent spurious oscillations near
   discontinuities:

$$
\bar{s}_L = \phi(s_1,\, s_0),\qquad
\bar{s}_R = \phi(s_2,\, s_1)
$$

   where $\phi$ is the chosen limiter function (see below).

3. Reconstruct interface values:

$$
P_L = P_i     + \beta\,\bar{s}_L\,\delta l_L,\qquad
P_R = P_{i+1} - \beta\,\bar{s}_R\,\delta l_R
$$

   where $\delta l_L,\,\delta l_R$ are the distances from the cell
   centres to the interface and $\beta$ is a blending parameter
   controlled by the shock detector.

!!! note "Role of β"

    $\beta = 1$ recovers the full second-order MUSCL scheme.
    $\beta = 0$ reduces reconstruction to first order (donor-cell),
    which is used in the immediate vicinity of shocks.
    The shock-detection sensor (see below) smoothly modulates $\beta$.

At the first interior face the slope behind the first cell, $s_0$ or $s_2$,
takes the ghost cell of the boundary face. Behind a zero-gradient face (`0`,
`401`–`407`, `420`) the ghost value is that of the first cell: that slope is
zero, and every limiter gives zero with a zero argument, so the reconstruction
of the first interior face is of first order on the side of the first cell.
Behind a mapped face (`410`) the ghost cell is the first cell mirrored through
the face (see *Gradient Computation for Diffusive Fluxes*): it has the length
of the first cell across the face, which is then the spacing $\Delta l_0$ (or
$\Delta l_2$) of the slope behind the first cell, and its centre lies that far
from the centre of the first cell. Along a straight grid line the slope is then
that of the field for the quantities of the file; for the turbulence variables
and the soot, whose specific value is that of the first cell, the slope of
$(\rho\phi)$ is $\phi_1$ times that of the density.

!!! note "Choosing the reconstruction"

    `space-reconstruction = MUSCL` limits the slopes with the function
    named by `flux-limiter` (see [Flux Limiters](#flux-limiters)); if
    `flux-limiter` is omitted, van Leer is used and a warning is printed.
    `space-reconstruction = first-order` sets every slope to zero
    (donor-cell); a `flux-limiter` given with it is ignored, with a warning.

---

## Flux Limiters

Several limiter functions $\phi(a,b)$ are available.  All are expressed in
terms of the slope ratio $r = b\,/\,a$.

**MinMod**

$$
\phi(a,b) =
\begin{cases}
\operatorname{sign}(a)\,\min(|a|,\,|b|) & \text{if } a\,b > 0 \\
0 & \text{otherwise}
\end{cases}
$$

Most diffusive TVD limiter.  Strictly enforces monotonicity.

**Van Leer**

$$
\phi(a,b) = \frac{r + |r|}{1 + |r|}\;a
$$

Smooth, differentiable limiter.  Good balance of dissipation and accuracy.

**Van Albada**

$$
\phi(a,b) = \frac{r^2 + r}{1 + r^2}\;a
$$

Smooth, less dissipative than van Leer near local extrema.

**MC (Monotonized Central)**

$$
\phi(a,b) =
\max\!\Bigl(0,\;\min\!\bigl(2r,\;\tfrac{1}{2}(1+r),\;2\bigr)\Bigr)\;a
$$

Symmetric TVD limiter that permits the steepest gradients within the TVD
region.

**Superbee**

$$
\phi(a,b) =
\max\!\Bigl(0,\;\min(2r,\,1),\;\min(r,\,2)\Bigr)\;a
$$

Sharpest TVD limiter; excellent at capturing steep gradients but may
introduce mild staircase artefacts.

---

### TVD Diagram

All non-trivial limiters fall within Sweby's TVD region in the $(r,\,\phi)$ diagram:

| Limiter | Dissipation | Smoothness |
|---------|:-----------:|:----------:|
| MinMod | High | $C^0$ |
| Van Leer | Medium | $C^\infty$ |
| Van Albada | Medium-low | $C^\infty$ |
| MC | Low | $C^0$ |
| Superbee | Lowest | $C^0$ |

<figure>
  {% include "theory/images/sweby.svg" %}
</figure>

---

## Shock Detection

Near strong shocks the MUSCL reconstruction can produce oscillations even with TVD limiters. This is the case of the well-knwon shock carbuncle. MOSE uses a Jameson-type pressure sensor to detect shocks and smoothly degrade to first order in their vicinity.

The key `shock-detector` of `input.ini` turns on the sensor below (`Tramel`) or that of Chen (`Chen`); without it, HLLC+Tramel and HLLE++ take the sensor below, HLLC+Chen and AUSM+M that of Chen, and the other Riemann solvers none ($\beta = 1$).

### Sensor formulation

A $3\times 3\times 3$ stencil of pressures centred on cell $(i,j,k)$ is
used to compute second-order undivided differences in each coordinate
direction:

$$
\kappa_d =
\frac{|p_{+} - 2p_0 + p_{-}|}{p_{+} + 2p_0 + p_{-}},
\qquad d = 1, 2, 3
$$

The shock strength is the maximum over all directions:

$$
\sigma = \max(\kappa_1,\;\kappa_2,\;\kappa_3)
$$

### Blending function

A smooth blending is applied using a tuneable threshold $\Delta = 20$:

$$
\beta =
\begin{cases}
1 - \tanh\!\bigl(10\,(\sigma/\Delta)^3\bigr)
  & \text{if } \sigma < 1/\Delta \\[4pt]
0 & \text{otherwise}
\end{cases}
$$

- $\beta \to 1$ in smooth flow — full second-order reconstruction
- $\beta \to 0$ at shocks — first-order (no reconstruction)

This $\beta$ is passed directly to the MUSCL reconstruction step, where
it multiplies the limited slope.

!!! note "Solver-level integration"
      Some Riemann solvers also use the shock-detection flag internally to tune their numerical dissipation, providing an additional layer of robustness.

---

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
$\Delta l$ of the MUSCL reconstruction. The ghost cells behind a boundary face
take their metric from nodes extrapolated linearly along the index lines,
$\mathbf{N}_{-g} = 2\mathbf{N}_{-g+1} - \mathbf{N}_{-g+2}$ (faces between
blocks: the cell across; symmetry: the mirror image; the wedge: the rotated
cell): both layers have the edge across the face of the first cell, so their
size across the face never vanishes nor turns over however stretched the mesh,
and the first is centred at $2\mathbf{x}_f - \mathbf{x}_1$ ($\mathbf{x}_f$ the
mean of the nodes of the face, $\mathbf{x}_1$ that of the cell). Their edges
along the face are extrapolated as well, and turn a layer over where they change
sharply across the first cell, as where the index lines across the face spread
apart away from it: on a plane 2-D mesh where they grow 3 times or more the
first layer, from 5/3 times the second. Behind a chimera face the first layer
enters the viscous flux (*Chimera faces*, below). Behind the other faces no flux
reads the metric of either layer: walls take the
metric of the cell, or the stencil of positions, the other faces carry no
viscous flux that needs it, and the wedge rotates the cell. The solver reads
only the length across the face, in the MUSCL reconstruction, and that is the
length of the first cell whatever the edges along the face: their sign is not
checked. The cells and the first layer of ghost cells behind chimera faces are
checked at setup (*Checks of the cells at setup*, above). `diffusive-metric` in
`[MOSE-Numerics]` selects the face metric:

| `diffusive-metric` | Face metric | |
|---|---|---|
| `stencil` (default) | $M_f = \left(J_f^{T}\right)^{-1}$, $J_f^{T} = \left[D_a\mathbf{x},\ D_b\mathbf{x},\ D_c\mathbf{x}\right]$ | the same differences of the cell centres as of the field: exact for a linear field on any mesh (at the symmetry faces and the faces of type `200` for the linear fields that they allow), except at the chimera faces (`102`) |
| `inverse-mean` | $M_f = \left[\tfrac12\left(A_L + A_R\right)\right]^{-T}$ | the face-normal edge is the centre-to-centre vector: exact for a linear field on any stretching normal to the face |
| `mean` | $M_f = \tfrac12\left(M_L + M_R\right)$ | faster; biased on stretched meshes |

With `stencil`, the columns of $J_f^{T}$ are the differences of the cell
centres (the mean of the eight nodes) that $D_a\phi$, $D_b\phi$ and $D_c\phi$
take of the field. For a linear field $\phi = \phi_0 + \mathbf{g}\cdot\mathbf{x}$
every difference is $D\phi = J_f\,\mathbf{g}$, so the face gradient is
$\mathbf{g}$ on skewed, stretched and non-smooth meshes and across folds at the
faces between blocks: the difference approximations of the field and of the
coordinates must be the same [8]. Vinokur states the rule for the gradient of a
potential (p. 26) and writes the face gradient of the viscous terms on the same
values as here, Eq. (90)–(91), with metrics from averaged areas and volumes,
Eq. (93)–(94), as in the cross-derivative terms of CFL3D [12]; `stencil`
applies the rule to that stencil. Its component along the centre-to-centre
vector is the difference of the two cell values over their distance, as in the
face-tangent and node-averaging schemes [9, 10], with the derivatives along the
stencil's tangential directions taken between means of four cells, at the mean
of their centres (near, not in, the plane of the face); on a 2-D mesh of
parallelograms it is the diamond scheme with least-squares weights [11], which
are exactly 1/4 there. Where the tangential differences of the stencil lie in
the plane of the face, its split of the face area vector is the over-relaxed
one of [13]. Exactness for linear fields does not certify accuracy on curved,
highly stretched cells [14]; that case is not measured here. One effect of the
off-face position of the tangential derivatives is measured (test `MMSResiduals`):
on a wedge whose cells grow along $r$ by a ratio $q$, the axial viscous force of
Homann's flow converges at first order on the first rows next to the axis, with
`stencil` as with `inverse-mean`, as the derivative of $u_r$ along $x$ at the
faces $j$ is that of the mean of the two rows of cells, off the face by about
$\Delta r\,(q-1)/4$. `inverse-mean`
takes the tangential columns
from the mean edges of the two cells instead. On smooth meshes they approach
the centre differences as the cells shrink, but where the mesh is not smooth
(cells of alternating size, nodes moved at random, folds at the faces between
blocks) the face gradient of a linear field stays wrong by several per cent or
more at any resolution. At a fold between two blocks of uniform cells the
tangential columns of the two metrics coincide; there the error of
`inverse-mean` came from the ghost cell of the connection, whose geometry is
extrapolated straight across the fold. The column across the face is kept as $\tfrac12(a_L + a_R)$, which
is $\mathbf{x}_R - \mathbf{x}_L$, and so are the columns of the index
directions with one cell (1-D runs, the $k$ direction of 2-D meshes and
wedges). On 1-D runs and plane 2-D meshes the field does not vary along them,
and 1-D runs give the bits of `inverse-mean`. On a wedge the $k$ neighbours are
the cell turned by $\mp\delta$, so the velocity differs across them, and that
column, the mean edge of the cells, is the difference of the turned positions
divided by $\cos^2(\delta/2)$: the error of order $\delta^2$ of the azimuthal
derivatives (see *Axisymmetric wedge* below). The ghost cells in the stencil take the
position at which their value is right for a linear field:

- faces between blocks (`101`, `201`): the centre of the cell across, plus the
  translation of a periodic pair;
- walls (`301`, `302`, `503`–`505`, and `103`, `104`, coupled to another
  solver) and `400` faces: the position of the class
  of each quantity, $2\mathbf{x}_f - \mathbf{x}_1$ for the quantities with a
  wall value and $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ for the
  extrapolated ones (see *Ghost Cells of Walls and 400 Faces* below); a face
  with such a ghost cell in its stencil takes one metric per group of
  quantities (species and scalars, temperature, velocity, turbulence
  variables);
- faces whose ghost value is $3u_1 - 3u_2 + u_3$ (`408` and the other types not
  listed here): $3\mathbf{x}_1 - 3\mathbf{x}_2 + \mathbf{x}_3$
  ($2\mathbf{x}_1 - \mathbf{x}_2$, or $\mathbf{x}_1$, in a block with two
  cells, or one, across the face; there the value reads the ghost cell of the
  opposite face, so value and position agree only with three cells or more);
  chimera (`102`) faces take the same position, an approximation, since their
  values come from donor cells; the tangential differences of their faces take
  as neighbours of a ghost cell along the face the ghost cells of the
  neighbouring faces, and beyond the edges of the block face the row of its
  ghost cells extrapolated linearly (in the k direction of a 2-D mesh, the
  ghost cell turned by the wedge angle); the MUSCL reconstruction of the first
  interior face takes the distance $l_0 = \Delta l_1$ of the ghost cell instead,
  so on cells stretched across these faces the slope behind the face is scaled
  by a factor of order one, which the limiter bounds;
- mapped faces (`410`): the file gives the state of the cell across the face,
  taken as the ghost cell, the first cell mirrored through the face; the
  quantities of the file sit at its centre, $2\mathbf{x}_f - \mathbf{x}_1$, and
  the turbulence variables, which the file does not give (the values of the
  first cell per unit mass), at $\mathbf{x}_1$; a face with such a ghost cell in
  its stencil takes one metric per group of quantities, as next to walls. The
  soot, copied in the same way, takes the position of the species: no flux takes
  its gradient. Where the cells of the solver that wrote the file have another
  width across the face, its values sit at half that width from the face, and
  the derivatives across the first row carry the error of that difference;
- zero-gradient faces (`0`, `401`–`407`, `420`): the cell itself;
- symmetry (`300`): the mirror image of the cell across the face; on a face of
  zero area (an axis at $r = 0$) the cell itself, since the ghost value there is
  a copy of the cell, and the viscous flux of that face is zero;
- the $k$ faces of a wedge or of a 3-D sector (`200`): the source of the ghost
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
periodic face. For a field that breaks the condition the value of the ghost cell
(reflected, copied or rotated) is not the value of the field at any point, so no
position makes the stencil exact: the error is proportional to the part that the
condition forbids ($\mathbf{g}\cdot\mathbf{n}$, $G - RGR$, $\mathbf{g}\cdot\mathbf{T}$),
it stays on the faces next to the boundary, and it vanishes as the solution
satisfies the condition. Behind the other faces between blocks and the
zero-gradient faces the gradient is exact for every linear field; behind walls
and `400` faces for the linear fields of the next section; behind mapped faces
for the linear fields whose values in the file are those at the centre of the
ghost cell. Behind the faces whose ghost value is $3u_1 - 3u_2 + u_3$ it is
exact for every linear velocity, but for the mass fractions, the temperature and
the specific turbulence and scalar variables only where the quantities that the
rule extrapolates, $\rho_s$, $p$ and $(\rho\phi)$, are linear too (for $Y_s$ and
$\phi$ where the density is uniform, for $T$ where the composition is uniform
as well): their quotients $Y_s = \rho_s/\rho$, $T = p/(\rho R)$ and
$\phi = (\rho\phi)/\rho$ are not linear where the density varies. There the
derivative across the first row has a relative error of second order in the
cell size on uniform cells and of first order on cells that grow by a fixed
ratio across the face.

$M_f$ is computed once at setup and stored on every face of the blocks of each
rank, 27 values per cell (the mesh does not move); the cell centres come from
the nodes, which every rank holds, so no message is needed. Where a face or a
cell is degenerate, the differences of its stencil can stop spanning the space
and the metric is not defined there: the run stops at setup with a message that
names the block, the cell and the face, a declared limit of `stencil`
(`inverse-mean` does not need this metric, but a viscous run with it stops on a
3-D mesh at a boundary face of zero area, below). A wall on faces of zero area,
as on
the axis of a 3-D mesh around it (a polar line), stops too: its gradient is not
defined there, and the axis takes the symmetry condition (`300`). With `300` on
the axis the polar line runs, and the face gradients next to it are exact for
linear fields. With `inverse-mean` and `mean` a viscous run on a 3-D mesh with a
boundary face of zero area, as on a polar line with the axis exactly at $r = 0$,
stops at setup and names `stencil`: next to the axis
the viscous force of the first row of cells grows as the mesh is refined (a
linear field, whose force is zero).

With the axis nodes exactly at $r = 0$ the polar line runs with `stencil`, but
the solution next to the axis is of first order. On smooth fields the residual
of the two rows of cells next to the axis does not converge as the mesh is
refined: its viscous part grows as $1/h$ and its convective part does not
decrease. With `inverse-mean` and `mean` the viscous residual of the cells next
to a symmetry plane (`300`) that bounds a sector of a polar mesh does not
converge either: on smooth fields that the plane allows it stays of order one as
the mesh is refined, while with `stencil` those cells converge as the others.

In one dimension $M = 1/h$, so `mean` divides the difference across the face by
the mean of $1/h$ instead of the centre-to-centre distance. With a cell-size
ratio $r = h_R/h_L$ it overstates the face-normal gradient, and so the viscous,
conductive and diffusive fluxes, by

$$
\frac{(1+r)^2}{4r}:\qquad
+0.2\,\% \;(r=1.1),\quad +0.8\,\% \;(r=1.2),\quad +4.2\,\% \;(r=1.5),\quad +8.9\,\% \;(r=1.8)
$$

This bias does not decrease with grid refinement at a fixed growth ratio, and it
is present in the converged steady solution. `mean` is adequate where adjacent
cells grow by about 1.2 or less; boundary-layer meshes often grow faster near
the wall. With `inverse-mean` MOSE stores $A^{T}$ per cell (9 values), so each
face needs one $3\times 3$ inversion.

With `stencil`, a face between blocks takes one face gradient with its stencil
metric, as an interior face, and a symmetry face the gradient between the cell
and its mirror image. With `inverse-mean` and `mean`, block-connection faces
use the inverse mean, as chimera faces do with every option: across a large
cell-size jump at a block interface, the mean makes the normal gradient so steep
that the explicit diffusion becomes unstable.

At a wall face (301, 302, and 503–505 when they do not ablate) the gradient
takes the difference between the cell and the wall value across the face.
With `stencil` it also takes the differences of the wall values along the
wall, between the neighbouring wall faces of the same block face and type
(central, or one-sided where the wall ends or goes on in another block face,
whose faces it does not take: exact there for a linear wall value, first order
for a curved one), with the metric of the same differences of the positions:
the vector from the face centre (the mean of its four nodes) to the cell centre
across the face, the differences of the face centres along it. For a linear
field whose wall values are its values at the face centres the wall gradient
is exact: the temperature of a 302 wall that
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
neighbouring face centres. A boundary face that keeps only the normal distance,
as in [13] (Eq. 3.57–3.59) and in OpenFOAM (OpenFOAM-dev, commit dfa6c7e,
`fvPatch.C`), has the error of `inverse-mean` and `mean` given below.
With `inverse-mean` and `mean`, and for the quantities whose wall value is not
known along the wall (the temperature of 301 walls and of the surface models
of 503–505), the derivatives along the wall are zero and the metric is that of
the wall cell: on a cell skewed by θ the normal gradient is then wrong by
$(\partial\phi_w/\partial s)\tan\theta$, which refinement does not reduce.
The turbulence variables take no difference of their wall values along the
wall, with every option: their wall values are uniform along a smooth wall
(zero for $\tilde\nu$ of Spalart–Allmaras, for $k$ and for the Reynolds
stresses), but not those of $\omega$ of the $k$–$\omega$ models and of
SSG/LRR, $C\mu_l/d^2$ with $d$ the distance of the cell centre from
the wall, which follow the height of the first cell, nor that of $\tilde\nu$
on a rough wall, taken from the value of the cell. Their derivative along the
wall is left out, with the error above on skewed cells.

A 301 or 302 wall of a stationary face of a rotating frame moves in the frame
with $\mathbf{u}_w = -\boldsymbol{\Omega}\times(\mathbf{x}-\mathbf{x}_o)$
($\mathbf{x}_o$ the origin of the frame), a velocity linear in $\mathbf{x}$
with the constant gradient $-[\boldsymbol{\Omega}\times]$. With every option
the wall gradient takes the differences of $\mathbf{u}_w$ along the position
differences of its metric, the columns of $M^{-1}$: those of the wall stencil
above with `stencil`, the edges of the wall cell with `inverse-mean` and
`mean`. The part of $\mathbf{u}_w$ in the wall gradient is then exact with
`stencil` on cells of any shape, and with `inverse-mean` and `mean` on
parallelepiped cells, where the metric of the wall cell is exact; the rest is
that of a wall at rest. Without these differences the wall stress would be
wrong by the order of $\mu\,|\boldsymbol{\Omega}|$ with every metric and on
any mesh, orthogonal cells included, because the stress takes the derivatives
of the velocity along the wall: a gas at rest in the absolute frame, a rigid
rotation in the rotating frame, would feel a shear on the wall. The part of
$\mathbf{u}_w$ is exact because the values sit where the metric places them:
the value of a cell at the centre of the cell (the mean of its eight nodes),
the wall velocity at the centre of the wall face (the mean of its four nodes).
A change of either position has to be made in the metric too.

**Chimera faces (`102`).** The first layer of ghost cells behind a chimera face is the
ghost cell of the face, whose metric enters its viscous flux: a first layer that
turns over (index lines that spread apart sharply across the first cell, above)
or whose determinant is not that of the first cell in sign and to $10^{-12}$ in
size stops the run at setup, with its block, cell and face; the remedy is the
mesh, with index lines that change smoothly across the first two node layers at
that face. The ghost values of a chimera face (`102`) come from the donor cells and weights that the tool which
wrote `bc.txt` computed for its own ghost cells; the face metric of their
viscous flux takes these linear ones, so where the tool builds its ghost cells
otherwise (on cells stretched across the face) the normal gradient is taken over
a distance that differs by as much (no test with `102` faces).
The viscous flux of a chimera face is the mean of two fluxes: one with the
central differences along the face of the row of the cell, one with those of
the row of the ghost cells, filled from the donors. With a stress linear in the
gradient (laminar, Spalart–Allmaras) the mean is the flux of the mean gradient,
whose differences along the face are the mean of those of the two rows, as on
an interior face (Vinokur 1986, NASA CR-177416, eqs. 90–91, p. 23). Beyond the
edges of the block face the row of the ghost cells is extrapolated linearly, in
the variables whose differences the stencil takes: the difference there is the
one-sided one, exact for a linear field (a copy of the ghost cell would halve
it, and the face would take 3/4 of the derivative of a linear field along the
edge). In the k direction of a 2-D mesh the neighbours of a ghost cell are the
ghost cell turned by the wedge angle, as those of a cell are; across a block
one cell wide in another direction they are the ghost cell itself.

!!! note "Tangential derivatives"
    With `inverse-mean` and `mean` the tangential derivatives on a face are
    central differences over two cells, divided by the cell's own tangential
    edge. With stretching in a tangential direction they carry the same
    $(1+r)^2/(4r)$ factor, which neither option corrects. In a thin boundary
    layer this enters only through the streamwise variation of the shear on
    faces normal to the flow, so its effect is much smaller than that of the
    face-normal gradient. `stencil` divides them by the difference of the
    centres they span, and has no such factor.

---

## Ghost Cells of Walls and 400 Faces

Behind a viscous wall (`301`, `302`, `503`–`505`, and the walls `103` and `104`
coupled to another solver) and an extrapolation face (`400`) the two layers of ghost cells are filled once per stage, before the
fluxes, and the viscous fluxes, the MUSCL reconstruction and the source terms
of the turbulence models read the same values. Extrapolated are the quantities
whose differences the viscous fluxes and the sources take, the mass fractions
$Y_s$, the temperature $T$, the velocity $\mathbf{u}$ and the specific
turbulence and scalar variables $\phi = (\rho\phi)/\rho$, and the pressure
$p$; the ghost state follows from them, $\rho = p/(R(Y)\,T)$,
$\rho_s = Y_s\rho$, $(\rho\phi) = \rho\phi$. With $q_1$, $q_2$ the values in the
wall cell and in the next one across the face:

| Class | Quantities | First ghost cell | Position of its value |
|---|---|---|---|
| wall value | velocity at `103`, `104`, `301` and `302` (zero, or the wall velocity of a stationary face of a rotating frame); the turbulence variables at the walls ($k$, $R_{ij}$, $\tilde\nu$ on a smooth wall: zero; $\omega$, $\tilde\nu$ on a rough wall: the wall value of the model) | $2q_w - q_1$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| wall temperature | $T$ at `302` | $\max(2T_w - T_1,\ T_w/2)$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| extrapolated | $Y_s$, $p$; $T$ at `103`, `104`, `301`, `400`, `503`–`505`; the velocity at `400` and `503`–`505`; turbulence and scalars at `400`, scalars at the walls | $q_1 + \psi(q_1 - q_2)$ | $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ |

$\mathbf{x}_f$ is the centre of the face (the mean of its four nodes) and
$s = l_0/l_1$ the ratio of the distances of the MUSCL reconstruction at the
first interior face, $l_0 = (\Delta l_g + \Delta l_1)/2 = \Delta l_1$ (the ghost
cell has the edge across the face of the first cell, see *Gradient Computation
for Diffusive Fluxes*), $l_1 = (\Delta l_1 + \Delta l_2)/2$:
$s = 2\Delta l_1/(\Delta l_1 + \Delta l_2)$, 1 on uniform cells. $\psi \le s$ is
the largest value that keeps each species, the pressure, the temperature and
each scalar within $[q_1/2,\ 2q_1]$; the species are then normalised to sum
one; the velocity is not limited ($\psi = s$). The wall values of the turbulence
variables are those of the wall fluxes, from the same wall state (`302`: $T_w$;
`301`: the wall temperature estimated from the heat flux;
`503`–`505`, whose surface model sets the wall state at every call, and the
coupled walls `103` and `104`, whose viscous flux the other solver gives: the
state of the wall cell, as an adiabatic wall). With one cell across the face the
extrapolated quantities are copied (position $\mathbf{x}_1$). In an inviscid
run the walls carry no wall values and every quantity is extrapolated, as
behind a `400` face. The second layer takes $2q_w - q_2$,
$\max(2T_w - T_2,\ T_w/2)$ and $q_1 + 2\psi_2(q_1 - q_2)$; no flux reads it.

Properties, for every mesh and in 3-D, with the hypotheses written:

- **Positivity.** The ghost state has $\rho_s \ge 0$, $p > 0$, $T > 0$ and
  so $\rho > 0$; the extrapolation amplifies $q_1$ at most by $1 + s$. The
  temperature and the mass fractions are what the viscous fluxes take
  differences of, so they must be those of a gas. The extrapolation of every
  primitive as $3q_1 - 3q_2 + q_3$, which the other types keep (below), can be
  negative (density that rises steeply away from a hot wall) and amplifies by 3
  and 6 on the two layers.
- **Wall value.** For the first two classes without the floor, the value
  interpolated at the face centre, $(G_1 + q_1)/2$, is the wall value, for any
  profile (zero velocity, zero $k$, the model value of $\omega$). The floor
  $T_w/2$ acts only where $T_1 > 1.5\,T_w$, a first cell much hotter than a
  cooled wall that the mesh does not resolve; there the wall flux still takes
  $T_w$, and the ghost cell stays on the cold side of the wall.
- **Linear fields.** The rules are affine combinations of the values whose
  coefficients sum to one, so for a linear field the ghost value is the field at
  the position of the table: for the extrapolated class for every linear field,
  for the other two for the linear fields that take the wall value at the face
  centres (zero velocity on the wall, $T(\mathbf{x}_f) = T_w$ of each face).
  With these positions in the face metrics (`stencil`) the face gradients next
  to walls and `400` faces are exact for those fields, whatever the density:
  the quantities differenced are extrapolated themselves. So is the MUSCL
  reconstruction at the first interior face along a straight grid line: both
  positions are $l_0 = \Delta l_1$ from $\mathbf{x}_1$, the distance it assumes.
  The density gradient of the source terms takes the position of the
  temperature. At `301`, `400` and `503`–`505` the pressure, the composition
  and the temperature all sit at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$,
  so the density of the ghost cell is that of the field there, exact for linear
  $p$, $Y_s$ and $T$ where no limit acts, with no further hypothesis. At `302`
  the temperature sits at $2\mathbf{x}_f - \mathbf{x}_1$, $\Delta l_1$ from
  $\mathbf{x}_1$ along the edge across the face of the wall cell, and the
  pressure and the composition at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$,
  along the line from $\mathbf{x}_2$ to $\mathbf{x}_1$, $\Delta l_1$ from
  $\mathbf{x}_1$ on a straight grid line: there the two positions are the same
  point and the density is exact; elsewhere they differ by about $\Delta l_1$
  times the angle $\theta$ between the two directions, and the density is exact
  where the pressure and the composition do not vary across the wall cell (a
  boundary layer), its relative error otherwise of the order of
  $\Delta l_1\,\theta\,|\nabla p|/p$, and the same for the composition.
- **Fluxes and sources.** One set of ghost cells: the wall flux routines write
  no ghost values, and the source terms take the positions of the table.

Where a limit acts the ghost value is no longer the value of the field at its
position, and the gradients next to the face are first order there: it happens
only where a quantity changes by more than a factor two between the first two
cells, or where the first cell is 1.5 times hotter than a `302` wall. At the
end of a run one line reports the ghost fills in which either acted:

```
 Ghost cells of walls and 400 faces, fills with the extrapolation cut: <n>, with the 302 temperature at T_w/2: <m>
```

Declared limits:

- `503`–`505`: the velocity and the temperature are extrapolated (the surface
  may blow, its state is set at every call of its routine) and the wall values
  of the turbulence variables are those of an adiabatic wall at the state of
  the cell: exact for $k$, $R_{ij}$ and $\tilde\nu$ (zero), while that of
  $\omega$ differs, at uniform pressure, from the value at the state of the
  surface by the factor $\mu(T_1)T_1/(\mu(T_s)T_s)$ (about $(T_1/T_s)^{1.7}$
  with $\mu \propto T^{0.7}$): of order one but bounded, it enters only the
  source terms (the gradient of $\omega$ in the cross-diffusion of SST,
  $\nabla k\cdot\nabla\omega$) and the limited MUSCL reconstruction of the
  first face, while the wall flux keeps the state of the surface; with the
  blowing correction of $\omega$ (`-blowcorr`) the wall flux takes the corrected
  value, the ghost cell the uncorrected one;
- the temperature of a `301` wall is extrapolated (consistent, exact for every
  linear field) rather than set from the heat flux;
- the SRM grain (`502`) and the other types keep the extrapolation
  $3q_1 - 3q_2 + q_3$ of the primitives; in a RANS run the turbulence variables
  of a `502` grain take their wall values instead, $2\phi_w - \phi_l$ at
  $2\mathbf{x}_f - \mathbf{x}_1$ as behind a wall, with the wall values of the
  flux of the grain (the blowing wall while it burns, the wall at the flame
  temperature when it does not).

## Rotational Periodicity and Axisymmetric Wedge (Type 200)

Type 200 on the $k$ faces closes a mesh that turns about the $x$ axis: an
axisymmetric wedge (one cell in $k$) or a 3-D sector (several). The flow
beyond face 5 is the flow next to face 6 turned back by the angle $\delta$
of the wedge or sector, and the flow beyond face 6 the flow next to face 5
turned forward. With $g = +1$ on face 5 and $g = -1$ on face 6, the ghost
cell takes the state of its source cell (on a wedge the cell itself, on a
sector the cell next to the opposite $k$ face of the same block, in the same
column) with the velocity turned about $x$ by $\varphi = -g\,\delta$:

$$
v_g = v\cos\varphi - w\sin\varphi = v\cos\delta + g\,w\sin\delta, \qquad
w_g = v\sin\varphi + w\cos\varphi = -g\,v\sin\delta + w\cos\delta .
$$

$u$ and the scalars do not change under a turn about $x$, and are copied.
The ghost cells and the fluxes of these faces take $\cos\delta$ and
$\sin\delta$ computed once, with $\delta$, when the mesh is read. The record
of a face of type 200 names no other block: a block closed by type 200 holds
the whole wedge or sector in $k$, a column of cells with type 200 on one of
its $k$ faces has it on the other, and a record of type 200 that does not fit
this stops the run at setup (see the note on type 200 in
[Boundary Conditions](../user/boundary-conditions.md#euler-boundaries)).

MOSE finds the wedge and the sector from the nodes (the rule is in the note
*Mesh type* of [Using MOSE](../user/using.md#case-directory-structure)): the
last $k$ node plane of the block must be the first one turned or translated to
within 1 % of the move of the nodes per cell in $k$, measured along the $k$
lines: the largest, over the nodes of the first $k$ node plane, of the sum of the
moves of the node across the cells in $k$ (each its largest component), over the
number of those cells (the thickness of a 2-D block). Along straight $k$ lines it
is the move from the first plane to the last per cell, and on a sector the
length of the arc per cell in the norm of the largest component (between 0.71
and 1 of the arc), which the chord from plane to plane understates past half a
turn (4.2 times at 270 degrees).

### Axisymmetric wedge

A 2-D axisymmetric case is a wedge of the 3-D mesh: one cell in $k$, whose two
$k$ faces lie in the meridional planes $\theta = \mp\delta/2$ and carry the
boundary type `200`. MOSE integrates the 3-D equations on the wedge; there are
no axisymmetric source terms. The ghost cell across each wedge face is the cell
itself rotated by $\mp\delta$ about the $x$ axis, so every face flux follows
from the state of the cell:

- inviscid: the pressure and the momentum carried by the cell state rotated to
  the face (no Riemann problem). In the radial momentum the two faces give the
  $+p/r$ and the centrifugal terms of the axisymmetric equations;
- viscous: the flux of a face between the cell and its rotated image, with the
  stencil of an interior face (the $i$ and $j$ neighbours of the image rotated
  in the same way) and, with `diffusive-metric = stencil`, the face metric of
  the same differences of the positions, the row of the image the row of the
  cell turned (exact for the linear fields that the rotation allows, see the 3-D
  sector below); with `inverse-mean` and `mean`, the inverse-mean face metric of
  the cell and its image. Both are computed once at setup (they depend on the
  mesh only).

On a planar 2-D mesh (faces 5 and 6 parallel, no wedge angle) type `200` makes
the image the cell itself, and the fluxes of the two faces cancel.

For an axisymmetric field the mass, energy and scalar fluxes of the two faces
are equal and opposite, and only the momentum fluxes are applied. Without swirl
the viscous traction of each face is the normal stress

$$
\tau_{\theta\theta} = 2\mu\,\frac{u_r}{r} + \lambda\,\nabla\!\cdot\mathbf{u},
\qquad \lambda = -\tfrac23\,\mu ,
$$

with $\mu$ including the eddy viscosity, and the two faces together give the
$-\tau_{\theta\theta}/r$ of the radial momentum equation [5, 6]; with swirl
they also carry $\tau_{r\theta}$ and $\tau_{x\theta}$. Near the axis
$\tau_{\theta\theta}/r$ is as large as the radial stress terms it balances: in
Homann's stagnation-point flow [7] $\tau_{rr} = \tau_{\theta\theta}$, and
$\tau_{rr}/r$ and $\tau_{\theta\theta}/r$ cancel exactly while each grows like
$1/r$ towards the axis.

The azimuthal derivatives of the $i$ and $j$ faces are differences across the
rotated ghost cells divided by the $k$ edge of the cell, which is a chord: they
carry a relative error of order $\delta^2$ that the mesh does not reduce
($\delta^2/12$ of the hoop stress $2\mu u_r/r^2$ in the radial balance,
$2.5\times10^{-5}$ for $\delta = 1^\circ$). Keep the wedge thin.

Put the axis nodes at a small radius, as MOSE's meshes do ($10^{-6}\,\Delta r$):
a viscous run on a wedge whose axis is exactly at $r = 0$ stops at setup, with
every `diffusive-metric`. On such an axis the faces have zero area and no
normal, so their ghost cells copy the cells. With `inverse-mean` and `mean` the
wedge faces and the faces $i$ of the first row of cells take the derivative in
$r$ across the axis with the ghost cell one cell away and would see half of
$\partial u_r/\partial r$: an error of the hoop stress of those cells that
refinement does not reduce ($1/6$ of $2\mu u_r/r^2$ in Homann's flow, on 32 and
64 cells). With `stencil` the copy sits at the cell and the faces take the
difference on one side of it: in Homann's flow, linear in $r$, the first row
gets the hoop stress of a small axis radius, so the error of $1/6$ is gone; a
field that is not linear in $r$ stays farther from the result of the small
radius (test `MMSResiduals`). The stop stays with every metric, a declared
limit: the axis at $r = 0$ is not treated completely (ghost cells beyond the
axis, not copies).

Next to the axis the hoop stress also stiffens the explicit step. In the radial
momentum of the first row of cells, at $r \approx \Delta r/2$, it is a damping
$(2\mu+\lambda)\,u_r/r^2 \approx \tfrac{16}{3}\,\mu\,u_r/\Delta r^2$ ($\mu$ with
the eddy viscosity), with couplings to the neighbours through $\lambda\,\nabla\cdot\mathbf{u}$.
The local time step includes them ([Time Integration](time-integration.md#local-time-step)),
so the `vnn` limits of the rule hold next to the axis too; without them, at a
Prandtl number of 10, the first rows broke down at the `vnn` the rule declares
stable (test `WedgeAxisStep`).

### 3-D sector

On a 3-D sector (type `200` with more than one cell in $k$) faces 5 and 6
make one interior face of the whole annulus, computed once, by the record of
face 6, for the inviscid and for the viscous flux. The image beyond face 6 is
the cell next to face 5 turned forward by $R$, the turn about $x$ by the sector
angle: the $i$ and $j$ differences of the image row are those of that row turned, the
metric of the image is that of its cell turned, $M R^{\mathsf T}$ (each row of
$M = \partial\xi/\partial\mathbf{x}$ is a vector), and the wall distance and
roughness of the face are the mean of the two cells, as on an interior face.
The cell next to face 6 takes the flux and the cell next to face 5 the same
flux with its momentum turned back, $R^{\mathsf T}F$: mass, energy, species,
turbulence variables, axial momentum and angular momentum about $x$ are
conserved. A translated sector (faces 5 and 6 parallel) is the same with
$R = I$. The viscous face of the wedge and of the sector takes, with
`diffusive-metric = stencil`, the metric of the same differences of the
positions as the gradient takes of the values: across the face the image, at
the position of its ghost cell, minus the cell; along the face the differences
of the row of the cell and of the row of the image, which are those of the row
of the source turned by $R$, positions and values (a ghost cell of the row of
the source, behind a wall too, turned with it). A linear field that the
periodicity allows, $\mathbf{u}(R\mathbf{x}) = R\,\mathbf{u}(\mathbf{x})$, is so
exact at these faces on any mesh; where a wall ghost cell enters the row of the
cell or of the source, the face takes one metric per group of quantities, as an
interior face. With `inverse-mean` and `mean` the face takes the inverse of the
mean of the inverse metrics of the cell and of its image, whose columns along
the face are the mean edges of the two cells, not the differences of the
centres: where the cell size grows along the face by a ratio $q$ the
derivatives along the face carry the factor $(1+q)^2/(4q)$ of the note
*Tangential derivatives* above, a relative error of $(q-1)^2/(4q)$
($8.3\times10^{-3}$ at $q = 1.2$), and in the first row next to an axis with
symmetry (`300`), whose ghost cell is the mirror image of the cell, $(3+q)/4$, a
relative error of $(q-1)/4$ ($5\times10^{-2}$ at $q = 1.2$); refinement at a fixed
ratio does not reduce them. For $u_r = U r$, whose viscous force is zero, the
radial viscous force of a cell of the wedge is, with `inverse-mean`, that
relative error times the hoop term $\tau_{\theta\theta}/r$ on top of the error
of order $\delta^2$ above, and with `stencil` that error alone,
$\sin^2(\delta/2)$ of the hoop term in every cell (test `MMSResiduals`).
A viscous run on a sector whose axis nodes are
exactly at $r = 0$ stops at setup, as on the wedge, with every
`diffusive-metric`: behind faces of zero area the ghost cells copy the cells,
and with `inverse-mean` and `mean` the periodic faces of the first row of cells
would see half of $\partial u_r/\partial r$ ($1/6$ of the hoop term in Homann's
flow; with `stencil` the error of a small axis radius, next to the periodic faces
too). Put the axis nodes at a small radius.

## Faces Between Blocks

A face shared by two blocks (block connection `101`, periodic connection `201`)
has one flux, computed once and given to both cells with opposite signs, as an
interior face. Of the two boundary records of the face, the one of the cell on
the low-index side (whose face is 2, 4 or 6 while the other is 1, 3 or 5)
computes the inviscid and the viscous flux with the shock sensor of its cell,
as the interior face between cells $n$ and $n+1$ uses that of cell $n$, and with
the preconditioning reference velocity $U_r$ of both cells (the ghost cell
carries that of the other block). When the two faces have the same parity, as
where blocks meet with their index directions permuted, the record with the
smaller (block, face, $k$, $j$, $i$) computes it, so that the result does not
depend on the order of the records in `bc.txt`. The face therefore conserves
mass, momentum and energy exactly whatever the Riemann solver, the shock sensor
and the preconditioning, on any number of MPI ranks; between ranks the other
cell's share travels in one more message per Runge–Kutta stage.

The inviscid flux of a face between blocks is, bit for bit, that of the same
face in the unsplit mesh when three conditions hold: the index directions of
the two blocks agree across the face; the two blocks store the nodes of the
common face with the same bits (the ghost cell takes its geometry from the
nodes of the other block, so nodes that differ in the last digits give a
slightly different face); and `preconditioning-Ur-smooth = 0`. This is an
identity of the fluxes, not of the solution: in 2-D and 3-D a cell next to the
cut adds the flux of the face between blocks after those of its interior faces,
so its sum is taken in another order than in the unsplit mesh and the two
solutions differ at round-off. Only in 1-D, where a cell sums two face fluxes,
is the solution of an inviscid run on the split mesh the same bit for bit (the
viscous flux of a face between blocks has its own stencil, below).

With `preconditioning-Ur-smooth` greater than 0 the max-dilation of $U_r$ in a
cell next to a face between blocks also sees the cell across the face, through
the ghost cell, with its value at the start of the stage. In the unsplit mesh
the smoothing sees the new value of the neighbour, so with the smoothing on a
split mesh does not give the fluxes of the unsplit one bit for bit.

With `diffusive-metric = stencil` the viscous flux of a face between blocks is
that of an interior face: one face gradient, from the cell, the ghost cell and
their tangential neighbours (those of the ghost cell are cells of the other
block, `bc%Pg(:,3:6)`), with the stencil face metric. With `inverse-mean` and
`mean` it uses its own stencil (two gradients, each with the tangential
differences of one of the two cells, and the mean of their fluxes) and the
inverse-mean face metric.

The ghost cells and the neighbours along the face of those of the faces
between blocks (`bc%Pg`) are filled at every Runge–Kutta stage, before the
fluxes, from the state of that stage, in an order such that every value the
fluxes read is that of the stage, with any number of threads and ranks, under
five conditions:

1. the source block of every face between blocks has at least two cells across
   that face;
2. every block with a face whose ghost values are extrapolated
   ($3u_1 - 3u_2 + u_3$) has at least three cells across it;
3. `bc.txt` has one record for every cell face of the boundary;
4. there is no chimera (`102`) nor mapped (`410`) face;
5. the tangential steps $(d_{11}, d_{12})$ and $(d_{21}, d_{22})$ of every
   `101` and `201` record are unit steps along one index.

The first is checked at setup: a connection whose source block has one cell
across its face stops the run. The others are not checked: where one does not
hold, a ghost value can be that of the previous stage, or depend on the order
of the records or of the threads.

---

## References

1. B. van Leer, "Towards the ultimate conservative difference scheme.
   V. A second-order sequel to Godunov's method," *J. Comput. Phys.*,
   32(1), 1979.
2. P. K. Sweby, "High resolution schemes using flux limiters for
   hyperbolic conservation laws," *SIAM J. Numer. Anal.*, 21(5), 1984.
3. A. Jameson, W. Schmidt, E. Turkel, "Numerical solution of the Euler
   equations by finite volume methods using Runge–Kutta time-stepping
   schemes," AIAA-81-1259, 1981.
4. R. Tramel, R. Nichols, and P. Buning. "Addition of improved shock-capturing schemes to OVERFLOW 2.1." 19th AIAA Computational Fluid Dynamics. 2009. 3988.
5. R. Vandenhoeck, A. Lani, "Implicit high-order flux reconstruction solver
   for high-speed compressible flows," arXiv:1811.08666v1, 2018, Appendix A.2,
   eqs. (A.6)–(A.7) (axisymmetric source $(0, 0, p - \tau_{\theta\theta}, 0)^T$).
6. Eilmer 4, Gas Dynamics Toolkit, `src/eilmer/fvcell.d`, lines 1213–1228 (viscous
   source of axisymmetric flow), https://github.com/gdtk-uq/gdtk, commit 826496a;
   the code is described in N. N. Gibbons, K. A. Damm, P. A. Jacobs, R. J. Gollan,
   "Eilmer: an open-source multi-physics hypersonic flow solver," *Comput. Phys.
   Commun.*, 282, 108551, 2023, doi:10.1016/j.cpc.2022.108551.
7. F. Homann, "Der Einfluss grosser Zähigkeit bei der Strömung um den Zylinder
   und um die Kugel," *Z. Angew. Math. Mech.*, 16, 153–164, 1936 (not read; as
   cited by J. Magnaudet and M. Abbas, arXiv:2011.10755v1, 2020). The equation
   of the flow, f''' + 2ff'' − f'² + 1 = 0 with f(0) = f'(0) = 0, f'(∞) = 1,
   u_r = a r f'(η), u_x = −2(aν)^{1/2} f(η), as written by P. Rajamanickam,
   arXiv:2508.17982v4, Eq. (12), and by H. Guo, X. Zhuang, X. Meng,
   T. Rabczuk, arXiv:2010.05682v1, Eq. (7) with β₀ = 2, β = 1, whose Table 8
   gives f''(0) = 1.31194.
8. M. Vinokur, "An analysis of finite-difference and finite-volume formulations
   of conservation laws," NASA CR-177416, 1986: pp. 15–16, Eq. (60a) (the
   area vector of a face of a hexahedron as half the cross product of its
   diagonals; the shape of each face "precisely defined and consistently used
   by the two neighboring cells" so that the volumes sum to the total volume;
   the doubly ruled face); p. 26, Eq. (101) ("the two
   difference approximations must be the same"); p. 23, Eq. (90)–(91) (the
   transport terms on an auxiliary cell, with means of four cells along the
   face); p. 40, Eq. (153) (the difference across a wall). Journal version:
   *J. Comput. Phys.*, 81, 1989.
9. B. Diskin, J. L. Thomas, E. J. Nielsen, H. Nishikawa, J. A. White,
   "Comparison of Node-Centered and Cell-Centered Unstructured Finite-Volume
   Discretizations: Viscous Fluxes," *AIAA J.*, 48(7), 1326–1338, 2010,
   doi:10.2514/1.44940, Sec. III.A.
10. H. Nishikawa, "Beyond Interface Gradient: A General Principle for
    Constructing Diffusion Schemes," AIAA Paper 2010-5093, 2010, App. C.
11. Y. Coudière, J.-P. Vila, P. Villedieu, "Convergence rate of a finite volume
    scheme for a two dimensional convection-diffusion problem," *ESAIM: M2AN*,
    33(3), 493–516, 1999, Eq. (9) and Lemma 6.2.
12. CFL3D, `source/cfl3d/libs/gfluxv1.F` (cross-derivative terms of the viscous
    fluxes), https://github.com/nasa/CFL3D, commit 776d516; S. L. Krist,
    R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual (Version 5.0)*,
    NASA/TM-1998-208444, 1998, App. C.2, p. 249 (the volume at a cell interface
    as the mean of its neighbours, for the thin-layer terms; the manual, p. 27,
    leaves out the cross-derivative terms, which gfluxv1.F adds).
13. H. Jasak, *Error Analysis and Estimation for the Finite Volume Method with
    Applications to Fluid Flows*, PhD thesis, Imperial College, University of
    London, 1996, Eq. (3.32).
14. D. J. Mavriplis, "Revisiting the Least-squares Procedure for Gradient
    Reconstruction on Unstructured Meshes," NASA/CR-2003-212683, 2003, p. 14.
