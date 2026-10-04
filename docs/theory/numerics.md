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

### Rotational periodicity (type 200)

Type 200 on the $k$ faces closes a mesh that turns about the $x$ axis: an
axisymmetric wedge (one cell in $k$) or a 3-D sector (several). The flow
beyond face 5 is the flow next to face 6 turned back by the angle $\delta$
of the wedge or sector, and the flow beyond face 6 the flow next to face 5
turned forward. With $g = +1$ on face 5 and $g = -1$ on face 6, the ghost
cell takes the state of its source cell (on a wedge the cell itself, on a
sector the cell next to the opposite $k$ face) with the velocity turned
about $x$ by $\varphi = -g\,\delta$:

$$
v_g = v\cos\varphi - w\sin\varphi = v\cos\delta + g\,w\sin\delta, \qquad
w_g = v\sin\varphi + w\cos\varphi = -g\,v\sin\delta + w\cos\delta .
$$

$u$ and the scalars do not change under a turn about $x$, and are copied.

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
Cartesian coordinates $(x, y, z)$ with a face metric tensor $M_f$. Each cell
stores its metric $M = A^{-T}$, where the rows of $A$ are the cell edge vectors
$\partial\mathbf{x}/\partial\xi$, $\partial\mathbf{x}/\partial\eta$,
$\partial\mathbf{x}/\partial\zeta$; the edge lengths also give the physical
spacings $\Delta l$ of the MUSCL reconstruction. The face metric combines the
two cells $L$ and $R$ that share the face, selected by `diffusive-metric` in
`[MOSE-Numerics]`:

| `diffusive-metric` | Face metric | |
|---|---|---|
| `inverse-mean` (default) | $M_f = \left[\tfrac12\left(A_L + A_R\right)\right]^{-T}$ | the face-normal edge is the centre-to-centre vector: exact for a linear field on any stretching |
| `mean` | $M_f = \tfrac12\left(M_L + M_R\right)$ | faster; biased on stretched meshes |

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

Block-connection and chimera faces always use the inverse mean: across a large
cell-size jump at a block interface, the mean makes the normal gradient so steep
that the explicit diffusion becomes unstable.

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
    The tangential derivatives on a face are central differences over two cells,
    divided by the cell's own tangential edge. With stretching in a tangential
    direction they carry the same $(1+r)^2/(4r)$ factor, which neither option
    corrects. In a thin boundary layer this enters only through the streamwise
    variation of the shear on faces normal to the flow, so its effect is much
    smaller than that of the face-normal gradient.

---

## Axisymmetric Wedge (BC 200)

A 2-D axisymmetric case is a wedge of the 3-D mesh: one cell in $k$, whose two
$k$ faces lie in the meridional planes $\theta = \mp\delta/2$ and carry the
boundary type `200`. MOSE integrates the 3-D equations on the wedge; there are
no axisymmetric source terms. The ghost cell across each wedge face is the cell
itself rotated by $\pm\delta$ about the $x$ axis, so every face flux follows
from the state of the cell:

- inviscid: the pressure and the momentum carried by the cell state rotated to
  the face (no Riemann problem). In the radial momentum the two faces give the
  $+p/r$ and the centrifugal terms of the axisymmetric equations;
- viscous: the flux of a face between the cell and its rotated image, with the
  stencil of an interior face (the $i$ and $j$ neighbours of the image rotated
  in the same way) and the face metric of the block connections.

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

Next to the axis the hoop stress also stiffens the explicit step. In the radial
momentum of the first row of cells, at $r \approx \Delta r/2$, it is a damping
$(2\mu+\lambda)\,u_r/r^2 \approx \tfrac{16}{3}\,\mu\,u_r/\Delta r^2$ ($\mu$ with
the eddy viscosity), with couplings to the neighbours through $\lambda\,\nabla\cdot\mathbf{u}$.
The local time step includes them ([Time Integration](time-integration.md#local-time-step)),
so the `vnn` limits of the rule hold next to the axis too; without them, at a
Prandtl number of 10, the first rows broke down at the `vnn` the rule declares
stable (test `WedgeAxisStep`).

The 3-D rotational periodicity (type `200` with more than one cell in $k$)
computes the inviscid flux only. On a planar 2-D mesh (faces 5 and 6 parallel,
no wedge angle) type `200` makes the image the cell itself, and the fluxes of
the two faces cancel.

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
the ghost cell, with its value at the start of the stage; before the flux of
these faces was computed once, the ghost held $U_r = 0$ and the smoothing
ignored that neighbour. Neither is the unsplit mesh, where the smoothing sees
the new value of the neighbour, so with the smoothing on a split mesh does not
give the fluxes of the unsplit one bit for bit.

The viscous flux of a face between blocks uses its own stencil (the two
tangential gradients of the two cells, averaged) and the inverse-mean face
metric, not those of an interior face.

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
