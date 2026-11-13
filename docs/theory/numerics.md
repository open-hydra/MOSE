# Spatial Discretization

MOSE uses a cell-centred finite volume method (FVM) on structured
multi-block grids.  This page describes the discretization framework,
the MUSCL reconstruction that provides second-order accuracy, the
slope limiters that keep each reconstructed face value between the values
of the two cells of the face on any mesh (van Albada excepted: see
[Flux Limiters](#flux-limiters)), and the shock-detection
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

Rotationally periodic faces join the two $k$ faces of one block, turned about the $x$
axis or translated (see
[Rotational Periodicity and Axisymmetric Wedge](numerics-wedge.md#rotational-periodicity-and-axisymmetric-wedge)).

### Cell geometry

A cell is the image of the unit cube by the trilinear map of its eight nodes. Each face is the bilinear surface through its four
nodes, the same for the two cells that share it, and its area vector, the one the fluxes take, is half the cross product of its
diagonals. The volume of the cell is the integral of the Jacobian of the map.

At set-up MOSE checks every cell and stops the run on a cell whose Jacobian is zero, folded at a corner or turned against its
block. The volume formula and the checks are in [Mesh Geometry and Block Interfaces](numerics-geometry.md).

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

3. Bound the limited slopes so that the value of each cell at both its faces
   lies between the value of the cell and that of the neighbour across the
   face, on a line of cells of any lengths: the change over half the cell is at
   most the jump to each neighbour [15],

$$
|\bar{s}_L|\,\delta l_L \le \min\bigl(|P_{i+1}-P_i|,\ |P_i-P_{i-1}|\bigr),\qquad
|\bar{s}_R|\,\delta l_R \le \min\bigl(|P_{i+2}-P_{i+1}|,\ |P_{i+1}-P_i|\bigr)
$$

   With $\phi(a,b) = \varphi(r)\,a$, $r = b/a$, the limiter takes, for $r > 0$,
   $\varphi \le \Delta l_a/\delta l$ and $\varphi \le r\,\Delta l_b/\delta l$,
   where $\Delta l_a$ and $\Delta l_b$ are the spacings of the slopes $a$ and $b$
   and $\delta l$ is the distance from the centre of the cell to the face: for
   the cell of length $h_i$ before the face, $\varphi \le (h_i + h_{i+1})/h_i$
   and $\varphi \le r\,(h_{i-1} + h_i)/h_i$. On a line of equal cells the bounds
   are $\varphi \le 2$ and $\varphi \le 2r$, the edges of the TVD region of Sweby's
   diagram [16], and the slope is the limiter's own; for
   data linear along the line $r = 1$, where no bound acts. The two faces of a
   cell take the same slope, except at the boundary faces below, which take a
   slope of their own. Minmod never exceeds the bounds and does not take them.
   The faces of symmetry planes and walls, whose ghost cell is the mirror
   image of the cell, and of the supersonic inflows (supersonic inlet, choked nozzle), whose
   ghost cell holds the inflow state and whose slopes are differences of index,
   take the bounds of a line of equal cells, 2 and 2: the limiter's own slope,
   which keeps their face value between the cell and the ghost cell.

4. Reconstruct interface values:

$$
P_L = P_i     + \beta\,\bar{s}_L\,\delta l_L,\qquad
P_R = P_{i+1} - \beta\,\bar{s}_R\,\delta l_R
$$

   where $\delta l_L,\,\delta l_R$ are half the lengths of cells $i$ and $i{+}1$
   along the line (the length of a cell along a direction is that of its mean
   edge in that direction), the spacings $\Delta l$ are the means of the lengths
   of two consecutive cells, and $\beta$ is a blending parameter
   controlled by the shock detector.

!!! note "Role of β"

    $\beta = 1$ recovers the full second-order MUSCL scheme.
    $\beta = 0$ reduces reconstruction to first order (donor-cell),
    which is used in the immediate vicinity of shocks.
    The shock-detection sensor (see below) sets $\beta$.

The reconstruction is one-dimensional along the grid lines: it treats the cell
centres of a line and the faces between them as points of a straight line,
spaced by the cell lengths. Where the centre of a face lies off the chord
between the two cell centres by a fraction of the cell size that refinement
does not reduce, the offset alone makes the face value of a linear field wrong
by that fraction of the cell size times the gradient across the chord: there
the reconstruction is of first order.

On grids whose cells are skewed the other
way of their neighbours the face centres lie off the chords on alternate sides,
the errors of the two faces of a cell add up, and the truncation error of the
cell does not converge. The error of the solution converges nonetheless: on
the zigzag grids of the
[manufactured solutions on skewed grids](../vv/2D-mms-skewed.md), whose face
centres lie up to 0.3 cell sizes off the chords at every refinement, the
observed $L_2$ orders on the finest pair are 0.97 to 1.10.

Diskin and Thomas
report the same behaviour for node-centred schemes on perturbed mixed-element
grids, where their edge-based flux integration is of first order: truncation
errors that do not converge, discretization errors that converge at first
order [18, Sec. VII.C]. On the same grids five of their six cell-centred
schemes, which reconstruct the face values from the value and the gradient of
each cell and take the flux at the face centre, keep the second order
[18, Sec. III, VII.C and X].

At the first interior face the slope behind the first cell, $s_0$ or $s_2$,
takes the ghost cell of the boundary face. Behind a zero-gradient face (unused
faces, and the inlets and outlets but the velocity inlet) the ghost value is that of the first cell: that slope is
zero, and every limiter gives zero with a zero argument, so the reconstruction
of the first interior face is of first order on the side of the first cell.

Behind a face with an assigned state the ghost cell is the first cell mirrored through
the face (see *Gradient Computation for Diffusive Fluxes*): it has the length
of the first cell across the face, which is then the spacing $\Delta l_0$ (or
$\Delta l_2$) of the slope behind the first cell, and its centre lies that far
from the centre of the first cell. Along a straight grid line the slope is then
that of the field for the quantities of the file; for the turbulence variables
and the soot, whose specific value is that of the first cell, the slope of
$(\rho\phi)$ is $\phi_1$ times that of the density.

!!! note "Choosing the reconstruction"

    MUSCL limits the slopes with one of the [flux limiters](#flux-limiters), van Leer
    by default; the first-order reconstruction sets every slope to zero (donor-cell).

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

Smooth, less dissipative than van Leer near local extrema. Unlike the others
it is not zero where the two slopes have opposite signs ($r < 0$): there the
bounds of step 3 do not apply, and the face value can lie outside the values
of the two cells of the face, on any mesh.

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

On a line of equal cells van Leer, MC, superbee and minmod lie in Sweby's TVD
region in the $(r,\,\phi)$ diagram; on other lines the bounds of the
reconstruction (step 3) keep them in the region that Zeng gives for lines of
cells of any lengths [16], and the forward Euler step of one-dimensional
advection is then TVD for Courant numbers up to 1/2 [15]:

| Limiter | Dissipation | Smoothness |
|---------|:-----------:|:----------:|
| MinMod | High | $C^0$ |
| Van Leer | Medium | $C^\infty$ |
| Van Albada | Medium-low | $C^\infty$ |
| MC | Low | $C^0$ |
| Superbee | Lowest | $C^0$ |

<figure>
  --8<-- "docs/theory/images/sweby.svg"
</figure>

---

## Shock Detection

Near strong shocks the MUSCL reconstruction can produce oscillations even with TVD limiters. MOSE uses a Jameson-type pressure sensor to detect shocks and degrade to first order in their vicinity.

Two sensors are available, the one below (Tramel) and that of Chen. HLLC+Tramel and HLLE++ take the sensor below, HLLC+Chen and AUSM+M that of Chen, and the other Riemann solvers none ($\beta = 1$) unless a sensor is chosen.

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

The Tramel sensor gives $\beta$ with the threshold $\Delta = 20$, a fixed constant:

$$
\beta =
\begin{cases}
1 - \tanh\!\bigl(10\,(\sigma/\Delta)^3\bigr)
  & \text{if } \sigma < 1/\Delta \\[4pt]
0 & \text{otherwise}
\end{cases}
$$

Where the first case applies $\sigma/\Delta < 1/\Delta^2$, so the tanh term is below $2\cdot10^{-7}$:

- $\beta = 1$, to within $2\cdot10^{-7}$, where $\sigma < 1/\Delta$ — full second-order reconstruction
- $\beta = 0$ elsewhere — first-order (no reconstruction)

This $\beta$ is passed directly to the MUSCL reconstruction step, where
it multiplies the limited slope.

!!! note "Solver-level integration"
      Some Riemann solvers also use the shock-detection flag internally to tune their numerical dissipation, providing an additional layer of robustness.

---

## Gradient Computation for Diffusive Fluxes

The velocity, temperature and species gradients of the viscous flux are computed with a 10-point stencil at each cell interface: the
two cells that share the face, and four cells along each of the two tangential directions. The gradient in computational space is
mapped to Cartesian coordinates with a face metric tensor $M_f$.

The default face metric takes the same differences of the cell centres as of the field, so the gradient of a linear field is exact on
skewed, stretched and non-smooth meshes. The alternative averages the metrics of the two cells: it is faster, and biased on stretched
meshes. The metric of the ghost cells, the exceptions at some boundaries and the evidence are in
[Face Gradients and Metrics](numerics-gradients.md).

## Ghost Cells of Walls and Extrapolation Faces

Behind a viscous wall and an extrapolation face the two layers of ghost cells are filled once per stage, before the fluxes, and the
viscous fluxes, the MUSCL reconstruction and the source terms of the turbulence models read the same values. What is extrapolated
are the quantities whose differences those terms take (the mass fractions, the temperature, the velocity, the specific turbulence and
scalar variables, the pressure); the ghost state follows from them. The rules for each boundary type are in
[Ghost Cells at Walls and Extrapolation Faces](numerics-ghost-cells.md).

## Rotational Periodicity and Axisymmetric Wedge

Rotational periodicity on the $k$ faces closes a mesh that turns about the $x$ axis: an axisymmetric wedge (one cell in $k$) or a 3-D
sector (several). The flow beyond one $k$ face is the flow next to the other turned by the angle $\delta$ of the wedge or sector.

MOSE integrates the 3-D equations on the wedge, with no axisymmetric source terms: the pressure and centrifugal terms come from the
fluxes of its two faces. The rotation of the ghost cells, the viscous terms of the wedge and the sector are in
[Axisymmetric Wedge and 3-D Sector](numerics-wedge.md).

## Faces Between Blocks

A face shared by two blocks (a block connection or a periodic connection) has one flux, computed once and given to both cells with
opposite signs, as an interior face. Which cell computes it, the shock sensor and the preconditioning velocity it takes, and when
the result equals that of the unsplit mesh, are in [Mesh Geometry and Block Interfaces](numerics-geometry.md#faces-between-blocks).

!!! info "Going deeper"
    - [Mesh Geometry and Block Interfaces](numerics-geometry.md): cell volume, checks of the cells, faces between blocks.
    - [Face Gradients and Metrics](numerics-gradients.md): the stencil, the face metric, ghost-cell metrics, chimera faces.
    - [Ghost Cells at Walls and Extrapolation Faces](numerics-ghost-cells.md): the extrapolation rules of each boundary type.
    - [Axisymmetric Wedge and 3-D Sector](numerics-wedge.md): rotational periodicity.

## References

[1] B. van Leer, "Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method," *J. Comput. Phys.*, 32(1), 1979.

[2] P. K. Sweby, "High resolution schemes using flux limiters for hyperbolic conservation laws," *SIAM J. Numer. Anal.*, 21(5), 1984.

[3] A. Jameson, W. Schmidt, E. Turkel, "Numerical solution of the Euler equations by finite volume methods using Runge–Kutta time-stepping schemes," AIAA-81-1259, 1981.

[4] R. Tramel, R. Nichols, and P. Buning. "Addition of improved shock-capturing schemes to OVERFLOW 2.1." 19th AIAA Computational Fluid Dynamics. 2009. 3988.

[15] M. Berger, M. J. Aftosmis, S. M. Murman, "Analysis of slope limiters on irregular grids," AIAA Paper 2005-0490, 43rd AIAA Aerospace Sciences Meeting, Reno, 2005, §III.A, Eq. (13)–(19): the TVD condition on the limited slope of a cell of a non-uniform grid, at most $\min\{2(u_{i+1}-u_i)/h_i,\ 2(u_i-u_{i-1})/h_i\}$ for monotone increasing data (Eq. 19, a maximum value as p. 8 reads it), with the forward Euler step of the upwind scheme for Courant numbers up to 1/2 (Eq. 17).

[16] X. Zeng, "A general approach to enhance slope limiters on non-uniform rectilinear grids," arXiv:1301.0967v2, 2014, Eq. (2.13), (3.11) and §3.3 (Sweby's diagram on non-uniform grids).

[18] B. Diskin, J. L. Thomas, "Comparison of node-centered and cell-centered unstructured finite-volume discretizations: inviscid fluxes," 48th AIAA Aerospace Sciences Meeting, Orlando, FL, 2010 (NASA NTRS 20100003409).
