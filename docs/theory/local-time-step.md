# Local Time Step: Derivation and Evidence

This page backs the [local time step](time-integration.md#local-time-step) rule and the limits of its table with derivations and
evidence. Read the rule and the recommended values there first; come here to see why the limits hold, and where they stop holding.

## Why the convective limits add up

Without dimensional splitting the convective limits of the directions add up: an explicit scheme
on square 2-D cells needs a condition of the type $\tau(|a_x|/\Delta x + |a_y|/\Delta y) \le C_0$, not
the largest of the two terms ([8], §9.3, which states the form without a proof). For forward
Euler the form follows, with the diffusive rates added as well, from the positivity condition (6.1)
of [8]: with first-order upwind convection and central diffusion on a Cartesian mesh the
diagonal entry of a row is $-\sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr)$ and the others are not negative, so
$\tau \sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr) \le 1$: the [local time step](time-integration.md#local-time-step) rule with CFL = 1 and VNN = 1/2.

## Where the constants come from

The constants of Heun and RK3 in the [table of limits](time-integration.md#local-time-step) do not come from that source: they are the von Neumann limits
of the linearised MUSCL–upwind advection with central diffusion, proven for MOSE's stencil in 1-D and
carried to 2-D and 3-D, for any split of the Courant number, by the convexity of the symbol; they hold
for the van Leer slope (the default), van Albada and MC, not for minmod and superbee. The source warns
that such an eigenvalue criterion can mislead for the non-normal operators that boundaries give, as an
inflow boundary with upwind differences ([8], §3.3, pp. 26–27): the limits are those of the interior.

## Proof of the limits

The
van Leer slope (the default; van Albada and MC linearise the same way) linearised in smooth flow is
the mean of the two one-sided slopes (Fromm's scheme), with the symbol
$g(\theta) = \tfrac12[(1-\cos\theta)^2 + i\sin\theta(3-\cos\theta)]$, $z = -\nu g$; the constants
hold for these three limiters, not for minmod and superbee, which take the smaller or the larger of
the two slopes and have no linearisation where they are equal. Heun is stable if and only if
$\nu \le 1$: $|R(z)|^2 - 1 = \tfrac{\nu}{16}(1-c)^2[\nu(\nu a - 2b)^2 + 4\nu b^2 - 16]$ with
$c = \cos\theta$, $a = 5 - c^2$, $b = 1 - c$, increasing in $\nu$ and equal to
$\tfrac{\nu}{16}(1-c)^4(1+c)(c-3)$ at $\nu = 1$.

RK3 is stable if and only if
$\nu \le \nu^\ast = 1.1757678815978$, a root of a polynomial of degree 24 (a cylindrical
decomposition in exact arithmetic), and with the first-order flux if and only if
$\nu \le 1.2563726633$, half the real stability interval $2.5127453$; diffusion gives VNN $\le 1/2$
and $0.6282$ (the table rounds them down). The curve of $-g$ bounds a convex set $K$
($x'y'' - y'x'' = -\tfrac12(1-\cos\theta)(\cos^2\theta + \cos\theta + 4) \le 0$ along it), so a sum
over directions is a point of $\nu K$, inside the stability region by the maximum principle: the
limits of 1-D hold in 2-D and 3-D for any split of the Courant number, with no need for a convex
stability region.

With diffusion the additive rule is stable for Heun with CFL = 1 and VNN = 1/2,
and exact at every Péclet number (the checkerboard mode reaches $z = -2$); for RK3 the region is not
convex, and the convex hull of $1.17K$ and $[-2.48, 0]$ (and of $1.1757K$ and $[-2.512, 0]$) is
shown inside it with a certificate in exact rational arithmetic.

## Skewed cells in 3-D

In 3-D, on a cell whose three directions are skewed against each other, the viscous rate is
multiplied by $F = (2 + m + 1/m)/4$: with $a$, $b$, $c$ the cosines of the angles between the three
$\nabla\xi_d$, where $a^2 + b^2 + c^2 + 2abc > 1$, $m + 1$ is the largest eigenvalue of the
$3\times3$ matrix of the cosines ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$).

There the spectral radius of the viscous operator exceeds $4D_{\max}\sum_d|\nabla\xi_d|^2$ by a
factor below $9/8$, which $F$ bounds (exactly on rhombohedra); elsewhere, orthogonal cells and cells
with one direction orthogonal to the other two included, $F = 1$, and in 2-D the rate holds as it
is.

## The axisymmetric wedge

On the 2-D axisymmetric wedge (rotational periodicity or symmetry on the $k$ faces) the stress on the wedge faces adds to the
rows of $u_r$ and $u_\theta$ of every cell terms that grow as $1/r^2$ next to the axis: in the row of
$u_r$ the diagonal $-\tfrac{4}{3}\nu_e/r^2$ and couplings $\tfrac{2}{3}\nu_e|\nabla\xi_d\cdot\mathbf{e}_r|/(2r)$
to the $u_r$ and $\tfrac{2}{3}\nu_e|\nabla\xi_d\cdot\mathbf{e}_x|/(2r)$ to the $u_x$ of the neighbours $d\pm1$, in the
row of $u_\theta$ the diagonal $-\nu_e/r^2$ and couplings $\nu_e|\nabla\xi_d\cdot\mathbf{e}_r|/(2r)$, with
$\nu_e = (\mu_\ell + \mu_t)/\rho$ and $r$ the distance of the cell centre from the axis (no coupling in the row of $u_x$). The diffusive rate takes the increase of the row sum they give,

$$
s_h = \nu_e \max\!\left( \frac{4/3}{r^2} + \frac{2}{3r}\sum_d \bigl(|\nabla\xi_d\cdot\mathbf{e}_x| + |\nabla\xi_d\cdot\mathbf{e}_r|\bigr),\;
\frac{1}{r^2} + \frac{1}{r}\sum_d |\nabla\xi_d\cdot\mathbf{e}_r| \right),
\qquad
\frac{1}{\Delta t_i} \mathrel{+}= \frac{s_h}{4\,\text{VNN}} ,
$$

so that $\Delta t\,(4D_{\max}\sum_d|\nabla\xi_d|^2 + s_h) = 4\,\text{VNN}$. The stability limits of the table are proven away
from the axis on uniform meshes (in 2-D for any skew and aspect ratio, in 3-D with the factor $F$ above) and, next to the
axis, for the scalars; for the momentum next to the axis they rest on evidence (see [What supports the limits](#what-supports-the-limits)).

Next to the axis of square cells $s_h = 8\,\nu/\Delta r^2$; far from it $s_h$ vanishes as $1/r$. The wedge stress, and with it $s_h$, applies
with either sign of the wedge angle (a right-handed mesh with the axis on the last $j$ face has a negative one). The step
takes $s_h$ on every 2-D mesh with a wedge angle, whatever the boundary of its $k$ faces: the viscous flux of rotationally
periodic faces and that of symmetry faces, whose mirror image of a cell is, without swirl, the cell turned by the wedge
angle, both carry this stress; unused $k$ faces carry none, and there $s_h$ only makes the step smaller than it needs to
be.

## What supports the limits

- Away from the axis, on a uniform mesh of equal cells (parallelograms in 2-D, parallelepipeds in 3-D) with frozen
  coefficients and no boundaries, the viscous operator of MOSE has a real symmetric Fourier symbol. With the face gradient of
  MOSE (the difference across the face, the mean of the central differences of the two cells along it, mapped with the face
  metric, which on such a mesh is the cell metric for every choice of face metric), a mode $e^{i\boldsymbol\theta\cdot\mathbf c}$
  gives $\rho\,\dot{\hat{\mathbf u}} = -S\hat{\mathbf u}$ with $S = \mu\,\mathrm{tr}(B)\,I + (\mu+\lambda)B$, $B = G\Theta G^T$,
  the columns of $G$ the rows $\nabla\xi_d$ of the cell metric, $\Theta_{dd} = 4\sin^2(\theta_d/2)$,
  $\Theta_{de} = \sin\theta_d\sin\theta_e$; a scalar of diffusivity $D$ has the symbol $D\,\mathrm{tr}(\Theta\Gamma)$,
  $\Gamma_{de} = \nabla\xi_d\cdot\nabla\xi_e$ (the cross metric terms of skewed cells included).
  $\Theta = \mathbf c\mathbf c^T + \mathrm{diag}\bigl((1-\cos\theta_d)^2\bigr)$, $c_d = \sin\theta_d$, is positive
  semidefinite: the eigenvalues are real and not positive, the operator is normal, and the eigenvalue criterion is von
  Neumann's ([8], §2.4, eq. (2.10), p. 18, and §3.3, pp. 25–28).
- In 2-D $\Theta \le 4I$: for any skew angle and aspect ratio the spectral radius is exactly $4D\sum_d|\nabla\xi_d|^2$ for a
  scalar and $4\nu\bigl(\sum_d|\nabla\xi_d|^2 + \lambda_{\max}(\Gamma)/3\bigr) \le \tfrac{16}{3}\nu\sum_d|\nabla\xi_d|^2$ for the
  momentum, so the rate $4D_{\max}\sum_d|\nabla\xi_d|^2$ bounds it, with equality for the scalar of largest diffusivity (the
  energy in air).

  In 3-D the sharp constant is $9/2$: with $a$, $b$, $c$ the cosines of the angles between the three
  $\nabla\xi_d$, the rate bounds the spectral radius if and only if $a^2 + b^2 + c^2 + 2abc \le 1$ (true for orthogonal
  cells, for cells with one direction orthogonal to the other two, and for a rhombohedron up to the edge angle
  $\arccos(-1/3) = 109.47^\circ$); elsewhere the radius exceeds it by a factor below $9/8$, and the local step multiplies the
  rate by $F = (2 + m + 1/m)/4$, where $m + 1 \in (2, 3)$ is the largest eigenvalue of the $3\times3$ matrix of the cosines
  ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$): a proven bound of the factor, exact on rhombohedra
  (1.04 at $115^\circ$, 1.11 at $119.5^\circ$), below $9/8$ on every cell, and 1 where the rate already holds.

  The real
  stability intervals of Heun, $[-2, 0]$, and of RK3, $[-2.5127, 0]$, then give VNN = 1/2 and 0.62 ($2.5127/4 = 0.628$).
- Next to the axis the operator is not translation invariant. For the scalars (energy, species, turbulence variables) it stays
  self-adjoint in the inner product weighted by the cell volumes: the ghost cells of the wedge copy the scalars and the wedge
  faces do not update them, so on a mesh of equal rectangular cells in $(x, r)$ each face couples two cells symmetrically, and
  the row sums are at most $4D\sum_d|\nabla\xi_d|^2$ (up to a relative $\varepsilon/6$ for an axis at $\varepsilon\,\Delta r$,
  covered by $s_h$ for diffusivities below $0.6\,\nu_e/\varepsilon$).

  A self-adjoint operator whose row sums are within the
  rate of each cell is stable with the local step, in the norm weighted by $V/\Delta t$, and the limits of the table hold for
  the scalars on every row. The momentum operator is not self-adjoint in any inner product with cell weights: the cross terms
  $(\lambda+\mu)\,\partial^2/\partial x\,\partial r$ of the faces in $x$ and in $r$ are weighted by the radius of the cell,
  $(j-\tfrac12)\Delta r$, and of the face, $j\Delta r$.

  In the first row the rows of $u_r$ and $u_\theta$ are diagonally
  dominant (their diagonals carry the wedge terms $-\tfrac43\nu_e/r^2$ and $-\nu_e/r^2$); the row of $u_x$ is not: next to the
  axis (doubled weight of the upper face, mirror image across the axis) it couples to $u_r$ of four neighbours with
  $\nu/(3\Delta x\,\Delta r)$ each, four times the cross coupling of the plane operator, a deficit of
  $\tfrac43\nu/(\Delta x\,\Delta r)$; its row sum stays below the rate of the cell.

  No proof covers the momentum next to the
  axis. The support is evidence: on the operator of the code for cells with $dx/dr$ from 0.01 to 128 and wedges of 0.25° to 5°,
  its symmetric part is negative semidefinite, its eigenvalues are real to $6\cdot10^{-4}$, and the largest stable VNN with the
  local step is 0.500 to 0.577 with Heun and 0.628 to 0.725 with RK3, approaching the limits of the table, 1/2 and 0.62
  (0.628), from above for very long or very short cells; on a diffusion-dominated 1° wedge of the tests
  (32 × 32 cells, 2000 local steps) the largest stable VNN is about 0.5 with Heun (0.508 to 0.519)
  and 0.63 to 0.65 with RK3 for $dx/dr$ = 1, 4 and 8 at Prandtl numbers 0.72 and 0.9 (and 10 at $dx/dr$ = 4), and the
  test runs 2000 steps at the limits and at the recommended values. The eigenvalue criterion can mislead for non-normal operators
  ([8], pp. 26–27).

## Positivity with the limiter

The positivity limits of its last column hold on a mesh of equal rectangular cells in smooth flow, the one of the
diffusion with a constant diffusivity.

These limits are linear (von Neumann, smooth solutions). With the limiter, positivity (a maximum principle:
in more than one dimension a second-order scheme is not TVD, Goodman et al. 1985 as cited in [9], p. 15 of the
report) is guaranteed only for smaller steps, for the
advection of a scalar with the limited flux: on a line of equal cells the van Leer slope of MOSE is the flux limiter
$\psi(\theta) = \theta/(1+\theta)$, and with the bounds of the reconstruction
([Spatial Discretization](numerics.md#muscl-reconstruction), step 3) the limited flux of every limiter but van Albada
satisfies $0 \le \psi \le 1$ and $0 \le \psi/\theta \le \mu$ with $\mu = 1$ on a line of cells of any lengths.

With
Heun and RK3 (SSP coefficient 1) positivity then holds for $(1+\mu)\,\nu \le 1$: on any mesh, when at every stage the
Courant number of the flow leaving each cell (the step times the volume flow out of the cell through its faces, over
its volume) is at most 1/2. Central diffusion of a scalar is positive for $\text{VNN} \le 1/2$ on a mesh of equal
rectangular cells with a constant diffusivity, where the face gradient has no cross metric terms (the condition (6.1) of
[8] above).

On other rectangular cells the bound can be lower for a cell that advances with the step of the rule (a
local step, or a level of the multirate; with the global step of a time-accurate run, the smallest of the mesh,
$\text{VNN} \le 1/2$ keeps every cell positive in one dimension with a constant diffusivity): the rule takes the width
and the largest diffusivity of the cell itself, while the flux of a face takes the distance between the two
centres (with the default stencil face metric) and the diffusivity of the mean state of the two cells,
which can exceed that of the cell; in one dimension, with a constant diffusivity and no convection, the positivity of
every stage asks $\text{VNN}\,h\,(1/\delta_- + 1/\delta_+) \le 1$ (the forward Euler step; Heun and RK3 inherit it as
convex combinations of such steps), with $h$ the width of the cell and $\delta_\mp$ the distances to
the centres of its two neighbours: $\text{VNN} \le 1/2$ on equal cells and on widths in a geometric progression, and
towards 1/4 in a cell much larger than both neighbours.

On skewed cells the face gradient also takes the central
differences of the two cells along the face, whose cross metric terms can give corner neighbours of a
cell a negative weight: on a mesh of equal parallelograms two corner neighbours of every cell have one, at every skew
angle, and no step makes the diffusion positive there. With the additive rule
both hold for CFL ≤ 1/2 and VNN ≤ 1/2 in smooth flow on a mesh of equal rectangular cells with a constant diffusivity,
where the flow leaving a cell through its faces is about the
convective rate of the rule (at captured shocks see below) ([8], §5.2, §6.2 and §6.4; experiments reported there give
positive results up to $\nu = 1$ in 1D and about $2/3$ in 2D).

The bound is lower in the first cell next to a wall of a
RANS run, whose ghost cell holds $-k_1\rho_G$ (and $-\tilde\nu_1\rho_G$; the wall value of $k$, and of $\tilde\nu$ on
a smooth wall, is zero), $\rho_G$ the density of the ghost cell: there it is $1/(2 + \rho_G/\rho_1)$, and in smooth
flow CFL ≤ 1/2 keeps it while the Mach number of the flow normal to the wall is below $2\rho_1/\rho_G$. Next to faces whose ghost
cells extrapolate the cells, $3q_1 - 3q_2 + q_3$ (the velocity inlet, and the flow and soot variables of a solid
propellant grain that blows its products in through a wall), positivity is not guaranteed.

## Shocks and the HLLC flux

The local time step takes the wave speeds of the cell alone. With a stationary shock captured on one face, the fastest
wave that the HLLC flux lets into a cell on the slow side of the shock exceeds the cell's own $|u_n| + a$ by up to 1.47
(a stretched cell crossed obliquely by a strong shock; 1.34 for a normal shock aligned with the cells, 1.25 for the sum
over the two directions of a square cell): with CFL = 0.8 the face Courant number stays below 1.256, the first-order
limit of RK3 that applies where the limiter has removed the slope.

The local time step takes the wave speeds of the cell alone, $|\mathbf{v}\cdot\nabla\xi_d| + a\,|\nabla\xi_d|$. The HLLC
flux lets waves enter a cell through a face with the speeds of Batten et al. (SIAM J. Sci. Comput. 18, 1553–1570, 1997, doi:10.1137/S1064827593260140, eq. 51), which
include the Roe-averaged state of the face: on the slow side of a shock the fastest wave entering a cell is faster than the
cell's own $|u_n| + a$, so the Courant number of that face exceeds the CFL number.

For a gas with $\gamma = 1.4$ and a
stationary shock captured on one face, the ratio is 1.27 at Mach 3 and at most 1.34 for a normal shock aligned with the cells,
at most 1.25 for the sum over the two directions of a square cell, the quantity the time step rule uses, and up to 1.47 in one
direction of a stretched cell whose face is inclined to a strong shock (1.32 at Mach 3, 1.45 at Mach 10).

With CFL = 0.8 the
face Courant number is thus at most 1.08, 1.00 and 1.18 respectively: below the linear stability limit of RK3 with the
first-order upwind flux, 1.256, the limit that applies where the limiter has removed the slope (the limit of the limited scheme
itself is not derived), and below the limit with the smooth MUSCL reconstruction, 1.176, except on stretched cells crossed
obliquely by shocks stronger than about Mach 50 (1.177 in the limit of a strong shock).

Lower $\gamma$ raises the ratios (1.49
instead of 1.34 for $\gamma = 1.2$). With strong shocks keep CFL ≤ 0.8 for RK3. The positivity of the first-order Euler step
with the HLLC flux, guaranteed by CFL ≤ 1/2 in smooth flow, also holds at shocks on square cells (the sum of the incoming waves
is at most 1.90 times the cell's own rate), but not on stretched cells crossed obliquely by a strong shock, which need
CFL ≤ 0.41 at Mach 10; with MUSCL it is observed, not proven.

## Reading the recommended values

The recommended values are in [Local time step](time-integration.md#local-time-step).

These are
values of the scheme, not of a case: a case that runs stable beyond them owes it to its mesh (a non-uniform spacing around the
cell that sets the step lowers the spectral radius of its diffusion), and its threshold is not a limit to carry elsewhere. The
viscous rate takes the largest diffusivity, that of the energy, and adds up the directions: a VNN tuned on a rule that took
the momentum diffusivity and the largest direction alone gives here a smaller step (by $\gamma/\text{Pr}_t$, 1.65 with the
default $\text{Pr}_t = 0.85$, in a
stretched cell where the turbulent viscosity sets it), and can be raised within the limits of the table.

## References

[8] W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam, 2000.

[9] E. M. Constantinescu, A. Sandu, "Multirate timestepping methods for hyperbolic conservation laws," *J. Sci. Comput.*, 33(3), 239–278, 2007, doi:10.1007/s10915-007-9151-y; the page cited is that of the authors' report, Virginia Tech, https://hdl.handle.net/10919/19488.
