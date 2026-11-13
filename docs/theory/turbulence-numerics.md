# Turbulence: Numerical Treatment

This page details how the turbulence equations are treated numerically: the ghost cells and the cell gradients that the source terms
read, the wall distance, and the advance of the stiff source terms in steady and time-accurate runs. The models are in
[Turbulence Modelling](turbulence.md).

## General RANS Features

### Wall ghost cells of the turbulence variables

Behind a viscous wall the ghost cell of each turbulence variable takes the value that puts the wall
value on the face,

$$
Z_\text{ghost} = 2\,Z_\text{wall} - Z_\text{cell}, \qquad Z = \tilde\nu,\ k,\ \omega ,
$$

with the wall values of the models in [Turbulence Modelling](turbulence.md), here of the specific variables $Z = \tilde\nu, k, \omega$
(the ghost value of $\rho Z$ is that times the density of the ghost cell). It is the rule of CFL3D,
«the actual wall boundary conditions for the turbulence quantities (unlike the primitive variables)
are applied at ghost cells» (S. L. Krist, R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual
(Version 5.0)*, NASA/TM-1998-208444, 1998, appendix H, eq. H-186).

MOSE sets these
values once per stage, from the wall values of the wall flux routines (the same wall state), with
the other ghost values of the wall ([Numerics](numerics.md), *Ghost Cells of Walls and Extrapolation Faces*):
the fluxes through the first interior faces and the source terms read the same ghost cells, at every
stage. The source terms read
them in the gradients of the wall cells: the $c_{b2}$ term of SA and its term
$-(\nu+\tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu/\sigma$, $F_1$ and the cross-diffusion of SST, the
cross-diffusion of Wilcox 2006.

With $\tilde\nu = a\,y$ next to a wall the $c_{b2}$ term of the wall
cell sees $\partial\tilde\nu/\partial y = a$ (manufactured-solution tests); a ghost cell that copied the cell
would give $a/2$. A solid propellant grain that blows its products in through a wall takes the
same rule for its turbulence variables, with the wall values of its flux: those of the blowing wall while
it burns, of the wall at the flame temperature when it does not; a grain modelled as an inflow is a
mass-flux inlet.

### Cell gradients of the source terms

The source terms take, in every cell, the central differences of the values along the three index
directions and map them to Cartesian components with a cell metric. With
the stencil metric it is $M_s = (J_c^{T})^{-1}$ with
$J_c^{T}(:,d) = (\mathbf{x}_{c+e_d} - \mathbf{x}_{c-e_d})/2$, the same central differences of the
cell centres (the mean of the eight nodes): the gradient of a linear field is exact on stretched,
skewed and non-smooth cells (along an index direction with one cell, the edge of the cell and no
difference).

On a mesh of equal parallelograms these gradients are the reference gradient of [14]
(Eq. 1), second order on smooth structured meshes. In a wall cell, for a field that is not linear,
the reflected ghost value $2\phi_\text{wall} - \phi_\text{cell}$ makes the gradient first order, as
the Green–Gauss gradient is in the boundary cells of [14] (Eq. 13).
The ghost cells sit where their value is right: those of the fluxes, at the positions
of the face metric ([Numerics](numerics.md)); behind a wall the position of the class of the
quantity, one metric per group: the velocity, the density (whose ghost value follows the temperature:
exact where the pressure and the composition do not vary across the wall cell), the turbulence
variables (a metric of their own only next to gas–surface interaction walls, where the velocity is
extrapolated and they take their wall value).

Behind a zero-gradient face (unused faces, and the inlets
and outlets but the velocity inlet), whose ghost cell copies the
cell, the position is the cell itself, so the difference across the boundary cell is one-sided. The
strain gradient of SA-RC takes, behind the faces between blocks and of a turned sector, the position
of the cell whose velocity gradient the ghost cell holds, and the cell itself behind the other faces
(see [SA-RC](turbulence-variants.md#sa-rc-spalartshur-rotation-curvature-correction)).

The metrics are computed once at setup and serve the source terms of SA, SA-RC, SST, Wilcox 2006 and
SSG/LRR and the eddy viscosity of the output. With the mean metric it is that of the cell:
the central difference across a cell next to a jump of size $r$ is then wrong by $(r-1)/4$, on a
skewed wall cell a copied density takes a normal gradient it does not have (manufactured-solution tests),
and next to a zero-gradient face the normal component is half the one-sided difference.

### Blowing / wall-suction correction

A wall value of $\omega$ for walls that inject or suck gas, for SST and Wilcox 2006, is
implemented but not offered as a model choice.

### Wall distance

The distance $y$ of the models (SA, the SST blending functions, the wall
value of $\omega$) is the minimum distance from the cell centre, the mean of
its eight nodes, to the wall surface, as the NASA Turbulence Modeling
Resource requires (ref. 13). It is not the distance to the nearest wall node
or face centre, which differs from it wherever the grid lines are not normal
to the wall; next to a convex corner or at the end of a wall it is the
distance to the corner or the end.

The distance to the centre of the nearest wall face is wrong by a fixed
fraction in the first cell of a skewed mesh. With the grid lines at an angle
$\theta$ to the wall normal, the centre of the first cell, at the height
$d_1$, lies $d_1\tan\theta$ along the wall from the centre of its wall face,
the nearest face centre while $d_1\tan\theta$ is at most half the face, at
the distance $d_1/\cos\theta$. The relative error,
$1/\cos\theta - 1$, does not fall with refinement: 1.5 % at 10°, 6.4 % at 20°,
15.5 % at 30°, 41 % at 45°. It enters squared the destruction of SA,
$c_{w1} f_w (\tilde\nu/d)^2$ (a factor $\cos^2\theta$: 0.97, 0.88, 0.75,
0.50), and the arguments of $F_1$ and $F_2$ and the wall value
$\omega_w = 60\nu/(\beta_1 d_1^2)$ of SST. That distance is exact only on
meshes orthogonal to the wall, with the cell centres above the face centres.

- Each face of a wall boundary (viscous walls, coupled walls, a grain that blows through a wall, gas–surface
  interaction walls) is the bilinear
  surface through its four nodes, the face of the cells
  ([Cell volume](numerics-geometry.md#cell-volume)):
  $\mathbf{x}(s,t) = \mathbf{x}_1 + s\mathbf{a} + t\mathbf{b} + st\mathbf{T}$,
  $\mathbf{a} = \mathbf{x}_2 - \mathbf{x}_1$, $\mathbf{b} = \mathbf{x}_4 - \mathbf{x}_1$,
  with the twist $\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$
  (corners in cyclic order). The distance to it is its minimum over the face,
  exact to the rounding of the size of the face and of the distance: a face plane
  to the resolution of its coordinates (every face of a 2-D mesh) takes the four
  triangles between its centre and its sides, any other face a certified minimum
  ([Code Structure](../development/structure.md#wall-distance) gives the
  algorithm). The four triangles of a twisted face would be wrong by up to
  $(|\mathbf{T}|/16)/y$, 6 % on a unit face twisted by 0.1 at the height 0.1.
- On a 2-D axisymmetric wedge the distance is measured in the
  meridional plane $(x, r)$, $r = \sqrt{y^2+z^2}$, from the cell centre to the
  wall curve. This is the distance to the wall surface of revolution, and it
  does not depend on the angle of the wedge. A wall on the faces of the wedge
  (faces 5–6) stops the run. It is the one exception to the single geometric
  model of the cells: for the wall distance the wall of the wedge is the
  surface of revolution of its meridional curve, for the volumes and the fluxes
  the cells are polyhedra with plane faces; the two differ by
  $r(1 - \cos(\delta/2)) \approx r\delta^2/8$.
- The distance and the face that gives it do not depend on the threads or on
  the ranks.
- Across periodic boundaries the wall goes on: the distance is the minimum over
  the periodic images of the wall too. A periodic connection gives
  the translation from its source face to its own face, and on a 3-D mesh
  without sectors so does a block closed in $k$ by rotational periodicity (its last node
  plane is the first one translated). The images are the integer combinations of
  these translations (along $x$ on an axisymmetric wedge). The two sides of a
  periodic boundary match only to the round-off of the coordinates in the mesh
  file: translations that
  differ by at most $10^{-2}$ of the thinnest cell across the periodic
  boundaries are the same, as in the type of the mesh. The thickness of a
  cell across a periodic face is the distance from the centre of the face to
  the centre of the opposite face of the cell, along the normal of the face;
  for a block closed in $k$ by rotational periodicity, the move of its nodes per cell in
  $k$. The round-off of a coordinate grows with its size, so a tolerance
  relative to the size of the coordinates would depend on where the mesh
  lies; relative to the cells it does not. On a 3-D mesh whose
  sectors are closed by the rotational periodicity (on faces 5–6)
  the images are also the turns by multiples of the sector angle about the
  $x$ axis, up to half a turn either way. A channel periodic in $x$ and turned by
  30 degrees about $z$ gives translations that differ by 1.1e-4 of its cells with
  five significant digits, by 8.8e-4 with four, and by 0.15 with two, which stops
  the run. Translations that do not form a lattice stop the run: one that is not
  an integer combination of the others, more than three independent ones, or one
  not along $x$ on a wedge or with sectors. Only a bounded set of images can hold
  the nearest wall, and only those are searched: a mesh without periodic
  boundaries costs nothing more.
- Each cell takes the sand-grain roughness $k_s$ of the face that gives the
  minimum distance; between faces at the same distance, the first one in the
  boundary file.

---

## Numerical treatment of the source terms

The turbulence source terms are added to the residual and advanced with the same explicit
Runge–Kutta stages as the mean flow. Their destruction parts are stiff near walls: for the
$k$–$\omega$ models $\omega$ grows like $1/y^2$ towards the wall, and an explicit stage on
$\dot q = -D(q)$ gives $q - \Delta t\,D = q\,(1 - \Delta t\,D/q)$, negative for $\Delta t\,D/q > 1$. Here
$q$ is the conserved turbulence variable; the destruction rate $D/q$ is at most $d$, the derivative
of the destruction ($d/2$ for $\tilde\nu$ and $\omega$, whose destruction is quadratic):

| Model | Equation | $d = \partial D/\partial q$ |
|-------|----------|------------------------------|
| SST / Wilcox | $k$ | $\beta^\ast\,\omega$ |
| SST / Wilcox | $\omega$ | $2\,\beta\,\omega$ |
| SA | $\tilde\nu$ | $2\,c_{w1}\,f_w\,\tilde\nu/y^2$ |
| SSG/LRR | $R_{ij}$ | $C_1\,C_\mu\,\omega + \tfrac12 C_1^\ast\max(P_{kk},0)/k$ (dissipation and return to isotropy) |
| SSG/LRR | $\omega$ | $2\,\beta_\omega\,\omega$ |

($y$ the wall distance; $f_w$ is held fixed.) Steady and time-accurate runs treat this
stiffness in two different ways.

### Steady runs: point-implicit factor

With the point-implicit factor (the default) the increment of each turbulence equation is divided
by $1 + \Delta t\,d$, with the local time step $\Delta t$ of the cell:

$$
q^{n+1} = q^n + \frac{\Delta t}{1 + \Delta t\,d}\,\bigl( F(q^n) + S(q^n) \bigr),
$$

where $F$ stands for the fluxes and $S = P - D$ for the source (in each Runge–Kutta stage, with the
residual of the stage). This is the point-implicit form of Bussing and Murman,
$D\,\Delta q = \Delta t\,(F + S)$ with $D = 1 - \Delta t\,\partial S/\partial q$ approximated by its
destruction part $1 + \Delta t\,d$, in its first-order form, «suitable for steady-state
calculations» (H. C. Yee, J. L. Shinn, NASA TM-89415, 1986, eq. 3.1).

Like any approximation of the
implicit operator of a delta form, it leaves the steady state where the residual vanishes
(T. H. Pulliam, *Time Accuracy and the Use of Implicit Methods*, AIAA Paper 93-3360, eq. 8): the
converged solution is that of $F + S = 0$, whatever the time step, and the factor only lets
steady runs keep the mean-flow CFL near walls. The residuals of $\tilde\nu$, $k$ and $\omega$ printed during
the run (the change of the variable over the step) are scaled by the same factor, so near walls they
are smaller than without it; the convergence check uses the density residual and is
not affected.

With the [energy coupling](turbulence-variants.md#energy-coupling) the energy of the mean flow has
the sink $-S_k$ of the $k$ equation, whose derivative with respect to $\rho k$ is $+\beta^\ast\omega$, so
the factor is that of the two equations together: the implicit operator of the $k$ and energy rows
is $D = \begin{pmatrix} 1 + \Delta t\,\beta^\ast\omega & 0 \\ -\Delta t\,\beta^\ast\omega & 1 \end{pmatrix}$, and the
energy takes its own increment plus $\Delta t\,\beta^\ast\omega$ times that of $k$,

$$
\Delta(\rho E_m) = \Delta t\,\bigl(F_E - S_k\bigr) + \Delta t\,\beta^\ast\omega\,\Delta(\rho k), \qquad
\Delta(\rho k) = \frac{\Delta t}{1 + \Delta t\,\beta^\ast\omega}\,\bigl(F_k + S_k\bigr).
$$

The two increments add up to $\Delta t\,(F_E + F_k)$, the fluxes of $E_m + k$ alone, at every stage, and
the converged solution, $F_k + S_k = 0$ and $F_E = S_k$, does not depend on the time step.

The Reynolds-stress model SSG/LRR, always in this energy form, takes one factor $1/(1 + \Delta t\,d)$ for its six
stress equations, with the rate $d$ of the table, and its own for $\omega$; its energy row is treated in the same way,
with $\rho k = \tfrac12(\rho R_{11} + \rho R_{22} + \rho R_{33})$ and $S_k$ half the trace of the stress sources.

### Time-accurate runs: explicit sources

In time-accurate runs the factor is not applied. Since $(F + S)/(1 + \Delta t\,d) = (F + S)(1 - \Delta t\,d) + O(\Delta t^2)$, the factor
would make the turbulence equations, and the mean flow through $\mu_t$, first-order accurate in
time. With the $k$–$\omega$ models the error of the factor would not even vanish when the mesh and the
step are refined together, as in an explicit run, whose global step goes as $h^2$ in the wall cells: in
the first cell $\omega \approx 6\nu/(\beta_1 y^2)$ gives $\Delta t\,d \approx 24.7$ VNN at every level (7.4 at
VNN = 0.3, a factor of 0.12), and the next two cells stay of order one.

With Spalart–Allmaras
$\Delta t\,d \approx 1.64\,f_w\,y_1^+$ at VNN = 0.3 ($y_1^+$ the height of the first cell in wall units): it
vanishes as the mesh is refined, but is of order one at the usual heights.
The sources are explicit, like the fluxes, and the local time step of each cell
([Time Integration](time-integration.md#local-time-step)) includes a bound of their stiffness:

$$
\frac{1}{\Delta t_i} = \bigl(\text{convective and diffusive rates}\bigr) + \hat d_i ,
\qquad \hat d_i \ge d_i = \max_q \frac{\partial D_q}{\partial q} ,
$$

with $\hat d$ a bound of the rates $d$ of the table that the state of the cell at the start of the step gives
alone, without gradients or ghost cells:

| Model | $\hat d$ | $\hat d / d$ |
|-------|----------|--------------|
| SST | $\max(\beta^\ast, 2\max(\beta_1, \beta_2))\,\omega$ | $\beta = F_1\beta_1 + (1 - F_1)\beta_2$, and $F_1$ takes gradients: $\beta_2/\beta_1 = 1.10$ where $F_1 = 1$ (next to walls), 1 where $F_1 = 0$ |
| Wilcox 2006 | $\max(\beta^\ast, 2\beta_0)\,\omega$ | $\beta = \beta_0 f_\beta$, and $f_\beta \le 1$ takes the velocity gradient: $1/f_\beta$, 1 in 2-D ($\chi_\omega = 0$), at most $1/0.85$ |
| SA | $2\,c_{w1}\,(1 + c_{w3}^6)^{1/6}\,\tilde\nu/y^2$ | $f_w = g\,[(1 + c_{w3}^6)/(g^6 + c_{w3}^6)]^{1/6}$ grows with $g \ge 0$ towards $(1 + c_{w3}^6)^{1/6} = 2.005$: $2.005/f_w$, 1 where $r$ reaches its cap of 10, about 2 where $r \approx 1$ (an equilibrium boundary layer, where the destruction is slow against the diffusion) |

($y$ the wall distance of the source, moved below a rough wall.) For the Reynolds-stress model SSG/LRR,
$\hat d = \max\bigl(\max(C_1^\text{LRR}, C_1^\text{SSG})\,C_\mu,\ 2\max(\beta_\omega^\text{LRR}, \beta_\omega^\text{SSG})\bigr)\,\omega
= 0.1656\,\omega$: with $\omega$ held, the dissipation and the return to isotropy of the stresses,
$-\tfrac23\rho\varepsilon\,\delta_{ij} - C_1\rho\varepsilon\,a_{ij}$ with $\varepsilon = C_\mu k\omega$, decay the trace at
$C_\mu\omega$ and the deviatoric part at $C_1 C_\mu\omega$, and $\omega$ at $2\beta_\omega\omega$; the bound leaves out the part
$\tfrac12 C_1^\ast P_{kk}/k$ of the return to isotropy, which needs the velocity gradient.

Then
$\Delta t\,d \le \Delta t\,\hat d \le 1$: a forward-Euler stage cannot make
$\tilde\nu$, $k$ or $\omega$ negative through their destruction, and neither can Heun and RK3, whose
stages are convex combinations of forward-Euler stages (W. Hundsdorfer, *Numerical Solution of
Advection–Diffusion–Reaction Equations*, lecture notes, 2000, §6.2, Theorem 6.4; C.-W. Shu,
S. Osher, *J. Comput. Phys.* 77, 1988, Proposition 2.1 of its report version, ICASE Report 87-33).
The destruction of $\tilde\nu$ and of $\omega$ is quadratic, $D = (d/2)\,q$, which leaves them a
margin of 2.

Where the destruction is faster than the transport, as in the wall cells of fine
meshes, above all with SST, the step is smaller than the convective and diffusive limits.

The bound covers the destruction only. Three sinks stay outside it: the cross-diffusion of SST,
$2(1-F_1)\,\rho\,\sigma_{\omega 2}\,\nabla k\cdot\nabla\omega/\omega$, negative where $\nabla
k\cdot\nabla\omega < 0$ (Wilcox 2006 switches its cross-diffusion off there); with the compressible
production, the dilatation part of the production of SST, $-\tfrac23\rho k\,\nabla\cdot\mathbf{u}$
and $-\tfrac23\gamma\rho\omega\,\nabla\cdot\mathbf{u}$, negative where the flow expands; and the
production of SA-R and SA-RC, which can be negative («the SA-RC production term can be negative»,
NASA Turbulence Modeling Resource). Like the transport, they can make a stage negative.

The convex combination keeps the bound of the forward-Euler stage only where that bound holds on the
state each stage starts from: Shu and Osher assume that the forward-Euler step keeps the property,
under one step restriction, for every state it is applied to (ICASE Report 87-33, (2.7) and
Proposition 2.1), while here the restriction depends on the state through $d$, and the rates are
those of the state at the start of the step.

A stage from a state where the rate is $d_s$ takes
$q_s - \Delta t\,(d_s/2)\,q_s$ through a quadratic destruction, positive while $\Delta t\,d_s < 2$; with
$\Delta t\,\hat d \le 1$ this holds as long as $d_s < 2\,\hat d$. The margin of 2 of $\tilde\nu$ and $\omega$
therefore covers a rate that grows within the step up to twice the bound at the start (at least twice the rate),
and no more; $k$, whose destruction is linear, keeps only the margin by which its rate $\beta^\ast\omega$ lies below
$\hat d$.

An
impulsive start of a $k$–$\omega$ model from a state far from the wall value of $\omega$ goes beyond
that: during the first step the wall condition raises $\omega$ in the wall cells by orders of
magnitude, and with the full step the explicit destruction overshoots ($\omega < 0$, then the floor
and a very large $\mu_t$). Start such runs with the
[CFL ramp](time-integration.md#global-time-step).

MOSE counts the cell updates in which a stage made $\tilde\nu$, $k$ or $\omega$ negative before the
realizability floor ([Running MOSE](../user/using.md#logging-and-diagnostics)). In a steady run the count covers the whole run: the start, and the fronts of $k$ that the fluxes
push below zero where $k$ is about zero, which can last to the last iteration; it is a total, which
convergence does not bring back to zero. The transport keeps the turbulence variables non-negative,
with every limiter but van Albada, when at every stage the Courant number of the flow leaving each
cell is at most 1/2 ([Time Integration](time-integration.md#local-time-step)); a front can go below
zero with a larger step (a CFL above 0.5, such as the 0.8 suggested with RK3, or the faster waves at captured shocks),
through the sources, in the first cell next to a wall, whose bound is lower, and next to faces whose
ghost cells extrapolate the cells (the velocity inlet, the solid propellant grain), and with residual smoothing of
every equation, which mixes the updates of neighbouring cells.

The residuals of the residual history are the change of
each variable over a step ([Stopping criteria](../user/using.md#stopping-criteria)): a variable
that the floor holds at its value in a cell from one step to the next adds nothing there to its
residual, so these fronts do not show in the residual history, and the count is their trace.

### Relation to under-relaxation

Point-implicit damping is a *physically weighted, automatic* form of
under-relaxation: the effective relaxation factor $1/(1+\Delta t\,d)$ is
large (≈1) where the source is mild and small only where the destruction is
stiff, and it scales with the local time step.  A constant global
under-relaxation factor on the turbulence update achieves a similar
stabilising effect but is blunter — it slows convergence everywhere, must be
hand-tuned per case/grid, and (if applied to the whole update rather than the
increment) can bias the transient.  Point-implicit treatment is therefore the
preferred technique; global under-relaxation remains a valid, simpler
fallback when a model's destruction Jacobian is not readily available.

---

## References

[13] NASA Langley Research Center, Turbulence Modeling Resource, "The Spalart–Allmaras Turbulence Model", https://turbmodels.larc.nasa.gov/spalart.html (TMBWG/turbmodels commit 20a39f5, consulted 2026-09-29): minimum distance to the wall.

[14] A. Syrakos, S. Varchanis, Y. Dimakopoulos, A. Goulas, J. Tsamopoulos, "A critical analysis of some popular methods for the discretisation of the gradient operator in finite volume methods," *Phys. Fluids*, 29, 127103, 2017, doi:10.1063/1.4997682.
