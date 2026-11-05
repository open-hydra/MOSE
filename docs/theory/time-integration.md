# Time Integration

MOSE advances the semi-discrete equations in time using an explicit
strong-stability-preserving (SSP) Runge–Kutta scheme.  This page
describes the time-stepping algorithm, the CFL and VNN stability
conditions, Strang operator splitting for chemistry coupling, implicit
residual smoothing (IRS) for convergence acceleration, the
multigrid strategy, and the multirate (local time stepping) integrator
for time-accurate runs.

---

## Explicit Time Integration Schemes

MOSE provides three explicit time integration schemes, all formulated in
the unified **Shu–Osher convex-combination** form.  The active scheme is
selected at run time via `n_RK` (1, 2, or 3), set by the key `time-scheme`
of `[MOSE-Numerics]` (`euler`, `RK2` or `RK3`):

| Scheme | Order | SSP coefficient |
|--------|:-----:|:---------:|
| Forward Euler | 1 | 1 |
| 2-stage SSP Runge–Kutta (Heun) | 2 | 1 |
| 3-stage SSP Runge–Kutta (Shu–Osher) | 3 | 1 |

`time-scheme` is a required key: a run without it stops at set-up with
`[ERROR] Required parameter not set: time-scheme`. The Default column of the
[input parameters](../user/registry.md) lists `euler` for it, a placeholder
that a required key never uses.

All three are **strong-stability preserving (SSP)**: each stage is a
convex combination of forward-Euler steps, so the positivity and TVD
properties of the spatial discretization under a forward-Euler step
$\Delta t$ hold for the same $\Delta t$ (SSP coefficient 1). These are not
the linear stability limits of the local time step: both are given in
[Stability Conditions](#stability-conditions).

---

### Forward Euler (`n_RK = 1`)

$$
\mathbf{U}^{n+1} = \mathbf{U}^n + \Delta t\;\mathcal{L}(\mathbf{U}^n)
$$

First-order accurate, one residual evaluation per step.  Useful for
initial testing or very stiff problems where a minimal time step is
already imposed by chemistry or acoustics.

---

### 2-stage SSP Runge–Kutta (`n_RK = 2`)

Also known as **Heun's method** or SSP-RK2:

$$
\begin{aligned}
\mathbf{U}^{(1)} &= \mathbf{U}^n + \Delta t\;\mathcal{L}(\mathbf{U}^n)\\[4pt]
\mathbf{U}^{n+1} &= \tfrac{1}{2}\,\mathbf{U}^n
  + \tfrac{1}{2}\bigl(\mathbf{U}^{(1)}
  + \Delta t\;\mathcal{L}(\mathbf{U}^{(1)})\bigr)
\end{aligned}
$$

Second-order accurate, two residual evaluations per step.

---

### 3-stage SSP Runge–Kutta (`n_RK = 3`)

The **Shu–Osher RK3** method:

$$
\begin{aligned}
\mathbf{U}^{(1)} &= \mathbf{U}^n + \Delta t\;\mathcal{L}(\mathbf{U}^n)\\[4pt]
\mathbf{U}^{(2)} &= \tfrac{3}{4}\,\mathbf{U}^n
  + \tfrac{1}{4}\bigl(\mathbf{U}^{(1)}
  + \Delta t\;\mathcal{L}(\mathbf{U}^{(1)})\bigr)\\[4pt]
\mathbf{U}^{n+1} &= \tfrac{1}{3}\,\mathbf{U}^n
  + \tfrac{2}{3}\bigl(\mathbf{U}^{(2)}
  + \Delta t\;\mathcal{L}(\mathbf{U}^{(2)})\bigr)
\end{aligned}
$$

Third-order accurate, three residual evaluations per step.  TVD under
the standard CFL limit when paired with a TVD spatial discretization.

---

### Unified coefficient form

All three schemes share the same kernel.  At stage $k$ the update is

$$
\mathbf{U}^{(k)} = \mathbf{U}^n
  + c_k\bigl(\mathbf{U}^{(k-1)} - \mathbf{U}^n + \mathbf{R}^\ast\bigr)
$$

where the stage coefficients $c_k$ are stored in the matrix
`RKcoeff(n_rk, stage)`:

$$
\mathbf{C} =
\begin{pmatrix}
1 & 0 & 0 \\[2pt]
1 & \tfrac{1}{2} & 0 \\[2pt]
1 & \tfrac{1}{4} & \tfrac{2}{3}
\end{pmatrix}
$$

Row index = `n_RK`; column index = stage.

### Stage times and time-dependent boundary data

Stage $k$ evaluates the residual of $\mathbf{U}^{(k-1)}$, which approximates
the solution at the stage time $t^n + \theta_k\,\Delta t$:

| Scheme | $\theta_1$ | $\theta_2$ | $\theta_3$ |
|--------|:---:|:---:|:---:|
| Forward Euler | 0 | | |
| SSP-RK2 (Heun) | 0 | 1 | |
| SSP-RK3 | 0 | 1 | ½ |

In time-accurate runs the time-dependent boundary data (inlet `402`,
assigned state `410`) are evaluated at the stage time; with Strang splitting,
inside each half step. This is the conventional imposition analysed by
Carpenter et al. [7]. On a fixed mesh, with $\Delta t \to 0$, both schemes keep
their order in time (second for Heun, third for SSP-RK3). Refining mesh and
step together at a fixed CFL number, the boundary error lowers the global
order of the third-order scheme towards two (in the fourth-order experiments
of [7], table I.a, the rate settles at 2.5 in the L2 norm); Heun stays second
order. Both statements follow from the boundary error derived below. The
experiments of [7] show the same loss with a weak (SAT) imposition of the
data, and its appendix A shows that data entering as a source term give the
same scheme; the characteristic boundary fluxes of MOSE are outside that
analysis and are checked by the tests (`TimeIntegration`/`TimeControl`).

[7] gives two ways to keep the full order. One imposes no boundary data at
the intermediate stages, only at the end of the step: the scheme keeps its
order, also for nonlinear systems, but in the experiments of [7] the CFL
number had to be lower by at least a factor of two, and often much more
(text after table I.a). The other imposes intermediate values built from the
data and its time derivatives (equations (22)–(27) of [7]) at the same time
step; it is general for linear problems only, and for nonlinear systems [7]
shows it for the third-order Runge–Kutta scheme of Heun (stage times 0, 1/3,
2/3; section 3), not for SSP-RK3. MOSE uses neither: tabulated data have no
time derivatives, and the step would shrink.

#### Boundary error of the stage-time data

??? note "Derivation for Heun and SSP-RK3"
    The derivation of [7], section 2, carried out for the two schemes of MOSE
    ([7] does it for the classical fourth-order scheme and notes that it can be
    done for any Runge–Kutta scheme). Model problem ([7], equations (1)–(4)):
    $u_t + u_x = 0$ for $x > 0$, $u(0,t) = g(t)$, semi-discretised as
    $\mathrm{d}V/\mathrm{d}t = (1/h)\,DV$ for the interior entries of $V$, with the
    boundary entry $v_0 = g$. Split $D$ by its first row and column,

    $$
    D = \begin{pmatrix} d_{00} & r^{\mathsf T} \\ c & B \end{pmatrix},
    $$

    so that $c$ couples the interior entries to the boundary one, and let
    $\lambda = \Delta t/h$. Each stage computes every entry with $D$ and then sets
    the boundary entry to the data at its stage time; that change is $G^k e_0$, as
    in [7], equations (9)–(16). With the exact solution in $V^n$,
    $\lambda^k (D^k V^n)_0 = \Delta t^k g^{(k)}(t^n) + \dots$ ([7], equation (18)),
    and the scheme without the boundary changes is accurate to its own order.

    *Heun* (stage times 0, 1):
    $V^{(1)} = V^n + \lambda D V^n + G^0 e_0$,
    $V^{n+1} = \tfrac12 V^n + \tfrac12\bigl(V^{(1)} + \lambda D V^{(1)}\bigr) + G^1 e_0$, with
    $G^0 = g(t^n + \Delta t) - g - \Delta t\,g' = \tfrac12 \Delta t^2 g'' + O(\Delta t^3)$.
    The interior entries of $V^{n+1}$ differ from those of the scheme without the
    boundary changes by

    $$
    \tfrac{\lambda}{2}\,G^0\,c = \tfrac{\lambda}{4}\,\Delta t^2\,g''\,c + O(\Delta t^3) .
    $$

    *SSP-RK3* (stage times 0, 1, ½):
    $V^{(1)} = V^n + \lambda D V^n + G^0 e_0$,
    $V^{(2)} = \tfrac34 V^n + \tfrac14\bigl(V^{(1)} + \lambda D V^{(1)}\bigr) + G^1 e_0$,
    $V^{n+1} = \tfrac13 V^n + \tfrac23\bigl(V^{(2)} + \lambda D V^{(2)}\bigr) + G^2 e_0$, with
    $G^0$ as above and, from the boundary entry of $V^{(2)}$ before the change,

    $$
    G^1 = g(t^n + \tfrac12\Delta t) - g - \tfrac12\Delta t\,g' - \tfrac14\Delta t^2 g'' - \tfrac14 G^0 - \tfrac{\lambda}{4}\,d_{00}\,G^0 + O(\Delta t^3)
    = -\tfrac14\Delta t^2 g'' - \tfrac{\lambda}{4}\,d_{00}\,G^0 + O(\Delta t^3) .
    $$

    The interior entries of $V^{n+1}$ differ by

    $$
    \bigl(\tfrac{\lambda}{3}G^0 + \tfrac{2\lambda}{3}G^1\bigr)\,c + \tfrac{\lambda^2}{6}\,G^0\,(d_{00}\,c + Bc) :
    $$

    the terms in $\lambda\,\Delta t^2$ cancel ($\tfrac16 - \tfrac16$), and so do those in
    $d_{00}\,c$, which leaves

    $$
    \tfrac{\lambda^2}{12}\,\Delta t^2\,g''\,Bc + O(\Delta t^3) .
    $$

Both errors sit next to the boundary, where $c$ and $Bc$ are not zero. At a
fixed CFL number they are $O(\Delta t^2)$ per step, first order there per unit
of time, as equations (20)–(21) of [7] for the fourth-order scheme, and the
global order is two ([7], abstract; it cannot exceed two, after equation
(21)): SSP-RK3 drops from three to two, and two is the order of Heun. On a
fixed mesh they are $\Delta t^3 g''c/(4h)$ and $\Delta t^4 g''Bc/(12h^2)$ per step,
$O(\Delta t^2)$ and $O(\Delta t^3)$ per unit of time: the orders of the schemes. The
derivation is linear and one-dimensional, like that of [7].

### State update

At each RK stage $k$ the full update sequence is:

1. Convert primitives to conservative: $\mathbf{U} = \text{prim2cons}(\mathbf{P})$
2. Scale residual: $\mathbf{R}^\ast = -\mathbf{R}\,\Delta t / V$
3. *(Optional)* Apply implicit residual smoothing to $\mathbf{R}^\ast$
4. RK combination using the coefficient above
5. Recover primitives: $\mathbf{P}^{(k)} = \text{cons2prim}(\mathbf{U}^{(k)})$
6. Enforce realizability (positive densities and pressures)

An alternative **primitive-variable update** path is available, where
the residual is scaled and applied directly in primitive space.  This
path is faster per step but less conservative.

### Reference velocity of the preconditioning

With `integration-variables = prec` the update of a stage, the local time step (below) and the flux `HLLC-PC` take in
every cell the reference velocity $U_r$ of the Weiss–Smith preconditioning. It is computed after the update of each
stage from the new state of the cell, and the ghost cells of the faces between blocks receive it with the state at the
start of the next stage. By default $U_r = \varepsilon a$ for $|\mathbf v| < \varepsilon a$, $a$ for
$|\mathbf v| \ge a$ and $f\,|\mathbf v|$ otherwise ($\varepsilon$ = `preconditioning-eps-min`, $f$ =
`preconditioning-Ur-factor`), then raised to the viscous floor

$$
U_{\text{visc}} = \max\Bigl( \frac{\mu + \mu_t}{\rho\,\Delta_{\min}},\;
\frac{k + \mu_t\,c_p/\text{Pr}_t}{\rho\,c_p\,\Delta_{\min}} \Bigr),
$$

with $\Delta_{\min}$ the smallest width of the cell, and to `preconditioning-Ur-min`
(`preconditioning-Uref` and `preconditioning-Mach` select the other modes). With a RANS model $\mu_t$ is the eddy
viscosity of the model with the velocity gradient of the cell: SST takes $\mu_t = \rho a_1 k/\max(a_1\omega, S F_2)$
and Wilcox 2006 $\mu_t = \rho k/\hat\omega$, $\hat\omega = \max\bigl(\omega, C_{\lim}\sqrt{2\bar S_{ij}\bar
S_{ij}/\beta^*}\bigr)$, so the strain rate lowers it; Spalart–Allmaras does not use the gradient. The gradient is
that of the state the stage starts from, computed after the ghost cells of the stage are filled, with the cell metric of
the source terms, as the turbulence sources take it; the setup computes it from the initial state. After the update the
ghost cells still hold the start of the stage: a gradient computed there would take, in a cell next to a face between
blocks, the cell across the face at the start of the stage and the cells of its own block after the update, and $U_r$
would depend on how the mesh is divided into blocks.

---

## Stability Conditions

### Local time step

Each cell $i$ gets a local time step that combines convection and diffusion in one additive rule
([8], §6 and §9.3; [5], eq. 26–28):

$$
\frac{1}{\Delta t_i} =
\frac{1}{\text{CFL}} \sum_{d=1}^{n_d} \Bigl( \bigl|\mathbf{v}\cdot\nabla\xi_d\bigr| + a\,\bigl|\nabla\xi_d\bigr| \Bigr)
\;+\; \frac{D_{\max}}{\text{VNN}} \sum_{d=1}^{n_d} \bigl|\nabla\xi_d\bigr|^2 ,
$$

where $\nabla\xi_d$ is the row $d$ of the cell metric (its modulus is the inverse of the width of the
cell across its faces of direction $d$), $n_d$ the number of computational directions (one cell in
$k$ gives $n_d = 2$), $a$ the speed of sound and $D_{\max}$ the largest diffusivity of the equations:

$$
D_{\max} = \max\!\left(
\frac{4}{3}\,\frac{\mu_\ell + \mu_t}{\rho},\;
\frac{k_\ell + \mu_t\,c_p/\text{Pr}_t}{\rho\,c_v},\;
D_s,\;
\frac{\mu_\ell + \rho\tilde\nu}{\sigma\,\rho}
\right).
$$

The conductivity acts on the internal energy $\rho c_v T$, so the energy diffusivity is
$\gamma\nu/\text{Pr} \approx 1.94\,\nu$ in air, larger than the momentum one; $D_s$ are the species
diffusivities (more than one species only); the last term enters with Spalart–Allmaras
($\sigma = 2/3$). The $k$ and $\omega$ diffusivities of SST and Wilcox 2006,
$(\mu_\ell + \sigma\mu_t)/\rho$ with $\sigma \le 1$, stay below the momentum one. Euler runs have
no diffusive term.

In 3-D, on a cell whose three directions are skewed against each other, the viscous rate is
multiplied by $F = (2 + m + 1/m)/4$: with $a$, $b$, $c$ the cosines of the angles between the three
$\nabla\xi_d$, where $a^2 + b^2 + c^2 + 2abc > 1$, $m + 1$ is the largest eigenvalue of the
$3\times3$ matrix of the cosines ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$).
There the spectral radius of the viscous operator exceeds $4D_{\max}\sum_d|\nabla\xi_d|^2$ by a
factor below $9/8$, which $F$ bounds (exactly on rhombohedra); elsewhere, orthogonal cells and cells
with one direction orthogonal to the other two included, $F = 1$, and in 2-D the rate holds as it
is.

Without dimensional splitting the convective limits of the directions add up: an explicit scheme
on square 2-D cells needs a condition of the type $\tau(|a_x|/\Delta x + |a_y|/\Delta y) \le C_0$, not
the largest of the two terms ([8], §9.3, which states the form without a proof). For forward
Euler the form follows, with the diffusive rates added as well, from the positivity condition (6.1)
of [8]: with first-order upwind convection and central diffusion on a Cartesian mesh the
diagonal entry of a row is $-\sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr)$ and the others are not negative, so
$\tau \sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr) \le 1$: the rule above with CFL = 1 and VNN = 1/2. The
constants of Heun and RK3 in the table do not come from that source: they are the von Neumann limits
of the linearised MUSCL–upwind advection with central diffusion, proven for MOSE's stencil in 1-D and
carried to 2-D and 3-D, for any split of the Courant number, by the convexity of the symbol; they hold
for the van Leer slope (the default), van Albada and MC, not for minmod and superbee. The source warns
that such an eigenvalue criterion can mislead for the non-normal operators that boundaries give, as an
inflow boundary with upwind differences ([8], §3.3, pp. 26–27): the limits are those of the interior.

??? note "Proof of the limits"
    The
    van Leer slope (the default; van Albada and MC linearise the same way) linearised in smooth flow is
    the mean of the two one-sided slopes (Fromm's scheme), with the symbol
    $g(\theta) = \tfrac12[(1-\cos\theta)^2 + i\sin\theta(3-\cos\theta)]$, $z = -\nu g$; the constants
    hold for these three limiters, not for minmod and superbee, which take the smaller or the larger of
    the two slopes and have no linearisation where they are equal. Heun is stable if and only if
    $\nu \le 1$: $|R(z)|^2 - 1 = \tfrac{\nu}{16}(1-c)^2[\nu(\nu a - 2b)^2 + 4\nu b^2 - 16]$ with
    $c = \cos\theta$, $a = 5 - c^2$, $b = 1 - c$, increasing in $\nu$ and equal to
    $\tfrac{\nu}{16}(1-c)^4(1+c)(c-3)$ at $\nu = 1$. RK3 is stable if and only if
    $\nu \le \nu^\ast = 1.1757678815978$, a root of a polynomial of degree 24 (a cylindrical
    decomposition in exact arithmetic), and with the first-order flux if and only if
    $\nu \le 1.2563726633$, half the real stability interval $2.5127453$; diffusion gives VNN $\le 1/2$
    and $0.6282$ (the table rounds them down). The curve of $-g$ bounds a convex set $K$
    ($x'y'' - y'x'' = -\tfrac12(1-\cos\theta)(\cos^2\theta + \cos\theta + 4) \le 0$ along it), so a sum
    over directions is a point of $\nu K$, inside the stability region by the maximum principle: the
    limits of 1-D hold in 2-D and 3-D for any split of the Courant number, with no need for a convex
    stability region. With diffusion the additive rule is stable for Heun with CFL = 1 and VNN = 1/2,
    and exact at every Péclet number (the checkerboard mode reaches $z = -2$); for RK3 the region is not
    convex, and the convex hull of $1.17K$ and $[-2.48, 0]$ (and of $1.1757K$ and $[-2.512, 0]$) is
    shown inside it with a certificate in exact rational arithmetic.

| Scheme | real stability interval | CFL ≤, MUSCL, smooth flow | CFL ≤, first order (limiter at extrema and shocks) | VNN ≤ | CFL and VNN ≤ for positivity with the limiter (equal rectangular cells, smooth flow) |
|---|---|---|---|---|---|
| Heun (`RK2`) | $[-2, 0]$ | 1 | 1 | 1/2 | 1/2 |
| `RK3` | $[-2.5127, 0]$ | 1.17 | 1.25 | 0.62 | 1/2 |
| explicit Euler (`euler`) | $[-2, 0]$ | unstable | 1 | 1/2 | 1/2 |

On the 2-D axisymmetric wedge (type `200` or `300` on the $k$ faces) the stress on the wedge faces adds to the
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
axis, for the scalars; for the momentum next to the axis they rest on evidence (test `WedgeAxisStep`).
The positivity limits of its last column hold on a mesh of equal rectangular cells in smooth flow, the one of the
diffusion with a constant diffusivity (below).

??? note "What supports the limits of the table"
    - Away from the axis, on a uniform mesh of equal cells (parallelograms in 2-D, parallelepipeds in 3-D) with frozen
      coefficients and no boundaries, the viscous operator of MOSE has a real symmetric Fourier symbol. With the face gradient of
      MOSE (the difference across the face, the mean of the central differences of the two cells along it, mapped with the face
      metric, which on such a mesh is the cell metric for every `diffusive-metric`), a mode $e^{i\boldsymbol\theta\cdot\mathbf c}$
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
      energy in air). In 3-D the sharp constant is $9/2$: with $a$, $b$, $c$ the cosines of the angles between the three
      $\nabla\xi_d$, the rate bounds the spectral radius if and only if $a^2 + b^2 + c^2 + 2abc \le 1$ (true for orthogonal
      cells, for cells with one direction orthogonal to the other two, and for a rhombohedron up to the edge angle
      $\arccos(-1/3) = 109.47^\circ$); elsewhere the radius exceeds it by a factor below $9/8$, and the local step multiplies the
      rate by $F = (2 + m + 1/m)/4$, where $m + 1 \in (2, 3)$ is the largest eigenvalue of the $3\times3$ matrix of the cosines
      ($m$ the largest root of $m^3 - (a^2+b^2+c^2)\,m - 2abc = 0$): a proven bound of the factor, exact on rhombohedra
      (1.04 at $115^\circ$, 1.11 at $119.5^\circ$), below $9/8$ on every cell, and 1 where the rate already holds. The real
      stability intervals of Heun, $[-2, 0]$, and of RK3, $[-2.5127, 0]$, then give VNN = 1/2 and 0.62 ($2.5127/4 = 0.628$).
    - Next to the axis the operator is not translation invariant. For the scalars (energy, species, turbulence variables) it stays
      self-adjoint in the inner product weighted by the cell volumes: the ghost cells of the wedge copy the scalars and the wedge
      faces do not update them, so on a mesh of equal rectangular cells in $(x, r)$ each face couples two cells symmetrically, and
      the row sums are at most $4D\sum_d|\nabla\xi_d|^2$ (up to a relative $\varepsilon/6$ for an axis at $\varepsilon\,\Delta r$,
      covered by $s_h$ for diffusivities below $0.6\,\nu_e/\varepsilon$). A self-adjoint operator whose row sums are within the
      rate of each cell is stable with the local step, in the norm weighted by $V/\Delta t$, and the limits of the table hold for
      the scalars on every row. The momentum operator is not self-adjoint in any inner product with cell weights: the cross terms
      $(\lambda+\mu)\,\partial^2/\partial x\,\partial r$ of the faces in $x$ and in $r$ are weighted by the radius of the cell,
      $(j-\tfrac12)\Delta r$, and of the face, $j\Delta r$. In the first row the rows of $u_r$ and $u_\theta$ are diagonally
      dominant (their diagonals carry the wedge terms $-\tfrac43\nu_e/r^2$ and $-\nu_e/r^2$); the row of $u_x$ is not: next to the
      axis (doubled weight of the upper face, mirror image across the axis) it couples to $u_r$ of four neighbours with
      $\nu/(3\Delta x\,\Delta r)$ each, four times the cross coupling of the plane operator, a deficit of
      $\tfrac43\nu/(\Delta x\,\Delta r)$; its row sum stays below the rate of the cell. No proof covers the momentum next to the
      axis. The support is evidence: on the operator of the code for cells with $dx/dr$ from 0.01 to 128 and wedges of 0.25° to 5°,
      its symmetric part is negative semidefinite, its eigenvalues are real to $6\cdot10^{-4}$, and the largest stable VNN with the
      local step is 0.500 to 0.577 with Heun and 0.628 to 0.725 with RK3, approaching the limits of the table, 1/2 and 0.62
      (0.628), from above for very long or very short cells; on the diffusion-dominated 1° wedge of the test `WedgeAxisStep`
      (32 × 32 cells, 2000 local steps) the largest stable VNN is about 0.5 with Heun (0.508 to 0.519)
      and 0.63 to 0.65 with RK3 for $dx/dr$ = 1, 4 and 8 at Prandtl numbers 0.72 and 0.9 (and 10 at $dx/dr$ = 4), and the
      test runs 2000 steps at the limits and at the recommended values. The eigenvalue criterion can mislead for non-normal operators
      ([8], pp. 26–27).

Next to the axis of square cells $s_h = 8\,\nu/\Delta r^2$; far from it $s_h$ vanishes as $1/r$. The wedge stress, and with it $s_h$, applies
with either sign of the wedge angle (a right-handed mesh with the axis on the last $j$ face has a negative one). Test
`WedgeAxisStep`. The step takes $s_h$ on every 2-D mesh with a wedge angle, whatever the type of its $k$ faces: the viscous
flux of faces of type `200` and that of symmetry faces (`300`), whose mirror image of a cell is, without swirl, the cell
turned by the wedge angle, both carry this stress; null (`0`) $k$ faces carry none, and there $s_h$ only makes the step
smaller than it needs to be.

These limits are linear (von Neumann, smooth solutions). With the limiter, positivity (a maximum principle:
in more than one dimension a second-order scheme is not TVD, Goodman et al. 1985 as cited in [9], p. 15 of the
report) is guaranteed only for smaller steps, for the
advection of a scalar with the limited flux: on a line of equal cells the van Leer slope of MOSE is the flux limiter
$\psi(\theta) = \theta/(1+\theta)$, and with the bounds of the reconstruction
([Spatial Discretization](numerics.md#muscl-reconstruction), step 3) the limited flux of every limiter but van Albada
satisfies $0 \le \psi \le 1$ and $0 \le \psi/\theta \le \mu$ with $\mu = 1$ on a line of cells of any lengths. With
Heun and RK3 (SSP coefficient 1) positivity then holds for $(1+\mu)\,\nu \le 1$: on any mesh, when at every stage the
Courant number of the flow leaving each cell (the step times the volume flow out of the cell through its faces, over
its volume) is at most 1/2. Central diffusion of a scalar is positive for $\text{VNN} \le 1/2$ on a mesh of equal
rectangular cells with a constant diffusivity, where the face gradient has no cross metric terms (the condition (6.1) of
[8] above). On other rectangular cells the bound can be lower for a cell that advances with the step of the rule (a
local step, or a level of the multirate; with the global step of a time-accurate run, the smallest of the mesh,
$\text{VNN} \le 1/2$ keeps every cell positive in one dimension with a constant diffusivity): the rule takes the width
and the largest diffusivity of the cell itself (`Mod_dt`), while the flux of a face takes the distance between the two
centres (with
`diffusive-metric = stencil`, the default) and the diffusivity of the mean state of the two cells (`Lib_Diffusive`),
which can exceed that of the cell; in one dimension, with a constant diffusivity and no convection, the positivity of
every stage asks $\text{VNN}\,h\,(1/\delta_- + 1/\delta_+) \le 1$ (the forward Euler step; Heun and RK3 inherit it as
convex combinations of such steps), with $h$ the width of the cell and $\delta_\mp$ the distances to
the centres of its two neighbours: $\text{VNN} \le 1/2$ on equal cells and on widths in a geometric progression, and
towards 1/4 in a cell much larger than both neighbours. On skewed cells the face gradient also takes the central
differences of the two cells along the face (`Lib_Diffusive`), whose cross metric terms can give corner neighbours of a
cell a negative weight: on a mesh of equal parallelograms two corner neighbours of every cell have one, at every skew
angle, and no step makes the diffusion positive there. With the additive rule
both hold for CFL ≤ 1/2 and VNN ≤ 1/2 in smooth flow on a mesh of equal rectangular cells with a constant diffusivity,
where the flow leaving a cell through its faces is about the
convective rate of the rule (at captured shocks see below) ([8], §5.2, §6.2 and §6.4; experiments reported there give
positive results up to $\nu = 1$ in 1D and about $2/3$ in 2D). The bound is lower in the first cell next to a wall of a
RANS run, whose ghost cell holds $-k_1\rho_G$ (and $-\tilde\nu_1\rho_G$; the wall value of $k$, and of $\tilde\nu$ on
a smooth wall, is zero), $\rho_G$ the density of the ghost cell: there it is $1/(2 + \rho_G/\rho_1)$, and in smooth
flow CFL ≤ 1/2 keeps it while the Mach number of the flow normal to the wall is below $2\rho_1/\rho_G$. Next to faces whose ghost
cells extrapolate the cells, $3q_1 - 3q_2 + q_3$ (type `408`, and the flow and soot variables of type `502` in its mode 1),
positivity is not guaranteed. With moving shocks and Heun keep `cfl` at
0.5–0.7. With preconditioning,
$|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$ is replaced by the preconditioned eigenvalue
$\bigl(u_n (1-\alpha) + c^\ast\bigr)|\nabla\xi_d|$.

The local time step takes the wave speeds of the cell alone. With a stationary shock captured on one face, the fastest
wave that the HLLC flux lets into a cell on the slow side of the shock exceeds the cell's own $|u_n| + a$ by up to 1.47
(a stretched cell crossed obliquely by a strong shock; 1.34 for a normal shock aligned with the cells, 1.25 for the sum
over the two directions of a square cell): with CFL = 0.8 the face Courant number stays below 1.256, the first-order
limit of RK3 that applies where the limiter has removed the slope. With strong shocks keep CFL ≤ 0.8 for RK3. The local
step changes the path of a steady run, not the steady solution it converges to.

??? note "Courant numbers of the HLLC flux at a captured shock"
    The local time step takes the wave speeds of the cell alone, $|\mathbf{v}\cdot\nabla\xi_d| + a\,|\nabla\xi_d|$. The HLLC
    flux lets waves enter a cell through a face with the speeds of Batten et al. (SIAM J. Sci. Comput. 18, 1553–1570, 1997, doi:10.1137/S1064827593260140, eq. 51), which
    include the Roe-averaged state of the face: on the slow side of a shock the fastest wave entering a cell is faster than the
    cell's own $|u_n| + a$, so the Courant number of that face exceeds the CFL number. For a gas with $\gamma = 1.4$ and a
    stationary shock captured on one face, the ratio is 1.27 at Mach 3 and at most 1.34 for a normal shock aligned with the cells,
    at most 1.25 for the sum over the two directions of a square cell, the quantity the time step rule uses, and up to 1.47 in one
    direction of a stretched cell whose face is inclined to a strong shock (1.32 at Mach 3, 1.45 at Mach 10). With CFL = 0.8 the
    face Courant number is thus at most 1.08, 1.00 and 1.18 respectively: below the linear stability limit of RK3 with the
    first-order upwind flux, 1.256, the limit that applies where the limiter has removed the slope (the limit of the limited scheme
    itself is not derived), and below the limit with the smooth MUSCL reconstruction, 1.176, except on stretched cells crossed
    obliquely by shocks stronger than about Mach 50 (1.177 in the limit of a strong shock). Lower $\gamma$ raises the ratios (1.49
    instead of 1.34 for $\gamma = 1.2$). With strong shocks keep CFL ≤ 0.8 for RK3. The positivity of the first-order Euler step
    with the HLLC flux, guaranteed by CFL ≤ 1/2 in smooth flow, also holds at shocks on square cells (the sum of the incoming waves
    is at most 1.90 times the cell's own rate), but not on stretched cells crossed obliquely by a strong shock, which need
    CFL ≤ 0.41 at Mach 10; with MUSCL it is observed, not proven.

**Recommended values.** For time-accurate runs with moving shocks and RANS models: with RK3, `cfl = 0.8` (the face Courant
numbers of the HLLC flux above stay within the first-order limit of the scheme, which applies at the shock: 1.47 × 0.8 =
1.18 < 1.256) and `vnn = 0.5`, 20 % below the limit 0.628 of the scheme (up to that limit, for example at `vnn = 0.6`,
4.5 % below it, the linearised scheme is stable); with Heun, `cfl` from 0.5 to 0.7 and `vnn = 0.4`, below its limit 0.5. These are
values of the scheme, not of a case: a case that runs stable beyond them owes it to its mesh (a non-uniform spacing around the
cell that sets the step lowers the spectral radius of its diffusion), and its threshold is not a limit to carry elsewhere. The
viscous rate takes the largest diffusivity, that of the energy, and adds up the directions: a `vnn` tuned on a rule that took
the momentum diffusivity and the largest direction alone gives here a smaller step (by $\gamma/\text{Pr}_t$, 1.65 with the
default $\text{Pr}_t = 0.85$, in a
stretched cell where the turbulent viscosity sets it), and can be raised within the limits of the table. The CFL number of a
run with local time steps, as a steady run, is a parameter of the pseudo-time, not a Courant number of the physical time: it is
not a value to carry to a time-accurate run.

### Global time step

In a time-accurate run all cells advance with the smallest local step:

$$
\Delta t = \min_i \Delta t_i .
$$

In time-accurate RANS runs the turbulence sources are explicit, and each local step also includes
a bound $\hat d_i$ of their stiffness, the derivative $d_i$ of the destruction terms, from the state of the cell at
the start of the step alone: $1/\Delta t_i$ gains the term $\hat d_i$, so that
$\Delta t\,d_i \le \Delta t\,\hat d_i \le 1$ and no stage makes $\tilde\nu$, $k$ or $\omega$ negative through their
destruction
([Turbulence](turbulence.md#time-accurate-runs-explicit-sources)). The Reynolds-stress model SSG/LRR, which the key
`turbulence` does not select, has no such bound: its $\hat d_i$ is zero.

!!! tip "CFL ramp-up"

    `cfl-rise-threshold` (`[MOSE-Numerics]`, 0 by default: no ramp) ramps the time step
    linearly over the first iterations: at iteration $n$ below $N$ = `cfl-rise-threshold`
    the local step of every cell is multiplied by $n/N$, with $n$ the iteration number of
    the grid level (1 at its first iteration). This improves robustness when the initial
    condition is far from the solution, and a time-accurate start of a $k$–$\omega$ model
    far from the wall value of $\omega$ needs it
    ([Turbulence](turbulence.md#time-accurate-runs-explicit-sources)).

---

## Strang Operator Splitting

For reacting flows the chemistry source term is integrated separately
from the fluid transport using **Strang splitting**, which preserves
second-order temporal accuracy:

$$
\mathbf{U}^{n+1} =
\mathcal{L}_\text{RK}\!\bigl(\tfrac{1}{2}\Delta t\bigr)
\;\circ\;
\mathcal{L}_\text{chem}(\Delta t)
\;\circ\;
\mathcal{L}_\text{RK}\!\bigl(\tfrac{1}{2}\Delta t\bigr)
\;\mathbf{U}^n
$$

1. Advance half a time step with the RK spatial operator
   (convection + diffusion)
2. Integrate the full chemical ODE system over $\Delta t$ in each cell
   (species densities + temperature)
3. Advance the remaining half step with the RK operator

The chemistry ODE is solved by the stiff integrators in the **OSLO** library.

---

## Implicit Residual Smoothing (IRS)

IRS increases the effective CFL limit of the explicit RK scheme by
smoothing the residual with a Laplacian operator before the state
update.

### Smoothing equation

In each coordinate direction $d$ the smoothed residual $\mathbf{R}^\ast$
satisfies, with $\varepsilon$ = `irs-beta` of `[MOSE-Numerics]`,

$$
\mathbf{R}^\ast_i =
\frac{\mathbf{R}_i
  + \varepsilon\,(\mathbf{R}^\ast_{i-1} + \mathbf{R}^\ast_{i+1})}
{1 + 2\varepsilon}
$$

This implicit tridiagonal system is approximated with **2 Jacobi
sweeps**, which is sufficient for practical convergence.

### Acceleration effect

The smoothing amplifies the stable CFL limit by approximately
$1 + 2\varepsilon$. The cases of the repository that use IRS take
$\varepsilon = 0.5$ (0.3 in the frozen rocket nozzle).

### Steady runs only

The smoothed residual is not the residual of the space discretization, so a
run with IRS does not follow the transient of the discretized equations: IRS
accelerates the convergence of steady (pseudo-time) runs. It is switched on
by `irs = true` or by `irs-beta > 0` alone, and MOSE refuses it at set-up
when `time-accurate = true`. In unsteady computations residual smoothing
belongs to the pseudo-time iterations of dual time stepping [5], which MOSE
does not implement.

The change does not vanish as $\Delta t \to 0$: the smoothing turns the
semi-discrete system into $(I - \varepsilon\,\delta^2)\,d\mathbf{U}/dt =
-\mathbf{R}/V$ in each direction, which slows a Fourier mode of wavenumber
$k$ by the factor $1 + 4\varepsilon\sin^2(k\Delta x/2)$: a change of order
$\varepsilon (k\Delta x)^2$ for well-resolved waves and of order one at
discontinuities. Residual smoothing has also been used in time-accurate
computations with a constant time step: Jorgenson and Chima raised the step
limit of a four-stage Runge–Kutta scheme by a factor of six, with negligible
differences in the unsteady results of periodic rotor–stator flows [6]; MOSE
refuses it because its time-accurate targets include moving shocks.

### Implementation

1. Extrapolate residuals to ghost cells (boundary conditions)
2. For each of the three coordinate directions:
    - Perform 2 Jacobi iterations of the 3-point stencil
    - Pass the smoothed residual to the next direction
3. Use the smoothed residual in the RK stage

By default (`irs-variables = flow`) only the flow residuals (species,
momentum, energy) are smoothed; the turbulence equations then advance
with the same local time step without smoothing, and above the explicit
limit they become the unstable part. `irs-variables = all` smooths every
residual (turbulence, soot and passive scalars too) in the same way. It is
meant for Spalart–Allmaras: the k-omega models diverge with it (the SST
and Wilcox 2006 flat plates of the tests at IRS 0.5).
Smoothing changes the path to the steady state, not the steady state
itself (a zero residual stays zero).

IRS is on with `irs = true` or `irs-beta > 0`, at the start and at every
reading of `input.ini` during the run. Its arrays are allocated at the start
only: a reading that switches it on in a run started without it is an error of
the reading, which `ini-error` decides ([Input](../user/input.md)). On a run
started with IRS, `irs-beta` and `irs-variables` can change during the run, and
IRS goes off when the values in force give neither `irs = true` nor
`irs-beta > 0`: `irs-beta = 0` switches it off unless `irs = true` is in force,
and `irs = false` leaves it on while `irs-beta > 0`.

---

## Multigrid Acceleration

MOSE supports **geometric multigrid** with a $2\!:\!1$ coarsening ratio
in each spatial direction, to start the fine-grid iterations from a
coarse-grid solution.

### Grid hierarchy

- **Number of levels**: the key `levels` of `[MOSE-Multigrid]` (`MGL` below).
- **Coarsening rule**: each coarse cell corresponds to $2^d$ fine cells
  ($d$ = number of active spatial dimensions).
- **Compatibility**: the fine-grid dimensions must be divisible by
  $2^{\text{MGL}-1}$ in each active direction.

!!! warning "Dimensionality and axis convention"
    Multigrid works for **3-D, 2-D and 1-D** grids, but the reduced-dimension
    cases must use a fixed axis layout, because only the trailing index
    directions are treated as collapsible (a dimension equal to $1$ is held
    fixed and not coarsened):

    | Case | Active directions | Must be singleton ($=1$) |
    |------|-------------------|--------------------------|
    | 3-D  | $i,\,j,\,k$        | —                        |
    | 2-D  | $i,\,j$ (the $i\!-\!j$ plane) | $k$            |
    | 1-D  | $i$               | $j$ **and** $k$          |

    In other words, a **1-D** case must be laid out **along $i$** ($n_j = n_k = 1$),
    and a **2-D** case must lie in the **$i\!-\!j$ plane** ($n_k = 1$). A mesh
    whose non-trivial extent is on $j$ or $k$ alone is *not* supported and will
    coarsen incorrectly.

!!! warning "Starting-mesh divisibility requirement"
    The supplied (finest) mesh must already satisfy the coarsening requirement:
    **every active direction must be divisible by $2^{\text{MGL}-1}$**, while
    singleton directions ($=1$) are exempt. For example, with `MGL = 5`
    ($2^{4} = 16$) a 1-D mesh needs $n_i$ a multiple of $16$ (e.g. $n_i = 4000$,
    $n_j = n_k = 1$). If this is not met, grid generation / setup will reject the
    mesh. Choose $n_i$ (and $n_j$ in 2-D) accordingly **before** generating the
    grid — the hierarchy is built from this mesh, not the other way around.

### Level schedule

The levels run in sequence, coarsest first: level $n$ runs `level<n>-iter`
iterations (or until `res-threshold` / `time-threshold`, which act on every
level), then its solution is prolongated onto level $n-1$. `level<n>-iter`
(and `iter-threshold` without multigrid) is re-read with `input.ini` every
`ini-diter` iterations: lowering it below the iterations already done ends
the current level at the iteration after that re-read.

Every field write of a coarse level also prolongates its solution to the
finest level and writes it (`field-prolongated`, `wall-prolongated`). On a
large mesh that write costs more than the coarse iterations between two
writes; `prolongated-diter`, a key taken at every reading of `input.ini` (column
During the run of the [Parameter Registry](../user/registry.md)), spaces it to at
least that many iterations of the level. The end of a level always writes it, so the one
written at the end of level 2 is the exact start of level 1.

### Restriction (fine → coarse)

Conservative volume averaging:

$$
\mathbf{U}_\text{coarse} =
\frac{1}{V_\text{coarse}}
\sum_{\text{fine}\,\in\,\text{coarse}}
\mathbf{U}_\text{fine}\;V_\text{fine}
$$

- 3-D: 8-cell average
- 2-D: 4-cell average (in the $i\!-\!j$ plane)
- 1-D: 2-cell average (along $i$)

Primitive variables on the coarse grid are recovered from the averaged
conservative state via Newton–Raphson (`cons2prim`), using the
fine-grid temperature as the initial guess.

### Prolongation (coarse → fine)

The coarse-grid solution is transferred to the next finer grid with
**cubic interpolation** (3-D), **biquadratic interpolation** (2-D, in the
$i\!-\!j$ plane), or **linear interpolation** (1-D, along $i$).

3-D interpolation weights:

$$
a_1 = \tfrac{27}{64},\quad
a_2 = \tfrac{9}{64},\quad
a_3 = \tfrac{3}{64},\quad
a_4 = \tfrac{1}{64}
$$

2-D interpolation weights:

$$
a_1 = \tfrac{9}{16},\quad
a_2 = \tfrac{3}{16},\quad
a_3 = \tfrac{1}{16}
$$

Each fine cell receives a weighted contribution from its parent coarse
cell and the surrounding neighbours, with the weight determined by the
fine cell's position within the coarse cell.

### Level sequence

MOSE runs the levels once each, from the coarsest to the finest (grid
sequencing). At set-up the initial solution is restricted to the coarse
levels; each coarse level $N$ runs `levelN-iter` iterations (for example
`level2-iter`), and its solution is prolongated to the next finer level,
down to level 1, which runs `level1-iter` iterations: with more than one
level, `iter-threshold` does not set the iterations of level 1. The residual
criterion of a steady run and `time-threshold` also end a level. The
iteration count starts again from zero on each level, and each level advances
its own time from the initial time of the run (zero in a new run). There are
no V or W cycles, no coarse-grid correction and no dual time stepping: the
coarse levels only provide the starting solution of the finer ones.

---

## Multirate Time Stepping

In a time-accurate run with the global step every cell advances with the
smallest admissible step of the mesh, $\Delta t = \min_i \Delta t_i$.  When a
few cells (near walls, in a throat, in a refined band) set that minimum, most
of the work is spent advancing the other cells with steps far below their own
limit.  The **multirate** integrator (`[MOSE-Multirate] enabled = .true.`)
lets each cell advance with a step close to its own limit, keeping, on a fixed
mesh, the second order in time, and the exact conservation of the single-rate
scheme.

### Levels and macro step

Each cell $c$ gets a **level** $\ell_c \in \{0,\dots,L\}$ and advances with
sub-steps $h_\ell = H/2^{\ell}$, where $H$ is the **macro step** and $L$ the
depth of the map, its deepest level.  The level
comes from the admissible local step $\Delta t^{\text{cap}}_c$ of the cell (the
[local time step](#local-time-step) above) and a safety factor $\sigma \ge 1$
(`safety`, 2 by default, see [Non-local CFL condition](#non-local-cfl-condition); $L_{\max}$ = `max-level`),
for the step $H$ for which the map is built:

$$
\ell_c = \min\Bigl(L_{\max},\ \max\bigl(0,\ \bigl\lceil \log_2(\sigma H/\Delta t^{\text{cap}}_c) - 10^{-12} \bigr\rceil\bigr)\Bigr),
$$

then two rules are imposed, iterated until nothing changes:

- **nesting**: the levels of two cells sharing a face differ by at most one;
- **equal levels across block connections** (types 101 and 201): the two
  sides of a connection always advance together, so a connection is never a
  face between levels: its one flux per stage, which the record that owns the
  face computes for both cells, needs the two sides at the same stage times
  ([Blocks and boundary conditions](#blocks-and-boundary-conditions)), the
  correction stays within a process, and such a map is nested, so Appendix A
  holds for it.  SeisSol keeps the elements of a dynamic rupture face in one time
  cluster in the same way [18].

The macro step is then

$$
H = \min\Bigl(H_{\text{target}},\ \min_c \Delta t^{\text{cap}}_c\,2^{\ell_c}/\sigma\Bigr),
$$

so every cell advances with $h_{\ell_c} \le \Delta t^{\text{cap}}_c/\sigma$.  The
$10^{-12}$ in the levels, about 4500 units in the last place of a double, keeps a cell at the
limit of a level from being pushed one level finer by the rounding of
$\sigma H/\Delta t^{\text{cap}}_c$ and of its logarithm, which stays below $10^{-13}$; the other comparisons
of steps of the integrator use the same tolerance, as Castro does on the landing of its sub-steps [16].

**Local steps.**  The levels, the depth of the map and $H$ come from the
admissible local steps, so each of them must be a positive finite number.  At
the start of every macro step, whatever the keys of `[MOSE-Multirate]`, a cell
whose $\Delta t^{\text{cap}}_c$ is zero, negative, infinite or not a number stops
the run before the map is built: a negative step would make $H$ negative, and
the halving of $H$ while it exceeds a negative $H_{\text{cap}}$ in the redo of
[Non-local CFL condition](#non-local-cfl-condition) would never end; a step
that is not a number can make $H$ not a number, which the floor of the redo
never stops.  Such a step comes from a state that is not physical: a
pressure and a density of opposite signs make the temperature $p/(\rho R)$
negative and the speed of sound not a number; a density and a pressure both
negative make the diffusivity $\mu/\rho$ of a viscous run negative, and the
step negative where the diffusion rate outweighs the convective one; a negative
eddy viscosity can do the same where it makes the diffusivity of every equation
negative and the diffusion rate outweighs the convective rate and the bound of
the stiffness of the turbulence sources; a zero pressure at rest makes the step
of an inviscid run infinite, and a density so small that the diffusion rate
overflows makes the step zero.  A state that is not physical can also give a
positive finite step (a density and a pressure both negative in an inviscid
run): the check does not see it, and the run meets it at the first stage, in
the checks of the faces and of the cells.  The run prints the number of such
cells and, for the first one in the order of the blocks and of $k$, $j$, $i$,
its local step and its primitive variables, and every process stops with the
others (exit code 1).  NaN and infinities are found by their bits, so the test
holds also in the RELEASE build, whose fast floating-point mode can fold a
comparison with a NaN; the DEBUG build, which traps invalid operations,
divisions by zero and overflows, stops with SIGFPE at the first of them, before
the check.

**Depth of the map.**  The map of the adaptive step is built for
$H_{\hat L} = 2^{\hat L}\,\Delta t_{\min}/\sigma$ (at most $H_{\text{target}}$), with
$\Delta t_{\min} = \min_c \Delta t^{\text{cap}}_c$ and $\hat L$ the largest depth up to
$L_{\max}$ at which the cell with the largest admissible step,
$\Delta t^{\text{cap}}_{\max}$, keeps level 0 in the rule above (its own term,
without the interface margin below).  The test is monotone in $\hat L$ and holds for
$\hat L = 0$, so $\hat L$ is found from $L_{\max}$ down; without a target
$\hat L = \min(L_{\max}, \lfloor K + 10^{-12} \rfloor)$, with
$K = \log_2(\Delta t^{\text{cap}}_{\max}/\Delta t_{\min})$ the octaves spanned by
the local steps.  The cell with $\Delta t^{\text{cap}}_{\max}$ has level 0, so
$H_{\hat L} \le \Delta t^{\text{cap}}_{\max}/\sigma$ (to $10^{-12}$), and $H \le H_{\hat L}$: the
target bounds $H$ when it sets $H_{\hat L}$, otherwise the cell with $\Delta t_{\min}$
has level $\hat L$, which nesting and the connections do not raise, nor the margin while
$\kappa \le \sigma$ (the default $\sigma = 2$ with cfl $\le 1$), and
its term in $H$ is $H_{\hat L}$, so the depth of the map is $L = \hat L$.  The macro step
therefore does not exceed the admissible step of the slowest cells divided by
$\sigma$, and the rule leaves no coarse level empty; with $\kappa > \sigma$ the margin can raise that cell, and a coarse
level can stay empty.  Nesting, the margin and
the connections only raise levels, so in a strongly stretched mesh they can
raise every cell of the slowest octave as well; the run shows it as a first count
0 in its line «cells per level».  With a target at or below
$\Delta t^{\text{cap}}_{\max}/\sigma$ every depth passes the test and
$\hat L = L_{\max}$: the target already keeps $H$ within the step of the slowest
cells, and the map is the one built for $\min(2^{L_{\max}}\Delta t_{\min}/\sigma,
H_{\text{target}})$, whose depth $L$ is the one the target gives, at most $L_{\max}$.
On a uniform mesh ($K < 1$) the map has one level and $H =
\Delta t_{\min}/\sigma$.  The depth is that of the map of the adaptive step, also
the one that `levels = fixed` keeps; where $H$ is given, with `fixed-macro-step`
and `levels = recompute` and in the redo of `nonlocal-cfl = redo`, which keeps
the $H$ of its attempt, the levels follow the cells up to $L_{\max}$: for a given
$H > 2^{\lfloor K\rfloor}\Delta t_{\min}/\sigma$ a cut at $\lfloor K \rfloor$
would take the cells with $\Delta t_{\min}$ beyond their limit.  The depth
depends only on the state, as the map, and the largest admissible step is that
of all the MPI processes.

The default $L_{\max} = 10$ is a choice, not a bound of the scheme.  The ideal gain, the cell updates of the single-rate
step over those of the multirate step in the same simulated time, does not fall as the depth grows and stops growing at
$\lfloor K \rfloor$, while every level adds its stages, with their exchanges and synchronisations between the processes,
also when it holds few cells; and with $L_{\max}$ above 10 the floor of the redo is not within the step of a
single-rate run ([Non-local CFL condition](#non-local-cfl-condition)).  Which depth costs least depends on how many
cells lie in each octave of the local steps and on the cost of a stage.

The map is computed at the first macro step and kept (`levels = fixed`) or
recomputed every `recompute-every` macro steps (`levels = recompute`), at the
iterations $1, 1 + k, \dots$ of the run counter.  With `nonlocal-cfl = redo` a
rejected macro step rebuilds the map with levels that only rise
([Non-local CFL condition](#non-local-cfl-condition)); with `recompute-every`
$> 1$ the raised map stays in use until the next rebuild.  One iteration of
MOSE is one macro step.

**Which map for which run.**  For a fixed map the macro step is a one-step
method of order 2 in time: halving $H$ halves every $h_{\ell_c}$ and the error
drops by four.  If instead the map is rebuilt for each $H$, with every cell at
the coarsest level within its limit, halving $H$ lowers by one the level of
every cell finer than level 0, since every ratio that sets a level
($\sigma H/\Delta t^{\text{cap}}_c$, and $\kappa H/\Delta t^{\text{cap}}_v$ with
the interface margin) is proportional to $H$, so its step $h_{\ell_c}$ stays
the same: its error does not drop, and the measured order would be that of the
change of map, not of the scheme.  A study of the temporal order therefore uses
`levels = fixed`: the map of the adaptive step of the first macro step, kept, also with
`fixed-macro-step`.  A production run needs the opposite: as the state changes,
the map for the $H$ actually taken keeps every cell within its limit at the
coarsest admissible level, while a map kept from another $H$ either wastes
steps or exceeds the limits; so `levels = recompute` builds the map for the
$H$ taken, also with `fixed-macro-step`.  With `levels = fixed` the map keeps
every cell within its limit only up to its own step, and it never rises: with an $H$
above that step the finest cells start the first macro step beyond $\Delta t^{\text{cap}}/\sigma$, and every later one
until their admissible step grows by the ratio of $H$ to the step of the map; the run reports it with a warning, once,
and its run line counts the macro steps above the limit, with the largest $h/(\Delta t^{\text{cap}}/\sigma)$, as
with `levels = recompute` where $L_{\max}$ levels are not enough.

The macro step lands exactly on the output instants $t_k = k\,\Delta t_{out}$
($\Delta t_{out}$ = `sol-dtime`, $k$ integer) and on `time-threshold`: when
$t_n + H$ reaches the next of them, or ends within $10^{-12} H$ of it, $H$
becomes $t_k - t_n$ and the time after the step is $t_k$ itself.  The map of
the step is kept: every $h_{\ell_c}$ only shrinks (up to the $10^{-12}$ of a
stretched step), so every cell stays within its limit, and the step is the
same scheme with a shorter $H$, conservative as any other.  A macro step whose
$H$ a redo halves does not land: a later macro step reaches the instant.  The
instants depend only on the time, and the map only on the state and on the run
counter, so a restart from a solution written with its counters (`vtk raw` or
`tecplot ascii`, both checked) goes on as the uninterrupted run (with `recompute-every` > 1, when it falls on
an iteration where the map is rebuilt).

### Non-local CFL condition

The CFL condition is not local: within a step a wave crosses the faces of a
cell and reaches its neighbours, which advance with their own steps.
Gnedin, Semenov and Kravtsov [13, section 2] state two requirements for a
local time stepping scheme:

- **(3)** the step of a cell satisfies its own admissible step and the
  admissible steps of all its neighbours;
- **(4)** if the admissible step of a cell at a time $t_0$ is $\Delta t$, all
  its immediate neighbours end their step at or before $t_0 + \Delta t$.

(3) is necessary but not sufficient [13, fig. 5]: when the admissible step of
a cell drops during the steps of its neighbours (a shock reaching it), a
neighbour that keeps its step delays the wave.  The source asks for (4) at any
time $t_0$ and imposes it in its code, where the tests show no delay; it gives
no proof that (4) is sufficient.

**(3) on the state at the start of the macro step.**  Let $c$ and $v$ share a
face in a block.  Nesting gives $\ell_c \ge \ell_v - 1$ and the levels give
$h_{\ell_v} \le \Delta t^{\text{cap}}_v/\sigma$, so

$$
h_{\ell_c} = \frac{H}{2^{\ell_c}} \le \frac{2H}{2^{\ell_v}} = 2\,h_{\ell_v}
\le \frac{2}{\sigma}\,\Delta t^{\text{cap}}_v ,
$$

and with $\sigma \ge 2$ every cell advances within the admissible step of each
of its face neighbours.  The two cells of a block connection have the same
level, so there $h_{\ell_c} = h_{\ell_p} \le \Delta t^{\text{cap}}_p/\sigma$.
Where the convective limit rules, the step of $c$ measured with the limit of
$v$ has a Courant number of at most $2\,\text{cfl}/\sigma$, which exceeds the
cfl of the limit of $v$ when $\sigma < 2$.  With the interface margin (below: `interface-margin`, and always with
`nonlocal-cfl = redo`) the levels and $H$ also keep $h_{\ell_c} \le \Delta t^{\text{cap}}_v/\kappa$, $\kappa \ge 1.1$,
for every face neighbour $v$ in the block, so (3) holds there for every $\sigma \ge 1$.  So $\sigma \ge 2$ is what gives
(3) without the margin; `safety` takes every value of at least 1.  Without the margin and the certificate of (4)
(`nonlocal-cfl = off`) the Sod tube `mr1d-sod`, with `max-level` 3, the levels recomputed at every macro step and cfl
0.8, breaks down with $\sigma = 1$, after the step taken reaches about twice the step admitted on the new state, and
completes with $\sigma = 2$ (that ratio stays below 0.96) and with $\sigma = 1$ under `nonlocal-cfl = redo`.  With
$L_{\max} = 0$ every cell advances with $H \le \min_c \Delta t^{\text{cap}}_c/\sigma$, (3) holds with $\sigma = 1$
and the macro step is the single-rate step.  The default is $\sigma = 2$, the smallest value that gives (3) without
the margin: the term $\sigma H/\Delta t^{\text{cap}}_c = 2^{L}\Delta t_{\min}/\Delta t^{\text{cap}}_c$ of the levels of
the adaptive step does not depend on $\sigma$, while $H$ is inversely proportional to it, so a larger $\sigma$ costs that
much more work per unit of simulated time.  With the margin, a $\sigma$ below $\kappa$ lets the margin set the levels of
every cell next to a neighbour with a step as small as its own, and a $\sigma$ near 1 makes the cell that sets $H$ start at
its limit, so that every fall of its admissible step within the macro step rejects the attempt (see the floor of
$\kappa$ below): more attempts, not a failure.

**Interface margin.**  On a face between levels the coarse cell $c$ next to
the finer cell $v$ advances with a Courant number, measured with the limit of
$v$, up to $2\,\text{cfl}/\sigma$.  In one dimension, for the advection
equation, the forward-Euler MUSCL step is positive and TVD for a Courant
number up to $1/(1+\theta)$: with
$u^-_{j+1/2} = \bar u_j + \tfrac12\phi(r_j)\,\Delta_+\bar u_j$ and
$r_j = \Delta_-\bar u_j/\Delta_+\bar u_j$ (`Lib_Reconstruction`),

$$
u^-_{j+1/2} - u^-_{j-1/2} = \Delta_-\bar u_j\,\Bigl[1 + \frac{\phi(r_j)}{2r_j}
- \frac{\phi(r_{j-1})}{2}\Bigr],
$$

so that, with $B_j$ the bracket, in $[0, 1+\theta]$, and $\nu\,(1+\theta) \le 1$,
the step $u^{n+1}_j = (1-\nu B_j)\,\bar u_j + \nu B_j\,\bar u_{j-1}$ is a
convex combination and
$u^{n+1}_j - u^{n+1}_{j-1} = (1-\nu B_j)\,\Delta_-\bar u_j + \nu B_{j-1}\,\Delta_-\bar u_{j-1}$
has non-negative coefficients, so the total variation does not grow (the
argument of Hundsdorfer [8], §6.4, p. 54, where the limiter condition (5.4) of his
§5.2, p. 39, bounds the bracket by $1+\mu$):
$\theta = 1$ for van Leer, MC and superbee ($\phi$ and $\phi/r$ in $[0, 2]$),
$\theta = 1/2$ for minmod ($[0, 1]$), $\theta = \sqrt2/2$ for van Albada
($[-0.207, 1.207]$), $\theta = 0$ with the first-order reconstruction.  These bounds hold on cells of equal width. `Lib_Reconstruction`, called by `Lib_Convective`, limits the slopes divided by the
distances between the cell centres and bounds the limited slope so that the value of a cell at each of its faces lies between
the value of the cell and that of the neighbour across the face
([Spatial Discretization](numerics.md#muscl-reconstruction)), so on cells of different widths the face value is
$u^-_{j+1/2} = \bar u_j + \hat\phi_j\,w^+_j\,\Delta_+\bar u_j$, with $w^+_j = \Delta x_j/(\Delta x_j + \Delta x_{j+1})$,
$w^-_j = \Delta x_j/(\Delta x_{j-1} + \Delta x_j)$, $r'_j = r_j\,(\Delta x_j + \Delta x_{j+1})/(\Delta x_{j-1} + \Delta x_j)$ and
$\hat\phi_j = \min(\phi(r'_j),\ 1/w^+_j,\ r'_j/w^-_j)$ for $r'_j > 0$ and $\hat\phi_j = \phi(r'_j)$ for $r'_j \le 0$, which is 0
for every limiter but van Albada.  For a limiter whose $\phi$ and $\phi/r$ lie in $[0, 2\theta]$ (all but van Albada) the
bracket lies in $[1 - \min(2\theta w^+_{j-1}, 1),\ 1 + \min(2\theta w^-_j, 1)] \subset [0, 2]$ on every mesh (with van Albada
in $[w^-_j - 0.207,\ w^-_j + 1.207]$), and the step is a convex combination when $\nu_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$
([Appendix B](#appendix-b-total-variation-at-a-face-between-levels), Lemma 1).  With one time step for all the cells this
follows from $\nu_j\,(1+\theta) \le 1$ on every cell.  Hoang
et al. [11, Theorem 4.2] prove a total-variation bound for a one-dimensional
scalar law, discretised by piecewise-linear discontinuous Galerkin elements
with the second-order strong-stability-preserving Runge–Kutta scheme in
two regions joined by one interface cell, with their limiter condition (4.14)
on the reconstructed values, for a solution regular at the interface and under
the hypothesis (4.15) that the coarse step meets this bound on the first fine
cell.  With `interface-margin = .true.`, and
always with `nonlocal-cfl = redo`, the levels and $H$ also satisfy
$h_{\ell_c} \le \Delta t^{\text{cap}}_v/\kappa$ for every face neighbour $v$ in
the block, $\kappa = \max(\text{cfl}\,(1+\theta),\ \kappa_0)$, $\kappa_0 = 1.1$, with the cfl
in force (a reading of `input.ini` during the run can change cfl within the bound below, and $\kappa$
follows it from the next macro step): the level of a cell comes from
$\max(\sigma H/\Delta t^{\text{cap}}_c,\ \kappa H/\min_v \Delta t^{\text{cap}}_v)$
and $H_{\text{cap}}$ takes the term $\min_v \Delta t^{\text{cap}}_v\,2^{\ell_c}/\kappa$.
With $\sigma \ge \kappa$ this adds nothing on the faces between equal levels or
towards coarser cells, and nesting and equal levels across connections only
raise levels.  It raises a ring of cells on the coarse side of a face between
levels where the fine cell runs above $1/\kappa$ of its limit.
The coarse cells run with $\text{cfl}/\sigma$, within $1/(1+\theta)$ when
$\text{cfl} \le \sigma/(1+\theta)$: with the default $\sigma = 2$ and $\theta \le 1$, $\sigma/(1+\theta) \ge 1$, so every
cfl up to 1 meets it for every limiter (with $L_{\max} > 0$ and the redo of `nonlocal-cfl = redo` the set-up and every
reading of `input.ini` reject a cfl above 1: the bounds of cfl and vnn, below); a $\sigma$ below 2 meets it only for
$\text{cfl} \le \sigma/(1+\theta)$, a property and not a check of the set-up.  This is the hypothesis of the theorem, not the theorem: the scheme is
a finite-volume MUSCL scheme, in several dimensions, for systems, with several
levels.  At the faces between levels the monotonicity is proven for linear
two-cell fluxes, the first-order upwind flux of a given velocity field without
discrete divergence and the central diffusion
([Appendix A](#appendix-a-order-and-stability-of-the-multirate-macro-step),
Theorem 2).  For the advection equation in one dimension, with the limited reconstruction of `Lib_Reconstruction`, a limiter
whose $\phi$ and $\phi/r$ lie in $[0, 2\theta]$ (van Albada, whose $\phi$ is negative for $-1 < r < 0$, is outside) and two
levels, on cells of any widths, one macro step raises the total variation at most by a term local to the face between levels,
a multiple of the Courant number of the coarse cell of the face times the differences of the solution next to the face. Summed
over a run, this bounds the total variation of a solution that is Lipschitz
continuous next to the face, as in Hoang et al. [11, Theorem 4.2], without their hypothesis (4.15). The derivation asks every
cell to satisfy $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, with $c_j$ its own Courant number.  With
$\text{cfl}/\sigma \le 1/2$ (cfl $\le 1$ and the default $\sigma = 2$) every cell has $c_j \le 1/2$, and with it this condition, on every mesh
and for every such limiter: also, with minmod, at the coarse cell downstream of a face from a finer to a coarser cell, where
it reads $\lambda\,(1 + r/(r+1)) \le 1$ with $r$ the ratio of the coarse to the fine width, which $c\,(1+\theta) \le 1$ alone
does not give ([Appendix B](#appendix-b-total-variation-at-a-face-between-levels)). Several levels, several dimensions, systems and the
positivity at the faces between levels are outside this result: for them, and for the limited reconstruction outside it, the
stability at the faces between levels rests on the tests of [Multirate Verification](../vv/multirate.md).

**The floor of $\kappa$.**  $\text{cfl}\,(1+\theta)$ is the hypothesis above; $\kappa_0$ is the margin of condition (4)
below.  With $\kappa \le 1$ and $\sigma = 2$ the margin would raise no cell above the level that nesting gives it: the cell that sets $H$, next to a coarser cell,
would start with $\hat h_v = \Delta t^{\text{cap}}_v$ exactly, any fall of its admissible step at a later start of a
sub-step shared with that cell would reject the attempt, and the next macro step, whose map comes from the state alone, would
meet the same limit.  With $1 < \kappa \le \sigma$ (the default $\sigma = 2$ with cfl $\le 1$) a rejection for a violation of
(4) needs a fall of the admissible step of a cell of more
than $1 - 1/\kappa$ within the macro step, on every mesh and state, for every limiter and number of dimensions and of
processes ([Appendix A, Proposition](#proposition-the-margin-on-condition-4)); a front is rejected only when
it makes the step fall by more than that, and then only where the pair of cells runs that close to its limit: a fall by more
than $1 - 1/\kappa$ is necessary for a rejection, not sufficient.  That proposition does not fix the value.  $\kappa_0 = 1.1$ is a practice taken from Castro, whose
check of the step tolerates the same quantity, the fall of the admissible step within the step, up to $1/1.1$ (`change_max`
[16]); it is not a quantity that the scheme gives.  PETSc and DOPRI5 take a step of 0.9 of the one they estimate on the same
principle, that a rejection be rare [14, 17], with another quantity, the error of the step.  Castro accepts a step
above its limit by up to that factor; here the margin is built into the map, and condition (4) stays imposed as stated.  The
margin raises a cell only where the admissible step jumps by more than $\sigma/\kappa$ across one of its faces, and
the rule of a cell reads the local steps of its neighbours, not their levels, so it does not propagate beyond those
neighbours.  How many cells it raises depends on the mesh and on the state: the run prints the cells whose level the margin term raises
before nesting.  Below
$\text{cfl}\,(1+\theta) = 1.1$ the margin does not depend on cfl and no key turns the floor off; a `safety` above 2 gives more
margin at the price of a shorter $H$.

**(4) with `nonlocal-cfl`.**  A state exists only at the start of the
sub-steps, so (4) is checked there: for every cell $v$ and every start $\tau$
of one of its sub-steps,

$$
\hat h_v(\tau) \le \Delta t^{\text{cap}}_v(\tau), \qquad
\hat h_v = 2\,h_{\ell_v} \text{ if a face neighbour one level coarser starts a
sub-step at } \tau, \quad \hat h_v = h_{\ell_v} \text{ otherwise}
$$

(finer neighbours start and end within the sub-step of $v$, nesting limits the
jumps to one, connections have none).  At $t_n$ this is (3).  The certificate
evaluates $\Delta t^{\text{cap}}_v(\tau)$ with the local-step formula of the
single-rate step, cell by cell, in the first stage of every sub-step after the
turbulence sources, on the start state of the sub-step; a value that is not a
number is a violation, and $\hat h_v$ above $\Delta t^{\text{cap}}_v$ by less than
$10^{-12}$ of it is not, the tolerance of the formula of the levels: a cell that
starts at its limit gives a ratio of 1 within some units in the last place, which
no map and no $H$ would change.  Two cells that share only an edge and are two levels
apart satisfy neither (3) nor (4) ($h_c = 4h_d \le 2\Delta t^{\text{cap}}_d$
with $\sigma = 2$): a wave reaches a cell only through its faces, along which
(4) holds pair by pair.

- `off`: (4) is not imposed; the run reports $r_{post}$ and $r_{nb}$ (below).
- `count`: the certificate counts the violations per level and the macro steps
  that would be rejected, printed at the end of a run that has some; the results
  do not change.
- `redo` (the default): a macro step whose certificate fails is rejected.  The state goes back
  to $t_n$, the map is rebuilt with the same $H$ from the smallest admissible
  step that every cell showed in the attempt, so the levels only rise, up to
  $L_{\max}$ (the map can become deeper than that of the first attempt), and $H$
  is halved while a cell would need a level above $L_{\max}$ or a value was not a
  number.  An accepted macro step satisfies (4) at every start of a sub-step.
  With $L > 0$ an attempt is rejected also when a cell or a face takes a state
  that is not physical: a cell whose update gives a conservative variable that
  is not finite or a total density that is not positive (tested before the
  conversion to primitive variables, which divides by the density), or a
  primitive variable that is not finite, a negative pressure or a temperature
  outside the tables of the gas (after the floor of the partial densities, the
  test of the faces); an internal face of a block or a face between levels
  whose reconstructed left or right state has a negative pressure or a
  temperature outside the tables (the faces of block connections and of
  boundaries have no such test, as in the single-rate step).  The cell keeps the
  start state of its sub-step and the face gives no flux; from the next stage
  on, the process that took the state computes no internal face and updates no
  cell (its exchanges go on, so the processes stay together), and the attempt is
  rejected on every process.  The state goes back to $t_n$ and the attempt is
  redone on the same map and local steps with $H$ halved once: the admissible
  steps that the certificate sees after such a state come from states that are
  not the solution, so they do not lower the local steps.  A face that is not
  physical at the first stage of a level at $t_n$ sees only the state at $t_n$,
  which a shorter $H$ does not change: there the run stops with the lines of the
  face, as the single-rate step does.  In a macro step accepted after attempts
  that the redo can redo (`levels = recompute`, no `fixed-macro-step`, $L > 0$)
  every cell ends with finite values, partial densities not below the floor, a
  pressure that is not negative and a temperature within the tables.
  The attempts end: the levels are bounded by $L_{\max}$, and as $H$ halves the states
  of the step tend to those at $t_n$, where (3) holds with the margins (the
  local step is continuous in the state where density and pressure are
  positive); every rejection for a state that is not physical halves $H$, so a
  macro step has at most 11 of them; below $2^{-10}$ of the first $H$ the run
  stops with a message that gives the rejections of the macro step by cause and,
  when its last attempt took a state that is not physical, the first one of that
  attempt in the order of the stages.  With a map of depth $L \le 10$, the default cap, the
  first $H$ is at most $2^{L}\Delta t_{\min}/\sigma$, so this floor is at most
  $\Delta t_{\min}/\sigma$: every cell then advances with at most the step of a single-rate
  run with the same margin.  Their number is bounded too, by a bound
  that grows with the mesh: $H$ takes at most 11 values, and with the same $H$
  the local steps only decrease, so the levels of the rebuilt maps only rise,
  and every rejection that keeps $H$, but the one of `force-redo`, raises a cell
  by at least one level; over the $N$ cells of the mesh a macro step has at
  most $1 + 11\,(N L_{\max} + 2)$ attempts (the term 1 is the attempt rejected by
  `force-redo`).  $H$ is halved, the factor of the retry of Castro when it has no suggested
  step [16], and not set to $H_{\text{cap}}$, so that it takes these 11 values, on which the
  bound rests; the macro step can then come out up to twice shorter than $H_{\text{cap}}$.
  `max-attempts` (21 by default) limits them: when that many
  attempts of a macro step are rejected (the attempt rejected by `force-redo` is
  not counted), the run stops with the same message, which names the limit in
  place of the floor and, for a last attempt without a state that is not
  physical, gives its violations of (4) or its local step that is not a
  positive number.  The floor and the limit are rules that stop a run with a
  message, not conditions of the scheme: the bound above ends the attempts
  without them.  The floor is 10 halvings of $H$, as many rejections as PETSc
  allows a step by default [14] and as many failures of the nonlinear solver or
  of the inequality constraints as end a step of the ARKODE integrators of
  SUNDIALS [15]; there every rejection shortens the step.  21 is one attempt for
  each of the 11 values of $H$ that the floor allows and 10 more for the
  rejections that keep $H$ and raise levels, which those integrators do not
  have.  The limit is checked after the floor, so with `max-attempts` $\ge 11$ a
  macro step whose rejections all halve $H$ ends as without it.
  With `fixed-macro-step` or `levels = fixed` the certificate only counts: there
  (4) is not imposed, and, as with one level, the stops of the single-rate step
  stay: a cell whose update gives a value that is not finite or a negative
  pressure, or a face that is not physical, stops the run.  With `recompute-every` $= k > 1$ the map raised by a redo
  stays in use until the next rebuild of the counter.  With $L = 0$ there is one
  level and the check is the own limit at $t_n$, true by construction.  When a
  macro step was redone, the run prints at its end the attempts per macro step,
  the macro steps redone, the halvings of $H$, the macro steps whose map a
  redo made deeper than that of their first attempt and, when there were some,
  the attempts rejected for a state that is not physical (with a cell, with a
  face); the counts of the violations of the certificate leave out these
  attempts; the stage times and the counts
  of negative turbulence values include the rejected attempts.  The test device `force-redo` rejects the first attempt of
  every macro step that the certificate accepts and redoes it with the same $H$
  and map: the state goes back to $t_n$ to the bit, so the results do not
  change.

**The bounds of cfl and vnn.**  At $t_n$ the map runs every cell with $h \le \Delta t^{\text{cap}}/\sigma$; with
`nonlocal-cfl = redo` the certificate accepts, at every later start $\tau$ of a sub-step, $\hat h_v$ up to
$\Delta t^{\text{cap}}_v(\tau)$ itself, without $\sigma$ (to the $10^{-12}$ of its test), so within the macro step a cell
can advance with the whole step that cfl and vnn admit on its state at $\tau$.  On a uniform mesh of equal cells, for the
frozen-coefficient advection–diffusion model with the first-order upwind flux and the central diffusion, the
checkerboard datum $u_j = (-1)^j$ gives a zero limited slope with every reconstruction of MOSE (first order, minmod,
superbee, van Leer, van Albada, MC), so on it every limited scheme is the first-order scheme, and the datum is an
eigenvector of Heun with the factor $R(z) = 1 + z + z^2/2$, $z = -2\,\text{cfl}\,a - 4\,\text{vnn}\,b$, where $a$ and $b$
($a + b = 1$) are the convective and diffusive shares of the local step of a cell at its limit; $|R(z)| > 1$ exactly when
$\text{cfl}\,a + 2\,\text{vnn}\,b > 1$.  Hence a sub-step with the whole admissible step is linearly stable on every cell
only if cfl $\le 1$ and vnn $\le 1/2$ (pure convection: $R(-2\,\text{cfl}) = 1 - 2\,\text{cfl} + 2\,\text{cfl}^2 > 1$ for
cfl $> 1$; pure diffusion: $R(-4\,\text{vnn}) = 1 - 4\,\text{vnn} + 8\,\text{vnn}^2 > 1$ for vnn $> 1/2$), for every
limiter; with a larger value it is not, on a cell where that number sets the step (with cfl $= 1$ and vnn $> 1/2$ on every
cell with $b > 0$; with cfl $< 1$ where $b > (1 - \text{cfl})/(2\,\text{vnn} - \text{cfl})$).  These bounds are also
sufficient for the first-order reconstruction and for the limiters with a linearisation (van Leer, MC, van Albada)
([Local time step](#local-time-step)), and are the bound of the other slopes for minmod and superbee, for which no linear
limit is proven: their one-sided linearisations are the fully upwind slope, whose Heun limit on the checkerboard mode is
1/2, and the central slope, which Heun does not integrate stably at any Courant number; what is proven for them is the
bound of the total variation of the forward-Euler step, 2/3 for minmod and 1/2 for superbee, on equal cells in one
dimension (interface margin, above).  In two and three dimensions, with the additive rule of the local step, the bounds
are those of the checkerboard in every direction and are the same numbers (on cells skewed in all three pairs of
directions the viscous rate takes a factor that raises the diffusive threshold on those cells).  For the Euler and
Navier–Stokes equations linearised about a uniform state the convective bound is attained by the fastest wave with the
Roe, Rusanov and HLL-type fluxes in one dimension, in several dimensions with Rusanov on every state and with Roe and HLLC
on a state at rest; the diffusive bound is attained by the scalar that owns the largest diffusivity (the energy with
gases of Prandtl number near 0.7, the species, the turbulence variables) in the limit of a cell where the diffusion sets
the step.  The bound vnn = 1/2 is proven on a uniform mesh of equal cells
of any skew and aspect ratio (in 3-D with the factor that the local step takes on cells skewed in all three pairs of
directions), with frozen coefficients and away from boundaries, for the scalars exactly and for the momentum
conservatively ([Local time step](#local-time-step)); on cells of varying size or shape, next to boundaries and at the
axis of the wedge the local-step rule applies the same bound cell by cell, and its support there is the evidence of
that section.  So with $L_{\max} > 0$ and `nonlocal-cfl = redo`, where the certificate redoes the attempts
(`levels = recompute` and no `fixed-macro-step`), the set-up
and every reading of `input.ini` reject a cfl above 1 and, in Navier–Stokes runs, a vnn above 1/2: these bounds are
necessary for every sub-step that the certificate accepts to be linearly stable there.  A larger value would not show at
$t_n$, where the cells run within $1/\sigma$ of their step.  That a run with a larger value fails is not proven: an
unstable sub-step occurs only where the admissible step of a cell falls within the macro step by more than
$\sigma/\text{cfl}$ or $\sigma/(2\,\text{vnn})$, for the sub-steps of that macro step that remain.  With `nonlocal-cfl =
off` and `count`, and with `redo` and `levels = fixed` or a `fixed-macro-step`, where the certificate only counts,
nothing bounds the sub-steps after the start of the macro step, where every
cell runs within cfl/$\sigma$ and vnn/$\sigma$, and the set-up checks neither bound.  The two conditions are distinct: at
the start of the macro step every cell runs within $\Delta t^{\text{cap}}/\sigma$, so every cell, the coarse cells of a
face between levels included, has on its own width a Courant number of at most $\text{cfl}/\sigma$, 1/2 with cfl $\le 1$
and the default $\sigma = 2$, within the bound $1/(1+\theta)$ of the total variation above for every limiter (measured
with the limit of the finer neighbour, the coarse step reaches $2\,\text{cfl}/\sigma$; the interface margin, on with
`redo` or `interface-margin`, keeps it within $1/(1+\theta)$); within the macro step, with `redo`, the certificate keeps
every sub-step within the linear stability limit of Heun, on the state where the sub-step starts.  For a cell next to a wall where the diffusion sets the step, with 1 µs its step
at the limit of Heun (a diffusion number $h\,D_{\max}\sum_d|\nabla\xi_d|^2$ of 1/2 on its state):

| | step at the start | diffusion number at the start | the certificate rejects above a diffusion number of |
|---|---|---|---|
| single rate, vnn 0.5 | 1 µs | 0.5, the limit | — (the step is computed again at every step) |
| multirate, vnn 0.5 | 0.5 µs | 0.25 | 0.5, the limit |
| multirate, vnn 0.98, which the set-up rejects with the redo | 0.98 µs | 0.49 | 0.98: between 0.5 and 0.98 an accepted sub-step is not linearly stable |

For the convection, cfl $\le 1$ puts the threshold of the certificate at the linear limit of Heun in the same way.  With
$L_{\max} = 0$ the macro step is the single-rate step, whose limits the set-up does not check
([Time step](../user/input.md#time-step-cfl-and-vnn)).

**Measures without the certificate.**  With `nonlocal-cfl = off` the levels and
$H$ come from the state at the start of the macro step and stay fixed during
it: if the admissible step of a cell drops within the macro step, its
neighbours keep their steps.  The run measures (4) on the state reached at
$t_{n+1} = t_n + H$, whose admissible steps the next macro step computes
anyway:

$$
r_{post} = \max_c \frac{h_{\ell_c}}{\Delta t^{\text{cap}}_c(t_{n+1})},
\qquad
r_{nb} = \max_c\, \max_{v \in N(c)} \frac{h_{\ell_c}}{\Delta t^{\text{cap}}_v(t_{n+1})},
$$

with $N(c)$ the face neighbours of $c$ in its block; across a block
connection $h_{\ell_c} = h_{\ell_p}$, so the ratio with the partner $p$ is the
$r_{post}$ of $p$.  The last sub-step of $c$ ends at $t_{n+1}$: if the
admissible step of $v$ does not fall below $\Delta t^{\text{cap}}_v(t_{n+1})$
during that sub-step, $r_{nb} \le 1$ gives (4) for it, and $r_{post} \le 1$
gives the own limit of $c$; a value above 1 flags a possible violation.  The
measure is partial: it sees only the state at the end of the macro step (not
the earlier sub-steps) and only the face neighbours in the block (not the
cells that share an edge or a vertex).  Every `shell-diter`
iterations the run prints both values of the last macro step and, when a cell
is above 1 or above $1/\text{cfl}$ (a Courant number above 1 where the
convective limit rules), the number of such cells, and at the end, for each of the two that went above 1
in a macro step, its largest value over the run with the number of macro steps
that had a cell above 1 (also for $r_{post}$ or $r_{nb}$ together).

### Recursive Heun scheme

With $\mathbf{U}^0$ the state of a cell at the start of its sub-step, the
macro step is the recursion below.  It is the *Prop. RK2* scheme of
Jeanmasson, Mary and Mieussens [10] with two differences: the coarse corrector
reads the final values of the finer cells (their Table 16 reads the predictor
of the fine cells), and the recursion over several levels, with bands of one
cell, is ours (they write the scheme for two levels).  Both keep the
conservation, which needs only that the two sides of a face receive the same
flux with the same weights.

Both differences keep the second order on a fixed mesh
([Appendix A](#appendix-a-order-and-stability-of-the-multirate-macro-step),
Theorem 1).  Every stage vector is first-order consistent on the stencil it
is read from (the coarse cells through the interpolant, the finer cells at
their final values, the cells of the level at the start state or at the
predictor), and every cell receives each face flux with the weights of Heun,
sum $H$ and first moment $H^2/2$, also on the corrected side of a face between
levels, where the coarse cell receives $h_m A = \tfrac{h_m}{4}\sum_j G_j$ at
the four fine stage times: the two differences enter only the $O(H^3)$ term,
so the macro step is a one-step method of order 2, uniformly at the faces
between levels, for every linear semi-discretisation and for $f \in C^2$ (the
limiter and the switches of the turbulence models break $C^2$ where they
act).  This is what `mr-order`, `mr-deep` and `mr-bc` measure: fixed mesh,
fixed map, $H \to 0$.  At a fixed Courant number ($H \propto \Delta x$) the
local error of a macro step is $O(H^2)$ in the cells next to a face between
levels and $O(H^3)$ elsewhere (first order locally, as the conventional
boundary data of [7]; Appendix A, Remark 2).  When the map is rebuilt
(`levels = recompute`) and a face between levels travels with a wave at the
speed of the wave, the errors of order $H^2$ that each macro step leaves next
to the face stay on the same part of the wave and can add up.
[Appendix A, Corollary](#corollary-every-sequence-of-maps) bounds the error
for every sequence of nested maps, this one included: on a fixed mesh by
$O(H^2)$, under (H2)–(H4) of Theorem 1; at a fixed Courant number, for the
linear two-cell fluxes of Theorem 2, by $O(H)$ in the maximum norm and by
$O(H^2)$ in the volume-weighted $L^1$ norm (with faces between levels of
bounded total area).  Outside that linear setting no
bound is proven at a fixed Courant number: the test `mr-moving-levels` of
[Multirate Verification](../vv/multirate.md) measures the $L^1$ order of this
case on a simple wave of the Euler equations (second order); its global
$L^\infty$ order, near 1.3 as for the single-rate step with
the same limiter, has its largest error at the maximum of the wave and does
not isolate the cells next to the moving faces, where the test reports,
without deciding on it, a lower $L^\infty$ order (1.67 on its last pair,
against 2.02 to 2.03 for the single-rate step on the same cells).  The
stability of the two differences is proven for the first-order upwind flux
of a velocity field without discrete divergence and the central diffusion:
the macro step is monotone (non-negative matrix,
maximum principle, conservation, norms 1, 2 and $\infty$ at most 1) for every
nested map and every number of levels when every cell advances within 0.86 of
its forward-Euler positivity step in the setting of Theorem 2 (linear two-cell fluxes, Appendix A): on a mesh of
equal rectangular cells with constant coefficients, where the diffusive flux of a scalar is the two-cell flux of
Theorem 2, the additive local-step rule
starts every cell within $\max(\text{cfl}, 2\,\text{vnn})/\sigma$ of that step, 0.4 with $\sigma = 2$, cfl 0.8 and
vnn 0.4, and at most 1/2 with cfl $\le 1$, vnn $\le 1/2$ and the default $\sigma = 2$.  The bounds that the set-up
enforces with the redo of `nonlocal-cfl = redo` are necessary conditions of linear stability
([Non-local CFL condition](#non-local-cfl-condition)), not of positivity: on skewed cells the face gradient also takes
the central differences of the two cells along the face (`Lib_Diffusive`), whose cross metric terms can give corner
neighbours a negative weight (on equal parallelograms they do, so that no step makes that diffusion positive), and
Theorem 2 does not cover it.  With a limited
reconstruction the coarse corrector is not a convex combination of
forward-Euler steps (Appendix A, Remark 4): for the advection equation in one
dimension, with two levels and with three levels and a band of one cell, Appendix B bounds the growth of the total variation
by terms local to the faces between levels, and elsewhere the stability at the
faces between levels rests on the tests of
[Multirate Verification](../vv/multirate.md).

```text
ADVANCE(l, tau):  U0 := U on the cells of level l
                  STAGE(l, 1, tau)                       predictor
                  if l < L:  ADVANCE(l+1, tau);  ADVANCE(l+1, tau + h_l/2)
                  STAGE(l, 2, tau + h_l)                 corrector
STAGE(l, k, t):   coarser cells read by the stage <- interpolant at time t
                  ghost cells read by the stage at time t
                  viscous runs: viscous variables of the cells and ghost cells read by the stage
                  fluxes of the faces with a cell of level l (direction by direction, the convective
                  then the diffusive faces, as the single-rate step), boundary fluxes
                  turbulence models: turbulence sources of the cells of level l
                  U := U0 + c_k (U - U0 + h_l L(U)),  c = (1, 1/2),  on the cells of level l only
```

A macro step with all levels $0\dots L$ populated has $2^{L+2}-2$ stages.
When a stage of level $\ell$ runs, every finer cell is already at the stage
time (finer levels complete their two sub-steps before the coarser
corrector), and every coarser cell $m<\ell$ that the stage reads receives the
**linear interpolant in conservative variables** between its start state and
its predictor (Hoang, Ju, Leng and Wang [11], eq. (3.2), with the predictor of
their App. B.1, eq. (B.4): the full Heun stage of the coarse cell with its true
neighbours, as stage 1 above, not the step 1 of their section 3, which uses the
coarse cell alone; the same values as [10], Tables 13 and 16):

$$
\mathbf{U}_c(t) = \mathbf{U}^0_c + \eta\,(\mathbf{U}^{(1)}_c - \mathbf{U}^0_c),
\qquad \eta = \frac{t - \tau_m}{h_m},
$$

with exact copies at $\eta = 0$ and $\eta = 1$.  Stage times are dyadic
fractions of $H$ and are handled as integers, so $\eta$ is exact.  Boundary
conditions are evaluated at the stage time.

The halo can hold cells two or more levels coarser (bands of one cell, below).
The interpolant has an error $O(h_m^2)$ for every $\eta$ in $[0, 1]$, and the
equations (B.3)–(B.6) of [11] hold for any integer ratio of the steps, so these
cells are read with the same local order.  They exchange no flux with the
cells that read them: only face fluxes enter the balance of a cell, and with
the nesting and the equal levels across connections a face joins levels at
most one apart, so the correction below keeps the conservation.  For linear
two-cell fluxes (first-order upwind, central diffusion) Theorem 2 of
Appendix A covers them (every nested map, bands of one cell included); with the limited reconstruction
Appendix B bounds the total variation for the advection equation in one dimension with three levels and a band of one
cell (Theorem 3b), and their effect on the stability is tested by `mr-deep` (bands of one cell, 2000
macro steps).

The coarser cells read by a stage (the **halo**) are the union of the true
stencils of the routines the stage calls: the four cells along the axis of
every convective face and the cells read through the ghost cells of the
boundary records.  With the laminar Navier–Stokes equations the halo also
holds the ten cells of every diffusive face (the two cells of the face and
their neighbours $\pm 1$ in the two other directions, diagonals included)
and the four tangential neighbours of the partner cell of every block
connection of the level, which the viscous connection flux reads through the
values `Pg(:,3:6)` refreshed at every stage after the ghost cells.  The
viscous fluxes take the viscous variables of these cells (mass fractions,
velocity, temperature and turbulence variables per unit mass) from an array
that the single-rate stage fills once for every cell after the ghost cells; a
stage of level $\ell$ fills it, after the interpolation and the ghost cells, on
the cells of the level, their face neighbours and their diagonal neighbours in
every pair of directions, ghost cells across one face included: the coarser
cells at the stage time, the finer ones at their current values, so that no
flux reads a value computed from an older state.  The
turbulence sources of a cell read its neighbours $\pm 1$ along the axes and its
own ghost cells, already in that union.  Bands of one cell are therefore
allowed.

### Conservative correction

On a face between a coarse cell (level $m$) and a fine cell (level $m+1$) the
two sides would integrate different fluxes.  The coarse side is corrected so
that it receives exactly the impulse of the fine side (the refluxing principle
of Berger and Colella [12], eq. (2), stated there for a one-step scheme; the
form for the two Heun stages is derived here, with the algebra of [11],
eq. (3.5)): at its stage 1 the coarse cell uses and saves its own flux $G_s$; each
of the four fine stages adds $\tfrac14$ of its flux to an accumulator $A$
(weight $b\,h_{m+1}/h_m = \tfrac12\cdot\tfrac12$); at its stage 2 the coarse cell
uses $2A - G_s$.  Over its step the coarse cell then receives
$h_m A = \sum_{\text{fine stages}} \tfrac12 h_{m+1} G$, the same as the fine cell,
and mass, momentum and energy are conserved to round-off.  The convective and
the diffusive flux of a face are corrected each with its own accumulator, so
every cell receives its contributions in the order of the single-rate step.

### Blocks and boundary conditions

The level map spans all the blocks.  Because the two sides of a block
connection (types `101` and `201`) always have the same level, a connection
is never a face between levels: at the stages of that level the record that
owns the face computes its one flux for both cells, as the single-rate step
does, from ghost cells that every stage of the level refreshes from the
partner block at the stage time.  The cells of
the partner block read through those ghost cells are part of the halo, so a
coarser one takes its interpolant.  Equal levels cost a few cells: a
connection on the edge between two bands raises the cells of the coarser side
next to it, and the nesting rule follows.

Every boundary record is evaluated at the stages of the level of its cell,
with `domain%time` equal to the stage time
$t_n + \tau$ of that level.  The time-dependent data of the inlet `402` and of
the assigned state `410` are therefore read at the stage times of each level,
the conventional imposition of [7] applied level by level; with one level the
times coincide with those of the single-rate Heun step (section
[Stage times](#stage-times-and-time-dependent-boundary-data)).

### Turbulence models

With the Spalart–Allmaras (`SA`) or SST model the stage of a level ends, as the
single-rate stage, with the turbulence sources of its cells, evaluated on the
state of the stage; they read the ghost cells of the fluxes of the stage (at
the walls those of the wall values, one set for the fluxes and the sources).  In
time-accurate runs the sources are explicit and the admissible local step
includes a bound $\hat d_c$ of the rate of their destruction, $D_q/q$, on the state of the cell at the start of the
macro step alone (for the Spalart–Allmaras model twice the largest value of that rate, for SST 1.84 times the rate of
$k$ and at least twice that of $\omega$) ([turbulence sources](turbulence.md#time-accurate-runs-explicit-sources)):
$1/\Delta t^{\text{cap}}_c$ gains $\hat d_c$.  The levels and the macro step
inherit it, so every cell advances with
$h_{\ell_c}\,D_q/q \le h_{\ell_c}\,\hat d_c \le 1/\sigma$, and each stage of a level is a forward-Euler
step on its cells that keeps $\tilde\nu$, $k$ and $\omega$ positive against
their destruction (the coarser cells read by the stage take a convex
combination of two states, the interpolant with $0 \le \eta \le 1$).  The
correction on the faces between levels is not a convex combination, so the
positivity there rests on the step margin and is checked by the tests
(`mr-positivity`, and the count of negative values that a RANS run prints at
its end when a stage made one negative).  With one level the stage calls the block routines of the
single-rate step, so the macro step is the single-rate step to the bit; with
more levels it calls their per-cell copies on the cells of the level.

### Several processes (MPI)

With MPI every block belongs to one process.  The level map is the one of a
single process: at every sweep of the two rules the levels of the partner
cells of the connections between processes are exchanged, and the sweeps stop
when no process changes its map; $H$ and $L$ are global.  Every process runs
the same recursion, so the stages of a level are simultaneous everywhere.

A stage of level $\ell$ exchanges only what it reads from other processes:
the ghost cells of the connection records whose partner block is on another
process and whose ghost positions the stage reads, in viscous runs the
tangential ghost values of the connection records of its cells, and the share
of the flux of the connection faces whose owner record is on the other
process.  The process
of the partner block sends them at the stage time, as the single-rate
exchange does: cells of level $\ell$ or finer as they are, coarser cells with
their interpolant, which that process evaluates as part of its own halo of
level $\ell$.  The two processes agree on these lists, record by record and
level by level, every time the lists are built (first the tangential values,
whose sources the partner's process may read from a third process, then the
ghost cells).  Because both sides of a connection have the same level, no
face between levels crosses a connection and the conservative correction
stays within a process.  The MPI calls are made by the master thread only
(`MPI_THREAD_FUNNELED`), between barriers of the OpenMP team.

### Scope

The multirate integrator runs time-accurate Heun (`time-scheme = RK2`) with
conservative integration variables.  At set-up or at its first macro step it stops with an error
(non-zero exit code and a message `[MOSE-Multirate] <item>`) for the
configurations outside its scope: other time schemes, primitive or
preconditioned variables, implicit residual smoothing, the shock detector
(`shock-detector`, also implied by `HLLC+Tramel`, `HLLC+Chen`, `HLLE++`, `AUSM+M`), chemistry,
soot, the rotating frame, multigrid sequences, steady (local time step) runs,
the CFL ramp, chimera, coupled (103), SRM-grain boundaries in mode 1 (502, a wall that blows the products in; a
grain in mode 2, an inflow, runs as the mass-flux inlet 404) and gas-surface-interaction (503-505) boundaries,
rotational periodicity in 3-D (the 2-D
axisymmetric wedge is supported, with levels too), and
programs that pass a non-empty per-stage external function, boundary types without a multirate
boundary flux, with $L_{\max} > 0$ and the redo of `nonlocal-cfl = redo` (`levels = recompute`, no `fixed-macro-step`) a
cfl above 1 or, in Navier–Stokes runs, a vnn above 1/2
(see [Non-local CFL condition](#non-local-cfl-condition)), and a restart
(`newrun = false`) with `levels = fixed` (the map is not written with the
solution) or with `fixed-macro-step` (a multiple of the smallest local step of
the first macro step of the run).  The boundary types `409`, `501` and `506`,
which MOSE does not apply, stop every run earlier, at the reading of `bc.txt`.
With `levels = fixed` and a `fixed-macro-step` above the step of the map the run goes on with the finest cells beyond
their limit, which it reports (see [Levels and macro step](#levels-and-macro-step)).  With `nonlocal-cfl = redo` and
`fixed-macro-step` or `levels = fixed` the certificate of condition (4) only
counts; with `recompute-every` $> 1$ a map raised by a redo stays in use until
the next rebuild.  It supports the Euler and the Navier–Stokes equations,
laminar or with the Spalart–Allmaras model (`SA`, and `SA-R`) or the SST model, with `k-energy-coupling` and either
form of `sst-production`, on one or several MPI processes, each with OpenMP threads.  The production of `SA-R` and the
dilatation part of the compressible production of SST can be a sink outside the bound of the stiffness that the step
takes ([turbulence sources](turbulence.md#time-accurate-runs-explicit-sources)), as in the single-rate step, and no
multirate test measures their accuracy (the conservation tests run them with levels).  `Wilcox2006`, whose sources the
stages do not evaluate, is rejected, and with $L_{\max} > 0$ so are `SA-RC` and `SST-RC`, whose rotation and curvature
terms are a field of the whole domain that the cell-wise sources of a level do not compute (with $L_{\max} = 0$ the
stages call the block routine of the single-rate stage, which computes them);
`SAcomp` the set-up of MOSE refuses for both integrators.  The keys are listed in
the [Parameter Registry](../user/registry.md#mose-multirate); the
verification tests in [Multirate Verification](../vv/multirate.md).

The stability statements above are proven only in the linear setting of
Theorem 2 of [Appendix A](#appendix-a-order-and-stability-of-the-multirate-macro-step)
(first-order upwind flux, central diffusion) and, for the limited
reconstruction, in the setting of
[Appendix B](#appendix-b-total-variation-at-a-face-between-levels) (the
advection equation in one dimension, two levels, and three levels with a band of one cell, a bound of the total
variation).  The limited reconstruction in several dimensions, with more than three
levels or with three levels in another arrangement, HLLC on the Euler and Navier–Stokes equations, the turbulence sources,
the terms of the axisymmetric wedge (faces 200) and shocks that cross faces
between levels are outside them: for them the stability rests on the tests of
[Multirate Verification](../vv/multirate.md).  Likewise the second order in
time is proven on a fixed mesh (Appendix A, Theorem 1, and for every sequence
of maps its Corollary (a)); with the mesh and the macro step refined together at a fixed Courant number the local error next to a face between levels is of first order (Appendix A, Remark 2), and the error is bounded only for linear two-cell fluxes (Corollary (b): $O(H)$ in the maximum norm and $O(H^2)$ in the volume-weighted $L^1$ norm); outside them it is measured
([Recursive Heun scheme](#recursive-heun-scheme)).

!!! note "Stability of Heun with a limiter"
    Every cell advances within its own local step divided by $\sigma$.  In
    one dimension, for the advection equation with a limited reconstruction,
    the forward-Euler step is a convex combination of the values of the cell
    and of its upwind neighbour, so it keeps a local maximum principle and the
    positivity, for a Courant number up to $1/(1+\theta)$, where $\theta$
    bounds the face increment of the reconstruction relative to the upwind
    difference: $\theta = 1$ for the limiters whose slope can reach twice the
    smaller difference (van Leer, MC, superbee) and $\theta = 1/2$ for minmod.
    On cells of different widths the bounds of the slopes of `Lib_Reconstruction` keep the property on every mesh, for
    $\nu_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, which $\nu_j\,(1+\theta) \le 1$ gives with one time step for all the cells
    (see the interface margin in [Non-local CFL condition](#non-local-cfl-condition)).
    This is the bound $(1+\mu)\,\nu \le 1$ of Hundsdorfer [8] (§6.4, p. 54, with
    $\theta$ the $\mu$ of the limiter condition (5.4) of his §5.2, p. 39),
    stated there for the explicit trapezoidal rule, which is Heun: Heun keeps
    the property of its forward-Euler stages (Hundsdorfer [8], §6.2, Theorem 6.4
    and the example after it, pp. 51–52; Shu and Osher [1], Proposition 2.1 and Remark 2.1).  In several
    dimensions no second-order scheme is TVD
    ([Local time step](#local-time-step)): the forward-Euler stage is a convex
    combination of one-dimensional steps with the Courant number
    $\sum_d \nu_d$ (sum of the directional Courant numbers), so
    $\sum_d \nu_d \le 1/(1+\theta)$ keeps a local maximum principle.  The
    local step adds up the directions, so a cell at its limit has
    $\sum_d \nu_d \le \text{CFL}/\sigma$, 0.4 with $\sigma = 2$ and CFL 0.8, on
    every cell, isotropic ones included.  On a face between levels the coarse
    step measured with the limit of the fine cell can reach
    $2\,\text{CFL}/\sigma$, above $1/(1+\theta)$, without the interface
    margin; with it (always with `nonlocal-cfl = redo`, the default) it stays
    within $1/(1+\theta)$ on the state at the start of the macro step, where the
    levels are built ([Non-local CFL condition](#non-local-cfl-condition)),
    the bound that Hoang et al. [11] (Theorem 4.2, condition (4.15)) assume on
    the first fine cell.  Their theorem, for a one-dimensional scalar law
    discretised by piecewise-linear discontinuous Galerkin elements with the
    second-order strong-stability-preserving Runge–Kutta scheme in two
    regions, with a solution regular at the interface, is not a foundation of
    the stability at the faces between levels of a finite-volume scheme for
    systems, in several dimensions and with several levels: there the
    monotonicity is proven for linear two-cell fluxes (first-order upwind,
    central diffusion), for every cell within 0.86 of its forward-Euler
    positivity step (Appendix A, Theorem 2).  With a limited reconstruction the
    coarse corrector is not a convex combination of forward-Euler steps
    (Appendix A, Remark 4); for the advection equation in one dimension with two
    levels Appendix B bounds the total variation, without the condition (4.15),
    and elsewhere the stability there rests on the tests.

---

## Appendix A: Order and stability of the multirate macro step

This appendix proves the two properties of the macro step that the section
[Recursive Heun scheme](#recursive-heun-scheme) states: the second order on a
fixed mesh (Theorem 1) and the monotonicity for linear two-cell fluxes
(Theorem 2), with their hypotheses and their range, and what they give for a
map that changes from one macro step to the next (Corollary).

**Setting.**  The semi-discretisation is $y' = f(y) = \sum_k f_k(y)$, with $y$
the cell averages of the conservative variables (and of the turbulence
variables) and $k$ running over the faces (convective and diffusive fluxes),
the boundary records and the turbulence sources; $f_k$ depends only on the
cells of its stencil $S_k$.  The level map is nested ($|\ell_c - \ell_v| \le 1$
across every face, equal levels across connections), $h_\ell = H/2^\ell$, and
every stage is $U := U^0 + c_k\,(U - U^0 + h_\ell L(U))$ with $c = (1, \tfrac12)$:
Heun in the form of Shu and Osher [1], $U^{(1)} = U^0 + h\,L(Y^{(1)})$ and
$U^1 = U^0 + \tfrac h2\,[L(Y^{(1)}) + L(Y^{(2)})]$, where $Y^{(1)}$ and $Y^{(2)}$
are the stage vectors, the states that the whole stencil reads.  On a face
between a coarse level $m$ and a fine level $m+1$ the coarse side uses and
saves $G_s$ at its stage 1, the four fine stages accumulate
$A = \tfrac14\sum_{j=1}^4 G_j$, and the coarse side uses $2A - G_s$ at its
stage 2 without evaluating the face
([Conservative correction](#conservative-correction)).

Two facts carry both proofs.

**(F1) Weights and first moments.**  Every stage is linear in the evaluations
of $L$, and the correction replaces a flux by a linear combination of fluxes,
so a macro step updates a cell $c$ as
$U^{n+1}_c = U^n_c + \sum_k \sum_j w_{kj}\, f_{k,c}(Y_{kj})$, with weights
$w_{kj} > 0$ and stage times $t_{kj} \in [0, H]$ counted from $t_n$.  For every
pair of a cell and a face, $\sum_j w_{kj} = H$ and $\sum_j w_{kj}\, t_{kj} = H^2/2$:

- a face without correction, seen by a cell of level $\ell$: the $2^\ell$ Heun
  sub-steps put the weight $h_\ell/2$ at the times $m h_\ell$ and $(m+1) h_\ell$,
  $m = 0, \dots, 2^\ell - 1$, so $\sum w = 2^\ell h_\ell = H$ and
  $\sum w t = \sum_m \tfrac{h_\ell^2}{2}(2m+1) = \tfrac{h_\ell^2}{2}\,4^\ell = \tfrac{H^2}{2}$;
- the coarse side of a corrected face: on each of its sub-steps
  $[\tau, \tau + h]$ it receives $\tfrac h2\,(G_s + 2A - G_s) = hA$, the weight
  $h/4$ at the four fine stage times $\tau$, $\tau + h/2$, $\tau + h/2$ and
  $\tau + h$: $\sum w = h$ and $\sum w t = h\tau + h^2/2$, the moments of a face
  without correction on the same sub-step;
- the fine side of a corrected face: its own four stages, as a face without
  correction.

**(F2) Consistency of the stage vectors.**  Every evaluation reads, on its whole
stencil, $Y_{kj} = y^n + t_{kj}\, f(y^n) + O(H^2)$, with a constant that depends
on the mesh and the map but not on the stage time.  By induction over the
recursion: the cells of the level are read at the start of their sub-step or
at the predictor, $y(\tau) + O(H^2)$ and $y(\tau + h_\ell) + O(H^2)$ (each
sub-step is first-order consistent, and at most $2^L$ of them precede $\tau$);
the finer cells are at the stage time, because the finer levels complete their
sub-steps before the coarser corrector; a coarser cell of level $m$ is read
through the interpolant,
$\hat U_m(t) = U_m(\tau_m) + (t - \tau_m)\, f_m(Y(\tau_m)) = y_m(t) + O(h_m^2)$ for
$t \in [\tau_m, \tau_m + h_m]$, also when it is two or more levels coarser; the
ghost cells are set at the stage time, from boundary data of that time and
from cells for which (F2) holds.

### Theorem 1 (order on a fixed mesh)

*Hypotheses.*  (H1) The mesh and the level map are fixed.  (H2) $f \in C^2$ in a
neighbourhood of the trajectory (no hypothesis for a linear $f$).  (H3) The
recursion, the interpolant and the correction are those above: steps in the
ratio 2 between levels, nesting, equal levels across connections, bands of one
cell allowed, the coarse corrector reading the final values of the finer
cells.  (H4) Time-dependent boundary data are evaluated at the stage times.

*Statement.*  The macro step $\Psi_H$ satisfies $\Psi_H(y) - \varphi_H(y) = O(H^3)$,
with $\varphi_H$ the exact flow of $y' = f(y)$ over $H$ and a constant that
depends on $f$, the mesh and the map, not on where a cell lies with respect to
the faces between levels.  The error of the macro steps up to a time $T$ is
$O(H^2)$, uniformly over the cells, the cells next to the faces between levels
included.

*Proof.*  By (F2), $f_k(Y_{kj}) = f_k(y) + t_{kj}\, f_k'(y)\, f(y) + O(H^2)$, where
the Jacobian $f_k'$ involves only the components of $S_k$, on which (F2)
holds.  With (F1),

$$
\Psi_H(y) - y = \sum_k \Bigl(\sum_j w_{kj}\Bigr) f_k(y)
+ \sum_k \Bigl(\sum_j w_{kj}\, t_{kj}\Bigr) f_k'(y)\, f(y) + O(H^3)
= H f + \frac{H^2}{2}\, f' f + O(H^3),
$$

the expansion of $\varphi_H(y) - y$; for a linear $f = Ly$ the same count gives
$\Psi_H = I + HL + \tfrac{H^2}{2}L^2 + O(H^3)$.  The macro step is a composition
of $2^{L+2} - 2$ maps $Y \mapsto Y + h\,L(Y)$ and of linear combinations with
bounded coefficients, so it is Lipschitz with a constant $1 + CH$, $C$
proportional to $\|f'\|$.  Then $|e_{n+1}| \le (1 + CH)\,|e_n| + DH^3$ and
$|e_n| \le DH^2\,(e^{CT} - 1)/C$. ∎

*Remark 1 (the two differences from [10]).*  (i) The coarse corrector reads the
finer cells at their final values, $y_F(\tau + h) + O(h^3)$, where the scheme of
[10] reads their predictor, $y_F(\tau + h) + O(h^2)$: both satisfy (F2), and the
choice enters only the $O(H^3)$ term.  With two-cell fluxes the two choices
give the same scheme, since the face between levels is not evaluated again at
the coarse stage 2.  (ii) With several levels and bands of one cell, (F1) and
(F2) hold level by level, and a halo cell two or more levels coarser changes
only the constant, through $O(h_m^2)$.  Neither difference adds an order
condition.

*Remark 2 (range).*  The constant of Theorem 1 contains $\|f'\| \sim 1/\Delta x$:
the theorem holds for $H \to 0$ on a fixed mesh, which is what `mr-order`,
`mr-deep` and `mr-bc` measure.  Where the limiter switches, and at the
switches of HLLC and of the turbulence models, (H2) fails locally, and the
local error of Heun itself is $O(H^2)$ there.  At a fixed Courant number
($H \propto \Delta x$) the theorem gives no bound (the Corollary below gives
one for linear two-cell fluxes).  At the coarse stage 2 the
coarse cells are at their forward-Euler predictor, with the error
$-\tfrac{h^2}{2}L^2 y$, and the finer cells at their final values, with an
error $O(h^3)$.  A flux that reads both sides of a face between levels (the
MUSCL slope of the coarse cell next to the face, the HLLC flux, the diffusive
flux) sees a jump $O(h^2)$ in the errors of its stencil, and its error,
$O(h^2)$, differs by as much from that of the next face: the error of the
state is $h \cdot O(h^2)/\Delta x = O(h^2)$ per macro step in the two cells of
the face, with opposite signs.  The local error is then $O(H^2)$ next to a face
between levels and $O(H^3)$ elsewhere: first order locally, as the
conventional boundary data of [7].

### Theorem 2 (monotonicity for linear two-cell fluxes)

*Hypotheses.*  $dU_i/dt = \sum_j \beta_{ij}\,(U_j - U_i)$ with $\beta_{ij} \ge 0$, the
flux of every face being a linear function of its two cells alone: the
first-order upwind flux ($\beta_{ij} = |a_f| A_f/V_i$ for the upwind neighbour,
with $a_f$ the normal velocity through the face of area $A_f$, and 0 for the
downwind one) and the central diffusion ($\beta_{ij} = D A_f/(\delta_f V_i)$, with
$\delta_f$ the distance between the two cell centres), on any mesh and in any
dimension (the diffusive flux of MOSE is this two-cell flux for a scalar on orthogonal cells; on skewed cells it also
reads the neighbours along the face, `Lib_Diffusive`).  With the upwind flux the velocity field has no discrete divergence
($\sum_f a_f A_f = 0$ over the faces of every cell, $a_f$ outward); this form
is then the conservative scheme of the face fluxes $|a_f| A_f U_j$, also in the
macro step, where every cell sends its outflow at its own stage values (on a
corrected face too: $2A - G_s = |a_f| A_f\,U_i^{(1)}$ when the coarse cell $i$
is upwind), and the conservation in the statement uses it.  The load of cell
$i$ with its step $h_i = H/2^{\ell_i}$ is
$x_i = h_i \sum_j \beta_{ij}$; the forward-Euler step
$(1 - x_i)\,U_i + h_i\sum_j \beta_{ij}U_j$ is positive for $x_i \le 1$, and on a
Cartesian mesh $x_i = h_i\sum_d\bigl(|a_d|/\Delta x_d + 2D/\Delta x_d^2\bigr)$,
the ratio of $h_i$ to the local step with CFL = 1 and VNN = 1/2
([Local time step](#local-time-step)).  Any nested map, any number of levels,
and:

- (C1) $x_i \le 1$ for every cell;
- (C2) $x_i \le 2\,c(x_j)$, with $c(x) = 1 - x + x^2/2 - x^3/8$, on every face
  where a finer cell $j$ feeds a coarser cell $i$ ($\beta_{ij} > 0$,
  $\ell_j = \ell_i + 1$).

*Statement.*  The matrix $M$ of the macro step has non-negative entries,
$M\mathbf 1 = \mathbf 1$ and $\mathbf 1^T V M = \mathbf 1^T V$, with $V$ the diagonal
of the cell volumes.  Hence $\|M\|_\infty = 1$ (maximum principle),
$\|M\|_{1,V} = 1$ (conservation and positivity) and $\|M\|_{2,V} \le 1$, and the
same holds for every product of macro steps, with different maps.

*Proof.*  By induction over the recursion, every stage vector, interpolant and
end-of-sub-step state is a combination of the values at $t_n$ with
non-negative coefficients that sum to 1; such combinations are preserved by
composition, and the only coefficient that can become negative arises in the
coarse corrector.  Take a cell $i$ on a sub-step $[\tau, \tau + h]$, $h = h_i$.
Stage 1 gives $U_i^{(1)} = (1 - x_i)\,U_i^0 + h\sum_j \beta_{ij}\,\hat U_j(\tau)$, a
convex combination by (C1), $\hat U_j$ being the neighbour at the stage time
(its state, or its interpolant, itself a convex combination).  On a corrected
face with a finer neighbour $j$ the four fine stages read the coarse cell at
$\eta = 0, \tfrac12, \tfrac12, 1$, whose mean is $(U_i^0 + U_i^{(1)})/2$, so
$2A - G_s = \beta_{ij}\,[\,2\bar U_j - U_j(\tau) - U_i^{(1)}\,]$, with $\bar U_j$ the
mean of the four fine stage vectors.  Stage 2 then gives

$$
U_i^1 = \Bigl(1 - \frac{x_i}{2}\Bigr)U_i^0 - \frac{x_i}{2}\,U_i^{(1)}
+ \frac h2 \sum_{j\ \text{not finer}} \beta_{ij}\,\bigl[\hat U_j(\tau) + \hat U_j(\tau + h)\bigr]
+ h \sum_{j\ \text{finer}} \beta_{ij}\,\bar U_j .
$$

With $U_i^{(1)}$ substituted, $U_i^0$ has the coefficient $1 - x_i + x_i^2/2 > 0$;
a neighbour that is not finer contributes
$\tfrac h2\,\beta_{ij}(1 - x_i)\,\hat U_j(\tau) + \tfrac h2\,\beta_{ij}\,\hat U_j(\tau + h)$,
non-negative by (C1); a finer neighbour contributes $h\beta_{ij}\,(c_j - x_i/2)$
on $U_j(\tau)$, with $c_j$ the coefficient of $U_j(\tau)$ in $\bar U_j$, plus
non-negative terms.  The four fine stage vectors carry $U_j(\tau)$ with the
coefficients 1, $1 - x_j$ (the predictor), at least $1 - x_j + x_j^2/2$ (the
state after the first fine sub-step: the diagonal of Heun,
$I + hL + h^2L^2/2$, plus non-negative products) and at least
$(1 - x_j)(1 - x_j + x_j^2/2)$ (the predictor of the second fine sub-step), so
$c_j \ge c(x_j)$, and (C2) makes the coefficient non-negative.  A cell of a
band of one cell is coarse on one side and fine on the other, and falls in the
cases above.  The rows sum to 1 because a constant is an exact solution (every
flux of a constant vanishes, and with it $A$, $G_s$ and the correction); the
columns weighted by the volumes sum to 1 because every face gives its two
sides the same impulse with opposite signs (in the conservative form of the
hypotheses), also through the correction ($h_m A$ on either side).
$M \ge 0$ with $M\mathbf 1 = \mathbf 1$ gives
$\|M\|_\infty = 1$; with $\mathbf 1^T V M = \mathbf 1^T V$ it gives
$\|M\|_{1,V} = 1$; and since the rows are convex,
$(MU)_i^2 \le \sum_j M_{ij}U_j^2$ (Jensen), so
$\sum_i V_i\,(MU)_i^2 \le \sum_j V_j U_j^2$. ∎

*Remark 3 (threshold).*  $c$ decreases on $[0, 1]$, with $c(1) = 3/8$, and
$2c(x) \ge x$ exactly when $x^3/4 - x^2 + 3x - 2 \le 0$, that is for
$x \le x^* = 0.8603$: $x_i \le 0.86$ for every cell implies (C1) and (C2).  The
threshold is exact for the pure upwind flux with the finer cell upstream.
With the additive local-step rule and $h_c \le \Delta t_c/\sigma$, writing
$r_c = \sum_d |a_d|/\Delta x_d$ and $r_v = \sum_d D/\Delta x_d^2$,

$$
x = h\,(r_c + 2r_v) \le \frac{r_c + 2r_v}{\sigma\,(r_c/\text{cfl} + r_v/\text{vnn})}
\le \frac{\max(\text{cfl},\, 2\,\text{vnn})}{\sigma},
$$

0.4 with $\sigma = 2$, cfl 0.8 and vnn 0.4.  A single-rate Heun step is
positive up to $x = 1$ (Hundsdorfer [8], §6.1, pp. 49–50): the faces between
levels cost at most 14 % of the positivity step.  On a mesh of equal rectangular cells with constant coefficients, with
cfl $\le 1$ and, in Navier–Stokes
runs, vnn $\le 1/2$, the bounds that the set-up enforces when $L_{\max} > 0$ with the redo of `nonlocal-cfl = redo`
([Non-local CFL condition](#non-local-cfl-condition)), and the default $\sigma = 2$, every cell starts the macro step with
$x \le 1/2$ (in Euler runs $r_v = 0$ and $x \le \text{cfl}/\sigma$), below $x^*$: (C1) and (C2) hold for every map that
the level rule builds, whatever the limiter of the run (the loads are those of the first-order fluxes of Theorem 2), for
every $\sigma \ge 2$, with or without the interface margin; with a smaller
$\sigma$ they hold while $\max(\text{cfl}, 2\,\text{vnn})/\sigma \le x^*$.  Within the macro
step the certificate of (4) (with `redo`) accepts, on the state at every start $\tau$ of a sub-step, a cell with $x \le 1$, the load of
(C1), and the finer
cell of a face between levels, at a start that it shares with the coarser one ($\hat h = 2h$), with $x \le 1/2$, where
$2c(1/2) = 1.22 > 1$, the load of (C2).  Theorem 2 takes the coefficients fixed within the macro step; when they change,
these bounds at the starts of the sub-steps are what the certificate checks, not a theorem.  The bound
$x \le \max(\text{cfl}, 2\,\text{vnn})/\sigma$ is written for a mesh of equal rectangular cells with constant
coefficients; on other rectangular cells the load of the diffusion can exceed $2r_v$
([Stability Conditions](#stability-conditions)), and on skewed cells the diffusion of MOSE is not a two-cell flux and
is outside Theorem 2.

*Remark 4 (what Theorem 2 does not cover).*  With a limited reconstruction the
coarse corrector is
$U_i^1 = \tfrac12 U_i^0 + \tfrac12\,[\,U_i^{(1)} + h\,\tilde L_i(Y^{(2)})\,]$, where
$\tilde L$ has, on a face between levels, the extrapolation $2A - G_s$ (twice
the mean of the fine fluxes minus the flux at $\tau$) in place of a flux: it
is not a convex combination of forward-Euler steps.  For a monotone flux
$G = aU$, $2A - G_s = a\,(2\bar U_j - U_j(\tau))$ leaves the range of the values
when the fine value decreases fast.  In the linear two-cell case this is the
term $-\tfrac{x_i}{2}\,h\beta_{ij}\,U_j(\tau)$ that $c_j$ absorbs above; with a
limiter the corresponding coefficient depends on the limiter at different
times, and the argument by convex combinations, by which Heun inherits the
properties of its forward-Euler stages, does not go through at the faces
between levels.  Outside Theorem 2 are also HLLC on the Euler and
Navier–Stokes equations, the turbulence sources, the terms of the
axisymmetric wedge (faces 200, which are not two-cell fluxes) and shocks that
cross faces between levels.

### Corollary (every sequence of maps)

*Hypotheses.*  The macro steps $n = 0, \dots, N-1$ take the steps $H_n \le H$,
with $\sum_n H_n = T$, and nested maps $M_n$ with levels in
$\{0, \dots, L\}$, $L$ fixed.  $M_n$ and $H_n$ may depend on the numerical
solution $y_n$: the map rebuilt at every macro step, the levels raised by a
redo, a step shortened to land on an output instant.  The error
$e_n = y_n - y(t_n)$ is taken against the exact solution of $y' = f(y)$ on the
same mesh, as in Theorem 1.

*Statement.*  (a) On a fixed mesh, under (H2)–(H4) and for $H$ small enough,
$|e_N| \le DH^2\,(e^{CT} - 1)/C$, with $D$ and $C$ the largest constants of
Theorem 1 over the nested maps.  (b) In the setting of Theorem 2, with every
$M_n$ within (C1) and (C2) for its step $H_n$ (as the maps of the level rule
are, on a mesh of equal rectangular cells with constant coefficients, with cfl $\le 1$, vnn $\le 1/2$ and
$\sigma \ge 2$, Remark 3), let the
second and third time derivatives of the exact
solution be bounded on $[0, T]$ by $K_2$ and $K_3$ independently of the mesh
(since $y^{(k)}(t) = e^{tL} L^k y(0)$ and $e^{tL}$ is non-negative with unit
row sums, bounds on $L^2 y(0)$ and $L^3 y(0)$ suffice), and let $|B|$ bound the
volume of the cells within $2^{L+2}$ cells of a face between levels of every
$M_n$.  Then

$$
\|e_N\|_\infty \le c_L K_2\, T H, \qquad
\|e_N\|_{1,V} \le \Bigl(c_L K_2\,\frac{|B|}{H} + \frac{K_3\,|\Omega|}{6}\Bigr)\, T H^2,
$$

with $c_L$ that depends only on $L$ and $|\Omega|$ the volume of the domain.
At a fixed Courant number ($H \propto \Delta x$), with faces between levels of
bounded total area on a quasi-uniform family of meshes, $|B| = O(\Delta x)$:
the error is $O(H)$ in the maximum norm and $O(H^2)$ in the volume-weighted
$L^1$ norm, for every sequence of maps, a face between levels that travels
with a wave included.

*Proof.*  In both parts
$e_{n+1} = [\Psi_{H_n}^{M_n}(y_n) - \Psi_{H_n}^{M_n}(y(t_n))] + \tau_n$, with
$\tau_n = \Psi_{H_n}^{M_n}(y(t_n)) - y(t_{n+1})$ the local error of the map
$M_n$ on the exact solution: every bound below holds for every admissible map,
so the choice of $M_n$ from the numerical solution does not enter.  (a) On a
fixed mesh the nested maps are finitely many, so the largest constants exist,
$|e_{n+1}| \le (1 + CH_n)\,|e_n| + DH_n^3$, and
$|e_N| \le DH^2 \sum_n H_n\, e^{C(T - t_{n+1})} \le DH^2\,(e^{CT} - 1)/C$; for
$H$ small enough, by induction, the numerical solution stays in the
neighbourhood of (H2).  (b) The macro step is linear, $y_{n+1} = M_n y_n$, so
$e_N = \sum_n M_{N-1} \cdots M_{n+1}\, \tau_n$, and every product has norm 1 in
$\|\cdot\|_\infty$ and in $\|\cdot\|_{1,V}$ (Theorem 2).  Take every value that
the macro step computes (start state, predictor, interpolant, end of a
sub-step) with its deviation from the exact solution at its time.  A predictor
is a convex combination, by (C1), of values read at $\tau$, and adds the
remainder $y_i(\tau + h) - y_i(\tau) - h\,y_i'(\tau)$, at most $h^2 K_2/2$; an
interpolant adds at most $h_m^2 K_2/8$; a corrector adds at most twice the
largest deviation of the values it reads (the weights of each face sum to $h$,
and $h \sum_j \beta_{ij} \le 1$) and a quadrature error at most $h^2 K_2$ (the
weights of each face have the moments of (F1), and
$|\beta_{ij}\,(y_j'' - y_i'')| \le 2\beta_{ij} K_2$).  Over the at most
$2^{L+2} - 2$ stages of a macro step every deviation, $\tau_n$ included, stays
within $c_L K_2 H_n^2$.  The value of a cell after the macro step depends only
on the cells within $2^{L+2}$ cells of it (each stage reads face neighbours);
if they all have its level $\ell$, it is the value after $2^\ell$ single-rate
Heun steps of $h_\ell$, whose rows there are non-negative with unit sums
(by (C1), since $I + hL + \tfrac{h^2}{2}L^2 = \tfrac12 I + \tfrac12 (I + hL)^2$),
each with the local error $y(\tau + h) - (I + hL + \tfrac{h^2}{2}L^2)\,y(\tau)$,
at most $h^3 K_3/6$: $|\tau_{n,i}| \le 2^\ell h_\ell^3 K_3/6 \le K_3 H_n^3/6$.
With $\sum_n H_n^2 \le HT$ and $\sum_n H_n^3 \le H^2 T$ the bounds follow. ∎

Part (b) holds in the setting of Theorem 2 only: outside it (Remark 4) the
bound of the products of macro steps is missing, and at a fixed Courant number
no bound is proven.  Part (a) holds on a fixed mesh, with constants that grow
with $1/\Delta x$ (Remark 2).  In the maximum norm the order 1 is what the
errors of order $H^2$ of the $T/H$ macro steps give when they add up on the
same cells, as next to a face between levels that travels with a wave at its
speed.

### Proposition (the margin on condition (4))

*Setting.*  $\Delta t_c$ is the admissible step of cell $c$ for the attempt: the local step at $t_n$, a positive finite
number ([Local steps](#levels-and-macro-step)), lowered by a redo to the smallest value that the cell showed;
$\Delta t_c(\tau)$ its value at the start $\tau$ of a sub-step, as the certificate evaluates it.  $v \sim c$ runs over the
face neighbours of $c$ in its block, every one of them on the process of $c$.  $\varepsilon = 10^{-12}$ and

$$
H_{\text{cap}} = \min_c\ 2^{\ell_c} \min\Bigl(\frac{\Delta t_c}{\sigma},\ \min_{v \sim c} \frac{\Delta t_v}{\kappa}\Bigr)
$$

over all the processes, for a nested map with equal levels across connections.

*Hypotheses.*  (M1) $\sigma \ge 2$, $L_{\max} > 0$, the margin on (`nonlocal-cfl = redo` or `interface-margin`).
(M2) $\kappa = \max(\text{cfl}\,(1+\theta), \kappa_0)$ with $1 < \kappa_0 \le 2$ and $\text{cfl}\,(1+\theta) \le \sigma$ (from
cfl $\le 1$, which the set-up and the readings check with the redo of `nonlocal-cfl = redo`, with $\theta \le 1$ and $\sigma \ge 2$), so $1 < \kappa \le \sigma$.  (M3) $H \le (1 + \varepsilon)\, H_{\text{cap}}$: true of every
attempt without `fixed-macro-step`, since the redo halves $H$ while it is above that bound of the map it has rebuilt, a target and
a halving only shorten $H$, and a landing lengthens it by at most $\varepsilon H$ plus half a unit in the last place of the time
(a relative $1.1\cdot10^{-16}\,t/H$, above $\varepsilon$ when $H < 10^{-4}\,t$); (a) then holds with that larger $\varepsilon$, without
consequence next to $1 - 1/\kappa_0$.  No hypothesis on the mesh, the state, the limiter, the fluxes, the number of dimensions, of blocks or of
processes.

*Statement.*

(a) If $\Delta t_v(\tau) \ge \Delta t_v/\kappa$ for every cell $v$ and every start $\tau$ of its sub-steps in the attempt,
the certificate finds no violation.  So every rejection of a violation of (4) comes from a cell whose admissible step fell by more
than $1 - 1/\kappa \ge 1 - 1/\kappa_0$ within the macro step.

(b) The margin raises a cell $c$ only where a face neighbour has $\Delta t_c/\Delta t_v > \sigma/\kappa$.  Where no face
of the mesh has that jump, the levels and $H$ are those of the rule without the margin, to the bit.

(c) Before nesting and connections the level of a cell is at most the largest level that its own term gives to the cell or to
one of its face neighbours.  So the margin does not propagate through the rule, and the deepest level does not change.

(d) In the adaptive step without a target, on a map built for $H_{\text{try}} = 2^{L}\Delta t_{\min}/\sigma$,
$2^{-\varepsilon} H_{\text{try}} \le H_{\text{cap}} \le H_{\text{try}}$: the first $H$ does not depend on $\kappa$, to a
relative $7\cdot10^{-13}$.

*Proof.*  (a) Take $v$ and a start $\tau$ of a sub-step.  If a face neighbour $g$ of $v$ in the block is one level coarser and
starts a sub-step at $\tau$, then $\hat h_v = h_{\ell_g} = H/2^{\ell_g}$, and $v \sim g$ gives in $H_{\text{cap}}$ the
term $2^{\ell_g} \Delta t_v/\kappa$, so $\hat h_v \le (1 + \varepsilon) \Delta t_v/\kappa$ by (M3).  Otherwise $\hat h_v =
h_{\ell_v} \le (1 + \varepsilon) \Delta t_v/\sigma \le (1 + \varepsilon) \Delta t_v/\kappa$ by (M2).  A neighbour across a connection has
the level of $v$, so it is never the coarser one, and a finer neighbour does not enter $\hat h_v$.  By hypothesis the
last member is at most $(1 + \varepsilon) \Delta t_v(\tau)$, the test of the certificate.  (b) The rule takes for $c$ the larger
of $\sigma H/\Delta t_c$ and $\kappa H/\Delta t_v$, $v \sim c$; the second exceeds the first only if
$\Delta t_c/\Delta t_v > \sigma/\kappa$.  Without that jump every cell has the first in the rule, and in
$H_{\text{cap}}$ the term of $v$ is not below that of $c$.  (c) $\kappa H/\Delta t_v \le \sigma H/\Delta t_v$ by (M2), and the
formula of the levels is monotone: the margin term of $c$ asks at most the level of the own term of $v$; it reads the local
steps of the neighbours, not their levels.  (d) By (c) the cell with $\Delta t_{\min}$ keeps level $L$, the rules cannot
take it higher (their fixed point does not exceed the largest level of the rule), and its term is $H_{\text{try}}$, so
$H_{\text{cap}} \le H_{\text{try}}$.  Every level satisfies $\log_2 x \le \ell_c + \varepsilon$ for both terms $x$ of the rule, so every
term of $H_{\text{cap}}$ is at least $2^{-\varepsilon} H_{\text{try}}$; no level is cut at $L_{\max}$, since $x \le 2^{L}$ by
(c).  ∎

*The extremes.*  cfl $\to 0$ gives $\kappa = \kappa_0$, the same relative fall.  At $\text{cfl}\,(1+\theta) = \kappa_0$ the two
terms of $\kappa$ meet, and at $\text{cfl}\,(1+\theta) = \sigma$, which (M1) and (M2) allow only with $\sigma = 2$, $\theta = 1$ and
cfl 1, the margin is the term of the cell itself.  With $\sigma > 2$ a coarse neighbour of a cell at the limit of its level is raised only if $\kappa > \sigma/2$;
(a) holds as it stands.  With $\sigma < \kappa$, outside (M1) and (M2), a cell at the limit of its own term has, at a start
that no coarser neighbour shares, $\hat h_v = \Delta t_v/\sigma$, and a fall of its admissible step by more than
$1 - 1/\sigma$ rejects the attempt: with $\sigma = 1$, every fall.  A cell with no neighbour in its block has no margin term.  At $L_{\max}$ a level can be cut, after a redo,
and (a) holds in that case too, because it rests on the term of $H_{\text{cap}}$ and not on the level.  With $L_{\max} = 0$ the margin
is off and there is no face between levels.  At $\kappa = 1$, $\sigma = 2$ (a) leaves no room: a cell at the limit of its level with a
coarser neighbour has $\hat h_v = \Delta t_v$.

*What it does not say.*  It does not say that (4) is sufficient (the source gives no proof either [13]); that a front is
not rejected (a fall by more than $1 - 1/\kappa$ is necessary for a rejection, not sufficient; where it happens, (4) asks for it); how much the raised cells cost, which depends on how many
faces have the jump of (b); or why $\kappa_0$ is 1.1, a value taken from another code (section
[Non-local CFL condition](#non-local-cfl-condition)).  With `fixed-macro-step` (M3) fails, and there the certificate only counts.
A map with the margin is a nested map, so Theorem 1, Theorem 2, the Corollary and Appendix B hold for it as for any other.

## Appendix B: Total variation at a face between levels

This appendix proves the bound of the total variation that the section
[Non-local CFL condition](#non-local-cfl-condition) states for the limited reconstruction, in a restricted setting: the
advection equation in one dimension, two levels and one face between them (Theorem 3), and three levels with a band of one
cell (Theorem 3b), with the reconstruction of `Lib_Reconstruction` on cells of different widths.  It is a bound, not the TVD
property: one macro step can raise the total variation, by terms local to the faces between levels (Remark 6).

**Setting.**

- (B1) $u_t + a\,u_x = 0$, $a > 0$ constant, on cells $I_j$ of widths $\Delta x_j$, with cell averages $u_j$ and a constant
  state on each side far from the face between levels; $d_j = u_j - u_{j-1}$ is the difference across the face $j - \tfrac12$
  and $TV(u) = \sum_j |d_j|$.  The weights $w^+_j = \Delta x_j/(\Delta x_j + \Delta x_{j+1})$ and
  $w^-_j = \Delta x_j/(\Delta x_{j-1} + \Delta x_j)$ satisfy $w^+_{j-1} + w^-_j = 1$ and are $\tfrac12$ on cells of equal
  width.  The face between levels lies between the cells 0 and 1, the cell 1 downstream; $\rho = \Delta x_0/\Delta x_1$ is the
  ratio of the upstream to the downstream width there, and $r$ the ratio of the coarse to the fine width ($\rho = r$ with the
  coarse cell upstream, $\rho = 1/r$ with the fine cell upstream).
- (B2) The reconstruction of `Lib_Reconstruction`, which `Lib_Convective` calls with the distances between the cell centres
  $l^0_j = (\Delta x_{j-1} + \Delta x_j)/2$ and $l^1_j = (\Delta x_j + \Delta x_{j+1})/2$ and the half width $\Delta x_j/2$, and the
  upwind flux $a\,q_{j+1/2}$: $q_{j+1/2} = u_j + s_j$ with $s_j = \hat\phi_j\,w^+_j\,\Delta_+u_j$,
  $r'_j = (\Delta_-u_j/\Delta_+u_j)\,(l^1_j/l^0_j)$, $\hat\phi_j = \min(\phi(r'_j),\ 1/w^+_j,\ r'_j/w^-_j)$ for $r'_j > 0$ and
  $\hat\phi_j = 0$ for $r'_j \le 0$, and $s_j = 0$ when $\Delta_+u_j = 0$.  The limiter receives the bounds
  $l^1_j/(\Delta x_j/2) = 1/w^+_j$ and $l^0_j/(\Delta x_j/2) = 1/w^-_j$, which keep the face value between $u_j$ and $u_{j+1}$;
  with minmod and the first-order reconstruction they never act.  The term $10^{-30}$ that `Lib_Limiters` adds to the
  denominator of the ratio of the slopes is left out.  On cells of equal width the bounds are 2 and $2r'_j$, which no limiter of
  (B3) exceeds, and this is the form of the interface margin above.
- (B3) The limiter: $\phi(r) = 0$ for $r \le 0$, and $0 \le \phi(r) \le 2\theta$, $0 \le \phi(r)/r \le 2\theta$ for $r > 0$,
  with the $\theta$ of the interface margin: 0 with the first-order reconstruction, $\tfrac12$ for minmod, 1 for van Leer, MC
  and superbee.  Van Albada, whose $\phi$ is negative for $-1 < r < 0$, is outside.
- (B4) Two levels with the macro step of MOSE: the coarse cells advance with $H$ and the fine cells with $h = H/2$, with the
  Courant numbers $\lambda_j = aH/\Delta x_j$ and $\mu_j = ah/\Delta x_j$.  The coarse stage 1 at $t_n$; the two Heun
  sub-steps of the fine cells, whose four stages read the coarse cells through the interpolant at the fractions
  $0, \tfrac12, \tfrac12, 1$ of $H$; the coarse stage 2, which reads the fine cells at their final values.  On the face between
  levels the coarse cell uses and saves its flux $G_s$ at its stage 1, each fine stage adds a quarter of its flux to $A$, and
  the coarse cell uses $2A - G_s$ at its stage 2 ([Conservative correction](#conservative-correction)).
- (B5) Every cell satisfies $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, with $c_j$ its own Courant number ($\lambda_j$ or
  $\mu_j$); on cells of equal width this is $c_j\,(1+\theta) \le 1$ (Remark 5).

A forward-Euler stage on a state $w$ (the values that the stage reads on its whole stencil) updates the cells of its level as
$w'_j = w_j - c_j\,(q_{j+1/2} - q_{j-1/2})$, with the face values of $w$.

**Lemma 1 (incremental form of a stage).**  For every state $w$ and cell $j$ there are
$\alpha_j \in [0, \min(2\theta w^+_j, 1)]$ and $\beta_j \in [0, \min(2\theta w^-_j, 1)]$ with
$s_j = \alpha_j\,(w_{j+1} - w_j) = \beta_j\,(w_j - w_{j-1})$, so that

$$
q_{j+1/2} - q_{j-1/2} = B_j\,(w_j - w_{j-1}), \qquad
B_j = 1 + \beta_j - \alpha_{j-1} \in [\,1 - \min(2\theta w^+_{j-1}, 1),\ 1 + \min(2\theta w^-_j, 1)\,] \subset [0, 2].
$$

*Proof.*  If $\Delta_+w_j = 0$ or $r'_j \le 0$, then $s_j = 0$: take $\alpha_j = \beta_j = 0$.  Otherwise take
$\alpha_j = w^+_j\,\hat\phi_j$ and, since $(l^1_j/l^0_j)\,w^+_j = w^-_j$, $\beta_j = w^-_j\,\hat\phi_j/r'_j$.  By (B2) and (B3),
$0 \le \hat\phi_j \le \phi(r'_j) \le 2\theta$ and $0 \le \hat\phi_j/r'_j \le \phi(r'_j)/r'_j \le 2\theta$, and by (B2)
$\hat\phi_j \le 1/w^+_j$ and $\hat\phi_j/r'_j \le 1/w^-_j$: so $\alpha_j \le \min(2\theta w^+_j, 1)$ and
$\beta_j \le \min(2\theta w^-_j, 1)$.  Then $q_{j+1/2} - q_{j-1/2} = (w_j - w_{j-1}) + s_j - s_{j-1}$, with
$s_j = \beta_j\,(w_j - w_{j-1})$ and $s_{j-1} = \alpha_{j-1}\,(w_j - w_{j-1})$, and $\alpha_{j-1} \le 1$, $\beta_j \le 1$ give
$0 \le B_j \le 2$. ∎

So $B_j \ge 0$ at every face, the face between levels included, on cells of any widths.  The bound of $\beta_j$ does not depend
on the width of the downstream cell: the slope of a cell that leaves through a face is written with its upstream difference,
whatever lies downstream, and the two sides of the face between levels can be treated apart.

**Lemma 2 (a region with a ghost upstream).**  Let $J = \{j_a, j_a + 1, \dots\}$ be consecutive cells that end in a constant
state or in an outflow cell $j_b$ (one whose slope reads a cell outside $J$), $g = w_{j_a - 1}$ the value upstream of $J$ (the
ghost), and $TV_J(w; g) = |w_{j_a} - g| + \sum_{j \in J,\ j > j_a} |w_j - w_{j-1}|$.  If (B5) holds in $J$, a forward-Euler stage on
$J$ with the ghost fixed gives

$$
TV_J(w'; g) \le TV_J(w; g) - c_{j_b}B_{j_b}\,|w_{j_b} - w_{j_b - 1}|,
$$

without the last term when $J$ ends in a constant state.

*Proof.*  By Lemma 1, $w'_{j_a} - g = (1 - c_{j_a}B_{j_a})(w_{j_a} - g)$ and
$d'_j = (1 - c_jB_j)\,d_j + c_{j-1}B_{j-1}\,d_{j-1}$ for $j > j_a$.  By Lemma 1 and (B5),
$0 \le c_jB_j \le c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$ in $J$.  Taking absolute values and summing, each
$c_{j-1}B_{j-1}|d_{j-1}|$ cancels the $-c_{j-1}B_{j-1}|d_{j-1}|$ of the row before; the outflow term remains. ∎

**Lemma 3 (a Heun sub-step).**  (a) Let $\tilde w = \tfrac12 w + \tfrac12 E_2(E_1(w))$, with $E_k$ a forward-Euler stage on $J$
with the ghost $g_k$, under the hypotheses of Lemma 2.  Then
$TV_J(\tilde w; g_2) \le TV_J(w; g_1) + |g_2 - g_1|$.
(b) If $J$ starts at a constant state and ends in an outflow cell $j_b$, then
$TV_J(\tilde w) \le TV_J(w) - D$, with
$D = \tfrac12\,c_{j_b}\,(B^{(1)}_{j_b}|d^{(1)}_{j_b}| + B^{(2)}_{j_b}|d^{(2)}_{j_b}|) \ge |\tilde w_{j_b} - w_{j_b}|$, the
superscripts naming the two stages.

*Proof.*  $TV_J(\cdot\,; g)$ is convex.  (a) $TV_J(\tilde w; g_2) \le \tfrac12\,TV_J(w; g_2) + \tfrac12\,TV_J(E_2E_1w; g_2)$, with
$TV_J(w; g_2) \le TV_J(w; g_1) + |g_2 - g_1|$ and, by Lemma 2 twice,
$TV_J(E_2E_1w; g_2) \le TV_J(E_1w; g_1) + |g_2 - g_1| \le TV_J(w; g_1) + |g_2 - g_1|$.  (b) Lemma 2 twice, and
the convexity; $\tilde w_{j_b} - w_{j_b} = -\tfrac12\,c_{j_b}(B^{(1)}_{j_b}d^{(1)}_{j_b} + B^{(2)}_{j_b}d^{(2)}_{j_b})$. ∎

### Theorem 3 (total variation at a face between levels)

*Statement.*  Under (B1)–(B5), one macro step with one face between levels gives

$$
TV(u^{n+1}) \le TV(u^n) + R,
$$

with $R$ as follows.

- *Coarse cell upstream* (cells $j \le 0$ coarse, $j \ge 1$ fine).  Let $p_0$ be the predictor of the cell 0,
  $d^p_0 = p_0 - p_{-1}$, $X = B^0_0\,d_0$ and $Y = B^1_0\,d^p_0$ the increments of Lemma 1 of the cell 0 at its two stages
  without the correction, and $\delta_0$ the change that the correction makes to the cell 0.  Then

  $$
  R = \tfrac12\lambda_0\,(|X| - |Y| + |X - Y|) + 2|\delta_0|
  \le \lambda_0\,\bigl[(1 + \tfrac32\theta)\,|d_0| + \theta\,|d^p_0|\bigr]
  \le \lambda_0\,\bigl[(1 + \tfrac52\theta)\,|d_0| + \theta\,|d_{-1}|\bigr],
  $$

- *Fine cell upstream* (cells $j \le 0$ fine, $j \ge 1$ coarse).  Let $\Phi = 2A - G_s$ and $\Phi_H$ be the face values that
  the cell 1 uses at its stage 2 with and without the correction (without it, from the fine cells at their final values),
  $d^{(s)}_0$ the difference of the last fine cell at the fine stage $s$ ($d^{(1)}_0 = d_0$) and $d^{n+1}_0$ at the end, and
  $m = \mu_0\,(1 + \theta) \le 1$.  Then

  $$
  R = \lambda_1\,|\Phi - \Phi_H| \le \lambda_1\,\bigl[(m + 3\theta)\,|d_0| + (m + 4\theta)\,|d_{-1}|
  + (\tfrac34 m + 3\theta)\,|d_{-2}| + (\tfrac14 m + \tfrac54\theta)\,|d_{-3}| + \tfrac14\theta\,|d_{-4}|\bigr],
  $$

The bounds of $R$ take the two cells upstream of the face of equal width (otherwise $\theta$ becomes $\min(2\theta w^-_0, 1)$ in
them).

*Proof, coarse cell upstream.*

1. *Regions.*  The coarse region $C = \{\dots, -1, 0\}$ reads the fine cells only in the slope $s_0$, which leaves through the
   face $\tfrac12$: written as $s_0 = \beta_0\,(w_0 - w_{-1})$ with $\beta_0 \in [0, \min(2\theta w^-_0, 1)]$, whatever the
   fine cell, it makes $C$ a region with the outflow cell 0, and $B_j \ge 0$ in $C$.  The fine region $F = \{1, 2, \dots\}$ reads
   the cells $-1$ and 0 through the interpolant: its ghost is $g(\tau) = u_0 + \tau\,(p_0 - u_0)$ at the fraction $\tau$ of $H$,
   and $w_{-1}$ enters only $\alpha_0 \in [0, \min(2\theta w^+_0, 1)]$, $w^+_0 = \rho/(\rho + 1)$; $B_j \ge 0$ for $j \ge 1$
   (Lemma 1).  $TV(u) = TV_C(u) + TV_F(u_F; u_0)$.
2. *Fine region.*  The four fine stages have the ghosts $g(0)$, $g(\tfrac12)$ and $g(\tfrac12)$, $g(1)$.  Lemma 3(a) on the
   two fine sub-steps (entry cell 1, $c = \mu_1$) gives $TV_F(u^{n+1}_F; p_0) \le TV_F(u_F; u_0) + |p_0 - u_0|$, the path
   of the ghost being monotone, of length $|p_0 - u_0|$; then
   $TV_F(u^{n+1}_F; u^{n+1}_0) \le TV_F(u^{n+1}_F; p_0) + |u^{n+1}_0 - p_0|$.
3. *Coarse region without the correction.*  The cells $j \le -1$ take the coarse Heun step, whose faces read no fine cell.
   With $v_0$ the Heun value of the cell 0 with its own flux at both stages, Lemma 3(b) gives
   $TV_C(\dots, u^{n+1}_{-1}, v_0) \le TV_C(u) - \tfrac12\lambda_0\,(|X| + |Y|)$.
4. *Correction.*  The cell 0 receives on the face $\tfrac12$ the mean of the four fine face values in place of the mean of its
   own two.  With $q_{1/2} = w_0 + \beta_0\,(w_0 - w_{-1})$ on each state, the parts $w_0$ have the same mean,
   $\tfrac14\,(u_0 + 2g(\tfrac12) + p_0) = \tfrac12\,(u_0 + p_0)$, and the interpolated differences are
   $(1 - \tau)\,d_0 + \tau\,d^p_0$, so that

   $$
   \delta_0 = -\lambda_0\,\bigl[(-\tfrac14\beta^{(1)} + \tfrac18\beta^{(2)} + \tfrac18\beta^{(3)})\,d_0
   + (\tfrac18\beta^{(2)} + \tfrac18\beta^{(3)} + \tfrac14\beta^{(4)} - \tfrac12\beta^{c})\,d^p_0\bigr],
   $$

   with $\beta^{(s)}$ that of the cell 0 at the fine stage $s$ and $\beta^c$ at its own stage 2: $|\delta_0| \le
   \lambda_0\,\theta\,(\tfrac14|d_0| + \tfrac12|d^p_0|)$ with the cells $-1$ and 0 of equal width.  Only the last difference of
   $C$ changes: $TV_C(u^{n+1}) \le TV_C(\dots, u^{n+1}_{-1}, v_0) + |\delta_0|$.
5. *Sum.*  With $|p_0 - u_0| = \lambda_0|X|$, $v_0 - p_0 = \tfrac12\lambda_0\,(X - Y)$ and
   $|u^{n+1}_0 - p_0| \le |v_0 - p_0| + |\delta_0|$, steps 2–4 give
   $TV(u^{n+1}) - TV(u^n) \le \lambda_0|X| + \tfrac12\lambda_0|X - Y| - \tfrac12\lambda_0(|X| + |Y|) + 2|\delta_0| = R$.
   Then $\tfrac12(|X| - |Y| + |X - Y|) \le |X| \le (1 + \theta)\,|d_0|$, and $|d^p_0| \le |d_0| + |d_{-1}|$ by Lemma 2 on the
   coarse cells, whose coefficients lie in $[0, 1]$. ∎

*Proof, fine cell upstream.*

1. *Regions.*  The fine region $F = \{\dots, -1, 0\}$ reads the cell 1 only in the slope $s_0 = \beta_0\,(w_0 - w_{-1})$,
   $\beta_0 \in [0, \min(2\theta w^-_0, 1)]$: it is a region with the outflow cell 0, with $B_j \ge 0$.  The coarse region
   $G = \{1, 2, \dots\}$ has the fine cell 0 as its ghost, and $B_1 = 1 + \beta_1 - \alpha_0$ lies in
   $[1 - \min(2\theta\rho/(\rho + 1), 1),\ 1 + \min(2\theta/(\rho + 1), 1)] \subset [0, 2]$; $B_j \ge 0$ for $j \ge 2$.  $TV(u) = TV_F(u_F) + TV_G(u_G; u_0)$.
2. *Fine region.*  Lemma 3(b) on the two fine sub-steps gives $TV_F(u^{n+1}_F) \le TV_F(u_F) - D_1 - D_2$, with
   $|u^{n+1}_0 - u_0| \le D_1 + D_2$.
3. *Coarse region without the correction.*  The coarse Heun step reads the ghost at $u_0$ at its stage 1 (face value $G_s$) and
   at $u^{n+1}_0$ at its stage 2 (face value $\Phi_H$); Lemma 3(a), with the entry cell 1 and $c = \lambda_1$, gives
   $TV_G(u^H; u^{n+1}_0) \le TV_G(u_G; u_0) + |u^{n+1}_0 - u_0|$.  Here (B5) of the cell 1 is used, since $B_1$ reaches
   $1 + \min(2\theta/(\rho + 1), 1)$.
4. *Correction.*  It changes the cell 1 only, by $\tfrac12\lambda_1\,(\Phi - \Phi_H)$, and the cell 1 enters two terms of
   $TV_G$: $TV_G(u^{n+1}_G; u^{n+1}_0) \le TV_G(u^H; u^{n+1}_0) + \lambda_1|\Phi - \Phi_H|$.
5. *Sum.*  Steps 2–4 give $TV(u^{n+1}) \le TV(u^n) + R$.
6. *The form of $R$.*  With $A = \tfrac14\sum_s q_{1/2}(S_s)$ over the states $S_s$ of the four fine stages, $G_s = q_{1/2}(S_1)$
   and $q_{1/2}(S) = w_0 + \beta_0\,(w_0 - w_{-1})$ on each state, $\Phi - \Phi_H = \tfrac12\,[q_{1/2}(S_2) + q_{1/2}(S_3) +
   q_{1/2}(S_4)] - \tfrac12\,q_{1/2}(S_1) - q_{1/2}(S^{n+1})$, with $S^{n+1}$ the fine cells at their final values.  The parts
   $w_0$, through the Heun stages of the fine cell 0, give $\tfrac12\mu_0\,(B^{(4)}_0d^{(4)}_0 - B^{(1)}_0d^{(1)}_0)$, and the
   slopes give $\tfrac12\,(\beta^{(2)}d^{(2)}_0 + \beta^{(3)}d^{(3)}_0 + \beta^{(4)}d^{(4)}_0) - \tfrac12\beta^{(1)}d^{(1)}_0 -
   \beta^{n+1}d^{n+1}_0$; with $B \le 1 + \theta$ and $\beta \le \theta$ on cells of equal width,
   $R \le \lambda_1\,[\tfrac12 m\,(|d^{(1)}_0| + |d^{(4)}_0|) + \theta\,(\tfrac12\sum_s |d^{(s)}_0| + |d^{n+1}_0|)]$.
7. *Locality.*  In $F$ a stage gives $|d'_j| \le |d_j| + |d_{j-1}|$ (Lemma 2, coefficients in $[0, 1]$): on the coefficients of
   $|d_0|, |d_{-1}|, \dots$ a stage adds the vector shifted by one cell upstream.  From $[1]$ this gives $[1, 1]$ after the
   first stage, $[1, 1, \tfrac12]$ at the end of the first sub-step, $[1, 2, \tfrac32, \tfrac12]$ after the third stage and
   $[1, 2, 2, 1, \tfrac14]$ at the end, which in step 6 give the bound of the statement. ∎

### Corollary (a solution Lipschitz continuous next to the face)

Let $\Delta x_*$ be the largest distance between the centres of two neighbouring cells among those whose differences $R$
reads (with cells of equal width within each level, $\Delta x_* = \Delta x_0$).  If $u$ is Lipschitz continuous with the constant
$K$, the cell averages have $|u_j - u_{j-1}| \le K\,(\Delta x_{j-1} + \Delta x_j)/2$.  If up to $T = NH$ every difference that
$R$ reads is at most $K\Delta x_*$, Theorem 3 summed over the $N$ macro steps gives, with the coarse cell upstream (the
first bound of $R$),

$$
TV(u^N) \le TV(u^0) + (1 + \tfrac52\theta)\,a K T\,\frac{\Delta x_*}{\Delta x_0},
$$

and with the fine cell upstream

$$
TV(u^N) \le TV(u^0) + (1 + 3\theta)\,a K T\,\frac{\Delta x_*}{\Delta x_1} .
$$

This is the form $TV(\bar u^n) \le TV(\bar u^0) + CT$ of Theorem 4.2 of Hoang et al. [11], here for the finite-volume scheme of
MOSE, and like theirs it asks the solution to be regular next to the face.  Each estimate is for one macro step and carries no
memory from one macro step to the next: it holds for every position of the face at every macro step, that is for every
sequence of two-level maps, with the regularity next to the face at every macro step.  Of the step, the proof uses only (B5), the own Courant number of every cell with the weights of the widths of (B1): not the coarse step measured on the first fine cell, which
the hypothesis (4.15) of [11] bounds and which in MOSE the interface margin bounds.

*Remark 5 (the step condition).*  With one time step for all the cells of a level, $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$
follows from $c\,(1 + \theta) \le 1$, and at a face from a coarser to a finer cell (B5) asks nothing more.  At a face from a
finer to a coarser cell with $r > 1$ it asks, of the coarse cell downstream, $\lambda_1\,(1 + \min(2\theta r/(r + 1), 1)) \le 1$:
with $\theta = 1$ (van Leer, MC, superbee) this is $\lambda_1\,(1 + \theta) \le 1$; with $\theta = \tfrac12$ (minmod) it is
$\lambda_1\,(1 + r/(r + 1)) \le 1$, which $\lambda_1\,(1 + \theta) \le 1$ does not give.  Where the convective limit rules,
$\lambda_1 \le \text{cfl}/\sigma$, so with minmod $\text{cfl} \le \sigma/(1 + r/(r + 1))$ suffices, and $\text{cfl} \le \sigma/2$
for every $r$: with $\sigma = 2$, every $r$ with cfl 1, the largest cfl that the set-up accepts when $L_{\max} > 0$, so
the set-up gives (B5) for every $r$ and every limiter.  It is a hypothesis of the proof, not a condition shown to be
necessary.

*Remark 6 (what the bound does not give).*  The total variation can grow at a face between levels, at first order too: with
the coarse cell upstream, $\lambda_0 = \mu_1 = 1$, cells of widths 1 and $\tfrac12$, $u_j = 0$ for $j \le -1$ and $u_j = 1$
for $j \ge 0$, the macro step gives 0.5 to the cell 0 and 0.375, 0.75, 0.875, 1 to the fine cells: the values stay in
$[0, 1]$, as Theorem 2 says, and the total variation goes from 1 to 1.25.  Away from the face between levels, on cells of any
widths, the stages of one level do not raise the total variation (Lemma 2).  More than three levels, three levels in another arrangement, several dimensions, systems, nonlinear fluxes, and the
positivity and the maximum principle at the faces between levels are outside this appendix.

### Theorem 3b (three levels, a band of one cell)

*Statement.*  Under (B1)–(B3) and (B5), let the map have three levels, with steps $H$, $H/2$ and $H/4$,
and one cell of the middle level (the band) between the coarse and the fine cells, in either order along the flow; the macro step is
the recursion above, every coarser cell read by a stage is taken at its interpolant at the stage time (the coarse cells also by the
stages of the fine cells, two levels finer), and the correction acts on both faces between levels, each with its own accumulator, the
band being the fine side of one face and, once per sub-step, the coarse side of the other.  Then one macro step gives
$TV(u^{n+1}) \le TV(u^n) + R_A + R_B$, with $R_A$ the remainder $R$ of Theorem 3 at the face between the coarse cells and the band
(the band in place of the fine cell or of the fine region) and $R_B$ the sum, over the two sub-steps of the band, of the remainder $R$
of Theorem 3 at the face between the band and the fine cells, with the band as the coarse cell and its own Courant number
$\mu_b = aH/(2\Delta x_b)$; with the fine cells upstream, $R_A$, the remainder of the downstream face, also carries one half of the remainder of the
upstream face in the second sub-step of the band, the correction that the band receives there.  With the coarse cells upstream every
difference read by $R_A + R_B$ lies at $t_n$ in the cells $-1$, $0$ and the band, and with $\theta = 1$ and a width ratio of 2,
$R_A + R_B \le \lambda_0\,(3.5\,|d_0| + |d_{-1}|) + \mu_b\,(7.34\,|d_b| + 5.34\,\lambda_0\,|d_0|)$, so that
$TV(u^N) \le TV(u^0) + 17.2\,aKT\,\Delta x_*/\Delta x_0$ when those differences stay below $K\Delta x_*$ up to $T = NH$.  With the fine
cells upstream the remainders read the differences of the last fine cell and of the band at the stages of the band; with $\theta = 1$
and a width ratio of 2, if every difference they read stays below $K\Delta x_*$, $TV(u^N) \le TV(u^0) + 13\,aKT\,\Delta x_*/\Delta x_2$,
$\Delta x_2$ the width of the coarse cells.  (B5) is required of the band with its own Courant number and, with the fine cells
upstream, of the first coarse cell after the band with $\min(2\theta \cdot 2/3, 1)$ (Remark 5).  As Theorem 3, this is a bound and not the
TVD property: the total variation can grow at a face between levels (Remark 6).

*Outline of the proof.*  Lemmas 1–3 hold as they stand.  The one new step is Lemma 3 written for a region of one cell: a band $b$ with
Courant number $c$ that takes a Heun sub-step with the ghosts $g_1$ and $g_2$, with $X$ and $Y$ the increments of Lemma 1 at its two
stages and $\tilde v_b$ the value of the sub-step, has $|\tilde v_b - g_2| \le |w_b - g_1| + |g_2 - g_1| - \tfrac12 c\,(|X| + |Y|)$.  The recursion is
nested: within each sub-step of the band, the band and the fine cells form the configuration of Theorem 3, with the band as the
coarse cell; on the scale of the macro step, the band and the fine cells together take the place of the fine region of Theorem 3
towards the coarse cells, with the same path of the ghost.  With the fine cells upstream the band reads the fine cell at its current
value and the coarse cells read the band at its final value; the correction that the band receives in its second sub-step reappears
in the remainder of the downstream face.  Each face between levels then gives its own remainder.

---

## References

1. C.-W. Shu, S. Osher, "Efficient implementation of essentially
   non-oscillatory shock-capturing schemes," *J. Comput. Phys.*,
   77(2), 1988.
2. S. Gottlieb, C.-W. Shu, E. Tadmor, "Strong stability-preserving
   high-order time discretization methods," *SIAM Rev.*, 43(1), 2001.
3. A. Brandt, "Multi-level adaptive solutions to boundary-value
   problems," *Math. Comp.*, 31(138), 1977.
4. A. Jameson, "Solution of the Euler equations for two-dimensional
   transonic flow by a multigrid method," *Appl. Math. Comput.*,
   13(3–4), 1983.
5. A. Arnone, M.-S. Liou, L. A. Povinelli, "Multigrid time-accurate
   integration of Navier–Stokes equations," NASA TM-106373
   (ICOMP-93-37), 1993.
6. P. C. E. Jorgenson, R. V. Chima, "An unconditionally stable
   Runge–Kutta method for unsteady flows," NASA TM-101347, 1989.
7. M. H. Carpenter, D. Gottlieb, S. Abarbanel, W.-S. Don, "The theoretical
   accuracy of Runge–Kutta time discretizations for the initial boundary
   value problem: a study of the boundary error," *SIAM J. Sci. Comput.*,
   16(6), 1241–1252, 1995, doi:10.1137/0916072 (ICASE Report 93-83, 1993).
8. W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction
   Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam,
   2000.
9. E. M. Constantinescu, A. Sandu, "Multirate timestepping methods for
   hyperbolic conservation laws," *J. Sci. Comput.*, 33(3), 239–278, 2007,
   doi:10.1007/s10915-007-9151-y; the page cited is that of the authors'
   report, Virginia Tech, https://hdl.handle.net/10919/19488.
10. G. Jeanmasson, I. Mary, L. Mieussens, "On some explicit local time
   stepping finite volume schemes for CFD," *J. Comput. Phys.*, 397, 2019,
   108818, doi:10.1016/j.jcp.2019.07.017 (Prop. RK2: section on the
   second-order scheme and its conservative correction).
11. T.-T.-P. Hoang, L. Ju, W. Leng, Z. Wang, "High order explicit local
   time stepping methods for hyperbolic conservation laws," *Math.
   Comp.*, 89(324), 2020, 1807–1842, doi:10.1090/mcom/3507
   (arXiv:1905.09705).
12. M. J. Berger, P. Colella, "Local adaptive mesh refinement for shock
   hydrodynamics," *J. Comput. Phys.*, 82(1), 1989.
13. N. Y. Gnedin, V. A. Semenov, A. V. Kravtsov, "Enforcing the
   Courant–Friedrichs–Lewy condition in explicitly conservative local time
   stepping schemes," *J. Comput. Phys.*, 359, 2018, 93–105,
   doi:10.1016/j.jcp.2018.01.008 (arXiv:1801.03108; section 2,
   conditions (3) and (4)).
14. PETSc, `TSSetMaxStepRejections` (option `-ts_max_step_rejections`) in
   `src/ts/interface/ts.c` and the default `max_reject = 10` of `TSCreate` in
   `src/ts/interface/tscreate.c`; a step rejected more often fails with
   `TS_DIVERGED_STEP_REJECTED` (`src/ts/impls/explicit/rk/rk.c`),
   https://github.com/petsc/petsc, commit febf3a0.
15. SUNDIALS, ARKODE: the defaults `MAXNCF` and `MAXCONSTRFAILS` (10) and
   `MAXNEF` (7) of `src/arkode/arkode_impl.h`, the failures of the nonlinear
   solver, of the inequality constraints and of the error test at which a step
   ends (`src/arkode/arkode.c`), set by `ARKodeSetMaxConvFails`,
   `ARKodeSetMaxNumConstrFails` and `ARKodeSetMaxErrTestFails`
   (`src/arkode/arkode_io.c`), https://github.com/LLNL/sundials, commit 312fc0f.
16. Castro: `change_max` (1.1 by default, `Source/driver/_cpp_parameters`) in the check before and after the advance of
   `Source/driver/Castro_advance.cpp`, the fall of the time step tolerated within a step; `retry_subcycle_factor` (0.5) in
   `Source/driver/Castro_advance_ctu.cpp`, with `dt_cutoff` (1e-12) on the landing of the sub-steps,
   https://github.com/AMReX-Astro/Castro, commit 371ba40.
17. E. Hairer, G. Wanner, the code DOPRI5, `WORK(2)`, the safety factor in step size prediction, 0.9 by default; and PETSc
   [14], the default `safety = 0.9` of `TSAdaptCreate` in `src/ts/adapt/interface/tsadapt.c`.
18. SeisSol, `docs/local-timestepping.rst` (the maximum difference property), https://github.com/SeisSol/SeisSol, commit a526776.
