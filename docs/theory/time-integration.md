# Time Integration

MOSE advances the semi-discrete equations in time using an explicit
strong-stability-preserving (SSP) Runge–Kutta scheme.  This page
describes the time-stepping algorithm, the CFL and VNN stability
conditions, Strang operator splitting for chemistry coupling, implicit
residual smoothing (IRS) for convergence acceleration, and the
multigrid strategy.

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

| Scheme | real stability interval | CFL ≤, MUSCL, smooth flow | CFL ≤, first order (limiter at extrema and shocks) | VNN ≤ | CFL and VNN ≤ for positivity with the limiter |
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

so that $\Delta t\,(4D_{\max}\sum_d|\nabla\xi_d|^2 + s_h) = 4\,\text{VNN}$. The limits of the table are proven away
from the axis on uniform meshes (in 2-D for any skew and aspect ratio, in 3-D with the factor $F$ above) and, next to the
axis, for the scalars; for the momentum next to the axis they rest on evidence (test `WedgeAxisStep`).

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
scalar advection of each direction with the limited flux: the van Leer slope of MOSE is the flux limiter $\psi(\theta) =
\theta/(1+\theta)$, which satisfies $0 \le \psi \le 1$ and $0 \le \psi/\theta \le \mu$ with $\mu = 1$, and with
Heun and RK3 (SSP coefficient 1) positivity holds for $(1+\mu)\,\nu \le 1$, that is
$\sum_d \nu_d \le 1/2$; central diffusion is positive for $\text{VNN} \le 1/2$; with the additive rule
both hold for CFL ≤ 1/2 and VNN ≤ 1/2 ([8], §5.2, §6.2 and §6.4; experiments reported there give
positive results up to $\nu = 1$ in 1D and about $2/3$ in 2D). With moving shocks and Heun keep `cfl` at
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

IRS is switched on only at setup (`irs = true` or `irs-beta > 0`), where its
arrays are allocated; on a run started with IRS `irs-beta` and `irs-variables`
can change at runtime and `irs = false` switches it off; `irs = true`, or
`irs-beta` above 0, added to a run started without IRS is ignored with a warning.

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
