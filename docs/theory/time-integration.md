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
selected at run time via `n_RK` (1, 2, or 3):

| Scheme | Order | SSP coefficient |
|--------|:-----:|:---------:|
| Forward Euler | 1 | 1 |
| 2-stage SSP Runge–Kutta (Heun) | 2 | 1 |
| 3-stage SSP Runge–Kutta (Shu–Osher) | 3 | 1 |

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

The **Shu–Osher RK3** method — the default integrator:

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
Carpenter et al. [5]. On a fixed mesh, with $\Delta t \to 0$, both schemes keep
their order in time (second for Heun, third for SSP-RK3). Refining mesh and
step together at a fixed CFL number, the boundary error lowers the global
order of the third-order scheme towards two (in the fourth-order experiments
of [5], table I.a, the rate settles at 2.5 in the L2 norm), also when the
boundary data enter through the boundary flux as in MOSE; Heun stays second
order. Keeping the full order would need intermediate boundary values built
from the time derivatives of the data [5], which tabulated data do not
provide.

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

---

## Stability Conditions

### Local time step

Each cell $i$ gets a local time step that combines convection and diffusion in one additive rule
(W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes,
Thomas Stieltjes Institute, 2000, §6 and §9.3; A. Arnone, M.-S. Liou, L. A. Povinelli, *Multigrid
Time-Accurate Integration of Navier–Stokes Equations*, NASA TM-106373, 1993, eq. 26–28):

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
the largest of the two terms (Hundsdorfer, §9.3, which states the form without a proof). For forward
Euler the form follows, with the diffusive rates added as well, from the positivity condition (6.1)
of Hundsdorfer: with first-order upwind convection and central diffusion on a Cartesian mesh the
diagonal entry of a row is $-\sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr)$ and the others are not negative, so
$\tau \sum_d \bigl(|a_d|/h_d + 2D/h_d^2\bigr) \le 1$: the rule above with CFL = 1 and VNN = 1/2. The
constants of Heun and RK3 in the table do not come from that source: they are a von Neumann analysis
of the linearised MUSCL–upwind advection with central diffusion, made for MOSE, and at the worst
cell Péclet number the additive law is exact. The source warns that such an eigenvalue criterion can
mislead for the non-normal operators that boundaries give, as an inflow boundary with upwind
differences (Hundsdorfer, §3.3, pp. 26–27): the limits are those of the interior.

| Scheme | CFL ≤ | VNN ≤ |
|---|---|---|
| Heun (`RK2`) | 1 | 1/2 |
| `RK3` | 1.17 | 0.62 |
| explicit Euler (`euler`) | 1 with first-order reconstruction; unstable with MUSCL | 1/2 |

These limits are linear (von Neumann, smooth solutions). With the limiter, positivity and the TVD property
are guaranteed only for smaller steps: the van Leer slope of MOSE is the flux limiter $\psi(\theta) =
\theta/(1+\theta)$, which satisfies $0 \le \psi \le 1$ and $0 \le \psi/\theta \le \mu$ with $\mu = 1$, and with
Heun and RK3 (SSP coefficient 1) positivity holds for $(1+\mu)\,\nu \le 1$, that is
$\sum_d \nu_d \le 1/2$; central diffusion is positive for $\text{VNN} \le 1/2$; with the additive rule
both hold for CFL ≤ 1/2 and VNN ≤ 1/2 (Hundsdorfer, §5.2, §6.2 and §6.4; experiments reported there give
positive results up to $\nu = 1$ in 1D and about $2/3$ in 2D). With moving shocks and Heun keep `cfl` at
0.5–0.7. With preconditioning,
$|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$ is replaced by the preconditioned eigenvalue
$\bigl(u_n (1-\alpha) + c^\ast\bigr)|\nabla\xi_d|$.

### Global time step

In a time-accurate run all cells advance with the smallest local step:

$$
\Delta t = \min_i \Delta t_i .
$$

In time-accurate RANS runs the turbulence sources are explicit, and each local step also includes
their stiffness, the derivative $d_i$ of the destruction terms evaluated at the start of the step:
$1/\Delta t_i$ gains the term $d_i$, so that $\Delta t\,d_i \le 1$ and no stage makes $\tilde\nu$, $k$ or
$\omega$ negative through their destruction
([Turbulence](turbulence.md#time-accurate-runs-explicit-sources)).

!!! tip "CFL ramp-up"

    During the initial transient MOSE supports a linear CFL ramp from a
    low starting value to the target CFL over a user-specified number of
    iterations (`rampa_iter`).  This improves robustness when the initial
    condition is far from the steady state.

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
satisfies

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
$1 + 2\varepsilon$.  Typical values of the smoothing parameter $\varepsilon$
range from 0.1 to 0.2.

### Implementation

1. Extrapolate residuals to ghost cells (boundary conditions)
2. For each of the three coordinate directions:
    - Perform 2 Jacobi iterations of the 3-point stencil
    - Pass the smoothed residual to the next direction
3. Use the smoothed residual in the RK stage

---

## Multigrid Acceleration

MOSE supports **geometric multigrid** with a $2\!:\!1$ coarsening ratio
in each spatial direction to accelerate convergence to steady state.

### Grid hierarchy

- **Number of levels**: user-specified (`MGL`).
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

The coarse-grid correction is transferred to the fine grid with
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

### Cycle structure

- **V-cycle** (default): one pre-smoothing sweep, restriction,
  coarse-grid solve, prolongation, one post-smoothing sweep.
- **W-cycle**: supported through recursive calling of the multigrid
  driver.
- Multigrid can also be used for **transient** problems to accelerate
  the convergence of each physical time step.

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
5. M. H. Carpenter, D. Gottlieb, S. Abarbanel, W.-S. Don, "The theoretical
   accuracy of Runge–Kutta time discretizations for the initial boundary
   value problem: a study of the boundary error," *SIAM J. Sci. Comput.*,
   16(6), 1241–1252, 1995, doi:10.1137/0916072 (ICASE Report 93-83, 1993).
