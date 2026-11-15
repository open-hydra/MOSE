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
the unified **Shu–Osher convex-combination** form:

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

### Forward Euler

$$
\mathbf{U}^{n+1} = \mathbf{U}^n + \Delta t\;\mathcal{L}(\mathbf{U}^n)
$$

First-order accurate, one residual evaluation per step.  Useful for
initial testing or very stiff problems where a minimal time step is
already imposed by chemistry or acoustics.

---

### 2-stage SSP Runge–Kutta

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

### 3-stage SSP Runge–Kutta

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

with the stage coefficients $c_k$:

$$
\mathbf{C} =
\begin{pmatrix}
1 & 0 & 0 \\[2pt]
1 & \tfrac{1}{2} & 0 \\[2pt]
1 & \tfrac{1}{4} & \tfrac{2}{3}
\end{pmatrix}
$$

Row = number of stages of the scheme; column = stage.

### Stage times and time-dependent boundary data

Stage $k$ evaluates the residual of $\mathbf{U}^{(k-1)}$, which approximates
the solution at the stage time $t^n + \theta_k\,\Delta t$:

| Scheme | $\theta_1$ | $\theta_2$ | $\theta_3$ |
|--------|:---:|:---:|:---:|
| Forward Euler | 0 | | |
| SSP-RK2 (Heun) | 0 | 1 | |
| SSP-RK3 | 0 | 1 | ½ |

In time-accurate runs the time-dependent boundary data (tabulated inlet
conditions, assigned boundary states) are evaluated at the stage time; with Strang splitting,
inside each half step. This is the conventional imposition analysed by
Carpenter et al. [7]. On a fixed mesh, with $\Delta t \to 0$, both schemes keep
their order in time (second for Heun, third for SSP-RK3). Refining mesh and
step together at a fixed CFL number, the boundary error lowers the global
order of the third-order scheme towards two (in the fourth-order experiments
of [7], table I.a, the rate settles at 2.5 in the L2 norm); Heun stays second
order.

Both statements follow from the boundary error derived below. The
experiments of [7] show the same loss with a weak (SAT) imposition of the
data, and its appendix A shows that data entering as a source term give the
same scheme; the characteristic boundary fluxes of MOSE are outside that
analysis and are checked by the test suite.

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

Steps 2–5 can advance one of three sets of variables:

- Conservative (the default): the steps above.
- Primitive: the residual of each cell is turned into the increment of the primitive variables (partial densities,
  velocity, pressure) by the chain rule, linearised about the state at the start of the stage, and added to them. The
  mass, momentum and energy of a cell then change by the residual only to first order in the increment: a step does not
  conserve them exactly.
- Preconditioned: as the primitive update, with the Weiss–Smith preconditioning (below). It changes the transient, so
  it serves steady runs only.

### Preconditioned update

The preconditioned update solves, at each stage and in every cell,

$$
\Gamma\,\Delta\mathbf Q = -\frac{\Delta t}{V}\,\mathbf R,\qquad
\mathbf Q = (p,\ u,\ v,\ w,\ T,\ Y_1, \dots, Y_{N-1}),
$$

where $\mathbf R$ is the residual of the conservative variables and $\Gamma$ is the matrix of Weiss and Smith [19]: the
Jacobian $\partial\mathbf U/\partial\mathbf Q$ in which $\rho_p = (\partial\rho/\partial p)_T = 1/(RT)$ is replaced by

$$
\Theta = \frac{1}{U_r^2} - \frac{\rho_T}{\rho\,c_p} = \frac{1}{U_r^2} + \frac{1}{c_p T},
$$

with $U_r$ the reference velocity of the cell (below). The rows of $\Gamma\,\Delta\mathbf Q$ are the changes of the
total density, $\Delta\rho^\ast = \Theta\,\Delta p + \rho_T\,\Delta T + \sum_i (\partial\rho/\partial Y_i)\,\Delta Y_i$,
of the momentum, $\mathbf v\,\Delta\rho^\ast + \rho\,\Delta\mathbf v$, of the energy,
$H\,\Delta\rho^\ast + \rho\,(c_p\,\Delta T + \sum_i h_i\,\Delta Y_i + \mathbf v\cdot\Delta\mathbf v) - \Delta p$, and of the
species, $Y_i\,\Delta\rho^\ast + \rho\,\Delta Y_i$. MOSE solves this system in closed form, then converts $\Delta\mathbf Q$ to
the primitive variables with the physical equation of state.

Where $U_r = a$, $\Theta = \rho_p$ (for a thermally perfect
gas $1/a^2 + 1/(c_pT) = 1/(RT)$) and the update is the primitive one. Below $a$ the pressure moves faster than the
density: the acoustic waves of the preconditioned system travel at about $U_r$, close to the flow speed, so at low Mach
number all the waves of a cell have speeds of one order and the local time step is set by the flow, not by the speed of
sound. The turbulence, soot and passive-scalar variables are not preconditioned: they take their residual as in the
primitive update.

Preconditioning changes the path of a steady run, not the steady state of a given flux. Only the HLLC-PC flux takes
$U_r$ ([Riemann solvers](riemann-solvers.md#hllc-pc)), and it is the one that keeps the low-Mach accuracy: its
dissipation scales with $U_r$ instead of $a$, so the pressure fluctuations of a flow at Mach number $M$ are of order
$M^2$, the incompressible limit, instead of order $M$. With another flux the preconditioned update meets a dissipation
scaled with $a$ and a time step scaled with $U_r$: on the [low-Mach cylinder](../vv/2D-low-mach-cylinder.md) at Mach
0.01 the preconditioned update with the HLLC flux diverges. The preconditioned update is paired with HLLC-PC.

### Reference velocity of the preconditioning

The preconditioned update of a stage, the local time step (below) and the HLLC-PC flux take in every cell the reference
velocity $U_r$ of the Weiss–Smith preconditioning, the same value in the three. It is computed after the update of each
stage from the new state of the cell, and the ghost cells of the faces between blocks receive it with the state at the
start of the next stage. It is defined in one of three ways, with a floor fraction $\varepsilon$ of the speed of sound and
a factor $f$ on the flow speed:

| Mode | $U_r$ | Viscous floor |
|---|---|---|
| Uniform, with a reference speed $U_\text{ref}$ | $U_\text{ref}$ where $\lvert\mathbf v\rvert \le U_\text{ref}$, $a$ elsewhere | no |
| Target Mach number $M_t$ | the $U_r$ whose preconditioned Mach number $\lvert\mathbf v\rvert(1-\alpha)/\sqrt{\alpha^2\lvert\mathbf v\rvert^2+U_r^2}$ equals $M_t$, kept within $[\varepsilon a, a]$; $a$ where the flow is faster than the target | yes |
| Weiss–Smith (default) | $\varepsilon a$ for $\lvert\mathbf v\rvert < \varepsilon a$, $a$ for $\lvert\mathbf v\rvert \ge a$, $f\,\lvert\mathbf v\rvert$ otherwise | yes |

$\alpha = \tfrac12(1 - U_r^2/a^2)$. A lower bound $U_{\min}$ can then raise $U_r$ in every mode, and $n$ smoothing passes
give each cell the largest $U_r$ of itself and its six face neighbours, $n$ times.
$\varepsilon$ sets how far the preconditioning reaches: at Mach numbers below $\varepsilon$ the floor $\varepsilon a$
holds $U_r$ above the flow speed, and the dissipation of HLLC-PC grows back towards that of a speed of sound; the
default, 0.05, suits flows down to about Mach 0.05. For a flow at Mach number $M$, take $\varepsilon$ about $M/10$.
The viscous floor is

$$
U_{\text{visc}} = \max\Bigl( \frac{\mu + \mu_t}{\rho\,\Delta_{\min}},\;
\frac{k + \mu_t\,c_p/\text{Pr}_t}{\rho\,c_p\,\Delta_{\min}} \Bigr),
$$

with $\Delta_{\min}$ the smallest width of the cell. With a RANS model $\mu_t$ is the eddy
viscosity of the model with the velocity gradient of the cell: SST takes $\mu_t = \rho a_1 k/\max(a_1\omega, S F_2)$
and Wilcox 2006 $\mu_t = \rho k/\hat\omega$, $\hat\omega = \max\bigl(\omega, C_{\lim}\sqrt{2\bar S_{ij}\bar
S_{ij}/\beta^*}\bigr)$, so the strain rate lowers it; Spalart–Allmaras does not use the gradient. The gradient is
that of the state the stage starts from, computed after the ghost cells of the stage are filled, with the cell metric of
the source terms, as the turbulence sources take it; the setup computes it from the initial state.

After the update the
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

On a skewed 3-D cell the viscous rate carries an extra factor, at most 9/8
([Skewed cells in 3-D](local-time-step.md#skewed-cells-in-3-d)), and on the 2-D axisymmetric wedge a term $s_h$ is added
([The axisymmetric wedge](local-time-step.md#the-axisymmetric-wedge)).

The rule adds up the convective and the diffusive limits of the directions, as an explicit scheme without dimensional
splitting needs. Its constants are the stability limits of the time scheme:

| Scheme | real stability interval | CFL ≤, MUSCL, smooth flow | CFL ≤, first order (limiter at extrema and shocks) | VNN ≤ | CFL and VNN ≤ for positivity with the limiter (equal rectangular cells, smooth flow) |
|---|---|---|---|---|---|
| Heun | $[-2, 0]$ | 1 | 1 | 1/2 | 1/2 |
| SSP-RK3 | $[-2.5127, 0]$ | 1.17 | 1.25 | 0.62 | 1/2 |
| forward Euler | $[-2, 0]$ | unstable | 1 | 1/2 | 1/2 |

The limits are linear (von Neumann, smooth solutions). With the limiter, positivity is guaranteed only for smaller steps:
both CFL and VNN at most 1/2 in smooth flow (last column).

With moving shocks and Heun, keep CFL at 0.5–0.7; with strong shocks and RK3, keep CFL ≤ 0.8. The local step changes the path
of a steady run, not the steady solution it converges to.

With preconditioning, $|\mathbf{v}\cdot\nabla\xi_d| + a|\nabla\xi_d|$ is replaced by the preconditioned eigenvalue
$\bigl(u_n (1-\alpha) + c^\ast\bigr)|\nabla\xi_d|$, with $u_n$ the velocity normal to the $d$ faces,
$\alpha = \tfrac12(1 - U_r^2/a^2)$ and $c^\ast = \sqrt{\alpha^2 u_n^2 + U_r^2}$.

**Recommended values.** For time-accurate runs with moving shocks and RANS models: with RK3, CFL = 0.8 (the face Courant
numbers of the HLLC flux at a captured shock, [derived here](local-time-step.md#shocks-and-the-hllc-flux), stay within the
first-order limit of the scheme, which applies at the shock: 1.47 × 0.8 = 1.18 < 1.256) and VNN = 0.5, 20 % below the limit 0.628 of the scheme (up to that limit, for example at VNN = 0.6,
4.5 % below it, the linearised scheme is stable); with Heun, CFL from 0.5 to 0.7 and VNN = 0.4, below its limit 0.5.

The CFL number of a
run with local time steps, as a steady run, is a parameter of the pseudo-time, not a Courant number of the physical time: it is
not a value to carry to a time-accurate run.

!!! info "Going deeper"
    [Local Time Step: Derivation and Evidence](local-time-step.md) gives where the constants come from, the proof of the
    limits, the evidence next to the axis of a wedge, positivity with the limiter, and the Courant numbers at captured shocks.

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
([Turbulence](turbulence-numerics.md#time-accurate-runs-explicit-sources)), the Reynolds-stress model SSG/LRR included (its bound
leaves out the part of the return to isotropy that the production drives).

!!! tip "CFL ramp-up"

    The time step can be ramped linearly over the first $N$ iterations (no ramp by default): at
    iteration $n$ below $N$ the local step of every cell is multiplied by $n/N$, with $n$ the iteration number of
    the grid level (1 at its first iteration). This improves robustness when the initial
    condition is far from the solution, and a time-accurate start of a $k$–$\omega$ model
    far from the wall value of $\omega$ needs it
    ([Turbulence](turbulence-numerics.md#time-accurate-runs-explicit-sources)).

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
satisfies, with the smoothing coefficient $\varepsilon$,

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
accelerates the convergence of steady (pseudo-time) runs, and MOSE applies it to
steady runs only. In unsteady computations residual smoothing
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

By default only the flow residuals (species, momentum, energy) are
smoothed; the turbulence equations then advance with the same local time
step without smoothing, and above the explicit limit they become the
unstable part. Optionally every residual (turbulence, soot and passive
scalars too) is smoothed in the same way. That suits Spalart–Allmaras: the
k-omega models diverge with it (the SST and Wilcox 2006 flat plates of the
tests at $\varepsilon = 0.5$).
Smoothing changes the path to the steady state, not the steady state
itself (a zero residual stays zero).

---

## Multigrid Acceleration

MOSE supports **geometric multigrid** with a $2\!:\!1$ coarsening ratio
in each spatial direction, to start the fine-grid iterations from a
coarse-grid solution.

### Grid hierarchy

- **Number of levels**: MGL.
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
    singleton directions ($=1$) are exempt. For example, with MGL = 5
    ($2^{4} = 16$) a 1-D mesh needs $n_i$ a multiple of $16$ (e.g. $n_i = 4000$,
    $n_j = n_k = 1$). If this is not met, grid generation / setup will reject the
    mesh. Choose $n_i$ (and $n_j$ in 2-D) accordingly **before** generating the
    grid — the hierarchy is built from this mesh, not the other way around.

### Level schedule

The levels run in sequence, coarsest first: level $n$ runs a given number of
iterations (or until the residual or time criterion of the run, which act on
every level), then its solution is prolongated onto level $n-1$. The output of
the coarse levels, prolongated to the finest grid, is described in the
[User Guide](../user/using.md#output).

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
conservative state via Newton–Raphson, using the
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
levels; each coarse level runs its own number of iterations, and its
solution is prolongated to the next finer level, down to level 1. The
residual criterion of a steady run and the final time also end a level. The
iteration count starts again from zero on each level, and each level advances
its own time from the initial time of the run (zero in a new run).

There are
no V or W cycles, no coarse-grid correction and no dual time stepping: the
coarse levels only provide the starting solution of the finer ones.

---

## Multirate Time Stepping

In a time-accurate run with the global step every cell advances with the
smallest admissible step of the mesh, $\Delta t = \min_i \Delta t_i$.  When a
few cells (near walls, in a throat, in a refined band) set that minimum, most
of the work is spent advancing the other cells with steps far below their own
limit.  The **multirate** integrator lets each cell advance with a step close
to its own limit, keeping, on a fixed mesh, the second order in time, and the
exact conservation of the single-rate scheme.

The idea and the guarantees are in [Multirate Time Stepping](multirate.md); the construction in detail is in
[Multirate Method](multirate-method.md) and the proofs of its order and stability in [Multirate Analysis](multirate-analysis.md).
The input keys are in the [User Guide](../user/input.md#multirate-time-stepping) and its verification in
[Multirate Verification](../vv/multirate.md).

---

## References

[1] C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.*, 77(2), 1988.

[2] S. Gottlieb, C.-W. Shu, E. Tadmor, "Strong stability-preserving high-order time discretization methods," *SIAM Rev.*, 43(1), 2001.

[3] A. Brandt, "Multi-level adaptive solutions to boundary-value problems," *Math. Comp.*, 31(138), 1977.

[4] A. Jameson, "Solution of the Euler equations for two-dimensional transonic flow by a multigrid method," *Appl. Math. Comput.*, 13(3–4), 1983.

[5] A. Arnone, M.-S. Liou, L. A. Povinelli, "Multigrid time-accurate integration of Navier–Stokes equations," NASA TM-106373 (ICOMP-93-37), 1993.

[6] P. C. E. Jorgenson, R. V. Chima, "An unconditionally stable Runge–Kutta method for unsteady flows," NASA TM-101347, 1989.

[7] M. H. Carpenter, D. Gottlieb, S. Abarbanel, W.-S. Don, "The theoretical accuracy of Runge–Kutta time discretizations for the initial boundary value problem: a study of the boundary error," *SIAM J. Sci. Comput.*, 16(6), 1241–1252, 1995, doi:10.1137/0916072 (ICASE Report 93-83, 1993).

[8] W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam, 2000.

[19] J. M. Weiss, W. A. Smith, "Preconditioning applied to variable and constant density flows," *AIAA J.*, 33(11), 2050–2057, 1995, doi:10.2514/3.12946.
