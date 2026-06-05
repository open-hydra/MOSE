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

| Scheme | Order | CFL limit |
|--------|:-----:|:---------:|
| Forward Euler | 1 | 1 |
| 2-stage SSP Runge–Kutta (Heun) | 2 | 1 |
| 3-stage SSP Runge–Kutta (Shu–Osher) | 3 | 1 |

All three are **strong-stability preserving (SSP)**: each stage is a
convex combination of forward-Euler steps, so positivity and TVD
properties of the spatial discretization are inherited under the
standard CFL limit.

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

### CFL condition

The convective time-step limit in direction $d$ for cell $i$ is

$$
\Delta t_{\text{CFL},\,i}^{(d)} =
\frac{\Delta x_d}{\bigl|\mathbf{v}\!\cdot\!\hat{\mathbf{e}}_d\bigr| + a}
\;\times\;\text{CFL}
$$

### VNN condition

When a turbulence model is active, the viscous (von Neumann) stability
limit is

$$
\Delta t_{\text{VNN},\,i}^{(d)} =
\frac{\rho\,(\Delta x_d)^2}{\mu_\ell + \mu_t}\;\times\;\text{VNN}
$$

### Global time step

The actual time step is the global minimum:

$$
\Delta t = \min_{i,\,d}\!\bigl[\min\!\bigl(
  \Delta t_{\text{CFL},\,i}^{(d)},\;
  \Delta t_{\text{VNN},\,i}^{(d)}\bigr)\bigr]
$$

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
  For 2-D grids ($k_\max = 1$), coarsening is applied only in the $i$
  and $j$ directions.
- **Compatibility**: the fine-grid dimensions must be divisible by
  $2^{\text{MGL}-1}$ in each active direction.

### Restriction (fine → coarse)

Conservative volume averaging:

$$
\mathbf{U}_\text{coarse} =
\frac{1}{V_\text{coarse}}
\sum_{\text{fine}\,\in\,\text{coarse}}
\mathbf{U}_\text{fine}\;V_\text{fine}
$$

- 3-D: 8-cell average
- 2-D: 4-cell average

Primitive variables on the coarse grid are recovered from the averaged
conservative state via Newton–Raphson (`cons2prim`), using the
fine-grid temperature as the initial guess.

### Prolongation (coarse → fine)

The coarse-grid correction is transferred to the fine grid with
**cubic interpolation** (3-D) or **biquadratic interpolation** (2-D).

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
