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

**MC (Monotonicity-Centred)**

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
1 - \tanh\!\bigl(10\,(\sigma\,\Delta)^3\bigr)
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
Cartesian coordinates $(x, y, z)$ using the face metric tensor $M_{3\times 3}$, which is computed from the grid geometry and stored at each face.  The same metric is used to compute the physical spacing $\Delta l$ for the MUSCL reconstruction.

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
