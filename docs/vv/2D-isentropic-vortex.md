# 2D Isentropic Vortex

A smooth (C-infinity) vortex superimposed on a uniform mean flow is an **exact
solution of the compressible Euler equations**: it is convected without change
of shape.  On a periodic box the vortex returns to its starting point after one
pass, so the exact solution at that time equals the initial condition.  Because
the field is infinitely smooth, the scheme's **observed order of accuracy** can
be measured directly by grid refinement — this is the primary order-verification
case for the 2D solver.

!!! note
    This is *not* the [Gresho vortex](2D-gresho-vortex.md): the Gresho profile has
    derivative kinks (its azimuthal velocity is only piecewise linear), which cap the
    observed order below two, so it probes low-Mach robustness rather than accuracy.
    The isentropic vortex is smooth precisely so that a second-order scheme shows a
    clean second-order slope.

**Reference**: C.-W. Shu, _Essentially Non-Oscillatory and Weighted ENO Schemes_,
ICASE Report 97-65 (1998).

## Problem setup

The vortex of strength $\beta$ is centered at $(x_0, y_0) = (5, 5)$ on a periodic
domain $[0,10]^2$, superimposed on a mean flow $(u_\infty, v_\infty) = (1, 1)$:

$$
\begin{aligned}
r^2 &= (x-x_0)^2 + (y-y_0)^2, \\
\delta u &= -\tfrac{\beta}{2\pi}\,(y-y_0)\,e^{(1-r^2)/2}, \qquad
\delta v = \tfrac{\beta}{2\pi}\,(x-x_0)\,e^{(1-r^2)/2}, \\
\delta T &= -\tfrac{(\gamma-1)\beta^2}{8\gamma\pi^2}\,e^{1-r^2}, \qquad
T = 1 + \delta T, \\
\rho &= T^{1/(\gamma-1)}, \qquad p = \rho^{\gamma} \quad (\text{entropy} = 1).
\end{aligned}
$$

| Parameter | Value |
|---|---|
| Domain | $[0, 10]^2$, periodic |
| Vortex strength $\beta$ | 5 |
| Mean flow $(u_\infty, v_\infty)$ | $(1, 1)$ |
| Final time | $t = 10$ (one full period) |
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |

The mean-flow Mach number is $M_\infty = \sqrt{u_\infty^2+v_\infty^2}/a \approx 1.2$;
the flow stays smooth (no shocks), so the exact solution remains the translated
initial vortex.

## Results and verification

The case is run on three successively refined grids ($32^2$, $64^2$, $128^2$,
refinement ratio 2) and the normalised $L_2$ density error is measured against
the exact solution.  The observed order of accuracy is
$p = \log_2\!\big(e(N)/e(2N)\big)$.

| Grid | Normalised $L_2$ density error | Observed order |
|---|---|---|
| $32^2$  | 2.56 % | — |
| $64^2$  | 0.56 % | 2.18 |
| $128^2$ | 0.17 % | 1.75 |

The error decreases at close to second order (mean $p \approx 2$), confirming the
MUSCL/RK3 scheme is second-order accurate on smooth flow. The slight drop on the
finest pair is the Van Leer limiter mildly clipping the smooth vortex core.

<figure>
  {% include "vv/images/IsentropicVortex.svg" %}
</figure>

## Limiter choice and observed order

The case fixes the **Van Leer** limiter on purpose. Only the *smooth* limiters
(Van Leer, Van Albada) recover a clean second-order slope here; the more
aggressive **MC** and **Superbee** do not — which is expected limiter behaviour,
**not** a loss of accuracy of the base scheme.

Two effects are at play:

- **Godunov barrier.** The vortex core is a *smooth
  extremum*, where the upwind/downwind slopes have opposite signs
  ($r = \Delta_\text{up}/\Delta_\text{down} < 0$). Every TVD limiter returns
  $\phi = 0$ there, so all of them locally drop to first order at the core. This
  is fundamental to TVD reconstruction and unavoidable.
- **Behaviour in the smooth, monotone regions.**
    - *Van Leer* $\phi=\dfrac{r+|r|}{1+|r|}$ and *Van Albada*
      $\phi=\dfrac{r^2+r}{1+r^2}$ are **smooth ($C^1$ at $r=1$)** and
      non-compressive: away from the tiny core the reconstruction is cleanly
      second order, so the global $L_2$ error converges at $\approx 2$.
    - *MC* $\phi=\max\!\big(0,\min(2r,\tfrac{1+r}{2},2)\big)$ is
      **piecewise-linear (kinked)** and steeper, so it clips the smooth extremum
      harder — the observed order dips (a marginal, coarse-grid effect that
      recovers toward 2 as the core is better resolved).
    - *Superbee* $\phi=\max\!\big(0,\min(2r,1),\min(r,2)\big)$ sits on the
      **maximally compressive edge** of the TVD region. It is anti-diffusive by
      design and artificially steepens smooth gradients into staircases
      ("terracing"), so its error barely converges — Superbee is a known poor
      choice for smooth-accuracy tests (it targets keeping contact
      discontinuities sharp, the opposite regime).
