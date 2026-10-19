# Oblique Shock

A Mach 4.0 supersonic flow encountering a wedge ramp at 33.0° deflection angle. This configuration generates a single oblique shock whose properties are governed by the oblique shock relations from inviscid shock theory. The test verifies that numerical solutions across different Riemann solvers accurately satisfy these theoretical shock jump conditions.

**Reference**: standard oblique shock relations (_e.g._ Anderson, _Modern Compressible Flow_, McGraw-Hill).

---

## Problem setup

Uniform Mach 4 flow enters from the left at unit pressure and temperature. The flow encounters a sharp ramp deflection angle $\delta$. The oblique shock is uniquely determined by the shock angle $\beta$ through the theta-beta-Mach relation (Ames Research Staff, _Equations, Tables, and Charts for Compressible Flow_, NACA Report 1135, p. 621, eqs. (139a)–(139b)):

$$\tan(\delta) = 2\cot(\beta) \frac{M_1^2 \sin^2(\beta) - 1}{M_1^2 \left(\gamma + \cos(2\beta)\right) + 2}$$

Post-shock properties (pressure ratio, temperature ratio, density ratio, exit Mach) are uniform. 

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+Tramel, HLLE, SLAU** |

## Grid structure

The mesh is composed of two structured blocks:

- **Block 1** (pre-shock): $20 \times 218$ cells over $[0, 0.1] \times [0, 0.15]$
- **Block 2** (post-shock): $38 \times 218$ cells over the parallelogram between $x = 0.1$ and $x = 0.3$, $0.15$ high, whose lower side is the ramp

!!! note
    The grid is built to trigger the unstable behaviour of some Riemann solvers.

Boundary conditions:

- Inlet: Riemann invariant extrapolation with Mach 4 state
- Wall: slip wall (symmetry)
- Upper boundary: zero-gradient extrapolation
- Exit: zero-gradient extrapolation

## Results and verification

Comparison of numerical shock jump conditions across all Riemann solvers:

<figure>
  {% include "vv/images/OS-fields.svg" %}
</figure>

Theoretically, a uniform flow should exist behind the oblique shock wave. However, most Riemann solvers except HLLE exhibit instabilities that produce an undesirable zebra-like oscillatory pattern in the post-shock region. The HLLC+ solver, which is explicitly designed to suppress such numerical artifacts, delivers the cleanest solution.

The transversal density profile normalized with the analytical solution is plotted along a vertical line in the post-shock region:

<figure>
  {% include "vv/images/OS-slice.svg" %}
</figure>

The HLLE solver produces a quite uniform flow with density matching the analytical solution very closely. In contrast, HLLC exhibits errors up to 10% compared to the theoretical value, while HLLC+ remains slightly unstable with oscillations of approximately 2% around the analytical density.
