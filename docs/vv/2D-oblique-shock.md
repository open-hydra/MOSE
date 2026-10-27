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

- Inlet (`405`): supersonic inflow at the Mach 4 state
- Wall: slip wall (symmetry)
- Upper boundary: zero-gradient extrapolation
- Exit: zero-gradient extrapolation

## Results and verification

`verify.py` averages $p_2/p_1$, $T_2/T_1$, $\rho_2/\rho_1$ and $M_2$ over a window of block 2 well behind the shock ($i$ = 14 to 34, $j$ = 5 to 49, counted from 0) and passes when every ratio is within 1.5 % of the oblique-shock relations. With HLLE, HLLC+Tramel and SLAU the four ratios are within 0.80 %. HLLC fails, with $T_2/T_1$ +1.94 % and $M_2$ −3.21 % ($p_2/p_1$ +0.42 %): its average over the striations described below is not the post-shock state. The density of both blocks, scaled by the pre-shock density:

<figure>
  {% include "vv/images/OS-fields.svg" %}
</figure>

Behind the oblique shock the exact flow is uniform. HLLC, SLAU and, much less, HLLC+Tramel show density striations along the grid lines behind the shock (a zebra-like pattern); HLLE does not.

The density along the vertical grid line $i$ = 30 of block 2 (cell centres at $x$ = 0.234), scaled by the exact post-shock density (the pre-shock density times the exact ratio):

<figure>
  {% include "vv/images/OS-slice.svg" %}
</figure>

In $0.09 \le y \le 0.12$ on that line the density departs from the exact post-shock value by up to 0.69 % with HLLE, 2.7 % with HLLC+Tramel, 9.2 % with SLAU and 12.0 % with HLLC.
