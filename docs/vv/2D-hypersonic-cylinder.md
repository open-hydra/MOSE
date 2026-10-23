# Hypersonic Cylinder

A Mach 8.1 supersonic flow around a cylindrical body, generating a detached bow shock upstream of the cylinder stagnation point. The case shows how four Riemann solvers behave in a high-speed flow prone to numerical instabilities such as the carbuncle phenomenon, and compares the total pressure and total temperature next to the stagnation point with the normal-shock relations, since the bow shock is approximately normal at the stagnation line. It has no pass criterion.

**Reference**: normal shock relations (_e.g._ Anderson, _Modern Compressible Flow_, McGraw-Hill).

---

## Problem setup

Uniform Mach 8.1 flow impinges on a circular cylinder. At the stagnation point on the cylinder surface, a detached bow shock forms approximately normal to the free-stream direction. The post-shock flow decelerates isentropically to the stagnation conditions of the cylinder. 

The analytical solution at the stagnation line is governed by normal shock relations:

$$M_2 = \sqrt{\frac{M_1^2 + 2/(\gamma-1)}{2\gamma M_1^2/(\gamma-1) - 1}}$$

$$\frac{p_2}{p_1} = 1 + \frac{2\gamma}{\gamma+1}(M_1^2 - 1)$$

$$\frac{T_2}{T_1} = \frac{p_2/p_1}{\rho_2/\rho_1}$$

The reference values of the total pressure $p_0$ and total temperature $T_0$ behind the shock follow from the post-shock state and Mach number through the isentropic relations.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| Time stepping | local (steady), 5000 iterations |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+Tramel, HLLE, SLAU** |

## Grid structure

The mesh is composed of a structured grid around the cylindrical body with 160 cells in the circumferential direction and 40 cells in the radial direction, extending from the cylinder surface to a far-field boundary located at 1.5 cylinder radii.

!!! note
    The grid is built to trigger the unstable behaviour of some Riemann solvers. In fact, there is not any refinement in the shock layer or the employment of cells with a high aspect ratio to avoid the carbuncle phenomenon.

Boundary conditions:

- Outer boundary (`405`): supersonic inlet, $M = 8.1$, $T = 288.15$ K, $p = 101\,325$ Pa
- Cylinder surface (`300`): slip wall
- Outflow ends (`400`): extrapolation
- Faces 5 and 6 (`200`): the closure of the planar 2-D mesh

The test reads no probe: `verify.py` takes the post-shock total quantities in the field cell next to the cylinder surface and to the stagnation streamline (indices 80 and 0 of the field, counted from 0).

## Results and verification

MOSE solutions compared across four Riemann solvers: contours of the density of the flow field divided by its value in the cell with indices 0 and 0 (no colour scale), and the density along the cylinder wall in kg/m³:

<figure>
  {% include "vv/images/HC-fields.svg" %}
</figure>

<figure>
  {% include "vv/images/HC-slice.svg" %}
</figure>

All Riemann solvers capture the fundamental shock structure, but exhibit varying levels of accuracy in preserving a smooth post-shock region. In particular, the HLLC solver shows significant oscillations in the shock layer, while HLLC+Tramel and HLLE provide smoother solutions. The SLAU solver demonstrates an intermediate behavior with moderate oscillations.

The density profile along the cylinder wall confirms the overall behavior previously observed.

Total pressure and total temperature in the cell that `verify.py` reads, after the 5000 iterations of the case, against
the normal-shock relations (GNU 12.2 RELEASE build, the build of the figures):

| Solver | $p_0$ [MPa] | Error $p_0$ [%] | $T_0$ [K] | Error $T_0$ [%] |
|---|---|---|---|---|
| HLLC | 9.26915 | +7.70 | 4119.72 | +1.24 |
| HLLC+Tramel | 8.76383 | +1.83 | 4099.98 | +0.76 |
| HLLE | 8.47606 | −1.51 | 4053.06 | −0.39 |
| SLAU | 8.31264 | −3.41 | 4068.34 | −0.02 |
| **Analytical** | **8.60638** | **—** | **4069.11** | **—** |

`verify.py` prints these values and draws the figures. They are not converged values: the
solution next to the stagnation point does not settle with any of the four solvers, and after the 5000 iterations the
density residual is 3.5 % (HLLE), 6.0 % (SLAU), 8.0 % (HLLC) and 16.6 % (HLLC+Tramel) of its first value. The state
of the cell changes with the iteration and, with HLLC and SLAU, with the build: with gfortran 13 TESTING the error of
$p_0$ is +6.35 % (HLLC) and −2.17 % (SLAU), with ifx 2024.1 RELEASE +6.13 % and −0.58 %; with HLLC+Tramel and HLLE the
three builds agree within 0.009 % ($p_0$) and 0.002 % ($T_0$).
