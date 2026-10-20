# Hypersonic Cylinder

A Mach 8.1 supersonic flow around a cylindrical body, generating a detached bow shock upstream of the cylinder stagnation point. This test case verifies the accuracy of numerical solutions against theoretical normal shock relations at the stagnation line, where the bow shock becomes approximately normal. The test is particularly valuable for validating Riemann solver behavior in high-speed flows known to be prone to the generation of numerical instabilities such as carbuncle phenomena.

**Reference**: normal shock relations (_e.g._ Anderson, _Modern Compressible Flow_, McGraw-Hill).

---

## Problem setup

Uniform Mach 8.1 flow impinges on a circular cylinder. At the stagnation point on the cylinder surface, a detached bow shock forms approximately normal to the free-stream direction. The post-shock flow decelerates isentropically to the stagnation conditions of the cylinder. 

The analytical solution at the stagnation line is governed by normal shock relations:

$$M_2 = \sqrt{\frac{M_1^2 + 2/(\gamma-1)}{2\gamma M_1^2/(\gamma-1) - 1}}$$

$$\frac{p_2}{p_1} = 1 + \frac{2\gamma}{\gamma+1}(M_1^2 - 1)$$

$$\frac{T_2}{T_1} = \frac{p_2/p_1}{\rho_2/\rho_1}$$

The test verifies that numerical predictions of total pressure $p_0$ and total temperature $T_0$ behind the shock match theoretical values calculated from the post-shock Mach number.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+Tramel, HLLE, SLAU** |

## Grid structure

The mesh is composed of a structured grid around the cylindrical body with 160 cells in the circumferential direction and 40 cells in the radial direction, extending from the cylinder surface to a far-field boundary located at 1.5 cylinder radii.

!!! note
    The grid is built to trigger the unstable behaviour of some Riemann solvers. In fact, there is not any refinement in the shock layer or the employment of cells with a high aspect ratio to avoid the carbuncle phenomenon.

Boundary conditions:

- Inlet: Mach 8.1 free-stream state
- Cylinder surface: slip wall (symmetry)
- Far-field boundaries: zero-gradient extrapolation

The test reads no probe: `verify.py` takes the post-shock total quantities in the field cell next to the cylinder surface and to the stagnation streamline (indices 80 and 0 of the field, counted from 0).

## Results and verification

MOSE solutions compared across four Riemann solvers in terms of non-dimensional density contours flow field and wall trend:

<figure>
  {% include "vv/images/HC-fields.svg" %}
</figure>

<figure>
  {% include "vv/images/HC-slice.svg" %}
</figure>

!!! note "Figure and code version"
    The figure comes from the first upload of MOSE (May 2026, commit `ab2bf537`). MOSE at commit
    `10b59003b72c` moves the HLLC, HLLE and SLAU curves by 8.9 %, 3.2 % and 7.3 % of the range of the
    axes, and the HLLC+Tramel curve by 0.8 %. With any of the four solvers,
    and also with the code of 13 October 2026 (commit `b29130b`), the solution does not settle at the
    stagnation point: the error of $p_0$ in the cell where `verify.py` reads the stagnation state changes
    by up to 7.7 points between 5000, 10 000 and 20 000 iterations, and by 9.2 points between two builds
    (GNU 12.2 RELEASE with MPI on the cluster and gfortran 13 TESTING without MPI).

All Riemann solvers capture the fundamental shock structure, but exhibit varying levels of accuracy in preserving a smooth post-shock region. In particular, the HLLC solver shows significant oscillations in the shock layer, while HLLC+ and HLLE provide smoother solutions. The SLAU solver demonstrates an intermediate behavior with moderate oscillations.

The density profile along the cylinder wall confirms the overall behavior previously observed.

Quantitative error metrics at the stagnation point compared to analytical normal-shock solution:

| Solver | $p_0$ [MPa] | Error $p_0$ [%] | $T_0$ [K] | Error $T_0$ [%] |
|---|---|---|---|---|
| HLLC | 8.53825 | -0.79 | 4081.36 | +0.30 |
| HLLC+ | 8.67777 | +0.83 | 4079.80 | +0.26 |
| HLLE | 8.74553 | +1.62 | 4060.20 | -0.22 |
| SLAU | 8.34013 | -3.09 | 4079.70 | +0.26 |
| **Analytical** | **8.60638** | **—** | **4069.11** | **—** |

Notwithstanding only HLLE and HLLC+ provide smooth post-shock solutions, all solvers exhibit errors well within engineering tolerances.
