# Woodward-Colella Supersonic Step

A Mach 3.0 supersonic flow over a forward step geometry. This is a classic test case for compressible flow solvers that exercises the scheme's ability to capture oblique shocks, expansion fans, and their interactions in a multi-dimensional setting. The solution features a reflected shock, expansion fan, and contact surface, providing a stringent validation of shock-capturing capabilities.

**Reference**: P. Woodward & P. Colella, _J. Comput. Phys._ **54** (1984) 115–173.

---

## Problem setup

A uniform Mach 3 flow enters from the left boundary and faces a forward step obstacle at the lower boundary. This provides a classical shock-deflection pattern: 
- Primary oblique shock deflecting around the step corner
- Second shock reflecting from the upper boundary
- Expansion fan originating from the step corner
- Contact surface separating regions of different entropy

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[0, 3] \times [0, 1]$ |
| Final time | $t = 4.0$ |
| Mach number | 3.0 |
| Time scheme | RK3 |
| CFL | 0.8 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLE |

## Grid structure

The mesh is composed of two structured blocks:

- **Block 1** (step region): $48 \times 16$ cells over $[0, 0.6] \times [0, 0.2]$
- **Block 2** (downstream): $240 \times 64$ cells over $[0, 3.0] \times [0.2, 1.0]$

Boundary conditions:

- Inlet: Riemann invariant extrapolation with Mach 3 state
- Step wall: slip wall (symmetry)
- Upper boundary: symmetry
- Exit: zero-gradient extrapolation

## Results and verification

MOSE solutions compared against OpenFOAM (rhoCentralFoam) reference on the baseline mesh:

<figure>
  {% include "vv/images/WC-fields.svg" %}
</figure>

Density profile cut at $y = 0.5$ demonstrating sharp resolution of shock structures:

<figure>
  {% include "vv/images/WC-slice.svg" %}
</figure>

The solution accurately captures all major flow features including the oblique shock angle, shock interaction with the upper boundary, and the expansion fan structure. Pressure and density jumps are well-resolved with minimal oscillation.
