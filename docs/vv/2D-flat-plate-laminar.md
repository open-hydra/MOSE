# 2D Flat Plate – Laminar Boundary Layer

Incompressible laminar flow over a flat plate at low Mach number. This configuration establishes a developing Blasius boundary layer whose properties are well-characterized by classical similarity theory. Verification is performed by comparison with the exact Blasius solution for skin friction coefficient and velocity profiles.

**Reference**: Blasius similarity solution for laminar boundary layer over a flat plate (see _e.g._ Anderson, _Introduction to Flight_, McGraw-Hill).

<figure>
  {% include "vv/images/LP-field.svg" %}
</figure>

---

## Problem setup

Uniform incompressible flow at Mach 0.2 enters from the left over a semi-infinite flat plate. The leading edge induces a boundary layer that develops downstream following the laminar Blasius solution. Post-stagnation, a no-slip adiabatic wall condition is imposed.

The Blasius equation governs the self-similar boundary layer profile:

$$f''' + f \cdot f'' = 0$$

with boundary conditions $f(0) = f'(0) = 0$ and $f'(\eta \to \infty) = 1$, where $\eta = y\sqrt{U_\infty / 2\nu x}$ is the similarity variable.

The skin friction coefficient follows:

$$C_f = \frac{\sqrt{2} \cdot f''(0)}{\sqrt{Re_x}}$$

where $Re_x = U_\infty x / \nu$ is the local Reynolds number.

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 0.2 |
| Temperature | 300 K |
| Pressure | 0.97250 Pa |
| Viscosity | 1.0 × 10⁻³ Pa·s |
| Re/m | 80.65 m⁻¹ |

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK2 |
| CFL | 0.5 |
| VNN | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |

## Grid structure

The mesh is a 2D structured grid (65 × 65 cells) spanning the physical domain with fine resolution near the leading edge and near-wall region to capture boundary layer development accurately.

Boundary conditions:

- **Inlet**: Riemann invariant boundary condition with freestream state
- **Adiabatic wall** (x ≥ 0): No-slip, no-heat-flux
- **Symmetry** (x < 0): Symmetry boundary condition upstream
- **Outlets**: Atmospheric ambient conditions

## Results and verification

Comparison of skin friction coefficient along the plate surface versus the analytical Blasius distribution:

<figure>
  {% include "vv/images/LP-Cf.svg" %}
</figure>

The solution demonstrates excellent agreement with Blasius theory away from the leading edge. Leading-edge errors (x < 0.05 m) are expected due to numerical stagnation effects and coarse initial mesh resolution; Blasius theory assumes an infinitesimal leading edge, whereas practical CFD has finite spatial resolution.

**Velocity profile**

Comparison of the normalized velocity profile $u/U_\infty$ versus the Blasius similarity coordinate $\eta$ at x = 0.29 m:

<figure>
  {% include "vv/images/LP-slice.svg" %}
</figure>

The numerical velocity profile shows excellent collapse onto the Blasius self-similar solution, confirming that the boundary layer develops according to classical theory.

