# 2D Flat Plate – Turbulent Boundary Layer

Incompressible turbulent flow over a flat plate at low Mach number. This configuration validates MOSE's turbulent boundary layer modeling against well-established compressible Navier–Stokes solvers from the NASA suite (CFL3D, FUN3D). Multiple turbulence closure models are tested: Spalart–Allmaras (SA), Shear Stress Transport (SST), and Wilcox k-ω.

**Reference**: [NASA Turbulence Modeling Resource](https://www.nasa.gov/nasa-turbulence-modeling-resource/)

---

## Problem setup

Uniform incompressible flow at Mach 0.2 enters from the left over a flat plate. The boundary layer rapidly transitions from laminar to turbulent flow, developing a fully turbulent boundary layer downstream. A symmetry condition is imposed upstream of the leading edge (pre-stagnation region), and a no-slip adiabatic wall condition is applied over the plate.

The skin friction coefficient is the primary validation metric:

$$C_f = \frac{\tau_w}{0.5 \rho_\infty U_\infty^2}$$

where $\tau_w = \mu \frac{\partial u}{\partial y}\big|_{\text{wall}}$ is the wall shear stress.

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 0.2 |
| Temperature | 555 K |
| Pressure | 1.02828 Pa |
| Viscosity | 1.1859 × 10⁻⁵ Pa·s |
| Re/m | ~110,000 m⁻¹ |

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 1.0 |
| VNN | 1.0 |
| IRS $\beta$ | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |

## Grid structure

The mesh is a 2D structured grid spanning the physical domain with fine resolution near the leading edge and near-wall region to capture boundary layer development accurately.

Boundary conditions:

- **Inlet**: Riemann invariant boundary condition with freestream state
- **Adiabatic wall** (x ≥ 0): No-slip, no-heat-flux
- **Symmetry** (x < 0): Symmetry boundary condition upstream
- **Outlets**: Atmospheric ambient conditions

## Turbulence models

1. **Spalart–Allmaras (SA)**: One-equation model; computationally efficient for wall-bounded flows
2. **Shear Stress Transport (SST)**: Two-equation model; blends k-ω near-wall and k-ε in the freestream
3. **Wilcox k-ω**: Traditional two-equation model; excellent wall-layer resolution

## Results and verification

For each turbulence model, MOSE computes the wall shear stress distribution along the flat plate. The skin friction coefficient is extracted and compared pointwise against reference solutions.

<figure>
  {% include "vv/images/TP-cf-sa.svg" %}
  Spalart–Allmaras
</figure>

<figure>
  {% include "vv/images/TP-cf-sst.svg" %}
  Shear Stress Transport (SST)
</figure>

<figure>
  {% include "vv/images/TP-cf-wilcox.svg" %}
   Wilcox k-ω
</figure>



