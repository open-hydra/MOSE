# Theoretical Guide

This section provides the theoretical foundations for the physical
models and numerical methods implemented in MOSE.  The material spans
the governing equations of compressible multi-species flow, the
finite-volume spatial discretization, Riemann solvers, time-integration
schemes, turbulence closures, and the thermochemistry framework.

## How this guide is organised

The pages of this guide explain what each part of MOSE does and why: the idea, the main equations and the values the solver
uses. Where a topic has a longer derivation, a proof or a body of evidence behind it, that material is in a separate **deep
dive**, linked from the page that needs it. Read the overview first; go to a deep dive when you want to check or extend the
argument.

| Deep dive | Backs | What it holds |
|-----------|-------|---------------|
| [Mesh Geometry and Block Interfaces](numerics-geometry.md) | [Spatial Discretization](numerics.md) | Cell volume, the checks of the cells, the flux of a face between blocks |
| [Face Gradients and Metrics](numerics-gradients.md) | [Spatial Discretization](numerics.md) | The stencil and the face metric of the viscous fluxes, ghost-cell metrics, chimera faces |
| [Ghost Cells at Walls and Extrapolation Faces](numerics-ghost-cells.md) | [Spatial Discretization](numerics.md) | The extrapolation rules of each boundary type |
| [Axisymmetric Wedge and 3-D Sector](numerics-wedge.md) | [Spatial Discretization](numerics.md) | Rotational periodicity, the wedge terms, the sector |
| [Local Time Step](local-time-step.md) | [Time Integration](time-integration.md#local-time-step) | Where the limits of the time step come from, their proof, and the evidence near the axis, with the limiter and at shocks |
| [Multirate Method](multirate-method.md) | [Multirate Time Stepping](multirate.md) | The construction of levels and macro step, the non-local CFL condition, the recursion, boundaries, turbulence and MPI |
| [Multirate Analysis](multirate-analysis.md) | [Multirate Time Stepping](multirate.md) | Proofs of the order, stability and monotonicity |
| [Turbulence Model Variants](turbulence-variants.md) | [Turbulence Modelling](turbulence.md) | SA variants; SST compressible production, dilatation dissipation, rotation-curvature correction, energy coupling |
| [Turbulence Numerical Treatment](turbulence-numerics.md) | [Turbulence Modelling](turbulence.md) | Wall ghost cells, cell gradients, wall distance, advance of the source terms |
| [Surface Models in Detail](gsi-models.md) | [Gas–Surface Interaction](gsi.md) | Energy and mass balances of pyrolysis, solid propellant, surface reactions and the pyrolysing ablator |

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Compressible Navier–Stokes for multi-species ideal-gas mixtures: conservative form, state vectors, convective and viscous fluxes, source terms, boundary conditions

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretization__

    ---

    Cell-centred finite volume method, MUSCL reconstruction, five flux limiters, Jameson shock detection

    [:octicons-arrow-right-24: Spatial discretization](numerics.md)

-   :material-rotate-3d:{ .lg .middle } __Single Rotating Frame__

    ---

    Coriolis and centrifugal forces, relative velocity formulation, turbomachinery applications, Spalart–Shur correction

    [:octicons-arrow-right-24: Rotating frame](rotating-frame.md)

-   :material-waves:{ .lg .middle } __Riemann Solvers__

    ---

    14 approximate and exact solvers: AUSM family, HLL family, Godunov, Lax–Friedrichs, SLAU

    [:octicons-arrow-right-24: Riemann solvers](riemann-solvers.md)

-   :material-arrow-expand-horizontal:{ .lg .middle } __Boundary Conditions__

    ---

    The state imposed by outlets, inlets at total conditions and the choked nozzle, and the directions of one cell

    [:octicons-arrow-right-24: Boundary conditions](boundary-conditions.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    SSP Runge–Kutta, CFL/VNN stability, Strang splitting, implicit residual smoothing, multigrid, multirate local time stepping

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-weather-windy:{ .lg .middle } __Turbulence Modelling__

    ---

    Spalart–Allmaras (SA, SA-R, SA-RC), Menter SST k–ω (SST, SST-RC), Wilcox 2006 k–ω

    [:octicons-arrow-right-24: Turbulence models](turbulence.md)

-   :material-thermometer:{ .lg .middle } __Thermodynamic and Transport Properties__

    ---

    Equation of state, mixture rules, property evaluation, and Wilke transport models

    [:octicons-arrow-right-24: Thermodynamic models](thermo.md)

-   :material-chart-bell-curve:{ .lg .middle } __Finite-Rate Kinetics__

    ---

    Arrhenius reactions, three-body collisions, Lindemann falloff, Troe formulation

    [:octicons-arrow-right-24: Kinetics models](kinetics.md)

-   :material-scale-balance:{ .lg .middle } __Chemical Equilibrium__

    ---

    NASA CEA algorithm for constant-UV equilibrium problems

    [:octicons-arrow-right-24: Equilibrium solver](equilibrium.md)

-   :material-fire:{ .lg .middle } __Gas-Surface Interaction__

    ---

    Ablating walls: melting, pyrolysis (HTPB, HDPE, PP), heterogeneous surface reactions, surface energy balance and wall temperature solution

    [:octicons-arrow-right-24: Gas-surface interaction](gsi.md)

</div>

---

## Overview

The MOSE solver advances the compressible Navier–Stokes equations in
conservative form on structured multi-block grids.  The numerical
pipeline can be summarised as follows:

| Stage | Method | Page |
|-------|--------|------|
| **Governing system** | Multi-species Euler / Navier–Stokes | [Governing Equations](governing-equations.md) |
| **Spatial discretization** | Cell-centred FVM + MUSCL-limited reconstruction | [Spatial Discretization](numerics.md) |
| **Interface fluxes** | 14 Riemann solvers (AUSM, HLL, Roe, SLAU, Godunov, Lax–Friedrichs) | [Riemann Solvers](riemann-solvers.md) |
| **Boundaries** | Outlets, inlets at total conditions, choked nozzle | [Boundary Conditions](boundary-conditions.md) |
| **Time marching** | 3-stage SSP RK3, IRS, multigrid, multirate local time stepping | [Time Integration](time-integration.md), [Multirate](multirate.md) |
| **Turbulence closure** | SA, SST $k$–$\omega$, Wilcox 2006 $k$–$\omega$ | [Turbulence Modelling](turbulence.md) |
| **Thermochemistry** | Ideal-gas mixture, NASA polynomials | [Thermodynamics](thermo.md) |
| **Chemistry** | Finite-rate (Arrhenius/Troe) or equilibrium (CEA) | [Kinetics](kinetics.md), [Equilibrium](equilibrium.md) |
| **Ablating walls** | Melting, pyrolysis and surface-reaction boundary conditions | [Gas-Surface Interaction](gsi.md) |

---
