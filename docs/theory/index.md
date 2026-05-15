# Theoretical Guide

This section provides the theoretical foundations for the physical
models and numerical methods implemented in MOSE.  The material spans
the governing equations of compressible multi-species flow, the
finite-volume spatial discretization, Riemann solvers, time-integration
schemes, turbulence closures, and the thermochemistry framework.

<div class="grid cards" markdown>

-   :material-math-integral:{ .lg .middle } __Governing Equations__

    ---

    Compressible Navier–Stokes for multi-species ideal-gas mixtures: conservative form, state vectors, convective and viscous fluxes, source terms, boundary conditions

    [:octicons-arrow-right-24: Governing equations](governing-equations.md)

-   :material-grid:{ .lg .middle } __Spatial Discretization__

    ---

    Cell-centred finite volume method, MUSCL reconstruction, seven flux limiters, Jameson shock detection

    [:octicons-arrow-right-24: Spatial discretization](numerics.md)

-   :material-rotate-3d:{ .lg .middle } __Single Rotating Frame__

    ---

    Coriolis and centrifugal forces, relative velocity formulation, turbomachinery applications, Spalart–Shur correction

    [:octicons-arrow-right-24: Rotating frame](rotating-frame.md)

-   :material-waves:{ .lg .middle } __Riemann Solvers__

    ---

    14 approximate and exact solvers: AUSM family, HLL family, Godunov, Lax–Friedrichs, SLAU

    [:octicons-arrow-right-24: Riemann solvers](riemann-solvers.md)

-   :material-timer-outline:{ .lg .middle } __Time Integration__

    ---

    SSP Runge–Kutta, CFL/VNN stability, Strang splitting, implicit residual smoothing, multigrid

    [:octicons-arrow-right-24: Time integration](time-integration.md)

-   :material-weather-windy:{ .lg .middle } __Turbulence Modelling__

    ---

    Spalart–Allmaras (+ SA-RC, SAcomp, SAR), Menter SST k–ω, Wilcox 2006 k–ω

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
| **Interface fluxes** | 14 Riemann solvers (AUSM, HLL, Godunov, LF, SLAU) | [Riemann Solvers](riemann-solvers.md) |
| **Time marching** | 3-stage SSP RK3, IRS, multigrid | [Time Integration](time-integration.md) |
| **Turbulence closure** | SA, SST $k$–$\omega$, Wilcox 2006 $k$–$\omega$ | [Turbulence Modelling](turbulence.md) |
| **Thermochemistry** | Ideal-gas mixture, NASA polynomials | [Thermodynamics](thermo.md) |
| **Chemistry** | Finite-rate (Arrhenius/Troe) or equilibrium (CEA) | [Kinetics](kinetics.md), [Equilibrium](equilibrium.md) |

---
