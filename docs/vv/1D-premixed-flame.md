# Premixed Laminar Flame

A one-dimensional, freely-propagating **premixed CH$_4$/air laminar flame** at
stoichiometry ($\phi = 1.0$). Unburnt reactants enter at the laminar flame speed
so that the flame is stationary in the computational frame; the case exercises the
**finite-rate chemistry**, **multicomponent transport**, and **compressible
reacting-flow** coupling together, and verifies that MOSE reproduces the flame
structure of an independent low-Mach reference.

Because no analytical solution exists, MOSE is validated against a **Cantera 1-D
freely-propagating flame** solved with the same kinetic mechanism and transport
data. The reference fixes the flame speed, the adiabatic flame temperature, and
the full species structure, so any discrepancy is numerical.

---

## Problem setup

A stoichiometric CH$_4$/air mixture ($Y_{\mathrm{CH_4}} = 0.0552$,
$Y_{\mathrm{O_2}} = 0.2201$, $Y_{\mathrm{N_2}} = 0.7247$) enters at $T_u = 300$ K
and $p = 101325$ Pa. The inflow velocity is set to the laminar flame speed,
$u_n = S_L = 0.2860$ m/s, anchoring the flame in the domain. The initial field is
a sharp unburnt/burnt discontinuity that relaxes to the steady flame structure.

| State | $T$ [K] | $u$ [m/s] | composition |
|---|---:|---:|---|
| Unburnt (reactants) | 300  | 0.286 | CH$_4$/O$_2$/N$_2$ at $\phi = 1$ |
| Burnt (products)    | 2153 | 2.428 | CO$_2$, H$_2$O, CO, OH, N$_2$, … |

Mass conservation across the flame fixes the burnt-gas velocity through the
density (temperature) ratio, $u_b/u_u = \rho_u/\rho_b = T_b/T_u \approx 7.2$.

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier–Stokes (laminar) + finite-rate chemistry |
| Mechanism | ZK (25 species) |
| Domain | $[0,\,0.30]$ m |
| Grid | $4000 \times 1$ |
| Multigrid | 5 levels (FMG start-up) |
| Time scheme | RK2 (time-accurate), Strang splitting for chemistry |
| ODE solver | ROS4 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| Riemann solver | HLLC (Batten) |
| Schmidt / Prandtl | 0.7 / 0.7 |

The inflow is a subsonic inlet and the outflow a fixed-pressure outlet; the spanwise faces use symmetry, reducing the problem to 1-D along $i$.

## Results and verification

### Density and temperature

Density (left axis) and temperature (right axis) against the progress variable
$c = (T - T_u)/(T_b - T_u)$, which aligns the two flames without any spatial
shift. Solid lines are MOSE; open symbols are Cantera.

<figure>
  {% include "vv/images/flame-thermo.svg" %}
</figure>

### Species structure

Major species ($Y_{\mathrm{CH_4}}, Y_{\mathrm{O_2}}, Y_{\mathrm{CO_2}},
Y_{\mathrm{H_2O}}$, left axis) and minor species ($Y_{\mathrm{OH}},
Y_{\mathrm{CO}}$, right axis) against $c$.

<figure>
  {% include "vv/images/flame-species.svg" %}
</figure>

The reactant inflow speed, the burnt-gas state, and the radical (OH) peak all
match the Cantera reference closely:

| Quantity | MOSE | Cantera |
|---|---:|---:|
| $S_L$ (unburnt $u$) [m/s] | 0.287 | 0.286 |
| $T_\text{burnt}$ [K]      | 2193  | 2149  |
| $u_\text{burnt}$ [m/s]    | 2.04  | 2.06  |
| $Y_{\mathrm{OH}}$ peak    | $5.16\times10^{-3}$ | $5.15\times10^{-3}$ |

!!! focus "Partially non-reflective boundaries"
    This case uses **partially non-reflective** inlet and outlet
    conditions, controlled by the reflection factor `rf`: the inflow runs at `rf = 0.5` and the outflow at `rf = 0.001` (nearly fully draining). The `rf` relaxes the characteristic entering the domain toward the imposed target, so outgoing acoustic waves leave instead of reflecting. This keeps the flame anchored and the burnt-gas field clean; a residual low-amplitude ripple remains and averages out, with the mean velocity following the Cantera profile.
