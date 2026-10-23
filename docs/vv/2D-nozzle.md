# Rocket Nozzle – Multi-Chemistry

A 2D axisymmetric Euler simulation of a rocket nozzle operating on RP-1/LOX propellants mocking the Falcon 9 Merlin1D. Three configurations are tested: **frozen composition** (fixed species frozen at chamber state), **equilibrium** (species adjust instantaneously to maintain thermochemical equilibrium), and **finite-rate** (species composition evolves according to reaction kinetics). Area-weighted quasi-1D averages extracted from 2D CFD solutions are compared against 1D CEA (Chemical Equilibrium with Applications) references at discrete area-ratio stations. The tests verify multi-species solver accuracy across different thermochemical assumptions.

**Reference**: CEA (Chemical Equilibrium with Applications), NASA Glenn Research Center.

<figure>
  {% include "vv/images/Merlin-field.svg" %}
  Contours of Mach number from the 2D frozen-flow simulation of the Merlin nozzle.
</figure>

---

## Problem setup

Combustion products of RP-1 burned with liquid oxygen at an oxidiser-to-fuel ratio of 2.35 are expanded from a stagnation pressure of 108 bar and a stagnation temperature of 3644 K through the nozzle contour.

### Grid structure

A single structured block with 120 axial × 40 radial cells covers the full nozzle from the converging inlet to the nozzle exit. The lower boundary is the axis of symmetry and the upper boundary follows the Merlin ideal contour.

Boundary conditions (same for all configurations):

- Inlet (`401`, $T_0 = 3644$ K, $p_0 = 108$ bar): CEA chamber equilibrium state (subsonic inflow) of the RP-1/LOX combustion products, with 13 species (frozen, equilibrium) or 19 (finite-rate: the 17 of the CORIA mechanism, and COOH and HCOOH, which it leaves inert)
- Exit: zero-gradient supersonic extrapolation
- Axis: symmetry
- Wall: symmetry (inviscid)

## Test configurations

| Parameter | Value |
|---|---|
| Time scheme | RK3 (frozen), RK2 (equilibrium, finite-rate) |
| CFL | 0.9 |
| IRS beta | 0.3 (frozen), 0.5 (equilibrium, finite-rate) |
| Space reconstruction | MUSCL |
| Flux limiter | vanleer |
| Riemann solver | SLAU |

- **Frozen composition** Species mass fractions are fixed at the chamber-equilibrium state and do not react during the expansion.  
- **Equilibrium chemistry** Species composition adjusts instantaneously at each point to maintain thermochemical equilibrium.  
- **Finite-rate chemistry** Species composition evolves according to reaction kinetics using the CORIA reaction mechanism (17 species, 44 reactions) with NASA9 thermodynamic database.


## Results and verification

The 2D averages from each configuration are compared against 1D CEA references. The frozen-composition solver reproduces isentropic expansion closely, validating the multi-species Euler solver. The equilibrium configuration is verified to ensure that thermochemistry updates are correctly integrated with the flow solver across the nozzle expansion. Finite-rate chemistry results are compared to confirm that the solver can handle non-equilibrium effects, although the high-temperature, high-pressure conditions in the nozzle lead to near-equilibrium behavior. Overall, the tests demonstrate accurate multi-species flow predictions across different thermochemical assumptions in a rocket nozzle geometry.

<figure>
  {% include "vv/images/Merlin-XY.svg" %}
  Area-weighted 1D averages of Mach number, pressure ratio, and temperature ratio compared against CEA reference stations.
</figure>

!!! note "Figure and code version"
    The curves are drawn from the solutions stored in `reference/field.tec` of each configuration, which
    come from the first upload of MOSE (May 2026, commit `ab2bf537`; the equilibrium one after 2000
    iterations). MOSE at commit `10b59003b72c` moves the equilibrium curve of the fourth panel by 2.4 % of
    the range of its axes, and every other curve by less than 1 %. 1.7 % of it is already there with the
    code and the inputs of the first upload, run for 2000 or 10 000 iterations: they do not give the
    stored equilibrium solution again. The other 0.7 % comes from the changes between the code of 13 October
    2026 (commit `b29130b`) and `10b59003b72c`.
