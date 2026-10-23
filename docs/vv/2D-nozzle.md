# Rocket Nozzle – Multi-Chemistry

A 2D axisymmetric Euler simulation of a rocket nozzle operating on RP-1/LOX propellants mocking the Falcon 9 Merlin1D. Three configurations are tested: **frozen composition** (fixed species frozen at chamber state), **equilibrium** (species adjust instantaneously to maintain thermochemical equilibrium), and **finite-rate** (species composition evolves according to reaction kinetics). Mass-flow-weighted quasi-1D averages of the 2D solutions are compared with 1D CEA (Chemical Equilibrium with Applications) references at the area ratios of the CEA stations. `verify.py` prints and draws this comparison; it has no pass criterion.

**Reference**: CEA (Chemical Equilibrium with Applications), NASA Glenn Research Center.

<figure>
  {% include "vv/images/Merlin-field.svg" %}
  Contours of Mach number from the 2D finite-rate simulation of the Merlin nozzle.
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

`verify.py` reads the solution of each run, averages Mach number, pressure, temperature, ratio of specific heats and molar mass
over each cross-section with the mass flow as weight, and takes the cell nearest to the area ratio of each CEA station (it prints
the area ratio of that cell: 1.883 for the subsonic station at 2, 1.944 for the supersonic one at 2, 3.932, 7.919 and 15.898 for
4, 8 and 16). The frozen run stops when its residual has fallen by $10^{-10}$ from the first iteration (6211 iterations); the
equilibrium run does not converge (after its 10 000 iterations the residual is 3 times its first value); the finite-rate run
stops at its 2000 iterations, with the residual at 1.7 % of its first value.

In the frozen and equilibrium runs temperature, ratio of specific heats and molar mass stay within 2.6 % of the CEA reference
of the same model at the six stations; Mach number and pressure differ more, as the list below gives. Part of the Mach
difference comes from the cells
taken, which are not at the area ratios of the stations: at the subsonic station the area ratio 1.883 instead of 2 alone
raises the quasi-1D isentropic Mach number by 7 % (from 0.312 to 0.334 with the ratio of specific heats of the frozen run,
1.214).

- **Frozen** against CEA frozen: temperature within −0.46 % and +2.57 % at the six stations, ratio of specific heats within 0.08 %,
  the molar mass of the chamber mixture; Mach number 19.9 % above CEA at the subsonic station and between −3.8 % and +1.9 % at
  the others; pressure 7.5 % to 17.5 % above CEA in the diverging part.
- **Equilibrium** against CEA equilibrium: temperature within −0.33 % and +1.56 %, ratio of specific heats within 0.03 %, molar
  mass within 0.19 %; Mach number 16.7 % above CEA at the subsonic station and between −5.5 % and −0.8 % at the others; pressure
  7.0 % to 16.4 % above CEA in the diverging part.
- **Finite-rate**: close to equilibrium: molar mass within 0.45 % of CEA equilibrium (and up to 4.1 % above CEA frozen),
  temperature within −1.57 % and +1.36 % of CEA equilibrium.

<figure>
  {% include "vv/images/Merlin-XY.svg" %}
  Mass-flow-weighted 1D averages of Mach number, temperature ratio T/T0, ratio of specific heats and molar mass of the three
  configurations, against the CEA stations (frozen and equilibrium).
</figure>
