# Rocket Nozzle – Multi-Chemistry

A 2D axisymmetric Euler simulation of a rocket nozzle operating on RP-1/LOX propellants mocking the Falcon 9 Merlin1D. Three configurations are tested: **frozen composition** (fixed species frozen at chamber state), **equilibrium** (species adjust instantaneously to maintain thermochemical equilibrium), and **finite-rate** (species composition evolves according to reaction kinetics). Mass-flow-weighted quasi-1D averages of the 2D solutions are compared with 1D CEA (Chemical Equilibrium with Applications) references at the area ratios of the CEA stations. `verify.py` prints and draws this comparison; it has no pass criterion. The finite-rate configuration has a check of its own, `finite-rate/verify.py`: it stops with an error when pressure, density, the two velocity components or the ratio of specific heats of the run differ by more than a relative $10^{-3}$ from those of the stored solution `reference/field.tec`, in the cell on the axis that contains the throat of the contour, x = 0.

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

`verify.py` reads the solution of each run and averages Mach number, pressure, temperature, ratio of specific heats and molar
mass over each column of cells, with the mass flow as weight. It places each CEA station where the wall of the grid has the
area ratio of the station, the area being taken relative to the narrowest section, at the wall node of smallest radius, and
interpolates the averages linearly in x between the centres of the two columns around the station. The station at 16 is the
nozzle exit: the grid reaches an area ratio of 15.993 there, and the solution at the exit is that of the last column, whose
state the exit face takes. The frozen run stops when its residual has fallen to $10^{-10}$ of its value at the first iteration (6196
iterations); the equilibrium run does not converge (after its 10 000 iterations the residual is 3 times its first value); the
finite-rate run stops at its 2000 iterations, with the residual at 1.7 % of its first value.

The Mach number is the velocity over the frozen speed of sound, $\sqrt{\gamma p/\rho}$ with $\gamma$ the ratio of specific heats
of the mixture at fixed composition, which MOSE writes for every configuration, for the three runs and for both CEA references.
CEA gives the Mach number of its equilibrium stations with the equilibrium speed of sound; `verify.py`
divides their velocity by the frozen speed of sound of the same state, from the specific heat with frozen reactions of
`reference/CEA.out` (at the equilibrium throat this Mach number is 0.974, where CEA prints 1).

In the frozen and equilibrium runs temperature, ratio of specific heats and molar mass stay within 1.8 % of the CEA reference
of the same model at the six stations; Mach number and pressure differ more, as the list below gives. The averages at the
stations change by up to 0.21 % with the order of the interpolation (linear, quadratic or cubic), and the Mach numbers of CEA are
printed with three decimals, 0.16 % at the subsonic station.

- **Frozen** against CEA frozen: temperature from −0.31 % to +1.79 % at the six stations, ratio of specific heats within 0.06 %,
  the molar mass of the chamber mixture; Mach number 13.3 % above CEA at the subsonic station and from −2.6 % to −0.6 % at
  the others; pressure 6.2 % to 12.8 % above CEA in the diverging part.
- **Equilibrium** against CEA equilibrium: temperature from −0.22 % to +1.05 %, ratio of specific heats within 0.03 %, molar
  mass within 0.15 %; Mach number 13.3 % above CEA at the subsonic station and from −2.5 % to −0.5 % at the others; pressure
  5.8 % to 11.9 % above CEA in the diverging part.
- **Finite-rate**: close to equilibrium: molar mass within 0.45 % of CEA equilibrium (and up to 4.1 % above CEA frozen),
  temperature from −1.57 % to +0.80 % of CEA equilibrium.

<figure>
  {% include "vv/images/Merlin-XY.svg" %}
  Mass-flow-weighted 1D averages of Mach number (with the frozen speed of sound), temperature ratio T/T0, ratio of specific
  heats and molar mass of the three configurations, against the CEA stations (frozen and equilibrium).
</figure>
