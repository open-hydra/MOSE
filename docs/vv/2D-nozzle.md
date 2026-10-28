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
iterations). The residual is the norm, over the cells, of the change of the density in one iteration, and the local time step
grows linearly over the first 100 iterations (`cfl-rise-threshold = 100`), so the first iteration takes a step a hundred times
smaller. The residual of the equilibrium run peaks at iteration 130, falls by a factor of 59 and then oscillates about 3.3 times
its first value (between 2.1 and 7.3 times from iteration 1000 on) up to its 10 000 iterations: the solution keeps changing, most
in the three rows of cells next to the axis, and in a run continued to 20 000 iterations the averages of the comparison below
change by less than a relative $2 \cdot 10^{-4}$ between iterations 10 000 and 20 000. The finite-rate run stops at its 2000
iterations, with the residual at 1.7 % of its first value: with the Strang splitting of the chemistry the residual measures the
change of the second half step of the transport alone, which does not vanish when the solution is steady; in a run continued to
10 000 iterations the residual stays at 1.7 % of the first one from iteration 4000 on, and the averages at the stations change by
less than a relative $10^{-7}$ between iterations 2000 and 10 000.

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

CEA computes a quasi-one-dimensional flow: one uniform state on each cross-section, moving along the axis, with the entropy and
the total enthalpy of the chamber (Gordon and McBride, NASA RP-1311, 1994, chapter 6). The frozen and equilibrium solutions keep
the mass flow, the total enthalpy and the entropy along the nozzle, up to the discretisation error, but their cross-sections are
not uniform: the flow is inclined to the axis, by the wall angle at the wall (−41° at the subsonic station, from 25° to 10° in the
divergent part), and the pressure changes from the axis to the wall (at the exit of the frozen run, 26 kPa on the axis and
100 kPa at the wall). The uniform state with the mass flow, the total enthalpy and the entropy of a cross-section is the CEA state
at the area ratio divided by the discharge coefficient $C_D = 0.996$, the ratio of the mass flow of the run to that of the
quasi-one-dimensional flow through the same throat; NASA SP-8120 (*Liquid Rocket Engine Nozzles*, 1976, figure 3) gives for
potential flow about 0.9945 when the radius of the wall arc upstream of the throat is 1.5 times the throat radius, as on this
contour (about 0.991 when it is the throat radius, about 0.998 when it is three times). In both runs this state differs from that
CEA state by at most 0.8 % in pressure and 0.33 % in Mach number at the five stations away from the throat. The mass-flow-weighted
averages differ more because they average a non-uniform cross-section. At the subsonic station the mass flux along the velocity
exceeds the axial one by 10 % (mass-flow-weighted cosine of the flow angle 0.91), which raises the averaged Mach number by 12 %,
and the change across the section adds 1.6 %. In the divergent part the flow angle and the change of pressure across the section
raise the averaged pressure by 6.5 % to 13.5 % and lower the averaged Mach number by 0.6 % to 2.9 %. The discharge coefficient
lowers the Mach number by 0.5 % at the subsonic station and the pressure by 0.5 % to 0.7 % in the divergent part, where it raises
the Mach number by 0.1 % to 0.2 %; the discretisation error, the digits that CEA prints and, at the exit, the half cell between
the centre of the last cells and the exit plane make up the rest of the differences listed above (at most 0.8 % each). At the
geometric throat the sonic line is curved and the cross-section holds subsonic and supersonic flow. Temperature, ratio of specific
heats and molar mass depend on pressure along the isentrope and differ less from CEA (at fixed composition,
$dT/T = ((\gamma - 1)/\gamma)\, dp/p$). `verify.py` does not compute these parts of the differences: they come from an analysis of
the solutions of the frozen and equilibrium runs.

<figure>
  {% include "vv/images/Merlin-XY.svg" %}
  Mass-flow-weighted 1D averages of Mach number (with the frozen speed of sound), temperature ratio T/T0, ratio of specific
  heats and molar mass of the three configurations, against the CEA stations (frozen and equilibrium).
</figure>
