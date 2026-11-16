# HyShot II Scramjet Combustor

Reacting flow in the HyShot II supersonic combustor at flight condition **Case 2** ($h \approx 28$ km), computed in three dimensions with finite-rate hydrogen/air chemistry on an overset (chimera) grid. This is the most demanding case in the suite: it couples, in a single run, everything the earlier cases exercise one at a time - a supersonic turbulent boundary layer, a shock train, sonic transverse fuel injection into crossflow, multicomponent diffusion, and stiff finite-rate kinetics.

The validation metric is the **wall pressure along the between-injector line**, measured in the HyShot II ground-test campaign and reported with its uncertainty in Karl's thesis, against which a Spalart-Allmaras RANS computed with DLR's TAU code is also available. The comparison is therefore a validation (against measurement) *and* a code-to-code check (against another SA RANS of the same configuration).

**References**: Karl, *Numerical Investigation of a Generic Scramjet Configuration*, PhD thesis, TU Dresden (2011) - the experiment with error bars and the TAU S-A solution; Chapuis et al., *A computational study of the HyShot II combustor performance*, Proc. Combust. Inst. **34** (2013) 2101-2109 - the same measurements in the paper's coordinates, used here only to audit the unit conversion.

<figure>
  <img src="../images/HyShot-field.png" alt="HyShot II combustor: temperature, pressure and hydrogen mass fraction">
  <figcaption>Reacting flow in the computed domain: static pressure, static temperature, and hydrogen mass fraction.</figcaption>
</figure>

---

## Problem setup

### Geometry

The combustor is a constant-area duct of rectangular section, fed by the intake, with a spanwise row of wall portholes injecting hydrogen transversely, and a diverging nozzle at the exit. Only **one half injector pitch** is computed: $y = 0$ is the symmetry plane through a porthole axis and $y = W$ the plane midway to its neighbour, so the row is represented by two symmetry planes rather than by repeating the geometry.

| Parameter | Value | In combustor coordinates |
|---|---|---|
| Combustor leading edge | $x = 350$ mm | $0$ |
| Inflow plane (domain inlet) | $x = 360$ mm | $10$ mm |
| Injector axis | $x = 408$ mm | $58$ mm |
| End of the constant-area duct | $x = 650$ mm | $300$ mm |
| Nozzle exit (domain outlet) | $x = 760$ mm | $410$ mm |
| Duct height $H$ | 9.8 mm | |
| Spanwise half-pitch $W$ | 9.375 mm | injector pitch 18.75 mm |
| Porthole diameter | 2 mm | pipe 2 mm below the wall |
| Nozzle ceiling angle | 12$^\circ$ | bottom wall stays flat |

The porthole is modelled as the half-disk on the $y \ge 0$ side of the symmetry plane, extruded 2 mm below the wall so the jet issues from a short pipe rather than from a prescribed surface state.

### Inflow

The combustor entrance state is **not uniform**, and this matters: the case is run with $p$, $T$ and $u$ profiles imposed point-by-point across the 9.8 mm duct height, rather than with a single freestream state. The profiles are **digitised from the figures of Karl's thesis**.

The digitisation deliberately stops short of both walls. The near-wall region is subsonic and the profiles there are far too steep to read reliably off a printed figure, so rather than extract numbers that cannot be trusted, those points are left out and padded with the sonic state.

| Quantity | Core ($z = 5\ldots8$ mm) | Mass-flux averaged |
|---|---|---|
| Mach number | 2.648 | 2.526 |
| Static temperature | 1303 K | 1366 K |
| Static pressure | 120.0 kPa | 133.6 kPa |
| Velocity | 1860 m/s | 1812 m/s |
| Mass flow per unit width | - | 5.837 kg/(s m) |

Composition is frozen air, $Y_{\mathrm{N_2}} = 0.7686$, $Y_{\mathrm{O_2}} = 0.2314$.

### Fuel injection

Hydrogen is injected through the porthole as a choked pipe flow, prescribed as a mass flux with a stagnation temperature:

| Parameter | Value |
|---|---|
| Species | pure H$_2$ ($Y_{\mathrm{H_2}} = 1$) |
| Stagnation temperature | 300 K |
| Mass flux | 297.13 kg/(m$^2$ s) |
| H$_2$ mass flow (computed half-port) | $4.664\times10^{-4}$ kg/s |
| Air mass flow (computed half-pitch) | $5.47\times10^{-2}$ kg/s |
| Global equivalence ratio | $\boldsymbol{\phi = 0.29}$ |

$\phi$ follows from the two boundary conditions and the stoichiometric H$_2$/air mass ratio of 0.029185; an independent mass-flow audit integrated over the block interfaces of the converged field returns $\phi = 0.2926$, with mass closure across the injection region to 0.12%.

### Boundary conditions

| Boundary | `bc.txt` | Condition |
|---|---|---|
| Inlet ($x = 360$ mm) | `405` | supersonic inflow along the face normal, profiled $M(z)$, $p(z)$, $T(z)$ |
| Injector pipe inlet ($z = -2$ mm) | `403` | mass-flux inlet, $g = 297.134$ kg/(m$^2$ s), $T_0 = 300$ K, pure H$_2$ |
| Bottom wall and ceiling | `302` | isothermal no-slip, $T_w = 300$ K |
| Pipe wall, below the wall plane | `301` | adiabatic no-slip, heat flux $q_w = 0$ |
| $y = 0$ and $y = W$ | `300` | symmetry |
| Outlet ($x = 760$ mm) | `406` | outlet with `p_back` = 0 and `rel_fac` = 1: where the gas leaves, the face takes the state of its boundary cell at any Mach number; where the gas of that cell moves into the domain, the face is a slip wall |
| Block-to-block | `101`, `102` | conformal connection and chimera overlap |

The pipe wall upstream of the wall plane is adiabatic and the part above it is a connection, switched by a range-based composite patch.

---

## Results

<figure>
  --8<-- "docs/vv/images/HyShot-pw.svg"
  <figcaption>Wall pressure along the between-injector line.</figcaption>
</figure>

The measured pressure rises through a shock train from about 120 kPa at the entrance to nearly 290 kPa at the end of the constant-area duct, then collapses to about 50 kPa through the nozzle. MOSE reproduces the whole structure: the level and slope of the rise, the position of the collapse at $x = 650$ mm, and the nozzle expansion.

### Quantitative agreement

Both solutions are interpolated onto the 25 measurement stations and compared against the measured value and its band. The band is wide - its half-width averages 12% of the measurement - which is what the shock train does to a wall-pressure measurement.

| Region | Stations | Experiment | MOSE | TAU S-A |
|---|---|---|---|---|
| Upstream of the injector, $x < 0.408$ m | 1 | 118.4 kPa | 122.4 kPa (+3.4%) | 126.5 kPa (+6.8%) |
| Combustor, $0.41 \le x \le 0.65$ m | 19 | 211.7 kPa | 208.6 kPa (**-1.5%**) | 217.8 kPa (+2.9%) |
| Nozzle, $x > 0.65$ m | 5 | 122.7 kPa | 115.8 kPa (-5.6%) | 118.8 kPa (-3.1%) |

| Pointwise, all 25 stations | MOSE | TAU S-A |
|---|---|---|
| Mean signed deviation | **+0.4%** | +4.8% |
| Mean absolute deviation | **9.6%** | 11.3% |
| RMS deviation | **13.2%** | 17.8% |
| Points inside the experimental band | 17 / 25 | 18 / 25 |

MOSE carries **no systematic bias** in the mean combustor pressure level - the quantity a combustor design actually depends on - and it is closer to the measurement than the reference TAU solution on every aggregate metric except the raw count of points inside the band, where the two are equivalent.

The residual pointwise scatter is dominated by **shock-train phase**, not by level. Both computations produce oscillations of the right amplitude at very nearly the right wavelength, but a small streamwise offset between computed and measured shock positions turns into a large pointwise error whenever a measurement station happens to fall between a computed peak and trough. That the two independent RANS solutions disagree with **each other** by 9.8% RMS - less than either differs from the measurement (13.2% and 17.8%) - confirms the scatter is phase, not physics: a shock train in a duct is exquisitely sensitive to the upstream boundary layer and to the effective back pressure, and neither is measured.

---
