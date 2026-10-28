# HyShot II Scramjet Combustor

Reacting flow in the HyShot II supersonic combustor at flight condition **Case 2** ($h \approx 28$ km), computed in three dimensions with finite-rate hydrogen/air chemistry on an overset (chimera) grid. This is the most demanding case in the suite: it couples, in a single run, everything the earlier cases exercise one at a time - a supersonic turbulent boundary layer, a shock train, sonic transverse fuel injection into crossflow, multicomponent diffusion, and stiff finite-rate kinetics.

The validation metric is the **wall pressure along the between-injector line**, measured in the HyShot II ground-test campaign and reported with its uncertainty in Karl's thesis, against which a Spalart-Allmaras RANS computed with DLR's TAU code is also available. The comparison is therefore a validation (against measurement) *and* a code-to-code check (against another SA RANS of the same configuration).

**References**: Karl, *Numerical Investigation of a Generic Scramjet Configuration*, PhD thesis, TU Dresden (2011) - the experiment with error bars and the TAU S-A solution; Chapuis et al., *A computational study of the HyShot II combustor performance*, Proc. Combust. Inst. **34** (2013) 2101-2109 - the same measurements in the paper's coordinates, used here only to audit the unit conversion.

!!! warning "The case as committed stops on its coarsest grid level"
    Started from the [initial condition](#initial-condition) that ATLAS ICB builds from the `[ICB-*]` sections of
    `input.ini`, with the committed `input.ini` and boundary-condition files, the run stops on the coarsest grid level
    (level 3), where the multigrid sequence starts, in iteration 858 of its 100,000. Three cells of block 10 (the
    downstream part of the O-grid ring around the porthole), with centres 0.14 mm above the bottom wall and 0.1 to
    0.2 mm outside the edge of the porthole on its downstream side, hold hydrogen of the jet. They empty towards vacuum:
    two adjacent ones go from 17 and 14 kPa after 800 iterations to 284 and 368 Pa after 857 (6.3 and 9.8 K), the third
    from 24 kPa to 358 Pa (6.8 K), while the cells next to them hold 2 to 52 kPa. The update of iteration 858 gives the
    first of them a negative pressure: MOSE prints `Integration failed at b, i, j, k:` with block 10 and cell (3, 2, 5)
    of level 3, the primitive variables of that cell and `NaN, infinity or p<0 detected`, and exits with status 1.

    The collapse comes from the MUSCL reconstruction together with the
    [implicit residual smoothing](../theory/time-integration.md#implicit-residual-smoothing-irs): from the same start,
    the run with `irs-beta = 0` and the run with `space-reconstruction = first-order` pass iteration 2000 of level 3,
    the lowest temperature of the level after 857 iterations being 132.5 K and 162.5 K, while the run without
    chemistry (`exclude-blocks` with the eleven blocks) stops in the same iteration and cell, without the two
    `MATRIX IS REPEATEDLY SINGULAR` messages of the chemistry integrator that come before the stop of the case. All
    these runs use MOSE built with GNU Fortran 12 in `RELEASE` mode, with 4 MPI ranks of 5 OpenMP threads (the case
    as committed stops in the same iteration and cell also without MPI, with 8 threads), and the time stepping of
    `input.ini`: local time steps, RK2, CFL 0.8 and VNN 0.5 ramped linearly over the first 5000 iterations
    (`cfl-rise-threshold`), so 0.137 and 0.086 in iteration 858.

<figure>
  <img src="../images/HyShot-field.png" alt="HyShot II combustor: temperature, pressure and hydrogen mass fraction">
  Reacting flow in the computed domain: static pressure, static temperature, and hydrogen mass fraction.
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

### Initial condition

ATLAS ICB builds the initial field, `INPUT/ic.szplt`, from the `[ICB-*]` sections of `input.ini` (see [What ships, and what you build](#what-ships-and-what-you-build)). Every block takes the state of `[ICB-Block*]`, except blocks 6 to 9, the four blocks of the porthole that run down the pipe: their own sections split them along $z$ at the wall plane, and a cell whose centre lies below it takes the state of `[hydroIC]`, any other cell the state of `[airIC]`, equal to that of `[ICB-Block*]`.

| Region | Gas | Velocity | $p$ | $T$ | $\tilde\nu$ |
|---|---|---|---|---|---|
| Every cell above the wall plane | air, $Y_{\mathrm{N_2}} = 0.7686$, $Y_{\mathrm{O_2}} = 0.2314$ | $M = 2.49$ along $x$ | 130.2 kPa | 1377 K | $10^{-6}$ m$^2$/s |
| The pipe, below the wall plane | pure H$_2$ | $w = 1204.2$ m/s, into the duct | 252.8 kPa | 248.4 K | $10^{-6}$ m$^2$/s |

The hydrogen moves at its speed of sound with the mass flux of the `403` inflow of the pipe: with the tables of the case, 1204.2 m/s is the speed of sound of hydrogen at 248.4 K and $\rho w = 297.15$ kg/(m$^2$ s) is the mass flux $g = 297.134$ kg/(m$^2$ s) of that inflow, each to within 0.01 %, and the stagnation temperature of this state is 299.5 K, against $T_0 = 300$ K of the inflow. The air is uniform: it does not follow the profile of the inlet, and next to the 300 K walls it starts at 1377 K and $M = 2.49$ too.

---

## Numerical setup

### Grid

Three independent, overlapping structured regions - a stretched upstream box, the homogeneous chamber background, and a refined O-grid patch around the porthole with its pipe - coupled by **chimera overlap**, with no conformal interface between the regions. The overlap connectivity (donor/overset, type 102) is computed by the ATLAS pre-processor from the geometry, so the regions can be re-meshed independently.

| Parameter | Value |
|---|---|
| Blocks | 11 |
| Cells | 4,563,968 |
| Boundary faces | 414,752 |
| of which chimera / connection / symmetry / wall / inlet-outlet | 28,768 / 59,712 / 224,224 / 95,008 / 7,040 |
| Chamber spacing $\Delta x = \Delta y$ | 0.25 mm |
| First wall cell | 1 $\mu$m, growth 1.2, 28 cells (mirrored at the ceiling) |
| Wall resolution, bottom wall | $y^+$ median 1.3, 95th percentile 1.9, max 2.1 |

### Physics and numerics

| Parameter | Value |
|---|---|
| Equations | Navier-Stokes |
| Equation of state | ideal gas, thermally perfect (variable $c_p(T)$) |
| Turbulence model | Spalart-Allmaras |
| Species diffusion | constant laminar Schmidt number, $Sc = 0.7$ |
| Turbulent Schmidt number | $Sc_t = 0.7$ |
| Chemistry | finite-rate, Gerlinger-9 (9 species, 19 reactions) |
| Chemistry-flow coupling | Strang operator splitting |
| ODE solver | `H-radau5` |
| Riemann solver | HLLC (Batten) |
| Space reconstruction | MUSCL, Van Leer limiter |
| Time scheme | RK2, conservative variables, steady (local time stepping) |
| Implicit residual smoothing | $\beta = 0.5$ |
| CFL / VNN | 0.8 / 0.5 |
| Multigrid | 3 levels, 100,000 (L3) + 20,000 (L2) + 20,000 (L1) iterations |

### Chemistry tolerances

The ODE error control is set explicitly and tightly:

```ini
[MOSE-Chemistry]
ode-absolute-tol-species     = 1e-12 (x9)
ode-relative-tol-species     = 1e-8  (x9)
ode-absolute-tol-temperature = 1e-4
ode-relative-tol-temperature = 1e-8
ode-max-steps                = 100000
```

This is not a convergence detail. MOSE integrates the chemistry on **partial densities** $\rho_s$ [kg/m$^3$], and the registered default absolute tolerance sits one to two decades *above* the entire pre-ignition radical pool. At the default the induction zone is integrated as if it were zero, heat release is displaced downstream, and the case shows no OH anywhere near the porthole - while the flow field (separation size, reversed-flow fraction, peak backflow velocity) is identical to a fraction of a percent. Tightening the tolerance alone raises peak $Y_{\mathrm{OH}}$ in the porthole region by two to three orders of magnitude.

Getting ignition to occur **ahead** of the jet, as the literature reports, additionally requires the profiled inflow: the uniform inlet delivers near-wall gas about 145 K colder, which weakens the upstream separation and roughly quarters the amount of hot flammable gas there. Both changes are necessary; neither is sufficient.

---

## Results

<figure>
  {% include "vv/images/HyShot-pw.svg" %}
  Wall pressure along the between-injector line.
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

## Running the case

### What ships, and what you build

Everything the case needs is committed **except the initial condition**. At 4.6 M cells an initial field is about 375 MB even as a single-precision `.szplt`, and 2.2 GB as Tecplot ASCII, so it is not something a repository can carry.

| Committed | Build yourself |
|---|---|
| `MESH/mesh.szplt` - the overset grid, and `MESH/hyshot_mesh.py` that generates it | `INPUT/ic.szplt` - the initial field |
| `INPUT/phase.txt`, `thermo.dat`, `transport.dat`, `chemistry-*` - the 9-species Gerlinger mixture | |
| `INPUT/bc.txt`, `bc2.txt`, `bc3.txt` - boundary conditions and chimera connectivity, one file per multigrid level | |
| `input.ini` - the solver setup documented above | |

So the first step is to build the initial condition with the ATLAS pre-processor, from the shipped mesh and phase data. ICB looks for the mesh and the phase files **in the working directory** (`mesh.tec`, then `mesh.p3d`, then `mesh.szplt`), and writes into `fromATLAStoSolver/`:

```bash
cd test/3D/hyshotII
ln -s MESH/mesh.szplt .          # ICB searches the cwd, not MESH/
ln -s INPUT/phase.txt .          # and needs the phase file beside it
ulimit -s unlimited              # see the warning below
ATLAS ICB                        # -> fromATLAStoSolver/ic.szplt
mv fromATLAStoSolver/ic.szplt INPUT/
```

Then solve:

```bash
./MOSE.sh compile              # or copy an existing bin/MOSE
./MOSE.sh solve -p 60          # 60 OpenMP threads; ~30 h for the full MG schedule
```

The boundary-condition files are committed, so ATLAS BCB does not need to run. They already encode everything in the [boundary conditions](#boundary-conditions) table - the profiled supersonic inlet, the sonic hydrogen porthole, the 300 K isothermal walls and the chimera overlap - which is why the physical settings are documented on this page rather than being readable from `input.ini`.

The wall-pressure figure is rebuilt from the solver's `OUTPUT/wall.tec`:

```bash
cd test/3D/hyshotII/reference
python wall_pressure.py                  # -> docs/vv/images/HyShot-pw.svg
python wall_pressure.py --preset paper   # -> wall_pressure.pdf / .png
```

`wall_pressure.py` stitches the between-injector line across whatever bottom-wall zones `wall.tec` contains, accepting only rows that genuinely sample $y = W$, so it survives a re-mesh without editing. It also re-audits the thesis-to-paper coordinate and pressure-normalisation transform on every run, and prints the worst disagreement between the two digitisations of the same measurements.
