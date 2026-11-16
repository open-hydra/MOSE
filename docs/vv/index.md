# Verification & Validation

This section shows that MOSE solves its equations correctly (verification) and that those equations describe real flows (validation). Each page states what is checked, against which reference, and with what result.

The cases are grouped from the numerics to the applications. If you are new to the section, start with **Order of accuracy**, which shows that the scheme converges at its design order, then move to the group that matches your flow.

## Order of accuracy

Exact or manufactured solutions refined on successive grids: the error must fall at the design order of the scheme.

| Case | What it verifies | Reference |
|---|---|---|
| [Isentropic Vortex](2D-isentropic-vortex.md) | Order of the inviscid scheme on a smooth flow | Exact solution (Shu 1998) |
| [Viscous Manufactured Solution](2D-mms.md) | Order of the viscous fluxes, gradients and metrics | Manufactured solution (Roache 2002) |
| [Manufactured Solutions on Skewed Grids](2D-mms-skewed.md) | The same on skewed, non-smooth two-block grids with walls | Manufactured solution (Salari and Knupp 2000) |
| [3D Manufactured Solutions on Skewed Grids](3D-mms.md) | The same in three dimensions | Manufactured solution (Salari and Knupp 2000) |

## Shock capturing

Discontinuities and strong waves against exact solutions.

| Case | What it verifies | Reference |
|---|---|---|
| [1D Euler Problems](1D.md) | Shocks, contacts, rarefactions, near-vacuum, strong shocks, shock–entropy-wave interaction (Sod, Einfeldt, Noh, Toro, Shu–Osher) | Exact Riemann solutions; grid convergence |
| [Oblique Shock](2D-oblique-shock.md) | Mach 4 oblique shock, post-shock oscillations of four Riemann solvers | Oblique shock relations |
| [Woodward–Colella Step](2D-woodward-colella.md) | Mach 3 step: shock reflections and interactions | OpenFOAM |
| [Hypersonic Cylinder](2D-hypersonic-cylinder.md) | Mach 8.1 bow shock, carbuncle | Normal shock relations |

## Low Mach number

The accuracy of a compressible solver as the Mach number goes to zero.

| Case | What it verifies | Reference |
|---|---|---|
| [Gresho Vortex](2D-gresho-vortex.md) | Low-Mach dissipation of three Riemann solvers, Mach 0.1 to 0.001 | Steady exact solution |
| [Low-Mach Cylinder](2D-low-mach-cylinder.md) | Preconditioning: Mach-independent solution at Mach 0.01 and 0.001 | Potential flow |

## Wall-bounded and turbulent flows

Boundary layers, turbulence models and shock–boundary-layer interaction.

| Case | What it verifies | Reference |
|---|---|---|
| [Laminar Flat Plate](2D-flat-plate-laminar.md) | Laminar boundary layer | Blasius solution |
| [Turbulent Flat Plate of the NASA TMR](2D-flat-plate-turbulent-tmr.md) | SA, SST, Wilcox 2006 and SSG/LRR-ω on nested grids | CFL3D, FUN3D, TAU (NASA TMR) |
| [Planar Nozzle on Grid Families](2D-planar-nozzle-grids.md) | Turbulent nozzle flow independent of grid direction and smoothness | Mirror symmetry; quasi-1D flow |
| [Shock Wave–Boundary Layer Interaction](2D-swbli.md) | Mach 5 separation and reattachment with SA and SST | OpenFOAM on the same grid; experiment, Wind-US, SU2 |

## Multicomponent and reacting flows

Species diffusion, chemistry and gas–surface interaction.

| Case | What it verifies | Reference |
|---|---|---|
| [Unity-Lewis Diffusion](1D-diffusion.md) | Species and thermal diffusion with $Le = 1$ | Cantera |
| [Multicomponent Diffusion](1D-multicomponent-diffusion.md) | Mixture-averaged diffusion | Cantera |
| [Premixed Laminar Flame](1D-premixed-flame.md) | Flame structure with finite-rate chemistry | Cantera |
| [Reactive Shock Tube](1D-fer14.md) | Shock-induced H$_2$/O$_2$ chemistry | Ferrer et al. (2014) |
| [Rocket Nozzle](2D-nozzle.md) | Frozen, equilibrium and finite-rate expansion | NASA CEA |
| [Ablating Wall](2D-ablating-wall.md) | Melting and pyrolysing walls: surface mass and energy balance | Closed-form surface balance |

## Applications

| Case | What it verifies | Reference |
|---|---|---|
| [HyShot II Scramjet Combustor](3D-hyshotII.md) | 3D reacting scramjet: fuel injection, shock train, chemistry, overset grids | Experiment; TAU |

## Solver features

| Case | What it verifies | Reference |
|---|---|---|
| [Rotating Frame](rotating-frame.md) | Coriolis and centrifugal source terms | Analytical balance |
| [Multirate Time Stepping](multirate/index.md) | Local time stepping: exactness, conservation, temporal order, accuracy on real flows | Single-rate solver; exact and published solutions |

## Running the cases

Each case lives in a folder under `test/`, whose `README.md` gives its setup and how to run it. Most folders hold a `verify.py` that reads the solver output, compares it with the reference and draws the figures of these pages; `test/TESTING-PLAN.md` lists every test and what it checks.
