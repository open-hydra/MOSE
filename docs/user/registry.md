# Input Parameters

MOSE reads `input.ini` again every `ini-diter` iterations. The column *During the run* says what such a reading does with each key:

- **taken**: the run takes the value written.
- **fixed**: the run keeps the value it read at its start; another value is an error of the reading, which `ini-error` decides.
- **start only**: the key serves the start of the run; a reading leaves it as it was.

## MOSE-Parameters

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| newrun | true | true , false |  no | start only | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | taken | Residual convergence threshold, relative to the first iteration (steady runs only) |
| time-threshold | 1e30 | > 0 |  no | taken | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 |  no | taken | Maximum number of iterations |

## MOSE-IO

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| ic-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | start only | Initial condition (INPUT/ic.*) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | fixed | Solution (OUTPUT/field.*) format |
| decomposition-map | none |  |  no | start only | Block-to-rank map in INPUT/ (none: MOSE assigns the blocks) |
| sol-diter | 1000000000 | > 0 |  no | taken | Solution output iter frequency |
| sol-dtime | 1e30 | > 0 |  no | taken | Solution output time frequency |
| sol-overwrite | true | true, false |  no | taken | Overwrite solution files |
| sol-variables | thermo |  |  no | fixed | Solution variables to write |
| wall-variables | mech thermal |  |  no | fixed | Wall variables to write |
| res-diter | 1 | > 0 |  no | taken | Residual history iter frequency |
| shell-diter | 1 | > 0 |  no | taken | Shell update iter frequency |
| ini-diter | 10000 | > 0 |  no | taken | input.ini update iter frequency |
| ini-error | warn | warn, stop |  no | taken | input.ini read again with an error: warn and keep the values held before the reading, or stop (at set-up every error stops) |
| timer-diter | 0 | >= 0 |  no | taken | Wall-clock timing report iter frequency (0 = off) |

## MOSE-Probes

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| probe1 | probe-placeholder |  |  no | fixed | Probe file name |

## probe-placeholder

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| variables | none |  |  no | fixed | Probe variables to write |
| dtime | 1e30 | > 0 |  no | fixed | Probe output time frequency |
| diter | 1000000000 | > 0 |  no | fixed | Probe output iter frequency |
| index-position | 0 0 0 0 | >= 0 |  no | fixed | Probe location by index |
| position | 0.0 0.0 0.0 |  |  no | fixed | Probe location by coordinates |

## MOSE-Numerics

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| time-scheme | euler | euler, RK2, RK3 | yes | fixed | Time integration solver |
| cfl | 0.5 | > 0 | yes | taken | CFL number of the convective spectral radii, summed over the directions (RK2 <= 1, RK3 <= 1.17; at most 1 with the multirate integrator) |
| vnn | 0.3 | > 0 |  no | taken | Diffusion number of the largest diffusivity, summed over the directions (RK2 <= 0.5, RK3 <= 0.62; at most 0.5 with the multirate integrator) |
| cfl-rise-threshold | 0 | >= 0 |  no | taken | CFL rise threshold |
| time-accurate | .false. | logical | yes | fixed | Time accurate switch |
| integration-variables | cons | cons, prim, prec |  no | fixed | Integration variables: cons, prim (linearised primitive update) or prec (Weiss-Smith preconditioning; steady runs, with HLLC-PC) |
| irs | .false. | logical |  no | taken | Implicit Residual Smoothing (steady runs only; a run started without it cannot switch it on) |
| irs-beta | 0.0 | >= 0 |  no | taken | IRS beta parameter (> 0 also enables IRS; steady runs only) |
| irs-variables | flow | flow, all |  no | taken | Residuals smoothed by IRS: flow (species, momentum, energy) or all (also the RANS, soot and passive-scalar equations; for Spalart-Allmaras: the k-omega models diverge with all) |
| preconditioning-Uref | 0.0 | >= 0 |  no | taken | Uniform Ur [m/s], a where the flow speed exceeds it, no viscous floor; 0: not used |
| preconditioning-Mach | 0.0 | >= 0 |  no | taken | Ur that gives this preconditioned Mach number, a where the flow is faster; 0: not used |
| preconditioning-eps-min | 0.05 | in (0, 1) |  no | taken | Floor of Ur as a fraction of a: keep it below the Mach number of the flow (about M/10) |
| preconditioning-Ur-min | 0.0 | >= 0 |  no | taken | Floor on Ur [m/s], in every mode |
| preconditioning-Ur-factor | 1.0 | > 0 |  no | taken | Ur = factor times the flow speed in the default (Weiss-Smith) mode |
| preconditioning-Ur-smooth | 0 | >= 0 |  no | taken | Passes of max-smoothing of Ur over the face neighbours |
| space-reconstruction |  | MUSCL, first-order | yes | fixed | Space reconstruction method |
| flux-limiter |  | vanalbada, minmod, superbee, vanleer, mc |  no | fixed | Flux limiter for MUSCL reconstruction (van Leer if omitted; ignored with first-order) |
| diffusive-metric | stencil | stencil, mean |  no | fixed | Face metric of the viscous fluxes: stencil (exact for a linear field on any mesh) or mean (mean of the two cell metrics) |
| shock-detector |  | Tramel, Chen |  no | fixed | Shock detector method |
| riemann-solver | HLLC | HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC, HLLE, HLLE++, SLAU, SLAU2, LMRoe, MiczekRoe, AUSM+, AUSM+M, LLF, Rusanov, exact |  no | fixed | Riemann solver |
| riemann-options-Mco | 0.0 | >= 0 |  no | taken | Low-Mach acoustic-dissipation cutoff Mach (floor) |

## MOSE-Multigrid

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| level1-iter | 0 | >= 0 |  no | taken | Iterations for multigrid level 1 |
| level2-iter | 0 | >= 0 |  no | taken | Iterations for multigrid level 2 |
| prolongated-diter | 0 | >= 0 |  no | taken | Coarse levels: least iterations between two writes of field/wall-prolongated; the end of a level always writes them (0 = every field write) |

## MOSE-Multirate

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| enabled | .false. | logical |  no | fixed | Multirate local time stepping: recursive Heun (SSP-RK2) with a level per cell, time-accurate runs only |
| max-level | 10 | in [0, 50] |  no | fixed | Deepest multirate level allowed: the finest cells advance with H/2^L (fewer levels when the local steps span fewer octaves) |
| safety | 2.0 | >= 1 |  no | fixed | Safety factor sigma: every cell advances with a sub-step of at most dtcap/sigma (dtcap = local CFL/VNN step) |
| macro-step-target | 0.0 | >= 0 |  no | fixed | Target macro step H [s]; 0 = set by the level map |
| levels | recompute | fixed, recompute |  no | fixed | Level map: computed at the first macro step and kept (fixed) or recomputed from the local steps (recompute) |
| recompute-every | 1 | > 0 |  no | fixed | Macro steps between two recomputations of the level map (levels = recompute) |
| fixed-macro-step | 0.0 | >= 0 |  no | fixed | Fixed macro step as a multiple of the dt_min of the first macro step; 0 = adaptive (verification runs) |
| conservative-correction | .true. | logical |  no | fixed | Conservative flux correction at the faces between levels (.false. only as a positive control) |
| trace | .false. | logical |  no | fixed | Print the stage sequence, the level map and the halo lists of the first macro step (standard output) |
| nonlocal-cfl | redo | off, count, redo |  no | fixed | Check of condition (4) of Gnedin et al. at every sub-step: off, count (counts only) or redo (a failing macro step is redone with finer levels or a smaller H) |
| interface-margin | .false. | logical |  no | fixed | Faces between levels: the coarse step stays within dtcap/kappa of the fine cell (always on with nonlocal-cfl = redo) |
| force-redo | .false. | logical |  no | fixed | Test device: the first attempt of every macro step is rejected and redone; the results must not change |
| max-attempts | 21 | >= 1 |  no | fixed | nonlocal-cfl = redo: most attempts of one macro step before the run stops |

## MOSE-Physics

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| equations |  | euler , navier-stokes | yes | fixed | Gas dynamics equations |
| turbulence | none | SA, SA-R, SA-RC, SAcomp, SST, SST-RC, Wilcox2006, SSGLRR, SSGLRR-SD, SSGLRR-w2019, SSGLRR-w2019-SD, none |  no | fixed | RANS turbulence model |
| chemistry | frozen | frozen, finite-rate, equilibrium |  no | fixed | Chemistry model |
| soot-generation | none | LL91, LIN, none |  no | fixed | Soot generation model (under development: not validated) |
| rotational-frame | none | rigid-body, none |  no | fixed | Rotational frame model |
| schmidt | 0.0 | >= 0 |  no | taken | Laminar Schmidt number |
| prandtl | 0.0 | >= 0 |  no | taken | Laminar Prandtl number |

## MOSE-Chemistry

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| exclude-blocks | none |  |  no | fixed | Blocks to exclude from chemistry (finite-rate and equilibrium, on every multigrid level) |
| ode-solver | H-radau5 | H-radau5, H-rodas, H-sdirk4, H-dopri5, radau2a, lobatto3c, gauss, radau1a, ros2, ros3, ros4, rodas3, rodas4, sdirk2a, sdirk2b, sdirk3a, sdirk4a, sdirk4b, dodesol, cvode |  no | fixed | ODE solver for chemistry |
| ode-max-steps | 100000 | > 0 |  no | fixed | Maximum ODE integration steps |
| ode-analytical-jacobian | .false. | logical |  no | fixed | Use the mechanism analytical Jacobian instead of finite differences |
| ode-relative-tol-species | 1e-5 | > 0 |  no | fixed | ODE relative tolerance for species |
| ode-relative-tol-temperature | 1e-5 | > 0 |  no | fixed | ODE relative tolerance for temperature |
| ode-absolute-tol-species | 1e-5 | > 0 |  no | fixed | ODE absolute tolerance for species |
| ode-absolute-tol-temperature | 1e-5 | > 0 |  no | fixed | ODE absolute tolerance for temperature |

## MOSE-Turbulence

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| Prt | 0.85 | > 0 |  no | taken | Turbulent Prandtl number |
| Sct | 0.90 | > 0 |  no | taken | Turbulent Schmidt number |
| point-implicit | .true. | logical |  no | taken | Point-implicit factor 1/(1 + dt d) on the increment of the turbulence equations, steady runs only (time-accurate runs: explicit sources, stiffness in the time step) |
| omega-wall-bc | practical | practical, asymptotic |  no | fixed | Wall value of omega: practical (Menter, 60 nu/(beta1 y^2)) or asymptotic (6 nu/(beta y^2), for wall-resolved meshes) |
| sst-production | incompressible | incompressible, compressible |  no | fixed | SST production term: incompressible (mu_t S^2) or compressible (tau_ij du_i/dx_j, as OpenFOAM kOmegaSST) |
| sst-dilatation-dissipation | none | none, wilcox, zeman |  no | fixed | Compressibility correction of the SST dissipation (Rumsey 2009): none, wilcox (free shear layers) or zeman (boundary layers) |
| k-energy-coupling | .false. | logical |  no | fixed | Two-equation models: k in the total energy and in the stress tensor (full NASA TMR form), .false. = TMR m/e variants; SSGLRR is always in this form |

## MOSE-Rotating-Frame

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| omega | 0.0 | >= 0 |  no | taken | Angular speed [rad/s] |
| axis | 0.0 0.0 1.0 |  |  no | taken | Rotation axis direction (3 components) |
| origin | 0.0 0.0 0.0 |  |  no | taken | Point on the rotation axis [m] |
