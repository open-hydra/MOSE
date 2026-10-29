# Input Parameters

MOSE reads `input.ini` again during the run, every `ini-diter` iterations. The column During the run says what such a reading does with each key: *taken*, the run takes the value written; *fixed*, the run keeps the value that it read at its start, and another value written is an error of the reading, which `ini-error` decides; *start only*, the key serves the start of the run, and a reading leaves it as it was.

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
| cfl | 0.5 | > 0 | yes | taken | CFL number of the convective spectral radii summed over the directions (RK2 <= 1, RK3 <= 1.17) |
| vnn | 0.3 | > 0 |  no | taken | Diffusion number of the largest diffusivity summed over the directions (RK2 <= 0.5, RK3 <= 0.62) |
| cfl-rise-threshold | 0 | >= 0 |  no | taken | CFL rise threshold |
| time-accurate | .false. | logical | yes | fixed | Time accurate switch |
| integration-variables | cons | cons, prim, prec |  no | fixed | Integration variables (cons/prim/prec) |
| irs | .false. | logical |  no | taken | Implicit Residual Smoothing (steady runs only; a run started without it cannot switch it on) |
| irs-beta | 0.0 | >= 0 |  no | taken | IRS beta parameter (> 0 also enables IRS; steady runs only) |
| irs-variables | flow | flow, all |  no | taken | Residuals smoothed by IRS: flow (species, momentum, energy) or all (also the RANS, soot and passive-scalar equations; for Spalart-Allmaras: the k-omega models diverge with all) |
| preconditioning-Uref | 0.0 | >= 0 |  no | taken | Reference velocity for preconditioning |
| preconditioning-Mach | 0.0 | >= 0 |  no | taken | Mach target for preconditioning |
| preconditioning-eps-min | 0.05 | in (0, 1) |  no | taken | Low-Mach cutoff for Ur |
| preconditioning-Ur-min | 0.0 | >= 0 |  no | taken | Floor on Ur [m/s] |
| preconditioning-Ur-factor | 1.0 | > 0 |  no | taken | Weiss-Smith multiplicative factor |
| preconditioning-Ur-smooth | 0 | >= 0 |  no | taken | Ur max-smoothing passes |
| space-reconstruction |  | MUSCL, first-order | yes | fixed | Space reconstruction method |
| flux-limiter |  | vanalbada, minmod, superbee, vanleer, mc |  no | fixed | Flux limiter for MUSCL reconstruction (van Leer if omitted; ignored with first-order) |
| diffusive-metric | stencil | stencil, mean |  no | fixed | Face metric of the viscous fluxes: stencil takes the face Jacobian from the same differences of the cell centres as the face gradient, exact for a linear field on any mesh (at the symmetry faces and the k faces of type 200 for the linear fields that they allow) except at the chimera faces (102), which take the inverse of the mean of the inverse metrics of the two cells with every value, and at the symmetry faces with one cell across them, the walls 502 and the walls 503-505 while they ablate, which take the metric of the cell; mean takes the mean of the two metrics (the faces between blocks and of type 200 the inverse of the mean of the inverse metrics), overstates the gradient by (1+r)^2/(4r) at a cell-size ratio r and does not converge next to a symmetry plane (300) that bounds a sector of a polar mesh |
| shock-detector |  | Tramel, Chen |  no | fixed | Shock detector method |
| riemann-solver | HLLC | HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC, HLLE, HLLE++, SLAU, SLAU2, LMRoe, MiczekRoe, AUSM+, AUSM+M, LLF, Rusanov, exact |  no | fixed | Riemann solver |
| riemann-options-Mco | 0.0 | >= 0 |  no | taken | Low-Mach acoustic-dissipation cutoff Mach (floor) |

## MOSE-Multigrid

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| level1-iter | 0 | >= 0 |  no | taken | Iterations for multigrid level 1 |
| level2-iter | 0 | >= 0 |  no | taken | Iterations for multigrid level 2 |
| prolongated-diter | 0 | >= 0 |  no | taken | Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write) |

## MOSE-Physics

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| equations |  | euler , navier-stokes | yes | fixed | Gas dynamics equations |
| turbulence | none | SA, SA-R, SA-RC, SAcomp, SST, SST-RC, Wilcox2006, none |  no | fixed | RANS turbulence model |
| chemistry | frozen | frozen, finite-rate, equilibrium |  no | fixed | Chemistry model |
| soot-generation | none | LL91, LIN, none |  no | fixed | Soot generation model |
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
| omega-wall-bc | practical | practical, asymptotic |  no | fixed | Omega wall condition for k-omega models: Menter practical (60*nu/(beta1*y^2), beta1 of SST, for both models) or asymptotic (6*nu/(beta*y^2), exact y->0 limit with the beta of the model, beta1 for SST and beta0 for Wilcox2006, for wall-resolved meshes) |
| sst-production | incompressible | incompressible, compressible |  no | fixed | SST production: incompressible (P = mu_t*S^2; NASA TMR SST-2003m, or SST-2003s with k-energy-coupling) or compressible (P = tau_ij*du_i/dx_j as OpenFOAM kOmegaSST arranges it: the limiter on mu_t*(S^2 - 2/3*divU^2) only, then -2/3*rho*k*divU added, and -2/3*gamma*rho*omega*divU in the omega equation; TMR SST-2003e limits the whole P and takes gamma/nu_t times it in the omega equation) |
| k-energy-coupling | .false. | logical |  no | fixed | Two-equation models (SST, Wilcox2006): k in the total energy, the full compressible RANS form of NASA TMR: -2/3*rho*k*delta_ij in the stress tensor of the momentum and energy equations, and the k-equation source (production - dissipation) taken from the mean-flow energy. .false. = TMR m/e variants |

## MOSE-Rotating-Frame

| Parameter | Default | Allowed | Required | During the run | Description |
|-----------|---------|---------|----------|----------------|-------------|
| omega | 0.0 | >= 0 |  no | taken | Angular speed [rad/s] |
| axis | 0.0 0.0 1.0 |  |  no | taken | Rotation axis direction (3 components) |
| origin | 0.0 0.0 0.0 |  |  no | taken | Point on the rotation axis [m] |
