# Input Parameters


## MOSE-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false |  no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | Residual convergence threshold |
| time-threshold | 1e30 | > 0 |  no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 |  no | Maximum number of iterations |

## MOSE-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ini-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | Initial condition (INPUT/ic.*) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | Solution (OUTPUT/field.*) format |
| sol-diter | 1000000000 | > 0 |  no | Solution output iter frequency |
| sol-dtime | 1e30 | > 0 |  no | Solution output time frequency |
| sol-overwrite | true | true, false |  no | Overwrite solution files |
| sol-variables | thermo |  |  no | Solution variables to write |
| wall-variables | mech thermal |  |  no | Wall variables to write |
| res-diter | 1 | > 0 |  no | Residual history iter frequency |
| shell-diter | 1 | > 0 |  no | Shell update iter frequency |
| ini-diter | 10000 | > 0 |  no | input.ini update iter frequency |
| timer-diter | 0 | >= 0 |  no | Wall-clock timing report iter frequency (0 = off) |

## MOSE-Probes

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| probe1 | probe-placeholder |  |  no | Probe file name |

## probe-placeholder

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| variables | none |  |  no | Probe variables to write |
| dtime | 1e30 | > 0 |  no | Probe output time frequency |
| diter | 1000000000 | > 0 |  no | Probe output iter frequency |
| index-position | 0 0 0 0 | >= 0 |  no | Probe location by index |
| position | 0.0 0.0 0.0 |  |  no | Probe location by coordinates |

## MOSE-Numerics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| time-scheme | euler | euler, RK2, RK3 | yes | Time integration solver |
| cfl | 0.5 | > 0 | yes | CFL number |
| vnn | 0.3 | > 0 |  no | VNN parameter |
| cfl-rise-threshold | 0 | >= 0 |  no | CFL rise threshold |
| time-accurate | .false. | logical | yes | Time accurate switch |
| integration-variables | cons | cons ,  prim , prec |  no | Integration variables (cons/prim/prec) |
| irs | .false. | logical |  no | Implicit Residual Smoothing |
| irs-beta | 0.0 | >= 0 |  no | IRS beta parameter |
| preconditioning-Uref | 0.0 | >= 0 |  no | Reference velocity for preconditioning |
| preconditioning-Mach | 0.0 | >= 0 |  no | Mach target for preconditioning |
| preconditioning-eps-min | 0.05 | in (0, 1) |  no | Low-Mach cutoff for Ur |
| preconditioning-Ur-min | 0.0 | >= 0 |  no | Floor on Ur [m/s] |
| preconditioning-Ur-factor | 1.0 | > 0 |  no | Weiss-Smith multiplicative factor |
| preconditioning-Ur-smooth | 0 | >= 0 |  no | Ur max-smoothing passes |
| space-reconstruction |  | MUSCL, first-order | yes | Space reconstruction method |
| flux-limiter |  | vanalbada, minmod, superbee, vanleer, mc |  no | Flux limiter for space reconstruction |
| shock-detector |  | Tramel, Chen |  no | Shock detector method |
| riemann-solver | HLLC | HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC, HLLE, HLLE++, SLAU, SLAU2, LMRoe, MiczekRoe, AUSM+, AUSM+M, LLF, Rusanov, exact |  no | Riemann solver |
| riemann-options-Mco | 0.0 | >= 0 |  no | Low-Mach acoustic-dissipation cutoff Mach (floor) |

## MOSE-Multigrid

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| level1-iter | 0 | >= 0 |  no | Iterations for multigrid level 1 |
| level2-iter | 0 | >= 0 |  no | Iterations for multigrid level 2 |

## MOSE-Physics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| equations |  | euler , navier-stokes | yes | Gas dynamics equations |
| turbulence | none | SA, SA-R, SA-RC, SAcomp, SST, Wilcox2006, none |  no | RANS turbulence model |
| chemistry | frozen | frozen, finite-rate, equilibrium |  no | Chemistry model |
| soot-generation | none | LL91, LIN, none |  no | Soot generation model |
| rotational-frame | none | rigid-body, none |  no | Rotational frame model |
| schmidt | 0.0 | >= 0 |  no | Laminar Schmidt number |
| prandtl | 0.0 | >= 0 |  no | Laminar Prandtl number |

## MOSE-Chemistry

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| exclude-blocks | none |  |  no | Blocks to exclude from chemistry |
| ode-solver | H-radau5 | H-radau5, sdirk4b, ros4 |  no | ODE solver for chemistry |
| ode-max-steps | 100000 | > 0 |  no | Maximum ODE integration steps |
| ode-analytical-jacobian | .false. | logical |  no | Use the mechanism analytical Jacobian instead of finite differences |
| ode-relative-tol-species | 1e-5 | > 0 |  no | ODE relative tolerance for species |
| ode-relative-tol-temperature | 1e-5 | > 0 |  no | ODE relative tolerance for temperature |
| ode-absolute-tol-species | 1e-5 | > 0 |  no | ODE absolute tolerance for species |
| ode-absolute-tol-temperature | 1e-5 | > 0 |  no | ODE absolute tolerance for temperature |

## MOSE-Turbulence

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| Prt | 0.85 | > 0 |  no | Turbulent Prandtl number |
| Sct | 0.90 | > 0 |  no | Turbulent Schmidt number |
| point-implicit | .true. | logical |  no | Point-implicit (Patankar) treatment of turbulence destruction source terms |
| omega-wall-bc | practical | practical, asymptotic |  no | Omega wall condition for k-omega models: Menter practical (60*nu/(beta1*y^2)) or asymptotic (6*nu/(beta1*y^2), exact y->0 limit, for wall-resolved meshes) |

## MOSE-Rotating-Frame

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| omega | 0.0 | >= 0 |  no | Angular speed [rad/s] |
| axis | 0.0 0.0 1.0 |  |  no | Rotation axis direction (3 components) |
| origin | 0.0 0.0 0.0 |  |  no | Point on the rotation axis [m] |
