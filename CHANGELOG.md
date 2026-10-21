# Changelog

All notable changes to MOSE are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- The key `diffusive-metric` of `[MOSE-Numerics]` chooses the face metric of the viscous fluxes at interior faces:
  `inverse-mean`, the default, the inverse of the mean of the inverse metrics of the two cells, whose face-normal edge is the
  vector between the cell centres; or `mean`, the mean of the two metrics, as before.
- Probes can record `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde`
  (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model).

### Changed

- The interior faces of the viscous fluxes take by default the inverse of the mean of the inverse metrics. The mean of the two
  metrics divided the difference across a face by the mean of 1/h instead of the distance between the cell centres: it
  overstated the face-normal gradient, and with it the viscous, conductive and diffusive fluxes, by (1 + r)²/(4r) at a cell-size
  ratio r (4.2 % at r = 1.5, 8.9 % at r = 1.8), at any resolution. Viscous runs change wherever neighbouring cells differ in
  size; `diffusive-metric = mean` keeps the previous metric.
- Set-up errors stop the run before the initial condition is read, on every MPI rank: the root rank found them only after
  reading the initial condition, and stopped by itself. Errors of the transport data are reported there, such as a mixture of
  more than one species with `schmidt` = 0, its default, without `INPUT/diffusion.dat`, which failed later, in the first step of
  the run, when the diffusive flux read the table that was not there; a single species needs no such table, as it has no species
  diffusion.
- A restart appends to the probe files, and creates those that are missing; it stopped when one was missing.
- Each rank keeps the state of the remote chimera donors of its own receivers only, and none of the remote sources of manifold
  faces (`501`), which it never reads; the chimera exchange sends each donor cell once per rank.
- The diagnostic output of a steady run with a single species needs `sol-variables` to name `thermo` (the default) or
  `transport`; otherwise the run stops at its first solution output with a message that says so, where it named only the
  size of a buffer.
- The local time step adds its limits over the directions and over convection and diffusion:
  1/Δt = Σ_d (|u·∇ξ_d| + c|∇ξ_d|)/cfl + D_max Σ_d |∇ξ_d|²/vnn, where ∇ξ_d are the rows of the cell metric (the inverse widths of
  the cell) and D_max is the largest diffusivity of the cell: momentum, 4/3 (μ + μ_t)/ρ; energy, (k + μ_t c_p/Pr_t)/(ρc_v),
  that is γν/Pr = 1.94ν in air; species, with more than one; Spalart–Allmaras, (μ + ρν̃)/(σρ). In 3-D, on a cell whose three
  directions are skewed against each other (a² + b² + c² + 2abc > 1 for the cosines a, b, c between the ∇ξ_d), the viscous part
  takes the factor (2 + m + 1/m)/4, below 9/8, with m + 1 the largest eigenvalue of the matrix of the cosines: there the sum
  alone does not bound the spectral radius of the viscous operator. It took the smallest of the direction limits, with the edge
  length of the cell, and of a viscous limit on (μ + μ_t)/ρ, which is not stable where neither
  one direction nor one process dominates: at cfl 0.8, sound waves in a box of square cells (with Heun the run stopped on a
  negative pressure; with RK3 the pressure between 0.56 and 1.48 times its initial value) and an entropy wave on a skewed mesh
  (stopped on a negative pressure, with Heun and with RK3), and heat conduction at the default vnn 0.3, which is 0.58 on the
  energy diffusivity, beyond the limit 1/2 of Heun (stopped on a negative pressure). With the same cfl and vnn the step is
  smaller, down to one half on square 2-D cells and by 1/1.94 in cells set by heat conduction: steady runs need more iterations
  (3.4 to 3.9 times as many in diffusion-dominated cases, at the same cost per iteration; 555 instead of 385 in the Prandtl–Meyer
  test, 8000 instead of 4000 in the inert configuration of the ablating-wall test), and time-accurate runs advance by a shorter
  step (0.59 of the global step in a turbulent nozzle case with Spalart–Allmaras, RK3, cfl 0.8, vnn 0.3). Steady solutions do not
  move. The defaults of cfl and vnn do not change; the user guide gives the stable values and how to recover the step of before.
- In time-accurate RANS runs the turbulence sources are explicit, and the local time step includes the stiffness of their
  destruction, so that Δt·d ≤ 1 and no Runge–Kutta stage can make ν̃, k or ω negative through it. The point-implicit factor
  (`point-implicit`, true by default) now applies to steady runs only; a time-accurate run with it prints a warning. The factor
  divided the whole source, production included, by 1 + Δt·d, which made the turbulence equations, and the mean flow through
  μ_t, first-order accurate in time (order 1.00 for ν̃, 0.90 for k and 1.02 for ω with Heun and RK3; now about 2 with Heun and
  3 with RK3). Where the destruction is faster than the transport the local step of those cells is smaller than before; with
  Spalart–Allmaras the global step of a turbulent nozzle case did not change. Steady runs do not change. A time-accurate run with
  SST or Wilcox 2006 and no CFL ramp (`cfl-rise-threshold`) prints a warning, and so does one with a rotation–curvature (-RC)
  correction, whose term is evaluated without its time derivative; at the end of a RANS run one line gives the number of cell
  updates in which a stage made ν̃, k or ω negative before the realizability floor.
- `res-threshold` ends steady runs only; time-accurate runs end on `time-threshold` or `iter-threshold`. The density residual
  that `res-threshold` reads says nothing about the convergence of a transient: a time-accurate run started from a steady state
  stopped after one step. A time-accurate case that ended on `res-threshold`, with the key or with its default, now runs to
  `time-threshold` or `iter-threshold`.

### Removed

- The compressibility correction of Spalart–Allmaras (`turbulence = SAcomp`): it is refused at set-up. It took c from
  c² f(c) = 25 Ω ν_t / a, a number against a speed, so the factor it gave the production depended on the units: 0.99 with
  lengths in metres and 0.22 with lengths in millimetres for Ω = 1e3 1/s, ν_t = 1e-3 m²/s and a = 340 m/s; at Ω = 1e5 1/s it
  cut the production by 78 % whatever the Mach number. No case or test of the repository used it.

### Fixed

- Faces between blocks and chimera faces take the inverse of the mean of the inverse metrics: across a large size jump at a block
  connection the mean of the two metrics made the normal gradient far too steep, and the explicit diffusion unstable on the large
  cell.
- With MPI the tangential neighbours of the source cell of a block connection, which on the edge of a face can be ghost cells
  filled from another rank, are taken after the ghost exchange: they were one stage old, and MPI runs drifted from the serial run
  at block edges and junctions.
- With MPI the output variables of a block are computed on the rank that owns it: on the root rank the eddy viscosity of a
  Spalart–Allmaras run with `sol-variables = transport` read the metrics of the blocks of the other ranks, which the root rank
  does not allocate, and the velocity gradients at block edges read ghost cells that were never filled.
- At chimera faces (`102`) the viscous flux takes the tangential gradient of the donor side, and in RANS runs the chimera ghost
  cells take the wall distance and roughness of the receiver instead of zero.
- Beyond the edges of a chimera block face the viscous stencil takes the row of the ghost cells extrapolated linearly, and in
  the k direction of a 2-D axisymmetric wedge the ghost cell turned by the wedge angle, as the cells do: the ghost cell copied
  there gave the faces at the edges 3/4 of the tangential derivative of a linear field, and on a wedge every chimera face half
  of the v/r of the divergence.
- With a turbulence model the diagnostic output writes the residuals of the turbulence variables and declares as many variables
  as it writes: it declared nres + 2 variables and named and filled seven.
- An unknown probe variable, or a probe file that cannot be opened, stops every rank with a message: the variable was left
  undefined, and the run wrote to a unit never opened.
- A degenerate cell (a zero determinant of the metric) and a rejected interface state stop the run with a message that locates
  them: the corner nodes of the cell; the block, the face, the two states of the interface and its stencil. The degenerate cell
  was announced with a fallback that never ran.
- The diagnostic output of a steady run names and fills the soot residuals (`rho*Yp`, `rho*Np`) between ρE and the turbulence
  variables: with a soot model it named fewer variables than it wrote, and with a turbulence model too it stopped on a run-time
  error in builds with bounds checking.
- Where two faces of type `410` of a block meet, the ghost cell of their common edge is filled one face direction at a time:
  both faces wrote it in the same parallel loop and read it back for the viscous stencil, a data race whose result could depend
  on the threads. Chimera faces (`102`), which had the same race, write no ghost cell beyond the edges of a block face.
- The exact Riemann solver (`riemann-solver = exact`) takes for the face state the wave of each side, shock or rarefaction, that
  the last Newton iteration computed. It chose the wave with the velocity of the contact updated after that iteration, and where
  the two fell on the two sides of the normal velocity of the cell it read the speed of the other wave, which that iteration had
  not computed: the value of an earlier iteration, or an undefined one, in a DEBUG build a signalling NaN that stopped the run
  with a floating-point exception (the coarse Sod tube of `RiemannSmoke`).
- Time-dependent boundary data (the total pressure of `402`, the assigned state of `410`) are evaluated at the time of each
  Runge–Kutta stage, t^n + θ_k Δt with θ = (0), (0, 1) and (0, 1, 1/2) for forward Euler, Heun and SSP-RK3. They were
  evaluated at t^(n+1) in every stage, which made such runs first-order accurate in time (in a test, observed order 1.05 with
  Heun and 1.09 with RK3; now 2.02 and 3.29). Runs without `402` or `410` do not change.
- A time table of `402` or `410` before its first instant gives its first value. The value was undefined, so a steady run (at
  time −1) with `402`, or a time-accurate run that starts before its table, ended with NaN and wrote no solution; such a run
  now gives the field of a `401` with the same p0.
- The Spalart–Allmaras equation for ρν̃ has the term −(1/σ)(ν + ν̃)∇ρ·∇ν̃ of its standard conservation form (Allmaras,
  Johnson and Spalart, ICCFD7-1902, 2012, eq. 9), with which it equals the original equation, written without the density.
  Without the term the right-hand side was off by its size (0.154 and 0.160 of the scale on 32×32 and 64×64 cells of a
  constructed field, not converging); now it converges with second order. Results change where density and ν̃ vary together,
  as in boundary layers at high Mach number and at shocks; fields of uniform density do not change.
- In steady RANS runs the point-implicit factor 1/(1 + Δt·d) multiplies the whole increment of the turbulence equations,
  (F + S)/(1 + Δt·d), the point-implicit form of Yee and Shinn (NASA TM-89415, 1986, eq. 3.1), so that the converged solution
  solves F + S = 0 whatever the local time step. It divided the source alone, so that the steady state, wall shear and heat flux
  included, depended on the local time step. Steady RANS solutions computed with the default `point-implicit = .true.` change
  accordingly. The residuals of ν̃, k and ω printed during a steady run are scaled by the same factor.
- In 3-D the rotation tensor of the turbulence models, W_ij = (∂u_i/∂x_j − ∂u_j/∂x_i)/2, has W_23 and W_32 with their sign.
  They had the opposite one, so the matrix was not a tensor, and the models that take its components, not only its magnitude,
  depended on how the flow is oriented in space wherever ω_x ≠ 0: the term r̃ of the rotation–curvature correction (SA-RC and
  SST-RC) and the vortex-stretching function f_β of Wilcox 2006 (in a test, the source of SA-RC changed by 1.01 times its
  largest value when the flow was turned). 3-D runs with SA-RC, SST-RC or Wilcox 2006 change; SA, SST and 2-D runs do not.
- In a rotating frame the rotation–curvature correction (SA-RC and SST-RC) takes the rotation of the frame into r̃, the term
  (ε_imn S_jn + ε_jmn S_in)Ω′_m of its definition. It was left out, so a uniform shear in a frame turning about the normal to its
  plane had r̃ = 0 (f_r1 = 1.5 instead of 4.65 in a test). Runs with that correction in a rotating frame change wherever the
  strain is not zero; other runs, and runs in a fixed frame, do not.
- Where strain and vorticity are both below 1e-30 1/s, in a fluid at rest in a fixed frame or in a uniform flow, the
  rotation–curvature correction does not act: its factor 2r*/(1 + r*), r* = S/Ω, computed as 2S/(S + Ω), takes 1, the value of
  simple shear, and r̃ is 0. With the lower bound of Ω the factor was close to 0 there, so next to a wall at rest SA-RC turned
  the production of ν̃ into its negative; SST-RC multiplies μ_t S², which vanishes at rest, so it did not change there. Where
  only the vorticity vanishes the factor takes its limit 2, as before; elsewhere f_r1 changes by round-off only.
- In 3-D the vortex-stretching function f_β of Wilcox 2006 takes the whole tensor Ŝ = S − (1/2)(∇·u)I: its six off-diagonal
  components were never assigned, so χ_ω, f_β and the destruction of ω took values that the program does not define. 3-D runs
  with Wilcox 2006 change wherever vorticity and off-diagonal strain are not zero; 2-D runs do not.
- With `omega-wall-bc = asymptotic` the wall value of ω of Wilcox 2006 is 6ν/(β₀d²), with the β₀ = 0.0708 of its own ω
  equation. It took 6/β₁ = 80, the coefficient of SST, so the wall value was 5.9 % low, as if the wall were moved by 2.9 % of
  the first cell, an error that vanishes as the first cell shrinks. SST, and the default `omega-wall-bc = practical` (Menter's
  condition, 800ν/d², for both models), do not change. Among the cases of MOSE the turbulent flat plate with Wilcox 2006 takes
  `asymptotic`, and its results change.
- The QCR2000 stress takes the laminar and the eddy viscosity in the order in which the solver passes them: the quadratic
  correction acted on the laminar stress, and the eddy-viscosity stress went in uncorrected. `turbulence` does not accept
  QCR2000 today, as the guide now says, so no run changes.
