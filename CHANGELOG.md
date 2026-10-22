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
- Builds with `USE_MPI=ON` register a `needs-mpi` test tier (`ctest -L needs-mpi`). Each of its tests writes its own small
  multi-block case, runs the solver through the MPI launcher that CMake finds (`MPIEXEC_EXECUTABLE`, at most four ranks) and
  compares output files byte by byte. The CI build has no MPI and does not run this tier.
- At start-up the root process reports how far the wall distance can be from exact on twisted 3-D wall faces: the largest
  bound, min(|T|/16, |T_n|/2), over the distance of the cell next to the face, with that cell (T is the twist of the corners
  of the face, T_n its part along the normal of the face). It is a warning above 1 %. 2-D and axisymmetric meshes have plane
  wall faces and print nothing. The distances do not change.

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
- A face between blocks (`101`, `201`) whose block across has fewer than two cells across that face, at any grid level of a
  multigrid sequence, stops the run at set-up with a message that names the record, its block and face, the block across and
  the level. The second ghost cell of such a face is a ghost cell of the block across, which another record fills in the same
  loop: a release build read past its node array and ran on with no message, and a build with bounds checking stopped with an
  index error.
- The two records of a face between blocks (`101`, `201`) must point at each other: a record that points at a cell and face
  with no connection record, at its own cell and face, or at a record that does not point back at it stops the run at set-up
  with a message. Such records were accepted, and each side computed the face with its own ghost cells.
- `rel_fac` no longer acts on the inlets at total conditions (`401`, `402`) and on the ambient boundary (`407`); it stays in
  their data lines. When the equation of the face state had no real root it moved p0 towards the pressure of the cell, and
  with `rel_fac = 1` it repeated the same computation 1000 times and then set the discriminant to zero. Now that the gas enters
  only when the reservoir can push it in (see Fixed), that case is not reached, and p0 is used as given. Results change at
  most at round-off.

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
- The fluxes of the boundary faces no longer depend on the order of the records in `bc.txt`, and each MPI rank walks only the
  records of its own blocks. The boundary fluxes were computed in one worksharing loop per range of records, in file order, on
  every rank: the face contributions of a corner cell were summed in the order of the file, and with the records out of order
  and more than one thread two threads could update the same corner cell inside one loop, a data race whose result is
  unspecified; every rank also read every boundary record of the mesh at every stage. The records of the local blocks are now
  grouped by block and face, keeping the file order within a face, with one loop per face of a block. A `bc.txt` sorted by
  block and face, as the guide lists the records, gives the same solution bit for bit; any other order gives the solution of
  the sorted file.
- A face between two blocks (block cut `101`, periodic connection `201`) is computed once, and its flux is given with opposite
  signs to the cells on its two sides. Each block computed the face for its own cell with its own inputs, the shock sensor of
  its own cell and, with HLLC-PC, a reference velocity of zero in the ghost cell. The two shares then differed: mass, momentum
  and energy were not conserved at the face (with HLLC+Tramel the total mass of a periodic shock tube in two blocks drifted by
  4.3e-3 of its value in 40 steps), and a mesh split into blocks did not give the solution of the unsplit mesh (up to 1.8e-1
  apart). The record that computes the flux is that of the cell on the low-index side of the face; between faces of the same
  parity (index directions permuted across the cut), the record with the smaller (block, face, k, j, i), so that the result
  does not depend on the order of the records. Between MPI ranks the share of the other cell travels in one more message per
  stage. With `preconditioning-Ur-smooth` > 0 the smoothing next to such a face now sees the cell across it, where it saw a
  zero. Results change where a shock sensor or the preconditioning acts at a face between blocks, and at round-off where the
  two blocks store the nodes of the common face with different digits or a periodic connection translates them.
- Every MPI exchange uses one tag, within 32767, the bound that the MPI standard guarantees on every library. The tags were the
  rank of the sender plus an offset per exchange, up to 30000 plus the rank for the residual smoothing, and the output gathers
  used the block number: beyond 2768 ranks, or 31767 blocks, a library that offers only the least bound may reject such tags as
  invalid, and an MPI error stops the run. Each receive names its source, so the rank in the tag told nothing more. Results do
  not change.
- In a planar case with one cell in k and symmetry (`300`) on all its k faces, round-off in w no longer grows until the run
  stops. The two mirror planes one cell apart give w the mode dw/dt = −2aw/Δz (a the speed of sound), which leaves the real
  stability interval of Heun when 2aΔt/Δz exceeds 2, and that of RK3 when it exceeds 2.51, as on a layer thin against the
  in-plane cells. Between such faces, when they are planes normal to z, the z momentum is no longer integrated and w is set
  to zero at set-up, with a warning. A layer whose k faces are not normal to z (planes inclined to z) keeps w, as the gas
  moves along them, across the layer: it is integrated as before, with a warning, and on a thin layer the mode can
  still stop the run. Null (`0`) and periodic (`200`) k faces and axisymmetric wedges are unchanged.
- On the 2-D axisymmetric wedge (one cell in k, type `200` on the k faces) the two k faces carry their viscous stress. Only
  their inviscid flux was computed, so the radial momentum lacked −τθθ/r, with τθθ = 2μu_r/r + λ∇·u, which near the axis is as
  large as the radial stress terms it balances, and with swirl the faces also lacked τrθ. On Homann's stagnation-point flow the
  radial viscous force was off by 14.45 and 28.67 times the largest exact force on 32×32 and 64×64 cells; it now converges with
  second order away from the axis. Results change for every viscous run on a wedge; planar meshes with type `200` on k are
  unchanged.
- On the 2-D axisymmetric wedge the local time step includes the stiffness that the stress of the wedge faces adds next to the
  axis, for either sign of the wedge angle. That stress damps the radial momentum of the first rows of cells off the axis,
  about (16/3)ν/Δr² in the first row of square cells, beyond the viscous limit of the step: at the vnn the rule declares
  stable (0.62 with RK3, 1/2 with Heun) a diffusion-dominated 1° wedge broke down in the first row. The rule adds s_h, the
  increase of the largest row sum of the viscous operator that the stress brings; nothing changes off the wedge. The local
  step is smaller next to the axis: in the wedge cases of the tests by up to 20.4 % and 25.5 % of the largest step, on 32×32
  and 64×64 cells.
- The source terms of the turbulence models in the cells next to a wall of type 301, 302 or 503–505 see the ghost values set from
  the wall values, ghost = 2 wall − cell (CFL3D User's Manual 5.0, NASA/TM-1998-208444, eq. H-186). The ghost update that
  precedes the source terms copied the cell over them, so with Spalart–Allmaras the c_b2 term of a wall cell kept a quarter of
  its value (∇ν̃ = a/2 instead of a for ν̃ = a y). In the fast tests the Spalart–Allmaras and SST fields next to a wall change.
  The walls coupled to another solver (103, 104), which have no such ghost values, keep the copy of the cell.
- The wall output takes the records of each face of each block whatever the order of `bc.txt`. With records not sorted by
  block and face, the wall faces took the records of other faces and the wall output wrote 0 for their wall values (320 of
  the 856 values of `wall.tec` in a four-block case, laminar and Spalart–Allmaras). A sorted file gives the same wall output
  as before.
- An ambient boundary (`407`) whose flow direction is given as `normal` lets the gas in along the face normal. It took the
  sentinel value of the token for an angle (tan(huge) = −4.96e-3), so the gas entered with a spurious transverse component
  (1.21e-3 of the normal velocity in a test duct) and the closure of the inflow used a wrong coefficient. Results change
  wherever a `407` record gives `normal`; records with angles are unchanged.
- A 1-D case (one cell in j and in k) with symmetry (`300`) on all its j faces, or on all its k faces, no longer fails on
  round-off in v or w: the two mirror planes one cell apart give v the mode dv/dt = −2av/Δy, outside the stability region of
  the explicit schemes at the step of the x cells (three such cases stopped with `NaN or p<0` after 14 and 15 steps). Such a
  velocity is now set to zero at set-up and not integrated, with a warning that names the faces, when those faces are planes
  normal to its axis (a tube of constant section); in a channel of variable section it is integrated as before, with a
  warning, and thin cells can still stop the run. Results change only in these 1-D cases.
- The wall distance of the turbulence models is the exact minimum distance from the cell centre to the wall surface, a wall face
  being the four triangles between its centre and its sides, and on a 2-D axisymmetric wedge the distance in the meridional
  plane. It was the distance to the nearest wall-face centre, too large wherever the grid lines are not normal to the wall and
  next to convex corners and the ends of walls (up to 10.9 %, 16.3 % and 15.9 % in three test geometries). Results change
  wherever the distance changes, and at round-off where the grid lines are normal to a plane wall. Equal distances go to the
  lowest face index, so the roughness a cell inherits is the same on every rank and thread. With a turbulence model, a wall on
  faces 5–6 of a wedge stops the run.
- The type of the mesh (plane, axisymmetric wedge, 3-D sector closed by type `200`, slab translated in k) is taken from every
  node of every block, to within 1 % of the move of the nodes per cell in k, no longer from two probe nodes of block 1. A plane
  2-D mesh whose j index runs towards the x axis was taken for a wedge (wall distance off by 27 % in a test), and so was a plane
  slab far from the x axis; the periodic images of a 3-D slab periodic in k were turned, a mesh with a wedge block and a plane
  block ran as one of the two, and a 3-D mesh with a thin block 1 was taken for a 2-D one. A mesh whose blocks disagree, a 2-D
  block whose k planes fit neither a turn nor a translation to that tolerance (too few digits, or k planes that are not
  parallel), and a 3-D block closed by type `200` whose last k node plane is neither the first one turned nor translated, now
  stop at set-up; the set-up summary names the type. Meshes recognised right keep their wedge angle to the bit.
- The wall distance is the minimum over the images of the wall across periodic boundaries too (periodic connections `201`, slabs
  closed in k by type `200`, rotational periodicity of 3-D sectors). It was measured to the wall inside the mesh only, so near
  the end of a period it could be several times too large (by 9.8 times the exact distance in a test channel). A mesh without
  periodic boundaries costs nothing more, and the images searched are bounded whatever the size of the mesh against its
  periods: every wall face and every cell centre is moved into one cell of the lattice of the periodic translations, reduced
  pairwise and three at a time. With a turbulence model, periodic translations that cannot be the images of one wall stop
  the run with a message: not integer combinations of each other, more than three independent ones, a number of images too
  large for an integer, or not along x on a wedge or with sectors. Translations that differ by at most 1 % of the thinnest
  cell across the periodic boundaries are the same, as the two sides of a periodic boundary match only to the round-off of
  the coordinates: a channel periodic in x, turned by 30 degrees and written with 4 significant digits, runs.
- On a 3-D mesh closed in k by type `200` (a sector turned about the x axis, or a slab whose last k node plane is the first
  one translated), the ghost cell beyond face 5 holds the state of the cell next to face 6 turned back by the sector angle
  (copied, on a translated slab): its w was the opposite of the turned one (face 6 was right). The second-order
  reconstruction of the faces next to face 5 reads that ghost cell, and so do the viscous fluxes of the cells next to it;
  a radial flow that is its own mirror image departed from it by 1.8e-2 of its scale in 20 steps. Results change on 3-D
  sectors wherever the velocity next to face 6, turned by the sector angle, has a z component, and on translated slabs
  wherever w is not zero.
- The inlets at total conditions (`401`, `402`) and the ambient boundary (`407`) let the gas in when the reservoir can push it
  in: when p0 is above p_stop = p (R/a)^(2γ/(γ−1)), the pressure that the gas of the boundary cell reaches when the outgoing
  wave brings it to rest (R the outgoing Riemann invariant). The choice followed the sign of the velocity of the cell: a tube at
  rest never started from a `401` or `402`, a `407` at rest took its outflow branch and computed an inflow outside its own
  domain (0.854 instead of 0.770 ρa at p0 = 2 p), and an inlet with p0 below p_stop let gas out. When the gas does not enter,
  a `401` or `402` face is a slip wall and a `407` face an outlet at p_back, as before for gas moving out; with the gas of the
  cell moving in, at min(p_back, p_stop), where the face velocity points out. Results change where the sign of the cell
  velocity and the criterion disagree.
- An inflow at total conditions whose face state would be supersonic across, built on an invariant that can no longer reach
  the face, takes the choked state, normal Mach 1. Past a pressure ratio (p0/p = 6.12 for a cell at rest at the reservoir
  temperature, normal inflow) the subsonic equations gave such a face: at p0 = 36 p normal Mach 1.85 and a mass flux of
  13.9 ρa, below the 20.8 ρa of the choked face. Results change for inlets at high pressure ratios, `401` and `402` included.
