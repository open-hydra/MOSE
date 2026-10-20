# Changelog

All notable changes to MOSE are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- The key `diffusive-metric` of `[MOSE-Numerics]` chooses the face metric of the viscous fluxes: `stencil`, the default, from the
  same differences of the cell centres that the face gradient takes of the field (see Changed); `inverse-mean`, the inverse of the
  mean of the inverse metrics of the two cells, whose face-normal edge is the vector between the cell centres; or `mean`, the mean
  of the two metrics, as before.
- Probes can record `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde`
  (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model).
- Builds with `USE_MPI=ON` register a `needs-mpi` test tier (`ctest -L needs-mpi`). Each of its tests writes its own small
  multi-block case, runs the solver through the MPI launcher that CMake finds (`MPIEXEC_EXECUTABLE`, at most four ranks) and
  compares output files byte by byte. The CI runs this tier in a second job, `test-mpi`, which builds MOSE with gfortran,
  the Open MPI of Ubuntu and `USE_MPI=ON`.
- MPI runs take the rank of every block from `INPUT/decomposition.map` (with the phase prefix of `bc.txt` in multi-phase runs)
  when the case has one: the decomposition map that the ATLAS mesh decomposer (MDB) writes next to the split mesh, one line per
  block with its cells along i, j and k, its number of cells and its rank, checked against the mesh. A map that does not fit the
  mesh or the run (a block line missing or one too many, blocks out of order, other cells, a rank the run does not have, a map
  made for another number of ranks, a value that is not an integer, a directory in place of the file) stops the run at set-up:
  every rank prints the same message, finalizes MPI and ends with exit code 1. The map is read in every build, so a map made
  for more ranks also stops a serial run. With fewer blocks than ranks some ranks idle, with a warning, as without a map.
  Without the file the assignment of the blocks and the outputs are those of before, bit for bit; with it the solution does not
  change, only the rank on which each block is computed.
- At set-up MOSE warns, with the line number, about every `key = value` of `input.ini` that the reader does not take in its
  section: such a key is not read, and its parameter keeps its default. The guide says that a key lost this way is an error in
  the input. Inputs whose keys are all read do not change.
- A restart goes on exactly as the run in one go would, from a solution in any format: every field solution stores the
  counters that decide what the run does and writes next (the global iteration count `CYCLE`; `MOSE_ITERATION`,
  `MOSE_OUTPUT_NUMBER`, `MOSE_GRID_LEVEL`, `MOSE_LAST_OUTPUT_TIME` and `MOSE_RES_REFERENCE`), as field data of the `.vts` files
  of a vtk solution and as data set auxiliary data of a Tecplot solution, and a restart of the same kind of run reads them back.
  A restart from such a solution keeps the CFL ramp, the numbering of the solutions, the multigrid level and the convergence
  test of the uninterrupted run, and drops the rows of `residual-history.dat`, of the probe files and of `torque.dat` written
  after the solution it reads: the run split in two parts writes the same files as the run in one go, byte for byte, within the
  limits that the guide lists. A restart from vtk raw solutions did not start at all, and one from Tecplot started its counters
  again: the numbered solutions of the second part replaced those of the first. A solution written before restarts as before.
- The manufactured-solution driver (`MOSE_mms`) takes its field from the environment variable `MOSE_MMS_FIELD`: `periodic2d`,
  the field of `test/2D/viscous/mms` and the default, `periodic3d`, and the wall fields `wall2d` and `wall3d`, with walls at
  y = 0 and y = 1 and the wall temperature varying along them; an unknown value stops the driver at start. Where the total
  boundary flux of an equation is zero for any state (every equation of the periodic fields, the mass with walls), the driver
  removes the mean of the discrete source once after set-up, so that the discrete steady problem has a solution; the mean and
  its ratio to the mean of |S| go to the log.

### Changed

- The interior faces of the viscous fluxes no longer take by default the mean of the two metrics, which divided the difference
  across a face by the mean of 1/h instead of the distance between the cell centres: it overstated the face-normal gradient, and
  with it the viscous, conductive and diffusive fluxes, by (1 + r)²/(4r) at a cell-size ratio r (4.2 % at r = 1.5, 8.9 % at
  r = 1.8), at any resolution. Viscous runs change wherever neighbouring cells differ in size; `diffusive-metric = mean` keeps the
  previous metric. The default is `stencil` (below).
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
- When `input.ini` is read again during a run (every `ini-diter` iterations, 10 000 by default), a value outside the rule of its
  key stops the run with the message of the rule, as it does at set-up; it was taken without a word (a negative `cfl` written
  during a steady run made the residual grow until the run stopped on a negative pressure nine iterations later). The keys that
  the file leaves out keep their values: the set-up no longer writes the names of its summary into the values of
  `integration-variables`, `riemann-solver` and `flux-limiter`, and the summary prints the same names as before.
- A key of `input.ini` whose value does not read as the type of its parameter stops the run, at set-up and when the file is read
  again during the run, with exit code 1 and a message that gives the file, the section, the key and the value as written: a
  number that is not one (a word, `NaN`, `Infinity`, a decimal comma such as `0,5`), an integer written with a period or an
  exponent (`2.0`, `1e3`: gfortran refused it, ifx cut the fraction, `2.7` read as 2), a logical that does not start with `T` or
  `F`, a list with more values than its parameter holds, a key with no value; so does `levels` of `[MOSE-Multigrid]`. Such a key
  was taken as one not given: its parameter kept the default, or the value of the previous reading, in silence.
- MOSE builds with FiNeR v2.3.0 instead of v2.0.6. FiNeR v2.1.0 and later skip a comment line also when it holds `=` or a
  bracket, and report the failures of their parser: with v2.0.6 such a line ended the reading of its section, and the keys
  after it kept their defaults, so such an `input.ini` now gives the values written after the line. The five libraries that
  FiNeR v2.0.7 and later no longer carry (PENF, FACE, FLAP, BeFoR64 and StringiFor) are submodules of MOSE under
  `lib/third_party/`, at their release tags (v2.0.1, v1.1.17, v2.5.3, v1.2.4 and v1.8.0), and CMake adds them before FiNeR
  from `FINER_DEPS_PATH` (`lib/third_party/` by default): a build needs neither FoBiS nor a download besides the submodules.
  After pulling this change into an existing clone, run `git submodule update --init --recursive`: CMake initializes a
  missing submodule, but does not update one that is checked out at an older commit, and stops with this command when
  FiNeR is older than v2.1.0. ifx 2023.0 stops with an internal compiler error on StringiFor v1.8.0 and is not supported:
  this version is tested with gfortran 12 and 13 and with ifx 2024.1 (see the installation guide).
- `res-threshold` is relative to the density residual of the first iteration of a steady run (of each grid level): the run ends
  when its residual has fallen by that factor, and the default, `1e-10`, is a fall of ten orders. The test compared the residual
  as it is, which scales with the level of density and with the square root of the number of cells: the same flow stopped at
  different iterations in different units, and a flow at low density could stop after its first iteration, unconverged. Every
  steady run that ends on `res-threshold` ends at another iteration (the PrandtlMeyer case of the test suite runs 761 iterations
  instead of 555); runs that end on `iter-threshold` do not change. `residual-history.dat` and the screen still report the
  residual as it is. A restart from a solution that holds no reference measures the threshold from its own first iteration, with a
  warning. An input written for the absolute test keeps its stop with `res-threshold` divided by the density residual of the first
  iteration (see the user guide, Stopping criteria).
- With `sol-format = vtk raw` the field solutions hold their data as raw binary appended to each file, as the wall and diagnostic
  outputs did: they were base64 inside the XML, in a time-accurate run, and in a steady one up to its first diagnostic output.
  The values are the same, bit for bit; base64 took a third more room.
- The `.vtm` files of the field and wall solutions list their block files relative to their own folder, `vtk/<name>BlockN.vts`,
  as VTK readers resolve them; they listed `../vtk//<name>BlockN.vts`, outside `OUTPUT/`. An initial condition `INPUT/ic.vtm`
  reads its block files where the `.vtm` lists them, so a solution of MOSE copied into `INPUT/` with its `vtk/` folder starts a
  run; one whose block files are at `INPUT/vtk/field<entry>.vts`, where they were read before, is still read there.
- The outputs of a steady run carry minus the global iteration count as their time in every format: the vtk `TIME` of the field,
  wall and diagnostic outputs was −1, and the Tecplot solutions of the coarse grid levels of a sequence carried 0. The Tecplot
  zones of a steady solution are static (no `STRANDID`), and those of a time-accurate one are in the strand of their block
  (`STRANDID` = *b*). Tecplot solutions written by earlier versions of MOSE hold the count of a steady run as a positive
  `SOLUTIONTIME`, with `STRANDID = 0`: a steady restart goes on from them as before, but a time-accurate restart takes them for a
  solution at time t = count; to start a transient from such a solution, use it as the initial condition (`INPUT/ic.tec`,
  `newrun = true`), which starts at t = 0. The text formats, `vtk ascii` and `tecplot ascii`, write 17 significant digits, the
  fewest that give back every 64-bit value.
- A Tecplot field solution ends with six `DATASETAUXDATA` records, the counters of the restart: a record of the ASCII format that
  the Tecplot Data Format Guide defines for data of the whole data set. The Python reader of ORION reads such files; a script
  that takes every number after a zone header must stop at them.
- The density residual of `residual-history.dat` and the torque of `torque.dat` are summed in a fixed order, the same for any
  number of OpenMP threads and MPI ranks, so their rows and the convergence test that uses them no longer change with the
  decomposition. At one thread and one rank the order changes, and the sums can change in their last bits.
- MOSE uses ORION main a978930 (PR #7 merged; was 8cda645), whose multi-block VTK reader and writer handle field data (`TIME`,
  `CYCLE` and named scalars), any number of blocks and surface blocks, and whose Tecplot reader and writer handle the solution
  time with its sign, strands and data set auxiliary data. The vtk wall files now hold the cell values of their faces.
- The viscous fluxes take by default (`diffusive-metric = stencil`) a face metric made of the same differences of the cell centres
  that the face gradient takes of the field, so that the gradient of a linear field is exact on any mesh at the interior faces and
  at the block connections, and at periodic pairs, symmetry faces and the k faces of a wedge or of a 3-D sector (type `200`) for
  the linear fields that they allow. At the faces of type `200` the differences along the face are those of the row of the cell
  and of the row of the image, the row of the source turned by the angle of the wedge or sector (positions and values), with one
  metric per group of quantities where a ghost cell of a wall enters that row; with `inverse-mean` and `mean` these faces take,
  as before, the inverse-mean metric of the cell and its image, whose derivatives along the face are wrong by (q - 1)^2/(4q) where
  the cell size grows by q along the face and by (q - 1)/4 on the first row next to an axis (8.4e-3 and 5.0e-2 of the hoop term of
  a linear field at q = 1.2, test `MMSResiduals`). `inverse-mean` builds the metric from the mean edges of the two cells, which
  are not those differences: on meshes with nodes moved at random, cells of alternating size, grid lines that fold at a face
  between blocks, points shared by three or five blocks or skewed cells on a symmetry plane, a linear velocity or temperature,
  whose viscous force and heat conduction vanish, gave up to 5.4e-2 of their face-flux scale in the patch test with it, at any
  resolution. The metric is computed at set-up and stored on the faces (27 values per cell); a face between blocks takes one
  gradient, as an interior face. Results change for viscous runs on meshes that are not Cartesian; on meshes stretched across a
  boundary layer, where the tangential derivatives of `mean` and `inverse-mean` carry the factor (1 + r)²/(4r); and at the faces
  between blocks with SST, Wilcox 2006 or QCR, whose fluxes are not linear in the gradient, also on Cartesian meshes. A face whose
  stencil does not span three directions (a degenerate cell, a face that collapses) stops the run at set-up with its block, cell
  and face, and the remedy `diffusive-metric = inverse-mean`.
- A viscous run on a 2-D axisymmetric wedge or on a 3-D sector closed by type `200` whose axis nodes are exactly at r = 0 stops
  at set-up, with every `diffusive-metric`, with a message that asks to put the axis nodes at a small radius (the meshes of MOSE
  put them at 1e-6 of the radial size of the first cell). On such an axis the faces have zero area and their ghost cells copy
  the cells. With `inverse-mean` and `mean` the faces of type `200` take the derivative in r across the axis with a ghost cell one
  cell away, and the hoop stress of the first row of cells is wrong, by 1/6 of the hoop term in Homann's stagnation-point flow on
  the wedge and on the sector, and refinement does not reduce it. With `stencil` the faces take the difference on one side of the
  copy: Homann's flow gets on the first row the hoop stress of a small axis radius (2.55e-5 against 2.58e-5 of the hoop term on 64
  cells), but such an axis is not treated completely (ghost cells beyond it, not copies), and the run stops with every metric. A
  wall on a face of zero area, as on the axis of a 3-D mesh
  around it, stops with its remedy: the symmetry condition (`300`) on the axis. With `inverse-mean` and `mean` a viscous run on a
  3-D mesh stops at set-up at any boundary face of zero area, as on such a polar line, with a message that asks for
  `diffusive-metric = stencil`: with them the viscous force of the first row of cells of a linear field, which is zero, grows as
  the mesh is refined, while `stencil` is exact there.
- In time-accurate RANS runs the local time step takes, for the stiffness of the turbulence sources, a bound that the state of
  the cell gives alone: Spalart–Allmaras 2 c_w1 (1 + c_w3⁶)^(1/6) ν̃/y², SST max(β*, 2 max(β₁, β₂)) ω, Wilcox 2006
  max(β*, 2 β₀) ω. Each step evaluated the turbulence sources once more, with a fill and an exchange of the ghost cells, only
  to take their derivative. The step of a cell where the stiffness sets it is shorter, by the ratio of the bound to the
  derivative (1.10 for SST where F₁ = 1, at most 1.18 for Wilcox 2006 in 3-D, about 2 for Spalart–Allmaras in an equilibrium
  boundary layer); the stages and the positivity of ν̃, k and ω do not change; steady runs do not change.
- The viscous fluxes compute the viscous variables of each cell once per stage, the face metric of the faces of type `200`
  once at set-up, the cosine and sine of the periodicity angle once, and read the stored face metrics in place: a viscous run
  holds one more array of the size of the state. The same operations act on the same values; with the `-fp-model fast=2` of
  the ifx release build the results can differ in the last bits.
- The wall output (`OUTPUT/wall.tec`) writes, on every face of every wall, the traction of the gas on the wall, the heat flux
  into the gas and the mass flux into the gas: it wrote them along the face normal of the solver, toward growing index, so two
  walls that mirror each other wrote opposite values by the parity of their faces. On faces 2, 4 and 6 `tauX`, `tauY`, `tauZ`
  and the `mdot` of the walls `503`–`505` change sign, on faces 1, 3 and 5 `qw`; `y+`, `pw`, `Tw` and the solution do not
  change. A script that undid the parity must drop that.
- The configuration refuses ifx 2024.1 at `-O2` or above without interprocedural optimization (a RELEASE build with
  `USE_IPO=OFF`, or whose toolchain cannot do IPO, and the TESTING build), with a message: that compiler builds a loop of the
  ghost metrics wrongly there (a minimal program in standard Fortran gives wrong values at `-O2` and `-O3` and the right ones at
  `-O0`, `-O1` and with `-ipo`), and MOSE stopped at set-up on the check of the ghost metrics, or ran with wrong ghost metrics.
  The RELEASE build with IPO, the default, and the DEBUG build configure as before, with the same compile and link lines.

### Removed

- The compressibility correction of Spalart–Allmaras (`turbulence = SAcomp`): it is refused at set-up. It took c from
  c² f(c) = 25 Ω ν_t / a, a number against a speed, so the factor it gave the production depended on the units: 0.99 with
  lengths in metres and 0.22 with lengths in millimetres for Ω = 1e3 1/s, ν_t = 1e-3 m²/s and a = 340 m/s; at Ω = 1e5 1/s it
  cut the production by 78 % whatever the Mach number. No case or test of the repository used it.
- Implicit residual smoothing (`irs = true` or `irs-beta > 0`) in time-accurate runs: it is refused at set-up; steady runs keep
  it. The smoothed residual is not the residual of the space discretization, so a run with it did not follow the transient of
  the discretized equations: in one step from the same state the time derivative changed, by an amount that did not depend on
  the step, and so did the final density of a Sod run. The change does not vanish as the step does. MOSE has no dual time
  stepping, in whose pseudo time residual smoothing belongs.

### Fixed

- Faces between blocks take the inverse of the mean of the inverse metrics with `diffusive-metric = inverse-mean` or `mean`, and
  chimera faces with every value of the key: across a large size jump at a block connection the mean of the two metrics made the
  normal gradient far too steep, and the explicit diffusion unstable on the large cell.
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
- With `k-energy-coupling = .true.` (SST, Wilcox 2006) the energy of the mean flow takes, in steady runs, the point-implicit
  increment of the energy and k equations together: besides the sink of k with the factor of the k row, the share 1 − f_k of
  the fluxes of k, f_k = 1/(1 + Δt·β*ω). It took the sink alone: with the whole increment of k scaled by the factor (above)
  the converged energy, wall heat flux included, depended on the local time step (6.6e-6 of the wall heat flux at half the
  step on the plate of `PointImplicitSteady` with SST, 5.6e-6 with Wilcox 2006), and E_m + k was not conserved at each step.
  Now the two equations add up to the fluxes of E_m + k at every stage, and the converged solution does not depend on the
  step. Steady runs with the coupling change; time-accurate runs and runs without the coupling do not.
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
- The volume of a cell is the integral of det J of the trilinear cell, whose faces are the bilinear faces of their four nodes,
  those of the area vectors of the fluxes: V = 1/6 sum_d e_d · (S_d⁻ + S_d⁺), e_d the mean edge of direction d and S_d⁻, S_d⁺
  the area vectors of its two faces of that direction. It was a sum of five tetrahedra in absolute value, which cut each face
  along one diagonal while the cell next door cuts it along the other: on a twisted face the two cells missed or counted twice
  T · S_f/12 (T = c1 − c2 + c3 − c4, the twist of the face), so the volume of a cell was not the volume its faces enclose
  (6.7 % and 18.5 % too large on cubes with one twisted face, T_n/h = 1 and 5) and the volumes of a mesh did not sum to its
  volume (4.0e-3 short on a twisted mesh of 8³ cells). On cells with plane faces (2-D and axisymmetric meshes, 3-D meshes of
  plane faces) the two are equal in exact arithmetic: results change in their last bits. A cell whose volume has not the sign
  of its block, a cell turned over away from its centre, stops the run at set-up with its block and cell.
- The wall distance of the turbulence models is the exact minimum distance from the cell centre to the wall surface, a wall face
  being the bilinear surface through its four nodes, the face of the cells, and on a 2-D axisymmetric wedge the distance in the
  meridional plane. It was the distance to the nearest wall-face centre, too large wherever the grid lines are not normal to the
  wall and next to convex corners and the ends of walls (up to 10.9 %, 16.3 % and 15.9 % in three test geometries). A wall face
  plane to the resolution of its coordinates, like every wall face of a 2-D or axisymmetric mesh, is the four triangles between
  its centre and its sides. A twisted 3-D wall face, whose corners are not coplanar, lies up to |T|/16 from those triangles
  (T = c1 − c2 + c3 − c4, its twist; with them the distance was off by up to 4.9 % in a test): its distance is a certified
  minimum over the face, exact to the rounding of the size of the face. Results change wherever the distance changes,
  and at round-off where the grid lines are normal to a plane wall. Equal distances go to the lowest face index, so the
  roughness a cell inherits is the same on every rank and thread. With a turbulence model, a wall on faces 5–6 of a wedge
  stops the run.
- The type of the mesh (plane, axisymmetric wedge, 3-D sector closed by type `200`, slab translated in k) is taken from every
  node of every block, to within 1 % of the move of the nodes per cell in k, no longer from two probe nodes of block 1. A plane
  2-D mesh whose j index runs towards the x axis was taken for a wedge (wall distance off by 27 % in a test), and so was a plane
  slab far from the x axis; the periodic images of a 3-D slab periodic in k were turned, a mesh with a wedge block and a plane
  block ran as one of the two, and a 3-D mesh with a thin block 1 was taken for a 2-D one. A mesh whose blocks disagree, a 2-D
  block whose k planes fit neither a turn nor a translation to that tolerance (too few digits, or k planes that are not
  parallel), and a 3-D block closed by type `200` whose last k node plane is neither the first one turned nor translated, now
  stop at set-up; the set-up summary names the type. A 3-D block closed on itself over a full turn about the x axis, its
  last k node plane on the first one (to the round-off), is translated by zero, the periodicity by the identity, as before
  the check. The move of the nodes per cell in k is measured along the k lines, the sum of the moves of each node across
  its cells in k: on a sector past half a turn the tolerance is 1 % of its arc per cell in the norm of the largest
  component (between 0.71 and 1 of the arc), not of its chord per cell (4.2 times larger at 270 degrees). In 2-D a block
  whose probe node keeps its azimuth while its k = 1 node plane is not the k = 0 plane translated to 1 % (k planes that
  touch or cross, or a k = 1 plane scaled radially about the axis), and which the translation fits no better
  than the turn by the null angle, stops at set-up as a block whose k planes are not parallel: it passed for plane (and,
  with the planes crossed, stopped at its folded cells). Meshes recognised right keep their wedge angle to the bit.
- The wall distance is the minimum over the images of the wall across periodic boundaries too (periodic connections `201`, slabs
  closed in k by type `200`, rotational periodicity of 3-D sectors). It was measured to the wall inside the mesh only, so near
  the end of a period it could be several times too large (by 9.8 times the exact distance in a test channel). A mesh without
  periodic boundaries costs nothing more, and the images searched are bounded whatever the size of the mesh against its
  periods: every wall face and every cell centre is moved into one cell of the lattice of the periodic translations, reduced
  pairwise and three at a time. With a turbulence model, periodic translations that cannot be the images of one wall stop
  the run with a message: not integer combinations of each other, more than three independent ones, a number of images too
  large for an integer, or not along x on a wedge or with sectors. Translations that differ by at most 1 % of the thinnest
  cell across the periodic boundaries are the same, as the two sides of a periodic boundary match only to the round-off of
  the coordinates: a channel periodic in x, turned by 30 degrees and written with 4 significant digits, runs. For a slab
  closed by type `200` the cell is the move of its nodes per cell in k along its k lines, as for the mesh type: a block
  closed on itself over a full turn about x gives the translation zero, or of the round-off, which is no periodicity.
- On a 3-D mesh closed in k by type `200` (a sector turned about the x axis, or a slab whose last k node plane is the first
  one translated), the ghost cell beyond face 5 holds the state of the cell next to face 6 turned back by the sector angle
  (copied, on a translated slab): its w was the opposite of the turned one (face 6 was right). The second-order
  reconstruction of the faces next to face 5 reads that ghost cell, and so do the viscous fluxes of the cells next to it;
  a radial flow that is its own mirror image departed from it by 1.8e-2 of its scale in 20 steps. Results change on 3-D
  sectors wherever the velocity next to face 6, turned by the sector angle, has a z component, and on translated slabs
  wherever w is not zero.
- Type `200` joins the two k faces of one block: beyond face 5 lies the cell next to face 6 of the same block, in the same
  column (i, j), and beyond face 6 the cell next to face 5, as the record names no other block. A record of type `200` on a
  face 1 to 4, with its cell outside its block, or with no record of type `200` on the other k face of its block in its
  column now stops the run at set-up, on every rank, with a message that names the record and the remedy. A 3-D sector
  cut along k into blocks (in a test the rows of cells next to the cut and to the k faces were wrong by up to 31 times the
  scale of a step), a ring of 360 degrees in several blocks closed by type `200` between the last block and the first,
  type `200` on a face 1 to 4 of a 2-D wedge, and a wedge with type `200` on one k face only, whose two k faces took two
  different conditions, ran in silence with wrong fluxes; on a 3-D mesh type `200` on a face 1 to 4 stopped at the first
  step, on the rank of its block alone, and a record outside its block read outside the arrays. A sector is cut along i
  or j instead, and the blocks of a ring are joined by `101`, the junction at 2π included. Other meshes do not change.
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
- The inlets at total conditions (`401`, `402`) and the ambient boundary (`407`) let the gas in when the gas of the boundary
  cell moves in at 2a/(γ−1) or faster (R ≤ 0, normal Mach 5 or more at γ = 1.4): no pressure brings that gas to rest, so
  p_stop = 0 and every reservoir with p0 > 0 feeds the face, as the criterion of the guide says. A `401` or `402` face was a
  slip wall and a `407` face had no flux, in silence: in a 1-D tube with the gas entering at 6a from a reservoir at p0 = p and
  T0 = T the face now lets in 0.579 ρa instead of nothing. The face takes the inflow state on the invariant R, which links the
  cell to the face whatever its sign, when that state exists and is subsonic across, and the choked state, normal Mach 1,
  otherwise; no inflow state exists when the gas of the cell, even expanded to vacuum, leaves faster than the gas of the
  reservoir can follow. A reservoir much hotter than the cell can give a subsonic face (T0 = 25 T, p0 = 0.1 p and the gas at
  5.5a: normal Mach 0.83 and 0.0113 ρa), where a face choked whenever R ≤ 0 would let in 2.8 % more and jump by 8 % where R
  changes sign; the mass flux is continuous there. Results change only at faces where R ≤ 0.
- With `sol-format = vtk` the wall output sets the block sizes of its files, which it never set: depending on what they held, a
  run with a viscous wall could stop in the wall output or write wall files of no extent and one point. Each vtk wall file now
  has the extent and the points of its face. Tecplot wall files do not change.
- With more than one OpenMP thread per MPI rank, MPI is called only from the thread that initialised it, as
  `MPI_THREAD_FUNNELED`, the thread level MOSE requests, requires. Inside the parallel region of the Runge–Kutta stages the halo
  exchange and the exchange of the smoothed residuals were wrapped in `single` constructs, which any thread of the team may
  execute: with more than one thread per rank, MPI calls could come from threads other than the primary one. They are now
  `master` constructs, which the primary thread executes, each followed by the barrier that the `single` construct implied.
  Results do not change.
- MOSE builds and runs with gfortran 10 or later and an MPI library whose Fortran interface gives no explicit interface for the
  routines with choice buffers, such as Intel MPI 2021.12 with gfortran 13. The build stopped with errors in `Mod_MPI.f90` (the
  same routine called with buffers of different types and ranks), and a build forced through them crashed in
  `MPI_Reduce`: MOSE compiles with `-fno-underscoring`, under which the special addresses that the `mpi` module of the library
  declares in COMMON blocks (`MPI_IN_PLACE` among them) are not those of the library. With GNU 10 or later and `USE_MPI`,
  `Mod_MPI.f90` is compiled with `-fallow-argument-mismatch`, which turns those errors into warnings, and the three reductions
  that used `MPI_IN_PLACE` take a separate send buffer, on every compiler; the values they return are the same.
- With `chemistry = equilibrium` an MPI run on more than one rank no longer passes to the equilibrium the blocks that the rank
  does not hold: the explicit step computed the equilibrium on every block of the mesh, also on the blocks of the other ranks,
  whose state is not allocated on the rank after the partition. It now skips them, as the finite-rate chemistry already did, and
  on several ranks such a run gives the bits of one rank. Runs on one rank, and runs without equilibrium chemistry, do not
  change.
- In an index direction with one cell (the k direction of a 2-D case, the j direction of a 1-D case) with symmetry (`300`) on its
  faces, the metrics of the second layer of ghost cells no longer depend on the order of the OpenMP threads: the set-up filled
  every ghost layer in one loop, and the second layer of one face read the first layer of the other while another thread wrote
  it. Each layer of the symmetry faces now has a loop of its own, after the previous one. The results of the tests do not change.
- At start-up one thread writes the number of OpenMP threads, in MOSE and in the manufactured-solution driver (`MOSE_mms`);
  every thread wrote it, a data race without effect on the value.
- A release build with ifx and interprocedural optimisation finds the `llvm-ar` and `llvm-ranlib` of the compiler in the
  layouts of oneAPI 2023 and 2024 (`<bin>/compiler` from 2024 on, `<bin>/../bin-llvm` up to 2023, `<bin>`), also for an `ifx`
  reached through a symbolic link; `-DIPO_AR` and `-DIPO_RANLIB` give them for any other layout. With oneAPI up to 2023 IPO
  was switched off, and the fallback hid the release flags set after it, so the targets were compiled without
  `-ffast-math`, and without `-fomit-frame-pointer` until a second configure; with oneAPI 2024 every configure also added one
  more `-ipo` to the link flags.
- `space-reconstruction = MUSCL` without `flux-limiter` runs van Leer, as its warning says, and `first-order` ignores a
  `flux-limiter` given, with a warning; the set-up summary names the limiter in use. MUSCL without a limiter ran first order
  (zero slopes) while the summary said MUSCL, first-order with a limiter ran MUSCL with that limiter, and first-order without a
  limiter printed a spurious warning. Results change only for these two inputs.
- `shock-detector = Tramel` or `Chen`, the values that its rule admits, turns that shock sensor on. The set-up compared the key
  with `tramel` and `chen`, so a value of the rule never turned a sensor on and, written, it also kept off the sensor that
  HLLC+Tramel, HLLE++, HLLC+Chen and AUSM+M take without the key: the reconstruction stayed second order at shocks, and
  HLLC+Tramel was HLLC everywhere. Runs that leave the key out do not change; runs that write it change near shocks.
- A run with a chemistry model and `exclude-blocks` other than `none` no longer stops at set-up: the list of the blocks without
  chemistry was deallocated before it was ever allocated, an error that ends the run (exit code 2 with gfortran).
- `exclude-blocks` excludes its blocks from the equilibrium chemistry too, and on every multigrid level. The equilibrium step
  changed the composition and the pressure of every block, and with multigrid the excluded blocks got the chemistry on the
  coarse levels, with either model, and carried it to the fine one. Runs without `exclude-blocks` do not change.
- A word of `exclude-blocks` that is not a block number, or a number beyond the blocks of the mesh, stops the run at set-up with
  a message. A word that did not read emptied the list and the next word was written beyond it: a build with bounds checking
  stopped with an index error, and a release build went on with chemistry in every block; a number of four or more digits kept
  its first three, and a number beyond the mesh named no block, in silence.
- With `chemistry = equilibrium` and no element of the species read (no `INPUT/composition.txt`), the run stops at set-up with a
  message that names the file, before the initial condition is read. The equilibrium keeps the amount of each element, and
  without the elements the solver had no such constraint: a build with bounds checks stopped inside it. The guide of the
  equilibrium says which file it needs, and its format.
- The divergence check of the state update stops the run on any value that is not finite (NaN or infinity), or on a negative
  pressure, with any compiler and options: it tested `isnan`, which gfortran folds to false with `-ffast-math`, an option of
  the RELEASE build (a diverging run of such a build ended with a segmentation fault), and which misses an infinity. Its
  message is "NaN, infinity or p<0 detected". The user guide gives the exit code of a run that MOSE stops on an error, and the
  known limit that errors raised inside ORION end with exit code 0.
- A value that is not finite (NaN or infinity) in the initial condition, or in the solution that a restart reads, stops the run
  at set-up with exit code 1, in every build and on every rank, with the place of the value: the variable of a cell (by its
  position after the coordinates, with its name when the file gives the names) or the coordinate of a node, the block and the
  indices. A build that traps the invalid floating-point operations (the DEBUG build) ended with SIGFPE, exit code 136, at the
  first operation on the value. In the other builds a NaN velocity stopped the run at the divergence check of the state update
  ("Integration failed"); an infinite coordinate ended it on a bound of an array (gfortran with bounds checking, exit code 2)
  or on a degenerate cell ("Metric tensor det=0", ifx). Runs whose initial condition is finite do not change.
- An inlet `402` names its total-pressure table as written: a file name with a directory, quoted, or longer than 32 characters
  (up to 256). The property line was read as a list, so a slash ended the read and the values after the name were never read,
  and a long name was cut. A table that cannot be opened, holds no (time, value) pair or has decreasing times stops the run at
  set-up with the name of the file: a missing table was not reported, and a table with decreasing times was read without a word,
  the run ending with exit code 0. So does a time or a value of the table that is not finite (NaN, infinity), with the line and
  the column: every comparison with a NaN is false, so a NaN pressure closed the face of the inlet in silence (no gas through
  it) and the run ended with exit code 0, a NaN time passed the check of decreasing times, and a build that traps the invalid
  floating-point operations (the DEBUG build) stopped with SIGFPE.
- A property line of `bc.txt` with an empty value (two commas with nothing between them) for a value without a default stops
  the run at set-up, with the line of the file, the record, its type and the name of the value; the value was left undefined,
  and the run went on. The values with a default keep it, and the guide lists them.
- A line of `bc.txt` that is short of a value or holds a value that is not a number stops the run at set-up, with the line, the
  record and its type: read as a list from the file, a short line took its last values from the next line, and every record
  after it was read out of place. Blank lines between records are skipped, as before.
- A flow direction (`alpha`, `beta`) that is neither a number nor the token `normal` stops the run at set-up: a misspelt token
  (`Normal`, `nromal`) was taken as 0 rad without a word.
- A value of a property line of `bc.txt` that is not finite (NaN or infinity), required or with a default value, the flow
  directions included, stops the run at set-up with exit code 1, in every build, with the line of the file, the record, its type
  and the name of the value. The test of an empty value let a NaN through, as every comparison with a NaN is false: a NaN total
  pressure of an inlet `401` closed its face in silence and the run ended with exit code 0 (with `-ffast-math`, an option of the
  ifx RELEASE build, it was taken for an empty value), an infinite total temperature of an inlet `402` ran to the end with exit
  code 0, a NaN roughness of a wall `302` did so in a gnu TESTING build and stopped at the divergence check in the ifx RELEASE
  build, an infinite total pressure or a NaN direction stopped at the divergence check ("Integration failed"), and a build that
  traps the invalid floating-point operations (the DEBUG build) stopped with SIGFPE at the first operation on the value. The
  messages of the reader name the file it reads, which on a coarse multigrid level is the file of that level, not `bc.txt`.
- A time (`SOLUTIONTIME`) or a value of the file of an assigned state (`410`) that is not finite (NaN or infinity) stops the
  run at set-up with exit code 1, in every build, with the file, the zone, the line of its header, and the variable and the
  cell of the value. A NaN pressure ended the run at its first step: with a run-time error of an index out of bounds (gfortran
  with bounds checking), with a non-physical state at a face that named no file (ifx), or with SIGFPE (the DEBUG build); an
  infinite time was taken without a word, and the run ended with exit code 0.
- A file of an assigned state (`410`) that cannot be read stops the run at set-up with exit code 1, in every build, with the file
  and, for a zone, the zone and the line of its header: a file missing or with no `VARIABLES` line or no zone, a zone header
  whose `I`, `J`, `SOLUTIONTIME` or block number of the title does not read, values that are not numbers or that end before the
  values that `I` and `J` give, one value too many. So does a file whose zones do not come in increasing time, with the two
  times: the reader took the times in the order of the file, and the face took the state of the wrong time in silence. The reader
  did not check its reads: the values not read kept those of the zone before, or the zones after were read shifted, and the run
  went on; its stops were ERROR STOP.
- A zone header of the file of an assigned state (`410`) that gives no `I=` or no `J=` stops the run at set-up with exit code
  1, in every build, with the file, the line of the header, the zone and the form to write (`I=<n>, J=<n>`, in upper case with
  no blank around `=`). The reader takes the nodes of a zone only so written: a header with `I = 3` or `j=2` kept the default
  of 2 nodes in silence, the values of the zone were read for that size, and the run ended with exit code 0. Other forms are
  still not read.
- With the rotating frame, `stationary-face-<n>` entries no longer make the start-up hang: the count of the entries asked for
  `stationary-face-1` at every pass, so stationary walls could not be used. A `rotational-frame` other than `rigid-body`, its
  default `none` included, turns the frame off before its entries are read; with `none` the switch of the frame was left
  undefined.
- The user guide gives the data lines of the boundary types `420`, `502` and `506` as MOSE reads them: the fifth value of `420`
  is the mass flux `g` through the face, not a throat radius; the line of `502` ends with the mass fractions of the products;
  `506` was missing.
- The guide names the boundary files of the multigrid levels as MOSE reads them, `bc.txt`, `bc2.txt`, `bc3.txt`, … (it said
  `bc_x1.txt`, `bc_x2.txt`, …).
- Every rule of the parameter registry is checked, at set-up and when `input.ini` is read again: intervals such as the `in (0,
  1)` of `preconditioning-eps-min`, and the rule of a list of numbers on each of its values (`index-position` of a probe, `>=
  0`; `ode-relative-tol-species` and `ode-absolute-tol-species`, `> 0`). Only the bounds of single numbers and the lists of
  names were checked, and a value outside the other rules ran with no message. A rule of another form stops the run as an error
  of the registry. Runs whose values keep to the rules do not change.
- A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity) stops the run at set-up, in every
  build, with a message that names the key (the run went on to the divergence check, or stopped with the message of an axis of
  zero length, or with SIGFPE in a DEBUG build). So does a real parameter of `input.ini` read as an infinity (gfortran reads
  `1e999` as one), with a rule or without, at the start and when `input.ini` is read again: a bound such as `> 0` let it through.
- An `input.ini` that cannot be read stops the run with the reason, at set-up and when it is read again during the run, on every
  MPI rank and with exit code 1: a file that cannot be opened, that holds no line that is not blank, or that the INI reader
  rejects. A missing, unreadable or comments-only file stopped the run on the first required key without naming the file, an
  empty one ended the run inside the reader, and a file that could not be read again during the run left the values of the
  previous reading in place, without a word.
- With MPI, an error that every rank finds stops every rank together, and the message of every rank reaches the log: an error in
  a time table, in a line of `bc.txt`, in a value of `input.ini` read again, in the mesh type or in the periodic images of a
  wall. Every rank called `MPI_ABORT`, which makes only a best attempt to stop the processes: the first call could stop the
  others before they wrote their line, and at times no line reached the log.
  An error that rank 0 alone finds at set-up (the data tables, the rules of `input.ini`, the boundary conditions) stops every
  rank together too, once rank 0 has written it: with `MPI_ABORT` from rank 0 alone, a run with Intel MPI 2021.12 could end with
  exit code 141 (SIGPIPE) and without the last line of its log. An initial condition that one rank cannot read stops every rank
  as well: the ranks decide together on the error of the read, rank 0 writes the line of the set-up summary with the verdict of
  all, and the first rank that found the error writes it, with its rank. Rank 0 decided alone, so a run in which rank 1 could not
  read its copy of the file went on to its end; and with the lean read of a `.szplt` file, where each rank reads the variables of
  its own blocks after that decision, an error of that read was never looked at.
- With MPI, a rank that stops every rank with `MPI_ABORT` (an error that it finds alone, such as a folded or degenerate cell,
  the axis at r = 0, a wall face of zero area or a probe of an unknown variable) writes out its output and its error output
  first. Built with gfortran and run with Open MPI 4.1.4 as one process with its output on a file, as CTest runs the fast
  tests in an MPI build, the run ended without writing them out: the message, the output that the process still held before
  it and the warnings on the error output were lost, and the log of a stop at set-up held only the lines of Open MPI. With
  Intel MPI 2021.12 nothing was lost. Builds without MPI do not change.
- With MPI, rank 0 reads `input.ini`, at set-up and whenever it is read again during the run, and passes its text to every rank,
  which parses it: every rank takes the same values. Every rank read the file on its own, so a file rewritten while the ranks
  read it again could give them different values, and a different `sol-diter`, `res-diter` or `ini-diter` made the ranks call
  their collective operations in different orders: the run hung until the time limit of the job (two ranks in two folders, the
  value changed in the `input.ini` of rank 1 alone). A file that rank 0 cannot read stops every rank once rank 0 has written
  why; the `input.ini` in the folders of the other ranks is not read.
- The VTK outputs are well-formed XML: the variable names were written with the double quotes of the Tecplot lists
  (`Name=""rho(1)""`), not well-formed XML (Python's xml.etree, which parses with expat, stops on them), and a reader
  that takes the text between the first two quotes
  got empty names. The diagnostic vtk files hold their own columns only; they also held arrays with empty names and stale
  content of the field. The values written do not change.
- A restart from a vtk solution reads it: the `.vtm` listed its blocks outside `OUTPUT/`, the restart read the blocks of the
  initial condition and in base64 whatever the files held, and the time of the solution was dropped, so a time-accurate restart
  went on from t = 0. With multigrid the coarse levels of a restart took their time and iteration count from themselves, not
  from the solution read.
- A restart drops the rows that the run it continues wrote after the solution it reads, and a last row cut in the middle: the
  rows of those iterations were in `residual-history.dat` and in the probe files twice, and a row cut by hand went on with the
  first row of the restart. The files end with a newline before the restart appends to them.
- A restart appends to `torque.dat` of the rotating frame, after the rows of the run it continues: it replaced the file, which
  then held only the rows of the second part.
- Every row of `residual-history.dat` and of the probe files is in its file as soon as MOSE writes it, and the run closes these
  files and `torque.dat` at its end, before MPI is finalized: built with gfortran, the rows waited in the buffers of the files
  until the process ended normally, so a run killed or stopped at the time limit of a queue lost its last rows, and a restart
  from its last solution wrote files with a gap; an MPI run whose ranks the launcher ended in the shutdown of MPI (OpenMPI
  4.1.4 reports such a run as "exiting improperly") could leave probe files empty.
- A time-accurate restart from a steady Tecplot solution starts its transient at t = 0, as one from vtk does: it started at
  t = the iteration count of the solution. A steady restart from Tecplot goes on from the iteration count as before.
- A restart from a solution that holds no time (a vtk file without `TIME`, a Tecplot ASCII file without `SOLUTIONTIME`) stops at
  set-up with a message that says how to start from it, as an initial condition: it went on from t = 0, or from an iteration
  count of −huge.
- A restart from a solution whose time (Tecplot `SOLUTIONTIME`, vtk `TIME`) or one of the counters it holds is not finite (NaN or
  infinity) stops the run at set-up with exit code 1, in every build, with the file and the name of the value: the tests of the
  time took a NaN as false (another message, SIGFPE in a DEBUG build, or with ifx a run from a NaN time), and an infinite counter
  was taken in silence.
- An initial condition given as `vtk raw` or `vtk ascii` is read in the format its block files declare, and from the folder of
  its `.vtm`: a solution of MOSE copied into `INPUT/` stopped the run at set-up.
- Time-based probe rows fall on the multiples of `dtime`, from the first: the time counter of a probe had no initial value, so a
  new run could write time-based rows from its first step, or none, and a probe with `diter` only wrote time-based rows.
- A probe given by position records the cell whose centre is the nearest to it: it took the cell that has the nearest node as
  its upper corner, a ghost cell when that node lies on the first plane of nodes in j or k (as on the mid-plane of a 2-D case
  or a wedge), otherwise a neighbour of the cell that holds the probe unless that node is its upper corner. The rows of a probe
  change where the two cells differ; probes given by index do not change.
- The wall output leaves the residual and the ghost cells of the wall cells as they were: the diagnostic output, written after
  it, showed in the wall cells the wall shear force and heat flow instead of the residual. The solution does not change.
- When `input.ini` is read again during a run (every `ini-diter` iterations), a new `omega`, `axis` or `origin` of the rotating
  frame is taken by every term of the frame: the centrifugal source and the torque file took the new `omega`, while the Coriolis
  source, the velocity of the stationary walls, of the inlets and of the far field, and the absolute vorticity of the turbulence
  models kept the old one, and a new axis or origin did not enter at all. A new `rotational-frame`, or a
  `stationary-face-<k>` key that names another face or that the set-up does not have, stops the run with a message, on
  every MPI rank together: they were ignored without a word. The `stationary-face-<k>` keys that the file read again
  leaves out keep their faces, as every key left out keeps its value.
- A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity), written in `input.ini` during the
  run, stops the run when the file is read again, every rank at that reading, with a message that names the key; in a gnu TESTING
  build it went into every term of the frame and the run failed at the next iteration; a DEBUG build stopped with SIGFPE at that
  reading, an ifx RELEASE build with the message of an axis of zero length.
- A run that ends where a numbered or time-based solution or a probe row is due writes it, with its wall and diagnostic files,
  before the final solution, whose counters follow it: the last iteration wrote the final solution only, so the run in one go
  lacked its last numbered solution and probe rows, and a restart from a first part that ended there numbered the solutions
  after it one lower, shifted the time-based ones and left the state of that iteration in no file. The guide no longer lists
  this as a limit of the restart.
- With `diffusive-metric = stencil` the gradient at a wall face takes, along the wall, the differences of the wall values (the
  temperatures of the neighbouring `302` faces; zero for the velocity of a no-slip wall and for the turbulence variables) with
  the same differences of the positions, so that it is exact for a linear field whose wall values are its values at the face
  centres. It set the derivatives along the wall to zero and mapped the result with the metric of the wall cell, whose direction
  across the face is the wall normal only on parallelogram cells: where the wall temperature varies along a skewed wall the
  normal gradient was off by (dT_w/ds)·tan θ, at any resolution, and on wall cells that are not parallelograms the direction of
  the velocity gradient was wrong too (1.4e-2 to 8.8e-2 of the scale in the patch test; now at most 3.5e-13). Results change only
  there.
- With `diffusive-metric = stencil` the source terms of the turbulence models (Spalart–Allmaras and its variants, SST,
  Wilcox 2006, SSG/LRR) and the eddy viscosity of the output take cell gradients whose metric comes from the same differences of
  the positions as of the values, with each ghost value at the point where it is right. With the metric of the cell the gradient
  of a linear field was wrong next to a jump of cell size (5 % in the wall cell of a layer growing by 1.2), inside a stretched
  layer, behind walls on skewed cells (a copied density took a normal gradient it does not have) and next to zero-gradient faces
  (half the one-sided difference). Results change for RANS runs, mostly in the cells next to boundaries.
- With `diffusive-metric = stencil` a mesh with faces of zero area, such as an axis at r = 0 or the polar line of a 3-D mesh
  around its axis with symmetry on it, no longer stops at set-up on a not-a-number: such a face has a zero normal, as in the
  rest of the solver, its ghost cell is the cell itself, and its flux is zero.
- The ghost cells behind viscous walls (`301`, `302`, `503`–`505`) and extrapolation faces (`400`) are filled once per stage, for the fluxes and the source terms alike, from the quantities whose differences the fluxes take:
  with a wall value q_w the ghost takes 2q_w − q₁ (the velocity of the walls and the turbulence variables), the temperature of a
  `302` wall max(2T_w − T₁, T_w/2), and the other quantities an extrapolation kept within [q₁/2, 2q₁]; the density follows from
  p/(R T). The extrapolation 3q₁ − 3q₂ + q₃ of every primitive variable could turn negative and amplify (a negative ghost
  density next to a hot wall), put the wall value on the face only when the profile was linear over three cells, and the source
  terms read a second set of ghost values. A line at the end of a run counts the ghost fills in which the cut or the floor
  acted. Results change for every run with a wall or a `400` face, in the first rows of cells next to them.
- The ghost cells behind the boundary faces other than `101`, `201` and `300` take their nodes by a linear extrapolation along
  the index lines, so that both layers keep the edge across the face of the first cell. The cubic extrapolation made the second
  ghost layer vanish where the second cell is 1.5 times as wide as the first, and the first layer at twice: a run stopped on a
  determinant of exactly zero, or went on with a metric as large as 1/det, and at chimera faces (`102`) the viscous flux took an
  infinite, undefined or reversed value where the second cell is three times the first or more. A cell whose determinant is
  below 1e-12 of the product of its edges, or a ghost cell of the first layer behind a chimera face turned over or degenerate
  against the first cell, now stops the run at set-up with its block, cell, face and ratio, and the message gives the
  thresholds and the remedy (index lines that change smoothly across the first two node layers at the face). Behind the
  other faces no flux reads the metric of the ghost cells, and only a ghost cell of zero size stops: a mesh whose edges along
  such a face grow 3 times or more across the first cell runs. Every run changes in the last bits next to such faces, and by
  O(r − 1) where the cells are stretched across them.
- A cell folded at one of its corners, or turned against the orientation of its block, stops the run at set-up, naming block,
  cell, corner node, the Jacobian there and the bound; the coarse levels of a sequence of grids are checked too. Such a cell
  passed in silence: the volume of a cell was a sum of tetrahedra taken in absolute value, and the cell check looked at the
  determinant of the mean edges only. On a valid mesh nothing changes, bit for bit.
- The two k faces of a 3-D sector closed by type `200`, one interior face of the whole annulus, carry its viscous flux (viscous
  stress, heat flux, diffusion of species and turbulence variables), computed once between the cell next to face 6 and the
  turned image of the cell next to face 5. They carried the inviscid flux only, so the solution tended to that of a sector
  whose periodic planes are free of traction and adiabatic, an error whose share grows as the mesh is refined. Runs on the wedge do not change.
- On a 3-D block closed in k by type `200` with no angle (a translated sector), the ghost cells sit one period away, where their
  values belong: with `diffusive-metric = stencil` the metrics put them at the cell at the other end of the block, which scaled
  the derivative along k of the two layers of cells next to the k faces by 2/(2 − nk) (−1/3 with 8 cells), at any resolution.
- A wall `301` applies its heat flux `q_wall` into the fluid on every face, as the guide defines it. On faces 1, 3 and 5 it
  applied it out of the fluid, so a positive `q_wall` cooled the gas there, and the wall temperature it estimated followed the
  same wrong sign. Results change for every wall `301` with a non-zero `q_wall` on faces 1, 3 or 5.
- With the rotating frame, the wall output of a stationary wall uses the velocity of that wall, −(Ω × r) at the face centre, as
  the boundary fluxes do: the skin friction and the heat flux written for it were those of a wall rotating with the frame.
- On a wall `301` or `302` of a stationary face of the rotating frame (`stationary-face-<n>`), the wall gradient takes the
  derivatives of the wall velocity along the wall, from its gradient −[Ω×], with every `diffusive-metric`. They were zero, so
  the wall stress lacked μ times them, an error of order μ|Ω| on any mesh, orthogonal cells included: a gas at rest in the
  laboratory felt a shear μΩ on such a wall. The wall stress changes there, in the wall flux and in the wall output.
- On the same walls the energy flux of the wall face carries the work of the wall stress on the wall velocity, (τ·n)·w, as the
  faces between cells carry it with the velocity of their cells. The energy of the wall cells lacked it, by an amount per unit
  area of wall that does not depend on the mesh.
- A `stationary-face-<n>` whose face covers a record that would take the velocity of the wall at rest in part only (the
  inviscid flux of a symmetry `300`, of the inlets and outlets `401`–`406`, `408` and `420`, and of the walls `103`, `104` and
  `502`–`505`), whose face has no wall `301` or `302`, or that has no record in `bc.txt`, stops the run at set-up with a message
  that names the key, and for a record the record and its type: on a wall `503`–`505` such a face was at rest in the inviscid
  flux and turning with the frame in the viscous flux, a boundary condition of no wall, and a face with no record did nothing.
  The walls of a stationary face may share it with the types that never read whether their face is stationary, whose boundary
  condition is the same with or without the key (block connections, periodic and chimera faces, extrapolation `400`, `407`,
  `410`, `501`, `506` and unused faces), as a casing face that is a wall in part and a block connection in part.
- With the rotation–curvature correction of Spalart–Allmaras (`turbulence = SA-RC`) the ghost cells behind the faces between
  blocks (`101`) and of the periodic pairs (`201`) hold the velocity gradient of the cell across, from another rank by message,
  and on a 3-D sector closed by `200` the gradient of the turned image, so the result no longer depends on how the mesh is
  divided into blocks. They held a copy of the gradient of the cell, so next to such a face the strain-rate gradient took half
  of a one-sided difference. With SA-RC a mesh of one block changes next to its periodic pairs, next to the
  faces of a 3-D sector turned about x and, with swirl, on the 2-D wedge; without SA-RC nothing changes.
- With a RANS model, the ghost cells behind an SRM grain (`502`) take the turbulence values that its flux uses in the same
  stage: those of the blowing wall while the surface burns (ν̃ = 0, k = 0 and the wall value of ω), those of the impermeable
  wall at the flame temperature when it does not, as 2φ_w − φ₁ times the density of the ghost cell, in both ghost layers; and
  the wall distance counts the grain as a wall. Next to the grain the distance was that to the other walls, or 10⁸ with no other
  wall, so the destruction and the production of Spalart–Allmaras, the blending functions of SST and the wall value of ω at the
  grain were those of a flow without the grain there. Without RANS, or without a `502` boundary, nothing changes.
- With `integration-variables = prec` and SST or Wilcox 2006, the viscous floor of the reference velocity of each cell takes the
  eddy viscosity at the strain rate of the cell. The velocity gradient was computed only by SA-RC, so for the other models it
  stayed zero, and SST and Wilcox 2006, which limit their eddy viscosity with the strain rate, gave the eddy viscosity of zero
  strain, larger wherever the limiter acts: there the reference velocity was too large, and with it the local time step, the
  preconditioned update and the flux HLLC-PC (in a test the local step was 0.771 of the rule with SST and 0.782 with
  Wilcox 2006). Runs without `prec` (the default is `cons`) do not change.
- The ghost cells of the walls coupled to another solver (`103`, `104`) take the ratio of the extrapolation and, with
  `diffusive-metric = stencil`, the positions of the ghost values of a wall: they kept a ratio of 0, so the extrapolated
  quantities were copied, and the face metrics placed their ghost values where they are not, so the viscous flux of the faces
  between two wall cells was wrong for a linear field on skewed cells. Their velocity and turbulence variables take the wall
  values 2q_w − q₁ of an adiabatic wall, as the other walls: the turbulence variables were copied from the cell.
