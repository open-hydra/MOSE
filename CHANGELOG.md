# Changelog

All notable changes to MOSE are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- The key `turbulence` of `[MOSE-Physics]` selects the SSG/LRR-ω Reynolds-stress model of the NASA Turbulence Modeling
  Resource, in four variants: `SSGLRR` (TMR SSGLRR-RSM-w2012), `SSGLRR-SD` (simple diffusion), `SSGLRR-w2019` (with the
  length-scale correction of Eisfeld and Rumsey) and `SSGLRR-w2019-SD`. `Lib_SSGLRR` held the model, but the registry
  refused the name, so no run used it. The model exchanges energy with the mean flow as `k-energy-coupling` does; its
  sources take the point-implicit factor of steady runs (the energy row with them) and, in time-accurate runs, a bound of
  their stiffness in the time step; a rotating frame adds the Coriolis term to the stress equations; `omega-wall-bc`
  sets its wall value of omega, blowing walls have their own routine, and omega is kept positive. The boundaries carry the
  stress tensor: walls take the wall values in the wall stress, symmetry planes of any orientation reflect it, faces of
  type `200` turn it with the velocity, and inflow, outflow and extrapolation faces carry the Reynolds stress of the
  boundary cell in the momentum and energy fluxes. Tests: `RSMDecay` (validation tier, the four variants against the
  closed-form decay of isotropic turbulence on an axisymmetric wedge) and the flat plate of the NASA TMR
  (`FlatPlateTMR` and `FlatPlateTMRFull`, `test/2D/viscous/flat-plate-turbulent-tmr`): SSG/LRR on its three grids
  against CFL3D and TAU.
- The key `sst-dilatation-dissipation` of `[MOSE-Turbulence]`: a compressibility correction of the SST dissipation for
  compressible shear layers (Rumsey, NASA/TM-2009-215705, eqs. 40–43). The destruction coefficients become
  β*(1 + ξ* F(M_T)) and β − β* ξ* F(M_T), M_T = √(2k)/a the turbulent Mach number: `wilcox` (F = M_T² − M_T0² above
  M_T0 = 0.25, ξ* = 2; for jets and mixing layers) or `zeman` (F = 1 − exp(−((M_T − M_T0)/Λ)²) above M_T0 = 0.2,
  ξ* = 0.75, Λ = 0.66; for boundary layers). β − β* ξ* F is floored at 0 (Wilcox: M_T > 0.65), the point-implicit factors
  take the corrected coefficients. Default `none`: every run without the key is unchanged.
- The SRM grain (`502`) has two modes, the fifth field of its property line: `1`, a wall that blows the products in, as
  before; `2`, an inflow, applied as the mass-flux inlet `404` at the total temperature `T_af` along the face normal, with
  the mass flux of Saint-Robert's law at the pressure of the boundary cell and the turbulence values of the line: the
  turbulence values, ghost cells and wall distance of an inlet. Unlike an inlet `404` a grain of mode 2 is never closed by
  the gas of the boundary cell at rest or moving out: the grain burns whatever that gas does, as in a chamber at rest
  before the ignition.
- Extrapolation faces (`400`) report the gas that enters the domain through them, which they let in with the state of the
  boundary cell, as they take the flux of the cell state for either direction of the gas. The first time the gas enters
  through the `400` faces of one face of a block (on each grid level), the run prints a line with the iteration, the faces and
  the largest |w|/a among them (w the velocity of the gas of the cell along the normal out of the domain, a its speed of
  sound), and one more each time that largest value at a stage is more than 10 times that of the last line; at the end of a
  run in which the gas entered through them, the flux evaluations where it entered and their largest |w|/a. The flux of the
  faces does not change.
- The key `diffusive-metric` of `[MOSE-Numerics]` chooses the face metric of the viscous fluxes: `stencil`, the default, from the
  same differences of the cell centres that the face gradient takes of the field (see Changed); or `mean`, the mean of the two
  metrics, as before. Its first value, `inverse-mean`, is removed (see Removed).
- The key `irs-variables` of `[MOSE-Numerics]` chooses the residuals that implicit residual smoothing smooths: `flow`, the
  default, those of the species, momentum and energy, as before; or `all`, every residual, those of the turbulence, soot and
  passive-scalar equations too. With `flow` the turbulence equations advance with the local time step of the smoothed flow
  equations but without smoothing, and above their explicit limit they are the part that diverges: on a 3-D aerospike nozzle
  restarted from a converged field, IRS 0.5 at CFL 1.6 stopped after 884 iterations, the Spalart–Allmaras residual growing
  first, and ran stable with `all` (measured before the time step and face metric of this release). `all` is for
  Spalart–Allmaras: with the k-omega models it diverges (the SST and Wilcox 2006 flat plates of the tests, IRS 0.5: stop at
  iterations 28 and 95; with `flow` both run). Smoothing changes the path to the steady state, not the steady state.
- Multigrid runs write the wall quantities of a coarse level, `wall-level<n>`, with every field output of that level
  (`field-level<n>`) and at the end of the level, with the same solution time. The name was built, but the wall file was
  written on the finest level only: the friction and heat flux of a coarse level could be read only from `wall-prolongated`,
  the wall of its solution prolongated to the finest grid, not from the level itself. The other files do not change.
- The key `prolongated-diter` of `[MOSE-Multigrid]` spaces the writes of `field-prolongated` and `wall-prolongated` on the
  coarse levels: with a value above 0 they are written with a field output of the level only once that many iterations of the
  level have passed since the last one; the end of a level always writes them, so the file written at the end of level 2 is
  still the exact start of level 1. On a mesh of 23 million cells that file is 3 GB and its write took about 6.5 minutes,
  half the time of a coarse level at `sol-diter` 5000. The default, 0, writes them with every field output, as before; the
  key is read again at runtime.
- Probes can record `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde`
  (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model).
- Builds with `USE_MPI=ON` register a `needs-mpi` test tier (`ctest -L needs-mpi`). Each of its tests writes its own small
  multi-block case, runs the solver through the MPI launcher that CMake finds (`MPIEXEC_EXECUTABLE`, at most four ranks) and
  compares output files byte by byte. The CI runs this tier in a second job, `test-mpi`, which builds MOSE with gfortran,
  the Open MPI of Ubuntu and `USE_MPI=ON`.
- A `verification` test tier holds the full forms of the verification cases, whose short forms run in the `validation`
  tier: `MMSViscous3DFull`, the full form of `MMSViscous3D` with the grid with walls (18 h at one thread with a RELEASE
  build on one core of a shared cluster node; TIMEOUT three days). The PR CI leaves the label out, and
  `ctest -L verification` runs it on request: a plain `ctest` now runs it too, so the whole suite without it is
  `ctest -LE verification`.
- The turbulent flat plate of the NASA Turbulence Modeling Resource as a verification case, `test/2D/viscous/flat-plate-turbulent-tmr`,
  with its V&V page: Spalart–Allmaras, SST, Wilcox 2006 and SSG/LRR-ω on the 35 × 25, 69 × 49 and 137 × 97 grids of the
  TMR, with the inflow at the total conditions of the TMR and its free-stream turbulence. `FlatPlateTMR` (validation tier,
  about 150 s of CPU at one thread) restarts each model for 5000 iterations from its converged 35 × 25 solution: the skin
  friction at x = 0.970084071 and 1.90334 and the drag must stay within 1e-5 (relative) of the converged values, and the
  criteria of the page, computed again from the stored results of the full form, must hold. `FlatPlateTMRFull` (label
  `verification`, about 51 h at one thread, an estimate) runs the twelve cases to 300 000 iterations and checks the
  137 × 97 grid against CFL3D and the second code of the TMR on that grid (FUN3D; TAU for SSG/LRR-ω) and the limit of
  the three grids against CFL3D on the 545 × 385 grid.
- The key `decomposition-map` of `[MOSE-IO]` names a file of `INPUT/` that gives the MPI rank of every block (with the phase
  prefix of `bc.txt` in multi-phase runs), in the format of the decomposition map that the ATLAS mesh decomposer (MDB) writes
  next to the split mesh; with `none`, the default, MOSE assigns the blocks itself and reads no map. A map that does not fit the
  mesh or the run, a named file that is not there, or a value with a folder stops the run at set-up, on every rank, with exit
  code 1; a named map is read in every build, so a map made for more ranks also stops a serial run. With fewer blocks than ranks
  some ranks idle, with a warning, as without a map. The map changes only the rank of each block, not the solution.
- At set-up MOSE warns, with the line number, about every `key = value` of `input.ini` that the reader does not take in its
  section: such a key is not read, and its parameter keeps its default. The guide says that a key lost this way is an error in
  the input. Inputs whose keys are all read do not change.
- The key `ini-error` of `[MOSE-IO]` decides what an error of `input.ini` read again during the run does (a file that cannot be
  read, a value that does not read, is out of its rule or is not finite, a change of the model or of the stationary faces of
  the rotating frame, a value of its axis or origin that is not finite or does not read, a stationary face that does not read,
  an axis of zero length). With `warn`, the default, the reading is not taken: a warning gives the iteration and the error,
  which names the key, and the run goes on with every value it had before the reading, those of the keys of the same reading
  that were right too, on every MPI rank. With `stop` the run stops, on every rank, with the message and exit code 1. The
  value in force when the reading starts decides. At set-up every error stops the run, whatever the key says.
- A restart goes on exactly as the run in one go would, from a solution in any format: every field solution stores the counters
  of the run (the global iteration count `CYCLE`; `MOSE_ITERATION`, `MOSE_OUTPUT_NUMBER`, `MOSE_GRID_LEVEL`,
  `MOSE_LAST_OUTPUT_TIME` and `MOSE_RES_REFERENCE`), and a restart reads them back: the run split in two parts writes the same
  files as the run in one go, byte for byte, within the limits that the guide lists. A restart from vtk raw solutions did not
  start at all, and one from Tecplot started its counters again. A solution written before restarts as before.
- The manufactured-solution driver (`MOSE_mms`) takes its field from the environment variable `MOSE_MMS_FIELD`: `periodic2d`,
  the field of `test/2D/viscous/mms` and the default, `periodic3d`, and the wall fields `wall2d` and `wall3d`, with walls at
  y = 0 and y = 1 and the wall temperature varying along them; an unknown value stops the driver at start. Where the total
  boundary flux of an equation is zero for any state (every equation of the periodic fields, the mass with walls), the driver
  removes the mean of the discrete source once after set-up, so that the discrete steady problem has a solution; the mean and
  its ratio to the mean of |S| go to the log.
- The patch test of the viscous fluxes checks the cells next to inlets and outlets too (`ViscousPatch`, part
  `ViscousOpenPatch`, fast tier): the ghost cells of `401`–`408` and `420` must sit, in the metric of the face gradients, where
  their values belong. The solver does not change.
- A multirate (local time stepping) integrator for time-accurate runs, `[MOSE-Multirate] enabled = .true.`: every cell
  advances with a sub-step close to its own stability limit, the step H/2^l of the macro step H for a cell of level l,
  instead of the smallest step of the mesh. The scheme is the recursive Heun scheme of Jeanmasson, Mary and Mieussens
  (2019) with the interpolant of the coarse data of Hoang, Ju, Leng and Wang (2020) and a conservative correction of the
  fluxes at the faces between levels: mass, momentum and energy are conserved to round-off, and the macro step is of
  second order in time on a fixed mesh with a fixed level map (Appendix A of `docs/theory/time-integration.md`). The
  keys of `[MOSE-Multirate]` set the deepest level allowed (10 by default; the map takes fewer levels when the local
  steps span fewer octaves, so that the macro step stays within the step of the slowest cells divided by the safety
  factor), the safety factor, a target or a fixed macro step, the level map (rebuilt from the local steps, by default at
  every macro step, or kept) and the check of the non-local condition (4) of Gnedin, Semenov and Kravtsov (2018) at the
  start of every sub-step (`nonlocal-cfl = count` or `redo`): by default (`redo`) a macro step that violates it is redone
  with finer levels or a smaller macro step,
  and the coarse cell of a face between levels runs within 1/kappa of the limit of the finer cell, kappa = max(cfl (1 +
  theta), 1.1), so that, with kappa at most the safety factor (the default 2 with cfl up to 1), a rejection for a
  violation of condition (4) needs a fall of more than 1 - 1/kappa of the admissible
  step of a cell within the macro step (proven in Appendix A; the value 1.1 is a practice taken from Castro, the fall that its check of the step
  tolerates, not a value that the scheme gives; below cfl (1 + theta) = 1.1 no key turns it off),
  and a macro step in which a cell or a face takes a state that is not physical (a value that is not finite, a total
  density that is not positive, a negative pressure or a temperature outside the tables of the gas) is redone with half
  the macro step on the same levels when the map has more than one level; a face that is not physical already at the
  start of the macro step stops the run, as in the single-rate step. A macro step whose H falls below 2^-10 of its first
  attempt, or whose attempts are all rejected up to `max-attempts` (21 by default), stops the run with its rejections by
  cause, and at the start of every macro step a local step that is not a positive finite number (zero, negative,
  infinite or not a number) stops the run with the first such cell. The integrator runs the Euler and the Navier–Stokes
  equations, laminar or with the Spalart–Allmaras model (`SA`, `SA-R`) or the SST model (with `k-energy-coupling` and
  either `sst-production`), and with `max-level` 0 also `SA-RC` and `SST-RC`, with the SRM grain (`502`) in its
  mode 2, an inflow, on several blocks (the 2-D axisymmetric wedge included)
  and on several MPI processes, each with OpenMP threads, with the same output bits as on one process; the macro step
  lands on the output instants, and with `levels = recompute` a restart from a `vtk raw` or a `tecplot ascii` solution
  goes on bit for bit as the run in one go. The configurations outside its scope stop the run at set-up or at its first
  macro step with a message `[MOSE-Multirate] <item>` and a non-zero exit code; with `max-level` > 0 and
  `nonlocal-cfl = redo`, with `levels = recompute` and no `fixed-macro-step`, so do a `cfl` above 1 and, in a Navier–Stokes
  run, a `vnn` above 0.5, the linear stability limits of Heun on the checkerboard mode, with every limiter, up to which the
  check of condition (4) admits a sub-step within the macro step (with `off` and `count`, and where the certificate only
  counts, nothing bounds those sub-steps, and the set-up checks neither value). `safety` takes
  every value of at least 1 (2 by default). The keys of `[MOSE-Multirate]`
  keep the values of the set-up for the whole run: a reading of `input.ini` during the run that changes one, or that takes
  `cfl`, `vnn`, `time-accurate` or `cfl-rise-threshold` out of the scope of the integrator, is an error
  of the reading, which `ini-error` decides. In the `Timing` block and in the `MOSE Timing` and `MOSE Ranks` lines the
  exchange wait counts the wait for the messages of the exchanges of the integrator and their copy from the receive
  buffers, and the collective wait its reductions; the start of the exchanges (receives and sends posted, send buffers
  filled) is part of the compute and has its own line, `Exchange start`, in the `Timing` block of a run on several
  processes; the wall times per level are read with the clock of these timers, in
  builds with and without OpenMP, and a run on several processes ends with the minimum, mean and maximum over the
  processes, per level, of the wall time of the stages without their exchanges and of the time of the exchanges. With
  the multirate disabled, the default, the results do not change (the `Timing` block changes as the entry under Changed
  says).
- The tests of the multirate integrator (directory `test/multirate`, whose cases are built by
  `test/multirate/mr-common/gen_mr_cases.py`; `ctest -L multirate -LE verification` runs them without the longer forms,
  `ctest -L multirate` with them): a smoke test in the fast tier (`MultirateSmoke`); the validation groups
  `MultirateIdentity`, `MultirateConservation`, `MultirateOrder`, `MultirateOrderInletsMoving`, `MultirateShocks`,
  `MultirateInput`, `MultirateDaruTenaud` (the viscous shock tube of Daru and Tenaud) and `MultirateNozzleStartUp` (the
  start-up of a nozzle); in builds with `USE_MPI=ON`, `MultirateHybrid` in the `needs-mpi` tier; and, with the label
  `verification`, the full forms of the order checks (`MultirateOrderFull`) and the medium and full forms of the shock
  tube of Daru and Tenaud (`MultirateDaruTenaudMedium`, `MultirateDaruTenaudFull`). The V&V page `docs/vv/multirate.md`
  lists the cases, the checks and their thresholds.
- SSG/LRR-ω Reynolds-stress model (`turbulence = SSGLRR`, `-SD`, `-w2019`, `-w2019-SD`); tests `RSMDecay`, `FlatPlateTMR`.
- Multirate (local time stepping) integrator for time-accurate runs, `[MOSE-Multirate]`; tests in `test/multirate`.
- NASA TMR turbulent flat plate as a verification case (SA, SST, Wilcox 2006, SSG/LRR on three grids).
- SRM grain `502` mode 2: the grain as a mass-flux inlet.
- `diffusive-metric` (`stencil`, default, or `mean`): face metric of the viscous fluxes.
- `irs-variables` (`flow` or `all`): which residuals IRS smooths.
- `prolongated-diter`: spacing of the prolongated outputs of coarse levels.
- `decomposition-map`: MPI rank of each block, from the ATLAS MDB map.
- `ini-error` (`warn` or `stop`): what an error in a re-read `input.ini` does.
- Exact restarts: solutions store the run counters; a split run matches the run in one go bit for bit.
- Coarse multigrid levels write their own wall files (`wall-level<n>`).
- Warnings for keys of `input.ini` that are not read; reports of inflow through `400` faces.
- Probes can record `mi_tilde`.
- Test tiers `needs-mpi` (with a CI job) and `verification`; `MOSE_mms` wall fields; `ViscousOpenPatch`.
- Preconditioning tests: `PrecTransform` (unit), `Preconditioning` (low-Mach cylinder), `InputChecks`/`PrecScope`.

### Changed

- Docs built with Zensical; new home page; Overview merged into Getting Started.
- Theory Guide presents the methods only; input keys, files and messages moved to the User Guide (Models and methods).
- **Breaking:** the `p0_file` of inlet `402` is read in Pa, not bar.
- **Breaking:** the grain `502` line is `mdot_ref, n, p_ref, T_af, mode, massf(1:ns)`; `T_af` is the total temperature.
- **Breaking:** `res-threshold` is relative to the first-iteration residual; it ends steady runs only.
- **Breaking:** the mass flux of inlets `403`, `404`, `420` and the velocity of `408` are positive into the domain.
- **Breaking:** the wall output gives traction, heat and mass flux into the gas on every face (signs change on some faces).
- Local time step: sums the limits over directions, convection and diffusion (smaller steps, stable at the guide's values).
- Time-accurate RANS: explicit turbulence sources, with their stiffness in the time step; point-implicit only in steady runs.
- Viscous fluxes default to the `stencil` face metric, exact for linear fields on any mesh.
- `k-energy-coupling`: walls take k = 0 and open faces carry −2/3 ρk of the boundary cell.
- Choked nozzle `420`: design exit state below the normal-shock pressure; stricter checks of its record.
- `rel_fac` no longer acts on `401`, `402` and `407`.
- Re-reads of `input.ini`: rules checked, fixed keys refused, all ranks take rank 0's text.
- Invalid values in `input.ini` stop the run at set-up.
- Set-up errors stop every rank before the initial condition is read.
- Unsupported records (`409`, `501`, `506`) and inconsistent block connections stop the run at set-up.
- Residual and torque sums are independent of threads and ranks.
- Outputs:
  - `vtk raw` fields are raw binary.
  - `.vtm` files list their blocks correctly.
  - Steady outputs carry −iteration as their time.
  - Text formats write 17 digits.
- ORION, OSLO and FLINT come by release range (`cmake/SemVerDeps.cmake`), fetched at configure: ORION 2.1, OSLO 1.1, FLINT 2.0 or a later release of the same major. OSLO 1.1 has no `dodesol`; ifx with OpenMP refuses FATODE's Rosenbrock solvers (`ros2`, `ros3`, `ros4`, `rodas3`, `rodas4`; use `H-rodas`).
- FiNeR v2.3.0 with its dependencies as submodules (run `git submodule update --init --recursive`).
- CI:
  - Runs on Ubuntu 24.04 and on merge groups.
  - Checks out only the submodules it needs.
  - The fast tier holds only tests of about 10 s or less.
- ifx 2024.1 at `-O2` or above without IPO is refused at configure.

### Removed

- The old turbulent flat plate (`test/2D/viscous/flat-plate-turbulent`), replaced by the TMR case.
- `diffusive-metric = inverse-mean`.
- `turbulence = SAcomp` (unit-dependent).
- IRS in time-accurate runs.
- The submodules `lib/ORION` and `lib/FLINT`, and the CMake variables `ORION_PATH`, `FLINT_PATH` and `OSLO_PATH` (use `FETCHCONTENT_SOURCE_DIR_<NAME>`).

### Fixed

- `prec`: the update no longer clamps Ur at 1e-3 a; runs with Ur below it (small `eps-min`) diverged.
- Turbulence models:
  - SA conservative term −(1/σ)(ν + ν̃)∇ρ·∇ν̃.
  - Point-implicit form in steady runs (the steady state no longer depends on Δt), also with `k-energy-coupling`.
  - Sign of W₂₃ in 3-D, the frame rotation in SA-RC/SST-RC, RC at rest, and the full Ŝ of Wilcox 2006.
  - Wilcox 2006 wall ω and QCR2000 argument order.
  - Wall ghost values in the sources; SA-RC gradients across blocks.
  - `prec` with SST/Wilcox.
- Wall distance is exact (bilinear faces, periodic images, SRM grain); mesh type taken from all nodes; exact cell volumes.
- MUSCL keeps face values between neighbours on stretched cells (Berger–Aftosmis–Murman bound).
- Axisymmetric wedge and 3-D sectors (`200`):
  - Viscous stress on the k faces, and wedge stiffness in the time step.
  - Ghost-cell orientation; consistency checks of the records.
- Boundary conditions:
  - Inflow at total conditions: decided by the reservoir, choked when supersonic.
  - `407`/`405` with `normal`; `420` supersonic root.
  - `301` heat-flux sign on odd faces; `502` gas constant and turbulence ghosts.
  - Faces between blocks computed once (conservative).
  - Fluxes independent of the order of `bc.txt`.
- Ghost cells: wall and `400` fills bounded, linear node extrapolation, `410` values placed correctly, coupled walls `103`/`104`.
- Rotating frame:
  - Stationary faces no longer hang; stationary-wall velocity in the outputs; wall-velocity derivatives and work.
  - Re-read `omega`/`axis`/`origin` taken by every term.
- Time-dependent boundary data evaluated at each Runge–Kutta stage (was first order); tables before their first time.
- Chemistry:
  - Integrator failures stop the run with a message.
  - Missing mechanism or element files stop at set-up.
  - Equilibrium chemistry honours `exclude-blocks` and MPI.
  - Fixed the `exclude-blocks` crash.
- MPI:
  - Stale neighbours at block edges; outputs computed on the owning rank.
  - Tags within 32767; FUNNELED calls from the primary thread.
  - Intel MPI with gfortran.
  - Aborts flush the log; errors found by every rank stop every rank together.
- Restarts:
  - From vtk; probe and `torque.dat` appends; duplicate rows dropped; files flushed.
  - Last outputs written at the end; time handling of steady and timeless solutions.
- Input validation stops at set-up with a located message:
  - Non-finite values (initial condition, `bc.txt`, `410` files, restart counters, rotating frame).
  - Malformed `bc.txt` lines, `402` tables and `410` files.
  - Unknown probe variables; bad `exclude-blocks` words.
  - Degenerate and folded cells.
- Exact Riemann solver wave choice; `shock-detector` values; MUSCL without a limiter; IRS switch logic.
- Outputs:
  - Diagnostic variables named and counted correctly.
  - Well-formed VTK XML; vtk wall block sizes.
  - Wall output independent of `bc.txt` order and leaving the residual intact.
  - Probes by position take the nearest cell; time-based probe rows.
- Data races: `410` edges, chimera edges, symmetry ghost metrics, thread-count print.
- Chimera:
  - Donor-side tangential gradient.
  - Receiver wall distance.
  - Extrapolated ghost row.
  - Inverse-mean metric at block faces.
- Divergence check catches NaN and infinity with `-ffast-math`.
- ifx IPO tools found on oneAPI 2023/2024.
- Test scripts: rocket nozzle vs CEA, finite-rate nozzle, ramp channel, 3-D MMS figure, Shu–Osher GCI, non-OpenMP builds.
- Registry rules checked in full; unreadable `input.ini` reported.
- Docs corrected:
  - The preconditioned update, the Ur modes and the choice of `eps-min` documented.
  - Boundary types `420`, `502`, `506`; multigrid `bc<n>.txt` names.
  - SST production and SA variants; CMake defaults.
  - Loss of order on wedges and zigzag grids.
