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
- A [Versioning](docs/development/versioning.md) page: the releases MOSE requires of its libraries, and how a build gets them.
- Boundary condition `506`, the pyrolysing ablator: the surface reactions of `505` plus the pyrolysis gas blown at `phi` times the char mass flux, steady conduction.
  - Line `pyro_id, reac_id, phi, h_v0, q_rad, eps_wall`; composition from `INPUT/<prefix>pyrolysis.dat` (ATLAS GPB), read by FLINT 2.1.
  - With `phi = 0` and `h_v0 = 0` it is `505` bit for bit; a `506` line of an older ATLAS must be regenerated.
- Tests of the reacting walls: `AblatingWall` runs `505` and `506` on the gas `test/common/CHO-N2-viscous`; unit checks `GSIBradley`, `GSIAblator`.

### Changed

- Tests in three tiers by time: `light` (under 5 s), `medium` (5 s to 100 s), `heavy` (100 s and more), in place of `fast`, `validation`, `verification`; pull requests run `light` and `medium` in parallel, `heavy` runs locally.
- One test per check or part of a check (`TimeIntegration.TimeControl`, `Multirate.Order`, `FlatPlateTMR.SST`); `test/fast` is now `test/checks`; the list of tests in `test/TESTING-PLAN.md` is generated from the registrations.
- `MOSE_mms`, the driver of the manufactured-solution tests, is `test/unit/mms.f90` and is built in `build/test` with the other test programs, no longer in `bin/` (the scripts take `build/test/MOSE_mms` or `MOSE_MMS_BIN`).
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
- Unsupported records (`409`, `501`) and inconsistent block connections stop the run at set-up.
- Residual and torque sums are independent of threads and ranks.
- Outputs:
  - `vtk raw` fields are raw binary.
  - `.vtm` files list their blocks correctly.
  - Steady outputs carry −iteration as their time.
  - Text formats write 17 digits.
- ORION, OSLO and FLINT come by release range (`cmake/SemVerDeps.cmake`), fetched at configure: ORION 2.1, OSLO 1.1, FLINT 2.0 or a later release of the same major. OSLO 1.1 has no `dodesol`; ifx with OpenMP refuses FATODE's Rosenbrock solvers (`ros2`, `ros3`, `ros4`, `rodas3`, `rodas4`; use `H-rodas`).
- FiNeR v2.3.0 is the submodule `lib/FiNeR` (run `git submodule sync && git submodule update --init`); the configure
  fetches its libraries, and stops when FiNeR is not the pinned commit.
- A local ORION, OSLO or FLINT clone outside MOSE's ranges stops the configure (`-DMOSE_PIN_CHECK=WARNING` allows it).
- CI:
  - Runs on Ubuntu 24.04 and on merge groups.
  - Checks out only the submodules it needs.
  - The fast tier holds only tests of about 10 s or less.
- ifx 2024.1 at `-O2` or above without IPO is refused at configure.

### Removed

- The turbulent flat plate with boundary conditions of its own (`test/2D/viscous/flat-plate-turbulent`, SA, SST,
  Wilcox 2006 and SSG/LRR, outside CTest) and its V&V page: the flat plate of the NASA TMR, with the boundary conditions
  of the TMR, covers the four models on three grids of the same family.
- The value `inverse-mean` of `diffusive-metric` (`[MOSE-Numerics]`): the key takes `stencil`, the default, exact for a
  linear field on any mesh, or `mean`. A deck that names `inverse-mean` stops at set-up. The inverse-mean face metric stays
  where the code takes it whatever the key: the chimera faces, and with `mean` the faces between blocks and of type `200`.
- Implicit residual smoothing (`irs = true` or `irs-beta > 0`) in time-accurate runs: it is refused at set-up; steady runs keep
  it. The smoothed residual is not the residual of the space discretization, so a run with it did not follow the transient of
  the discretized equations: in one step from the same state the time derivative changed, by an amount that did not depend on
  the step, and so did the final density of a Sod run. The change does not vanish as the step does. MOSE has no dual time
  stepping, in whose pseudo time residual smoothing belongs.
- The submodules `lib/ORION` and `lib/FLINT`, the CMake variables `ORION_PATH` and `FLINT_PATH`, and `OSLO_PATH`, which
  `install.sh` handed to FLINT: `FETCHCONTENT_SOURCE_DIR_ORION`, `FETCHCONTENT_SOURCE_DIR_FLINT` and
  `FETCHCONTENT_SOURCE_DIR_OSLO` give a local checkout instead.
- The submodules of `lib/third_party` and the CMake variables `FINER_PATH` and `FINER_DEPS_PATH`:
  `FETCHCONTENT_SOURCE_DIR_FINER` (`--include-finer`) gives an external FiNeR.

### Fixed

- The wall state of a burning SRM grain (`502`) takes the gas constant of the composition the grain imposes: it took that
  of the gas of the boundary cell, so that its density, and with it the blowing velocity, its momentum and kinetic energy,
  were off by R(cell)/R(products) wherever the cell held another gas (a chamber not yet filled with the products, a
  multi-species run: 0.80 of the right blowing velocity for products of 371 J/(kg K) blown into N2). The mass flux and the
  enthalpy of the injected gas do not change; in a single-species run the gas constant of the cell is that of the
  products, and this change alone changes nothing.
- A property line of the SRM grain (`502`) stops the run at set-up when its mode is neither 1 nor 2, its `p_ref` is not
  positive, its `mdot_ref` is negative (or 0 in mode 2), and when `T_af` is not positive (the boundary-condition builder of
  ATLAS writes 0 for a grain without `eq-CEA-file`) or does not lie on the tables of the gas, from their lowest temperature
  above 0 K to one kelvin below their last row: the tables hold the enthalpy of `Tmax - 1` above it, so the gas entered
  with that enthalpy while the wall stayed at `T_af`, and a `T_af` of 0 made the grain an isothermal wall at 0 K (a slip
  wall in an inviscid run).
- The MUSCL reconstruction keeps the value of a cell at each of its faces between the value of the cell and that of the
  neighbour across the face on lines of cells of any lengths, with every limiter (van Albada excepted where the two slopes
  have opposite signs): the limited slope is bounded so that its change over half the cell is at most the jump to each
  neighbour, the condition of Berger, Aftosmis and Murman (AIAA 2005-0490, eq. 19). Where a cell was longer than its
  neighbour across a face, van Leer, MC and superbee could take the face value beyond that neighbour, by up to
  2c/(1 + c) - 1 of the jump for a cell c times longer: a variable about zero there, such as k at the edge of a boundary
  layer, became negative at the face, the cell could go below zero after a stage and be set to the realizability floor, and
  the positivity of the transport that the guide gave for `cfl = 0.5` held only on lines of equal cells. On lines of equal
  cells the slopes do not change, bit for bit in builds without fast math: there the bounds are 2 and 2r, which no limiter
  exceeds. Elsewhere a slope changes only where a bound acts, where the ratio of the two slopes is away from 1 by an amount
  that the ratio of the lengths sets, so data linear along the line keep their slopes; minmod and first order do not
  change. In a `RELEASE` build, compiled with `-ffast-math`, the compiler is free to evaluate the changed code in another
  order, and face values differ in their last digits also where no bound acts (measured with ifx 2024.1 on lines of equal
  cells: at most 4 units in the last place of the values of the two cells of the face). The
  guide gives the condition of positivity of the transport, a Courant number of the flow leaving each cell of at most 1/2,
  with a lower bound in the first cell next to a wall and none next to faces whose ghost cells extrapolate the cells
  (`408`, `502`).
- A failed integration of the finite-rate chemistry stops the run with exit code 1 and a message that gives the integrator,
  its return code and what the code means, the block (with the grid level in a multigrid run), the cell, the chemistry step
  and the temperature, pressure and density of the cell; every thread, and with MPI every rank, stops at the end of the
  chemistry of that step, before any flux reads the cell, and no solution of that step is written. MOSE did not read the
  return code of the integrator: the cell took the state where the integrator had stopped, or kept the one where it had
  started, and the run went on to its end with exit code 0 (with `H-radau5`, the default, only the line of RADAU5 told it).
  The message is written once, by the rank that holds the smallest block with a failed cell, for its cell of smallest k,
  j, i, so that it is the same with any number of threads and ranks. Runs whose integrations succeed give the same bits;
  the chemistry of each finite-rate step ends with two barriers of the threads and, with MPI, a reduction of one number
  over the ranks.
- With `chemistry = finite-rate` and no `INPUT/chemistry-info.txt`, the reaction mechanism, the run stops at set-up with a
  message that names the file, before the initial condition is read, on every rank. The run went on: without the file no
  rate table is read, and FLINT answers every evaluation of the source terms with its bail-out value, -1 in every
  component, which every integrator of OSLO but `cvode` takes as a rate, so that each chemistry step changed the partial
  densities and the pressure of every cell, with no message after the set-up, where the warning of FLINT and the set-up
  summary gave an empty mechanism name. Runs with `frozen`
  or `equilibrium` chemistry, which need no mechanism, and finite-rate runs with the file do not change. The tests of the
  finite-rate chemistry in `ExitCodes` and `MPIExitCodes` take the gas and the mechanism of `test/1D/Fer14`.
- The compressibility correction of Spalart–Allmaras (`turbulence = SAcomp`, Paciorri and Sabetta, J. Spacecraft Rockets
  40(3), 2003): the local convective Mach number M̃c is taken from M̃c² f₂(M̃c) = 25 ν_t |ω| / a², the paper's eq. 12, where the
  code divided by a instead of a²; the right side was a speed, so the factor on the production depended on the unit of
  length (0.99 with metres, 0.22 with millimetres for Ω = 1e3 1/s, ν_t = 1e-3 m²/s, a = 340 m/s) and cut the production
  wherever the shear was strong, whatever the Mach number (78 % at Ω = 1e5 1/s). With a² both sides are numbers; the same
  state gives a factor of 1.00. The set-up accepts the model again (it had been refused since the error was found); the
  Newton iteration for M̃c is capped at 100 steps. The constants (25 = 1/(4 τᵢ), τᵢ = 0.01), f₁, f₂, the eddy viscosity
  ν̃ f_v1, the vorticity and the local speed of sound follow the paper. No case or test of the repository used the model.
- The testing guide and the test plan say what the tests take and give: the fast tier holds parts of about 10 s of CPU time or
  less at one thread each, the rule of a row of the plan (the guide said "a few seconds each"); `MMSViscous3DFull` takes 18 h
  at one thread with a RELEASE build on one core of a shared cluster node and about 1 h 10 min with 18 or 32 threads (the
  plan, the guide and the comment of its `TIMEOUT` gave an estimate of about a day); the results of the cases outside CTest
  and of the isentropic vortex are those of the code of the repository (the plan named the commit of MOSE they came from).
  The plan gives the times of `AblatingWall`, `Noh87` and `PrandtlMeyer` without a cause: it put them down to the time step
  alone, while `Noh87` takes the same steps, each dearer, and `AblatingWall` also runs one more column. No test or
  registration changes.
- The face of a choked nozzle (`420`) below `p_sup` takes the supersonic design exit state of its record, the supersonic root
  of the mass flux `g` on the isentrope of the reservoir, at every exit Mach number. From an exit Mach number of about 2.8 in
  air, and of about 3.6 in the gas of the frozen rocket nozzle of the tests, the Newton iteration of the face reached the
  subsonic root, which let several times `g` into the domain (8.4 times at Mach 4 in air), or left the tables and stopped the
  run (an index out of bounds, or a segmentation fault without bounds checks). The state is computed once at set-up, with a
  Newton iteration kept inside the bracket of the supersonic root.
- Faces between blocks take the inverse of the mean of the inverse metrics with `diffusive-metric = mean`, and
  chimera faces with every value of the key: across a large size jump at a block connection the mean of the two metrics made the
  normal gradient far too steep, and the explicit diffusion unstable on the large cell.
- `level<n>-iter`, and `iter-threshold` without multigrid, take effect when `input.ini` is read again during the run: the
  iterations of each level were copied once at set-up, so a changed value did nothing. Lowering one below the iterations
  already done ends that level at the iteration after that reading.
- `irs = true`, or `irs-beta` above 0, added to `input.ini` during a run started without IRS is an error of the reading
  (with `ini-error = warn` the reading is not taken, with `stop` the run stops): the arrays of IRS exist only when it is on at
  set-up, and the run used them unallocated.
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
  QCR2000, as the guide says, so no run changes.
- The fluxes of the boundary faces no longer depend on the order of the records in `bc.txt`, and each MPI rank walks only the
  records of its own blocks. The boundary fluxes were computed in one worksharing loop per range of records, in file order, on
  every rank: the face contributions of a corner cell were summed in the order of the file, and with the records out of order
  and more than one thread two threads could update the same corner cell inside one loop, a data race whose result is
  unspecified; every rank also read every boundary record of the mesh at every stage. The records of the local blocks are now
  grouped by block and face, keeping the file order within a face, with one loop per face of a block. A `bc.txt` sorted by
  block and face, as the guide lists the records, gives the same solution bit for bit; any other order gives the solution of
  the sorted file.
- A face between two blocks (block cut `101`, periodic connection `201`) is computed once, and its flux is given with opposite
  signs to the cells on its two sides. Each block computed the face with its own inputs (the shock sensor of its own cell and,
  with HLLC-PC, a reference velocity of zero in the ghost cell): mass, momentum and energy were not conserved at the face, and
  a mesh split into blocks did not give the solution of the unsplit mesh. The record of the cell on the low-index side computes
  the flux (between faces of the same parity, the record with the smaller block, face, k, j, i), whatever the order of the
  records; between MPI ranks the share of the other cell travels in one more message per stage. With
  `preconditioning-Ur-smooth` > 0 the smoothing next to such a face sees the cell across it. Results change where a shock
  sensor or the preconditioning acts at a face between blocks, and at round-off where the two blocks store the nodes of the
  common face with different digits or a periodic connection translates them.
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
- A supersonic inlet (`405`) whose flow direction is given as `normal` lets the gas in along the normal of the face into the
  domain, on every face; with one angle `normal` and the other a number, that angle is the azimuth or the elevation of this
  normal. It took the angles of the normal toward growing index, atan(n_y/n_x) and atan(n_z/n_x), and built a velocity
  (cos α cos β, sin α cos β, sin β) whose x component is never negative: where the normal into the domain has a negative x
  component (face 2 of a block whose i grows along x, face 1 of one whose i grows along −x) the gas left through the inlet,
  0.463 ρa out of the face of a tube at Mach 2 instead of 1.968 ρa in, and where it has components along both y and z the gas
  entered off it, by 9.7° for (1, 1, 1)/√3. Results change where the angles taken, atan(n_y/n_x) for `alpha` and
  atan(n_z/n_x) for `beta`, are not the azimuth and the elevation of the normal into the domain, for the angles that the
  record gives as `normal`; where they are, as on a face whose normal into the domain is (1, 0, 0), at most the last digits
  change. Records with angles are unchanged.
- A 1-D case (one cell in j and in k) with symmetry (`300`) on all its j faces, or on all its k faces, no longer fails on
  round-off in v or w: the two mirror planes one cell apart give v the mode dv/dt = −2av/Δy, outside the stability region of
  the explicit schemes at the step of the x cells (three such cases stopped with `NaN or p<0` after 14 and 15 steps). Such a
  velocity is now set to zero at set-up and not integrated, with a warning that names the faces, when those faces are planes
  normal to its axis (a tube of constant section); in a channel of variable section it is integrated as before, with a
  warning, and thin cells can still stop the run. Results change only in these 1-D cases.
- The volume of a cell is the integral of det J of the trilinear cell, whose faces are the bilinear faces of the fluxes. It was a
  sum of five tetrahedra, which cut each face along a diagonal while the cell next door cuts it along the other: on a twisted
  face the volume of a cell was not the volume its faces enclose (up to 18.5 % too large on a cube with one twisted face), and
  the volumes of a mesh did not sum to its volume. On cells with plane faces (2-D and axisymmetric meshes, 3-D meshes of plane
  faces) results change in their last bits. A cell whose volume has not the sign of its block, a cell turned over away from its
  centre, stops the run at set-up with its block and cell, on every grid level of a multigrid run.
- The wall distance of the turbulence models is the exact minimum distance from the cell centre to the wall surface, a wall face
  being the bilinear surface through its four nodes, and on a 2-D axisymmetric wedge the distance in the meridional plane. It
  was the distance to the nearest wall-face centre, too large wherever the grid lines are not normal to the wall and next to
  convex corners and the ends of walls (up to 16.3 % in the test geometries). A twisted 3-D wall face takes a certified minimum
  over the face, exact to the rounding of its size. Results change wherever the distance changes, and at round-off where the
  grid lines are normal to a plane wall. Equal distances go to the lowest face index, so the roughness a cell inherits is the
  same on every rank and thread. With a turbulence model, a wall on faces 5–6 of a wedge stops the run. On a 3-D mesh with
  twisted wall faces the start-up takes longer: 26 % at one thread on a mesh of 1e6 cells whose 1e4 wall faces are all
  twisted (ifx 2024.1), within the spread of the measurements with 16 threads; meshes whose wall faces are plane take the
  search of before.
- The type of the mesh (plane, axisymmetric wedge, 3-D sector closed by type `200`, slab translated in k) is taken from every
  node of every block, to within 1 % of the move of the nodes per cell in k. Two probe nodes of block 1 took a plane 2-D mesh
  whose j index runs towards the x axis, and a plane slab far from the x axis, for a wedge (wall distance off by 27 % in a
  test), turned the periodic images of a 3-D slab periodic in k, ran a mesh with a wedge block and a plane block as one of the
  two, and took a 3-D mesh with a thin block 1 for a 2-D one. A mesh whose blocks disagree, a 2-D block whose k planes fit
  neither a turn nor a translation to that tolerance (too few digits, or k planes that are not parallel), and a 3-D block
  closed by type `200` whose last k node plane is neither the first one turned nor translated, stop at set-up; the set-up
  summary names the type. A block closed on itself over a full turn about the x axis is translated by zero, the periodicity by
  the identity. Meshes recognised right keep their wedge angle to the bit.
- The wall distance is the minimum over the images of the wall across periodic boundaries too (periodic connections `201`, slabs
  closed in k by type `200`, rotational periodicity of 3-D sectors): it was measured to the wall inside the mesh only, and near
  the end of a period it could be several times too large. A mesh without periodic boundaries costs nothing more, and the images
  searched are bounded whatever the size of the mesh. With a turbulence model, periodic translations that cannot be the images
  of one wall stop the run with a message. Translations that differ by at most 1 % of the thinnest cell across the periodic
  boundaries are the same, as the two sides of a periodic boundary match only to the round-off of the coordinates: a channel
  periodic in x, turned by 30 degrees and written with 4 significant digits, runs.
- On a 3-D mesh closed in k by type `200` (a sector turned about the x axis, or a slab whose last k node plane is the first
  one translated), the ghost cell beyond face 5 holds the state of the cell next to face 6 turned back by the sector angle
  (copied, on a translated slab): its w was the opposite of the turned one (face 6 was right). The second-order
  reconstruction of the faces next to face 5 reads that ghost cell, and so do the viscous fluxes of the cells next to it;
  a radial flow that is its own mirror image departed from it by 1.8e-2 of its scale in 20 steps. Results change on 3-D
  sectors wherever the velocity next to face 6, turned by the sector angle, has a z component, and on translated slabs
  wherever w is not zero.
- Type `200` joins the two k faces of one block: a record of type `200` on a face 1 to 4, with its cell outside its block, or
  with no record of type `200` on the other k face of its block in its column stops the run at set-up, on every rank, with a
  message that names the record and the remedy. A 3-D sector cut along k into blocks, a ring of 360 degrees in several blocks
  closed by type `200` between the last block and the first, type `200` on a face 1 to 4 of a 2-D wedge, and a wedge with
  type `200` on one k face only ran in silence with wrong fluxes, and a record outside its block read outside the arrays. Cut a
  sector along i or j, and join the blocks of a ring by `101`, the junction at 2π included. Declare the axis of a wedge or
  of a sector, on one of the faces 1 to 4, symmetry (`300`): the boundary-condition builder of ATLAS writes type `200` for
  `axisymmetric`. The mesh decomposer of ATLAS (MDB) cuts blocks along k unless `split-directions` says otherwise: give it
  `split-directions` without k for a block closed by type `200` (see the guide). Other meshes do not change.
- The inlets at total conditions (`401`, `402`) and the ambient boundary (`407`) let the gas in when the reservoir can push it
  in, when p0 is above p_stop, the pressure that the gas of the boundary cell reaches when the outgoing wave brings it to rest.
  The choice followed the sign of the velocity of the cell: a tube at rest never started from a `401` or `402`, a `407` at rest
  computed an inflow outside its own domain, and an inlet with p0 below p_stop let gas out. When the gas does not enter, a
  `401` or `402` face is a slip wall and a `407` face an outlet at p_back. Results change where the sign of the cell velocity
  and the criterion disagree.
- An inflow at total conditions whose face state would be supersonic across, built on an invariant that can no longer reach
  the face, takes the choked state, normal Mach 1. Past a pressure ratio (p0/p = 6.12 for a cell at rest at the reservoir
  temperature, normal inflow) the subsonic equations gave such a face: at p0 = 36 p normal Mach 1.85 and a mass flux of
  13.9 ρa, below the 20.8 ρa of the choked face. Results change for inlets at high pressure ratios, `401` and `402` included.
- The inlets at total conditions (`401`, `402`) and the ambient boundary (`407`) let the gas in when the gas of the boundary
  cell moves in at 2a/(γ−1) or faster (normal Mach 5 or more at γ = 1.4): no pressure brings that gas to rest, so every
  reservoir with p0 > 0 feeds the face. Such a face was a slip wall (`401`, `402`) or had no flux (`407`), in silence. The face
  takes the inflow state when it exists and is subsonic across, the choked state otherwise, and the mass flux is continuous
  where the outgoing invariant changes sign. Results change only at faces where the gas of the cell moves in that fast.
- With `sol-format = vtk` the wall output sets the block sizes of its files, which it never set: depending on what they held, a
  run with a viscous wall could stop in the wall output or write wall files of no extent and one point. Each vtk wall file now
  has the extent and the points of its face. Tecplot wall files do not change.
- With more than one OpenMP thread per MPI rank, MPI is called only from the thread that initialised it, as
  `MPI_THREAD_FUNNELED`, the thread level MOSE requests, requires. Inside the parallel region of the Runge–Kutta stages the halo
  exchange and the exchange of the smoothed residuals were wrapped in `single` constructs, which any thread of the team may
  execute: with more than one thread per rank, MPI calls could come from threads other than the primary one. They are now
  `master` constructs, which the primary thread executes, each followed by the barrier that the `single` construct implied
  (none after the block that posts the sends, which was `single nowait`). Results do not change.
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
  (up to 256); a slash ended the read of the property line, and a long name was cut. A table that cannot be opened, holds no
  (time, value) pair, has decreasing times or holds a time or a value that is not finite stops the run at set-up with the name
  of the file, and the line and the column of the value: a missing table was not reported, and the others ran in silence to
  exit code 0 (a NaN pressure closed the face of the inlet).
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
  and the name of the value: such values closed a face in silence, ran to the end with exit code 0, or stopped later at the
  divergence check. The messages of the reader name the file it reads, which on a coarse multigrid level is the file of that
  level, not `bc.txt`.
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
  not read.
- A file of an assigned state (`410`) without `p`, without one of `u`, `v` and `w`, or without any `rho(N)` stops the run at
  set-up with exit code 1, with the file and the variables it lacks: the file gives the whole state of the cell across the
  face. The ghost cells took the value of the cell for a missing `p`, `u`, `v` or `w`, and the run ended with exit code 0;
  without any `rho(N)` they had zero density, and the run ended at its first step with a run-time error (gfortran with bounds
  checking).
- A file of an assigned state (`410`) that leaves a cell of a face of type `410` without a value at one of its times stops
  the run at set-up with exit code 1, in every build, with the file, the block and the record: a block with no zone in the
  file, with no zone at one of its times, or whose zones at a time do not hold the cell, as when the blocks of the mesh were
  split or numbered anew after the file was written. The reader left such a cell with the zero state: the run ended at its
  first step with a run-time error that named no file (an index out of bounds, with gfortran with bounds checking and with
  ifx), or took that state at that time and ended with exit code 0.
- With the rotating frame, `stationary-face-<n>` entries no longer make the start-up hang: the count of the entries asked for
  `stationary-face-1` at every pass, so stationary walls could not be used. A `rotational-frame` other than `rigid-body`, its
  default `none` included, turns the frame off before its entries are read; with `none` the switch of the frame was left
  undefined.
- The user guide, the theory guide and the text of the key `diffusive-metric` describe these points as MOSE does them. With MPI
  some errors that one rank finds alone (at set-up in the mesh of its blocks or in its probes, on rank 0 in the diagnostic
  output), and two errors of `bc.txt` that every rank finds, stop the run through `MPI_ABORT`,
  whose best attempt to stop the processes can lose what the other ranks hold for the log (the guide said that every
  error found by every rank stops the ranks together). A binary Tecplot file in a build without TecIO stops the run inside
  ORION with `ERROR STOP` (the guide said that the run ended with exit code 0). MDB writes its decomposition map by default as
  `decomposition.map` in the folder where it runs (the guide said next to the split mesh). ATLAS writes `<name>-bc.txt` for a
  named phase, and MOSE reads `INPUT/bc.txt`. The outlet `406` with `p_back` = 0 and `rel_fac` = 1, what ATLAS writes for an
  `outlet` without `p` and `rf`, takes the state of the boundary cell whatever its Mach number, except where that gas moves
  into the domain (a slip wall). The flow angles give v = u tan(alpha) and w = u tan(beta), and the velocity of `405` the
  direction (cos(alpha) cos(beta), sin(alpha) cos(beta), sin(beta)). `rel_fac` has no effect on `403` and `420` either. With
  `diffusive-metric = stencil` the symmetry faces with one cell across them, the
  walls `502` and the walls `503`-`505` while they ablate take the metric of the cell, and with `mean` the faces of type `200`
  take the inverse-mean metric, as those between blocks. The coupled walls `103` and `104` share the declared limit of the
  wall values of the turbulence variables of `503`-`505`; the wall flux of the walls that blow takes C = 800 in the wall value
  of omega with either value of `omega-wall-bc`, and the guide described a blowing correction of omega that no model name
  that MOSE accepts selects. The residuals do not show the cells that the turbulence floor holds from one step to the next.
  The home page builds with `./install.sh build`. No run changes.
- The user guide gives the data lines of the boundary types `420`, `502` and `506` as MOSE reads them: the fifth value of `420`
  is the mass flux `g` through the face, not a throat radius; the line of `502` ends with the mass fractions of the products;
  `506` was missing.
- The guide names the boundary files of the multigrid levels as MOSE reads them, `bc.txt`, `bc2.txt`, `bc3.txt`, … (it said
  `bc_x1.txt`, `bc_x2.txt`, …).
- The text of the key `sst-production`, the turbulence guide and the SWBLI page describe the `compressible` production as MOSE
  computes it: the exact production with the dilatation terms arranged as in OpenFOAM's `kOmegaSST`, with the limiter on its
  deviatoric part only and -2/3 gamma rho omega div(u) in the omega equation. They named it NASA TMR's SST-2003e (SST-2003 with
  `k-energy-coupling`), which limits the whole production and takes gamma/nu_t times it in the omega equation. The guide also
  writes the Boussinesq stress with its -1/3 div(u) delta_ij part, and says that its -2/3 rho k delta_ij term enters the momentum
  and energy equations only with `k-energy-coupling = true` (it said that the two-equation models include it), that `SA` is the
  model without the f_t2 term (SA-noft2 in TMR naming), with the limiter of S~ of Note 1(c) of the TMR page and the cap r <= 10
  that the code applies (the guide gave neither), and that with `compressible` the rotation/curvature factor of SST-RC multiplies
  only the deviatoric part of the production. The SWBLI page gives the whole variant of its two models, and the models of the
  Wind-US and SU2 solutions as their sources give them: SPALART and SST, without variant, in the NPARC Alliance archive, and SA-neg
  and SST-2003m in the SU2 validation study of the case. No run changes.
- The comparison of the rocket nozzle with CEA (`test/2D/euler/rocket-nozzle/verify.py`, page `docs/vv/2D-nozzle.md`) reads the
  solutions of the runs of the three configurations, `OUTPUT/field.tec`, and stops with a message, exit status 1, when one is
  missing; it drew the stored solutions `reference/field.tec`. The mixture density is the sum of every species density: it summed
  12 of the 13 species of the frozen and equilibrium configurations and 18 of the 19 of the finite-rate one, so that the Mach
  number of the frozen run was 0.39 % low. The script prints the averages at the CEA stations and their differences from CEA; the
  page calls the averages weighted by the mass flow (it called them area-weighted) and names the run of the Mach contours,
  finite-rate, and the panels of the comparison.
- The comparison of the rocket nozzle with CEA (`test/2D/euler/rocket-nozzle/verify.py`, page `docs/vv/2D-nozzle.md`) takes
  the Mach number of the CEA equilibrium stations with the frozen speed of sound, as for the runs of MOSE and the CEA frozen
  stations, and the averages of the runs at the area ratios of the stations, interpolated linearly in x between two columns of
  cells, with the areas relative to the narrowest section of the grid. It divided the velocity of the equilibrium stations by
  the equilibrium speed of sound of CEA, so that one state had Mach number 1 for CEA and 0.974 in the definition of the curves
  of MOSE, and it took the column of cells nearest to each area ratio, relative to the middle of the narrowest column, 0.16 %
  wider than the narrowest section: at the subsonic station the column had an area ratio of 1.883 instead of 2. With the runs
  of the case, the Mach number of the equilibrium run at the subsonic station is 13.3 % above CEA instead of 16.7 %, and its
  pressure in the diverging part 5.8 % to 11.9 % above CEA instead of 7.0 % to 16.4 %.
- The comparison of the rocket nozzle with CEA (`test/2D/euler/rocket-nozzle/verify.py`) takes u, v, p, T, γ and R by their
  names in the solution files, and stops with a message, exit status 1, when one is missing. It took them at fixed places after
  the species densities, those of the files that the case writes: a file with the variables of another model between p and T
  (soot, passive scalars, turbulence), or one written with `sol-variables = transport`, which writes `mil` and `kl` where
  `thermo` writes `g` and `R`, gave it other quantities, with no message.
- The figure of the grids of the 3-D manufactured-solution page (`docs/vv/images/mms3d-grids.svg`, drawn by `verify.py` of
  `test/3D/viscous/mms`) shows the iface grid, whose panel was empty: `verify.py` drew the block faces on the plane z = 0, and the
  iface grids, kinked in z too, have none. It draws the bottom face of each block, the face that the period of the grid in z
  moves onto a face of a block, and stops with a message, exit status 1, when a block has no such face or more than one.
- The script of the ramp channel (`test/2D/euler/ramp-channel/verify.py`) finds the ORION Python package in `lib/ORION/src/python`
  of the repository, as the other `verify.py` scripts of the tests do: it imported ORION from the Python path alone, and stopped
  with `No module named 'ORION'` unless ORION was installed or on `PYTHONPATH`.
- The verification script of the finite-rate rocket nozzle (`test/2D/euler/rocket-nozzle/finite-rate/verify.py`) takes p, u, v,
  T and γ by their names in the solution files, and the mixture density as the sum of every species density. It took them at the
  positions they have in a file with 13 species, and the case writes 19: at the throat it compared four species densities and the
  sum of 12 of the 19 instead of p, ρ, u, v and γ, and two of those densities lie below the absolute tolerance of the comparison,
  1e-8, so that any value of them up to it passed; its figures read the stored solution `reference/field.tec` instead of the run,
  at wrong positions too. The case passes the check: at the throat p, ρ, u, v and γ of the run are within a relative 1e-3 of
  `reference/field.tec`.
- The Grid-Convergence Index of the Shu–Osher study (`test/1D/Shu-Osher/convergence.py`) is that of Roache: the safety factor
  1.25 times the Richardson estimate of the error of the finest grid, which is the normalised L1 error of that grid against the
  extrapolant. The script divided that estimate by r^p − 1 once more (1.02 for the grids of the case). The Shu–Osher page gives
  0.90 % (it gave 0.88 %). The criteria of the test, the observed order and the decrease of the errors, do not change.
- The tests that check the number of OpenMP threads of their runs pass in a MOSE built without OpenMP (`USE_OPENMP=OFF`, the
  default of the configuration), which writes `Serial execution` where it writes the number of threads with OpenMP; they failed
  there on that line. `BoundaryConditions` skips the part of `Extrapolation400` that compares 1 and 4 threads, `MPIBitIdentity`
  compares 1×1, 2×1 and 4×1 and does not run 1×4 and 2×2, and `MPIThreadFunnel`, which counts the MPI calls of the threads
  other than the one that initialised MPI, exits with 77 after its first run, the code its CTest definition declares as that of a
  skipped test (`SKIP_RETURN_CODE`); each says so in its output. With OpenMP they run as before.
- Every rule of the parameter registry is checked, at set-up and when `input.ini` is read again: intervals such as the `in (0,
  1)` of `preconditioning-eps-min`, and the rule of a list of numbers on each of its values (`index-position` of a probe, `>=
  0`; `ode-relative-tol-species` and `ode-absolute-tol-species`, `> 0`). Only the bounds of single numbers and the lists of
  names were checked, and a value outside the other rules ran with no message. A rule of another form stops the run as an error
  of the registry. Runs whose values keep to the rules do not change.
- A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity) stops the run at set-up, in every
  build, with a message that names the key (the run went on to the divergence check, or stopped with the message of an axis of
  zero length, or with SIGFPE in a DEBUG build). So does a real parameter of `input.ini` read as an infinity (gfortran reads
  `1e999` as one), with a rule or without, at the start, and read again during the run it is an error of the reading (see Added):
  a bound such as `> 0` let it through.
- An `input.ini` that cannot be read stops the run at set-up with the reason, on every MPI rank and with exit code 1, and read
  again during the run it is an error of the reading (see Added): a file that cannot be opened, that holds no line that is not blank, or that the INI reader
  rejects. A missing, unreadable or comments-only file stopped the run on the first required key without naming the file, an
  empty one ended the run inside the reader, and a file that could not be read again during the run left the values of the
  previous reading in place, without a word.
- With MPI, an error that every rank finds stops every rank together, and the message of every rank reaches the log: an error in
  a time table, in a line of `bc.txt`, in a value of `input.ini` read again, in the mesh type or in the periodic images of a
  wall. So does an error that rank 0 alone finds at set-up (the data tables, the rules of `input.ini`, the boundary conditions),
  once rank 0 has written it, and an initial condition that one rank cannot read, which the first rank that found the error
  writes, with its rank. Every rank called `MPI_ABORT`, which makes only a best attempt to stop the processes, so at times no
  line reached the log (and with Intel MPI 2021.12 a run could end with exit code 141); rank 0 decided alone on the read of the
  initial condition, so a run whose rank 1 could not read its copy went on to its end.
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
  why, or, read again during the run with `ini-error = warn`, is not taken on any rank; the `input.ini` in the folders of the
  other ranks is not read.
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
  `stationary-face-<k>` key that names another face or that the set-up does not have, is an error of the reading (see
  Added), on every MPI rank together: they were ignored without a word. The `stationary-face-<k>` keys that the file read again
  leaves out keep their faces, as every key left out keeps its value.
- A value of `axis` or `origin` of the rotating frame that is not finite (NaN or infinity), written in `input.ini` during the
  run, is an error of the reading when the file is read again (see Added), with a message that names the key; in a gnu TESTING
  build it went into every term of the frame and the run failed at the next iteration; a DEBUG build stopped with SIGFPE at that
  reading, an ifx RELEASE build with the message of an axis of zero length.
- A value of `axis`, `origin` or of a `stationary-face-<k>` key of the rotating frame that does not read as numbers, written in
  `input.ini` during the run, is an error of the reading when the file is read again (see Added), with a message that names the
  key and its value: it ended the run with an error of the run-time library at that reading.
- A run that ends where a numbered or time-based solution or a probe row is due writes it, with its wall and diagnostic files,
  before the final solution, whose counters follow it: the last iteration wrote the final solution only, so the run in one go
  lacked its last numbered solution and probe rows, and a restart from a first part that ended there numbered the solutions
  after it one lower, shifted the time-based ones and left the state of that iteration in no file. The guide no longer lists
  this as a limit of the restart.
- With `diffusive-metric = stencil` the gradient at a wall face takes, along the wall, the differences of the wall values with
  the same differences of the positions, so that it is exact for a linear field whose wall values are its values at the face
  centres. The derivatives along the wall were zero and the metric was that of the wall cell: where the wall temperature varies
  along a skewed wall the normal gradient was off by (dT_w/ds)·tan θ at any resolution, and on wall cells that are not
  parallelograms the velocity gradient was wrong too. Results change only there.
- With `diffusive-metric = stencil` the source terms of the turbulence models (Spalart–Allmaras and its variants, SST,
  Wilcox 2006, SSG/LRR) and the eddy viscosity of the output take cell gradients whose metric comes from the same differences of
  the positions as of the values, with each ghost value at the point where it is right. With the metric of the cell the gradient
  of a linear field was wrong next to a jump of cell size (5 % in the wall cell of a layer growing by 1.2), inside a stretched
  layer, behind walls on skewed cells (a copied density took a normal gradient it does not have) and next to zero-gradient faces
  (half the one-sided difference). Results change for RANS runs, mostly in the cells next to boundaries.
- With `diffusive-metric = stencil` a mesh with faces of zero area, such as an axis at r = 0 or the polar line of a 3-D mesh
  around its axis with symmetry on it, no longer stops at set-up on a not-a-number: such a face has a zero normal, as in the
  rest of the solver, its ghost cell is the cell itself, and its flux is zero.
- The ghost cells behind viscous walls (`301`, `302`, `503`–`505`) and extrapolation faces (`400`) are filled once per stage,
  for the fluxes and the source terms alike, from the quantities whose differences the fluxes take: with a wall value q_w the
  ghost takes 2q_w − q₁ (the velocity of the walls and the turbulence variables), the temperature of a `302` wall
  max(2T_w − T₁, T_w/2), and the other quantities an extrapolation kept within [q₁/2, 2q₁]; the density follows from
  p/(R T). The extrapolation 3q₁ − 3q₂ + q₃ of every primitive variable could turn negative and amplify (a negative ghost
  density next to a hot wall), put the wall value on the face only when the profile was linear over three cells, and the source
  terms read a second set of ghost values. A run in which the cut or the floor acted counts those ghost fills in a line at
  its end. Results change for every run with a wall or a `400` face, in the first rows of cells next to them.
- The ghost cells behind the boundary faces other than `101`, `201` and `300` take their nodes by a linear extrapolation along
  the index lines, so that both layers keep the edge across the face of the first cell. The cubic extrapolation made a ghost
  layer vanish where the cells grow 1.5 or 2 times across the face: a run stopped on a zero determinant, or went on with a
  metric as large as 1/det, and at chimera faces (`102`) the viscous flux could take an infinite or reversed value. A cell whose
  determinant is below 1e-12 of the product of its edges, or a ghost cell of the first layer behind a chimera face turned over
  or degenerate against the first cell, stops the run at set-up with its block, cell, face and ratio, and the remedy. Behind
  the other faces no flux reads the metric of the ghost cells. Every run changes in the last bits next to such faces, and by
  O(r − 1) where the cells are stretched across them.
- A cell folded at one of its corners, or turned against the orientation of its block, stops the run at set-up, naming block,
  cell, corner node, the Jacobian there and the bound; the coarse levels of a sequence of grids are checked too. Such a cell
  passed in silence: the volume of a cell was a sum of tetrahedra taken in absolute value, and the cell check looked at the
  determinant of the mean edges only. On a valid mesh nothing changes, bit for bit.
- The two k faces of a 3-D sector closed by type `200`, one interior face of the whole annulus, carry its viscous flux (viscous
  stress, heat flux, diffusion of species and turbulence variables), computed once between the cell next to face 6 and the
  turned image of the cell next to face 5. They carried the inviscid flux only, so the solution tended to that of a sector
  whose periodic planes are free of traction and adiabatic, an error whose share grows as the mesh is refined. Runs on the
  wedge do not change.
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
  that names the key and the record: such a face was a boundary condition of no wall, or did nothing. A stationary face may be
  shared with the types whose boundary condition does not read whether their face is stationary (the guide lists them), as a
  casing face that is a wall in part and a block connection in part.
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
- Next to a face of type `410` (assigned state from a file), with `diffusive-metric = stencil` (the default), the metrics of
  the faces between the cells of the first row and of the source terms of the turbulence models put the state of the file at
  the centre of the ghost cell of the face, the first cell mirrored through it, where the reconstruction already took it, and
  the turbulence variables, which the file does not give and the ghost cell takes from the cell next to the face, at the
  centre of that cell. Every quantity sat at 3x₁ − 3x₂ + x₃: the derivative of the copied ones across the first row was
  half the right one on uniform cells, and that of the quantities of the file was wrong where the grid lines across the face
  are not straight and uniform, at any resolution. The ghost cells take the turbulence variables and the soot of the cell per
  unit mass: they took them per unit volume with the density of the file, so that where the two densities differ a uniform
  field of these quantities did not stay uniform across the face, and the gas entering carried them scaled by the ratio of the
  densities. Viscous runs with a `410` face change next to it; with turbulence or soot, where the density of the file differs
  from that of the cell, runs with a `410` face change in the flux of these quantities through the face with every
  `diffusive-metric`. Runs without a `410` face do not change.
- The installation guide and the table of CMake options of `docs/development/structure.md` give the defaults that a `cmake`
  configuration of a new build directory takes for the options it does not set: `USE_TECIO` and `USE_SUNDIALS` are ON (options
  of ORION and of OSLO), `USE_OPENMP`, `USE_MPI` and `USE_CANTERA` OFF, while `install.sh build` turns TecIO and SUNDIALS on only
  with `--use-tecio` and `--use-sundials`. The table gave OFF for `USE_TECIO` and `USE_SUNDIALS`, and the installation guide said that SUNDIALS is built
  if `--use-sundials` is selected. The installation guide and the README say that every configuration needs C and C++ compilers:
  they said that these compilers were needed only for the optional components (TecIO, Cantera, SUNDIALS).
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
  - Surface reactions `505`: the O2 rate is the Nagle–Strickland-Constable form, 1 + k6 p in the denominator (was up to 4.7 times too high).
  - Surface reactions `505`: a gas without one of their eight species stops at set-up; no division by zero without O2.
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
