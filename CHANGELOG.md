# Changelog

All notable changes to MOSE are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

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
