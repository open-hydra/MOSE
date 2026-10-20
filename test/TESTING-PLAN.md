# MOSE Testing Plan

Reference plan for expanding the MOSE test suite to cover **physics**, **numerics**,
and **implementation**, plus a fast tier suitable for a `git push` pre-push hook
running on GitHub's CI runner.

> Status legend: ☐ to do · ◐ in progress · ☑ done

---

## 1. State of play

- **All driving is through CTest** ([`test/CMakeLists.txt`](CMakeLists.txt)), with
  tiers expressed as `LABELS`. The legacy `test/test.sh` runner has been removed —
  cases run via `ctest` (whole suite or `-L fast`) or directly per-case via each
  case's `MOSE.sh`.
- **The pre-push hook** ([`.githooks/pre-push`](../.githooks/pre-push)) runs the fast
  tier (`ctest -L fast`); **PR CI** ([`ctest.yml`](../.github/workflows/ctest.yml))
  runs everything except the `needs-*` heavy-dependency labels and, in a second job
  (`test-mpi`), the `needs-mpi` tier.
- **CI build is minimal** ([`.github/workflows/ctest.yml`](../.github/workflows/ctest.yml)):
  `USE_MPI=OFF`, `USE_TECIO=OFF`, `USE_SUNDIALS=OFF`, `USE_CANTERA=OFF`; the second
  job builds with `USE_MPI=ON` (gfortran and the Open MPI of Ubuntu).

### Hard constraint for anything on GitHub's runner
The minimal CI build means **no finite-rate chemistry, no equilibrium (Cantera/CEA),
no binary Tecplot (TecIO)**, and MPI only in the second job. Every flame / reactive /
equilibrium case is therefore off-limits for pre-push and PR CI — those belong in a
nightly/local tier.

---

## 2. Coverage map (exercised vs. available)

| Capability | Registered options | Tested today | Gap |
|---|---|---|---|
| Equations | euler, navier-stokes | both | — |
| Time scheme | euler, RK2, RK3 | RK3 only | RK2, forward-euler untested |
| Integration vars | cons, prim, prec | cons, prim (`NumericsSmoke`); prec with HLLC-PC (`Connections`/`ConnectionFlux`, `TimeIntegration`/`LocalTimeStep`) | no accuracy check of `prim` and `prec` |
| Space recon | first-order, MUSCL | both (first-order: Noh87, LimiterDefaults) | — |
| Limiters | vanleer, vanalbada, minmod, superbee, mc | vanleer (mostly) | 4 of 5 never smoke-tested |
| Riemann solvers | 15 | 13: 11 in `RiemannSmoke` (no NaN, ρ, p > 0), HLLC+Chen and HLLC-PC in `Connections`/`ConnectionFlux` | LMRoe only in the Gresho vortex, outside CTest; MiczekRoe in no test |
| Low-Mach (Mco) | floor | HLLC+Chen and HLLC-PC in `Connections`/`ConnectionFlux` (two blocks against one, conservation); the Gresho vortex is outside CTest | no CTest check of low-Mach accuracy |
| IRS / multigrid / CFL-ramp | yes | implicit only | no dedicated check |
| Turbulence | SA(+R/RC/comp), SST, Wilcox2006 | SA, SA-RC, SST, SST-RC and Wilcox2006 in checks of sources, right-hand sides and invariants (`TimeIntegration`/`TurbulenceSources`, `MMSResiduals`, `PointImplicitSteady`, `PlanarNozzleGrids`); the turbulent flat plate is outside CTest | no CTest check against a reference solution; SA-R untested; SAcomp refused at set-up |
| Chemistry | frozen, finite-rate, equilibrium | finite-rate and equilibrium runs that must end, stop at set-up or keep the bits of a run without chemistry (`InputChecks`/`ExitCodes`; with MPI, equilibrium in `MPIVertexStencil`); the reactive cases are outside CTest | no CTest check of reactive results against a reference |
| Multi-block | core feature | faces between blocks: one block against two and records sorted against reversed, bit for bit (`Connections`/`ConnectionFlux`, `Connections`/`ConnectionOrder`; with MPI `MPIConnectionFlux`, `MPIConnectionOrder`) | 3-D faces with permuted index directions, connections with multigrid and with RANS untested |
| Ablative wall BCs (503/504/505) | melting, pyrolysis, surface reactions | 503/504 vs closed form (`AblatingWall`) | 505 needs a CO/CO2 thermo table |
| Restart | core feature | every format bit for bit, counters included (`Restart`, `Restart`/`RestartCounters`, `MPIRestart`; `tecplot binary` with TecIO) | `MPIRestart` in `vtk raw` only |
| I/O formats | tec ascii/bin, vtk ascii/raw | ascii; `vtk raw` files are well-formed XML with the right names, the field, wall and diagnostic outputs of a steady run carry minus the global iteration count as their time in `vtk raw`, `vtk ascii` and Tecplot ASCII (`OutputFiles`/`VtkNames`), and `vtk raw` files hold raw appended data (`OutputFiles`/`VtkRawFormat`) | round-trip untested |
| Parallelism | OpenMP, MPI | 1 vs 4 threads (`OpenMPEquiv`); MPI ranks × threads bit-identical in the `needs-mpi` tier (`test/mpi/`) | MPI tier runs only in MPI builds: the second CI job, `test-mpi` (gfortran, Open MPI) |
| 3D | full | source terms only | no 3D NS verification |

Cases of [`docs/vv`](../docs/vv/index.md) count as tested in this map only when CTest runs them. For every case
folder: its CTest test, which fails when the criterion fails, or outside CTest; the reference data; the criterion
of the case:

| Case (`test/`) | Page (`docs/vv/`) | CTest | Reference | Criterion |
|---|---|---|---|---|
| `1D/Sod79` | `1D.md` | `Sod79` | analytic solution (Sod 1978) | `check.py` |
| `1D/Ein91` | `1D.md` | `Einfeldt91` | analytic solution (ExactPack) | `verify.py` |
| `1D/Noh87` | `1D.md` | `Noh87` | analytic solution (ExactPack) | `verify.py` |
| `1D/Tor99` | `1D.md` | `Toro99` | analytic solution (ExactPack) | `verify.py` |
| `1D/Shu-Osher` | `1D.md` | `ShuOsher`, `ShuOsherConvergence` | MOSE solution on N = 1600 (`reference/`) | `verify.py`, `convergence.py` (F8, P5) |
| `2D/euler/prandtl-meyer` | — | `PrandtlMeyer` | analytic Prandtl–Meyer expansion | `check.py` |
| `2D/euler/isentropic-vortex` | `2D-isentropic-vortex.md` | `IsentropicVortex` | exact solution | `verify.py` (N1) |
| `2D/viscous/mms` | `2D-mms.md` | `MMSViscous` | manufactured solution | `verify.py` (N2) |
| `2D/viscous/mms/skewed` | `2D-mms-skewed.md` | `MMSViscousSkewed` | manufactured solutions | `verify.py` (N2) |
| `3D/viscous/mms` | `3D-mms.md` | `MMSViscous3D` | manufactured solutions | `verify.py` (N2) |
| `2D/viscous/planar-nozzle` | `2D-planar-nozzle-grids.md` | `PlanarNozzleGrids` | mirror symmetry of O and H grids | `verify.py` (N8) |
| `2D/viscous/ablating-wall` | `2D-ablating-wall.md` | `AblatingWall` | closed form of the wall state | `run.sh` (P1) |
| `1D/Fer14` | `1D-fer14.md` | outside | digitized velocity and temperature profiles (`reference/*.csv`) | none: `verify.py` plots |
| `1D/diff1` | `1D-diffusion.md` | outside | Cantera (`reference/reference_unity_lewis.py`); MOSE fields in `reference/` | none: no verifier in the folder |
| `1D/diff2` | `1D-multicomponent-diffusion.md` | outside | Cantera (`reference/reference_multicomponent.py`); MOSE fields in `reference/` | none: no verifier in the folder |
| `1D/premixed-flame` | `1D-premixed-flame.md` | outside | Cantera free flame (`reference/flame_phi_1.0.csv`) | none: `verify.py` prints and plots |
| `2D/euler/chang-nozzle` | — | outside | Mach on the axis and at the wall (`reference/ref1.dat`, `ref2.dat`) | none: `verify.py` plots |
| `2D/euler/gresho-vortex` | `2D-gresho-vortex.md` | outside | the initial vortex, a steady solution | none: `verify.py` prints and plots |
| `2D/euler/hypersonic-cylinder` | `2D-hypersonic-cylinder.md` | outside | normal-shock relations at the stagnation line | none: `verify.py` prints and plots |
| `2D/euler/oblique-shock` | `2D-oblique-shock.md` | outside | oblique-shock relations | `verify.py`: p2/p1, T2/T1, ρ2/ρ1 and M2 within 1.5 % for each solver |
| `2D/euler/ramp-channel` | — | outside | oblique-shock relations | none: `verify.py` plots |
| `2D/euler/rocket-nozzle/frozen` | `2D-nozzle.md` | outside | CEA (`reference/CEA.out`); `reference/field.tec` | none: `../verify.py` plots |
| `2D/euler/rocket-nozzle/finite-rate` | `2D-nozzle.md` | outside | `reference/field.tec` | `verify.py`: p, ρ, u, v and γ at the throat within a relative 1e-3 |
| `2D/euler/rocket-nozzle/equilibrium` | `2D-nozzle.md` | outside | CEA (`reference/CEA.out`); `reference/field.tec` | none: `../verify.py` plots |
| `2D/euler/supersonic-forward-step` | `2D-woodward-colella.md` | outside | OpenFOAM solution (`reference/OpenFOAM/`) | none: `verify.py` prints L2 and L∞ differences |
| `2D/viscous/flat-plate-laminar` | `2D-flat-plate-laminar.md` | outside | Blasius solution, computed by `verify.py` | `verify.py` prints PASS for an RMS c_f error below 2.5 % at x ≥ 0.05 m; its exit code is 0 either way |
| `2D/viscous/flat-plate-turbulent/SA` | `2D-flat-plate-turbulent.md` | outside | c_f of CFL3D (`reference/cf.dat`) | `verify.py`: RMS c_f error below 5 % |
| `2D/viscous/flat-plate-turbulent/SST` | `2D-flat-plate-turbulent.md` | outside | c_f of FUN3D (`reference/cf.dat`) | `verify.py`: RMS c_f error below 8 % |
| `2D/viscous/flat-plate-turbulent/Wilcox2006` | `2D-flat-plate-turbulent.md` | outside | c_f of CFL3D (`reference/cf.dat`) | `verify.py`: RMS c_f error below 8 % |
| `2D/viscous/swbli/SA` | `2D-swbli.md` | outside | Schülein experiment; SU2, Wind-US, OpenFOAM (`reference/`) | none: no verifier; its initial condition is a binary Tecplot file (TecIO) |
| `2D/viscous/swbli/SST` | `2D-swbli.md` | outside | as SA | as SA |
| `3D/coriolis-channel` | `rotating-frame.md` | outside | — | none: no verifier |
| `3D/pressure-centrifugal` | `rotating-frame.md` | outside | closed-form isentropic solution of the page | none: no verifier |
| `3D/hyshotII` | `3D-hyshotII.md` | outside | wall pressure of the experiment and of TAU with SA (`reference/`) | none: no verifier; mesh and initial condition in binary Tecplot (TecIO) |

---

## 3. Proposed new tests

### Numerics (biggest gap: order-of-accuracy verification)
- ☑ **N1 — Isentropic vortex convection** (2D, smooth, periodic): convect one period,
  grid-refine 32²/64²/128², **measure observed order** (→2 with MUSCL). *Done:*
  `test/2D/euler/isentropic-vortex` (`IsentropicVortex`, validation tier). Normalised
  L2 density error 2.56 → 0.56 → 0.17 %, observed order 2.18 then 1.75 (≈2). Meshes/ICs/BCs
  pre-generated with GRIB + `build_ic.py` + ATLAS BCB (periodic). *Possible extension:*
  a Riemann × limiter sweep.
- ☑ **N2 — Method of Manufactured Solutions (NS)**: source-term forcing, grid refine,
  verify 2nd order on the *viscous* operator. Catches gradient/metric bugs N1 cannot.
  *Done:* `test/2D/viscous/mms` (`MMSViscous`, validation tier). A smooth periodic
  field is made an exact steady NS solution by an analytic source (derived with
  sympy, FD-checked) injected by a dedicated driver `MOSE_mms`
  ([`src/app/mms.f90`](../src/app/mms.f90), a build variant of `main.f90` whose
  external callback adds the source). Constant `mu=10` (flat transport table),
  `Pr=0.72`, `Re~5` (viscous ~20% of convective), mean `M~0.47` (high enough to
  avoid low-Mach accuracy degradation). On 32/64/128 every primitive (rho,u,v,p)
  converges at ~2nd order in both L2 and Linf (rho L2 0.059 → 0.013 → 0.0028 %,
  order ~2.2). CTest runs the fast 32/64 subset (single-threaded); the full
  32/64/128 study makes the V&V figure. Periodic meshes/ICs pre-generated
  (GRIB + `build_ic.py`); the periodic `bc.txt` (topology-only) is shared with N1.
  *On general grids:* `test/2D/viscous/mms/skewed` (`MMSViscousSkewed`) and
  `test/3D/viscous/mms` (`MMSViscous3D`), validation tier: the fields periodic2d,
  wall2d and periodic3d of the driver (`MOSE_MMS_FIELD`) on two-block grid families
  whose blocks meet with permuted indices, zigzag (every cell skewed, its
  neighbours the other way) and iface (a cell-size jump of 1.5 and kinked grid
  lines). CTest runs the short forms (2-D on N = 10, 20; 3-D zigzag on N = 8, 12):
  L2 errors decreasing on the finest pair, an observed L2 order of at least 0.7 on
  zigzag and 1.4 on iface (70 % of the order each family should show), with walls
  the same order for the heat flux integrated over the wall, and the totals of the
  equations with zero boundary flux within the smallest L2 error of each grid.
- ☑ **N3 — Riemann-solver smoke matrix**: coarse Sod, few steps, once per solver;
  assert positivity (ρ,p>0) and no NaN. *Done:* `test/fast/riemann-smoke` (11
  shock-capturing solvers; low-Mach family deferred to a Gresho smoke).
- ☑ **N4 — Limiter / scheme sweep**: all 5 limiters + 3 time schemes + recon +
  integration-variables; assert positivity. *Done:* `test/fast/numerics-smoke`.
  (Monotonicity assertion is a possible future tightening.)
- ☑ **N8 — Planar nozzle on grid families (RANS)**: a turbulent planar
  converging-diverging nozzle, Spalart–Allmaras and SST, on grid families that
  differ only in the direction of their grid lines (O, H, ZZ, SB2). *Done:*
  `test/2D/viscous/planar-nozzle` (`PlanarNozzleGrids`, validation tier). The short
  form runs 300 iterations on the coarse level from a quasi-1D field: on O and H,
  mirror images of themselves whose four blocks have different index
  orientations, the largest mirror asymmetry of every variable below 1e-10 of its
  largest magnitude; on ZZ (Spalart–Allmaras) a finite flow with positive density
  and pressure. The full form runs every family on three levels to convergence.

### Physics
- ☑ **Turbulent flat plate, SA / SST / Wilcox2006** — one case per model under
  `test/2D/viscous/flat-plate-turbulent/{SA,SST,Wilcox2006}/`, each vs NASA CFL3D/FUN3D
  reference Cf, outside CTest (200 000 iterations with SA, 100 000 with SST and
  Wilcox 2006; see §2). *No further turbulent flat-plate cases needed.*
- ☐ **P3 — Couette / Poiseuille** (2D viscous, analytic): clean fast laminar-NS check
  independent of Blasius subtleties.
- ☐ **P4 — 3D Euler verification**: 3D isentropic vortex or manufactured case (3D is
  currently source-term-only).
- ☑ **P5 — Shu–Osher** (1D): Mach-3 shock / sine interaction. *Done:*
  `test/1D/Shu-Osher` (hand-generated IC via `gen_ic.py`). Two CTest entries:
  `ShuOsher` (**fast**, single N=200 grid vs an N=1600 MOSE reference, 2.9 % L1) and
  `ShuOsherConvergence` (**validation**, 4-grid Richardson study N=200/400/800/1600 —
  the N=1600 solution is a study point — asserting a sane observed order (p ≈ 1.02)
  and monotone L1 decrease vs the Richardson extrapolant, 3.55 → 2.83 → 1.46 → 0.72 %,
  and producing the V&V figure).

### Implementation / infrastructure (also serve as fast tests)
- ☐ **I1 — Freestream preservation / GCL**: uniform flow on a stretched, curvilinear,
  multi-block mesh; assert residual stays at machine-zero.
- ◐ **I2 — Single-block vs multi-block equivalence**: same problem as 1 block and 2×2
  blocks; fields identical to round-off. Direct halo/interface-flux test. *Done for two
  blocks:* `test/fast/connections/connection-flux` (`Connections`/`ConnectionFlux`, a periodic shock tube as one and
  as two blocks, bit for bit, and conservation) and `test/fast/connections/connection-order`
  (`Connections`/`ConnectionOrder`, two 2-D blocks with permuted index directions, the `bc.txt` records
  sorted and reversed); with one block per MPI rank, `test/mpi/connection-flux`
  (`MPIConnectionFlux`) and `test/mpi/connection-order` (`MPIConnectionOrder`), needs-mpi
  tier. Four blocks at a vertex and 3-D faces with permuted directions are not covered.
- ☑ **I3 — Restart round-trip**: N steps vs (N/2 → write → restart → N/2); match to
  ~2e-15. *Done:* `test/fast/restart` (restart = `newrun=false`, reads OUTPUT/field.tec).
  With `sol-format = vtk raw` the same test asks for identical bits, iteration counters
  and output numbering included, and a restart from a solution without a time (Tecplot
  ASCII without SOLUTIONTIME, vtk raw without TIME) must stop at setup, and a solution of
  MOSE copied into INPUT/ must start a new run (vtk raw, vtk ascii, Tecplot ASCII); *MPI variant (needs-mpi tier):* `test/mpi/restart`
  (`MPIRestart`, 2 ranks × 2 threads). `test/fast/restart/restart-counters` (`Restart`/`RestartCounters`)
  covers the other counters, in every format (vtk raw, vtk ascii, Tecplot ASCII, Tecplot
  binary when MOSE is built with TecIO): steady runs (global iteration count, probes, diagnostic
  output), time-based output, the multigrid level, the global iteration count of a
  multigrid restart from Tecplot ASCII, a transient started from a steady solution (vtk raw
  and Tecplot ASCII), a steady restart from a vtk solution that holds TIME only, and
  the residual history and probe rows of a restart from a solution older than the end
  of the run that wrote it, also from history files cut by hand without a final newline,
  and the rows of `torque.dat` of the rotating frame in the same restart. A first part that
  ends on an output iteration (by iteration, also with `sol-overwrite = true` and with two
  grid levels, and by time) writes there the solution and the probe row due, and the restart
  in the same folder writes the files of the uninterrupted run. A first part killed with
  SIGKILL at its third numbered solution, as a queue kills a run at its time limit, leaves in
  `residual-history.dat` and in the probe files every row it wrote, and the restart from its
  last solution writes the files of the uninterrupted run.
- ☑ **I4 — Conservation invariant**: closed slip-wall box, shock-tube IC; total mass
  conserved to round-off (~2e-13). *Done:* `test/fast/conservation`. (Energy invariant
  is a possible future addition.)
- ☑ **I5 — Serial vs OpenMP equivalence**: 1 vs 4 threads bit-identical, also with chimera faces meeting at
  the edges of a block. *Done:*
  `test/fast/openmp-equiv`. *MPI variant (needs-mpi tier):* `test/mpi/boundary-records`
  (`MPIBoundaryRecords`) — a strip of four viscous blocks bit-identical on 1×1, 2×1,
  4×1 and 2×2 ranks × threads, and with the `bc.txt` records in reverse order (1 and 2
  ranks, with and without threads), the wall output included; the same with the
  Spalart–Allmaras model, ordered and reversed (1×1, 1×4, 2×2), and with SA-RC, whose
  velocity gradients cross the faces between ranks (2×1, 4×1, 2×2 against 1×1);
  `test/mpi/vertex-stencil` (`MPIVertexStencil`) — a
  2×2 box of viscous blocks with a vortex on the common vertex, bit-identical on 1×1,
  1×4, 2×1, 4×1 and 2×2 ranks × threads; with chemistry = equilibrium, bit-identical
  on 1×1 and 2×1; `test/mpi/bit-identity` (`MPIBitIdentity`) — every file a run writes
  in `OUTPUT/` (solutions, residual history, probe) bit-identical on 1×1, 1×4, 2×1, 4×1
  and 2×2 ranks × threads for five cases: a 2×2-block box with inlets, outlets, walls
  301 and 302 and symmetry planes, Spalart–Allmaras and SST steady (local time step,
  implicit residual smoothing) and Spalart–Allmaras time-accurate, and the points
  shared by 3 and by 5 blocks of `ViscousPatch`/`ViscousLinearPatch` with a field that changes at every
  stage. Its 25 solver runs take ~40 s of CPU time (67 s of wall time, load average ~178
  on 192 logical CPUs); main has no needs-mpi tier to compare with.
- ☑ **I6 — MPI thread level**: MOSE requests `MPI_THREAD_FUNNELED`, so only the thread
  that initialised MPI may call it. *Done (needs-mpi tier):* `test/mpi/thread-funnel`
  (`MPIThreadFunnel`) — a probe library preloaded into the solver counts no MPI call
  from other threads on 2×2, 2×4 and 2×2-with-IRS (steady) ranks × threads.
- ☑ **I8 — Block-to-rank assignment**: the optional `INPUT/decomposition.map`, in the
  format of the ATLAS mesh decomposer (MDB), gives the MPI rank of every block. *Done:*
  `test/fast/input-checks/decomposition-map` (`InputChecks`/`DecompositionMap`) — on one rank a map with every block
  on rank 0 gives the bits of the run without it, also the map of the halo-trade test
  of MDB kept byte for byte, and maps that do not fit the mesh or the run (a number
  written as a real and a directory included) stop it at setup with the fault named;
  *MPI variant (needs-mpi tier):* `test/mpi/decomposition-map` (`MPIDecompositionMap`)
  — on 2, 3 and 4 ranks every map (MOSE's own assignment, permutations, an idle rank
  0, 110 blocks, fewer blocks than ranks) gives the one-rank field, and a map for other
  ranks stops every rank with the same message.

---

## 4. Pre-push fast suite (GitHub-runner friendly)

Rules: minimal build only, tiny meshes (≤ ~1k cells), ≤ a few seconds each, assert
**hard invariants** (no digitized references, no plotting, no Cantera). Target **< 60 s** total (not met: see the count below).

| # | Test | Asserts | Status | ~Cost |
|---|---|---|---|---|
| F1 | Freestream preservation (stretched multi-block) | residual ≈ machine-zero | ☐ | <2 s |
| F2 | Sod coarse (existing) | shock/contact positions, positivity | ☑ `Sod79` | ~1 s |
| F3 | Conservation in closed box | mass conserved to round-off | ☑ `Conservation` | <1 s |
| F4 | Single- vs multi-block equivalence | identical field | ☑ `Connections`/`ConnectionFlux`, `Connections`/`ConnectionOrder` (I7, I7b) | see I7, I7b |
| F5 | Riemann smoke matrix; shock-detector = Tramel with HLLC and HLLC+Tramel | no NaN, ρ,p>0; the set-up summary names the shock sensor, and HLLC+Tramel gives the field.tec of the run without the key | ☑ `RiemannSmoke` | ~3 s |
| F5b | Numerics sweep (limiters/schemes/recon) | no NaN, ρ,p>0 | ☑ `NumericsSmoke` | ~5 s |
| F5c | Limiter defaults (MUSCL without limiter, first-order with limiter) | field bit-identical to the explicit setting; matching warnings and summary | ☑ `InputChecks`/`LimiterDefaults` | <1 s |
| F5d | Residual smoothing with time-accurate runs | refused at set-up, no solution; steady runs keep it | ☑ `InputChecks`/`IRSTimeAccurate` | <1 s |
| F5e | Time control (1-D duct, time-dependent inlet 402; 1-D tube of air at Mach 2 with an entropy wave, steady, `res-threshold` = 1e-6, at 1 bar and at 1 bar × 2⁻²⁴, and at 1 bar in two parts with a restart from `vtk raw` after 60 iterations) | order in time ≥ 1.8 (RK2, RK3); table held before its first instant; time-threshold reached; the steady runs end at the same iteration, with the density residual at most 1e-6 of that of the first iteration, the run in two parts as the run in one go | ☑ `TimeIntegration`/`TimeControl` | ~5 s |
| F6 | Restart round-trip; probes at a restart, the file of one of them removed | identical to uninterrupted run; the missing probe file written again, the other appended, with the probe rows of the uninterrupted run | ☑ `Restart`, `Restart`/`RestartCounters` | ~3 s, ~6 s |
| I11 | Wall output in vtk raw (2x2-block box of 8 × 8 cells, adiabatic no-slip walls (301) on its four sides, 5 RK3 steps) | the run ends; every wall file `wallB<b>F<f>.vts` has the extent of its face (WholeExtent and Piece Extent, from `fieldBlock<b>.vts`), its points and one value of each wall variable per cell, as raw appended data | ☑ `OutputFiles`/`WallOutput` | <1 s |
| I5 | Serial vs OpenMP equivalence: the Sod tube; a chimera case, a block overset on another whose four chimera faces meet at its edges (Navier–Stokes; it passes on the parent of its commit too: the race it guards shows under ThreadSanitizer with Archer, in three runs of three) | 1 vs 4 threads bit-identical | ☑ `OpenMPEquiv` | ~3 s |
| I12 | Names of the VTK outputs (steady viscous run on a 2 × 2-block box with no-slip walls, 5 iterations, `vtk raw`, laminar and with Spalart–Allmaras and `sol-variables = thermo transport`; the laminar run also in `vtk ascii` and in `tecplot ascii`) | every `.vts` that `field.vtm`, `wall.vtm` and `diagnostic.vtm` list parses as XML, and its data arrays are named as the variables that MOSE writes for the case (field, walls, diagnostic; the turbulence variable, the laminar viscosity and conductivity and the eddy viscosity with Spalart–Allmaras), with no other array; every `.vtm` lists its block files as `vtk/<file>.vts`, and every data array has the format of `sol-format`; the outputs of the steady run carry minus the iteration count as their time (`TIME`, `SOLUTIONTIME` without `STRANDID`) | ☑ `OutputFiles`/`VtkNames` | <1 s |
| I13 | `vtk raw` solutions as raw appended data (the Sod tube of `test/1D/Sod79`, 10 iterations, a numbered solution every 3, steady with the diagnostic output after each solution, and time-accurate) | every data array of every block file of every `.vtm` in `OUTPUT/` is `format="appended"`, with the raw appended data in the file | ☑ `OutputFiles`/`VtkRawFormat` | <1 s |
| F9 | Viscous stencil of chimera faces: two blocks of 8 x 6 cells (8 x 6 x 4 in 3-D) overlapping by 4 cells in x, each ghost cell on a cell of the other block; one step of a steady run with u = G y (2-D), u = G (y + z) (3-D), and a radial velocity on a 2-D wedge of 1 degree against the same wedge as one block; a control: the 2-D case with the y spacing of block 1 growing 4 times across its last cell, which turns over the first layer of extrapolated ghost cells behind the chimera face | y and z momentum residual zero within 1e-9 of μ G / H, also in the cells at the ends of the chimera faces (1/4 of it with the ghost cell copied beyond the edge); on the wedge the momentum residual of every cell that of the one block within 1e-9 of μ G / H (half of v/r at the chimera faces with the ghost cell as its own k neighbour); the control stops at set-up on the ghost cell behind the chimera face, whose viscous flux takes its metric | ☑ `Connections`/`ChimeraFace` | ~1 s |
| F7 | Symmetry (coarse forward-step/wedge) | top/bottom mirror preserved | ☐ | <3 s |
| F8 | Shu–Osher coarse (N=200) | L1 vs N=1600 ref < 5 %, finite | ☑ `ShuOsher` | ~2 s |
| P1 | Ablative wall BCs (503/504); the melting column also without `schmidt` (default 0: one species needs no `diffusion.dat`) | wall state vs closed form, both face orientations, inert-surface fallback | ☑ `AblatingWall` | ~12 s |
| I9 | Exit codes of the error stops: the Sod tube as a viscous mixture of four species without `schmidt` and without `diffusion.dat` ([thermo](../docs/theory/thermo.md)); the Sod tube in air with a probe of an unknown variable; a steady run with the soot model LL91, with and without Spalart–Allmaras (six species with C(gr)), and the same state with `turbulence = SAcomp`; the Sod tube in air, steady, with sol-variables = none; the Sod tube that diverges (CFL 20), with an initial condition of too few variables, with a rotating frame of zero axis, of NaN axis or of infinite origin, with chemistry = equilibrium and no `INPUT/composition.txt`, with a value outside the rule of the registry: an interval, a list of reals, a list of integers; the Sod tube with a NaN velocity in one cell and with an infinite coordinate of one node; the duct of BCInput with a NaN and an infinite pressure, and a NaN time, in the table of its inlet 402; the duct with a NaN and an infinite total pressure and a NaN alpha of an inlet 401, an infinite total temperature of its inlet 402 and a NaN roughness of a wall 302 in its bc.txt; the duct with its inlet a face 410 whose file holds a NaN pressure and an infinite SOLUTIONTIME, and the same file missing, cut in the values of its second zone, with a word among them, with I=two in the header of its second zone, with I = 3 or j=2 there (no I= or no J=), with one value more in its first zone and with its two zones in decreasing time; the Sod tube restarted from its solution with a NaN SOLUTIONTIME and with an infinite MOSE_LAST_OUTPUT_TIME; input.ini missing, empty, of comment lines only; the Sod tube with sol-dtime = NaN, Infinity and 0,5, sol-diter = seven and 2.5, sol-overwrite = yes and levels = 2.0 of [MOSE-Multigrid]; a valid run with finite-rate chemistry and exclude-blocks = 1; exclude-blocks with a word that is not a block number or with a block beyond the mesh, which stops at set-up; the tube filled with steam far from chemical equilibrium, its block excluded from chemistry = equilibrium, on one grid and on two multigrid levels; on 2 MPI ranks (`MPIExitCodes`, needs-mpi tier), cfl = -1, without input.ini, the NaN velocity with the block on rank 1 alone and the infinite coordinate, and cfl = -1 written during the run in the input.ini of rank 0 (two folders, MPMD), and the initial condition of rank 1 alone cut in the middle (two folders) | non-zero exit code and the diffusion-table error, no solution; exit status 1 and the name of the variable, no solution; exit status 0 and a diagnostic output that names every variable it writes, the local time step under dt; with SAcomp, exit code 1 to 128 at set-up, the message and no diagnostic output; exit status 1 and a message that names thermo or transport; exit code 1 to 128 and the error message (for input.ini, the reason it cannot be read, and no solution); for a key of input.ini whose value does not read as its type, exit code 1, the file, the section, the key and the value as written and the STOP line, no solution; for a value that is not finite, exit code 1, the place of the value and the STOP line, no solution; for the table, exit code 1, the line and the column, no solution; for bc.txt, exit code 1, the line, the record, its type and the name of the value, no solution; for the file of the face 410, exit code 1, the file, the zone, the line of its header, the variable and the cell, no solution, and for a file that cannot be read, exit code 1, the file and, for a zone, the zone and the line of its header (with the two times for zones out of time order, and the form of I and J to write for a header without I= or J=), no solution; for the restart, exit code 1, the file and the name of the value, the STOP line, the solution left as it was; exit code 0 and block 1 under the inactive blocks of the set-up summary; the field of the excluded block bit-identical to the run without chemistry, and different with the block not excluded; on 2 ranks, with cfl = -1 exit code 1 within 60 s, the message and the STOP line of both ranks, without input.ini exit code 1 to 128 within 60 s and the message, and for a value that is not finite exit code 1, the place written once and the STOP line of both ranks; with the input.ini of rank 0, which rank 0 passes to every rank, exit code 1 at iteration 20 and the rule, not a hang; with the initial condition of rank 1 cut, exit code 1, Initial conditions FAIL, the error of the read written once and the STOP line of both ranks | ☑ `InputChecks`/`ExitCodes` | ~12 s |
| F10 | Quasi-2D with symmetry k faces; 1-D with symmetry j and k faces; a channel of variable section in a layer between parallel planes inclined to z and a 1-D channel of variable section, steady subsonic flow | p, u, v uniform; round-off in w does not grow (Heun, RK3); in 1-D, round-off in v and w does not grow; in the layer and in the channel the velocity across is integrated, with the warning, and p within 0.03 p0 of the quasi-1-D flow | ☑ `TimeIntegration`/`Quasi2DSymmetry` | ~8 s |
| F11 | Local time step (square and skewed 2-D cells, cubic 3-D cells, heat-conduction and species-diffusion cells, 1-D at cfl 0.99; heat conduction on a periodic block of 8³ rhombohedra at 115°, Heun at vnn 0.5, 600 steps); with integration-variables = prec, a 2-D column in simple shear over a wall with SST and Wilcox 2006 | stable, uniform state kept, no NaN; the temperature wave of the rhombohedra decays; with prec the step of every cell (but the row next to the symmetry plane) equal to the rule with U_r from the viscous floor with the eddy viscosity of the model at the strain rate of the cell, within 1e-10, after one and after two steps | ☑ `TimeIntegration`/`LocalTimeStep` | ~7 s |
| F12 | Turbulence sources in time-accurate runs (wall channel, SA and SST, Heun and RK3; SST-RC, RK3, one step); Wilcox 2006 on a 3-D linear flow, divergence-free, with S^ zero on the diagonal (cube of 6³ cells), one step; Wilcox 2006 and SST with omega-wall-bc = asymptotic on the channel at rest, omega at the wall value 6 nu/(beta d²) of the model, one step | time order of nu~, k, omega ≥ 1.8 (Heun) and ≥ 2.7 (RK3); no negative stage value, and a first step equal to the rule with the bound of the stiffness that the state of each cell gives alone within 1e-6; the start-up warnings of k-omega and of -RC both printed (SST-RC); change of rho omega = dt (2 gamma rho S:S - beta_0 f_beta rho omega²) within 1e-6; in the wall cells, change of rho omega = -dt beta rho omega² (the destruction alone) within 1e-6 | ☑ `TimeIntegration`/`TurbulenceSources` | ~6 s |
| F13 | Spalart–Allmaras on constructed fields: variable density (periodic field, no wall); SA-RC and SST-RC on a 3-D flow with vorticity along x, y and z (cube of 8³ cells) and on the same flow turned (x to y, y to z, z to x), one step each; SA-RC on a uniform shear in a frame turning about z and on a uniform flow (8 × 8 cells), one step each; SA-RC and SA on air at rest next to a wall and on air flowing along it on cells stretched along the flow (8 × 8 cells), one step each; the viscous stress on the axisymmetric wedge (200): Homann's stagnation-point flow and a polynomial field with dilatation and swirl (1° wedge, axis included), and the axis at r = 0 on the wedge and on a 3-D sector of 4° (200), where viscous runs stop at setup, as on the polar line of a 3-D mesh with inverse-mean and mean; Spalart–Allmaras next to a wall (fluid at rest, nu~ = a y, chi < 1e-3), also on cells growing by 1.2 from the wall and on skewed cells (2-D, 3-D) with the density linear along the wall; SA-RC on the column with a shear flow whose strain rate varies along x, as one block, two blocks and two blocks with permuted indices | right-hand side of rho nu~ = standard conservation form, error < 0.02 and order > 1.5; change of rho nu~ (rho k and rho omega with SST-RC) the same in corresponding cells of the two 3-D runs within 1e-9 of its largest value; change of rho nu~ = dt rho c_b1 f_r1 S~ nu~ with the f_r1 of the specification within 1e-6 (shear), no NaN and no change (uniform); the same change of rho nu~ with SA-RC and SA within 1e-12 (rest; stream, off the wall row); viscous force vs exact: Homann, radial and axial error < 0.02 and order > 1.5, azimuthal zero; polynomial, radial and azimuthal error < 0.02 and order > 1.5, axial exact; right-hand side of rho nu~ below 0.05 of (cb2/sigma) rho a^2 in every cell, wall cells included (grad nu~ = a there); on the skewed cells the terms in grad rho (difference with the run at uniform temperature) below 0.05 of the same scale; SA-RC: two blocks, and two blocks with permuted indices, give the change of rho nu~ over the step of one block to 1e-9 of its largest value (a ghost cell that copies the cell: about 3e-2) | ☑ `MMSResiduals` | ~7 s |
| I7 | Faces between blocks (101, 201): periodic shock tube as one and as two blocks, HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC (low Mach) | two blocks = one block bit for bit; mass, momentum, energy conserved to 1e-13 | ☑ `Connections`/`ConnectionFlux` | ~2 s |
| I7b | Faces between blocks of the same parity (2-D, index directions permuted across the cut), HLLC+Tramel, Navier–Stokes; the same with the block below the cut one cell thick across it, and two cells thick with two multigrid levels | bc.txt records sorted = reversed, bit for bit; the thin blocks stop the run at set-up with a non-zero exit code and the message (the record, the block across, its cells, the level) | ☑ `Connections`/`ConnectionOrder` | <1 s |
| F15 | Wall output and residual (converged laminar hypersonic plate) | residual of ρu and ρe in the wall cells below 1e-6 of the wall shear force and heat flow of the face | ☑ `OutputFiles`/`WallOutputResidual` | ~2 s |
| I8 | Block-to-rank map `INPUT/decomposition.map` on one rank, also a map written by MDB | rank-0 map = no map bit for bit; misfit maps stop at setup with the fault named | ☑ `InputChecks`/`DecompositionMap` | ~6 s |
| I10 | input.ini read again during a run (2x2-block viscous box, 48 × 48 cells per block: ini-diter = 20, 40 iterations; an input.ini with the required keys only, and one that also writes integration-variables = prim, riemann-solver = HLLE and flux-limiter = minmod, ini-diter = 5, 10 iterations); a comment line with "=" before cfl in [MOSE-Numerics], cfl = 0.4, 3 iterations, and the same run without the line; the same box with 24 × 24 cells per block in a rigid-body frame about x with the wall of x = 0 at rest (stationary-face-1 = 1 1), 30 iterations, solutions every 7, omega, rotational-frame, the stationary face or an axis of NaN written during the run, and, with two stationary faces at the set-up (stationary-face-2 = 2 2), omega written with both keys, with the first alone and with none; on 2 MPI ranks (`MPIRuntimeInput`, needs-mpi tier), the 2x2-block box in two folders (MPMD), 40 iterations, sol-diter = 7 written during the run in the input.ini of rank 0 alone, then sol-diter = 7, res-diter = 7, ini-diter = 30 and cfl = -1 in the input.ini of rank 1 alone | cfl = -1 stops the run at setup and, written during the run, at the reading of iteration 20 with the message of the rule; with the required keys only, and with the three keys written, the run reaches iteration 10 with the residual history of the same run without readings, byte for byte; with the comment line the run writes the residual history of the run without it, byte for byte, and the setup names no key as not read; omega from 1e4 to 2e4 rad/s: from iteration 21 the run of the restart from iteration 20 with 2e4 at set-up (residual history and solution of iteration 30 byte for byte); rotational-frame = none, stationary-face-1 = 1 3, axis = NaN 0.0 0.0: stop at the reading of iteration 20 with the message; the keys left out keep their faces: with the first key alone and with none, iteration 30 and every output file of the run with both keys, byte for byte; on 2 ranks every rank takes the input.ini of rank 0, which rank 0 passes to every rank: with the change of rank 0 exit code 0 within 60 s at iteration 40 and the solutions of iterations 10, 20, 21, 28 and 35; with the changes of rank 1 exit code 0 within 60 s at iteration 40 and the files of the run without the change, byte for byte, not a hang | ☑ `InputChecks`/`RuntimeInput` | ~20 s (`MPIRuntimeInput` 30-63 s) |
| F17 | Inlets: the ambient inlet (407) with the flow direction `normal` (planar duct, air entering at 50 m/s); inlets at total conditions (1-D tube of air, one forward Euler step): 407 and 401 from rest at p0 = 2 p and 36 p, an inlet the reservoir cannot feed, a 407 that it cannot feed with the gas moving in at 0.5 a and 1.2 a, and with p_back above p_stop; 401 and 407 with the gas moving in at 2a/(γ−1) or faster (6 a at p0 = p, 5.5 a from a reservoir at T0 = 25 T and p0 = 0.1 p, 20 a), and the hot inlet with the invariant of the cell at ±1e-9 a | same field as the equivalent angles bit for bit, no transverse velocity (|v| <= 64 eps (|u| + c), rounding); mass flux of each inlet face = that of the face state (subsonic, choked at normal Mach 1, slip wall, outlet at min(p_back, p_stop)) within 1e-9 of rho a; the two fluxes at ±1e-9 a agree within 1e-9 of rho a | ☑ `BoundaryConditions`/`Inlets` | ~2 s |
| F18 | Wall distance (fluid at rest, nu~ uniform; convex corner of two walls, curved wall under a leaning mesh, nozzle lip on a 1° wedge; channel and pipe on a wedge periodic in x, 3-D sector closed by rotational periodicity, the channel with its periodic sides matched to 1e-8, 3-D slab periodic in k by type 200, each with the wall on part of the period; 3-D channel periodic in x and z with periods a hundredth of its height, wall on the whole period; the 2-D channel ten periods from the origin; a 2-D staggered row of plane blades 1.7 periods long; a 3-D block periodic along three nearly dependent translations; the 2-D channel turned by 30° with its coordinates rounded to 5, 4 and 2 significant digits; 3-D walls of twisted faces: a block under leaning cells, where in 9 cells the nearest face by the four triangles of each face is not the nearest face, and a block periodic in x) | distance read back from the Spalart–Allmaras source = exact distance to the wall and its periodic images within 1e-9 (1e-6 with the sides matched to 1e-8; for the turned channel, between the distances with the translation of each periodic face; on the twisted walls, the distance to the bilinear faces from the stationarity equations in rational arithmetic); with 2 digits, stop at set-up with the message on the periodic translations; no warning on twisted wall faces | ☑ `MeshGeometry`/`WallDistance` | ~5 s |
| F19 | Mesh type (plane slab with y falling along j, right-handed wedge with the axis on its last j face, a wedge of 5 degrees in two blocks with coordinates of 5 significant digits, a plane slab far from the x axis (y from 1 to 1.001, 1e-4 thick), a wedge of 1 degree with coordinates of 2 significant digits, a plane slab with its k planes tilted by 2 % of its thickness, 3-D sectors of 90 degrees with the periodic plane turned by half a cell more (100 cells in k) or a tenth (20 cells) but at two nodes, a wedge block with a plane block, a plane slab whose k = 1 plane touches the k = 0 plane along x = 0, 3-D slab periodic in k with type 200, 3-D rings of 2 × 2 × 36 cells closed by type 200 over the full turn, with the last k plane written as the first or computed with cos and sin, and the second with a wall on its inner face) and folded cells (a 2-D and a 3-D block with one interior node moved across the next cell; one cell of a real 3-D mesh) | plane and meridional wall distance (to the wall of the nodes as written) within 1e-9, and the type in the set-up summary; set-up stop with a non-zero exit code (2 digits, tilted slab, touching k planes, shifted sectors, mixed blocks, folded cells); uniform flow kept within 1e-12, and on the rings, named 3-D with no sector angle; on the ring with a wall the distance to the plane wall face of the sector of each cell within 1e-9 | ☑ `MeshGeometry`/`MeshType` | ~3 s |
| F21 | Rotational periodicity of a 3-D sector (annulus sector from -45 to 45 degrees, type 200 on the k faces, radial flow with uniform density and pressure, 20 steps of RK3) | the columns k and NK + 1 - k mirror images (rho, u, v, p equal, w opposite) within 1e-10 of their scales | ☑ `Connections`/`SectorRadial` | <1 s |
| N7 | Patch test of the viscous fluxes: linear velocity and linear temperature on 2-D meshes with nodes moved at random, a fold between two blocks, a fan of index lines at a wall (the edges along it growing 7/4 times across the first cell: the second layer of extrapolated ghost cells turned over; 4 times: both layers, which no flux there reads), points shared by 3 and 5 blocks, a symmetry plane under skewed cells, a 3-D pair of blocks; walls with a prescribed heat flux; a control on parallelogram cells (walls 301 with q = 0) | viscous force and heat conduction of every cell below 1e-11 of their scale, also at every stage of an RK3 step around the points shared by 3 and 5 blocks; inverse-mean and mean equal to a replica of their formulas within 1e-10 | ☑ `ViscousPatch`/`ViscousLinearPatch` | ~10 s (44 one-step runs) |
| N7b | The same next to walls of prescribed temperature (302) and next to walls coupled to another solver (103, 104: T linear along the wall, a velocity linear and zero on it, the stress on the wall left to the other solver): the wall temperature varying along the wall on skewed cells, also behind a wall whose two layers of extrapolated ghost cells are turned over (the edges along it growing 4 times across the first cell); a velocity linear and zero on a no-slip wall whose cells have nodes moved at random (2-D and 3-D), and the same on a wall that does not turn with a rotating frame, its velocity in the frame -Omega x r varying along it, with every face metric (inverse-mean and mean on parallelogram cells, at a wall 302 and at an adiabatic wall 301), and there the energy of every cell; a control on parallelogram cells (T_w uniform along the wall) | heat conduction and viscous force of the wall cells below 1e-11 of their scale; on the parallelogram cells of the frame the energy right-hand side of every cell equal to tau : grad u within 1e-11 of its scale (the work of the wall stress on the wall velocity included); next to the 302 wall also within the heat flux of 6 units in the last place of T_w, k ulp(T_w) A_w/(d V) each (1.9e-11 of the scale): the temperature of the cell, p/(rho R), holds 8 roundings from the T_w of the case, 11 with fast math (11 u T_w = 6.4 ulp at 300 K), so in every declared build (gfortran and ifx, DEBUG, TESTING, RELEASE) | ☑ `ViscousPatch`/`ViscousWallPatch` | ~14.5 s (63 one-step runs) |
| F24 | Boundary input (1-D duct): inlet 402 table under four file names (plain, with a directory, longer than 32 characters, quoted with a blank), a missing table, a table whose times decrease and one that starts after the run; property lines with an empty value, short of a value, with a word for a number or with a misspelt flow direction; an SRM grain boundary (502) with the Spalart-Allmaras model and with SST | the four runs bit-identical; the missing and the decreasing table stop the run at set-up with the name of the file; the late table gives the run of the same table with its first value at t = 0, bit for bit; the bad lines stop the run with their line, record and type; the 502 with Spalart-Allmaras runs 40 steps, and with SST (k and omega never below zero), and one step at rest from a uniform rho nu~ m = mu changes it next to the grain by (m/(sigma rho h^2)) (cb2 m - 2 mu) dt plus the source terms at the distance h/2 to 1e-6 of 2 mu (the wall value nu~ = 0 in the ghost cells; an extrapolated ghost cell gives -2 mu), and in the other cells by the source terms at the distance (i - 1/2) h to the grain within 1e-9 (the grain a wall for the wall distance) | ☑ `BoundaryConditions`/`BCInput` | ~2 s |
| F23 | Rotating frame (viscous Sod tube with walls 301 at its ends, rigid-body frame) with two `stationary-face-<n>` entries, the two ends; then one on a face of type 400 (no wall), and one on a block the mesh does not have; a stationary wall across the axis, gas at rest in the absolute frame; the same with two cells across, the face a wall 302 in part and a block connection (101) in part, and with a symmetry plane (300) in place of the connection | the run ends within 60 s with exit code 0; the set-up summary lists the two stationary patches; the two other faces stop the run at set-up with a non-zero exit code and the message; wall shear in the wall output below 1e-3 of mu Omega y_c / (dx/2); the face with the connection runs (exit code 0, patch listed, wall shear below 1e-3 of the scale), the one with the symmetry plane stops at set-up with the message | ☑ `BoundaryConditions`/`RotatingFrame` | ~1 s |
| F25 | Probes given by position (planar 2-D duct, a different pressure in every cell, probes at cell centres in the mid-plane) | each probe holds the values of its cell, never a neighbour or a ghost cell | ☑ `OutputFiles`/`ProbePosition` | <1 s |
| F27 | Viscous flux of the rotational periodicity of a 3-D sector (sector of 90 degrees against the sector of 180 degrees that holds it, viscous air at 10 Pa, a velocity periodic for a turn by 90 degrees with normal stress and shear on the periodic faces; a radial flow with the k faces 200 and 300; a block periodic in z closed by 200 with no angle, a translated sector, and by periodic pairs 201, with Spalart–Allmaras; one explicit step) | viscous force of every cell equal to that of the same cell in the larger sector within 1e-9 of its scale (momentum and energy); 200 against 300 within 0.02 of the normal stress away from the j faces; translated 200 against 201: the step of rho, momentum, energy and rho nu~ of every cell within 1e-9 of the largest value | ☑ `Connections`/`SectorViscous` | <1 s |
| F26 | Sign of the wall heat flux and of the wall stress (the same wall on the two faces of one index direction, mirror images: a planar 2-D channel with the walls on the j faces, a 3-D channel with the walls on the i faces and on the k faces; one step) | `301` at q_wall = 1e4 W/m², air at rest: the cells next to the two walls warm up alike within 1e-9, qw = q_wall on both walls; `302` at 400 K: qw the same on the two walls within 1e-9 and positive; `302` with the air sliding along x at 10 m/s: tauX the same on the two walls within 1e-9 and positive, tauY opposite and tauZ the same within 1e-9 of tauX | ☑ `BoundaryConditions`/`WallHeatSign` | <1 s |
| P6 | Quadratic constitutive relation QCR2000, a unit check (QCR2000 is not selectable from the input): the stress of a 3-D velocity gradient with the eddy viscosity alone and with the laminar viscosity alone, called as the solver calls it | tau_QCR · n of the definition and the laminar stress, within 1e-12 | ☑ `QCR2000Stress` | <1 s |
| N9 | Volume of a cell, a unit check of Compute_Metric_Tensor: ten cells (cubes with one twisted face, T_n/h = 1 and 5; a valid non-convex cell with six twisted faces; a collapsed edge with two twisted faces; a 2-D prism; wedge cells of 5°, one on the axis; the first cube and the prism moved by 1e3; the first cube mirrored, left-handed), a twisted mesh of 8³ cells, an inverted cell with det A > 0 at its centre and V < 0 | the integral of det J of the trilinear cell (exact, in rational arithmetic from the binary64 nodes) within 32 u; the volumes of the mesh sum to the volume enclosed by its boundary faces within 1e-13; the inverted cell raises the flag of degenerate or inverted cells | ☑ `CellVolume` | <1 s |
| N10 | Order of rounding of `VOLATILE` variables, on which the geometry of the cells relies: (q − c1) − (c2 − c1) with q = 1, c1 = 2⁵⁴, c2 = 3 from the command line, through `VOLATILE` variables (a pair, an array) and through plain ones | −4 through the `VOLATILE` variables, the order of the source (plain ones printed: −2 with `-ffast-math` is allowed) | ☑ `VolatileOrder` | <1 s |
| N11 | Wall distance to a bilinear face, a unit check of the certified minimum (Bilinear_Face_Distance) on 203 cases: twisted (twist 1 to 1e-6), plane, non-convex, degenerate (a collapsed side, nearly collapsed by 1e-6 to 1e-12, two collapsed sides, three coincident corners, a bowtie, a pinch, a needle of aspect 1e6), turned at random and moved by 1e3, points on the normal down to 1e-9, beyond sides and corners, at the focal point of a saddle, adverse cases; and of the fast-path test (Face_Fast_Test) on 108 faces | within 288 u S of the exact distance of the stationarity equations (rational arithmetic), certificate within 256 u S, no piece cut at the deepest level, at most 1e5 pieces; the fast path taken by every plane face (convex, thin tapered up to aspect 1e5, triangles of a collapsed side also far from the origin, parallelograms, 2-D extruded), not by twisted faces and darts, its bound of the gap at least the exact gap | ☑ `BilinearDistance` | <1 s |

Current `ctest -L fast`: **24 tests, ~160 s** (Sod79, Einfeldt91, Noh87, Toro99,
ShuOsher, PrandtlMeyer, RiemannSmoke, NumericsSmoke, InputChecks,
TimeIntegration, BoundaryConditions, Conservation, OpenMPEquiv, Connections,
Restart, OutputFiles, MMSResiduals, ViscousPatch, AblatingWall, MeshGeometry,
QCR2000Stress, CellVolume, VolatileOrder, BilinearDistance). MPI builds also register the `needs-mpi` tier
(`test/mpi/`; `MPIConnectionFlux`, `MPIConnectionOrder`: the two-block cases with
one block per rank; `MPIDecompositionMap`).

The time is CPU time at one thread, one test at a time, measured on one machine
with 192 logical CPUs at a load average of ~190, where the 12 fast tests of main
take ~41 s, and at ~100 for the times of the later changes.

These catch most regressions (metrics, interfaces, solver crashes, restart/IO,
conservation) with no reference data or heavy dependencies.

The 4-grid Shu–Osher Richardson study (`ShuOsherConvergence`) uses the N=1600
solution and matplotlib, so it lives in the **validation** tier rather than the
fast pre-push suite (~4 s for the three solved grids).

`PointImplicitSteady` ([`2D/viscous/point-implicit-steady`](2D/viscous/point-implicit-steady))
converges a short hypersonic plate (6 × 16 cells; Spalart–Allmaras, SST, and SST and
Wilcox 2006 with `k-energy-coupling = .true.`) with the point-implicit factor, then
continues it at half the (cfl, vnn) and with implicit residual smoothing: the wall
shear and heat flux must not move, to 1e-9. Its 14 000 iterations per model (about
12 s of CPU at one thread each, ~50 s for the four) put it in the **validation**
tier as well.

`WedgeAxisStep` ([`2D/viscous/wedge-axis-step`](2D/viscous/wedge-axis-step)) checks the local time step
next to the axis of a 1° wedge (diffusion-dominated, Prandtl number 10, a little swirl; square
cells, cells twice as long in x, and a mesh with the axis on the last j face): 2000 steps
at RK3 vnn 0.62 and Heun vnn 0.5 with no NaN; and, on cells four times as long in x with
the Prandtl number of air and with 0.9, 2000 steps at the recommended vnn (RK3 0.5, Heun
0.4). Its eight runs of 2000 iterations (~155 s of CPU, TIMEOUT 600 s) put it in the
**validation** tier too.

`MMSViscousSkewed` ([`2D/viscous/mms/skewed`](2D/viscous/mms/skewed)), `MMSViscous3D`
([`3D/viscous/mms`](3D/viscous/mms)) and `PlanarNozzleGrids`
([`2D/viscous/planar-nozzle`](2D/viscous/planar-nozzle)), the short forms of N2 on general
grids and of N8, take ~76 s, ~20 s and ~72 s of CPU time at one thread: they are in the
**validation** tier, with a TIMEOUT of 900 s as `MMSViscous`.

---

## 5. Infrastructure (do this first)

The current blocker is that only 2 cases are in CTest. Restructure before adding cases:

1. ◐ **Wire cases into [`test/CMakeLists.txt`](CMakeLists.txt)** with `LABELS` via the
   `mose_add_test()` / `mose_add_fast_script()` helpers. *Done:* the dependency-free
   cases with proper exit codes (Sod79, Einfeldt91, Noh87, Toro99, PrandtlMeyer) plus
   the new fast scripts. *Remaining:* the `verify.py` cases that don't yet exit non-zero
   on failure (most 2-D/viscous cases) and the `needs-*` heavy cases.
2. ☑ **Pre-push hook** → `ctest -L fast` ([`.githooks/pre-push`](../.githooks/pre-push)).
3. ☑ **PR CI** → `ctest --label-exclude "needs-cantera|needs-sundials|needs-mpi"`
   ([`ctest.yml`](../.github/workflows/ctest.yml)); a second job with `USE_MPI=ON`
   → `ctest -L needs-mpi`.
4. ☐ **Nightly / `workflow_dispatch`** → full suite with Cantera+Sundials+MPI enabled.
5. ☐ Standardize on one verifier name (currently both `check.py` and `verify.py`) and a
   common exit-code contract (0 = pass, non-zero = fail) — needed before the remaining
   `verify.py` cases can join tier 1.

---

## 6. Tiering summary

| Tier | Trigger | Selector | Contents |
|---|---|---|---|
| 0 — fast | pre-push hook | `ctest -L fast` | F1–F7 |
| 1 — PR CI | PR to `main` | `ctest -LE needs-*`; `ctest -L needs-mpi` (second job, MPI build) | fast + dependency-free validation; MPI tier |
| 2 — nightly | schedule / dispatch | full | flames, reactive, equilibrium, MPI, parallel-equivalence |

---

## 7. Suggested implementation order

1. CTest label restructure (#5.1) + pre-push/CI selector wiring (#5.2–5.3).
2. F1, F3, F4 (no reference data needed) — makes the pre-push hook meaningful.
3. F5–F7, then N3/N4 (cheap, high coverage of numerics paths).
4. N1/N2 order-of-accuracy (multi-grid, slower → tier 1).
5. P3–P5 physics validation (tier 1/2). (Turbulent flat plate SA/SST/Wilcox2006
   exist as cases outside CTest: see §2 and §3.)
