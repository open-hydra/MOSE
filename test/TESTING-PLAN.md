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
| Integration vars | cons, prim, prec | cons (+ prec in Gresho) | `prim` untested |
| Space recon | first-order, MUSCL | both (first-order: Noh87, LimiterDefaults) | — |
| Limiters | vanleer, vanalbada, minmod, superbee, mc | vanleer (mostly) | 4 of 5 never smoke-tested |
| Riemann solvers | ~14 | ~6 in V&V | half never executed in CI |
| Low-Mach (Mco) | floor | Gresho | ok |
| IRS / multigrid / CFL-ramp | yes | implicit only | no dedicated check |
| Turbulence | SA(+R/RC/comp), SST, Wilcox2006 | SA, SST, Wilcox2006 (flat plate) | RC/comp variants untested |
| Chemistry | frozen, finite-rate, equilibrium | all (needs Cantera/Sundials) | can't run in CI |
| Multi-block | core feature | faces between blocks: one block against two and records sorted against reversed, bit for bit (`Connections`/`ConnectionFlux`, `Connections`/`ConnectionOrder`; with MPI `MPIConnectionFlux`, `MPIConnectionOrder`) | 3-D faces with permuted index directions, connections with multigrid and with RANS untested |
| Ablative wall BCs (503/504/505) | melting, pyrolysis, surface reactions | 503/504 vs closed form (`AblatingWall`) | 505 needs a CO/CO2 thermo table |
| Restart | core feature | none | untested |
| I/O formats | tec ascii/bin, vtk ascii/raw | ascii; `vtk raw` files are well-formed XML with the right names (`OutputFiles`/`VtkNames`) | round-trip untested |
| Parallelism | OpenMP, MPI | 1 vs 4 threads (`OpenMPEquiv`); MPI ranks × threads bit-identical in the `needs-mpi` tier (`test/mpi/`) | MPI tier runs only in MPI builds, not on the CI runner |
| 3D | full | source terms only | no 3D NS verification |

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
- ☑ **N3 — Riemann-solver smoke matrix**: coarse Sod, few steps, once per solver;
  assert positivity (ρ,p>0) and no NaN. *Done:* `test/fast/riemann-smoke` (11
  shock-capturing solvers; low-Mach family deferred to a Gresho smoke).
- ☑ **N4 — Limiter / scheme sweep**: all 5 limiters + 3 time schemes + recon +
  integration-variables; assert positivity. *Done:* `test/fast/numerics-smoke`.
  (Monotonicity assertion is a possible future tightening.)

### Physics
- ☑ **Turbulent flat plate, SA / SST / Wilcox2006** — already covered, one case per
  model under `test/2D/viscous/flat-plate-turbulent/{SA,SST,Wilcox2006}/`, each
  vs NASA CFL3D/FUN3D reference Cf. *No further turbulent flat-plate tests needed.*
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
- ☑ **I4 — Conservation invariant**: closed slip-wall box, shock-tube IC; total mass
  conserved to round-off (~2e-13). *Done:* `test/fast/conservation`. (Energy invariant
  is a possible future addition.)
- ☑ **I5 — Serial vs OpenMP equivalence**: 1 vs 4 threads bit-identical, also with chimera faces meeting at
  the edges of a block. *Done:*
  `test/fast/openmp-equiv`. *MPI variant (needs-mpi tier):* `test/mpi/boundary-records`
  (`MPIBoundaryRecords`) — a strip of four viscous blocks bit-identical on 1×1, 2×1,
  4×1 and 2×2 ranks × threads, and with the `bc.txt` records in reverse order (1 and 2
  ranks, with and without threads), the wall output included; the same with the
  Spalart–Allmaras model, ordered and reversed (1×1, 1×4, 2×2);
  `test/mpi/vertex-stencil` (`MPIVertexStencil`) — a
  2×2 box of viscous blocks with a vortex on the common vertex, bit-identical on 1×1,
  1×4, 2×1, 4×1 and 2×2 ranks × threads; with chemistry = equilibrium, bit-identical
  on 1×1 and 2×1.
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
| F5e | Time control (1-D duct, time-dependent inlet 402) | order in time ≥ 1.8 (RK2, RK3); table held before its first instant; time-threshold reached | ☑ `TimeIntegration`/`TimeControl` | ~5 s |
| F6 | Restart round-trip; probes at a restart, the file of one of them removed | identical to uninterrupted run; the missing probe file written again, the other appended, with the probe rows of the uninterrupted run | ☑ `Restart` | ~1 s |
| I11 | Wall output in vtk raw (2x2-block box of 8 × 8 cells, adiabatic no-slip walls (301) on its four sides, 5 RK3 steps) | the run ends; every wall file `wallB<b>F<f>.vts` has the extent of its face (WholeExtent and Piece Extent, from `fieldBlock<b>.vts`) and its points, as raw appended data | ☑ `OutputFiles`/`WallOutput` | <1 s |
| I5 | Serial vs OpenMP equivalence: the Sod tube; a chimera case, a block overset on another whose four chimera faces meet at its edges (Navier–Stokes; it passes on the parent of its commit too: the race it guards shows under ThreadSanitizer with Archer, in three runs of three) | 1 vs 4 threads bit-identical | ☑ `OpenMPEquiv` | ~3 s |
| I12 | Names of the VTK outputs (steady viscous run on a 2 × 2-block box with no-slip walls, 5 iterations, `vtk raw`) | every `.vts` that `field.vtm`, `wall.vtm` and `diagnostic.vtm` list parses as XML, and its data arrays are named as the variables that MOSE writes for the case (field, walls, diagnostic), with no other array | ☑ `OutputFiles`/`VtkNames` | <1 s |
| F9 | Viscous stencil of chimera faces: two blocks of 8 x 6 cells (8 x 6 x 4 in 3-D) overlapping by 4 cells in x, each ghost cell on a cell of the other block; one step of a steady run with u = G y (2-D), u = G (y + z) (3-D), and a radial velocity on a 2-D wedge of 1 degree against the same wedge as one block | y and z momentum residual zero within 1e-9 of μ G / H, also in the cells at the ends of the chimera faces (1/4 of it with the ghost cell copied beyond the edge); on the wedge the momentum residual of every cell that of the one block within 1e-9 of μ G / H (half of v/r at the chimera faces with the ghost cell as its own k neighbour) | ☑ `Connections`/`ChimeraFace` | ~1 s |
| F7 | Symmetry (coarse forward-step/wedge) | top/bottom mirror preserved | ☐ | <3 s |
| F8 | Shu–Osher coarse (N=200) | L1 vs N=1600 ref < 5 %, finite | ☑ `ShuOsher` | ~2 s |
| P1 | Ablative wall BCs (503/504); the melting column also without `schmidt` (default 0: one species needs no `diffusion.dat`) | wall state vs closed form, both face orientations, inert-surface fallback | ☑ `AblatingWall` | ~12 s |
| I9 | Exit codes of the error stops: the Sod tube as a viscous mixture of four species without `schmidt` and without `diffusion.dat` ([thermo](../docs/theory/thermo.md)); the Sod tube in air with a probe of an unknown variable; a steady run with the soot model LL91, with and without Spalart–Allmaras (six species with C(gr)), and the same state with `turbulence = SAcomp`; the Sod tube in air, steady, with sol-variables = none; the Sod tube that diverges (CFL 20), with an initial condition of too few variables, with a rotating frame of zero axis, of NaN axis or of infinite origin, with chemistry = equilibrium and no `INPUT/composition.txt`, with a value outside the rule of the registry: an interval, a list of reals, a list of integers; the Sod tube with a NaN velocity in one cell and with an infinite coordinate of one node; the duct of BCInput with a NaN and an infinite pressure, and a NaN time, in the table of its inlet 402; the duct with a NaN and an infinite total pressure and a NaN alpha of an inlet 401, an infinite total temperature of its inlet 402 and a NaN roughness of a wall 302 in its bc.txt; the duct with its inlet a face 410 whose file holds a NaN pressure and an infinite SOLUTIONTIME, and the same file missing, cut in the values of its second zone, with a word among them, with I=two in the header of its second zone, with I = 3 or j=2 there (no I= or no J=), with one value more in its first zone and with its two zones in decreasing time; input.ini missing, empty, of comment lines only; the Sod tube with sol-dtime = NaN, Infinity and 0,5, sol-diter = seven and 2.5, sol-overwrite = yes and levels = 2.0 of [MOSE-Multigrid]; a valid run with finite-rate chemistry and exclude-blocks = 1; exclude-blocks with a word that is not a block number or with a block beyond the mesh, which stops at set-up; the tube filled with steam far from chemical equilibrium, its block excluded from chemistry = equilibrium, on one grid and on two multigrid levels; on 2 MPI ranks (`MPIExitCodes`, needs-mpi tier), cfl = -1, without input.ini, the NaN velocity with the block on rank 1 alone and the infinite coordinate, and cfl = -1 written during the run in the input.ini of rank 1 alone (two folders, MPMD), and the initial condition of rank 1 alone cut in the middle (two folders) | non-zero exit code and the diffusion-table error, no solution; exit status 1 and the name of the variable, no solution; exit status 0 and a diagnostic output that names every variable it writes, the local time step under dt; with SAcomp, exit code 1 to 128 at set-up, the message and no diagnostic output; exit status 1 and a message that names thermo or transport; exit code 1 to 128 and the error message (for input.ini, the reason it cannot be read, and no solution); for a key of input.ini whose value does not read as its type, exit code 1, the file, the section, the key and the value as written and the STOP line, no solution; for a value that is not finite, exit code 1, the place of the value and the STOP line, no solution; for the table, exit code 1, the line and the column, no solution; for bc.txt, exit code 1, the line, the record, its type and the name of the value, no solution; for the file of the face 410, exit code 1, the file, the zone, the line of its header, the variable and the cell, no solution, and for a file that cannot be read, exit code 1, the file and, for a zone, the zone and the line of its header (with the two times for zones out of time order, and the form of I and J to write for a header without I= or J=), no solution; exit code 0 and block 1 under the inactive blocks of the set-up summary; the field of the excluded block bit-identical to the run without chemistry, and different with the block not excluded; on 2 ranks, with cfl = -1 exit code 1 within 60 s, the message and the STOP line of both ranks, without input.ini exit code 1 to 128 within 60 s and the message, and for a value that is not finite exit code 1, the place written once and the STOP line of both ranks; with rank 1 alone, exit code 1 at iteration 20 and the rule, not a hang; with the initial condition of rank 1 cut, exit code 1, Initial conditions FAIL, the error of the read written once and the STOP line of both ranks | ☑ `InputChecks`/`ExitCodes` | ~18 s |
| F10 | Quasi-2D with symmetry k faces; 1-D with symmetry j and k faces; a channel of variable section in a layer between parallel planes inclined to z and a 1-D channel of variable section, steady subsonic flow | p, u, v uniform; round-off in w does not grow (Heun, RK3); in 1-D, round-off in v and w does not grow; in the layer and in the channel the velocity across is integrated, with the warning, and p within 0.03 p0 of the quasi-1-D flow | ☑ `TimeIntegration`/`Quasi2DSymmetry` | ~8 s |
| F11 | Local time step (square and skewed 2-D cells, cubic 3-D cells, heat-conduction and species-diffusion cells, 1-D at cfl 0.99; heat conduction on a periodic block of 8³ rhombohedra at 115°, Heun at vnn 0.5, 600 steps) | stable, uniform state kept, no NaN; the temperature wave of the rhombohedra decays | ☑ `TimeIntegration`/`LocalTimeStep` | ~9 s |
| F12 | Turbulence sources in time-accurate runs (wall channel, SA and SST, Heun and RK3; SST-RC, RK3, one step); Wilcox 2006 on a 3-D linear flow, divergence-free, with S^ zero on the diagonal (cube of 6³ cells), one step; Wilcox 2006 and SST with omega-wall-bc = asymptotic on the channel at rest, omega at the wall value 6 nu/(beta d²) of the model, one step | time order of nu~, k, omega ≥ 1.8 (Heun) and ≥ 2.7 (RK3); no negative stage value; the start-up warnings of k-omega and of -RC both printed (SST-RC); change of rho omega = dt (2 gamma rho S:S - beta_0 f_beta rho omega²) within 1e-6; in the wall cells, change of rho omega = -dt beta rho omega² (the destruction alone) within 1e-6 | ☑ `TimeIntegration`/`TurbulenceSources` | ~9 s |
| F13 | Spalart–Allmaras on constructed fields: variable density (periodic field, no wall); SA-RC and SST-RC on a 3-D flow with vorticity along x, y and z (cube of 8³ cells) and on the same flow turned (x to y, y to z, z to x), one step each; SA-RC on a uniform shear in a frame turning about z and on a uniform flow (8 × 8 cells), one step each; SA-RC and SA on air at rest next to a wall and on air flowing along it on cells stretched along the flow (8 × 8 cells), one step each; the viscous stress on the axisymmetric wedge (200): Homann's stagnation-point flow and a polynomial field with dilatation and swirl (1° wedge, axis included); Spalart–Allmaras next to a wall (fluid at rest, nu~ = a y, chi < 1e-3) | right-hand side of rho nu~ = standard conservation form, error < 0.02 and order > 1.5; change of rho nu~ (rho k and rho omega with SST-RC) the same in corresponding cells of the two 3-D runs within 1e-9 of its largest value; change of rho nu~ = dt rho c_b1 f_r1 S~ nu~ with the f_r1 of the specification within 1e-6 (shear), no NaN and no change (uniform); the same change of rho nu~ with SA-RC and SA within 1e-12 (rest; stream, off the wall row); viscous force vs exact: Homann, radial and axial error < 0.02 and order > 1.5, azimuthal zero; polynomial, radial and azimuthal error < 0.02 and order > 1.5, axial exact; right-hand side of rho nu~ below 0.05 of (cb2/sigma) rho a^2 in every cell, wall cells included (grad nu~ = a there) | ☑ `MMSResiduals` | ~7 s |
| I7 | Faces between blocks (101, 201): periodic shock tube as one and as two blocks, HLLC, HLLC+Tramel, HLLC+Chen, HLLC-PC (low Mach) | two blocks = one block bit for bit; mass, momentum, energy conserved to 1e-13 | ☑ `Connections`/`ConnectionFlux` | ~2 s |
| I7b | Faces between blocks of the same parity (2-D, index directions permuted across the cut), HLLC+Tramel, Navier–Stokes; the same with the block below the cut one cell thick across it, and two cells thick with two multigrid levels | bc.txt records sorted = reversed, bit for bit; the thin blocks stop the run at set-up with a non-zero exit code and the message (the record, the block across, its cells, the level) | ☑ `Connections`/`ConnectionOrder` | <1 s |
| I8 | Block-to-rank map `INPUT/decomposition.map` on one rank, also a map written by MDB | rank-0 map = no map bit for bit; misfit maps stop at setup with the fault named | ☑ `InputChecks`/`DecompositionMap` | ~6 s |
| I10 | input.ini read again during a run (2x2-block viscous box, 48 × 48 cells per block: ini-diter = 20, 40 iterations; an input.ini with the required keys only, and one that also writes integration-variables = prim, riemann-solver = HLLE and flux-limiter = minmod, ini-diter = 5, 10 iterations); a comment line with "=" before cfl in [MOSE-Numerics], cfl = 0.4, 3 iterations, and the same run without the line | cfl = -1 stops the run at setup and, written during the run, at the reading of iteration 20 with the message of the rule; with the required keys only, and with the three keys written, the run reaches iteration 10 with the residual history of the same run without readings, byte for byte; with the comment line the run writes the residual history of the run without it, byte for byte, and the setup names no key as not read | ☑ `InputChecks`/`RuntimeInput` | ~14 s |
| F17 | Inlets: the ambient inlet (407) with the flow direction `normal` (planar duct, air entering at 50 m/s); inlets at total conditions (1-D tube of air, one forward Euler step): 407 and 401 from rest at p0 = 2 p and 36 p, an inlet the reservoir cannot feed, a 407 that it cannot feed with the gas moving in at 0.5 a and 1.2 a, and with p_back above p_stop | same field as the equivalent angles bit for bit, no transverse velocity; mass flux of each inlet face = that of the face state (subsonic, choked at normal Mach 1, slip wall, outlet at min(p_back, p_stop)) within 1e-9 of rho a | ☑ `BoundaryConditions`/`Inlets` | ~2 s |
| F18 | Wall distance (fluid at rest, nu~ uniform; convex corner of two walls, curved wall under a leaning mesh, nozzle lip on a 1° wedge; channel and pipe on a wedge periodic in x, 3-D sector closed by rotational periodicity, the channel with its periodic sides matched to 1e-8, 3-D slab periodic in k by type 200, each with the wall on part of the period; 3-D channel periodic in x and z with periods a hundredth of its height, wall on the whole period; the 2-D channel ten periods from the origin; a 2-D staggered row of plane blades 1.7 periods long; a 3-D block periodic along three nearly dependent translations; the 2-D channel turned by 30° with its coordinates rounded to 5, 4 and 2 significant digits; 3-D wall of twisted faces) | distance read back from the Spalart–Allmaras source = exact distance to the wall and its periodic images within 1e-9 (1e-6 with the sides matched to 1e-8; for the turned channel, between the distances with the translation of each periodic face); with 2 digits, stop at set-up with the message on the periodic translations; start-up report of the twist bound of the wall faces | ☑ `MeshGeometry`/`WallDistance` | ~5 s |
| F19 | Mesh type (plane slab with y falling along j, right-handed wedge with the axis on its last j face, a wedge of 5 degrees in two blocks with coordinates of 5 significant digits, a plane slab far from the x axis (y from 1 to 1.001, 1e-4 thick), a wedge of 1 degree with coordinates of 2 significant digits, a plane slab with its k planes tilted by 2 % of its thickness, 3-D sectors of 90 degrees with the periodic plane turned by half a cell more (100 cells in k) or a tenth (20 cells) but at two nodes, a wedge block with a plane block, 3-D slab periodic in k with type 200) | plane and meridional wall distance (to the wall of the nodes as written) within 1e-9, and the type in the set-up summary; set-up stop with a non-zero exit code (2 digits, tilted slab, shifted sectors, mixed blocks); uniform flow kept within 1e-12 | ☑ `MeshGeometry`/`MeshType` | ~2 s |
| F21 | Rotational periodicity of a 3-D sector (annulus sector from -45 to 45 degrees, type 200 on the k faces, radial flow with uniform density and pressure, 20 steps of RK3) | the columns k and NK + 1 - k mirror images (rho, u, v, p equal, w opposite) within 1e-10 of their scales | ☑ `Connections`/`SectorRadial` | <1 s |
| F24 | Boundary input (1-D duct): inlet 402 table under four file names (plain, with a directory, longer than 32 characters, quoted with a blank), a missing table, a table whose times decrease and one that starts after the run; property lines with an empty value, short of a value, with a word for a number or with a misspelt flow direction | the four runs bit-identical; the missing and the decreasing table stop the run at set-up with the name of the file; the late table gives the run of the same table with its first value at t = 0, bit for bit; the bad lines stop the run with their line, record and type | ☑ `BoundaryConditions`/`BCInput` | ~3 s |
| F23 | Rotating frame (Sod tube, rigid-body frame) with two `stationary-face-<n>` entries | the run ends within 60 s with exit code 0; the set-up summary lists the two stationary patches | ☑ `BoundaryConditions`/`RotatingFrame` | <1 s |
| P6 | Quadratic constitutive relation QCR2000, a unit check (QCR2000 is not selectable from the input): the stress of a 3-D velocity gradient with the eddy viscosity alone and with the laminar viscosity alone, called as the solver calls it | tau_QCR · n of the definition and the laminar stress, within 1e-12 | ☑ `QCR2000Stress` | <1 s |

Current `ctest -L fast`: **20 tests, ~140 s** (Sod79, Einfeldt91, Noh87, Toro99,
ShuOsher, PrandtlMeyer, RiemannSmoke, NumericsSmoke, InputChecks,
TimeIntegration, BoundaryConditions, Conservation, OpenMPEquiv, Connections,
Restart, OutputFiles, MMSResiduals, AblatingWall, MeshGeometry, QCR2000Stress).

The time is CPU time at one thread, one test at a time, measured on one machine
with 192 logical CPUs at a load average of ~190, where the 12 fast tests of main
take ~41 s, and at ~100 for the times of the later changes.

These catch most regressions (metrics, interfaces, solver crashes, restart/IO,
conservation) with no reference data or heavy dependencies.

The 4-grid Shu–Osher Richardson study (`ShuOsherConvergence`) uses the N=1600
solution and matplotlib, so it lives in the **validation** tier rather than the
fast pre-push suite (~4 s for the three solved grids).

`PointImplicitSteady` ([`2D/viscous/point-implicit-steady`](2D/viscous/point-implicit-steady))
converges a short hypersonic plate (6 × 16 cells, Spalart–Allmaras and SST) with the
point-implicit factor, then continues it at half the (cfl, vnn) and with implicit
residual smoothing: the wall shear and heat flux must not move, to 1e-9. Its 14 000
iterations per model (~20–40 s) put it in the **validation** tier as well.

`WedgeAxisStep` ([`2D/viscous/wedge-axis-step`](2D/viscous/wedge-axis-step)) checks the local time step
next to the axis of a 1° wedge (diffusion-dominated, Prandtl number 10, a little swirl; square
cells, cells twice as long in x, and a mesh with the axis on the last j face): 2000 steps
at RK3 vnn 0.62 and Heun vnn 0.5 with no NaN; and, on cells four times as long in x with
the Prandtl number of air and with 0.9, 2000 steps at the recommended vnn (RK3 0.5, Heun
0.4). Its eight runs of 2000 iterations (~155 s of CPU, TIMEOUT 600 s) put it in the
**validation** tier too.

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
   already exist; remaining turbulence gap is the SA RC/comp variants.)
