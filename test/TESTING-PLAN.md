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
  runs everything except the `needs-*` heavy-dependency labels.
- **CI build is minimal** ([`.github/workflows/ctest.yml`](../.github/workflows/ctest.yml)):
  `USE_MPI=OFF`, `USE_TECIO=OFF`, `USE_SUNDIALS=OFF`, `USE_CANTERA=OFF`.

### Hard constraint for anything on GitHub's runner
The minimal CI build means **no finite-rate chemistry, no equilibrium (Cantera/CEA),
no MPI, no binary Tecplot (TecIO)**. Every flame / reactive / MPI / equilibrium case
is therefore off-limits for pre-push and PR CI — those belong in a nightly/local tier.

---

## 2. Coverage map (exercised vs. available)

| Capability | Registered options | Tested today | Gap |
|---|---|---|---|
| Equations | euler, navier-stokes | both | — |
| Time scheme | euler, RK2, RK3 | RK3 only | RK2, forward-euler untested |
| Integration vars | cons, prim, prec | cons (+ prec in Gresho) | `prim` untested |
| Space recon | first-order, MUSCL | MUSCL | first-order path untested |
| Limiters | vanleer, vanalbada, minmod, superbee, mc | vanleer (mostly) | 4 of 5 never smoke-tested |
| Riemann solvers | ~14 | ~6 in V&V | half never executed in CI |
| Low-Mach (Mco) | floor | Gresho | ok |
| IRS / multigrid / CFL-ramp | yes | implicit only | no dedicated check |
| Turbulence | SA(+R/RC/comp), SST, Wilcox2006 | SA, SST, Wilcox2006 (flat plate) | RC/comp variants untested |
| Chemistry | frozen, finite-rate, equilibrium | all (needs Cantera/Sundials) | can't run in CI |
| Multi-block | core feature | implicit only | **no interface correctness test** |
| Ablative wall BCs (503/504/505) | melting, pyrolysis, surface reactions | 503/504 vs closed form (`AblatingWall`) | 505 needs a CO/CO2 thermo table |
| Restart | core feature | none | untested |
| I/O formats | tec ascii/bin, vtk ascii/raw | ascii only | round-trip untested |
| Parallelism | OpenMP, MPI | runs, not verified | **no serial==parallel check** |
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
- ☐ **I2 — Single-block vs multi-block equivalence**: same problem as 1 block and 2×2
  blocks; fields identical to round-off. Direct halo/interface-flux test.
- ☑ **I3 — Restart round-trip**: N steps vs (N/2 → write → restart → N/2); match to
  ~2e-15. *Done:* `test/fast/restart` (restart = `newrun=false`, reads OUTPUT/field.tec).
- ☑ **I4 — Conservation invariant**: closed slip-wall box, shock-tube IC; total mass
  conserved to round-off (~2e-13). *Done:* `test/fast/conservation`. (Energy invariant
  is a possible future addition.)
- ☑ **I5 — Serial vs OpenMP equivalence**: 1 vs 4 threads bit-identical, also with chimera faces meeting at
  the edges of a block. *Done:*
  `test/fast/openmp-equiv`. (MPI variant → nightly.)

---

## 4. Pre-push fast suite (GitHub-runner friendly)

Rules: minimal build only, tiny meshes (≤ ~1k cells), ≤ a few seconds each, assert
**hard invariants** (no digitized references, no plotting, no Cantera). Target **< 60 s** total.

| # | Test | Asserts | Status | ~Cost |
|---|---|---|---|---|
| F1 | Freestream preservation (stretched multi-block) | residual ≈ machine-zero | ☐ | <2 s |
| F2 | Sod coarse (existing) | shock/contact positions, positivity | ☑ `Sod79` | ~1 s |
| F3 | Conservation in closed box | mass conserved to round-off | ☑ `Conservation` | <1 s |
| F4 | Single- vs multi-block equivalence | identical field | ☐ | <3 s |
| F5 | Riemann smoke matrix | no NaN, ρ,p>0 | ☑ `RiemannSmoke` | ~4 s |
| F5b | Numerics sweep (limiters/schemes/recon) | no NaN, ρ,p>0 | ☑ `NumericsSmoke` | ~5 s |
| F5e | Time control (1-D duct, time-dependent inlet 402) | order in time ≥ 1.8 (RK2, RK3); table held before its first instant; time-threshold reached | ☑ `TimeIntegration`/`TimeControl` | ~5 s |
| F6 | Restart round-trip; probes at a restart, the file of one of them removed | identical to uninterrupted run; the missing probe file written again, the other appended, with the probe rows of the uninterrupted run | ☑ `Restart` | ~1 s |
| I5 | Serial vs OpenMP equivalence: the Sod tube; a chimera case, a block overset on another whose four chimera faces meet at its edges (Navier–Stokes; it passes on the parent of its commit too: the race it guards shows under ThreadSanitizer with Archer, in three runs of three) | 1 vs 4 threads bit-identical | ☑ `OpenMPEquiv` | ~3 s |
| F9 | Viscous stencil of chimera faces: two blocks of 8 x 6 cells (8 x 6 x 4 in 3-D) overlapping by 4 cells in x, each ghost cell on a cell of the other block; one step of a steady run with u = G y (2-D), u = G (y + z) (3-D), and a radial velocity on a 2-D wedge of 1 degree against the same wedge as one block | y and z momentum residual zero within 1e-9 of μ G / H, also in the cells at the ends of the chimera faces (1/4 of it with the ghost cell copied beyond the edge); on the wedge the momentum residual of every cell that of the one block within 1e-9 of μ G / H (half of v/r at the chimera faces with the ghost cell as its own k neighbour) | ☑ `Connections`/`ChimeraFace` | ~1 s |
| F7 | Symmetry (coarse forward-step/wedge) | top/bottom mirror preserved | ☐ | <3 s |
| F8 | Shu–Osher coarse (N=200) | L1 vs N=1600 ref < 5 %, finite | ☑ `ShuOsher` | ~2 s |
| P1 | Ablative wall BCs (503/504); the melting column also without `schmidt` (default 0: one species needs no `diffusion.dat`) | wall state vs closed form, both face orientations, inert-surface fallback | ☑ `AblatingWall` | ~12 s |
| I9 | Exit codes of the error stops: the Sod tube as a viscous mixture of four species without `schmidt` and without `diffusion.dat` ([thermo](../docs/theory/thermo.md)); the Sod tube in air with a probe of an unknown variable; a steady run with the soot model LL91, with and without Spalart–Allmaras (six species with C(gr)); the Sod tube in air, steady, with sol-variables = none | non-zero exit code and the diffusion-table error, no solution; exit status 1 and the name of the variable, no solution; exit status 0 and a diagnostic output that names every variable it writes, the local time step under dt; exit status 1 and a message that names thermo or transport | ☑ `InputChecks`/`ExitCodes` | ~4 s |
| F11 | Local time step (square and skewed 2-D cells, cubic 3-D cells, heat-conduction and species-diffusion cells, 1-D at cfl 0.99; heat conduction on a periodic block of 8³ rhombohedra at 115°, Heun at vnn 0.5, 600 steps) | stable, uniform state kept, no NaN; the temperature wave of the rhombohedra decays | ☑ `TimeIntegration`/`LocalTimeStep` | ~9 s |
| F12 | Turbulence sources in time-accurate runs (wall channel, SA and SST, Heun and RK3; SST-RC, RK3, one step) | time order of nu~, k, omega ≥ 1.8 (Heun) and ≥ 2.7 (RK3); no negative stage value; the start-up warnings of k-omega and of -RC both printed (SST-RC) | ☑ `TimeIntegration`/`TurbulenceSources` | ~9 s |
| F13 | Spalart–Allmaras on constructed fields: variable density (periodic field, no wall); SA-RC and SST-RC on a 3-D flow with vorticity along x, y and z (cube of 8³ cells) and on the same flow turned (x to y, y to z, z to x), one step each; SA-RC on a uniform shear in a frame turning about z and on a uniform flow (8 × 8 cells), one step each; SA-RC and SA on air at rest next to a wall and on air flowing along it on cells stretched along the flow (8 × 8 cells), one step each | right-hand side of rho nu~ = standard conservation form, error < 0.02 and order > 1.5; change of rho nu~ (rho k and rho omega with SST-RC) the same in corresponding cells of the two 3-D runs within 1e-9 of its largest value; change of rho nu~ = dt rho c_b1 f_r1 S~ nu~ with the f_r1 of the specification within 1e-6 (shear), no NaN and no change (uniform); the same change of rho nu~ with SA-RC and SA within 1e-12 (rest; stream, off the wall row) | ☑ `MMSResiduals` | ~3 s |

Current `ctest -L fast`: **16 tests, ~76 s** (Sod79, Einfeldt91, Noh87, Toro99,
ShuOsher, PrandtlMeyer, RiemannSmoke, NumericsSmoke, TimeIntegration,
Conservation, OpenMPEquiv, Connections, Restart, InputChecks, MMSResiduals,
AblatingWall).

The time is CPU time at one thread, one test at a time, measured on one machine
(load average ~190 on 192 logical CPUs), where the 12 fast tests of main take ~41 s.

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
   ([`ctest.yml`](../.github/workflows/ctest.yml)).
4. ☐ **Nightly / `workflow_dispatch`** → full suite with Cantera+Sundials+MPI enabled.
5. ☐ Standardize on one verifier name (currently both `check.py` and `verify.py`) and a
   common exit-code contract (0 = pass, non-zero = fail) — needed before the remaining
   `verify.py` cases can join tier 1.

---

## 6. Tiering summary

| Tier | Trigger | Selector | Contents |
|---|---|---|---|
| 0 — fast | pre-push hook | `ctest -L fast` | F1–F7 |
| 1 — PR CI | PR to `main` | `ctest -LE needs-*` | fast + dependency-free validation |
| 2 — nightly | schedule / dispatch | full | flames, reactive, equilibrium, MPI, parallel-equivalence |

---

## 7. Suggested implementation order

1. CTest label restructure (#5.1) + pre-push/CI selector wiring (#5.2–5.3).
2. F1, F3, F4 (no reference data needed) — makes the pre-push hook meaningful.
3. F5–F7, then N3/N4 (cheap, high coverage of numerics paths).
4. N1/N2 order-of-accuracy (multi-grid, slower → tier 1).
5. P3–P5 physics validation (tier 1/2). (Turbulent flat plate SA/SST/Wilcox2006
   already exist; remaining turbulence gap is the SA RC/comp variants.)
