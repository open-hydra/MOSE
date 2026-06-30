# MOSE Testing Plan

Reference plan for expanding the MOSE test suite to cover **physics**, **numerics**,
and **implementation**, plus a fast tier suitable for a `git push` pre-push hook
running on GitHub's CI runner.

> Status legend: ☐ to do · ◐ in progress · ☑ done

---

## 1. State of play

- **Only 2 cases are wired into CTest** (`Sod79`, `PrandtlMeyer`, via `check.py` in
  [`test/CMakeLists.txt`](CMakeLists.txt)). The other ~18 cases ship a `verify.py`
  but run only through [`test/test.sh`](test.sh) — they are **not** in CTest/CI.
- **The pre-push hook** ([`.githooks/pre-push`](../.githooks/pre-push)) just calls
  `ctest`, so it inherits those same 2 tests.
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
| Restart | core feature | none | untested |
| I/O formats | tec ascii/bin, vtk ascii/raw | ascii only | round-trip untested |
| Parallelism | OpenMP, MPI | runs, not verified | **no serial==parallel check** |
| 3D | full | source terms only | no 3D NS verification |

---

## 3. Proposed new tests

### Numerics (biggest gap: order-of-accuracy verification)
- ☐ **N1 — Isentropic vortex convection** (2D, smooth, periodic): convect one period,
  grid-refine 32²/64²/128², **measure observed order** (→2 with MUSCL). Highest-value
  missing test. Extend with a Riemann × limiter sweep.
- ☐ **N2 — Method of Manufactured Solutions (NS)**: source-term forcing, grid refine,
  verify 2nd order on the *viscous* operator. Catches gradient/metric bugs N1 cannot.
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
- ☐ **P5 — Shu–Osher** (1D): smooth+shock interaction stressing high-order recon.

### Implementation / infrastructure (also serve as fast tests)
- ☐ **I1 — Freestream preservation / GCL**: uniform flow on a stretched, curvilinear,
  multi-block mesh; assert residual stays at machine-zero.
- ☐ **I2 — Single-block vs multi-block equivalence**: same problem as 1 block and 2×2
  blocks; fields identical to round-off. Direct halo/interface-flux test.
- ☐ **I3 — Restart round-trip**: N steps vs (N/2 → write → restart → N/2); identical.
- ☑ **I4 — Conservation invariant**: closed slip-wall box, shock-tube IC; total mass
  conserved to round-off (~2e-13). *Done:* `test/fast/conservation`. (Energy invariant
  is a possible future addition.)
- ☐ **I5 — Serial vs OpenMP equivalence**: 1 vs 4 threads identical (MPI variant → nightly).

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
| F6 | Restart round-trip | identical to uninterrupted run | ☐ | <3 s |
| F7 | Symmetry (coarse forward-step/wedge) | top/bottom mirror preserved | ☐ | <3 s |

Current `ctest -L fast`: **8 tests, ~27 s** (Sod79, Einfeldt91, Noh87, Toro99,
PrandtlMeyer, RiemannSmoke, NumericsSmoke, Conservation).

These catch most regressions (metrics, interfaces, solver crashes, restart/IO,
conservation) with no reference data or heavy dependencies.

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
