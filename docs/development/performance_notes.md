# MOSE Performance Notes

Running record of what has been tried for solver speed, what it was worth, and
how it was measured. **Numbers here are measured, not predicted** — anything not
yet measured is marked as such and lives in [Open items](#open-items).

If you optimise something, add a row. If you measure something and it turns out
not to help, add it to [Measured and rejected](#measured-and-rejected) — knowing
what does *not* pay is worth as much as knowing what does, and it stops the next
person repeating the experiment.

---

## How to measure (read before adding a row)

The protocol below exists because several early entries in this file were wrong,
and each was wrong for one of these reasons.

0. **Pick the host by its load, and run a control.** `s3` (sprop3) is the build
   host, but it is heavily shared — loads of 150–215 on 192 threads have been
   normal — and there the run-to-run spread reaches 12–20%, which swamps
   anything under 10%. **sprop2 is much quieter** (~20 of 96) and the same
   comparison there had a base spread of 1.5%. Its oneAPI is older (2023.1 vs
   2024.1), so *build* on s3 and *run* on sprop2; s3-built binaries execute
   there. Always include a case the change cannot possibly affect: a control
   showing a 5% "gain" is how you learn the measurement is worthless.
1. **Interleave the A/B.** Run `base, new, base, new`, not all of one then all of
   the other. The build hosts are shared — background load of 75–90 from other
   users is normal — and interleaving is what makes the comparison survive it.
   Report the **minimum** of at least 3 repetitions, not the mean.
2. **Check bit-identity of the output, not just the runtime.** A change that is
   both faster and value-preserving is a different thing from a change that is
   faster because it does less work. `cmp` the whole `OUTPUT/` directory.
3. **When results do differ, normalise to each variable's own scale.** Pointwise
   relative differences on trace radicals or on a machine-zero transverse
   velocity are meaningless. Compare `max|Δ|` against `max|value|` per variable.
4. **To prove two builds do the same work, force strict FP.** Build both with
   `-fp-model=precise -ffp-contract=off`. If the outputs are then bit-identical,
   the work is provably identical and any remaining time difference is code
   speed. This technique settled the `-ipo` chemistry question below.
5. **Time the solver with `timer-diter`, not with the job's wall clock.**
   `MOSE_Mod_Timers` is the measurement instrument: it reports the time inside
   the iteration loop, as the maximum over the ranks, separated from set-up and
   from output, plus the exchange and collective waits
   ([Timing](../user/using.md#timing)). It replaces the old
   `Time of operation` line, which was CPU time on rank 0 divided by the thread
   count — a figure that measured wall clock only when every thread was busy
   and no rank ever blocked. Entries below dated before that change quote it.
6. **Do not compare VTune totals across binaries with different inlining.**
   VTune's user-mode sampler instruments call/return, so a heavily inlined binary
   makes fewer calls and looks artificially cheap. Use wall clock for the verdict
   and VTune only for the breakdown *within* one binary. (Hardware sampling is
   Intel-only and unavailable on the AMD build hosts.)

Reference cases used below, all single-threaded:

| Case | What it stresses |
|---|---|
| `test/2D/euler/supersonic-forward-step` | convective path, no chemistry |
| `test/2D/viscous/swbli` | diffusive + turbulence; production-like |
| `test/1D/Fer14` | finite-rate chemistry via the `general` evaluator |
| `test/1D/premixed-flame` | finite-rate chemistry via a hard-coded mechanism (ZK) |

For SWBLI, turn the periodic dumps off (`sol-diter`/`ini-diter` past
`iter-threshold`) before timing anything. At the shipped cadence the case writes
~100 MB of ASCII per dump, which dominates the wall clock and dilutes whatever
you are trying to measure.

---

## Landed

### 1. Inline integer power in the Arrhenius rate loops

`FLINT_Lib_Chemistry_wdot::general` raised concentrations to reaction orders with
`coi(is)**nint(...)`. A run-time integer exponent lowers to a call into libgcc's
`__powidf2`, which was **4.8% of Fer14 self time** on its own — and, more
importantly, the call sat in the innermost loop and stopped the whole
rate-of-progress loop from optimising. Replaced with a `pure function ipow` that
expands orders 0–3 inline and falls back to the intrinsic otherwise, so results
are unchanged.

| Case | Before | After | |
|---|---|---|---|
| Fer14 | 13185 ms | 10536 ms | **−20.1%** |

Output **bit-identical**, 26/26 files. The gain is far larger than the 4.8% the
profile attributed to `__powidf2` — that 4.8% was only the call itself.

Only the Arrhenius branch was converted. The Troe and Lindemann branches use a
*real* exponent (`coi(is)**ni1_troe_tab(is,ir)`), which is libm `pow`, not
`__powidf2`; converting those would change their semantics and was left alone.

### 2. Interprocedural optimisation (`-ipo`), on by default

MOSE's hot path is a chain of small cross-module leaf functions — `rlimiter`,
`f_ss`, `f_cp`, `h`, `f_Rtot`, `co_rotot_Rtot`, `check_gas_state` — none of which
could be inlined, because nothing enabled cross-file optimisation. That chain is
**38.5% of runtime on the forward step** and **32.6% on SWBLI**.

| Case | Gain |
|---|---|
| SWBLI | **−10.3%** |
| Forward step (2D Euler) | **−9.6%** |
| premixed-flame (ZK) | **−18.6%** |
| Fer14 (with `ipow`) | **−10.6%** |

Combined with `ipow`, Fer14 goes 13185 → 9314 ms, **−29.4%**. `ctest` 14/14.
`Reconstruction` and `f_ss` disappear as symbols entirely; `Convective_Flux`
grows 2475 → 3530 bytes; the binary shrinks 7.07 → 5.65 MB.

**Correctness under parallelism** — 2 MPI ranks and 4 OpenMP threads both produce
output bit-identical to serial, same as the non-IPO build. Worth stating because
`ctest` has no MPI coverage at all: the `needs-mpi` label appears only in the
header comment of `test/CMakeLists.txt` and is not attached to any test.

**Numerics** — IPO changes results at round-off, growing from ~1e-14 relative at
t = 10 µs to ~3e-9 at 250 µs on Fer14, three orders below the ODE solver's own
1e-5 tolerance.

#### The order of these two matters

Before `ipow`, `-ipo` *cost* 7.4% on Fer14 while gaining ~10% everywhere else.
The reason is worth recording, because it is a trap that will recur:

- With chemistry frozen, the same binary pair flipped sign (−4.5%), so the
  regression was confined to the chemistry path.
- RHS evaluation counts were **identical** (11,728,090 in both), so it was not
  extra ODE work.
- The penalty tracked FP aggressiveness — strict −9%, `fast=1` +4%, `fast=2` +8%.

That ordering is the signature of **lost vectorisation**: `-ffast-math` vectorised
the Arrhenius loops, `-ipo` then inlined `f_kf` into them, the body grew past what
the vectoriser would accept, and it fell back to scalar. `ipow` shrinks that loop
enough that the problem disappears — with `ipow` in, full `-ipo` is the fastest
option on Fer14 by 10.6%. An earlier workaround that excluded
`Lib_Chemistry_wdot.f90` from IPO is therefore **not** needed, and was removed;
it cost 3.5% on premixed-flame.

#### Turning it off

```sh
./install.sh build --no-ipo          # or: cmake -DUSE_IPO=OFF
```

`USE_IPO` defaults to `ON` and **self-disables** rather than producing a broken
build: the flag is probed with a link step, and if either the probe or the
archiver lookup fails, the build falls back to a plain `RELEASE` with a warning.

Three things were needed to make `-ipo` work at all, all now handled in
`cmake/SetFortranFlags.cmake`:

- **The flag must be on the link line too.** IPO objects hold IR, not ELF; code
  generation happens at link time.
- **An IPO-aware archiver.** GNU `ar`/`ranlib` cannot index IR objects
  (*"archive has no index"*), and `xiar` merely forwards to `ar`. oneAPI ships
  `llvm-ar`/`llvm-ranlib` beside the compiler. Note that setting only the *cache*
  variables `CMAKE_AR`/`CMAKE_RANLIB` is not enough — `CMakeFindBinUtils` leaves
  plain variables of the same name in scope, and those shadow the cache when
  CMake expands the archive rules.
- **A probe that actually tests linking.** `SET_COMPILE_FLAG` passed candidate
  flags through `TRY_COMPILE(... COMPILE_DEFINITIONS)`, i.e. compile-only, so
  `-ipo` compiled fine, failed to link, and was silently rejected. It had
  therefore **never** been enabled on any Intel compiler, despite being listed in
  the release flags since the first commit. `SetCompileFlag.cmake` now takes a
  `LINK` keyword that also puts the flag on the test program's link line. It is
  opt-in per call site on purpose: `LINK_OPTIONS` is a list, so a multi-word flag
  like `-check bounds` would reach the linker as one quoted argument.

#### Why GNU `-flto` is not offered

`-flto` is the equivalent and is deliberately excluded. It is fragile with
mixed-language static libraries and MPI wrappers, and it is untested here: the
only GNU compiler on the build hosts is **gfortran 7.5.0** (2019), where Fortran
LTO is markedly less reliable than in modern GCC, and OSLO links the Intel-ABI
static library `lib/Intel-ODE/lib/intel64/*`. If a GNU build becomes a real
target, this needs measuring before any claim is made about it.

### 3. Kill the per-face `dl` array temporary in `Fluxes_blk`

`Convective_Flux` takes the four-cell length stencil as `dimension(-1:2)`, an
explicit-shape dummy, and it was being handed
`blk % dl(i-1:i+2,j,k) % c(1)` — a **stride-3 section of a derived type**. The
compiler cannot pass that by reference, so it materialized a contiguous
temporary on *every face*, with a runtime-length copy loop plus an
`_intel_fast_memcpy` call for the long-vector branch. Replaced by four explicit
scalar stores into a `dl4(-1:2)` subroutine local (already per-thread —
`Fluxes_blk` runs inside an orphaned parallel region, so no `private` clause is
needed).

Static evidence: `_intel_fast_memcpy` call sites in `Mod_Fluxes.s` drop
**72 → 57**, i.e. 5 per direction × 3 directions.

| Case | Before | After | min of 6 |
|---|---|---|---|
| Forward step (2D Euler) | 19.111 s | 18.252 s | **−4.5%** |
| SWBLI (100 it, dumps off) | 97.648 s | 89.987 s | **−7.8%** |

Output **bit-identical** on both cases, all files. `ctest` 14/14.

Measurement note: this was taken with the build host at load 156-172 of 192
threads, where run-to-run spread (12–20%) exceeds the effect. Six interleaved
repetitions; 10 of 12 pairs favour the new code, and the two that do not are one
per case. Minimum, pairwise win count and median ratio (-4.3% fstep, -6.3% SWBLI)
all agree on sign and rough size, and the static memcpy count corroborates. Do
not try to reproduce this from a single run.

This is the same defect class as dima-MOSE `e31af5e` ("replace 10-element array
literals in `Fluxes_blk` with explicit stores"). That commit itself does **not**
port — its literals exist only to feed the `rho_c`/`T_c` thermo cache, which
this tree does not have. The `dl` section is where the equivalent cost lives
here.

### 4. Single-species fast path for the Wilke mixing rule (FLINT)

For a single species the Wilke mixing rule is the **identity**, and the code was
computing it anyway — two `dsqrt` plus an O(ns²) loop per call, once per face.
With `ns = 1` the only molecular-weight ratio is `wm_tab(1)/wm_tab(1) = 1`, so
`Mi_Mj_pow_m025(1,1) = 1` and `inv_sqrt8_1p(1,1) = 1/sqrt(16) = 1/4` — both
exact — hence `co_fiij` returns `fi(1,1) = (1+1)²/4 = 1`, the mole fraction
`Xi(1)` is 1, `lam_den(1)` is 1, and the mixture value reduces to the
pure-species table interpolation. Guarded in `co_k_mi_lam_Wilke_expr` and
`f_laminarViscosity` (`lib/FLINT/src/lib/Lib_ThermoTransport.f90`) — FLINT
`a4eebbd`, "add ns=1 shortcut for viscous properties".

| Case | Before | After | min of 6 |
|---|---|---|---|
| SWBLI (100 it, SA) | 77.858 s | 69.490 s | **−10.7%** |
| Fer14 (control, multi-species) | 11.349 s | 11.326 s | ±0 |

All six interleaved pairs landed between −10.5% and −12.0%; paired mean −11.1%.
`ctest` 14/14.

**The control is the point.** Fer14 is multi-species, so the branch is never
taken, and its output is **bit-identical, 26/26 files** — which is what proves
the guard is doing what it claims and not quietly changing the general path.

**Numerics.** Not bit-identical on the single-species case, by construction:
`ratio = sqrt(mi)/sqrt(mi+1d-20)` and `Xi(1)` each land within a rounding of 1
rather than on it. Measured after 100 iterations: `residual-history.dat`
identical, worst per-variable relative difference in the field **4.5e-14**
(`v`; `p` 1.5e-14, `T` 1.4e-14, `rho` 2.1e-14) — round-off, ~200 ulp of double
accumulated over 300 RK stages.

Applies to any single-species viscous run (SWBLI, both flat plates). No effect
on multi-species cases, and no memory cost.

### 5. Kill the array temporaries inside the diffusive kernels

Two temporaries, both **inside** `Lib_Diffusive` rather than at a call boundary,
both paid on every viscous face:

- **`Gradient = matmul ( Gradient, M )`** (`Diffusive_Flux`). `Gradient` appears
  on both sides, so the standard requires the RHS to be evaluated in full before
  the assignment and the compiler must materialize the whole `(nprim,3)` result.
  Because `nprim` is a run-time bound it does that with library calls: **3
  `_intel_fast_memset` to zero the temporary's columns plus 3
  `_intel_fast_memcpy` to copy it back**. Replaced by an explicit loop over `v`
  that holds the three row values in registers.
- **`Gradient(nu:nw,:)`** (`Compute_Diffusive_Flux`). Three *rows* of a
  `(nprim,3)` array is a strided section, so each of the three consumers —
  `Eddy_Viscosity`, `stress_vector`, `RANS_Diffusive_Flux` — built its own
  contiguous temporary, 3 memcpy apiece. `Gradient` is `intent(in)` and is never
  modified there, so one compile-time-shaped `VelGrad(3,3)` serves all three.

The same two patterns appear in the wall BCs (`Lib_BC_Fluxes_Wall_Temperature`,
`Lib_BC_Fluxes_Wall_Heat`) and were fixed identically.

Static evidence, counted in the **linked** binary, not per-file assembly:

| Symbol | Before | After |
|---|---|---|
| `diffusive_flux` | 6 | 0 |
| `compute_diffusive_flux` | 16 | 7 |
| `bc_wall_temperature` / `bc_wall_heat` (per-file) | 20 each | 8 each |

Solver time (`Time of operation`), SWBLI/SA, 100 it, single thread, 4
interleaved repetitions on quiet sprop2:

| Arm | min | median | vs base |
|---|---|---|---|
| base | 56.79 s | 56.95 s | — |
| internal faces only | 53.03 s | 53.16 s | **−6.6%** |
| + wall BCs | 52.54 s | 52.70 s | **−7.5%** |

`base > new > new2` held in **all 4** repetitions; min and median agree to
0.03 pp. On wall clock (which carries ~11 s of fixed setup/IO) the same change
is −5.5%, from a separate 6-repetition run where all 6 pairs favoured the new
build.

Output is **bit-identical**, 4/4 files, for both stages — the matmul reordering
does not even move round-off. `ctest` 14/14.

**The control is what makes this trustworthy.** `supersonic-forward-step` is
`equations = euler`, so `model = 0` and `Diffusive_Flux` is never called. It
measured **+0.15%** (2 of 5 pairs favouring the new build, i.e. noise) against
−5.5% on the target over the same interleaved run.

**Why the earlier survey missed these.** [Open item 1](#open-items) checked the
argument-setup basic block of the `Convective_Flux`/`Diffusive_Flux` *call sites*
in `Fluxes_blk` and correctly found zero memcpy there — that conclusion still
holds, and was re-confirmed here in the `-ipo` binary (the 105 `_intel_fast_memcpy`
sites inside `fluxes_blk` are all in versioned fallback paths guarded by a
run-time stride test, not on the hot path). The temporaries were one level down,
inside the callee. **When hunting for this defect class, attribute the calls to
source lines** — compile with `-S -g` and read the `.loc` directives — rather
than counting them per file.

---

## Measured and rejected

### Static dispatch for `rlimiter` — no gain on its own

Previously listed in this file as "expected 10–20%". **Measured: 1.00x.**

A standalone kernel reproducing `Reconstruction` exactly (elemental over `nprim`,
two limiter calls per element), built with the production flags, limiter chosen
from a run-time string so it cannot be constant-folded:

| Dispatch | no `-ipo` | with `-ipo` |
|---|---|---|
| procedure pointer (current) | 1.044 s / 2.317 s | 0.511 s / 1.120 s |
| `select case` → other module | 1.046 s / 2.317 s | 0.471 s / 1.050 s |
| `select case` → bodies co-located | 0.846 s / 1.875 s | 0.471 s / 1.050 s |

*(nprim = 13 and 29.)*

The cost was never the indirect jump — it was the **lost inlining**. A direct call
into another module is exactly as slow as the pointer. Only co-locating the bodies
helps (1.24x), and `-ipo` beats that (2.2x) without touching the source.

The profile agreed the target was real: `rlimiter_vanleer` is the single hottest
function in the Euler case at **10.5%** of total runtime, above `roe_averages`,
and above `Reconstruction` itself at 8.3%. The fix was a compiler flag, not a
refactor.

### Static dispatch for `Riemann` — not worth it

Previously listed as "expected 5–10%". One indirect call per face into a routine
of several hundred flops; the middle row of the table above is exactly that
pattern and shows 1.00x. Would touch 7 call sites for nothing. `-ipo` covers it.

---

## Parallel scaling

The results themselves, and the layout rule that follows from them, are in the
user guide: [Parallel Execution](../user/parallel.md). This section records only
what a developer needs — how those numbers were obtained, what the measurement
cost to get right, and which of the remaining losses are code and which are
hardware.

Reference measurement: SWBLI 6.8 M cells / 24 blocks, restart from binary
`szplt`, `res-diter` and `shell-diter` past the run length, 40 iterations,
`timer-diter = 10` and the **last** window scored, minimum of repeats. Nodes of 4 sockets ×
20 cores (Xeon Gold 6230, no hyperthreading, 3 of 6 memory channels populated
per socket).

### Report scaling with the clock held constant

Parallel efficiency is a ratio against a single-core run, and on a turbo-capable
machine that reference runs at a completely different frequency from a full
node — here 3.9 GHz against ~2.7 GHz, i.e. a clock ratio of **0.685**. Left
uncorrected, that alone accounts for about 28 of the 40 efficiency points
apparently lost at 80 cores. **An uncorrected curve measures the machine's power
management as much as the code**, and it is not comparable with anything
published from other hardware.

Disabling turbo needs root (`/sys/devices/system/cpu/intel_pstate/no_turbo`).
The same effect is available from user space: run MOSE on `N` cores and a
**ballast** — an L1-resident spin loop, no DRAM traffic, so memory conditions
are unchanged — on the other `80−N`. The package is then fully loaded for every
point including the single-core reference, so all of them share one clock. At
full node no ballast is needed or possible, and the corrected and uncorrected
runs are the same run.

Validation and its limits:

- Per-core P-states mean the ballast's *identity* matters: an AVX-512 ballast
  downclocks itself and frees package power, shifting absolute times by 2.7–5.5 %
  against a scalar one. It largely cancels in a ratio — the 20 → 40 core
  efficiency is 97.1 % with AVX ballast and 98.6 % with scalar, against an
  8–9 point correction — but it is the dominant uncertainty, and it grows as the
  ballast fraction shrinks: **±0.1 points at 4 cores, ±1.5 at 20, ±2.6 at 40,
  ±5 at 80**. Quote the conservative edge.
- MOSE stalls on memory and draws less power than either synthetic ballast, so a
  MOSE-like neighbour would leave more headroom and make the reference *faster* —
  which pushes the quoted efficiencies down, not up.
- The all-core clock differs between nominally identical nodes by up to 5.5 %,
  while their unloaded clocks agree to 0.3 %. **Every clock-pinned curve needs
  its reference measured on the same node, in the same job.**
- `raw_eff / clock_ratio` reproduces the ballasted efficiency at every point to
  better than 0.2 points, and the 1-core ratio reproduced across three jobs to
  0.25 %.

The one subset that needs no correction at all is inter-node scaling at a fixed
per-node layout (80 / 160 / 240 cores, all at 80 cores per node): every node is
equally loaded by construction, and it comes out at **98.7–100.9 % of ideal**.

### Traps that produced wrong answers here

Each of these was believed before it was checked:

1. **Confusing NUMA with turbo.** The socket-coverage effect — spreading ranks
   over all sockets is worth 21–27 % at fixed core count — is *mostly per-socket
   power headroom* below ~40 cores, not memory locality. Clock-pinned it is
   ~6 % at 40 cores and **0.01 %** at 16 (pure OpenMP vs pure MPI, which as
   measured differ by 25 %). Above half a node it is unambiguously real. Any
   claim about a placement mechanism has to be re-measured with every socket
   loaded before it is attributed to NUMA.
2. **Comparing across binaries or node sets.** Speed-ups and efficiencies must
   come from one binary, one node set, and one reference measured in the same
   job. Mixing them produced a 25 % "difference" that was entirely protocol.
3. **Too few repetitions.** `1×80` spans 10.9 % over 8 runs, but one job alone
   produced two runs 1.1 % apart. Configurations where one rank spans sockets
   need ≥ 3 repetitions; socket-contained layouts replicate to <1.5 %.
4. **Too short a window, or residuals left on.** With `res-diter = 1`, the
   global reduction every iteration changes the shape of the curve; with a short
   run, set-up dominates. Score the last timer window of a run of at least a few
   tens of iterations.
5. **Reading clocks from `/proc/cpuinfo`.** `sort -rn | head -N` samples the
   highest-reporting cores, which at low load are idle cores boosting, and
   `intel_pstate` reports the requested P-state rather than the delivered one.
   Infer the clock from a calibrated kernel's throughput instead.

### What is left, and whose fault it is

| Loss | Where it comes from | Actionable in the code? |
|---|---|---|
| ~12 points at full node (88 % vs 100 %) | shared memory bandwidth — 3 of 6 channels populated here, STREAM saturating near 16 threads/socket — plus a residual clock effect the ballast cannot remove | **No**, on this hardware. MOSE sits *above* the machine's own streaming-efficiency curve at the knee. Expect a higher plateau where all channels are populated |
| 0–42 % extra cells at awkward rank counts | block splitting when ranks do not divide blocks | Partly — a splitter that optimised surface area rather than balance alone would help |
| Pure OpenMP collapse beyond one socket | first touch defeated by transparent huge pages at this block size | `First_Touch_Block` in `Mod_Allocate_Data` is correct and works with 4 kB pages or larger blocks; it cannot express itself when a field is a handful of 2 MB pages |
| OOM at 160+ pure-MPI ranks | per-rank allocation scales with the *total* block count | **Yes** — this is a real defect, and it is what caps pure MPI |
| Serial solution dumps | `IO_Solution` gathers the whole field to rank 0 | **Yes** — `teciompi` is already linked but unused. The largest single gap against comparable codes |

---

## Profiles

VTune user-mode sampling, single thread, self time as % of total run.

**Forward step (2D Euler)** — `Convective_Flux` subtree ≈ 42%:

| Function | % |
|---|---|
| `rlimiter_vanleer` | 10.5 |
| `roe_averages` | 9.5 |
| `Reconstruction` | 8.3 |
| `Convective_Flux` | 6.5 |
| `riemann_HLLE` | 4.2 |
| `check_gas_state` | 3.1 |

**SWBLI, pre-`-ipo`** — the diffusive path is the bigger target, not the
convective one: `co_fiij` is the hottest function overall at 6.8%, and `co_fiij`
+ `Diffusive_Flux` + `Tangential_Gradient` + `Compute_Diffusive_Flux` +
`co_k_mi_lam_Wilke_expr` together are **26.6%**, against ~18% for the convective
path.

**SWBLI, post-`-ipo`** (SA, 60 it, 57.90 s CPU) — use this one for anything you
plan to act on. `co_rotot_Rtot` no longer appears as a symbol at all: `-ipo`
inlines it, so its cost now sits inside `Diffusive_Flux` and
`Tangential_Gradient`.

| Function | % | | Function | % |
|---|---|---|---|---|
| `Convective_Flux` | 7.7 | | `Compute_Diffusive_Flux` | 4.0 |
| `for_index` — *setup/IO, not solver* | 6.5 | | `Spalart_Source_Terms` | 3.8 |
| `Diffusive_Flux` | 6.5 | | `co_k_mi_lam_Wilke_expr` | 3.6 |
| `Tangential_Gradient` | 5.8 | | `Fill_Ghost_Cell` | 3.7 |
| `Fluxes_blk` | 5.3 | | `riemann_HLLC` | 2.1 |
| `co_fiij` | 5.1 | | `f_laminarViscosity` | 1.4 |
| `rlimiter_vanleer` | 4.7 | | `Compute_dt` | 0.8 |
| `roe_averages` | 4.5 | | | |

Grouped: convective path 19.0%, the three diffusive kernels 16.3%, face
Wilke/cp (`co_fiij` + `wilke_expr` + `f_cp_expr`) 9.7%. The last of those is
what [Landed §4](#4-single-species-fast-path-for-the-wilke-mixing-rule-flint)
removed for single-species runs.

**Read this profile with the run length in mind.** 60 iterations is short enough
that setup and I/O still register: `for_index` is string parsing in the Tecplot
reader and `MOSE_Setup`, not solver work, and `write` and `Fill_Ghost_Cell`
carry some of the same. Profile a longer run before treating any entry outside
the flux kernels as a solve-loop target.

**Fer14** — chemistry dominates: `general` alone is 44%.

---

## Open items

Predictions carried over from the original (April 2026) version of this file.
**None of these have been measured** — the estimates are the original author's and
should be treated as hypotheses, given that the two predictions in this file that
*were* checked (limiter and Riemann dispatch) both came out at 1.00x.

| # | Action | Effort | Predicted |
|---|---|---|---|
| 1 | ~~Eliminate `Prim`/`Res` array temporaries in `Convective_Flux` calls~~ | — | **no such temporaries exist** |
| 2 | Precompute `rho`, `T`, `Rgas` per cell (the dima-MOSE thermo cache) | medium | **assessed, 3–4%** |
| 3 | Fuse convective + diffusive loops per direction (cache locality) | small | 5–10% |
| 4 | AoS → SoA for the metric derived types (`dl%c(3)`, `N(3)`) | large | 5–10% |

**Item 1 is closed — the temporaries it predicted do not exist.** The `dl`
section was the only one, and it is [Landed §3](#3-kill-the-per-face-dl-array-temporary-in-fluxes_blk).
Checked directly in the `-S` output: the argument-setup basic block of each of
the three `Convective_Flux` call sites, and of all five `Diffusive_Flux` sites,
contains **zero** `_intel_fast_memcpy`. The arguments are pushed as computed
addresses. `blk % P(:,i-1:i+2,j,k)` and `blk % R(:,i:i+1,j,k)` *are* contiguous —
full first dimension, stride-1 range in the second, scalar in the rest — and ifx
proves it, so no copy is made even against an explicit-shape dummy.

Two corrections behind this, both mine, both instructive:

- The 57 surviving `_intel_fast_memcpy` sites in `Mod_Fluxes.s` are real but are
  **not** argument temporaries; they sit in loop bodies and remainder/peel
  variants. Counting call sites in a file says nothing about which call they
  belong to — check the basic block.
- `for_index` at 6.5% in the SWBLI profile is **not** solver work. Every caller
  is string handling: `tec_read_structured_multiblock`, `MOSE_Setup`,
  StringiFor, the input registry. A 60-iteration profile is short enough for
  setup and I/O to dominate an entry like that — the same short-window trap this
  file warns about in §1 of *How to measure*, applied to a profile rather than a
  timing.

The only genuinely non-contiguous section left in `Fluxes_blk` is the shock
detector's `blk % P(np,i-1:i+1,j-1:j+1,k-1:k+1)` — a 3×3×3 strided gather. It is
dead at runtime in every shipped case: `shock-detector` defaults to unset, so
`SD_id = 0` and neither `select case` branch is entered. Only worth touching if
a production case turns the detector on.

**Item 2 — the dima-MOSE thermo cache (`Mod_Thermo_Cache.f90` + Option B
`6b88309`) — assessed, deferred.** The cell cache is worth −18% there, but it
drifts at O(Δx²) because it face-averages cell-computed thermo. **Option B**
removes the drift by recomputing T/cp/mil/kl at the face — and that is exactly
the 9.7% face-thermo block in the profile above, so Option B declines to touch
the biggest thing the cache could address. What is left addressable here is
`Tangential_Gradient`'s thermo half, a slice of `Diffusive_Flux`,
`f_laminarViscosity` and `Compute_dt`: **a 5–6% ceiling, 3–4% realistic**, for
~24 B/cell of ghosted storage against a memory ceiling that already caps you
near 11 ranks/node at 20M cells.

Two things the assessment turned up:

- **The transport half of the cache would be a net loss in this tree.** Under
  Option B `kl_c` has no consumer at all, `cp_c` has one (`Lib_Newstate`), and
  `mil_c`'s consumers (`Compute_dt`, Spalart) evaluate twice per step while the
  refresh would run three times. If it is ever built, build `rho_c`/`Rg_c`/`T_c`
  only.
- The single-species Wilke observation, which came out of the same profile and
  turned out to be worth more than the whole cache — [Landed §4](#4-single-species-fast-path-for-the-wilke-mixing-rule-flint).

Caveat: SWBLI is single-species, and dima measured Option B on 7-species
ONERA-7, where `co_rotot_Rtot` is O(ns) and the cacheable slice is genuinely
larger. Their −12.6% is likely real *for that case*. **There is no
multi-species viscous reference case here** — Fer14 and premixed-flame are 1D
and chemistry-dominated. That gap blocks this item and the transposed-table
question both.

**Refinement (from the §5 work): most of the redundancy is over *cell* states,
and caching those is exact.** Count the temperature evaluations on one viscous
face: `Diffusive_Flux` does 2 (`T1`, `T2`), each `Tangential_Gradient` does 4 so
8 more, and `Compute_Diffusive_Flux` does 1 — **11 per face**, each a
`co_rotot_Rtot` plus a division. Ten of the eleven are of **cell-centred**
states (`blk % P` at the owner, neighbour, or a tangential neighbour); only the
last is of the face average `0.5*(Prim1+Prim2)`. In 3D there are ~3 faces per
cell, so a cell's `rho`/`T` is recomputed on the order of **30 times per RK
stage**.

This matters for the drift objection above: it applies to a cache that
*face-averages cell-computed thermo*, which is an O(Δx²) approximation. Reading
back a per-cell `T` for the ten cell-centred evaluations is **not** an
approximation — it is the same number, so that part is bit-identical by
construction, and only the eleventh evaluation has to stay at the face.

Two things still argue for caution, and neither is resolved:

- **The payoff scales with `ns`.** At `ns = 1` `co_rotot_Rtot` collapses to
  `rho = P(1)` with a constant `Rgas`, so all that is saved is the division.
  That is why SWBLI is the wrong case to judge this on, and it is consistent
  with dima's larger number on 7 species.
- **It trades flops for bandwidth, in the direction the machine is worst at.**
  `P` is already being read for the gradients, so computing `T` from it costs no
  extra traffic; a separate `T` array adds 8 B/cell of ghosted storage *and* a
  load per use. On nodes that are already DRAM-starved at high thread counts,
  that can be a loss at exactly the thread counts that matter. Measure it
  threaded, not just single-threaded.

### Resolved

- *"Priority 0 — stray `stop` in `Mod_Fluxes.f90` disables directions 2 and 3."*
  No `stop` exists in that file; fixed at some point since April 2026.

---

## Comparison against UCNS3D

Done by reading both trees, **not** by running a common case. That was the
original plan and it was dropped on purpose:

- UCNS3D is an *unstructured* mixed-element high-order FV/DG code (Fluent `.msh`
  input, WENO/DG reconstruction over large stencils); MOSE is a *structured
  multiblock curvilinear* second-order MUSCL code reading Tecplot multiblock.
  There is no mesh both codes accept, and no scheme both codes run.
- Even after converting the SWBLI mesh to hexahedra and forcing UCNS3D to
  second order, the measurement would be "structured direct addressing vs
  unstructured indirect addressing on a structured mesh". The sign of that is
  known in advance and nothing actionable follows from the magnitude.

What the reading did establish, per subsystem:

| Subsystem | Verdict |
|---|---|
| Halo exchange | **MOSE is ahead.** MOSE uses persistent requests (`MPI_SEND_INIT`/`RECV_INIT`/`STARTALL`) with `WAITALL`, and 4 barriers in the whole tree. UCNS3D uses `mpi_sendrecv` (60 sites) and **335 `mpi_barrier`**. Consistent with MOSE's measured 97.5–100% inter-node efficiency — there is nothing to copy here. |
| Threading | **Even.** MOSE opens **one** `!$omp parallel` region around the whole Strang/RK loop in `Mod_Explicit` and uses orphaned `!$omp do` inside, so it pays no repeated fork/join. UCNS3D does the same thing. Neither code uses a single SIMD directive. |
| Solution output | **UCNS3D is ahead, and this is MOSE's real gap.** UCNS3D writes collectively with MPI-IO (`mpi_file_set_view`, `mpi_file_write_all`, `mpi_type_create_indexed_block`). MOSE gathers the entire field to rank 0 and writes it serially — `IO_Solution` allocates the full `IOfield` on the root only, by design. That serializes every dump *and* puts the whole domain on one rank. Note `teciompi` is already linked but unused. |
| Kernel specialisation | **UCNS3D is ahead in kind, MOSE in one instance.** UCNS3D keeps hand-specialised `_ideal` clones of the hot flux kernels (`calculate_fluxeshi_convective_inner_cell_ideal`, `..._diffusive_inner_cell_ideal_noturb_3d`) beside the general ones, so the common physics path carries no run-time branching. MOSE has exactly one instance of that idea — [Landed §4](#4-single-species-fast-path-for-the-wilke-mixing-rule-flint) — and it was the single largest win in this file. That is the transferable lesson. |

The concrete speed defects found in MOSE during this comparison were **not**
things UCNS3D does better; they were Fortran array-temporary bugs
([Landed §5](#5-kill-the-array-temporaries-inside-the-diffusive-kernels)), and
UCNS3D has the same class of problem in its own kernels.

---

## Also worth knowing

- **`-xHost` is resolved on the build host.** The sprop machines are AMD
  (EPYC 9474F). If binaries built there are ever run on Intel hardware, the ISA
  selection is being made against the wrong CPU.
- **The release flag list contains `-O3` twice.** Harmless; an artefact of
  `SET_COMPILE_FLAG` appending to CMake's own default `RELEASE` flags.
