# Multigrid Cycles (FAS)

The full approximation scheme (FAS) multigrid cycle described in
[Time Integration](../theory/time-integration.md#multigrid-cycle-fas) accelerates the convergence of steady
computations. With grid sequencing alone, the coarser grid levels only provide the initial field of the finer ones;
with the cycle, they keep correcting the solution on the finest level until convergence. The cycle changes the path to
convergence, not the converged solution, which must remain that of the finest grid. It is therefore worthwhile only if
it reaches that solution in less wall-clock time than grid sequencing.

The cycle is assessed on the steady test cases of MOSE in three steps:

1. **Verification** – the cycle leaves the solver unchanged when it is switched off, gives reproducible results, and
   converges to the same discrete solution as grid sequencing.
2. **Performance** – the wall-clock time to convergence with the cycle and with grid sequencing, on thirteen steady
   cases ranging from inviscid nozzle flows to a turbulent shock–boundary-layer interaction.
3. **Settings** – how the default parameters of the cycle were chosen, and guidelines for adjusting them.

**Reference**: A. Jameson, "Solution of the Euler equations for two-dimensional transonic flow by a multigrid method,"
*Appl. Math. Comput.* 13 (1983), the scheme on which the cycle is based.

---

## Verification

| Property | Result |
|---|---|
| `cycle = none` with the other keys of the cycle set | output files identical bit for bit to those of the same run without them |
| Different numbers of MPI processes and OpenMP threads | identical bit for bit |
| Restart in the middle of a run with the cycle | identical bit for bit to the run without interruption |
| Cycle, or turbulence in it, started at a given iteration (`cycle-start`, `cycle-turbulence-start`) | identical bit for bit to switching it on at that iteration with a restart |
| Settings with which the cycle cannot work | the run stops at set-up, or when `input.ini` is re-read, with a message |
| Laminar flat plate, cycle at its defaults and grid sequencing (residual threshold 10⁻¹⁰, at most 100 000 iterations of the finest level) | both skin frictions within 2.5 % RMS of the Blasius solution and 0.05 % RMS from each other (at most 0.5 % required); a density residual of 10⁻⁷ after 14 394 iterations of the finest level with the cycle, 65 554 with grid sequencing (at most half required) |

The cycle converges to the same discrete solution as grid sequencing: run to convergence, the two methods give the
same final fields.

The checks and their acceptance criteria are listed in the
[test plan](https://github.com/open-hydra/MOSE/blob/main/test/TESTING-PLAN.md) (`TimeIntegration.Multigrid*`,
`FlatPlateLaminar.FAS`, and the cases with the cycle in `MPIBitIdentity` and `Restart.Counters`); the setup of the
laminar flat plate with the cycle, and how to run it, are in
[`test/2D/viscous/flat-plate-laminar-fas/README.md`](https://github.com/open-hydra/MOSE/blob/main/test/2D/viscous/flat-plate-laminar-fas/README.md).

---

## Performance

### Method

Each case is run from the same input, initial field, grid levels and iterations on the coarser levels, once with grid
sequencing and once with the cycle. Implicit residual smoothing (IRS) is a setting of the solver, not of the cycle,
and it changes the convergence of grid sequencing as well: where the input of the test uses it, grid sequencing is run
both with and without it, and the cycle is compared with the faster of the two, so that the reference is the fastest
grid sequencing and not one slowed down by IRS. Both runs are timed from their launch, coarser levels included, until
the residual on the finest level reaches

$$R^* = 2\max\left(R_\text{min}^\text{grid seq.},\, R_\text{min}^\text{FAS}\right),$$

twice the higher of their two final residuals, so that both runs reach it. The **speed-up** is the ratio of the two
times, $t_\text{grid seq.}/t_\text{FAS}$: above 1, the cycle is faster. A run with the cycle **stalls** when its
residual levels off above the final residual of grid sequencing, and **fails** when it stops on a non-physical state.

All runs were made on one workstation (Intel Core i9-13900K, 8 OpenMP threads, one run at a time). The time of the
same run varies by 5–10 % from one day to another.

### Test cases

| Case | Physics | Cells | Grid levels | Scheme, CFL | IRS β |
|---|---|---:|:---:|---|---|
| Shock–boundary-layer interaction, SA (`2D/viscous/swbli/SA`) | RANS SA, Mach 5 | 237 440 | 4 | RK3, 1.0 | 0.5 |
| Shock–boundary-layer interaction, SST (`2D/viscous/swbli/SST`) | RANS SST, Mach 5 | 237 440 | 4 | RK3, 0.5 | 0 |
| Laminar flat plate | laminar NS | 4 096 | 3 | RK2, 0.5 | 0 |
| Turbulent flat plate, SA | RANS SA | 13 056 | 3 | RK3, 0.8 | 0.5 |
| Turbulent flat plate, SST | RANS SST | 3 264 | 3 | RK3, 0.8 | 0.5 |
| Turbulent flat plate, Wilcox 2006 | RANS k–ω | 3 264 | 3 | RK3, 1.0 | 0.5 |
| Planar nozzle, SA, O grid level 2 | RANS SA | 55 296 | 3 | RK3, 0.8 | 0 |
| Manufactured solution 2D (`mms`) | laminar NS + source | 4 096 | 4 | RK3, 0.2 | 0 |
| Chang nozzle | Euler | 7 600 | 2 | RK3, 0.8 | 0.5 |
| Rocket nozzle, frozen | Euler | 4 800 | 3 | RK3, 0.9 | 0.3 |
| Prandtl–Meyer expansion | Euler | 2 500 | 2 | RK3, 1.0 | 0.5 |
| Pressure-centrifugal 3D | Euler, rotating frame | 4 000 | 2 | RK3, 2.0 | 0.5 |
| Hypersonic cylinder, Mach 8.1 | Euler | 6 400 | 3 | RK3, 0.5 | 0 |

The inputs are those of the tests, with the largest number of grid levels that the grid allows, up to four. The
iterations on the coarser levels are the same for both methods: those of the test where it defines them
(shock–boundary-layer interaction), and otherwise 2000 per level, a fixed value of this assessment, not tuned. The other
steady tests are not included: their grids cannot be coarsened (odd or too few cells), their levels are not nested,
their coarser levels do not run on their own, or they are not run to convergence.

### Results with the default settings

<figure>
  --8<-- "docs/vv/images/FAS-summary.svg"
  <figcaption>Speed-up of the cycle with the default settings over the faster grid-sequencing run (shock–boundary-layer
  interaction: time to a skin friction within 0.1 % of the converged solution). Red: the cycle is faster; blue: it is
  slower. Circles: the scheme of Jameson, where it converges. Right: cost of an iteration on the finest level.</figcaption>
</figure>

The largest gains are obtained on the shock–boundary-layer interaction, with a speed-up of 8.1 with Spalart–Allmaras
and 3.3 with SST ([below](#shockboundary-layer-interaction)). They are followed by the planar nozzle, 2.3 (38 min
instead of 87 min to the final residual of grid sequencing), the manufactured solution, 2.2, and the laminar flat
plate, 1.9. On the turbulent flat plates the gain is smaller: 1.22 with Spalart–Allmaras and 1.12 with SST.

The cycle is slower than grid sequencing on the flat plate with the Wilcox 2006 model (0.84) and on the small inviscid
cases, which converge within seconds to minutes: Prandtl–Meyer expansion 0.75, rocket nozzle 0.61,
pressure-centrifugal case 0.25 and Chang nozzle 0.19. On these cases the iterations saved do not make up for the
higher cost of each iteration ([Cost of an iteration](#cost-of-an-iteration)). On the hypersonic cylinder the
cycle fails with every setting tested ([Known limitations](#known-limitations)).

### Shock–boundary-layer interaction

The Mach 5 interaction of the [SWBLI page](2D-swbli.md) is the largest case of the assessment, with 237 440 cells and
four grid levels. Here the density residual of both methods levels off at about 4×10⁻³, so convergence is measured on
the skin friction instead: the time after which it stays within 0.1 % of the converged solution. Grid sequencing does
not reach this tolerance within its runs, and its time is extrapolated from its last 8000 iterations.

<figure>
  --8<-- "docs/vv/images/FAS-swbli_conv.svg"
  <figcaption>Shock–boundary-layer interaction: largest difference of the skin friction from the converged solution
  against wall-clock time. Dotted: extrapolation of grid sequencing; dots: the cycle with the default settings within
  0.1 %. (a) Spalart–Allmaras; (b) SST.</figcaption>
</figure>

With Spalart–Allmaras the cycle is within 0.1 % after 0.60 h. At the end of its run, after 2.7 h, grid sequencing
still differs by 2 % in the separation region, and the extrapolation gives 4.8 h: a speed-up of 8.1. With SST the
cycle is within 0.1 % after 3.27 h, while grid sequencing still differs by 1.5 % after 6.35 h and extrapolates to
10.7 h: a speed-up of 3.3.

<figure>
  --8<-- "docs/vv/images/FAS-swbli_bubble.svg"
  <figcaption>Shock–boundary-layer interaction, Spalart–Allmaras: skin friction across the separation bubble at about
  the same wall-clock time. The cycle after 0.60 h lies on the converged solution (dashed); grid sequencing after
  0.59 h places the separation 1 mm further downstream.</figcaption>
</figure>

### Cost of an iteration

With the default settings, an iteration on the finest level costs 1.9 to 2.4 times as much with the cycle as without
it on the larger cases (planar nozzle, shock–boundary-layer interaction), and up to 3.1 times on grids of a few
thousand cells with IRS, whose smoothing is repeated at every stage on the coarser levels. The cycle reduces the
wall-clock time only where it reduces the number of iterations by a larger factor.

---

## Settings

### Choice of the defaults

Each parameter of the cycle was first varied alone, starting from the original scheme of Jameson (`pre-iter = 1`,
`damp-prolongation = 1`), and then in combination with the others. No setting is the fastest everywhere: the scheme of
Jameson is the fastest on the laminar flat plate, but stalls on the planar nozzle and fails on the shock–boundary-layer
interaction with SST; two pre-smoothing iterations (`pre-iter = 2`) are the fastest on the planar nozzle, but stall on
the SST flat plate. Two pre-smoothing iterations with the prolongation damped by one half (`pre-iter = 2`,
`damp-prolongation = 0.5`) is the only setting that converges on every viscous case where the cycle is faster than grid
sequencing, and the fastest of those that do: it is the default.

The table gives, for each case, the speed-up with the default settings and with the fastest setting found.

| Case | `irs-beta` of the test | Defaults | Fastest setting | Its speed-up |
|---|:---:|:---:|---|:---:|
| Shock–boundary-layer interaction, SA / SST | 0.5 / 0 | 8.1 / 3.3 | the defaults (SA: Jameson and `damp-prolongation = 1` oscillate; SST: Jameson stops) | 8.1 / 3.3 |
| Laminar flat plate | 0 | 1.88 | Jameson | 4.75 |
| Manufactured solution | 0 | 2.16 | Jameson with `cycle = v` | 2.87 |
| Planar nozzle, SA | 0 | 2.27 | `pre-iter = 2`, `damp-prolongation = 1` | 2.59 |
| Turbulent flat plate, SA | 0.5 | 1.22 | Jameson with `cycle-turbulence = .true.` | 3.14 |
| Turbulent flat plate, SST | 0.5 | 1.12 | Jameson | 1.60 |
| Turbulent flat plate, Wilcox 2006 | 0.5 | 0.84 | Jameson with `irs-beta = 0` | 1.45 |
| Prandtl–Meyer expansion | 0.5 | 0.75 | Jameson with `irs-beta = 0` | 1.47 |
| Rocket nozzle | 0.3 | 0.61 | the defaults with `irs-beta = 0` | 1.26 |
| Chang nozzle | 0.5 | 0.19 | the defaults with `irs-beta = 0` | 0.74 |
| Pressure-centrifugal 3D | 0.5 | 0.25 | the defaults with `irs-beta = 0` and `cfl = 0.9` (2.0 in the test) | 0.79 |
| Hypersonic cylinder | 0 | fails | none: every setting fails | — |

*Defaults*: `pre-iter = 2`, `damp-prolongation = 0.5`. *Jameson*: `pre-iter = 1`, `damp-prolongation = 1`. IRS is not a
parameter of the cycle: unless the setting states otherwise, each run uses the `irs-beta` of the test input.

### Guidelines

The guidelines below summarise what was observed on these cases. They are not general rules: the best setting depends
on the case.

1. **Start from the defaults.** Where they gain little, change one parameter at a time, and return to the defaults if
   the residual levels off above that of grid sequencing.
2. **Prolongation damping.** Less damping (`damp-prolongation = 1`) can be faster, but stalls or stops the run more
   often.
3. **Turbulence in the cycle.** With Spalart–Allmaras, try `cycle-turbulence = .true.`; with the k–ω models (SST,
   Wilcox 2006), keep it off.
4. **Implicit residual smoothing.** In general, avoid it with the cycle (`irs-beta = 0`).
5. **Short runs and strong shocks.** On cases that converge within seconds or a few minutes the cycle gains little or
   nothing, and strong shocks can stop it: use `cycle = none`.
6. **Coarse-level iterations.** Keep `coarse-iter = 1`: more iterations on the coarser levels make each cycle more
   expensive and, in general, do not make it converge faster.
7. **Delayed start.** When a level starts far from its solution, start the cycle, or the turbulence in it, after the
   initial transient of the residual (`cycle-start`, `cycle-turbulence-start`).

---

## Known limitations

- **Strong shocks.** On the hypersonic cylinder at Mach 8.1, the cycle stops on a non-physical state at a face next to
  the bow shock, in the first iterations after the prolongation to a finer level, where grid sequencing runs without
  problems.
- **Turbulence in the cycle with the k–ω models.** With SST and Wilcox 2006, `cycle-turbulence = .true.` stops the run
  or brings no gain, also when it is switched on later.
- **Steady, non-reacting flows only.** The cycle stops the run at set-up with `time-accurate = true` or with
  chemistry; the other settings with which it cannot work stop the run too
  ([When the run stops](../user/input.md#when-the-run-stops), `TimeIntegration.MultigridChecks`).
