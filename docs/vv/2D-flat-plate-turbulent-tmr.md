# Turbulent Flat Plate of the NASA TMR

The zero-pressure-gradient flat plate of the NASA Turbulence Modeling Resource (TMR) is a verification
case: on a family of nested grids, the results of a correct implementation of a turbulence model should
approach, as the grid is refined, the same answer as the codes of the TMR, for the same flow and the same
boundary conditions. This page runs the case with the Spalart–Allmaras (SA), SST and Wilcox (2006) k-ω
models on three grids of the TMR family and compares each grid, and the limit of the three grids, with
the CFL3D and FUN3D results of the TMR. It is a verification of the implementation of the three models in
MOSE as the TMR intends it, not a validation against experiments: demonstrated for SA, observed for SST and
Wilcox 2006 (see below). The [turbulent flat plate](2D-flat-plate-turbulent.md) page runs the same flow on
two of these grids with boundary conditions of its own.

**Reference**: NASA Turbulence Modeling Resource, 2D Zero Pressure Gradient Flat Plate Verification Case:
[introduction](https://tmbwg.github.io/turbmodels/flatplate.html), [grids](https://tmbwg.github.io/turbmodels/flatplate_grids.html),
results with [SA](https://tmbwg.github.io/turbmodels/flatplate_sa.html), [SST-Vm](https://tmbwg.github.io/turbmodels/flatplate_sst.html)
and [Wilcox2006-klim-m](https://tmbwg.github.io/turbmodels/flatplate_w06.html). The grids (`MESH/`) and the data
of CFL3D and FUN3D (`reference/tmr/`) are the files of the repository
[TMBWG/turbmodels](https://github.com/TMBWG/turbmodels) at commit `20a39f5`, as published.

---

## Problem setup

A plate from $x = 0$ to $x = 2$ in a domain from $x = -0.333$ to $x = 2$ and from $y = 0$ to $y = 1$, at
$M = 0.2$ and $Re = 5\times10^6$ per unit length (the TMR). Here $T_\infty = 555$ K and $p_\infty = 10^5$ Pa,
with a constant viscosity of $1.1858686\times10^{-5}$ Pa·s, which gives $Re = 5.0\times10^6$ per metre.

| Boundary | TMR | This case |
|---|---|---|
| Inflow, $x = -0.333$ | total pressure and temperature, $p_0/p_\infty = 1.02828$, $T_0/T_\infty = 1.008$, one quantity from the interior | `401` at the same $p_0$ and $T_0$, flow normal to the face, the Riemann invariant from the interior |
| Outflow, $x = 2$ | static pressure $p_\infty$, the other quantities from the interior | `406` at $p_\infty$ |
| Top, $y = 1$ | far-field Riemann condition | `406` at $p_\infty$; where gas would enter, the face acts as a symmetry face |
| $y = 0$, $x < 0$ | symmetry | `300` |
| Plate, $0 \le x \le 2$ | adiabatic no-slip wall | `301`, $q = 0$ |

The free-stream turbulence is that of the TMR runs: $\tilde\nu_\infty = 3\,\nu_\infty$ for SA,
$k_\infty = 9\times10^{-9}\,a_\infty^2$ and $\omega_\infty = 10^{-6}\,\rho_\infty a_\infty^2/\mu_\infty$ for SST and
Wilcox 2006, here $\tilde\nu_\infty = 5.6667\times10^{-5}$ m² s⁻¹, $k_\infty = 2.00699\times10^{-3}$ m² s⁻² and
$\omega_\infty = 11\,806$ s⁻¹ (computed with $R = 287$ J/(kg·K); with the gas constant of MOSE,
286.9984 J/(kg·K), $\tilde\nu_\infty = 3.00002\,\nu_\infty$). The wall value of $\omega$ is the practical one,
$\omega = 10 \cdot 6\nu/(\beta_1 d^2) = 800\,\nu/d^2$ with $d$ the distance of the centre of the wall cell
(`omega-wall-bc = practical`), the value of the original SST reference that the TMR gives; for Wilcox 2006
the TMR notes that many codes use it. The initial field is uniform: the free stream with
$\tilde\nu = 10^{-6}$ m² s⁻¹ (SA) or $k = 2.25\times10^{-3}$ m² s⁻² and $\omega = 6250$ s⁻¹ (SST, Wilcox 2006).

**Gas.** Air as a perfect gas, $\gamma = 1.4$, $\mathrm{Pr} = 0.72$ (`prandtl`), $\mathrm{Pr}_t = 0.9$ (`Prt`):
`test/common/Air-viscous`, the files that ATLAS GPB writes at ATLAS commit `ebeaf4f` for an ideal gas `Air`
with `gamma = 1.4`, `R = 287.0`, `mil = 1.1858685985e-5` and `Pr = 0.72`.

**Differences from the TMR, with their effect.**

- *Top boundary.* MOSE has no far-field condition with a moving free stream (`407` is a reservoir at
  rest), so the top boundary is an outflow at $p_\infty$. Outside the boundary layer ($y = 0.1$) the
  velocity of MOSE on the 137 × 97 grid is $u/U_\infty$ = 0.99985 to 0.99987 at $x = 0.97008$ and 1.00003
  to 1.00012 at $x = 1.90334$ for the three models, where CFL3D on the 545 × 385 grid has 0.99852 and
  0.99996 (TMR, SA). The effect on $c_f$ is not measured; as an order of magnitude, $c_f$ referred to
  $U_\infty$ scales about as $u_e^2$, so the difference of 0.13 % in $u_e$ at $x = 0.97$ is worth about
  0.3 % of $c_f$ (an estimate, not a measurement). The TMR notes that variations of the boundary
  conditions at the inflow, top and outflow may work and give similar results, and that moving the top
  boundary down to $y = 0.48$ changed $c_f$ at $x = 0.97$ and the drag by less than 0.2 %.
- *Viscosity.* Constant, where the TMR uses Sutherland's law: with Sutherland's law of the TMR, a run of
  SA on the 69 × 49 grid gives a $c_f(0.97)$ higher by 0.09 %. With a constant viscosity the temperature
  scale does not enter the nondimensional flow, so $T_\infty = 555$ K instead of 540 °R does not change it.
- *Model variants.* The TMR data are SA, SST-Vm and Wilcox2006-klim-m; MOSE runs
    - SA without the $f_{t2}$ term (SA-noft2), which the TMR says gives essentially identical results for
      this case;
    - SST-2003m (the strain rate in $\mu_t$ and in $P = \mu_t S^2$, the limiter $10\beta^*\rho\omega k$,
      $\gamma_1 = 5/9$ and $\gamma_2 = 0.44$), where SST-Vm takes the vorticity, the limiter 20 on the
      production of $k$ alone and the constants of the original model; the TMR notes that $\gamma_1 = 5/9$
      is about 0.43 % higher than the original constant and that SST-V and SST give nearly identical
      results for this flat plate. The effect is not measured;
    - Wilcox2006-klim-e (the exact production $\tau_{ij}\,\partial u_i/\partial x_j$), where the TMR codes
      approximate it with $2\mu_t S_{ij}S_{ij}$, which the TMR calls very accurate for this flow: on the
      columns of cells nearest $x$ = 0.5, 0.97 and 1.9 of the final field of the 137 × 97 grid, where the
      production exceeds 1 % of its largest value, the two differ by at most $2.7\times10^{-4}$ of it.

    The margin $m$ of the criteria below is for these differences.

## Grids

Three grids of the TMR family `flatplate_clust2`, each every other node of the next finer one:
35 × 25, 69 × 49 and 137 × 97 nodes (816, 3264 and 13 056 cells), with 29, 57 and 113 nodes on the plate,
four, three and two levels below the finest grid of the TMR, 545 × 385. The 273 × 193 and 545 × 385
grids are not run.

## Numerical setup

| Parameter | SA | SST | Wilcox 2006 |
|---|---|---|---|
| Riemann solver | AUSM+ | SLAU | HLLC |
| CFL / VNN | 0.8 / 0.3 | 0.8 / 0.3 | 1.0 / 1.0 |
| Time scheme | RK3, local time step (steady) | RK3, local time step (steady) | RK3, local time step (steady) |
| IRS $\beta$ | 0.5 | 0.5 | 0.5 |
| Space reconstruction | MUSCL, van Leer | MUSCL, van Leer | MUSCL, van Leer |
| Iterations | 300 000 | 300 000 | 300 000 |

The point-implicit factor of the turbulence equations is on (the default). Each run writes a numbered
solution every 20 000 iterations. The nine runs were made with ifx 2024.1 (RELEASE build with IPO, OpenMP
and Intel MPI, one rank) with 2, 4 and 8 threads on the 35 × 25, 69 × 49 and 137 × 97 grids; the build
of the CI (gfortran, `TESTING`), restarted from the three converged 35 × 25 solutions, stays within
$1.3\times10^{-7}$ of their values in 5000 iterations (relative; the short form below).

## Criterion

The quantities are $c_f = |\tau_x|/q_\infty$ with $q_\infty = \gamma p_\infty M^2/2$
= 2800 Pa at the face centres of the plate; $c_f$ at $x = 0.970084071$ and $x = 1.90334$, by linear
interpolation between the two nearest face centres; and the drag coefficient $C_D = \sum c_f\,\Delta x / 2$
over the plate. The iterative error $e_{it}$ of a quantity is the largest minus the smallest value over the
last three numbered solutions and the final one.

On the three grids (1 the finest, 3 the coarsest, ratio 2), as in the TMR summary of the procedure of
Celik et al. (2008): $\varepsilon_{21} = f_2 - f_1$, $\varepsilon_{32} = f_3 - f_2$, the observed order
$p = \ln|\varepsilon_{32}/\varepsilon_{21}|/\ln 2$ (with a constant ratio the term $q(p)$ of the TMR
procedure is zero), the limit $f_{ext} = (2^p f_1 - f_2)/(2^p - 1)$, $GCI = 1.25\,|\varepsilon_{21}/f_1|/(2^p - 1)$
and the uncertainty $U = GCI\,|f_1| + e_{it}$, with $e_{it}$ the largest of the three grids. The formal
order 2 and the factor 3 replace $p$ and 1.25 for oscillatory or diverging convergence, for $p < 0.5$ or
$p > 4$, or when a difference of the grids is below $10\,e_{it}$; none of these occurs here. The TMR limits
its $GCI$ further for $p < 0.95$ and $p > 3.05$; with $p$ between 1.21 and 2.02 for the three quantities of
the three models, those limits do not apply.

- **(d)** $e_{it} \le GCI\,|f_1|/10$, for $c_f(0.97)$, $c_f(1.90)$ and $C_D$;
- **(1)** on each grid, $|c_f(0.97) - (c_\mathrm{CFL3D} + c_\mathrm{FUN3D})/2| \le |c_\mathrm{CFL3D} - c_\mathrm{FUN3D}|/2 + m\,c_f$,
  against CFL3D and FUN3D on the same grid;
- **(2)** $|f_{ext} - c_\mathrm{CFL3D,\,545\times385}| \le U + m\,c_{f,1}$, for $c_f(0.97)$,

with $m$ = 0.5 % for SA, whose variant the TMR calls essentially identical to that of its data, and 1 % for
SST and Wilcox 2006, whose variants differ from those of the TMR data with an effect not measured (see the
differences above). The case passes when, for each model, (d) holds for the three quantities, (1) holds on
the 137 × 97 grid and (2) holds; (1) on the two coarser grids is given for comparison and is not required.

## Results

$c_f(x = 0.970084071) \times 10^3$ on each grid, criterion (1): the distance of $c_f$ from the mean of CFL3D and
FUN3D against the band $|c_\mathrm{CFL3D} - c_\mathrm{FUN3D}|/2 + m\,c_f$.

| Model, grid | MOSE | CFL3D | FUN3D | $\lvert c_f - \text{mean}\rvert$, band | (1) |
|---|---|---|---|---|---|
| SA, 35 × 25 | 2.58510 | 2.79569 | 2.66738 | 0.1464 > 0.0771 | outside |
| SA, 69 × 49 | 2.67302 | 2.72835 | 2.69498 | 0.0386 > 0.0301 | outside |
| SA, 137 × 97 | 2.69807 | 2.71115 | 2.70215 | 0.0086 ≤ 0.0180 | inside |
| SST, 35 × 25 | 2.37528 | 2.55183 | 2.51562 | 0.1584 > 0.0419 | outside |
| SST, 69 × 49 | 2.55977 | 2.62625 | 2.60951 | 0.0581 > 0.0340 | outside |
| SST, 137 × 97 | 2.63669 | 2.66477 | 2.65845 | 0.0249 ≤ 0.0295 | inside |
| Wilcox 2006, 35 × 25 | 2.42178 | 2.56221 | 2.56755 | 0.1431 > 0.0269 | outside |
| Wilcox 2006, 69 × 49 | 2.60919 | 2.65134 | 2.66962 | 0.0513 > 0.0352 | outside |
| Wilcox 2006, 137 × 97 | 2.67867 | 2.69148 | 2.69406 | 0.0141 ≤ 0.0281 | inside |

Limit of the three grids, criterion (2), $c_f(0.97) \times 10^3$:

| | SA | SST | Wilcox 2006 |
|---|---|---|---|
| observed order $p$ | 1.81 | 1.26 | 1.43 |
| limit $f_{ext}$ | 2.70806 | 2.69170 | 2.71960 |
| $GCI$ | 0.46 % | 2.61 % | 1.91 % |
| $U$ | 0.0126 | 0.0688 | 0.0512 |
| CFL3D, 545 × 385 | 2.70562 | 2.69085 | 2.71510 |
| FUN3D, 545 × 385 | 2.70540 | 2.69055 | 2.71770 |
| $f_{ext}$ − CFL3D | +0.090 % | +0.031 % | +0.166 % |
| $\lvert\Delta\rvert / (U + m\,c_f)$ | 0.09 | 0.009 | 0.06 |
| (2) | met | met | met |

On the same three grids the observed orders of CFL3D and FUN3D (from the TMR data) are 1.97 and 1.94
with SA, 0.95 and 0.94 with SST-Vm, 1.15 and 2.06 with Wilcox2006-klim-m. $c_f(1.90)$ and $C_D$ converge
monotonically too, with $p$ from 1.21 to 2.02. Criterion (d) holds for the three quantities of the three
models: the largest $e_{it}$ is $6.2\times10^{-5}$ of the value ($c_f(1.90)$ of SA on the 137 × 97 grid).
The density residual falls from its largest value by 2.9, 4.4 and 5.2 orders (SA), 3.7, 6.6 and 5.2
(SST), 4.4, 6.5 and 5.9 (Wilcox 2006) on the three grids; on the 35 × 25 grid with SA it stays between
$3.9\times10^{-7}$ and $5.1\times10^{-7}$ over the last 2000 iterations, while $c_f(0.97)$ moves by
$2\times10^{-7}$ of its value over the last 40 000. (For CFL3D the TMR reports the density residual typically
driven below $10^{-13}$.)

Left: $c_f(0.97)$ against $h = N^{-1/2}$ ($N$ the cells), MOSE on the three grids with the limit and its
uncertainty $U$, CFL3D and FUN3D on the five grids of the TMR. Right: $c_f$ along the plate, MOSE on the
three grids, CFL3D and FUN3D on the 545 × 385 grid.

<figure>
  {% include "vv/images/flat-plate-tmr-sa.svg" %}
  Spalart–Allmaras
</figure>

<figure>
  {% include "vv/images/flat-plate-tmr-sst.svg" %}
  SST (TMR data: SST-Vm)
</figure>

<figure>
  {% include "vv/images/flat-plate-tmr-wilcox2006.svg" %}
  Wilcox 2006 (TMR data: Wilcox2006-klim-m)
</figure>

## What is demonstrated and what is observed

- **SA: demonstrated.** The limit of the three grids lies within $U + m\,c_f$ of CFL3D on the 545 × 385
  grid, at 0.09 of that bound, with an observed order of 1.81, near the formal 2, and a $GCI$ of 0.46 %.
- **SST and Wilcox 2006: observed.** Criterion (2) is met, but with observed orders of 1.26 and 1.43,
  below 2, and a $GCI$ of 2.6 % and 1.9 %: these three grids are not in the asymptotic range for these
  models, as for CFL3D and FUN3D with SST-Vm (0.95 and 0.94) and CFL3D with Wilcox2006-klim-m (1.15) on the
  same grids. The criterion passes because $U$ is large
  ($|\Delta|/(U + m\,c_f)$ = 0.009 and 0.06), not because the limit is precise.
- **Each grid.** Criterion (1) holds on the 137 × 97 grid for the three models and not on the two coarser
  grids, where CFL3D and FUN3D also differ from each other (by 4.7 % with SA on the 35 × 25 grid).

## Limits

- Three grids: the 273 × 193 and 545 × 385 grids of the TMR are not run with this definition.
- The seven wall spacings of the study of the effect of $y^+$ of the TMR are not run.
- Constant viscosity and the top boundary at $p_\infty$ (see the differences above).
- The variants of the models differ from those of the TMR data (see above); their effect on $c_f$ is not
  measured (for Wilcox 2006 the difference of the production is).
- On the 35 × 25 grid with SA the density residual stops falling at about $4\times10^{-7}$; $c_f$ does not
  move (see the results).

## Running the case

- **Short form**, CTest `FlatPlateTMR` (validation tier, about 71 s of CPU at one thread): `./MOSE.sh testfast`
  restarts each model from its converged 35 × 25 solution (`reference/field-35x25-<model>.tec`, the final
  solution of the runs of this page) for 5000 iterations. `verify.py` checks that each restart did its
  iterations (the global iteration count of its final solution and the rows of `residual-history.dat`) and
  that $c_f(0.97)$, $c_f(1.90)$ and $C_D$ stay within $\tau = 10^{-5}$ (relative) of the converged values; then
  it computes again the tables and the criteria of this page from the stored results of the full form
  (`reference/levels.json` and the nine final wall outputs `reference/wall-<model>-<grid>.tec`), writes the
  figures in `RUN/figures/` and fails if (d), (1) on the 137 × 97 grid or (2) fails.
- **What the short form checks.** It guards the steady state of this page on the 35 × 25 grid: a change of
  MOSE that moves that steady state starts the restart away from its own steady state, and the test fails
  when $c_f(0.97)$, $c_f(1.90)$ or $C_D$ moves by more than $\tau$ in the 5000 iterations; a change that
  moves them by less in that time goes unseen. With the asymptotic wall value of $\omega$ in place of the
  practical one (SST), $c_f(0.97)$ moves by $960\,\tau$ in 500 iterations and by $9500\,\tau$ in 5000; with
  the top boundary closed (`300` in place of `406`, SA), by $30\,\tau$ in 500 and $46\,\tau$ in 5000
  ($c_f(1.90)$ by $180\,\tau$). With the CFL number halved, which changes the path and not the steady
  state, and in the builds and thread counts of the tests, after the 5000 iterations $c_f$ and $C_D$ differ
  from the converged values by at most $0.013\,\tau$. It does not check the order or the limit of the grids:
  the full form does.
- **Full form**, CTest `FlatPlateTMRFull` (label `verification`, out of the PR CI): `./MOSE.sh test` runs the
  nine cases from the uniform field, each in its folder `RUN/<model>-<grid>`; `python3 -B verify.py` gives
  the tables, the criteria and the figures, and writes in `RUN/reference/` the files of `reference/` that
  the runs give. At one thread it takes about 24.5 h (an estimate); on a cluster node
  `./MOSE.sh -p <threads> test` takes a few hours. The runs write about 1.1 GB in `RUN/` (each its 16 field,
  diagnostic and wall solutions and a `residual-history.dat` of 45 to 51 MB), which git ignores.
- When a change of MOSE moves the steady state, run the full form, copy `RUN/reference/` to `reference/`,
  and take the tables and the figures of this page from the output of `verify.py`.

```bash
cd test/2D/viscous/flat-plate-turbulent-tmr
./MOSE.sh testfast                # short form
python3 -B verify.py
./MOSE.sh -p <threads> test       # full form (long)
python3 -B verify.py
```
