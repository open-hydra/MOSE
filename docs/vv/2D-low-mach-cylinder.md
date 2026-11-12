# Low-Mach Cylinder: Preconditioning

Inviscid flow past a circular cylinder at free-stream Mach numbers 0.01 and 0.001. At these speeds the flow is
incompressible to $O(M^2)$, and a compressible solver must reach that limit: its pressure field, scaled by the dynamic
pressure, must not depend on the Mach number. An upwind flux whose dissipation scales with the speed of sound misses it,
with pressure errors of order $M$ instead of $M^2$ (Guillard and Viozat, 1999). The case verifies the Weiss–Smith
preconditioning of MOSE, `integration-variables = prec` with the flux `HLLC-PC`
([Time Integration](../theory/time-integration.md#preconditioned-update),
[Riemann Solvers](../theory/riemann-solvers.md#hllc-pc)), against the same flow without it.

**References**: J. M. Weiss, W. A. Smith, *AIAA J.* 33(11), 1995; H. Guillard, C. Viozat, *Comput. Fluids* 28(1),
1999; potential flow past a circular cylinder.

<figure>
  --8<-- "docs/vv/images/LMC-cp-fields.svg"
  Pressure coefficient near the cylinder, flow from left to right: preconditioned at Mach 0.01 and 0.001, and without
  preconditioning at Mach 0.01 (its own colour scale).
</figure>

---

## Problem setup

| Parameter | Value |
|---|---|
| Cylinder radius $R$ | 0.5 m |
| Free stream | air of `test/common/Air` (cp = 1004.5 J/(kg K)), $p_\infty = 10^5$ Pa, $T_\infty = 300$ K |
| Mach number | 0.01 ($U_\infty$ = 3.47 m/s) and 0.001 |
| Inflow, upstream third of the far boundary | `401`, total conditions of the free stream, flow along $x$ |
| Outflow, rest of the far boundary | `406`, $p_\infty$ |
| Cylinder | slip wall (`300`) |
| Initial field | uniform free stream |

The exact solution of incompressible potential flow is

$$
C_p = \frac{p - p_\infty}{q_\infty} = 2\Bigl(\frac{R}{r}\Bigr)^2 \cos 2\theta - \Bigl(\frac{R}{r}\Bigr)^4,
\qquad q_\infty = p_0 - p_\infty,
$$

$\theta$ measured from the downstream axis; on the cylinder $C_p = 1 - 4\sin^2\theta$, from 1 at the stagnation points
to −3 at the shoulders.

### Grid

One block, an O-grid of 96 × 48 cells and one cell in $k$. The $i$ lines run around the cylinder, clockwise from the
downstream axis, where the two ends of the block meet in a block connection (`101`) on itself. The $j$ lines run radially
from $r = R$ to $r = 40R$, with cells growing geometrically from a first cell as wide as it is long.

<figure>
  --8<-- "docs/vv/images/LMC-grid.svg"
  The O-grid: the whole domain (every second grid line) and the region next to the cylinder.
</figure>

## Numerical setup

| Run | `integration-variables` | `riemann-solver` | Mach | `preconditioning-eps-min` |
|---|---|---|---|---|
| `prec-m2` | `prec` | `HLLC-PC` | 0.01 | 1e-3 |
| `prec-m3` | `prec` | `HLLC-PC` | 0.001 | 1e-4 |
| `cons-m2` | `cons` | `HLLC` | 0.01 | — |

All three runs are steady (`time-accurate = false`, local time steps), first order in space, with forward Euler at
CFL 0.8, for 3000 iterations. `preconditioning-eps-min` is a tenth of the Mach number, so that the floor
$\varepsilon a$ of the reference velocity stays below the flow speed, except near the stagnation points, and the two
preconditioned runs are the same problem in scaled variables.

## Results

### Pressure on the cylinder

<figure>
  --8<-- "docs/vv/images/LMC-cp-wall.svg"
  Pressure coefficient of the first row of cells against potential flow at the same cell centres. Left: the two
  preconditioned runs, on top of each other. Right: without preconditioning, against the preconditioned run.
</figure>

The preconditioned runs follow potential flow over the front of the cylinder, from the stagnation point to the
shoulders: the largest difference on the front half is 0.20, on the shoulders. Behind the shoulders they depart from it, as the wake below shows. Without
preconditioning, Cp reaches 7.5 at the front and −5.3 at the back: the pressure field of the Mach 0.01 run is off by
up to 6.5 dynamic pressures, the $O(1/M)$ error of a flux whose dissipation scales with the speed of sound.

### The incompressible limit

<figure>
  --8<-- "docs/vv/images/LMC-mach-independence.svg"
  Difference of the pressure coefficient between Mach 0.01 and Mach 0.001, preconditioned.
</figure>

Between Mach 0.01 and 0.001 the pressure coefficient changes by at most $3.0 \times 10^{-4}$ in the cells within
$1.5R$ of the axis, of the order of $M^2 = 10^{-4}$ of the faster run: the preconditioned scheme has an incompressible
limit. The difference is largest at the shoulders, where the flow is fastest.

### The wake

<figure>
  --8<-- "docs/vv/images/LMC-total-pressure.svg"
  Loss of total pressure, preconditioned run at Mach 0.01.
</figure>

Potential flow has no loss of total pressure. The dissipation of the first-order scheme on this grid makes one at the
back of the cylinder, up to 0.91 of the dynamic pressure, and the flow separates there into a wake that the Euler
equations would not hold. It sets the rear stagnation pressure, 0.3 instead of 1, and displaces the flow around the
cylinder, so the front half is compared with potential flow within a loose bound. The wake is the same at both Mach
numbers: it is an error of the discretisation, not of its low-Mach limit.

### Convergence

<figure>
  --8<-- "docs/vv/images/LMC-residuals.svg"
  Density residual of the three runs, relative to the first iteration.
</figure>

The density residual of the preconditioned runs falls by 1.8 × 10⁻³ (Mach 0.01) and 2.4 × 10⁻³ (Mach 0.001) in
3000 iterations and levels off there, held by a small unsteadiness of the wake. The curves of the two Mach numbers lie
on top of each other: in scaled variables the preconditioned runs are the same problem, and they converge alike. The
unpreconditioned run reaches a similar residual, on a solution that is wrong.

## What the test checks

`Preconditioning` (a medium test, about 37 s at one thread, `test/2D/euler/low-mach-cylinder/run.sh`):

| Check | Threshold | Measured |
|---|---|---|
| The preconditioned runs converge | residual below 1e-2 of the first | 1.8e-3, 2.4e-3 |
| Mach independence, cells with $r < 1.5R$ | $\lvert\Delta C_p\rvert$ < 2e-3 | 3.0e-4 |
| Front half of the cylinder against potential flow | < 0.3 | 0.20 |
| The unpreconditioned run differs from the preconditioned one | > 1 | 6.4 |

`InputChecks.PrecScope` checks that `prec` is refused in time-accurate runs and `HLLC-PC` without `prec`, and the
unit check `PrecTransform` the update itself, against the Weiss–Smith matrix.

## Limits

- **Second order.** With MUSCL reconstruction the density residual of the preconditioned runs stops falling at about
  a sixth of its first value, and the solution does not settle: the case is run at first order.
- **Other fluxes.** `prec` with `HLLC`, at Mach 0.01, diverges: the preconditioned update and time step meet a
  dissipation scaled with the speed of sound. Use `prec` with `HLLC-PC`.
- **`preconditioning-eps-min`.** The two Mach numbers give the same field because `eps-min` is scaled with them
  (here $M/10$). With one `eps-min` for both, the floor $\varepsilon a$ is another fraction of the flow speed at
  each Mach number, and the runs are different problems: at Mach 0.001, `eps-min` = 1e-3 lifts the reference
  velocity to the free-stream speed wherever the gas is slower.

## Reproduce

```bash
cd test/2D/euler/low-mach-cylinder
./run.sh                    # the three runs and the checks
python3 plot_vv.py          # the figures of this page and the numbers quoted
```
