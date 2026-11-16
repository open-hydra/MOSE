# Low-Mach Cylinder: Preconditioning

The results are discussed in the documentation, [Low-Mach Cylinder: Preconditioning](../../../../docs/vv/2D-low-mach-cylinder.md); this file holds the setup of the case and how to run it.

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

## Reproduce

```bash
cd test/2D/euler/low-mach-cylinder
./run.sh                    # the three runs and the checks
python3 plot_vv.py          # the figures of this page and the numbers quoted
```
