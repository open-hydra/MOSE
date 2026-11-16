# Unity-Lewis Diffusion

The results are discussed in the documentation, [Unity-Lewis Diffusion](../../../docs/vv/1D-diffusion.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier–Stokes (laminar) |
| Domain | $[0,\,0.05]$ m |
| Grid | $200 \times 1$ |
| Final time | $200$ ms (snapshots every $50$ ms) |
| Time scheme | RK2 (time-accurate) |
| CFL / VNN | 0.95 / 0.95 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |
| Laminar Schmidt $Sc$ | 1.0 |
| Laminar Prandtl $Pr_l$ | 1.0 |

The two ends in $x$ are **periodic** (`201`) and the $j$ faces use symmetry (`300`),
so the only error sources are the diffusive flux and the time integration.
