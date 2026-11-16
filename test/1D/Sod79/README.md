# Sod Shock Tube

The results are discussed in the documentation, [1D Euler Problems](../../../docs/vv/1D.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[0, 1]$ |
| Final time | $t = 0.2$ |
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | MC |
| Riemann solver | SLAU |
| Boundaries | extrapolation (`400`) at both ends and on the $j$ faces |
