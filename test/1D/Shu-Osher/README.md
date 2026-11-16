# Shu–Osher Shock Entropy-Wave Interaction

The results are discussed in the documentation, [1D Euler Problems](../../../docs/vv/1D.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[-5, 5]$ |
| Final time | $t = 1.8$ |
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |
| Boundaries | extrapolation (`400`) at both ends and on the $j$ faces |
