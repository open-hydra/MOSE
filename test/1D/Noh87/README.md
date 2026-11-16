# Noh Implosion

The results are discussed in the documentation, [1D Euler Problems](../../../docs/vv/1D.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[0, 1]$ |
| Final time | $t = 1.0$ |
| Time scheme | RK3 |
| CFL | 0.5 |
| Space reconstruction | first-order |
| Riemann solver | HLLE |
| Boundaries | extrapolation (`400`) at both ends and on the $j$ faces |
