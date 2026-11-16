# Finite-Rate Reactive Shock Tube

The results are discussed in the documentation, [Finite-Rate Reactive Shock Tube](../../../docs/vv/1D-fer14.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Euler |
| Chemistry | Finite-rate |
| Domain | $[0, 0.12]$ m |
| Grid | $400 \times 1$ |
| Final time | $250\,\mu s$ |
| Output times used for comparison | $170\,\mu s$, $190\,\mu s$, $230\,\mu s$ |
| Time scheme | RK2 |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC+Tramel |
