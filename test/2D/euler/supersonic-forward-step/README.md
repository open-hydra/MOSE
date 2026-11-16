# Woodward-Colella Supersonic Step

The results are discussed in the documentation, [Woodward-Colella Supersonic Step](../../../../docs/vv/2D-woodward-colella.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Domain | $[0, 3] \times [0, 1]$ |
| Final time | $t = 4.0$ |
| Mach number | 3.0 |
| Time scheme | RK3 |
| CFL | 0.8 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLE |
