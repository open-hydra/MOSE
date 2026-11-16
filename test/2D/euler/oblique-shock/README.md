# Oblique Shock

The results are discussed in the documentation, [Oblique Shock](../../../../docs/vv/2D-oblique-shock.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+Tramel, HLLE, SLAU** |
