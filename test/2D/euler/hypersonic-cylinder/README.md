# Hypersonic Cylinder

The results are discussed in the documentation, [Hypersonic Cylinder](../../../../docs/vv/2D-hypersonic-cylinder.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 |
| Time stepping | local (steady), 5000 iterations |
| CFL | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| **Riemann solvers tested** | **HLLC, HLLC+Tramel, HLLE, SLAU** |
