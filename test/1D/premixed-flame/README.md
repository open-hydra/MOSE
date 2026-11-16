# Premixed Laminar Flame

The results are discussed in the documentation, [Premixed Laminar Flame](../../../docs/vv/1D-premixed-flame.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier–Stokes (laminar) + finite-rate chemistry |
| Mechanism | ZK (25 species) |
| Domain | $[0,\,0.30]$ m |
| Grid | $4000 \times 1$ |
| Multigrid | 5 levels (FMG start-up) |
| Time scheme | RK2 (time-accurate), Strang splitting for chemistry |
| ODE solver | ROS4 |
| CFL | 0.9 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Albada |
| Riemann solver | HLLC (Batten) |
| Schmidt / Prandtl | 0.7 / 0.7 |

The inflow is a subsonic inlet and the outflow a fixed-pressure outlet; the spanwise faces use symmetry, reducing the problem to 1-D along $i$.
