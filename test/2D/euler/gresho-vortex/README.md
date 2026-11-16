# Gresho Vortex

The results are discussed in the documentation, [Gresho Vortex](../../../../docs/vv/2D-gresho-vortex.md); this file holds the setup of the case and how to run it.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 (time-accurate) |
| Final time | $t = 1$ |
| CFL | 0.6 |
| Space reconstruction | first order |
| Low-Mach cutoff | $M_{co} = 0.005$ ($\approx 5\,M_{\rm peak}$) |
| **Riemann solvers tested** | **AUSM+M** (`input.ini`), **HLLC** (baseline), **LMRoe** |

The cutoff $M_{co}$ is the shared low-Mach dissipation floor described in
[Low-Mach Cutoff Mach $M_{co}$](../../../../docs/theory/riemann-solvers.md#low-mach-cutoff-mach-m_co).
**HLLC** is the *control*: it carries no low-Mach treatment and shows what the
baseline upwind dissipation does to the vortex. **LMRoe** and **AUSM+M** are the
low-Mach-corrected schemes. HLLC and LMRoe run with `riemann-solver` changed in
`input.ini`.
