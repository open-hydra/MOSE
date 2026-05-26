# Finite-Rate Reactive Shock Tube

One-dimensional finite-rate reactive shock tube using an H2/O2/Ar mixture and a reduced finite-rate mechanism.
The test validates thermochemical coupling by comparing MOSE velocity and temperature profiles against reference curves at three times: 170, 190, and 230 microseconds.

### Initial conditions

| Region | Composition (mass fractions) | $T$ [K] | $u$ [m/s] | $p$ [Pa] |
|---|---|---:|---:|---:|
| Left ($x < 0.06$ m) | $Y_{H2}=0.0127$, $Y_{O2}=0.1013$, $Y_{Ar}=0.8861$ | 378.656 | 0.0 | 7173.0 |
| Right ($x \geq 0.06$ m) | $Y_{H2}=0.0127$, $Y_{O2}=0.1013$, $Y_{Ar}=0.8861$ | 748.472 | -487.34 | 35594.0 |

### Numerical setup

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
| Space reconstruction | MUSCL-SD |
| Flux limiter | Van Leer |
| Riemann solver | HLLC+ |

### Results and verification

<figure>
  {% include "vv/images/Fer14.svg" %}
</figure>

A perfect agreement is observed between MOSE and the reference data, demonstrating accurate thermochemical coupling in the solver. The velocity and temperature profiles at all three time points closely match the digitized reference curves from Ferrer et al. (2014).
