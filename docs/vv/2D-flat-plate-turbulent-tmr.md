# Turbulent Flat Plate of the NASA TMR

The zero-pressure-gradient flat plate of the NASA Turbulence Modeling Resource (TMR) is the standard check of a turbulence-model implementation. Every code that implements a model correctly should approach the same skin friction as the grid is refined. MOSE runs it with four models, Spalart–Allmaras (SA), SST, Wilcox (2006) $k$-$\omega$ and the SSG/LRR-$\omega$ Reynolds-stress model, and compares the result with the TMR codes: CFL3D and FUN3D, or CFL3D and TAU for SSG/LRR-$\omega$.

This is a verification of the implementation, not a validation against experiments.

**Reference**: NASA TMR, 2D Zero Pressure Gradient Flat Plate Verification Case: [introduction](https://tmbwg.github.io/turbmodels/flatplate.html), [grids](https://tmbwg.github.io/turbmodels/flatplate_grids.html), results with [SA](https://tmbwg.github.io/turbmodels/flatplate_sa.html), [SST-Vm](https://tmbwg.github.io/turbmodels/flatplate_sst.html), [Wilcox2006-klim-m](https://tmbwg.github.io/turbmodels/flatplate_w06.html) and [SSGLRR-RSM-w2012](https://tmbwg.github.io/turbmodels/flatplate_ssglrrrsm.html).

## Summary

| Model | Grid limit of $c_f(0.97)$ against CFL3D on its finest grid | Observed order | Verdict |
|---|---|---|---|
| SA | +0.09% | 1.81 | **Verified**: near second order, GCI 0.46% |
| SST | +0.03% | 1.26 | Agrees, but the grids are not fine enough to be conclusive |
| Wilcox 2006 | +0.17% | 1.43 | Agrees, but the grids are not fine enough to be conclusive |
| SSG/LRR-$\omega$ | −0.03% | 1.38 | Agrees, but the grids are not fine enough to be conclusive |

On its finest grid MOSE lies within the spread of the two TMR codes for all four models. The TMR codes themselves converge slowly on these grids with SST, Wilcox 2006 and SSG/LRR-$\omega$ (observed orders 0.56 to 2.06), so a low observed order is a property of the grids, not a defect of MOSE.

---

## The case

A flat plate at $M = 0.2$ and $Re = 5\times10^6$ per unit length, with an adiabatic wall, in the domain and with the free-stream turbulence of the TMR. MOSE runs three grids of the TMR family, 35 × 25, 69 × 49 and 137 × 97 nodes, each obtained by taking every other node of the next finer one. The TMR codes ran two further grids, up to 545 × 385.

### Differences from the TMR

MOSE does not reproduce the TMR setup exactly. Each difference is small, and a margin $m$ in the acceptance criteria accounts for it.

| Difference | Effect on $c_f$ |
|---|---|
| Top boundary: outflow at $p_\infty$ instead of a far-field condition | about 0.3% (estimated from the edge velocity; the TMR reports < 0.2% for a similar change) |
| Constant viscosity instead of Sutherland's law | 0.09% (measured with SA on the 69 × 49 grid) |
| SA: without the $f_{t2}$ term | none for this case, according to the TMR |
| SST: SST-2003m instead of SST-Vm (strain rate instead of vorticity, other limiter and constants) | not measured; the TMR finds the variants nearly identical here |
| Wilcox 2006: exact production instead of $2\mu_t S_{ij}S_{ij}$ | below $3\times10^{-4}$ of the production |
| SSG/LRR-$\omega$: second-order advection of the turbulence variables, where CFL3D and TAU used first order | changes how $c_f$ approaches its limit, not the limit |

The margin is $m$ = 0.5% for SA and 1% for the other three models.

---

## Acceptance criteria

The quantity compared is the skin friction $c_f$ at $x = 0.97$. Two criteria must hold for each model:

1. **Finest grid.** On the 137 × 97 grid, MOSE must lie inside the band spanned by CFL3D and the second TMR code on the same grid, widened by $m\,c_f$.
2. **Grid limit.** The limit extrapolated from the three grids (Richardson extrapolation, with the observed order $p$) must agree with CFL3D on its finest grid, 545 × 385, within the uncertainty of the extrapolation $U$ plus $m\,c_f$.

$U$ is the grid convergence index (GCI) of Celik et al. (2008), as summarised by the TMR, plus the iterative error. A third check makes sure the runs are converged: the iterative error must be less than a tenth of the GCI, for $c_f(0.97)$, $c_f(1.90)$ and the drag coefficient. The full procedure is in the [test README](https://github.com/open-hydra/MOSE/blob/main/test/2D/viscous/flat-plate-turbulent-tmr/README.md).

---

## Results

### Skin friction on each grid

$c_f(0.97) \times 10^3$. The last column compares the distance of MOSE from the mean of the two TMR codes with the half-width of the band.

| Model | Grid | MOSE | CFL3D | FUN3D / TAU | Distance vs half-width | Inside |
|---|---|---|---|---|---|---|
| SA | 35 × 25 | 2.58510 | 2.79569 | 2.66738 | 0.1464 vs 0.0771 | no |
| | 69 × 49 | 2.67302 | 2.72835 | 2.69498 | 0.0386 vs 0.0301 | no |
| | **137 × 97** | 2.69807 | 2.71115 | 2.70215 | 0.0086 vs 0.0180 | **yes** |
| SST | 35 × 25 | 2.37528 | 2.55183 | 2.51562 | 0.1584 vs 0.0419 | no |
| | 69 × 49 | 2.55977 | 2.62625 | 2.60951 | 0.0581 vs 0.0340 | no |
| | **137 × 97** | 2.63669 | 2.66477 | 2.65845 | 0.0249 vs 0.0295 | **yes** |
| Wilcox 2006 | 35 × 25 | 2.42178 | 2.56221 | 2.56755 | 0.1431 vs 0.0269 | no |
| | 69 × 49 | 2.60919 | 2.65134 | 2.66962 | 0.0513 vs 0.0352 | no |
| | **137 × 97** | 2.67867 | 2.69148 | 2.69406 | 0.0141 vs 0.0281 | **yes** |
| SSG/LRR-$\omega$ | 35 × 25 | 2.44614 | 2.62550 | 2.55474 (TAU) | 0.1440 vs 0.0598 | no |
| | 69 × 49 | 2.59577 | 2.65419 | 2.63056 (TAU) | 0.0466 vs 0.0378 | no |
| | **137 × 97** | 2.65319 | 2.67370 | 2.66745 (TAU) | 0.0174 vs 0.0297 | **yes** |

Only the finest grid is required to pass. On the coarser grids the TMR codes also disagree with each other, by 4.7% with SA on the 35 × 25 grid.

### Grid limit

$c_f(0.97) \times 10^3$:

| | SA | SST | Wilcox 2006 | SSG/LRR-$\omega$ |
|---|---|---|---|---|
| Observed order $p$ | 1.81 | 1.26 | 1.43 | 1.38 |
| Limit | 2.70806 | 2.69170 | 2.71960 | 2.68894 |
| GCI | 0.46% | 2.61% | 1.91% | 1.68% |
| Uncertainty $U$ | 0.0126 | 0.0688 | 0.0512 | 0.0447 |
| CFL3D, 545 × 385 | 2.70562 | 2.69085 | 2.71510 | 2.68987 |
| FUN3D / TAU, 545 × 385 | 2.70540 | 2.69055 | 2.71770 | 2.69044 (TAU) |
| Limit − CFL3D | +0.090% | +0.031% | +0.166% | −0.034% |
| Used fraction of the allowed band | 0.09 | 0.009 | 0.06 | 0.013 |

On the same three grids the observed orders of the TMR codes are:

| | SA | SST-Vm | Wilcox2006-klim-m | SSGLRR-RSM-w2012 |
|---|---|---|---|---|
| CFL3D | 1.97 | 0.95 | 1.15 | 0.56 |
| FUN3D / TAU | 1.94 | 0.94 | 2.06 | 1.04 (TAU) |

$c_f(1.90)$ and the drag coefficient also converge monotonically, with orders from 1.21 to 2.02. All runs are iteratively converged: the largest iterative error is $6.2\times10^{-5}$ of the value.

### Figures

Left: $c_f(0.97)$ against the grid size $h = N^{-1/2}$ ($N$ the number of cells), MOSE on its three grids with the limit and its uncertainty, and the TMR codes on their five grids. Right: $c_f$ along the plate, MOSE on its three grids and the TMR codes on the 545 × 385 grid.

<figure>
  --8<-- "docs/vv/images/flat-plate-tmr-sa.svg"
  <figcaption>Spalart–Allmaras.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/flat-plate-tmr-sst.svg"
  <figcaption>SST (TMR data: SST-Vm).</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/flat-plate-tmr-wilcox2006.svg"
  <figcaption>Wilcox 2006 (TMR data: Wilcox2006-klim-m).</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/flat-plate-tmr-ssglrr.svg"
  <figcaption>SSG/LRR-ω (TMR data: SSGLRR-RSM-w2012, CFL3D and TAU).</figcaption>
</figure>

---

## SSG/LRR-ω in more detail

### Velocity profile and Reynolds stresses

<figure>
  --8<-- "docs/vv/images/flat-plate-tmr-ssglrr-profile.svg"
  <figcaption>SSG/LRR-ω on the 137 × 97 grid at x = 0.951 (the cell column nearest x = 0.97). Left: u⁺ against CFL3D on the 545 × 385 grid. Right: the Reynolds stresses in wall units.</figcaption>
</figure>

- The velocity profile $u^+$ differs from CFL3D by 1.2% RMS over $2 < y^+ < 1000$.
- In the log layer ($50 < y^+ < 300$) the normal stresses are $\overline{u'u'} : \overline{v'v'} : \overline{w'w'} = 0.908 : 0.437 : 0.656$ of $k$. An eddy-viscosity model would give 2/3 to each: the anisotropy is what a Reynolds-stress model adds.
- The shear stress is $-\overline{u'v'}/k = 0.305$, close to the value 0.31 of Bradshaw's relation that SST uses for its constant $a_1$.

The TMR gives the Reynolds stresses of CFL3D and TAU only as contour plots, so they are not compared.

### Model variants

On the 69 × 49 grid:

| Variant | $c_f(0.97) \times 10^3$ | Difference |
|---|---|---|
| `SSGLRR`, practical $\omega$ wall value (as the TMR) | 2.59577 | – |
| `SSGLRR-SD` (simple diffusion) | 2.61266 | +0.65% |
| `SSGLRR-w2019` (length-scale correction) | 2.59577 | none: the correction is inactive in an attached boundary layer |
| `SSGLRR`, asymptotic $\omega$ wall value | 2.68068 | +3.27% |

The wall value of $\omega$ has by far the largest effect. Turning off the point-implicit treatment of the sources changes $c_f$ by less than $3\times10^{-7}$: it changes the path to the steady state, not the state.

---

## Conclusions

- **SA is verified.** The grid limit agrees with CFL3D within 0.09%, the observed order is 1.81 (close to the formal 2), and the GCI is 0.46%.
- **SST, Wilcox 2006 and SSG/LRR-$\omega$ agree with the TMR codes, but not conclusively.** Both criteria pass, but the observed orders (1.26 to 1.43) and GCIs (1.7 to 2.6%) show that the three grids are not in the asymptotic range for these models. The criterion passes because the uncertainty is large, not because the limit is precise. The TMR codes behave the same way on these grids.
- **On the finest grid** MOSE lies inside the band of the two TMR codes for every model. With SSG/LRR-$\omega$, MOSE is still 1.3% below its own limit on that grid, while CFL3D and TAU are 0.6% and 0.85% below theirs: MOSE starts further from the limit and approaches it faster.

## Limits

- Only three grids: the two finest grids of the TMR are not run.
- The effect of the wall spacing ($y^+$ study of the TMR) is not run.
- The SST and Wilcox 2006 variants differ from those of the TMR data, and their effect is not measured.
- The Reynolds stresses of SSG/LRR-$\omega$ are not compared with the TMR codes.
- On some grids the density residual stops falling (SA on 35 × 25 at about $4\times10^{-7}$; SSG/LRR-$\omega$ on the two finer grids at about $8\times10^{-7}$ and $5\times10^{-6}$), but $c_f$ no longer changes.
