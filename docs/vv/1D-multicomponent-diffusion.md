# Multicomponent Diffusion

A one-dimensional, multicomponent **diffusion-only** test: a quiescent CH$_4$ /
air mixture with a sharp temperature and composition front that relaxes purely by
molecular transport. With no convective forcing, the case isolates the
**species and thermal diffusion** closure of the [viscous
fluxes](../theory/governing-equations.md) — the per-species mixture-averaged
diffusion coefficients, the mixture conductivity, and the inter-species enthalpy
flux — from the rest of the solver.

This is the temperature- and composition-dependent **mixture-averaged** companion
to the [Unity-Lewis Diffusion](1D-diffusion.md) case (same geometry and states):
there every diffusivity collapses to $\mu/\rho$, whereas here each species carries
its own coefficient $D_k$.

Because there is no analytical solution for a real multicomponent mixture, MOSE
is validated against an **independent low-Mach mixture-averaged solver built on
Cantera** (GRI-3.0 thermodynamics and transport). The reference solves the *same*
physics with the *same* diffusion model on the *same* grid, so any discrepancy is
purely numerical.

---

## Problem setup

A periodic domain of length $L = 0.05$ m at uniform pressure $p = 101325$ Pa and
zero velocity. The initial temperature and species mass fractions follow a smooth
Gaussian front between an **oxidizer** state at the domain centre and a **fuel**
state in the far field,

$$
\phi(x) = \phi_0 + (\phi_f - \phi_0)\,f(x), \qquad
f(x) = 1 - \exp\!\left[-\frac{(x - x_0)^2}{d^2}\right],
$$

with $x_0 = 25\times10^{-3}$ m and $d = 2.5\times10^{-3}$ m. Subscript $0$ denotes
the oxidizer side ($f=0$ at the centre), subscript $f$ the fuel side
($f\to 1$ far field).

| State | $T$ [K] | $Y_{\mathrm{CH_4}}$ | $Y_{\mathrm{H_2O}}$ | $Y_{\mathrm{N_2}}$ | $Y_{\mathrm{O_2}}$ |
|---|---:|---:|---:|---:|---:|
| Oxidizer (centre, $f=0$) | 1350 | 0.000 | 0.100 | 0.758 | 0.142 |
| Fuel (far field, $f=1$)  |  320 | 0.214 | 0.000 | 0.591 | 0.195 |

Density follows from the ideal-gas EoS, $\rho = p\,W_{\rm mix}/(R_u T)$ with
$1/W_{\rm mix} = \sum_k Y_k/W_k$, and partial densities $\rho_k = Y_k\rho$. The
$\sim\!4{:}1$ temperature ratio gives a $\sim\!4{:}1$ density ratio across the
front — the key feature of this case (see below).

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier–Stokes (laminar) |
| Domain | $[0,\,0.05]$ m |
| Grid | $200 \times 1$ |
| Final time | $200$ ms (snapshots every $50$ ms) |
| Time scheme | RK2 (time-accurate) |
| CFL / VNN | 0.95 / 0.95 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |
| Laminar Schmidt $Sc$ | $\le 0$ (multicomponent) |
| Laminar Prandtl $Pr_l$ | $0$ (computed mixture conductivity) |

Both transverse boundaries are **periodic** and the spanwise faces use symmetry,
so the only error sources are the diffusive flux and the time integration.

!!! note "Selecting the model"
    Mixture-averaged diffusion is selected with `Sc <= 0` in the `[MOSE-Physics]`
    section (`Sc = -1` here); a positive `Sc` recovers the constant-Schmidt closure
    $D_k = \mu/(\rho\,Sc)$. Leaving `Prl = 0` keeps the temperature-dependent Wilke
    mixture conductivity. The binary-diffusion table `INPUT/diffusion.dat` (generated
    by the preprocessor) must be present.

## Mixture-averaged multicomponent diffusion

Each species $k$ is given its own composition- and temperature-dependent diffusion
coefficient through the **Curtiss–Hirschfelder mixture-averaged** rule,

$$
D_k = \frac{1 - X_k}{\displaystyle\sum_{j \ne k} X_j / \mathcal{D}_{kj}},
$$

where $X_k$ are mole fractions and $\mathcal{D}_{kj}(T,p)$ are the binary diffusion
coefficients, tabulated from kinetic theory at a reference pressure and rescaled as
$\mathcal{D}\propto 1/p$ at run time. The species mass flux uses this $D_k$ with the
mass-fraction gradient and is closed by the mass-conservation correction velocity so
that $\sum_k \mathbf{j}_k = 0$:

$$
\mathbf{j}_k = -\rho D_k \nabla Y_k
             + Y_k \sum_{l} \rho D_l \nabla Y_l .
$$

The diffusive energy flux carries the mixture conductivity together with the
inter-species enthalpy transport,

$$
\mathbf{q} = -\kappa\,\nabla T - \rho \sum_k h_k D_k \nabla Y_k ,
$$

with $\kappa$ the full Wilke mixture conductivity. The Cantera reference evaluates
$D_k$ with the **same $(1-X_k)$ rule** built from the same binary coefficients, so
the comparison measures discretization error only.

!!! warning "Mixture rule: $(1-X_k)$ vs $(1-Y_k)$"
    MOSE uses the mole-fraction numerator $(1-X_k)$. Cantera's `mix_diff_coeffs`
    and the Chemkin/Kee standard use the mass-fraction numerator $(1-Y_k)$; the two
    share the same denominator and binary data but differ by $\sim\!10\text{–}20\%$
    per species in light/heavy mixtures. A meaningful verification must therefore use
    the **same rule on both sides** — the reference here does. The $(1-Y_k)$ variant
    is shipped as commented reference code in the diffusion routine for users who
    prefer the Chemkin/Cantera convention.

## Dilatation velocity

The reference is **not** a zero-velocity diffusion model. With a $4{:}1$ density
contrast, the gas expands and contracts as heat diffuses, inducing a
**dilatational velocity** that advects the scalars — precisely what MOSE's
compressible solver generates. At constant thermodynamic pressure, continuity and
the ideal-gas EoS fix the velocity divergence (low-Mach constraint):

$$
\nabla\!\cdot\mathbf{u}
= \frac{D\ln T}{Dt} - \frac{D\ln W}{Dt}
= \frac{L_h - \sum_k h_k L_{Y_k}}{\rho\,c_p\,T}
+ \frac{W}{\rho}\sum_k \frac{L_{Y_k}}{W_k},
$$

where $L_{Y_k}$ and $L_h$ are the multicomponent diffusion operators for the
partial densities and the enthalpy. This relation is closure-independent (it follows
only from continuity and the EoS). In 1D periodic it is integrated directly for the
face velocity (the mean is removed, i.e. $p_0$ is held fixed), and the partial
densities are advected conservatively so that $\rho = \sum_k \rho Y_k$ keeps
continuity satisfied automatically.

!!! warning
    Omitting this velocity leaves a systematic $\sim\!10\%$ gap at the front centre
    (MOSE diffuses *more*); including it closes the gap to $\sim\!1\%$, confirming the
    discrepancy is the compressible transport physics.

## Results and verification

Temperature, CH$_4$ and N$_2$ mass-fraction profiles at $t = 50, 100, 150, 200$ ms.
Solid lines are MOSE; open circles are the Cantera low-Mach mixture-averaged
reference.

<figure>
  {% include "vv/images/multicomp_T.svg" %}
</figure>

<figure>
  {% include "vv/images/multicomp_CH4.svg" %}
</figure>

<figure>
  {% include "vv/images/multicomp_N2.svg" %}
</figure>

The reference markers sit on the MOSE curves at every time and every location —
the hot centre, the steep flanks, and the far field. The centre values track to
about $1\%$ throughout:

| $t$ [ms] | $T$ centre (MOSE / Cantera) [K] | $Y_{\mathrm{CH_4}}$ centre | $Y_{\mathrm{N_2}}$ centre |
|---:|---:|---:|---:|
|  50 | 689 / 692 | 0.134 / 0.135 | 0.655 / 0.654 |
| 100 | 595 / 598 | 0.155 / 0.155 | 0.639 / 0.638 |
| 150 | 549 / 551 | 0.165 / 0.165 | 0.630 / 0.630 |
| 200 | 521 / 522 | 0.171 / 0.171 | 0.625 / 0.625 |

Over the full transient the relative errors stay at or below $\sim\!1\%$ (Linf) and
$\sim\!0.3\%$ (L2) for all three quantities — species and thermal diffusion errors
are the same order, the signature of clean numerics rather than a model mismatch.

