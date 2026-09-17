# Unity-Lewis Diffusion

A one-dimensional, multicomponent **diffusion-only** test: a quiescent CH$_4$ /
air mixture with a sharp temperature and composition front that relaxes purely by
molecular transport. With no convective forcing, the case isolates the
**species and thermal diffusion** closure of the [viscous
fluxes](../theory/governing-equations.md) — diffusion coefficients, mixture
conductivity, and the inter-species enthalpy flux — from the rest of the solver.

This page uses a **unity Lewis number** ($Le = 1$) closure, which makes every
species and thermal diffusivity equal to $\mu/\rho$ and admits an exact
apples-to-apples Cantera reference. For the temperature- and composition-dependent
**mixture-averaged** closure, see [Multicomponent
Diffusion](1D-multicomponent-diffusion.md).

Because there is no analytical solution for a real multicomponent mixture, MOSE
is validated against an **independent low-Mach unity-Lewis solver built on
Cantera** (GRI-3.0 thermodynamics and transport). The reference solves the *same*
physics on the *same* grid, so any discrepancy is purely numerical.

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
| Laminar Schmidt $Sc$ | 1.0 |
| Laminar Prandtl $Pr_l$ | 1.0 |

Both transverse boundaries are **periodic** and the spanwise faces use symmetry,
so the only error sources are the diffusive flux and the time integration.

## Unity Lewis number and the reference solution

The transport closure is set to **unity Lewis number** ($Le = Sc/Pr = 1$) by
imposing $Sc = 1$ and the laminar override $Pr_l = 1$ (see [diffusive
fluxes](../theory/governing-equations.md)). This makes every species diffusivity and
the thermal diffusivity equal to $\mu/\rho$:

$$
D_k = \frac{\mu}{\rho\,Sc} = \frac{\mu}{\rho}, \qquad
\kappa = \frac{\mu\,c_p}{Pr_l} = \mu\,c_p .
$$

Under $Le = 1$ the species mass fractions and the specific enthalpy obey the
*same* scalar diffusion equation with coefficient $\mu$ — the standard collapse
$\kappa\nabla T - \sum_k h_k \mathbf{j}_k = \rho D\,\nabla h = \mu\,\nabla h$.
MOSE carries exactly this: its diffusive energy flux includes the inter-species
enthalpy term $\rho\sum_k D_k(\nabla Y_k)h_k$ alongside $\kappa\nabla T$.

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

where $L_{Y_k} = \partial_x(\mu\,\partial_x Y_k)$ and $L_h =
\partial_x(\mu\,\partial_x h)$ are the diffusion operators. In 1D periodic this is
integrated directly for the face velocity (the mean is removed, i.e. $p_0$ is
held fixed), and the partial densities are advected conservatively so that $\rho =
\sum_k \rho Y_k$ keeps continuity satisfied automatically.

!!! warning
    Omitting this velocity leaves a systematic $\sim\!10\%$ gap at the front centre
    (MOSE diffuses *more*); including it closes the gap to $\sim\!1\%$, confirming the
    discrepancy is the compressible transport physics.

## Results and verification

Temperature and CH$_4$ mass-fraction profiles at $t = 50, 100, 150, 200$ ms.
Solid lines are MOSE; open circles are the Cantera low-Mach reference.

<figure>
  {% include "vv/images/diffusion_T.svg" %}
</figure>

<figure>
  {% include "vv/images/diffusion_CH4.svg" %}
</figure>

The reference markers sit on the MOSE curves at every time and every location —
the hot centre, the steep flanks, and the far field. The centre values track to
about $1\%$ throughout:

| $t$ [ms] | $T$ centre (MOSE / Cantera) [K] | $Y_{\mathrm{CH_4}}$ centre (MOSE / Cantera) |
|---:|---:|---:|
| 50  | 742 / 742 | 0.123 / 0.124 |
| 200 | 556 / 555 | 0.164 / 0.164 |
