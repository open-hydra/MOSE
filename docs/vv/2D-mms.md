# 2D Viscous Manufactured Solution

The **Method of Manufactured Solutions (MMS)** is the reference technique for
verifying the *observed order of accuracy* of a discretisation, and in particular
of the **viscous** operator of the compressible Navier–Stokes equations. Rather
than looking for a physically meaningful exact solution (which rarely exists for
the full viscous system), one *manufactures* one: a smooth analytic field is
chosen, substituted into the governing equations, and the non-zero residual that
remains is added back as an **analytic source term**. That source makes the
chosen field an *exact steady solution of the discretised equations in the
continuum limit*, so the error of any grid can be measured directly against a
field that is known everywhere.

Because the manufactured field is infinitely smooth, refining the mesh and
tracking how the error decays exposes the order of accuracy of the **complete**
discretisation — convective and viscous fluxes, the velocity and temperature
gradients that feed the stress tensor and heat flux, and the mesh metrics.

!!! note
    This case complements the [isentropic vortex](2D-isentropic-vortex.md). The vortex is an *inviscid* exact solution and certifies the convective operator, but it cannot exercise the viscous terms. MMS forces the viscous fluxes and their gradients, so it catches gradient/metric bugs the vortex cannot see.

## Methodology

Let $\mathbf{U}$ be the vector of conserved variables and write the governing
system as $\partial_t\mathbf{U} + \nabla\!\cdot\mathbf{F}(\mathbf{U}) = 0$, where
$\mathbf{F} = \mathbf{F}_c - \mathbf{F}_v$ collects the convective and viscous
fluxes. Choosing a smooth manufactured state $\hat{\mathbf{U}}(x,y)$ and
inserting it leaves a residual

$$
\mathbf{S}(x,y) \;=\; \nabla\!\cdot\mathbf{F}(\hat{\mathbf{U}})
\;=\; \nabla\!\cdot\big(\mathbf{F}_c(\hat{\mathbf{U}}) - \mathbf{F}_v(\hat{\mathbf{U}})\big),
$$

which is added to the right-hand side as a spatially varying source. With that
source in place, $\hat{\mathbf{U}}$ is by construction the exact steady solution;
the computed field converges to it at the scheme's formal order.

In MOSE the source is:

1. **derived symbolically** with `sympy` for the full 2D compressible
   Navier–Stokes system (calorically perfect ideal gas, constant $\mu$ and $k$,
   Stokes' hypothesis $\lambda=-\tfrac23\mu$);
2. **cross-checked** against an independent finite-difference evaluation of the flux divergence $\nabla\!\cdot\mathbf{F}$ (agreement to $\sim10^{-10}$), so a derivation slip cannot silently corrupt the study;
3. **injected each Runge–Kutta stage** by a dedicated driver —
   a build variant of the standard solver whose external per-stage callback adds
   $-\mathbf{S}\,V$ to the residual of every cell.

## Manufactured field

A smooth, periodic field is prescribed on the unit square $[0,1]^2$ with
$\omega = 2\pi$:

$$
\begin{aligned}
\rho &= 1.0   + 0.1\,\sin(\omega x)\cos(\omega y), \\
u    &= 40.0  + 8.0\,\sin(\omega x)\sin(\omega y), \\
v    &= 30.0  + 8.0\,\cos(\omega x)\cos(\omega y), \\
p    &= 8.0\times10^{3} + 8.0\times10^{2}\,\cos(\omega x)\sin(\omega y).
\end{aligned}
$$

Being $C^\infty$ and periodic, the field is represented exactly at the domain edges by the periodic connectivity, so **no boundary error contaminates the order study** — the measured error is purely interior truncation error.

<figure>
  {% include "vv/images/mms-solution.svg" %}
  Manufactured solutions
</figure>

## Governing equations and source

The solver integrates the compressible Navier–Stokes equations

$$
\partial_t\!\begin{pmatrix}\rho\\ \rho\mathbf{u}\\ \rho E\end{pmatrix}
+\nabla\!\cdot\!\begin{pmatrix}\rho\mathbf{u}\\ \rho\mathbf{u}\otimes\mathbf{u}+p\mathbf{I}-\boldsymbol{\tau}\\ (\rho E+p)\mathbf{u}-\boldsymbol{\tau}\mathbf{u}-k\nabla T\end{pmatrix}
=\mathbf{S},
$$

with $\boldsymbol{\tau}=\mu\big(\nabla\mathbf{u}+\nabla\mathbf{u}^{\mathsf T}\big)+\lambda(\nabla\!\cdot\mathbf{u})\mathbf{I}$,
$\lambda=-\tfrac23\mu$, $p=\rho R T$ and $E=p/[(\gamma-1)\rho]+\tfrac12|\mathbf{u}|^2$.
The transport properties are held **constant** (a flat viscosity table, and
$k=\mu c_p/\mathrm{Pr}$), which keeps the manufactured source analytic and free of
$\mathrm d\mu/\mathrm dT$ terms. The single-species gas makes the mass-diffusion
terms vanish identically, so the source stresses exactly the stress tensor and
Fourier heat flux under test.

| Parameter | Value |
|---|---|
| Domain | $[0,1]^2$, periodic |
| Gas | calorically perfect air ($R = 287$, $c_p = 1004.5$, $\gamma = 1.4$) |
| Viscosity | constant, $\mu = 10$ (flat transport table) |
| Prandtl | $\mathrm{Pr} = 0.72\;\Rightarrow\;k = \mu c_p/\mathrm{Pr}$ |
| Reynolds number | $\approx 5$ (viscous terms $\sim 20\%$ of convective) |
| Mean Mach number | $\approx 0.47$ |
| Equations | Navier–Stokes (laminar) |
| Time scheme | RK3, local time stepping to steady state ($\|\text{res}\|\sim10^{-9}$) |
| CFL / VNN | 0.2 / 0.1 |
| Space reconstruction | MUSCL, Van Leer |
| Riemann solver | HLLC |

!!! note
    $\mu=10$ the Reynolds number is only $\approx 5$, so the viscous fluxes carry a large fraction of the flux balance and the viscous operator is genuinely stressed (rather than a negligible correction).

## Verification results

The case is run on three successively refined grids ($32^2$, $64^2$, $128^2$,
refinement ratio 2), each marched to a converged steady residual. The normalised
$L_2$ and $L_\infty$ errors are measured for every primitive variable against the
manufactured field, and the observed order of accuracy is
$p = \log_2\!\big(e(N)/e(2N)\big)$.

**Normalised $L_2$ error (%)**

| Variable | $32^2$ | $64^2$ | $128^2$ | order $32\!\to\!64$ | order $64\!\to\!128$ |
|---|---|---|---|---|---|
| $\rho$ | 0.0592 | 0.0128 | 0.00282 | 2.21 | 2.19 |
| $u$    | 0.0205 | 0.00311 | 0.000900 | 2.72 | 1.79 |
| $v$    | 0.0352 | 0.00718 | 0.00166 | 2.29 | 2.11 |
| $p$    | 0.0519 | 0.0104 | 0.00204 | 2.33 | 2.34 |

**Normalised $L_\infty$ error (%)**

| Variable | $32^2$ | $64^2$ | $128^2$ | order $32\!\to\!64$ | order $64\!\to\!128$ |
|---|---|---|---|---|---|
| $\rho$ | 0.147 | 0.0405 | 0.0100 | 1.86 | 2.01 |
| $u$    | 0.0498 | 0.00695 | 0.00194 | 2.84 | 1.84 |
| $v$    | 0.0718 | 0.0143 | 0.00313 | 2.33 | 2.19 |
| $p$    | 0.130 | 0.0358 | 0.00879 | 1.85 | 2.03 |

Every variable converges at close to second order in both norms, confirming that
the full viscous MUSCL/RK3 discretisation is second-order accurate.

<figure>
  {% include "vv/images/mms.svg" %}
</figure>

## References

- P. J. Roache, *Code Verification by the Method of Manufactured Solutions*,
  J. Fluids Eng. **124**(1), 4–10 (2002).
- K. Salari and P. Knupp, *Code Verification by the Method of Manufactured
  Solutions*, Sandia report SAND2000-1444 (2000).
- SU2 Verification & Validation, *Method of Manufactured Solutions for the
  compressible Navier–Stokes equations*,
  [su2code.github.io/vandv](https://su2code.github.io/vandv/MMS_FVM_Navier_Stokes/).
