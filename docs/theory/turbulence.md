# Turbulence Modelling

MOSE provides Reynolds-Averaged Navier–Stokes (RANS) closures for
turbulent flows, selected at run time via the `turbulence` key in the
`[MOSE-Physics]` section of the input file.  The eddy-viscosity models
described in detail below range from the one-equation Spalart–Allmaras
model (with several variants) to the two-equation Menter SST and
Wilcox 2006 $k$–$\omega$ models.  A full Reynolds-stress closure
(SSG–LRR) and an algebraic non-linear correction (QCR2000) are also
available — see [Other models](#other-models).

---

## Boussinesq Hypothesis

All eddy-viscosity models assume that the Reynolds stress tensor is
proportional to the mean strain rate:

$$
\tau_{ij}^R = 2\,\mu_t\,S_{ij} - \tfrac{2}{3}\,\rho\,k\,\delta_{ij}
$$

where $S_{ij} = \tfrac{1}{2}(\partial v_i / \partial x_j + \partial v_j / \partial x_i)$
is the mean strain-rate tensor, $\mu_t$ is the eddy (turbulent) viscosity,
and $k$ is the turbulence kinetic energy.  The isotropic term
$-\tfrac{2}{3}\rho k\,\delta_{ij}$ is included for two-equation models
and omitted for SA.

---

## Spalart–Allmaras (SA) — One-Equation Model

### Transport equation

The SA model solves a single transport equation for the modified
turbulent viscosity $\tilde\nu$:

$$
\frac{\partial(\rho\tilde\nu)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,\tilde\nu)
= \underbrace{c_{b1}\,\tilde{S}\,\rho\tilde\nu}_{\text{Production}}
+ \underbrace{\frac{1}{\sigma}\Bigl[
  \nabla\!\cdot\!\bigl((\mu + \rho\tilde\nu)\,\nabla\tilde\nu\bigr)
  + c_{b2}\,\rho\,|\nabla\tilde\nu|^2
\Bigr]}_{\text{Diffusion}}
- \underbrace{c_{w1}\,f_w\,\rho\!\left(\frac{\tilde\nu}{y}\right)^{\!2}}_{\text{Destruction}}
$$

### Model constants

| Constant | Value | Description |
|:--------:|:-----:|-------------|
| $c_{b1}$ | 0.1355 | Production coefficient |
| $c_{b2}$ | 0.622 | Diffusion coefficient |
| $\sigma$ | 2/3 | Turbulent Schmidt number |
| $\kappa$ | 0.41 | von Kármán constant |
| $c_{w1}$ | $c_{b1}/\kappa^2 + (1+c_{b2})/\sigma$ | Destruction coefficient |
| $c_{w2}$ | 0.3 | Destruction coefficient |
| $c_{w3}$ | 2.0 | Destruction coefficient |
| $c_{v1}$ | 7.1 | Damping-function constant |

### Auxiliary functions

$$
\chi = \frac{\tilde\nu}{\nu}, \qquad
f_{v1}(\chi) = \frac{\chi^3}{\chi^3 + c_{v1}^3}, \qquad
f_{v2}(\chi) = 1 - \frac{\chi}{1 + \chi\,f_{v1}(\chi)}
$$

Modified vorticity:

$$
\tilde{S} = \Omega + \frac{\tilde\nu}{\kappa^2 y^2}\,f_{v2}(\chi)
$$

where $\Omega = \sqrt{2\,W_{ij}\,W_{ij}}$ is the vorticity magnitude and
$y$ is the distance to the nearest wall.

Destruction function:

$$
r = \frac{\tilde\nu}{\kappa^2 y^2 \tilde{S}}, \qquad
g = r + c_{w2}(r^6 - r), \qquad
f_w = g\!\left(\frac{1 + c_{w3}^6}{g^6 + c_{w3}^6}\right)^{\!1/6}
$$

### Eddy viscosity

$$
\mu_t = \rho\,\tilde\nu\,f_{v1}(\chi)
$$

### Wall boundary condition

$$
\tilde\nu_\text{wall} = 0
$$

---

### SA Variants

#### SAcomp — Compressibility Correction

Paciorri–Sabetta correction that scales the production term with a
function of the turbulent stress ratio:

$$
S_\tau = \frac{\omega\,\tilde\nu\,f_{v1}}{a^2}
$$

Activated when compressibility effects on turbulence are significant
(high-speed boundary layers, mixing layers).

#### SAR — Rotation Correction

Adds a rotation term proportional to the difference between the strain
rate and the vorticity magnitude:

$$
P_\text{SAR} = P_\text{SA} + c_\text{rot}\,(\lVert S\rVert - \Omega)
\qquad (c_\text{rot} = 2.0)
$$

#### SA-RC — Spalart–Shur Rotation-Curvature Correction

A more sophisticated rotation/curvature correction that modulates the
production term via a correction factor $f_{r1}$:

$$
r^\ast = \frac{S}{\Omega}, \qquad
\tilde{r} = \frac{C_D}{D_A^2}
$$

where $C_D$ involves the material derivative interaction between the
strain and vorticity tensors:

$$
C_D = W_{ij}\,S_{jk}\!\left(\frac{DS_{ki}}{Dt} - \frac{DW_{ki}}{Dt}\right)
$$

The correction factor:

$$
f_{r1} = (1 + c_{r1})\,\frac{2\,r^\ast}{1 + r^\ast}\,
  \bigl(1 - c_{r3}\arctan(c_{r2}\,\tilde{r})\bigr) - c_{r1}
$$

bounded to $[0,\;1.25]$, with constants $c_{r1} = 1$, $c_{r2} = 2$,
$c_{r3} = 1$.  Production is multiplied by $\max(1,\,f_{r1})$.

---

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear Stress Transport model blends a $k$–$\omega$ formulation
(near walls) with a $k$–$\varepsilon$-like behaviour (in the
freestream) using blending functions.

### Transport equations

$$
\frac{\partial(\rho k)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,k)
= P_k - \beta^\ast\,\rho\,\omega\,k
+ \nabla\!\cdot\!\bigl[(\mu + \mu_t/\sigma_{k})\,\nabla k\bigr]
$$

$$
\frac{\partial(\rho\omega)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,\omega)
= \gamma\,\frac{\rho\,P_k}{k} - \beta\,\rho\,\omega^2
+ \nabla\!\cdot\!\bigl[(\mu + \mu_t/\sigma_{\omega})\,\nabla\omega\bigr]
+ 2(1 - F_1)\frac{\rho}{\sigma_{\omega 2}\,\omega}\,\nabla k\!\cdot\!\nabla\omega
$$

The last term is the **cross-diffusion** term, active only in the
freestream ($F_1 \to 0$).

### Production limiter

$$
P_k = \min\!\bigl(\mu_t\,S^2,\;10\,\beta^\ast\,\rho\,\omega\,k\bigr)
$$

This prevents unbounded growth of $k$ in stagnation regions.

!!! note "$\omega$-production uses the *limited* $P_k$"
    Per the NASA Turbulence-Modelling-Resource SST-2003 specification, the
    $\omega$-equation production is

    $$
    P_\omega = \frac{\alpha}{\nu_t}\,\tilde P_k
             = \frac{\gamma\,\rho}{\mu_t}\,\tilde P_k, \qquad
    \tilde P_k = \min\!\bigl(\mu_t S^2,\;10\,\beta^\ast\rho\,\omega\,k\bigr),
    $$

    i.e. it is built from the **same limited** production $\tilde P_k$ used in
    the $k$-equation, *not* the unlimited $\gamma\rho S^2$.  (The unlimited form
    that appears in the original 2003 paper is a typographical error, corrected
    on the TMR page.)  MOSE therefore applies the production limiter first and
    forms $P_\omega$ from the limited value; this matches the reference SST-2003
    codes (e.g. SU2's default `V2003`).

### Eddy viscosity

$$
\mu_t = \frac{\rho\,k\,a_1}{\max(a_1\,\omega,\;S\,F_2)}
$$

where $S = \sqrt{2\,S_{ij}\,S_{ij}}$ is the strain-rate magnitude and
$a_1 = 0.31$.

### Model constants

All blended coefficients $\phi$ are computed as
$\phi = F_1\,\phi_1 + (1 - F_1)\,\phi_2$.

| Constant | Set 1 ($\phi_1$) | Set 2 ($\phi_2$) |
|:--------:|:-----------------:|:-----------------:|
| $\sigma_k$ | 0.85 | 1.0 |
| $\sigma_\omega$ | 0.5 | 0.856 |
| $\beta$ | 0.075 | 0.0828 |
| $\gamma$ | 5/9 | 0.44 |

Universal: $\beta^\ast = 0.09$, $a_1 = 0.31$, $\kappa = 0.41$.

### Blending functions

$$
F_1 = \tanh\!\bigl(\arg_1^4\bigr), \qquad
\arg_1 = \min\!\left(
  \max\!\left(\frac{\sqrt{k}}{0.09\,\omega\,y},\;
              \frac{500\,\nu}{\omega\,y^2}\right),\;100\right)
$$

$$
F_2 = \tanh\!\bigl(\arg_2^2\bigr), \qquad
\arg_2 = \max\!\left(\frac{2\sqrt{k}}{0.09\,\omega\,y},\;
                      \frac{500\,\nu}{\omega\,y^2}\right)
$$

### Wall boundary conditions

$$
k_\text{wall} = 0, \qquad
\omega_\text{wall} = 10\,\frac{6\,\nu}{\beta_1\,y^2} = \frac{800\,\nu}{y^2}
$$

with $\beta_1 = 0.075$.  The factor of 10 over the analytical near-wall limit
$6\nu/(\beta_1 y^2)$ is Menter's recommended over-specification, which forces
the correct $\omega$ behaviour in the first off-wall cell.  Note that this
value grows as $1/y^2$ under grid refinement and is the origin of the
near-wall stiffness handled by the [point-implicit source
treatment](#numerical-treatment-of-the-source-terms).

### Energy coupling *(optional)*

When enabled, the turbulence kinetic energy production/destruction
contributes to the mean-flow energy equation:

$$
\frac{\partial(\rho e)}{\partial t} \mathrel{+}= S_k
$$

---

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

### Model constants

| Constant | Value |
|:--------:|:-----:|
| $\sigma_k$ | 0.6 |
| $\sigma_\omega$ | 0.5 |
| $\beta^\ast$ | 0.09 |
| $\beta_0$ | 0.0708 |
| $\gamma$ | 13/25 |
| $C_\text{lim}$ | 7/8 |
| $\sigma_d$ | 1/8 |

### Eddy viscosity

$$
\mu_t = \frac{\rho\,k}{\hat\omega}, \qquad
\hat\omega = \max\!\left(\omega,\;
  C_\text{lim}\,\frac{\sqrt{2\,S_{ij}\,S_{ij}}}{\beta^\ast}\right)
$$

The limiter $C_\text{lim}$ prevents excessive eddy viscosity in regions
where $\omega$ is small relative to the strain rate.

### Destruction with stress-limiter correction

The $\omega$-destruction coefficient is modified by Wilcox's
vortex-stretching function $f_\beta$:

$$
\beta = \beta_0\,f_\beta, \qquad
X_\omega = \frac{|W_{ij}\,W_{jk}\,S_{ki}|}{(\beta^\ast\omega)^3},
\qquad
f_\beta = \frac{1 + 85\,X_\omega}{1 + 100\,X_\omega}
$$

### Cross-diffusion *(conditional)*

Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega \le 0$:

$$
\text{CD} = \frac{\sigma_d\,\rho}{\omega}\,
  (\nabla k\!\cdot\!\nabla\omega)
\qquad \text{if } \nabla k\!\cdot\!\nabla\omega \le 0
$$

### Production limiter

More permissive than SST:

$$
P_k = \min\!\bigl(\mu_t\,S^2,\;20\,\beta^\ast\,\rho\,\omega\,k\bigr)
$$

---

## Other models

### SSG–LRR (Reynolds-stress model)

In addition to the eddy-viscosity closures above, MOSE provides the
**SSG–LRR** differential Reynolds-stress model (`turbulence = SSGLRR`),
which transports the six independent components of the Reynolds-stress
tensor together with a length-scale variable instead of assuming the
Boussinesq relation. It blends the Speziale–Sarkar–Gatski (SSG)
pressure–strain model away from walls with the Launder–Reece–Rodi (LRR)
model near walls. Being anisotropy-resolving, it captures secondary
flows and strong streamline-curvature effects that the two-equation
models miss, at the cost of the additional transport equations.

### QCR2000 (non-linear constitutive correction)

The **Quadratic Constitutive Relation** (Spalart, 2000) is an algebraic,
non-linear correction to the Boussinesq stress rather than a standalone
model. It is enabled as a suffix on a base model (e.g.
`turbulence = SA-QCR2000`) and adds a quadratic dependence on the
rotation tensor to the modelled stress, improving the prediction of
anisotropy-driven secondary flows (corner flows, square ducts) while
reusing the eddy viscosity of the underlying model.

---

## General RANS Features

### Procedure-pointer architecture

All turbulence models in MOSE are accessed through **function pointers**
defined in `Mod_RANS`:

| Pointer | Purpose |
|---------|---------|
| `Eddy_Viscosity` | Compute $\mu_t$ from model variables |
| `RANS_Diffusive_Flux` | Turbulent diffusion terms for $k$, $\omega$, $\tilde\nu$ |
| `Stress_Vector` | Viscous + Reynolds stress on a face |
| `RANS_Set_Wall_Values` | Set turbulence BC at walls |
| `RANS_Extrapolate_Wall` | Ghost-cell extrapolation for RANS variables |

This design allows switching models at run time without recompilation.

### Blowing / wall-suction correction

An optional correction (`blowing_corr`) modifies the wall boundary-layer
treatment for SST and Wilcox 2006 models in the presence of wall
injection or suction.

---

## Freestream and inlet values ($k$, $\omega$, $\tilde\nu$)

For the two-equation models the freestream/inlet values of $k$ and $\omega$
must be chosen consistently — the destruction of $k$ scales with $\omega$
($D_k = \beta^\ast\rho\,\omega\,k$), so an $\omega$ that is set too low
starves the $k$-equation of dissipation and lets $k$ grow without bound.
The recommended (NASA Turbulence-Modelling-Resource) freestream values are
built from a turbulence intensity $Tu$ and an eddy-to-molecular viscosity
ratio $\mu_t/\mu$:

$$
k_\infty = \tfrac{3}{2}\,\bigl(Tu\,U_\infty\bigr)^2, \qquad
\omega_\infty = \frac{\rho_\infty\,k_\infty}{\mu\,(\mu_t/\mu)}
             = \frac{k_\infty}{\nu_\infty\,(\mu_t/\mu)} .
$$

Typical verification values are $Tu \approx 0.04\%\!-\!1\%$ and
$\mu_t/\mu \approx 0.009\!-\!1$ (smaller $\mu_t/\mu$ ⇒ larger $\omega_\infty$
⇒ more near-wall dissipation and a more robust start-up).  For example, a
Mach-5 stream at $p=4000$ Pa, $T=68.3$ K ($\rho=0.204$ kg m⁻³,
$U_\infty=828$ m s⁻¹, $\nu_\infty=5.8\times10^{-5}$ m² s⁻¹) with
$Tu=0.04\%$, $\mu_t/\mu=0.009$ gives $k_\infty\approx0.16$ m² s⁻²,
$\omega_\infty\approx3\times10^{5}$ s⁻¹.  Values orders of magnitude below
this (e.g. $\omega_\infty=100$ s⁻¹) are a frequent cause of $k$ runaway.

The same values are set for both the initial condition (`[ICB-Block*]`) and
the inflow (`[inflow]`) in the input file.

---

## Numerical treatment of the source terms

### Point-implicit (Patankar) destruction

The turbulence source terms are added to the residual and advanced with the
same explicit Runge–Kutta step as the mean flow.  The **destruction** terms,
however, are stiff near walls — for $k$–$\omega$ models the wall value
$\omega_\text{wall}=6\nu/(\beta_1 y^2)$ grows like $1/y^2$, so on a fine
near-wall grid the explicit stability limit $\Delta t < 1/(\beta\,\omega)$ is
easily violated and $k$/$\omega$ diverge while the mean flow stays healthy
(typically at a multigrid fine-grid transition).

To remove this restriction MOSE treats the destruction terms
**point-implicitly**.  Writing the update for a conserved turbulence
variable $q=\rho\phi$ as $q^{n+1}=q^n+\Delta t\,S(q)$ and linearising only
the (stabilising) destruction part $D$,

$$
\Delta q = \Delta t\,\bigl[S^n - d\,\Delta q\bigr]
\;\Longrightarrow\;
\Delta q = \frac{\Delta t\,S^n}{1 + \Delta t\,d}, \qquad
d \equiv \frac{\partial D}{\partial q}\ge 0 .
$$

In practice the net source increment is simply divided by $(1+\Delta t\,d)$,
using the local time step $\Delta t$.  The destruction Jacobians are

| Model | Equation | $d = \partial D/\partial q$ |
|-------|----------|------------------------------|
| SST / Wilcox | $k$ | $\beta^\ast\,\omega$ |
| SST / Wilcox | $\omega$ | $2\,\beta\,\omega$ |
| SA | $\tilde\nu$ | $2\,c_{w1}\,f_w\,\tilde\nu/y^2$ |

Because the factor multiplies only the *increment*, it vanishes at
convergence ($\Delta q\to0$) and therefore **does not change the converged,
zero-residual solution** — it only enlarges the stable time step. This makes
SST/Wilcox run at the mean-flow CFL independent of near-wall spacing.

!!! warning "Effect on time-accurate (URANS) runs"
    For **steady** computations (local time-stepping) the treatment is exact:
    at convergence the residual is zero, so the factor has no effect on the
    result.

    For **time-accurate** runs it *does* enter the solution.  Dividing the
    destruction increment by $(1+\Delta t\,d)$ is a backward-Euler
    linearisation, so it formally reduces the turbulence **source** to
    first-order in time wherever $\Delta t\,d$ is not small — i.e. near walls,
    where $d\approx\beta\,\omega_\text{wall}\propto 1/y^2$ and
    $\Delta t\,d\gg 1$ even at a modest CFL.  Note that:

    - only the turbulence *destruction* is affected — the mean-flow equations
      and the turbulence production/diffusion/convection retain the full
      Runge–Kutta order;
    - the terms made implicit are exactly the *stiff, fast* ones, whose
      near-wall relaxation time is far shorter than any resolved URANS scale,
      so the physically relevant (resolved-scale) accuracy is essentially
      unchanged — the first-order error lives only in the unresolved fast
      transient.

    Disabling `point-implicit` for a URANS run does **not** recover accuracy
    for free: the explicit near-wall $\omega$ source then violates its
    stability limit and diverges.  To genuinely verify temporal convergence,
    set `point-implicit = .false.` **and** reduce $\Delta t$ until
    $\Delta t\,\beta\,\omega_\text{wall} < 1$ in the first off-wall cell
    (usually impractically small).  For production URANS the recommendation is
    to leave the treatment enabled.

The treatment is enabled by default and controlled by

```ini
[MOSE-Turbulence]
point-implicit = .true.   ; (default) point-implicit turbulence destruction
```

### Relation to under-relaxation

Point-implicit damping is a *physically weighted, automatic* form of
under-relaxation: the effective relaxation factor $1/(1+\Delta t\,d)$ is
large (≈1) where the source is mild and small only where the destruction is
stiff, and it scales with the local time step.  A constant global
under-relaxation factor on the turbulence update achieves a similar
stabilising effect but is blunter — it slows convergence everywhere, must be
hand-tuned per case/grid, and (if applied to the whole update rather than the
increment) can bias the transient.  Point-implicit treatment is therefore the
preferred technique; global under-relaxation remains a valid, simpler
fallback when a model's destruction Jacobian is not readily available.

---

## References

1. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for
   aerodynamic flows," AIAA-92-0439, 1992.
2. F. R. Menter, "Two-equation eddy-viscosity turbulence models for
   engineering applications," *AIAA J.*, 32(8), 1994.
3. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial
   experience with the SST turbulence model," in *Turbulence, Heat
   and Mass Transfer 4*, Begell House, 2003.
4. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries,
   2006.
5. P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models
   to rotation and curvature," *Aerosp. Sci. Technol.*, 1(5), 1997.
6. P. R. Spalart, "Strategies for turbulence modelling and simulations,"
   *Int. J. Heat Fluid Flow*, 21(3), 2000 (Quadratic Constitutive
   Relation, QCR2000).
7. C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain
   correlation of turbulence: an invariant dynamical systems approach,"
   *J. Fluid Mech.*, 227, 1991 (SSG model); B. E. Launder, G. J. Reece,
   W. Rodi, "Progress in the development of a Reynolds-stress turbulence
   closure," *J. Fluid Mech.*, 68(3), 1975 (LRR model).
