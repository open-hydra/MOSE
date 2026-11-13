# Turbulence Modelling

MOSE provides Reynolds-Averaged Navier–Stokes (RANS) closures for
turbulent flows.  The eddy-viscosity models
described in detail below range from the one-equation Spalart–Allmaras
model (with several variants) to the two-equation Menter SST and
Wilcox 2006 $k$–$\omega$ models; beyond them, the
[SSG/LRR-$\omega$ Reynolds-stress model](#ssglrr-reynolds-stress-model)
transports the Reynolds stresses themselves.  An algebraic non-linear
correction (QCR2000) is also available — see [Other models](#other-models).

---

## Boussinesq Hypothesis

All eddy-viscosity models assume that the Reynolds stress tensor is
proportional to the mean strain rate:

$$
\tau_{ij}^R = 2\,\mu_t\bigl(S_{ij} - \tfrac{1}{3}\,\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \tfrac{2}{3}\,\rho\,k\,\delta_{ij}
$$

where $S_{ij} = \tfrac{1}{2}(\partial v_i / \partial x_j + \partial v_j / \partial x_i)$
is the mean strain-rate tensor, $\mu_t$ is the eddy (turbulent) viscosity,
and $k$ is the turbulence kinetic energy.  The isotropic term
$-\tfrac{2}{3}\rho k\,\delta_{ij}$ is left out of the stress of the momentum
and energy equations unless the energy coupling is on, with which the
two-equation models (SST, Wilcox 2006) include it (see
[Energy coupling](turbulence-variants.md#energy-coupling)); SA has no $k$.

---

## Spalart–Allmaras (SA) — One-Equation Model

### Transport equation

The SA model solves a single transport equation for the modified
turbulent viscosity $\tilde\nu$. MOSE takes the model without the $f_{t2}$
term, SA-noft2 in the naming of the NASA Turbulence Modeling Resource, whose
standard SA has it, with $c_{t3} = 1.2$; SA-R and SA-RC apply their
corrections to this model (SA-noft2-R and SA-noft2-RC in that naming). The
equation is

$$
\frac{\partial(\rho\tilde\nu)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,\tilde\nu)
= \underbrace{c_{b1}\,\tilde{S}\,\rho\tilde\nu}_{\text{Production}}
+ \underbrace{\frac{1}{\sigma}\Bigl[
  \nabla\!\cdot\!\bigl((\mu + \rho\tilde\nu)\,\nabla\tilde\nu\bigr)
  + c_{b2}\,\rho\,|\nabla\tilde\nu|^2
  - (\nu + \tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu
\Bigr]}_{\text{Diffusion}}
- \underbrace{c_{w1}\,f_w\,\rho\!\left(\frac{\tilde\nu}{y}\right)^{\!2}}_{\text{Destruction}}
$$

This is the conservation form of the original equation, which has no density and holds for
compressible flows too: it is $\rho$ times the original equation plus $\tilde\nu$ times the
continuity equation (S. R. Allmaras, F. T. Johnson, P. R. Spalart, *Modifications and
Clarifications for the Implementation of the Spalart–Allmaras Turbulence Model*, ICCFD7-1902,
2012, eq. 9; NASA Turbulence Modeling Resource, Spalart–Allmaras page). Since
$\rho\,\nabla\cdot\bigl((\nu+\tilde\nu)\nabla\tilde\nu\bigr) = \nabla\cdot\bigl(\rho(\nu+\tilde\nu)\nabla\tilde\nu\bigr)
- (\nu+\tilde\nu)\nabla\rho\cdot\nabla\tilde\nu$, the last diffusion term is needed wherever density
and $\tilde\nu$ vary together (boundary layers at high Mach number, shocks): without it, as in
implementations that put the density inside the derivatives, the predictions of supersonic flows
change (NASA Turbulence Modeling Resource, Note 4).

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
\tilde{S} =
\begin{cases}
\Omega + \bar S, & \bar S \ge -c_2\,\Omega,\\[1ex]
\Omega + \dfrac{\Omega\,(c_2^2\,\Omega + c_3\,\bar S)}{(c_3 - 2c_2)\,\Omega - \bar S}, & \bar S < -c_2\,\Omega,
\end{cases}
\qquad
\bar S = \frac{\tilde\nu}{\kappa^2 y^2}\,f_{v2}(\chi), \quad c_2 = 0.7, \quad c_3 = 0.9,
$$

where $\Omega = \sqrt{2\,W_{ij}\,W_{ij}}$ is the vorticity magnitude, with the rotation tensor
$W_{ij} = \tfrac{1}{2}\left(\partial u_i/\partial x_j - \partial u_j/\partial x_i\right)$ (NASA Turbulence
Modeling Resource, Spalart–Allmaras page), and $y$ is the distance to the nearest wall.  This is
Note 1(c) of that page (from ICCFD7-1902), which replaces $\Omega + \bar S$ where $\bar S < -c_2\,\Omega$:
with it $\tilde S > 0.1\,\Omega$ wherever $\Omega > 0$, without clipping $\tilde S$.

Destruction function:

$$
r = \min\!\left(\frac{\tilde\nu}{\kappa^2 y^2 \tilde{S}},\ 10\right), \qquad
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

on the viscous walls, on the grain of a solid rocket motor that blows its products in through a wall,
and on the gas–surface interaction walls; the wall distance $y$ is measured from them and from the walls
coupled to another solver. A slip wall (in MOSE the same boundary as the symmetry plane) follows the
standard model: the ghost cell copies the cell, a zero normal derivative of $\tilde\nu$, and the
slip wall is not a wall for the distance («the free-slip wall points are not included in the search
when $d$ is computed», P. R. Spalart, S. R. Allmaras,
*La Recherche Aérospatiale* 1994(1), p. 21).

Allmaras, Johnson and Spalart (ICCFD7-1902, §3.3)
propose instead $\tilde\nu = 0$ on slip walls and their inclusion in the distance, as distinct
from symmetry planes; MOSE, whose one boundary serves both, does not make that change.

Negative values are clipped: after every stage MOSE keeps $\rho\tilde\nu \ge 10^{-10}$ (the
realizability floor) and counts the cell updates in which a stage made it negative
([Time-accurate runs](turbulence-numerics.md#time-accurate-runs-explicit-sources)). ICCFD7-1902 (§3.2) replaces the
clipping with the negative SA model, since clipping «prevents the convergence of discrete PDE
residuals»; MOSE does not implement SA-neg.

### Rough walls — SA-rough

A wall face with a sand-grain roughness height $k_s > 0$ (see the
[wall BCs](../user/boundary-conditions.md#viscous-walls)) switches on the
Boeing extension of Aupoix & Spalart (2003). Each cell takes the $k_s$ of its
nearest wall face, and the model changes in three places.

The wall is displaced below the surface by $d_0 = 0.03\,k_s$. Every $y$ above
becomes

$$
d = y + 0.03\,k_s
$$

The damping variable carries the roughness, and $f_{v2}$ is written without
$\chi$ alone, which no longer equals $\tilde\nu/\nu$:

$$
\chi = \frac{\tilde\nu}{\nu} + c_{R1}\,\frac{k_s}{d}, \qquad
f_{v2} = 1 - \frac{\tilde\nu}{\nu + \tilde\nu\,f_{v1}(\chi)}, \qquad c_{R1} = 0.5
$$

The wall condition becomes a Robin condition, so $\tilde\nu$ is nonzero at the wall:

$$
\left.\frac{\partial\tilde\nu}{\partial n}\right|_\text{wall} = \frac{\tilde\nu_\text{wall}}{d_0}
\quad\Longrightarrow\quad
\tilde\nu_\text{wall} = \tilde\nu_1\,\frac{d_0}{y_1 + d_0}
$$

where $\tilde\nu_1$ and $y_1$ are the value and wall distance of the boundary
cell. The eddy viscosity at the wall, $\mu_{t,\text{wall}} = \rho_\text{wall}\,\tilde\nu_\text{wall}\,f_{v1}(\chi_\text{wall})$
with $d = d_0$, then enters the wall shear stress and, on isothermal walls, the
wall heat flux through $\mu_t c_p/\mathrm{Pr}_t$.

With $k_s = 0$ every expression reduces to the smooth model, and MOSE takes the
smooth code path bit for bit.

!!! warning "Scope"
    Roughness is modelled for SA only. With SST, Wilcox 2006, SSG–LRR or a
    laminar run, MOSE prints a warning at start-up and treats every wall as
    smooth. Gas–surface interaction walls are always smooth.

---

### SA Variants

The model has three variants: a compressibility correction (SAcomp), a rotation correction (SA-R) and the Spalart–Shur
rotation-curvature correction (SA-RC). Their definitions are in
[Turbulence Model Variants and Corrections](turbulence-variants.md#sa-variants).

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear Stress Transport model blends a $k$–$\omega$ formulation
(near walls) with a $k$–$\varepsilon$-like behaviour (in the
freestream) using blending functions.

MOSE implements the **SST-2003** version of the model (Menter, Kuntz &
Langtry [3]; the original model is Menter 1994 [2]) as specified by the NASA
Turbulence Modeling Resource [8]: eddy viscosity built from the strain
invariant $S$, production limiter $10\,\beta^\ast\rho\,\omega k$ in both
equations, and the 2003 constants.  Two choices set the remaining
details, the form of the production and the energy coupling; the last two columns give the TMR name of each variant:

| Production | $P$ | Without energy coupling (default) | With energy coupling |
|---|---|---|---|
| incompressible (default) | $\mu_t S^2$ | SST-2003m | SST-2003s |
| compressible | $\tau_{ij}\,\partial u_i/\partial x_j$ (exact), arranged as in OpenFOAM's `kOmegaSST` | none | none |

The suffix records the approximations [8]: without the
[energy coupling](turbulence-variants.md#energy-coupling) the
$-\tfrac23\rho k\,\delta_{ij}$ part of $\tau_{ij}$ is left out of the
momentum and energy equations ("m", "e"); "m" and "s" approximate the production
by $\mu_t S^2$, which TMR describes as exact for incompressible flow and
"typically considered a very good approximation, except perhaps for very high
Mach number flows"; SST-2003e and SST-2003 take the exact production and limit
all of it.  With the compressible production the limiter acts on the deviatoric part of $P$
only, and the $\omega$ equation takes the dilatation term
$-\tfrac23\gamma\rho\omega\,\nabla\!\cdot\!\mathbf{v}$, as in OpenFOAM: these two
variants have no TMR name (see [below](turbulence-variants.md#compressible-production)). SST-RC adds the
[rotation/curvature correction](turbulence-variants.md#rotationcurvature-correction-sst-rc).

### Transport equations

$$
\frac{\partial(\rho k)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,k)
= P_k - \beta^\ast\,\rho\,\omega\,k
+ \nabla\!\cdot\!\bigl[(\mu + \sigma_{k}\,\mu_t)\,\nabla k\bigr]
$$

$$
\frac{\partial(\rho\omega)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,\omega)
= \frac{\gamma\,\rho}{\mu_t}\,P_k - \beta\,\rho\,\omega^2
+ \nabla\!\cdot\!\bigl[(\mu + \sigma_{\omega}\,\mu_t)\,\nabla\omega\bigr]
+ 2(1 - F_1)\,\frac{\rho\,\sigma_{\omega 2}}{\omega}\,\nabla k\!\cdot\!\nabla\omega
$$

The last term is the **cross-diffusion** term, active only in the
freestream ($F_1 \to 0$).

### Production limiter

$$
P_k = \min\!\bigl(\mu_t\,S^2,\;10\,\beta^\ast\,\rho\,\omega\,k\bigr)
$$

This prevents unbounded growth of $k$ in stagnation regions.

!!! note "$\omega$-production uses the *limited* $P_k$"
    Per the NASA Turbulence-Modelling-Resource SST-2003 specification [8], the
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

### SST Variants and Options

Beyond the model above, SST has four options: the compressible form of the production, a compressibility correction of the
dissipation, a rotation-curvature correction (SST-RC), and the coupling of $k$ with the energy of the mean flow, which the default
drops. Their definitions are in [Turbulence Model Variants and Corrections](turbulence-variants.md).

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
              \frac{500\,\nu}{\omega\,y^2}\right),\;
  \frac{4\,\rho\,\sigma_{\omega 2}\,k}{CD_{k\omega}\,y^2}\right), \qquad
CD_{k\omega} = \max\!\bigl(2\,\rho\,\sigma_{\omega 2}\,\omega^{-1}\,\nabla k\!\cdot\!\nabla\omega,\;10^{-10}\bigr)
$$

$$
F_2 = \tanh\!\bigl(\arg_2^2\bigr), \qquad
\arg_2 = \max\!\left(\frac{2\sqrt{k}}{0.09\,\omega\,y},\;
                      \frac{500\,\nu}{\omega\,y^2}\right)
$$

### Wall boundary conditions

$$
k_\text{wall} = 0, \qquad
\omega_\text{wall} = C\,\frac{\nu_w}{y_c^2}
$$

with $y_c$ the wall distance of the wall-adjacent cell centre and $\nu_w$ the
laminar viscosity at the wall.  $C$ takes one of two values:

- practical (default): $C = 10 \times 6/\beta_1 = 800$, Menter's
  over-specification of the near-wall limit, robust on meshes that do not
  resolve the viscous sublayer;
- asymptotic: $C = 6/\beta_1 = 80$, the viscous-sublayer solution
  $\omega = 6\nu/(\beta_1 y^2)$ evaluated at $y_c$, for wall-resolved meshes.

The wall flux of the walls that blow (pyrolysis and surface-reaction walls while they ablate, a solid
propellant grain while it burns) takes $C = 800$ with either value.

MOSE imposes $\omega_\text{wall}$ on the wall **face** (ghost value
$2\,\omega_\text{wall} - \omega_P$) and solves the $\omega$ equation in the
wall-adjacent cell.  Codes that instead **fix the wall-adjacent cell** to
$6\nu_w/(\beta_1 y_c^2)$ (OpenFOAM's `omegaWallFunction`) are not equivalent
even with the same $C$: since the sublayer solution is infinite at the wall, the
face value $80\,\nu_w/y_c^2$ yields $\omega \approx 6\nu/[\beta_1 (y + y_c)^2]$,
about a quarter of the sublayer value in the first cell.  In the
[SWBLI case](../vv/2D-swbli.md) that alone moves the separation point by about
1 mm, so a code-to-code comparison must apply the condition in the same way.

Either value grows as $1/y_c^2$ under grid refinement, which is the origin of the
near-wall stiffness handled by the [numerical treatment of the source
terms](turbulence-numerics.md#numerical-treatment-of-the-source-terms).

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

The Wilcox (2006) model [4, 12] with the exact production
$P = \tau_{ij}\,\partial u_i/\partial x_j$ and the optional TMR $k$-equation
limiter (below): TMR's Wilcox2006-klim, in its "e" form (Wilcox2006-klim-e) unless
the [energy coupling](turbulence-variants.md#energy-coupling) is on [8].

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
  C_\text{lim}\,\sqrt{\frac{2\,\bar S_{ij}\,\bar S_{ij}}{\beta^\ast}}\right), \qquad
\bar S_{ij} = S_{ij} - \tfrac{1}{3}\,\frac{\partial u_k}{\partial x_k}\,\delta_{ij}
$$

The limiter $C_\text{lim}$ prevents excessive eddy viscosity in regions
where $\omega$ is small relative to the strain rate.

### Destruction with stress-limiter correction

The $\omega$-destruction coefficient is modified by Wilcox's
vortex-stretching function $f_\beta$:

$$
\beta = \beta_0\,f_\beta, \qquad
X_\omega = \frac{|W_{ij}\,W_{jk}\,\hat S_{ki}|}{(\beta^\ast\omega)^3},
\qquad
f_\beta = \frac{1 + 85\,X_\omega}{1 + 100\,X_\omega}
$$

with $\hat S_{ki} = S_{ki} - \tfrac{1}{2}\,(\partial u_m/\partial x_m)\,\delta_{ki}$ (NASA Turbulence Modeling
Resource, Wilcox 2006 page): off the diagonal $\hat S = S$. On that page the formula of $\chi_\omega$ reads
$\Omega_{ij}\Omega_{ij}\hat S_{ki}$, whose indices do not close; MOSE takes $W_{ij}W_{jk}\hat S_{ki}$, the form of the Wilcox
1998 model on the same page ($\Omega_{ij}\Omega_{jk}S_{ki}$) and the only one for which the statement of the page that
$\chi_\omega$ is zero in 2-D holds. In 2-D flows the product $W_{ij}W_{jk}\hat S_{ki}$ is zero and $f_\beta = 1$.

### Wall boundary conditions

$k_\text{wall} = 0$. With the asymptotic condition the wall takes the near-wall limit of the $\omega$ equation of the
model, $\omega_\text{wall} = 6\nu/(\beta_0 y^2)$ with $\beta_0 = 0.0708$ (NASA Turbulence Modeling Resource, Wilcox 2006
page), that is $84.7\,\nu/y^2$, $y$ the distance of the centre of the wall cell from the wall. With the practical one, the default,
it takes Menter's condition of the [SST section](#wall-boundary-conditions), $800\,\nu/y^2$ with the $\beta_1$ of SST,
as the same page advises for this model. The wall flux of the walls that blow (pyrolysis and surface-reaction walls while
they ablate, a solid propellant grain while it burns) takes $800\,\nu/y^2$ with either value.

### Cross-diffusion *(conditional)*

Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega > 0$:

$$
\text{CD} = \frac{\sigma_d\,\rho}{\omega}\,
  (\nabla k\!\cdot\!\nabla\omega), \qquad
\sigma_d = \tfrac{1}{8} \ \text{if } \nabla k\!\cdot\!\nabla\omega > 0, \qquad
\sigma_d = 0 \ \text{if } \nabla k\!\cdot\!\nabla\omega \le 0
$$

### Production limiter

More permissive than SST, and in the $k$ equation only (the $\omega$ production
uses the unlimited $P$), as in TMR's Wilcox2006-klim:

$$
P_k = \min\!\bigl(\tau_{ij}\,\partial u_i/\partial x_j,\;20\,\beta^\ast\,\rho\,\omega\,k\bigr)
$$

---

## SSG/LRR-$\omega$ — Reynolds-Stress Model { #ssglrr-reynolds-stress-model }

The SSG/LRR-$\omega$ model of Eisfeld, Rumsey & Togiti [15] does not use the
Boussinesq hypothesis: it transports the six Reynolds stresses
$\rho R_{ij} = \overline{\rho\,u_i''u_j''}$ and the specific dissipation rate
$\omega$, seven equations in all.  The pressure–strain model and the $\omega$
equation blend the Launder–Reece–Rodi (LRR) coefficients near walls with the
Speziale–Sarkar–Gatski (SSG) ones away from them [7].  Being anisotropy-resolving, the
model captures secondary flows and the effects of streamline curvature and
rotation that the eddy-viscosity models miss, at the cost of the additional
equations.

MOSE follows the specification of the NASA Turbulence Modeling Resource [16],
in four variants:

| TMR name | Diffusion of the stresses | Length-scale correction |
|---|---|---|
| SSGLRR-RSM-w2012 | generalized gradient | no |
| SSGLRR-RSM-w2012-SD | simple (eddy viscosity) | no |
| SSGLRR-RSM-w2019 | generalized gradient | yes |
| SSGLRR-RSM-w2019-SD | simple (eddy viscosity) | yes |

### Transport equations

$$
\frac{\partial(\rho R_{ij})}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,R_{ij})
= \rho P_{ij} + \rho\,\Pi_{ij} - \rho\,\varepsilon_{ij} + \rho D_{ij}
$$

with $k = \tfrac12 R_{ii}$ and $\varepsilon = C_\mu\,k\,\omega$, $C_\mu = 0.09$.
The production is exact and the dissipation isotropic:

$$
\rho P_{ij} = -\rho R_{ik}\,\frac{\partial u_j}{\partial x_k}
              -\rho R_{jk}\,\frac{\partial u_i}{\partial x_k},
\qquad
\rho\,\varepsilon_{ij} = \tfrac23\,\rho\,\varepsilon\,\delta_{ij}
$$

The pressure–strain correlation is written with the anisotropy tensor
$a_{ij} = R_{ij}/k - \tfrac23\delta_{ij}$, the strain rate $S_{ij}$, its
deviator $S^\ast_{ij} = S_{ij} - \tfrac13 S_{kk}\delta_{ij}$ and the rotation
tensor $W_{ij} = \tfrac12(\partial u_i/\partial x_j - \partial u_j/\partial x_i)$:

$$
\begin{aligned}
\rho\,\Pi_{ij} ={}& -\bigl(C_1\,\rho\,\varepsilon + \tfrac12 C_1^\ast\,\rho P_{kk}\bigr)\,a_{ij}
 + C_2\,\rho\,\varepsilon\,\bigl(a_{ik}a_{kj} - \tfrac13 a_{kl}a_{kl}\,\delta_{ij}\bigr) \\
 & + \bigl(C_3 - C_3^\ast\sqrt{a_{kl}a_{kl}}\bigr)\,\rho k\,S^\ast_{ij}
 + C_4\,\rho k\,\bigl(a_{ik}S_{jk} + a_{jk}S_{ik} - \tfrac23 a_{kl}S_{kl}\,\delta_{ij}\bigr) \\
 & + C_5\,\rho k\,\bigl(a_{ik}W_{jk} + a_{jk}W_{ik}\bigr)
\end{aligned}
$$

The diffusion is the generalized gradient model of Daly & Harlow [17], or, with
the -SD suffix, a simple gradient model built on $\mu_t = \rho k/\omega$:

$$
\rho D_{ij} = \frac{\partial}{\partial x_k}\!\left[
  \left(\mu\,\delta_{kl} + D\,\frac{\rho R_{kl}}{C_\mu\,\omega}\right)
  \frac{\partial R_{ij}}{\partial x_l}\right]
\qquad\text{or}\qquad
\rho D_{ij} = \frac{\partial}{\partial x_k}\!\left[
  \left(\mu + \frac{D_\text{SD}}{C_\mu}\,\mu_t\right)
  \frac{\partial R_{ij}}{\partial x_k}\right]
$$

with $D_\text{SD} = 0.5\,C_\mu F_1 + \tfrac23\,0.22\,(1 - F_1)$.  The length-scale
equation is

$$
\frac{\partial(\rho\omega)}{\partial t}
+ \nabla\!\cdot\!(\rho\,\mathbf{v}\,\omega)
= \alpha_\omega\,\frac{\omega}{k}\,\frac{\rho P_{kk}}{2} - \beta_\omega\,\rho\,\omega^2
+ \nabla\!\cdot\!\Bigl[\Bigl(\mu + \sigma_\omega\frac{\rho k}{\omega}\Bigr)\nabla\omega\Bigr]
+ \sigma_d\,\frac{\rho}{\omega}\,\max\!\bigl(\nabla k\!\cdot\!\nabla\omega,\;0\bigr)
$$

### Model constants

Every coefficient $\phi$ is blended as
$\phi = F_1\,\phi^{(\omega)} + (1 - F_1)\,\phi^{(\varepsilon)}$ with Menter's function

$$
F_1 = \tanh\!\bigl(\zeta^4\bigr), \qquad
\zeta = \min\!\left[
  \max\!\left(\frac{\sqrt{k}}{C_\mu\,\omega\,d},\;\frac{500\,\mu}{\rho\,\omega\,d^2}\right),\;
  \frac{4\,\sigma_\omega^{(\varepsilon)}\rho k}{\text{CD}\,d^2}\right],
\qquad
\text{CD} = \sigma_d^{(\varepsilon)}\,\frac{\rho}{\omega}\,\max\!\bigl(\nabla k\!\cdot\!\nabla\omega,\;0\bigr)
$$

where $d$ is the distance to the nearest wall.

| Set | $\alpha_\omega$ | $\beta_\omega$ | $\sigma_\omega$ | $\sigma_d$ | $C_1$ | $C_1^\ast$ | $C_2$ | $C_3$ | $C_3^\ast$ | $C_4$ | $C_5$ | $D$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $\omega$ (LRR, near wall) | 0.5556 | 0.075 | 0.5 | 0 | 1.8 | 0 | 0 | 0.8 | 0 | $\tfrac{9C_2^\text{LRR}+6}{11}$ | $\tfrac{-7C_2^\text{LRR}+10}{11}$ | $0.75\,C_\mu$ |
| $\varepsilon$ (SSG, outer) | 0.44 | 0.0828 | 0.856 | 1.712 | 1.7 | 0.9 | 1.05 | 0.8 | 0.65 | 0.625 | 0.2 | 0.22 |

with $C_2^\text{LRR} = 0.52$.

### Length-scale correction (w2019)

The w2019 variants [18] multiply the destruction of $\omega$ by
$1 - F^\text{(LSC)}$, which removes it where the turbulent length scale exceeds
that of the log layer, near stagnation and reattachment points:

$$
F^\text{(LSC)} = \tfrac12\Bigl\{1 + \tanh\bigl[A\,(\chi - \chi_T)\bigr]\Bigr\},
\qquad
\chi = \max\!\bigl[(r - 1)\,r^2,\;0\bigr],
\qquad
r = \frac{\sqrt{k}}{C_\mu^{1/4}\,\kappa\,\omega\,d}
$$

with $A = 31$, $\chi_T = 1$ and $\kappa = 0.41$.  In an attached boundary layer
$r \le 1$ and the correction is inactive.

### Coupling with the mean flow

The momentum and energy equations carry the modelled stress directly,

$$
\tau_{ij} = 2\mu\bigl(S_{ij} - \tfrac13\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \rho R_{ij},
$$

and the turbulent heat and species fluxes are gradient models with
$\mu_t = \rho k/\omega$ and the turbulent Prandtl and Schmidt numbers $Pr_t$ and
$Sc_t$ (TMR uses $Pr_t = 0.9$).  The same $\mu_t$ enters the time-step limit.

The energy is treated as with the [energy coupling](turbulence-variants.md#energy-coupling),
which the Reynolds-stress model always has: MOSE advances the mean-flow energy
$E_m = e + \tfrac12 u_iu_i$, whose equation is the total-energy equation of
TMR [11] minus the $k$ equation (half the trace of the stress equations).  The
diffusion of $k$ cancels in the difference, and what is left is the exchange
with the turbulence,

$$
\frac{\partial(\rho E_m)}{\partial t} + \dots = \dots
- \bigl(\tfrac12\rho P_{kk} - \rho\,\varepsilon\bigr),
$$

so that $E_m + k$ is conserved: the work of the Reynolds stresses in the energy
flux contains the production of $k$, which is returned to the turbulence here,
and the dissipated $k$ comes back as heat.

### Boundary conditions

**Walls.**  $R_{ij} = 0$ and
$\omega_\text{wall} = C\,\nu_w/y_c^2$, imposed on the wall face exactly as for
[SST](#wall-boundary-conditions), with the same two values of $C$; the default,
practical ($C = 800$), is the condition of TMR.  The wall stress is the
molecular one.  On blowing (ablating) walls the stresses are zero as well.

**Symmetry planes.**  The ghost cell holds the mirror image of the stress
tensor, $\mathbf{Q}\,\mathbf{R}\,\mathbf{Q}$ with
$\mathbf{Q} = \mathbf{I} - 2\,\mathbf{n}\mathbf{n}^T$, and the face their mean,
which has no shear stress on the plane.  For a plane normal to an axis this is
the condition of TMR (the two shear stresses carrying that index vanish on the
plane, the other stresses have zero gradient); the plane may have any
orientation.

**Inlets and initial conditions.**  ATLAS gives the three normal stresses one
value (isotropic turbulence, $R_{11} = R_{22} = R_{33} = \tfrac23 k$) and sets
the shear stresses to $10^{-8}$; with $\omega$ this fixes the state.  At an
inlet the specific values are multiplied by the inlet density, as for the
two-equation models ([Initial Conditions](../user/initial-conditions.md#required-variables)
gives the input keys).  The freestream values follow
the [same rule](#freestream-and-inlet-values-k-omega-tildenu) as for the
two-equation models, with $R_{ii} = \tfrac23 k_\infty$.

**Inflow, outflow and extrapolation faces.**  They carry no diffusive flux, but
the Reynolds stress is not one: it acts on these faces as the pressure does, and
its flux in the momentum and energy equations, $\rho R_{ij}n_j$ and
$\rho R_{ij}n_j u_i$, is taken from the boundary cell.  Without it the pressure
jumps by $\rho R_{nn}$ in the cells next to an open boundary.

**Axisymmetric and rotationally periodic faces.**  The stress tensor turns
with the ghost cell, as the velocity does, in the inviscid flux of these faces
and in their viscous flux, whose
stress includes $-\rho R_{ij}$ of the face: on the two side faces of an
axisymmetric wedge the azimuthal normal stress $\rho R_{\theta\theta}$ thus
enters the radial momentum balance with the viscous hoop stress.

### Realizability

After every update the normal stresses are kept positive, the shear stresses
within $|R_{ij}| \le \sqrt{R_{ii}R_{jj}}$ (no sum), and $\omega$ positive.

### Rotating frame

In a [rotating frame](rotating-frame.md) of angular velocity $\boldsymbol\Omega$
the stress equations gain the Coriolis term

$$
\rho\,G_{ij} = -2\,\Omega_k\,\bigl(\epsilon_{ikm}\,\rho R_{jm} + \epsilon_{jkm}\,\rho R_{im}\bigr),
$$

which redistributes the stresses without producing $k$ ($G_{kk} = 0$), and the
pressure–strain model uses the absolute rotation tensor
$W_{ij} + \epsilon_{mji}\,\Omega_m$ [19, 20].

---

## Other models

### QCR2000 (non-linear constitutive correction)

The **Quadratic Constitutive Relation** (Spalart, 2000) is an algebraic,
non-linear correction to the Boussinesq stress rather than a standalone
model. It adds a quadratic dependence on the rotation tensor to the modelled
stress, improving the prediction of anisotropy-driven secondary flows
(corner flows, square ducts) while reusing the eddy viscosity of the
underlying model. MOSE does not offer it as a model choice; its stress is
checked against the definition by the tests.

---

## General RANS Features

The wall values of the turbulence variables are imposed through the ghost cells, which take $Z_\text{ghost} = 2Z_\text{wall} - Z_\text{cell}$.
The source terms take the central differences of the values along the three index directions, mapped with a cell metric. The wall
distance of the models is the minimum distance from the cell centre to the wall surface, as the NASA Turbulence Modeling Resource
requires. The details are in [Turbulence: Numerical Treatment](turbulence-numerics.md).

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
⇒ more near-wall dissipation and a more robust start-up).  The TMR SST
freestream used in the [SWBLI case](../vv/2D-swbli.md) is
$k_\infty = 9\times10^{-9}\,a_\infty^2$, $\omega_\infty = 10^{-6}\,a_\infty^2/\nu_\infty$
($\mu_t/\mu = 0.009$): at $p=4000$ Pa, $T=68.3$ K
($\rho_\infty=0.2041$ kg m⁻³, $a_\infty=165.7$ m s⁻¹,
$\nu_\infty=2.257\times10^{-5}$ m² s⁻¹ with Sutherland's law) this is
$k_\infty = 2.47\times10^{-4}$ m² s⁻², $\omega_\infty = 1216$ s⁻¹.
Values of $\omega_\infty$ orders of magnitude below the $\mu_t/\mu$ estimate are
a frequent cause of $k$ runaway.

**Specific and density-weighted values.**  The values above are specific, as TMR and
OpenFOAM give them; $\tilde\nu$ is the variable the Spalart–Allmaras model transports,
not the eddy viscosity $\mu_t = \rho\tilde\nu f_{v1}$.  MOSE solves for the
density-weighted variables $\rho\tilde\nu$, $\rho k$, $\rho\omega$: the specific values of an
inlet are multiplied by the density of the inlet state, and those of the initial condition by
the density of each cell.  The input keys and file variables are in
[Initial Conditions](../user/initial-conditions.md#required-variables).

---

## Numerical treatment of the source terms

The source terms are added to the residual and advanced with the same explicit Runge–Kutta stages as the mean flow. Their destruction
parts are stiff near walls. A steady run divides the increment of each turbulence equation by $1 + \Delta t\,d$ (the point-implicit
factor of Bussing and Murman). A time-accurate run keeps the sources explicit, since the factor would make the equations first-order
accurate in time, and bounds their stiffness in the time step. The derivations are in
[Turbulence: Numerical Treatment](turbulence-numerics.md#numerical-treatment-of-the-source-terms).

!!! info "Going deeper"
    - [Turbulence Model Variants and Corrections](turbulence-variants.md): SA variants; SST compressible production, dilatation dissipation, rotation-curvature correction and energy coupling.
    - [Turbulence: Numerical Treatment](turbulence-numerics.md): wall ghost cells, cell gradients, wall distance, and the advance of the source terms.

## References

[1] P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for aerodynamic flows," AIAA-92-0439, 1992.

[2] F. R. Menter, "Two-equation eddy-viscosity turbulence models for engineering applications," *AIAA J.*, 32(8), pp. 1598–1605, 1994, doi:10.2514/3.12149.

[3] F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial experience with the SST turbulence model," in *Turbulence, Heat and Mass Transfer 4*, Begell House, pp. 625–632, 2003.

[4] D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.

[5] P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models to rotation and curvature," *Aerosp. Sci. Technol.*, 1(5), 1997.

[6] P. R. Spalart, "Strategies for turbulence modelling and simulations," *Int. J. Heat Fluid Flow*, 21(3), 2000 (Quadratic Constitutive Relation, QCR2000).

[7] C. G. Speziale, S. Sarkar, T. B. Gatski, "Modelling the pressure–strain correlation of turbulence: an invariant dynamical systems approach," *J. Fluid Mech.*, 227, 1991 (SSG model); B. E. Launder, G. J. Reece, W. Rodi, "Progress in the development of a Reynolds-stress turbulence closure," *J. Fluid Mech.*, 68(3), 1975 (LRR model).

[8] NASA Langley Turbulence Modeling Resource, "The Menter Shear Stress Transport turbulence model," <https://tmbwg.github.io/turbmodels/sst.html> (formerly turbmodels.larc.nasa.gov/sst.html): definitions and naming of the SST variants (SST-2003, SST-2003m/s/e, SST-RC) and of Wilcox2006-klim.

[11] NASA Langley Turbulence Modeling Resource, "Implementing turbulence models into the compressible RANS equations," <https://tmbwg.github.io/turbmodels/implementrans.html>.

[12] D. C. Wilcox, "Formulation of the k-ω turbulence model revisited," *AIAA J.*, 46(11), pp. 2823–2838, 2008, doi:10.2514/1.36541.

[13] NASA Langley Research Center, Turbulence Modeling Resource, "The Spalart–Allmaras Turbulence Model", https://turbmodels.larc.nasa.gov/spalart.html (TMBWG/turbmodels commit 20a39f5, consulted 2026-09-29): minimum distance to the wall.

[15] B. Eisfeld, C. Rumsey, V. Togiti, "Verification and validation of a second-moment-closure model," *AIAA J.*, 54(5), pp. 1524–1541, 2016, doi:10.2514/1.J054718; erratum *AIAA J.*, 54(9), p. 2926, 2016, doi:10.2514/1.J055336.

[16] NASA Langley Turbulence Modeling Resource, "SSG/LRR full Reynolds stress model," <https://tmbwg.github.io/turbmodels/rsm-ssglrr.html>: definitions and naming of the variants (SSGLRR-RSM-w2012, -SD, -w2019).

[17] B. J. Daly, F. H. Harlow, "Transport equations in turbulence," *Phys. Fluids*, 13, pp. 2634–2649, 1970, doi:10.1063/1.1692845.

[18] B. Eisfeld, C. Rumsey, "Length-scale correction for Reynolds-stress modeling," *AIAA J.*, 58(4), pp. 1518–1528, 2020, doi:10.2514/1.J058858.

[19] B. E. Launder, D. P. Tselepidakis, B. A. Younis, "A second-moment closure study of rotating channel flow," *J. Fluid Mech.*, 183, pp. 63–75, 1987.

[20] C. G. Speziale, "Turbulence modeling in noninertial frames of reference," *Theor. Comput. Fluid Dyn.*, 1, pp. 3–19, 1989.
