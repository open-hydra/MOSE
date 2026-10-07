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
\tilde{S} = \Omega + \frac{\tilde\nu}{\kappa^2 y^2}\,f_{v2}(\chi)
$$

where $\Omega = \sqrt{2\,W_{ij}\,W_{ij}}$ is the vorticity magnitude, with the rotation tensor
$W_{ij} = \tfrac{1}{2}\left(\partial u_i/\partial x_j - \partial u_j/\partial x_i\right)$ (NASA Turbulence
Modeling Resource, Spalart–Allmaras page), and $y$ is the distance to the nearest wall.

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

on the viscous walls (types `301`, `302`, `503`–`505`); the wall distance $y$ is measured from
them and from the walls coupled to another solver (`103`, `104`). A slip wall (type `300`, which
is also the symmetry plane) follows the standard model: the ghost cell copies the cell, a zero
normal derivative of $\tilde\nu$, and the slip wall is not a wall for the distance («the free-slip
wall points are not included in the search when $d$ is computed», P. R. Spalart, S. R. Allmaras,
*La Recherche Aérospatiale* 1994(1), p. 21). Allmaras, Johnson and Spalart (ICCFD7-1902, §3.3)
propose instead $\tilde\nu = 0$ on slip walls and their inclusion in the distance, as distinct
from symmetry planes; MOSE, whose type `300` serves both, does not make that change.

Negative values are clipped: after every stage MOSE keeps $\rho\tilde\nu \ge 10^{-10}$ (the
realizability floor) and counts the cell updates in which a stage made it negative
([Time-accurate runs](#time-accurate-runs-explicit-sources)). ICCFD7-1902 (§3.2) replaces the
clipping with the negative SA model, since clipping «prevents the convergence of discrete PDE
residuals»; MOSE does not implement SA-neg.

### Rough walls — SA-rough

A wall face with a sand-grain roughness height $k_s > 0$ (the `k_rough` field of
the [wall BCs](../user/boundary-conditions.md#viscous-walls)) switches on the
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

#### SAcomp — Compressibility Correction

Paciorri–Sabetta correction that scales the production term with a
function of the turbulent stress ratio:

$$
S_\tau = \frac{\omega\,\tilde\nu\,f_{v1}}{a^2}
$$

Activated when compressibility effects on turbulence are significant
(high-speed boundary layers, mixing layers).

!!! warning "Stopped at set-up"
    A run with `turbulence = SAcomp` stops at set-up. The correction takes
    $c$ from $c^2 f(c) = 25\,\Omega\,\tilde\nu f_{v1}/a$, with
    $f(c) = 0.44/(1 + 14c^5) + 0.56$, and multiplies the production by
    $f_1(c)\,f(c)$, with $f_1(c) = 0.6/(1 + 9c^6) + 0.4$: the left side is a
    number and the right side a speed (the code divides by $a$, not by $a^2$
    as written above), so the factor depends on the units. For
    $\Omega = 10^3$ 1/s, $\nu_t = 10^{-3}$ m²/s and $a = 340$ m/s it is 0.99
    with lengths in metres and 0.22 with lengths in millimetres; and since
    $a$ only divides, at $\Omega = 10^5$ 1/s it cuts the production by 78 %
    whatever the Mach number. Until the correction is checked against its
    source (the article the code cites, doi:10.2514/2.3967) it is not
    available: use `SA` or another model.

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
\tilde{r} = \frac{2\,\omega_{ik}\,S_{jk}}{D^4}\left(\frac{DS_{ij}}{Dt}
  + (\varepsilon_{imn}S_{jn} + \varepsilon_{jmn}S_{in})\,\Omega'_m\right), \qquad
D^2 = \tfrac{1}{2}\,(S^2 + \Omega^2)
$$

where $\omega_{ik}$ is the rotation tensor of the absolute vorticity and $\Omega'$ the rotation rate
of the frame, zero in a fixed frame (NASA Turbulence Modeling Resource, SA-RC). MOSE computes the
factor $2r^\ast/(1 + r^\ast)$ of $f_{r1}$ below as $2S/(S + \Omega)$: where the vorticity vanishes it
takes its limit, 2. Where strain and vorticity are both below $10^{-30}\ \mathrm{s^{-1}}$, in a fluid
at rest in a fixed frame or a uniform flow, $r^\ast$ is not defined (the velocity gradient there is the
$10^{-40}$ that MOSE adds to the velocity differences, times the metric of the cell): MOSE takes the
value of simple shear, 1, with $\tilde r = 0$, so that the correction does not act and the production
there is that of SA.

The correction factor:

$$
f_{r1} = (1 + c_{r1})\,\frac{2\,r^\ast}{1 + r^\ast}\,
  \bigl(1 - c_{r3}\arctan(c_{r2}\,\tilde{r})\bigr) - c_{r1}
$$

with constants $c_{r1} = 1$, $c_{r2} = 12$, $c_{r3} = 1$ (NASA Turbulence Modeling Resource,
SA-RC). Production is multiplied by $f_{r1}$, which is not bounded and can be negative.

MOSE evaluates the Lagrangian derivative $DS_{ij}/Dt = \partial S_{ij}/\partial t + u_k\,\partial
S_{ij}/\partial x_k$ with its convective part only. Without the time term the correction holds for
steady flows only: «Without the time term, the model is not correct in a time-dependent flowfield»
(NASA Turbulence Modeling Resource, Spalart–Allmaras page, SA-RC). Time-accurate runs with the
correction print a warning at start-up.

MOSE takes $\partial S_{ij}/\partial x_k$ with central differences of the velocity gradients of the
cells, and the ghost cells of every face of a block hold a copy of the gradient of the cell
(`Compute_Velocity_Gradient`). Next to a face between blocks the central difference is then
$(S_{c+1} - S_c)/2$, half the one-sided difference: for $S_{ij}$ linear across the face, its
derivative across the face in the cell next to it, and the part of $DS_{ij}/Dt$ that comes from it,
is half the true value when the mesh is divided there, and whole when it is not, so the result
depends on how the mesh is divided into blocks (and among MPI ranks).
On a column with a shear flow whose strain rate varies across the face, the right-hand side of
$\rho\tilde\nu$ of the column divided in two blocks differs from that of the undivided column in
the two columns of cells next to the face, and agrees to round-off elsewhere. A run with SA-RC on a
mesh of more than one block therefore stops
at setup with a message; before, it ran with that dependence. Exchanging the velocity gradients
across the faces between blocks would remove the limit, and is left for a later version. At the
faces of the physical boundaries (walls, inflow, outflow) the copy is a choice of the model.

---

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear Stress Transport model blends a $k$–$\omega$ formulation
(near walls) with a $k$–$\varepsilon$-like behaviour (in the
freestream) using blending functions.

MOSE implements the **SST-2003** version of the model (Menter, Kuntz &
Langtry [3]; the original model is Menter 1994 [2]) as specified by the NASA
Turbulence Modeling Resource [8]: eddy viscosity built from the strain
invariant $S$, production limiter $10\,\beta^\ast\rho\,\omega k$ in both
equations, and the 2003 constants.  Two options select the remaining
details, and together they span TMR's naming of the variants:

| `sst-production` | Production $P$ | `k-energy-coupling = false` (default) | `k-energy-coupling = true` |
|---|---|---|---|
| `incompressible` (default) | $\mu_t S^2$ | SST-2003m | SST-2003s |
| `compressible` | $\tau_{ij}\,\partial u_i/\partial x_j$ (exact) | SST-2003e | SST-2003 |

The suffix records the approximations [8]: without the
[energy coupling](#energy-coupling-k-energy-coupling) the
$-\tfrac23\rho k\,\delta_{ij}$ part of $\tau_{ij}$ is left out of the
momentum and energy equations ("m", "e"); "m" and "s" approximate the production
by $\mu_t S^2$, which TMR describes as exact for incompressible flow and
"typically considered a very good approximation, except perhaps for very high
Mach number flows".  (With `compressible`, the limiter acts on the deviatoric
part of $P$ only, as in OpenFOAM, where TMR's SST-2003 limits all of it; see
[below](#compressible-production-sst-production).)  `turbulence = SST-RC` adds
the [rotation/curvature correction](#rotationcurvature-correction-sst-rc).

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

### Compressible production (`sst-production`)

The production above is the SST-2003m form $\mu_t S^2$, which drops the
dilatation terms of $P = \tau_{ij}\,\partial u_i/\partial x_j$.  With
`sst-production = compressible` (default `incompressible`) MOSE keeps them, with
$\tau_{ij} = \mu_t\bigl(2S_{ij} - \tfrac23\,\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \tfrac23\,\rho k\,\delta_{ij}$:

$$
P_k = \min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\rho k\,\nabla\!\cdot\!\mathbf{v},
$$

$$
P_\omega = \frac{\gamma\,\rho}{\mu_t}\,\min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\gamma\,\rho\,\omega\,\nabla\!\cdot\!\mathbf{v}.
$$

This is the exact SST-2003 production (TMR's SST-2003e, since the
$-\tfrac23\rho k$ term still stays out of the momentum and energy equations), with
the dilatation terms arranged as in OpenFOAM's `kOmegaSST` [9]: the deviatoric
part, $\mu_t[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2] = 2\mu_t S^{d}_{ij}S^{d}_{ij} \ge 0$,
stays under the limiter, and the dilatation terms are added outside it (TMR's
SST-2003 limits the whole of $P$).  In attached, nearly incompressible boundary
layers the two forms coincide; through a compression ($\nabla\!\cdot\!\mathbf{v}<0$:
shock feet, reattachment) the dilatation terms amplify $k$ more than $\omega$
($\gamma<1$), so $\mu_t$ grows.  In the Mach-5
[shock/boundary-layer interaction](../vv/2D-swbli.md) the incompressible form gives
a separation bubble about 16% longer and a slower recovery of $C_f$.

### Rotation/curvature correction (SST-RC)

`turbulence = SST-RC` sensitizes the model to streamline curvature and system
rotation with the Spalart–Shur correction as adapted to SST by Smirnov and
Menter [10] (TMR SST-RC [8]; on SST-2003, SST-2003RC).  The production $P$ of
both equations is multiplied by

$$
f_{r1} = \max\!\bigl[\min(f_\text{rot},\,1.25),\,0\bigr], \qquad
f_\text{rot} = (1 + c_{r1})\,\frac{2\,r^\ast}{1 + r^\ast}
  \bigl[1 - c_{r3}\arctan(c_{r2}\,\tilde r)\bigr] - c_{r1},
$$

$$
r^\ast = \frac{S}{\Omega}, \qquad
\tilde r = \frac{2\,\Omega_{ik}\,S_{jk}}{\Omega\,D^3}\,\frac{DS_{ij}}{Dt}, \qquad
D^2 = \max(S^2,\;0.09\,\omega^2),
$$

with $c_{r1} = 1$, $c_{r2} = 2$, $c_{r3} = 1$, $S = \sqrt{2S_{ij}S_{ij}}$,
$\Omega = \sqrt{2\Omega_{ij}\Omega_{ij}}$.  $f_{r1}$ multiplies $P$ before the
production limiter, so it reaches the $\omega$ equation through the limited
$P$.  In simple shear $r^\ast = 1$ and $\tilde r = 0$, so $f_{r1} = 1$ and the
model reduces to SST.

Implementation notes: $r^\ast$ and the numerator of $\tilde r$ are shared with
SA-RC (`Lib_SpalartShur`).  $DS_{ij}/Dt$ is its steady convective part
$\mathbf{u}\!\cdot\!\nabla S_{ij}$; in a [rotating frame](rotating-frame.md)
the frame enters through the absolute vorticity in $\Omega_{ik}$ and $\Omega$,
but the frame-rotation term $(\epsilon_{imn}S_{jn} + \epsilon_{jmn}S_{in})\,\Omega^\text{rot}_m$
of the published $\tilde r$ is not included.  $\Omega$ is floored at
$10^{-30}$ s⁻¹, so cells with no vorticity take the limit $f_{r1} = 1.25$
instead of an undefined value (their production is negligible there).

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
\omega_\text{wall} = C\,\frac{\nu_w}{y_c^2}
$$

with $y_c$ the wall distance of the wall-adjacent cell centre and $\nu_w$ the
laminar viscosity at the wall.  `omega-wall-bc` selects $C$:

- `practical` (default): $C = 10 \times 6/\beta_1 = 800$, Menter's
  over-specification of the near-wall limit, robust on meshes that do not
  resolve the viscous sublayer;
- `asymptotic`: $C = 6/\beta_1 = 80$, the viscous-sublayer solution
  $\omega = 6\nu/(\beta_1 y^2)$ evaluated at $y_c$, for wall-resolved meshes.

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
terms](#numerical-treatment-of-the-source-terms).

### Energy coupling (`k-energy-coupling`)

The compressible RANS equations [11] carry the turbulent kinetic energy in two
places: the total energy, $E = e + \tfrac12 u_iu_i + k$, and the Boussinesq
Reynolds stress,
$\tau_{ij} = 2\mu_t\bigl(S_{ij} - \tfrac13\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \tfrac23\rho k\,\delta_{ij}$.
By default MOSE drops both (the TMR "m"/"e" variants): its energy is that of the
mean flow, $E_m = e + \tfrac12 u_iu_i$, and the $-\tfrac23\rho k$ term is not in
the stress.  With `k-energy-coupling = true` (two-equation models: SST and
Wilcox 2006) MOSE keeps both:

- the stress of the momentum and energy equations includes
  $-\tfrac23\rho k\,\delta_{ij}$, so it also enters the work term
  $\nabla\!\cdot\!(\boldsymbol{\tau}\cdot\mathbf{v})$ of the energy equation;
- the mean-flow energy gives up what the turbulence receives,

$$
\frac{\partial(\rho E_m)}{\partial t} + \dots = \dots - \bigl(P_k - \beta^\ast\rho\,\omega\,k\bigr),
$$

with the same source as the $k$ equation (including its
[point-implicit](#steady-runs-point-implicit-factor) factor), so $E_m + k$ is
conserved.  This is the total-energy form of TMR [11] and Wilcox [4] written for
$E_m$: subtracting the $k$ equation from the equation for $E = E_m + k$ leaves
exactly this sink, and the diffusion of $k$ in the two equations cancels.  The
pressure follows from $e = E_m - \tfrac12 u_iu_i$, MOSE's conserved energy being
$E_m$ in either case: for an ideal gas this is TMR's
$p = (\gamma-1)\,\rho\,(E - \tfrac12 u_iu_i - k)$.  TMR notes that the $k$
terms are sometimes dropped "for non-supersonic speed flows"; they matter where
$k$ is not small against the kinetic energy of the mean flow.
With the default, MOSE leaves $k$ out of the energy and the pressure of the
fluxes, as OpenFOAM-dev does: its compressible solvers carry the kinetic energy
of the mean flow, $\tfrac12|\mathbf{U}|^2$, and no $k$ (OpenFOAM-dev, commit
10f026fd, `applications/modules/shockFluid/thermophysicalPredictor.C` and
`applications/modules/fluid/thermophysicalPredictor.C`).

---

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

The Wilcox (2006) model [4, 12] with the exact production
$P = \tau_{ij}\,\partial u_i/\partial x_j$ and the optional TMR $k$-equation
limiter (below): TMR's Wilcox2006-klim, in its "e" form (Wilcox2006-klim-e) unless
[`k-energy-coupling`](#energy-coupling-k-energy-coupling) is on [8].

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

$k_\text{wall} = 0$. With `omega-wall-bc = asymptotic` the wall takes the near-wall limit of the $\omega$ equation of the
model, $\omega_\text{wall} = 6\nu/(\beta_0 y^2)$ with $\beta_0 = 0.0708$ (NASA Turbulence Modeling Resource, Wilcox 2006
page), that is $84.7\,\nu/y^2$, $y$ the distance of the centre of the wall cell from the wall. With `practical`, the default,
it takes Menter's condition of the [SST section](#wall-boundary-conditions), $800\,\nu/y^2$ with the $\beta_1$ of SST,
as the same page advises for this model.

### Cross-diffusion *(conditional)*

Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega \le 0$:

$$
\text{CD} = \frac{\sigma_d\,\rho}{\omega}\,
  (\nabla k\!\cdot\!\nabla\omega)
\qquad \text{if } \nabla k\!\cdot\!\nabla\omega \le 0
$$

### Production limiter

More permissive than SST, and in the $k$ equation only (the $\omega$ production
uses the unlimited $P$), as in TMR's Wilcox2006-klim:

$$
P_k = \min\!\bigl(\tau_{ij}\,\partial u_i/\partial x_j,\;20\,\beta^\ast\,\rho\,\omega\,k\bigr)
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
model. It adds a quadratic dependence on the rotation tensor to the modelled
stress, improving the prediction of anisotropy-driven secondary flows
(corner flows, square ducts) while reusing the eddy viscosity of the
underlying model. The code takes it as a suffix of the model name
(`-QCR2000`), but `turbulence` does not accept such a name today (its
allowed values are listed in the [registry](../user/registry.md)), so no run
uses it; the test `QCR2000Stress` checks the stress against the definition.

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
| `RANS_Extrapolate_Wall` | Ghost-cell value of the RANS variables at walls, $2Z_\text{wall} - Z_\text{cell}$ |

This design allows switching models at run time without recompilation.

### Wall ghost cells of the turbulence variables

Behind a viscous wall the ghost cell of each turbulence variable takes the value that puts the wall
value on the face,

$$
Z_\text{ghost} = 2\,Z_\text{wall} - Z_\text{cell}, \qquad Z = \tilde\nu,\ k,\ \omega ,
$$

with the wall values of the models above (`RANS_Extrapolate_Wall`). It is the rule of CFL3D,
«the actual wall boundary conditions for the turbulence quantities (unlike the primitive variables)
are applied at ghost cells» (S. L. Krist, R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual
(Version 5.0)*, NASA/TM-1998-208444, 1998, appendix H, eq. H-186). The wall flux routines set these
values; the source terms read them in the gradients of the wall cells: the $c_{b2}$ term of SA and
its term $-(\nu+\tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu/\sigma$, $F_1$ and the cross-diffusion of SST,
the cross-diffusion of Wilcox 2006. They do so at every stage
and in the stiffness pass of time-accurate runs (`Turbulence_Stiffness`), which calls the wall flux
routines for this alone (`Wall_RANS_Ghosts`). With $\tilde\nu = a\,y$ next to a wall the $c_{b2}$
term of the wall cell sees $\partial\tilde\nu/\partial y = a$ (part 5 of `MMSResiduals`); a ghost cell
that copied the cell, as MOSE did before, gave $a/2$. The other primitive variables are copied into
the wall ghost cells, with the velocity reversed (`Ghost_Wall_Extrapolation`); the fluxes through the
first interior faces use the second-order extrapolation of `Fill_Ghost_Cell`.

### Cell gradients of the source terms

The source terms take, in every cell, the central differences of the values along the three index
directions and map them to Cartesian components with a cell metric. With
`diffusive-metric = stencil` the metric is $M_s = (J_c^{T})^{-1}$ with
$J_c^{T}(:,d) = (\mathbf{x}_{c+e_d} - \mathbf{x}_{c-e_d})/2$, the same central differences of the
cell centres (the mean of the eight nodes): the gradient of a linear field is exact on stretched,
skewed and non-smooth cells (along an index direction with one cell, the edge of the cell and no
difference). On a mesh of equal parallelograms these gradients are the reference gradient of [14]
(Eq. 1), second order on smooth structured meshes. In a wall cell, for a field that is not linear,
the reflected ghost value $2\phi_\text{wall} - \phi_\text{cell}$ makes the gradient first order, as
the Green–Gauss gradient is in the boundary cells of [14] (Eq. 13).
The ghost cells sit where their value is right when the sources are computed: behind a
wall the velocity and the turbulence variables, whose ghost values are $2\phi_\text{wall} -
\phi_\text{cell}$, at the reflection of the cell centre in the face centre; the densities, the pressure
and the other copied variables at the cell itself, which makes the difference across the wall cell
one-sided and exact for any linear field; behind the other faces the ghost cells of the fluxes, at
the positions of the face metric ([Numerics](numerics.md)): behind a zero-gradient face (0,
401-407, 420), whose ghost cell copies the cell, the cell itself, so the difference across the
boundary cell is one-sided there too. SA-RC copies the velocity gradient of the
cell into the ghost cells of every boundary face, so its strain gradient takes the cell itself there.
The metrics are computed once at setup (`Cell_Source_Metrics`) and chosen in one place
(`Source_Metrics`) by the source terms of SA, SA-RC, SST, Wilcox 2006 and SSG/LRR and by the eddy
viscosity of the output. With `inverse-mean` and `mean` the metric is that of the cell, as before:
the central difference across a cell next to a jump of size $r$ is then wrong by $(r-1)/4$, on a
skewed wall cell a copied density takes a normal gradient it does not have (part 5 of `MMSResiduals`),
and next to a zero-gradient face the normal component is half the one-sided difference.

### Blowing / wall-suction correction

An optional correction (`blowing_corr`) modifies the wall boundary-layer
treatment for SST and Wilcox 2006 models in the presence of wall
injection or suction.

### Wall distance

The distance $y$ of the models (SA, the SST blending functions, the wall
value of $\omega$) is the minimum distance from the cell centre, the mean of
its eight nodes, to the wall surface, as the NASA Turbulence Modeling
Resource requires (ref. 13). It is not the distance to the nearest wall node
or face centre, which differs from it wherever the grid lines are not normal
to the wall; next to a convex corner or at the end of a wall it is the
distance to the corner or the end.

The distance to the centre of the nearest wall face is wrong by a fixed
fraction in the first cell of a skewed mesh. With the grid lines at an angle
$\theta$ to the wall normal, the centre of the first cell, at the height
$d_1$, lies $d_1\tan\theta$ along the wall from the centre of its wall face,
the nearest face centre while $d_1\tan\theta$ is at most half the face, at
the distance $d_1/\cos\theta$. The relative error,
$1/\cos\theta - 1$, does not fall with refinement: 1.5 % at 10°, 6.4 % at 20°,
15.5 % at 30°, 41 % at 45°. It enters squared the destruction of SA,
$c_{w1} f_w (\tilde\nu/d)^2$ (a factor $\cos^2\theta$: 0.97, 0.88, 0.75,
0.50), and the arguments of $F_1$ and $F_2$ and the wall value
$\omega_w = 60\nu/(\beta_1 d_1^2)$ of SST. That distance is exact only on
meshes orthogonal to the wall, with the cell centres above the face centres.

- Each face of a wall boundary (types 301, 302, 503–505) is taken as the four
  triangles between its centre and its sides, which are the face itself when
  the face is plane and convex, as on every 2-D and axisymmetric mesh (on a
  plane face that is not convex, a pathological cell, the centre can fall
  outside the face, and the triangles are not the face). The distance to a
  triangle is the distance to its plane when the projection of the point falls
  inside it, otherwise the distance to the nearest side.
- On a 2-D axisymmetric wedge (BC 200) the distance is measured in the
  meridional plane $(x, r)$, $r = \sqrt{y^2+z^2}$, from the cell centre to the
  wall curve. This is the distance to the wall surface of revolution, and it
  does not depend on the angle of the wedge. A wall on the faces of the wedge
  (faces 5–6) stops the run.
- A k-d tree over the face centres, with the radius of each face, finds the
  nearest face: the exact distance is computed only for the faces whose centre
  lies within the best distance found so far plus their radius.
- Across periodic boundaries the wall goes on: the distance is the minimum over
  the periodic images of the wall too. A periodic connection (type 201) gives
  the translation from its source face to its own face, and on a 3-D mesh
  without sectors so does a block closed in $k$ by type `200` (its last node
  plane is the first one translated). The distinct translations are reduced
  to a basis $B_1 \ldots B_n$ of the lattice they generate, shortest first,
  and then pairwise (Lagrange: $B_j$ becomes $B_j - m B_i$, $m$ the integer
  nearest to $B_i\cdot B_j / B_i\cdot B_i$, while that shortens it) and, with
  three translations, by threes ($B_l$ becomes $B_l - m_i B_i - m_j B_j$, with
  $m_i, m_j$ the integers nearest to the least-squares coefficients of $B_l$ on
  $B_i, B_j$ or next to them, while that shortens it). Three nearly dependent
  translations can be reduced pair by pair (no pairwise step shortens them)
  and still have a short sum: $B_1 = (1, 0, 0)$,
  $B_2 = (-\tfrac12, \tfrac{\sqrt3}{2}, 0)$,
  $B_3 = (-\tfrac12, -\tfrac{\sqrt3}{2}, \varepsilon)$ give
  $B_1 + B_2 + B_3 = (0, 0, \varepsilon)$, which only the step by threes finds.
  Each step adds to a vector an integer combination of the others, so the
  lattice is the same, and the bound below holds for any basis of it: the
  reduction only keeps the number of images small. The
  images are the integer combinations $\sum_i k_i B_i$ of the basis (along $x$
  on an axisymmetric wedge). The two sides of a periodic boundary match only
  to the round-off of the coordinates in the mesh file: translations that
  differ by at most $10^{-2}$ of the thinnest cell across the periodic
  boundaries are the same, as in the type of the mesh. The thickness of a
  cell across a periodic face is the distance from the centre of the face to
  the centre of the opposite face of the cell, along the normal of the face;
  for a block closed in $k$ by type `200`, the move of its nodes per cell in
  $k$. The round-off of a coordinate grows with its size, so a tolerance
  relative to the size of the coordinates would depend on where the mesh
  lies; relative to the cells it does not. On a 3-D mesh whose
  sectors are closed by the rotational periodicity (type `200` on faces 5–6)
  the images are also the turns by multiples of the sector angle about the
  $x$ axis, up to half a turn either way.
- Only a bounded set of images can hold the nearest wall, whatever the size of
  the mesh against its periods. Let $B^*_i$ be the dual basis
  ($B^*_i\cdot B_j = \delta_{ij}$), $G_{ij} = B_i\cdot B_j$ the Gram matrix,
  and $w_i = B^*_i\cdot w$ the coordinates of a vector $w$. Every wall face
  is moved, whole, by the lattice vector that brings its centre $m$ to
  $|(m - a)_i| \le 1/2$, $a$ the centre of the first wall face, and so is the
  centre $c$ of every cell before its search: the periodic wall is the same
  set. Let $p$ be the point of a face nearest to $c$ over all the images, and
  $s$ the image that holds it: $s$ is a lattice vector closest to $v = c - p$
  (another one would bring a point of the same face nearer), so $u = v - s$ is
  not farther from 0 than from $\pm B_j$, $|u\cdot B_j| \le G_{jj}/2$, and
  since $u_i = \sum_j (G^{-1})_{ij}\, u\cdot B_j$,
  $|u_i| \le \beta_i = \tfrac12 \sum_j |(G^{-1})_{ij}|\, G_{jj}$ ($1/2$ on
  an orthogonal basis). A face is the four triangles between its centre $m$
  and its sides, so $|(p - m)_i| \le \rho_i$, $\rho_i$ the largest
  $|(q - m)_i|$ of a corner $q$ of a face, and
  $|v_i| \le |(c - a)_i| + |(m - a)_i| + |(p - m)_i| \le 1 + \rho_i$. The
  coefficients of $s$ are $k_i = v_i - u_i$: the images searched are
  $|k_i| \le \lfloor 1 + \rho_i + \beta_i \rfloor$, one on each side
  ($|k_i| \le 1$) when $\rho_i < 1/2$ on an orthogonal lattice. On a sector
  every translation lies along $x$, and the turns keep the coordinate along
  $x$, so the same bound holds for each turn. An image is searched only when
  the box around the moved wall, seen from the cell moved back by the image,
  is not farther than the distance found so far, so a mesh without periodic
  boundaries costs nothing more. Translations that do not form a lattice stop
  the run: one that is not an integer combination of the others, more than
  three independent ones, or one not along $x$ on a wedge or with sectors; so
  does a number of images too large for an integer.
- Each cell takes the sand-grain roughness $k_s$ of the face that gives the
  minimum distance; between faces at the same distance, the first one in the
  boundary file.

!!! warning "Twisted wall faces in 3-D"
    A wall face whose corners $\mathbf{x}_1 \ldots \mathbf{x}_4$ (in cyclic
    order) are not coplanar is the bilinear surface through them, with twist
    $\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$.
    The four triangles lie within $|\mathbf{T}|/16$ of it, so the distance is
    exact only to within $|\mathbf{T}|/16$: a relative error up to
    $(|\mathbf{T}|/16)/y$, small next to cells thicker than $|\mathbf{T}|/16$
    and of order one or more in thinner ones. On a unit face twisted by
    $t = |\mathbf{T}|$, at height $h$ above it, the largest error is 6 % for
    $t = h = 0.1$, 62 % for $t = 10^{-2}$ and $h = 10^{-3}$, and 62 times the
    distance for $t = 0.1$ and $h = 10^{-4}$. On a smooth wall $|\mathbf{T}|$
    falls as the square of the face size, so the error falls with refinement.
    2-D and axisymmetric meshes are not affected: their wall faces are plane,
    and the four triangles are the face. At start-up MOSE reports the largest
    gap bound of the wall faces over the distance of the cell next to the face,
    $\min(|\mathbf{T}|/16, |\mathbf{T}_n|/2)/y$ with $\mathbf{T}_n$ the part of
    the twist along the normal of the face (zero on a plane face; the second
    bound holds on faces convex in that plane), as a warning above 1 %.

Test `MeshGeometry`/`WallDistance`.

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
⇒ more near-wall dissipation and a more robust start-up).  The TMR SST
freestream used in the [SWBLI case](../vv/2D-swbli.md) is
$k_\infty = 9\times10^{-9}\,a_\infty^2$, $\omega_\infty = 10^{-6}\,a_\infty^2/\nu_\infty$
($\mu_t/\mu = 0.009$): at $p=4000$ Pa, $T=68.3$ K
($\rho_\infty=0.2041$ kg m⁻³, $a_\infty=165.7$ m s⁻¹,
$\nu_\infty=2.257\times10^{-5}$ m² s⁻¹ with Sutherland's law) this is
$k_\infty = 2.47\times10^{-4}$ m² s⁻², $\omega_\infty = 1216$ s⁻¹.
Values of $\omega_\infty$ orders of magnitude below the $\mu_t/\mu$ estimate are
a frequent cause of $k$ runaway.

**Inputs and files.**  The turbulence keys of ATLAS come in pairs, a specific one
and a density-weighted one, under the names of the
[NASA Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/): `nuTilde` $= \tilde\nu$ and `rhoNuTilde` $= \rho\tilde\nu$, `k` and
`rhoK` $= \rho k$, `omega` and `rhoOmega` $= \rho\omega$.  $\tilde\nu$ is the
variable the Spalart-Allmaras model transports, not the eddy viscosity
$\mu_t = \rho\tilde\nu f_{v1}$, which is not an input and is written as `mut`
(`sol-variables = transport`).  The specific values are those
TMR and OpenFOAM give, and the same ones serve the initial condition
(`[ICB-Block*]`) and the inflow (`[inflow]`).  MOSE solves for the
density-weighted variables $\rho\tilde\nu$, $\rho k$, $\rho\omega$:

- the initial condition: ATLAS ICB multiplies the specific keys `nuTilde`, `k`, `omega`
  by the density of each cell (it also takes the density-weighted keys `rhoNuTilde`,
  `rhoK`, `rhoOmega`, written as given);
- the inlet: the boundary takes the specific `nuTilde`, `k`, `omega` only, and MOSE
  multiplies the `bc.txt` values by the density of the inlet state;
- the IC, restart and solution files all carry the density-weighted variables
  under the names of the density-weighted keys: `rhoNuTilde`, `rhoK`, `rhoOmega` (divide
  by $\rho$ to read $\tilde\nu$, $k$ and $\omega$).  MOSE reads them by position
  and, in an ASCII file, checks the names
  ([Initial Conditions](../user/initial-conditions.md#required-variables)).

For the example above: `k = 2.4699e-4`, `omega = 1215.8`.

---

## Numerical treatment of the source terms

The turbulence source terms are added to the residual and advanced with the same explicit
Runge–Kutta stages as the mean flow. Their destruction parts are stiff near walls: for the
$k$–$\omega$ models $\omega$ grows like $1/y^2$ towards the wall, and an explicit stage on
$\dot q = -D(q)$ gives $q - \Delta t\,D = q\,(1 - \Delta t\,D/q)$, negative for $\Delta t\,D/q > 1$. Here
$q$ is the conserved turbulence variable; the destruction rate $D/q$ is at most $d$, the derivative
of the destruction ($d/2$ for $\tilde\nu$ and $\omega$, whose destruction is quadratic):

| Model | Equation | $d = \partial D/\partial q$ |
|-------|----------|------------------------------|
| SST / Wilcox | $k$ | $\beta^\ast\,\omega$ |
| SST / Wilcox | $\omega$ | $2\,\beta\,\omega$ |
| SA | $\tilde\nu$ | $2\,c_{w1}\,f_w\,\tilde\nu/y^2$ |

($y$ the wall distance; $f_w$ is held fixed.) Steady and time-accurate runs treat this
stiffness in two different ways.

### Steady runs: point-implicit factor

With `point-implicit = .true.` (the default) the increment of each turbulence equation is divided
by $1 + \Delta t\,d$, with the local time step $\Delta t$ of the cell:

$$
q^{n+1} = q^n + \frac{\Delta t}{1 + \Delta t\,d}\,\bigl( F(q^n) + S(q^n) \bigr),
$$

where $F$ stands for the fluxes and $S = P - D$ for the source (in each Runge–Kutta stage, with the
residual of the stage). This is the point-implicit form of Bussing and Murman,
$D\,\Delta q = \Delta t\,(F + S)$ with $D = 1 - \Delta t\,\partial S/\partial q$ approximated by its
destruction part $1 + \Delta t\,d$, in its first-order form, «suitable for steady-state
calculations» (H. C. Yee, J. L. Shinn, NASA TM-89415, 1986, eq. 3.1). Like any approximation of the
implicit operator of a delta form, it leaves the steady state where the residual vanishes
(T. H. Pulliam, *Time Accuracy and the Use of Implicit Methods*, AIAA Paper 93-3360, eq. 8): the
converged solution is that of $F + S = 0$, whatever the time step, and the factor only lets
steady runs keep the mean-flow CFL near walls. The residuals of $\tilde\nu$, $k$ and $\omega$ printed during
the run (the change of the variable over the step) are scaled by the same factor, so near walls they
are smaller than without it; the convergence check `res-threshold` uses the density residual and is
not affected.

Up to this version the factor divided the source alone,
$q^{n+1} = q^n + \Delta t\,\bigl(F + S/(1 + \Delta t\,d)\bigr)$, and the converged solution, wall shear
and heat flux included, depended on the local time step. Steady solutions computed with the default
change accordingly.

```ini
[MOSE-Turbulence]
point-implicit = .true.   ; (default) point-implicit factor on the turbulence increment, steady runs
```

### Time-accurate runs: explicit sources

In time-accurate runs the factor is not applied, whatever `point-implicit` says, and a warning at
start-up recalls it. Since $(F + S)/(1 + \Delta t\,d) = (F + S)(1 - \Delta t\,d) + O(\Delta t^2)$, the factor
would make the turbulence equations, and the mean flow through $\mu_t$, first-order accurate in
time. With the $k$–$\omega$ models the error of the factor would not even vanish when the mesh and the
step are refined together, as in an explicit run, whose global step goes as $h^2$ in the wall cells: in
the first cell $\omega \approx 6\nu/(\beta_1 y^2)$ gives $\Delta t\,d \approx 24.7$ `vnn` at every level (7.4 at
`vnn = 0.3`, a factor of 0.12), and the next two cells stay of order one. With Spalart–Allmaras
$\Delta t\,d \approx 1.64\,f_w\,y_1^+$ at `vnn = 0.3` ($y_1^+$ the height of the first cell in wall units): it
vanishes as the mesh is refined, but is of order one at the usual heights.
The sources are explicit, like the fluxes, and the local time step of each cell
([Time Integration](time-integration.md#local-time-step)) includes their stiffness:

$$
\frac{1}{\Delta t_i} = \bigl(\text{convective and diffusive rates}\bigr) + d_i ,
\qquad d_i = \max_q \frac{\partial D_q}{\partial q} ,
$$

with the rates $d$ of the table, evaluated on the state at the start of the step (one evaluation
of the sources before the first stage, with the wall ghost cells of the turbulence variables set
from the wall values as in the stages). Then $\Delta t\,d \le 1$: a forward-Euler stage cannot make
$\tilde\nu$, $k$ or $\omega$ negative through their destruction, and neither can Heun and RK3, whose
stages are convex combinations of forward-Euler stages (W. Hundsdorfer, *Numerical Solution of
Advection–Diffusion–Reaction Equations*, lecture notes, 2000, §6.2, Theorem 6.4; C.-W. Shu,
S. Osher, *J. Comput. Phys.* 77, 1988, Proposition 2.1 of its report version, ICASE Report 87-33).
The destruction of $\tilde\nu$ and of $\omega$ is quadratic, $D = (d/2)\,q$, which leaves them a
margin of 2. Where the destruction is faster than the transport, as in the wall cells of fine
meshes, above all with SST, the step is smaller than the convective and diffusive limits.

The bound covers the destruction only. Three sinks stay outside it: the cross-diffusion of SST,
$2(1-F_1)\,\rho\,\sigma_{\omega 2}\,\nabla k\cdot\nabla\omega/\omega$, negative where $\nabla
k\cdot\nabla\omega < 0$ (Wilcox 2006 switches its cross-diffusion off there); with `sst-production =
compressible`, the dilatation part of the production of SST, $-\tfrac23\rho k\,\nabla\cdot\mathbf{u}$
and $-\tfrac23\gamma\rho\omega\,\nabla\cdot\mathbf{u}$, negative where the flow expands; and the
production of SA-R and SA-RC, which can be negative («the SA-RC production term can be negative»,
NASA Turbulence Modeling Resource). Like the transport, they can make a stage negative, and the count
at the end of the run shows it.

The convex combination keeps the bound of the forward-Euler stage only where that bound holds on the
state each stage starts from: Shu and Osher assume that the forward-Euler step keeps the property,
under one step restriction, for every state it is applied to (ICASE Report 87-33, (2.7) and
Proposition 2.1), while here the restriction depends on the state through $d$, and the rates are
those of the state at the start of the step. A stage from a state where the rate is $d_s$ takes
$q_s - \Delta t\,(d_s/2)\,q_s$ through a quadratic destruction, positive while $\Delta t\,d_s < 2$; with
$\Delta t\,d \le 1$ this holds as long as $d_s < 2\,d$. The margin of 2 of $\tilde\nu$ and $\omega$ therefore
covers a rate that grows within the step up to twice its value at the start, and no more; $k$, whose
destruction is linear, keeps only the margin by which its rate $\beta^\ast\omega$ lies below $d$. An
impulsive start of a $k$–$\omega$ model from a state far from the wall value of $\omega$ goes beyond
that: during the first step the wall condition raises $\omega$ in the wall cells by orders of
magnitude, and with the full step the explicit destruction overshoots ($\omega < 0$, then the floor
and a very large $\mu_t$). Start such runs with the CFL ramp (`cfl-rise-threshold`).

At the end of a RANS run MOSE prints how many cell updates had a stage that made $\tilde\nu$, $k$ or
$\omega$ negative before the realizability floor:

```text
 Turbulence variables below zero after a stage (set to the floor): 0
```

In a steady run the count covers the whole run: the start, and the fronts of $k$ that the fluxes
push below zero where $k$ is about zero, which can last to the last iteration; it is a total, which
convergence does not bring back to zero.

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
   engineering applications," *AIAA J.*, 32(8), pp. 1598–1605, 1994,
   doi:10.2514/3.12149.
3. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial
   experience with the SST turbulence model," in *Turbulence, Heat
   and Mass Transfer 4*, Begell House, pp. 625–632, 2003.
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
8. NASA Langley Turbulence Modeling Resource, "The Menter Shear Stress
   Transport turbulence model," <https://tmbwg.github.io/turbmodels/sst.html>
   (formerly turbmodels.larc.nasa.gov/sst.html): definitions and naming of the
   SST variants (SST-2003, SST-2003m/s/e, SST-RC) and of Wilcox2006-klim.
9. OpenFOAM Foundation, OpenFOAM-10, `kOmegaSST` model
   (`src/MomentumTransportModels/momentumTransportModels/Base/kOmegaSST/kOmegaSSTBase.C`).
10. P. E. Smirnov, F. R. Menter, "Sensitization of the SST turbulence model to
    rotation and curvature by applying the Spalart–Shur correction term,"
    *J. Turbomach.*, 131(4), 041010, 2009, doi:10.1115/1.3070573.
11. NASA Langley Turbulence Modeling Resource, "Implementing turbulence models
    into the compressible RANS equations,"
    <https://tmbwg.github.io/turbmodels/implementrans.html>.
12. D. C. Wilcox, "Formulation of the k-ω turbulence model revisited,"
    *AIAA J.*, 46(11), pp. 2823–2838, 2008, doi:10.2514/1.36541.
13. NASA Langley Research Center, Turbulence Modeling Resource, "The
    Spalart–Allmaras Turbulence Model",
    https://turbmodels.larc.nasa.gov/spalart.html (TMBWG/turbmodels commit
    20a39f5, consulted 2026-09-29): minimum distance to the wall.
14. A. Syrakos, S. Varchanis, Y. Dimakopoulos, A. Goulas, J. Tsamopoulos, "A
    critical analysis of some popular methods for the discretisation of the
    gradient operator in finite volume methods," *Phys. Fluids*, 29, 127103,
    2017, doi:10.1063/1.4997682.
