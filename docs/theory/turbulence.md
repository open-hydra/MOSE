# Turbulence Modelling

MOSE provides Reynolds-Averaged Navier–Stokes (RANS) closures for
turbulent flows, selected at run time via the `turbulence` key in the
`[MOSE-Physics]` section of the input file.  The eddy-viscosity models
described in detail below range from the one-equation Spalart–Allmaras
model (with several variants) to the two-equation Menter SST and
Wilcox 2006 $k$–$\omega$ models.  A full Reynolds-stress closure
(SSG–LRR) and an algebraic non-linear correction (QCR2000) are in the
code, but `turbulence` does not accept them — see [Other models](#other-models).

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

on the viscous walls and on the grain of a solid rocket motor (types `301`, `302`, `502`,
`503`–`505`); the wall distance $y$ is measured from them and from the walls coupled to another
solver (`103`, `104`). A slip wall (type `300`, which is also the symmetry plane) follows the
standard model: the ghost cell copies the cell, a zero normal derivative of $\tilde\nu$, and the
slip wall is not a wall for the distance («the free-slip wall points are not included in the search
when $d$ is computed», P. R. Spalart, S. R. Allmaras,
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
    number and the right side a speed, so the factor depends on the unit of
    length. For $\Omega = 10^3$ 1/s, $\nu_t = 10^{-3}$ m²/s and $a = 340$ m/s
    it is 0.99 with lengths in metres and 0.22 with lengths in millimetres.
    The correction (the code cites doi:10.2514/2.3967) is not available: use
    `SA` or another model.

#### SA-R — Rotation Correction

`turbulence = SA-R` replaces $\tilde S$ in the production by $\tilde S + C_\text{rot}\min(0,\,S - \Omega)$,
with the strain-rate magnitude $S = \sqrt{2\,S_{ij}\,S_{ij}}$:

$$
P_\text{SA-R} = c_{b1}\,\bigl[\tilde{S} + C_\text{rot}\,\min(0,\,S - \Omega)\bigr]\,\rho\tilde\nu
\qquad (C_\text{rot} = 2)
$$

The production is reduced where the vorticity exceeds the strain rate, as in vortex cores, and
unchanged where it does not (NASA Turbulence Modeling Resource, Spalart–Allmaras page, SA-R, after
Dacles-Mariani et al., *AIAA J.* 33(9), 1995).

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

MOSE takes $\partial S_{ij}/\partial x_k$ with central differences of the velocity gradients
$G = \nabla\mathbf{u}$ of the cell and of its neighbours (`Compute_Velocity_Gradient`,
`Compute_RC_Terms`), mapped to Cartesian components with a cell metric. Next to a block face the
neighbour is the first ghost cell, and the difference is right when the ghost cell holds the value
of a definite cell and the metric takes the position of that same cell: for $S_{ij}$ linear in
space the difference is then exact (with `diffusive-metric = stencil`, the metric of the
differences of these positions: [Cell gradients of the source terms](#cell-gradients-of-the-source-terms)).

- Behind a face between blocks (`101`) or of a periodic pair (`201`) the ghost cell holds $G$ of the
  cell across, the source cell of the record (from another MPI rank by message), and the metric
  takes its centre, plus the translation of a periodic pair. $G$ is a tensor in Cartesian
  components: it is copied whatever the orientation of the two blocks. The cell next to the face
  then takes the central difference of the undivided mesh, and the result does not depend on how
  the mesh is divided into blocks or among ranks.
- Behind a face of rotational periodicity (`200`), the ghost cell holds $R\,G\,R^{\mathsf T}$, $G$ of
  the image cell that the velocity takes there, turned by the same rotation $R$ about $x$: a flow
  invariant under the turn, $\mathbf{u}(R\mathbf{x}) = R\,\mathbf{u}(\mathbf{x})$, has
  $G(R\mathbf{x}) = R\,G(\mathbf{x})\,R^{\mathsf T}$. On a 3-D sector the image is the cell next to the
  opposite face and the metric takes its position turned; on the 2-D wedge it is the cell itself,
  and along $k$ the metric takes the edge of the cell, as for the velocity. With a zero sector angle
  (a planar mesh, a translated sector) the ghost cell holds a copy of the cell, as below.
- Behind the other faces (walls, symmetry, inflow, outflow, chimera, mapped faces) the ghost cell
  holds a copy of the gradient of the cell and the metric takes the cell itself: the difference
  across the boundary cell is one-sided, a choice of the model.

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
laminar viscosity at the wall.  `omega-wall-bc` selects $C$:

- `practical` (default): $C = 10 \times 6/\beta_1 = 800$, Menter's
  over-specification of the near-wall limit, robust on meshes that do not
  resolve the viscous sublayer;
- `asymptotic`: $C = 6/\beta_1 = 80$, the viscous-sublayer solution
  $\omega = 6\nu/(\beta_1 y^2)$ evaluated at $y_c$, for wall-resolved meshes.

The wall flux of the walls that blow (`504` and `505` while they ablate, a `502` grain while it
burns) takes $C = 800$ with either value.

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

with the same source as the $k$ equation, so $E_m + k$ is conserved.  In steady
runs the [point-implicit](#steady-runs-point-implicit-factor) factor treats the
energy and $k$ equations as one system: the two increments add up to that of the
fluxes of $E_m + k$ at every stage, and the steady state does not depend on the
time step.  This is the total-energy form of TMR [11] and Wilcox [4] written for
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

$k_\text{wall} = 0$. With `omega-wall-bc = asymptotic` the wall takes the near-wall limit of the $\omega$ equation of the
model, $\omega_\text{wall} = 6\nu/(\beta_0 y^2)$ with $\beta_0 = 0.0708$ (NASA Turbulence Modeling Resource, Wilcox 2006
page), that is $84.7\,\nu/y^2$, $y$ the distance of the centre of the wall cell from the wall. With `practical`, the default,
it takes Menter's condition of the [SST section](#wall-boundary-conditions), $800\,\nu/y^2$ with the $\beta_1$ of SST,
as the same page advises for this model. The wall flux of the walls that blow (`504` and `505` while they ablate, a `502`
grain while it burns) takes $800\,\nu/y^2$ with either value.

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

## Other models

### SSG–LRR (Reynolds-stress model)

In addition to the eddy-viscosity closures above, the code has the
**SSG–LRR** differential Reynolds-stress model (`Lib_SSGLRR`),
which transports the six independent components of the Reynolds-stress
tensor together with a length-scale variable instead of assuming the
Boussinesq relation. It blends the Speziale–Sarkar–Gatski (SSG)
pressure–strain model away from walls with the Launder–Reece–Rodi (LRR)
model near walls. Being anisotropy-resolving, it captures secondary
flows and strong streamline-curvature effects that the two-equation
models miss, at the cost of the additional transport equations.
`turbulence` does not accept the name `SSGLRR` (its allowed values are listed
in the [registry](../user/registry.md)), so no run uses it.

### QCR2000 (non-linear constitutive correction)

The **Quadratic Constitutive Relation** (Spalart, 2000) is an algebraic,
non-linear correction to the Boussinesq stress rather than a standalone
model. It adds a quadratic dependence on the rotation tensor to the modelled
stress, improving the prediction of anisotropy-driven secondary flows
(corner flows, square ducts) while reusing the eddy viscosity of the
underlying model. The code takes it as a suffix of the model name
(`-QCR2000`), but `turbulence` does not accept such a name (its
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

with the wall values of the models above, here of the specific variables $Z = \tilde\nu, k, \omega$
(the ghost value of $\rho Z$ is that times the density of the ghost cell). It is the rule of CFL3D,
«the actual wall boundary conditions for the turbulence quantities (unlike the primitive variables)
are applied at ghost cells» (S. L. Krist, R. T. Biedron, C. L. Rumsey, *CFL3D User's Manual
(Version 5.0)*, NASA/TM-1998-208444, 1998, appendix H, eq. H-186). `Fill_Ghost_Cell` sets these
values once per stage, from the wall values of the wall flux routines (the same wall state), with
the other ghost values of the wall ([Numerics](numerics.md), *Ghost Cells of Walls and 400 Faces*):
the fluxes through the first interior faces and the source terms read the same ghost cells, at every
stage. The source terms read
them in the gradients of the wall cells: the $c_{b2}$ term of SA and its term
$-(\nu+\tilde\nu)\,\nabla\rho\cdot\nabla\tilde\nu/\sigma$, $F_1$ and the cross-diffusion of SST, the
cross-diffusion of Wilcox 2006. With $\tilde\nu = a\,y$ next to a wall the $c_{b2}$ term of the wall
cell sees $\partial\tilde\nu/\partial y = a$ (part 5 of `MMSResiduals`); a ghost cell that copied the cell
would give $a/2$. An SRM grain (`502`) takes the
same rule for its turbulence variables, with the wall values of its flux: those of the blowing wall while
it burns, of the wall at the flame temperature when it does not.

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
The ghost cells sit where their value is right: those of the fluxes, at the positions
of the face metric ([Numerics](numerics.md)); behind a wall the position of the class of the
quantity, one metric per group: the velocity, the density (whose ghost value follows the temperature:
exact where the pressure and the composition do not vary across the wall cell), the turbulence
variables (a metric of their own only next to `503`–`505`, where the velocity is extrapolated and they
take their wall value). Behind a zero-gradient face (0, 401-407, 420), whose ghost cell copies the
cell, the position is the cell itself, so the difference across the boundary cell is one-sided. The
strain gradient of SA-RC takes, behind the faces between blocks and of a turned sector, the position
of the cell whose velocity gradient the ghost cell holds, and the cell itself behind the other faces
(see SA-RC above).
The metrics are computed once at setup (`Cell_Source_Metrics`) and chosen in one place
(`Source_Metrics`) by the source terms of SA, SA-RC, SST, Wilcox 2006 and SSG/LRR and by the eddy
viscosity of the output. With `mean` the metric is that of the cell:
the central difference across a cell next to a jump of size $r$ is then wrong by $(r-1)/4$, on a
skewed wall cell a copied density takes a normal gradient it does not have (part 5 of `MMSResiduals`),
and next to a zero-gradient face the normal component is half the one-sided difference.

### Blowing / wall-suction correction

The code has a wall value of $\omega$ for walls that inject or suck gas
(`SST_Blowing_Correction`, used by SST and Wilcox 2006), which the suffix
`-blowcorr` of the model name selects; `turbulence` does not accept such a
name, so no run uses it.

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

- Each face of a wall boundary (types 103, 104, 301, 302, 502, 503–505) is the bilinear
  surface through its four nodes, the face of the cells
  ([Cell volume](numerics.md#cell-volume)):
  $\mathbf{x}(s,t) = \mathbf{x}_1 + s\mathbf{a} + t\mathbf{b} + st\mathbf{T}$,
  $\mathbf{a} = \mathbf{x}_2 - \mathbf{x}_1$, $\mathbf{b} = \mathbf{x}_4 - \mathbf{x}_1$,
  with the twist $\mathbf{T} = \mathbf{x}_1 - \mathbf{x}_2 + \mathbf{x}_3 - \mathbf{x}_4$
  (corners in cyclic order). The distance to it is its minimum over the face,
  exact to the rounding of the size of the face and of the distance: a face plane
  to the resolution of its coordinates (every face of a 2-D mesh) takes the four
  triangles between its centre and its sides, any other face a certified minimum
  ([Code Structure](../development/structure.md#wall-distance) gives the
  algorithm). The four triangles of a twisted face would be wrong by up to
  $(|\mathbf{T}|/16)/y$, 6 % on a unit face twisted by 0.1 at the height 0.1.
- On a 2-D axisymmetric wedge (BC 200) the distance is measured in the
  meridional plane $(x, r)$, $r = \sqrt{y^2+z^2}$, from the cell centre to the
  wall curve. This is the distance to the wall surface of revolution, and it
  does not depend on the angle of the wedge. A wall on the faces of the wedge
  (faces 5–6) stops the run. It is the one exception to the single geometric
  model of the cells: for the wall distance the wall of the wedge is the
  surface of revolution of its meridional curve, for the volumes and the fluxes
  the cells are polyhedra with plane faces; the two differ by
  $r(1 - \cos(\delta/2)) \approx r\delta^2/8$.
- The distance and the face that gives it do not depend on the threads or on
  the ranks.
- Across periodic boundaries the wall goes on: the distance is the minimum over
  the periodic images of the wall too. A periodic connection (type 201) gives
  the translation from its source face to its own face, and on a 3-D mesh
  without sectors so does a block closed in $k$ by type `200` (its last node
  plane is the first one translated). The images are the integer combinations of
  these translations (along $x$ on an axisymmetric wedge). The two sides of a
  periodic boundary match only to the round-off of the coordinates in the mesh
  file: translations that
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
  $x$ axis, up to half a turn either way. A channel periodic in $x$ and turned by
  30 degrees about $z$ gives translations that differ by 1.1e-4 of its cells with
  five significant digits, by 8.8e-4 with four, and by 0.15 with two, which stops
  the run. Translations that do not form a lattice stop the run: one that is not
  an integer combination of the others, more than three independent ones, or one
  not along $x$ on a wedge or with sectors. Only a bounded set of images can hold
  the nearest wall, and only those are searched: a mesh without periodic
  boundaries costs nothing more.
- Each cell takes the sand-grain roughness $k_s$ of the face that gives the
  minimum distance; between faces at the same distance, the first one in the
  boundary file.

Tests `MeshGeometry`/`WallDistance` and `BilinearDistance`.

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

With [`k-energy-coupling = true`](#energy-coupling-k-energy-coupling) the energy of the mean flow has
the sink $-S_k$ of the $k$ equation, whose derivative with respect to $\rho k$ is $+\beta^\ast\omega$, so
the factor is that of the two equations together: the implicit operator of the $k$ and energy rows
is $D = \begin{pmatrix} 1 + \Delta t\,\beta^\ast\omega & 0 \\ -\Delta t\,\beta^\ast\omega & 1 \end{pmatrix}$, and the
energy takes its own increment plus $\Delta t\,\beta^\ast\omega$ times that of $k$,

$$
\Delta(\rho E_m) = \Delta t\,\bigl(F_E - S_k\bigr) + \Delta t\,\beta^\ast\omega\,\Delta(\rho k), \qquad
\Delta(\rho k) = \frac{\Delta t}{1 + \Delta t\,\beta^\ast\omega}\,\bigl(F_k + S_k\bigr).
$$

The two increments add up to $\Delta t\,(F_E + F_k)$, the fluxes of $E_m + k$ alone, at every stage, and
the converged solution, $F_k + S_k = 0$ and $F_E = S_k$, does not depend on the time step.

```ini
[MOSE-Turbulence]
point-implicit = .true.   ; (default) point-implicit factor on the turbulence increment, steady runs
```

### Time-accurate runs: explicit sources

In time-accurate runs the factor is not applied, whatever `point-implicit` says; when `input.ini`
sets `point-implicit = .true.`, a warning at start-up recalls it (a run that leaves the key out
takes the default and gets no warning). Since $(F + S)/(1 + \Delta t\,d) = (F + S)(1 - \Delta t\,d) + O(\Delta t^2)$, the factor
would make the turbulence equations, and the mean flow through $\mu_t$, first-order accurate in
time. With the $k$–$\omega$ models the error of the factor would not even vanish when the mesh and the
step are refined together, as in an explicit run, whose global step goes as $h^2$ in the wall cells: in
the first cell $\omega \approx 6\nu/(\beta_1 y^2)$ gives $\Delta t\,d \approx 24.7$ `vnn` at every level (7.4 at
`vnn = 0.3`, a factor of 0.12), and the next two cells stay of order one. With Spalart–Allmaras
$\Delta t\,d \approx 1.64\,f_w\,y_1^+$ at `vnn = 0.3` ($y_1^+$ the height of the first cell in wall units): it
vanishes as the mesh is refined, but is of order one at the usual heights.
The sources are explicit, like the fluxes, and the local time step of each cell
([Time Integration](time-integration.md#local-time-step)) includes a bound of their stiffness:

$$
\frac{1}{\Delta t_i} = \bigl(\text{convective and diffusive rates}\bigr) + \hat d_i ,
\qquad \hat d_i \ge d_i = \max_q \frac{\partial D_q}{\partial q} ,
$$

with $\hat d$ a bound of the rates $d$ of the table that the state of the cell at the start of the step gives
alone, without gradients or ghost cells:

| Model | $\hat d$ | $\hat d / d$ |
|-------|----------|--------------|
| SST | $\max(\beta^\ast, 2\max(\beta_1, \beta_2))\,\omega$ | $\beta = F_1\beta_1 + (1 - F_1)\beta_2$, and $F_1$ takes gradients: $\beta_2/\beta_1 = 1.10$ where $F_1 = 1$ (next to walls), 1 where $F_1 = 0$ |
| Wilcox 2006 | $\max(\beta^\ast, 2\beta_0)\,\omega$ | $\beta = \beta_0 f_\beta$, and $f_\beta \le 1$ takes the velocity gradient: $1/f_\beta$, 1 in 2-D ($\chi_\omega = 0$), at most $1/0.85$ |
| SA | $2\,c_{w1}\,(1 + c_{w3}^6)^{1/6}\,\tilde\nu/y^2$ | $f_w = g\,[(1 + c_{w3}^6)/(g^6 + c_{w3}^6)]^{1/6}$ grows with $g \ge 0$ towards $(1 + c_{w3}^6)^{1/6} = 2.005$: $2.005/f_w$, 1 where $r$ reaches its cap of 10, about 2 where $r \approx 1$ (an equilibrium boundary layer, where the destruction is slow against the diffusion) |

($y$ the wall distance of the source, moved below a rough wall.) The Reynolds-stress model SSG/LRR,
which the key `turbulence` does not select, takes no bound: none was derived for its
sources, and its $\hat d$ is zero. Then
$\Delta t\,d \le \Delta t\,\hat d \le 1$: a forward-Euler stage cannot make
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
$\Delta t\,\hat d \le 1$ this holds as long as $d_s < 2\,\hat d$. The margin of 2 of $\tilde\nu$ and $\omega$
therefore covers a rate that grows within the step up to twice the bound at the start (at least twice the rate),
and no more; $k$, whose destruction is linear, keeps only the margin by which its rate $\beta^\ast\omega$ lies below
$\hat d$. An
impulsive start of a $k$–$\omega$ model from a state far from the wall value of $\omega$ goes beyond
that: during the first step the wall condition raises $\omega$ in the wall cells by orders of
magnitude, and with the full step the explicit destruction overshoots ($\omega < 0$, then the floor
and a very large $\mu_t$). Start such runs with the CFL ramp (`cfl-rise-threshold`).

At the end of a RANS run in which a stage made $\tilde\nu$, $k$ or $\omega$ negative before the
realizability floor, MOSE prints how many cell updates had such a stage; a run in which none did
prints no such line:

```text
 Turbulence variables below zero after a stage (set to the floor): <n>
```

In a steady run the count covers the whole run: the start, and the fronts of $k$ that the fluxes
push below zero where $k$ is about zero, which can last to the last iteration; it is a total, which
convergence does not bring back to zero. The transport keeps the turbulence variables non-negative,
with every limiter but van Albada, when at every stage the Courant number of the flow leaving each
cell is at most 1/2 ([Time Integration](time-integration.md#local-time-step)); a front can go below
zero with a larger step (a `cfl` above 0.5, such as the 0.8 suggested with RK3, or the faster waves at captured shocks),
through the sources, in the first cell next to a wall, whose bound is lower, and next to faces whose
ghost cells extrapolate the cells (types `408` and `502`), and with `irs-variables = all`, whose smoothing mixes the
updates of neighbouring cells. The residuals of `residual-history.dat` are the change of
each variable over a step ([Stopping criteria](../user/using.md#stopping-criteria)): a variable
that the floor holds at its value in a cell from one step to the next adds nothing there to its
residual, so these fronts do not show in the residual history, and the count is their trace.

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
