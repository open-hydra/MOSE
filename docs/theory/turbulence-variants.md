# Turbulence Model Variants and Corrections

This page details the variants and corrections that can be applied to the [turbulence models](turbulence.md): the Spalart–Allmaras
compressibility, rotation and rotation-curvature corrections, and for SST the compressible production, the dilatation dissipation,
the rotation-curvature correction and the coupling of $k$ with the mean-flow energy.

## SA Variants

### Compressibility Correction

The correction of Paciorri and Sabetta (J. Spacecraft Rockets 40(3), 2003,
doi:10.2514/2.3967) for the reduced growth of compressible mixing layers.
The production term is multiplied by $f_1(\tilde M_c)\,f_2(\tilde M_c)$, with

$$
f_1(M_c) = \frac{0.6}{1 + 9 M_c^6} + 0.4, \qquad
f_2(M_c) = \frac{0.44}{1 + 14 M_c^5} + 0.56
$$

(their eqs. 5 and 6: the growth-rate and shear-stress ratios of a
compressible mixing layer to the incompressible one, fitted to measurements
against the convective Mach number $M_c$). The local estimate $\tilde M_c$
solves their eq. 12,

$$
\tilde M_c^2\, f_2(\tilde M_c) = \frac{\nu_t\,|\omega|}{4\,\tau_i\,a^2}
= 25\,\frac{\tilde\nu f_{v1}\,|\omega|}{a^2},
$$

with $\tau_i = 0.01$ the shear-stress ratio of the incompressible model,
$\nu_t = \tilde\nu f_{v1}$ the eddy viscosity, $|\omega|$ the vorticity and
$a$ the local speed of sound; both sides are numbers. The code solves it by
Newton's method from $\tilde M_c = 1$ (the left side grows with $\tilde M_c$:
one root). At low $\tilde M_c$ the factor is 1; it falls to 0.22 at
$\tilde M_c \gtrsim 1.5$.

!!! note "Boundary layers"
    The paper fits mixing layers and warns that the correction alters the
    compressible boundary layer at high Mach number, where $\tilde M_c$ is
    large too; it gives no switch, and MOSE applies the correction in every
    cell. Before 2026-11 the code divided by $a$ instead of $a^2$: the factor
    then depended on the unit of length and cut the production wherever the
    shear was strong, whatever the Mach number.

### SA-R — Rotation Correction

SA-R replaces $\tilde S$ in the production by $\tilde S + C_\text{rot}\min(0,\,S - \Omega)$,
with the strain-rate magnitude $S = \sqrt{2\,S_{ij}\,S_{ij}}$:

$$
P_\text{SA-R} = c_{b1}\,\bigl[\tilde{S} + C_\text{rot}\,\min(0,\,S - \Omega)\bigr]\,\rho\tilde\nu
\qquad (C_\text{rot} = 2)
$$

The production is reduced where the vorticity exceeds the strain rate, as in vortex cores, and
unchanged where it does not (NASA Turbulence Modeling Resource, Spalart–Allmaras page, SA-R, after
Dacles-Mariani et al., *AIAA J.* 33(9), 1995).

### SA-RC — Spalart–Shur Rotation-Curvature Correction

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
correction are outside its validity.

MOSE takes $\partial S_{ij}/\partial x_k$ with central differences of the velocity gradients
$G = \nabla\mathbf{u}$ of the cell and of its neighbours, mapped to Cartesian components with a cell metric. Next to a block face the
neighbour is the first ghost cell, and the difference is right when the ghost cell holds the value
of a definite cell and the metric takes the position of that same cell: for $S_{ij}$ linear in
space the difference is then exact (with the stencil metric, the metric of the
differences of these positions: [Cell gradients of the source terms](turbulence-numerics.md#cell-gradients-of-the-source-terms)).

- Behind a face between blocks or of a periodic pair the ghost cell holds $G$ of the
  cell across (from another MPI rank by message), and the metric
  takes its centre, plus the translation of a periodic pair. $G$ is a tensor in Cartesian
  components: it is copied whatever the orientation of the two blocks. The cell next to the face
  then takes the central difference of the undivided mesh, and the result does not depend on how
  the mesh is divided into blocks or among ranks.
- Behind a face of rotational periodicity, the ghost cell holds $R\,G\,R^{\mathsf T}$, $G$ of
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

## Compressible production

The production of the [SST model](turbulence.md#production-limiter) is the SST-2003m form $\mu_t S^2$, which drops the
dilatation terms of $P = \tau_{ij}\,\partial u_i/\partial x_j$.  The
compressible production (not the default) keeps them, with
$\tau_{ij} = \mu_t\bigl(2S_{ij} - \tfrac23\,\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \tfrac23\,\rho k\,\delta_{ij}$:

$$
P_k = \min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\rho k\,\nabla\!\cdot\!\mathbf{v},
$$

$$
P_\omega = \frac{\gamma\,\rho}{\mu_t}\,\min\!\Bigl(\mu_t\bigl[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2\bigr],\;10\,\beta^\ast\rho\,\omega\,k\Bigr)
      - \tfrac23\,\gamma\,\rho\,\omega\,\nabla\!\cdot\!\mathbf{v}.
$$

This is the exact production, with the dilatation terms arranged as in
OpenFOAM's `kOmegaSST` [9]: the deviatoric
part, $\mu_t[S^2 - \tfrac23(\nabla\!\cdot\!\mathbf{v})^2] = 2\mu_t S^{d}_{ij}S^{d}_{ij} \ge 0$,
stays under the limiter, and $-\tfrac23\rho k\,\nabla\!\cdot\!\mathbf{v}$ ($k$ equation) and
$-\tfrac23\gamma\rho\omega\,\nabla\!\cdot\!\mathbf{v}$ ($\omega$ equation) are added outside it.  TMR's
SST-2003e and SST-2003 limit the whole of $P$ and take $\gamma/\nu_t$ times the
limited value as the production of $\omega$, whose term in $\rho k\,\nabla\!\cdot\!\mathbf{v}$, below the
limiter, is $-\tfrac23\gamma\rho(k/\nu_t)\,\nabla\!\cdot\!\mathbf{v}$.  The two
arrangements coincide where $\nabla\!\cdot\!\mathbf{v} = 0$, and elsewhere where the
production of neither reaches the limiter and $a_1\omega \ge S F_2$, so that
$k/\nu_t = \omega$.

Without the energy coupling (the default) the
$-\tfrac23\rho k\,\delta_{ij}$ part of $\tau_{ij}$ enters the production but not the
momentum and energy equations, as in TMR's "e" variants.  In attached, nearly
incompressible boundary layers the compressible and the incompressible productions
coincide; through a compression ($\nabla\!\cdot\!\mathbf{v}<0$:
shock feet, reattachment) the dilatation terms amplify $k$ more than $\omega$
($\gamma<1$), so $\mu_t$ grows.  In the Mach-5
[shock/boundary-layer interaction](../vv/2D-swbli.md) the incompressible form gives
a separation bubble about 16% longer and a slower recovery of $C_f$.

## Dilatation dissipation (`sst-dilatation-dissipation`)

A compressibility correction of the dissipation, for the reduced growth of
compressible shear layers (Rumsey, NASA/TM-2009-215705, eqs. 40–43, after
Wilcox and Zeman). The destruction coefficients of $k$ and $\omega$ become

$$
\beta^\ast_c = \beta^\ast\,[1 + \xi^\ast F(M_T)], \qquad
\beta_c = \beta - \beta^\ast\,\xi^\ast F(M_T), \qquad
M_T = \frac{\sqrt{2k}}{a},
$$

with $M_T$ the turbulent Mach number and $a$ the local speed of sound.
`sst-dilatation-dissipation` (default `none`) selects $F$:

| Value | $F(M_T)$ | Constants | Meant for |
|---|---|---|---|
| `wilcox` | $(M_T^2 - M_{T0}^2)\,\mathcal H(M_T - M_{T0})$ | $\xi^\ast = 2$, $M_{T0} = 0.25$ | jets and mixing layers (Sarkar/Zeman/Wilcox family) |
| `zeman` | $\bigl[1 - \exp(-((M_T - M_{T0})/\Lambda)^2)\bigr]\,\mathcal H(M_T - M_{T0})$ | $\xi^\ast = 0.75$, $M_{T0} = 0.2$, $\Lambda = 0.66$ | boundary layers (milder) |

Both raise the dissipation of $k$ and lower that of $\omega$ where $M_T$ is
above the threshold: less eddy viscosity, a shear layer that spreads more
slowly. Rumsey's survey: the Wilcox form gives the measured trend of mixing
layers and over-corrects high-Mach boundary layers (skin friction too low);
the Zeman form improves cold-wall hypersonic skin friction and is too weak
for free shear layers. In a boundary layer $M_T$ stays below 0.25 up to
about Mach 4 (his Fig. 7b), so the Wilcox form leaves it alone there.

MOSE floors $\beta_c$ at 0 (not in the papers): with the Wilcox constants
$\beta_c$ would turn negative for $M_T > 0.65$ and feed $\omega$. The
point-implicit factors of steady runs take $\beta^\ast_c$ and $\beta_c$.
With `none` the coefficients are $\beta^\ast$ and $\beta$ exactly. The
option acts on `SST` and `SST-RC`, not on `Wilcox2006`; the multirate integrator
carries it too (its cell-wise copy of the SST sources, checked by `verify_copies.py`).

## Rotation/curvature correction (SST-RC)

SST-RC sensitizes the model to streamline curvature and system
rotation with the Spalart–Shur correction as adapted to SST by Smirnov and
Menter [10] (TMR SST-RC [8], here on the SST-2003 variants of the model).  The production $P$ of
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
$P$; with the compressible production it multiplies the deviatoric part of
$P$, and the dilatation terms, added after the limiter, are not multiplied.  In
simple shear $r^\ast = 1$ and $\tilde r = 0$, so $f_{r1} = 1$ and the
model reduces to SST.

Implementation notes: $r^\ast$ and the numerator of $\tilde r$ are those of
SA-RC.  $DS_{ij}/Dt$ is its steady convective part
$\mathbf{u}\!\cdot\!\nabla S_{ij}$; in a [rotating frame](rotating-frame.md)
the frame enters through the absolute vorticity in $\Omega_{ik}$ and $\Omega$,
but the frame-rotation term $(\epsilon_{imn}S_{jn} + \epsilon_{jmn}S_{in})\,\Omega^\text{rot}_m$
of the published $\tilde r$ is not included.  $\Omega$ is floored at
$10^{-30}$ s⁻¹, so cells with no vorticity take the limit $f_{r1} = 1.25$
instead of an undefined value (their production is negligible there).

## Energy coupling

The compressible RANS equations [11] carry the turbulent kinetic energy in two
places: the total energy, $E = e + \tfrac12 u_iu_i + k$, and the Boussinesq
Reynolds stress,
$\tau_{ij} = 2\mu_t\bigl(S_{ij} - \tfrac13\nabla\!\cdot\!\mathbf{v}\,\delta_{ij}\bigr) - \tfrac23\rho k\,\delta_{ij}$.
By default MOSE drops both (the TMR "m"/"e" variants): its energy is that of the
mean flow, $E_m = e + \tfrac12 u_iu_i$, and the $-\tfrac23\rho k$ term is not in
the stress.  With the energy coupling (two-equation models: SST and
Wilcox 2006) MOSE keeps both:

- the stress of the momentum and energy equations includes
  $-\tfrac23\rho k\,\delta_{ij}$, so it also enters the work term
  $\nabla\!\cdot\!(\boldsymbol{\tau}\cdot\mathbf{v})$ of the energy equation;
- the mean-flow energy gives up what the turbulence receives,

$$
\frac{\partial(\rho E_m)}{\partial t} + \dots = \dots - \bigl(P_k - \beta^\ast\rho\,\omega\,k\bigr),
$$

with the same source as the $k$ equation, so $E_m + k$ is conserved.  In steady
runs the [point-implicit](turbulence-numerics.md#steady-runs-point-implicit-factor) factor treats the
energy and $k$ equations as one system: the two increments add up to that of the
fluxes of $E_m + k$ at every stage, and the steady state does not depend on the
time step.  This is the total-energy form of TMR [11] and Wilcox [4] written for
$E_m$: subtracting the $k$ equation from the equation for $E = E_m + k$ leaves
exactly this sink, and the diffusion of $k$ in the two equations cancels.

The
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

The $-\tfrac23\rho k\,\delta_{ij}$ part of the stress is a pressure in all but
name, and is applied wherever the pressure is: besides the interior faces, on the
inflow, outflow and extrapolation faces (from the boundary cell) and, on an
axisymmetric mesh, on the two side faces of the wedge, where it adds
$\tfrac23\rho k$ to the pressure in the radial momentum balance.  At a wall
$k = 0$ and the term vanishes.

---

## References

[4] D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries, 2006.

[8] NASA Langley Turbulence Modeling Resource, "The Menter Shear Stress Transport turbulence model," <https://tmbwg.github.io/turbmodels/sst.html> (formerly turbmodels.larc.nasa.gov/sst.html): definitions and naming of the SST variants (SST-2003, SST-2003m/s/e, SST-RC) and of Wilcox2006-klim.

[9] OpenFOAM Foundation, OpenFOAM-10, `kOmegaSST` model (`src/MomentumTransportModels/momentumTransportModels/Base/kOmegaSST/kOmegaSSTBase.C`).

[10] P. E. Smirnov, F. R. Menter, "Sensitization of the SST turbulence model to rotation and curvature by applying the Spalart–Shur correction term," *J. Turbomach.*, 131(4), 041010, 2009, doi:10.1115/1.3070573.

[11] NASA Langley Turbulence Modeling Resource, "Implementing turbulence models into the compressible RANS equations," <https://tmbwg.github.io/turbmodels/implementrans.html>.
