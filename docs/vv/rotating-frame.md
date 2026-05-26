# Rotating Frame Verification

This page collects the rotating-frame verification cases that exercise the Coriolis and centrifugal source terms in isolation. Both tests verify the source-term coupling in [the SRF theory page](../theory/rotating-frame.md).

## Coriolis channel (`test/3D/coriolis-channel`)

A straight rectangular channel carries an axial throughflow of air in a frame rotating about an axis aligned with $\hat{\boldsymbol{x}}$. The rotation axis is offset transversely so that $\boldsymbol{r}_{\perp}$ has both an antisymmetric component across the cross-section (in $y$) and a streamwise component (in $z$); both Coriolis and centrifugal source terms are therefore active, and the case exercises the full body-force coupling rather than isolating one term.

**Geometry and grid.** Single block covering $[0,0.01]\times[0,0.02]\times[0,0.1]$ m with $5\times20\times30$ cells (uniform spacing).

**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0.005,\,0.010,\,0)$ m. With $\boldsymbol{r} = (x-0.005,\,y-0.01,\,z)$, the perpendicular component is $\boldsymbol{r}_{\perp} = (0,\,y-0.01,\,z)$: antisymmetric in $y$ across the cross-section but strictly positive in $z$ along the channel.

**Boundary conditions (from `test/3D/coriolis-channel/input.ini`).**

| Face | Location | Preset | Type |
|---|---|---|---|
| 1–4 | lateral sides | 11 | outflow |
| 5 | inlet, $z=0$ | 4 | inflow (subsonic with total conditions) |
| 6 | outlet, $z=0.1$ m | 11 | outflow |

The inlet section sets `T = 300`, `g = 12.25`, `yAir = 1.0`. The numeric type IDs above are the user-side preset codes from `input.ini`; faces 1–4 and face 6 share the same preset (11). The inflow boundary lays on the rotation axis and injects the flow in the normal direction following the rotation.

**Solver.** Euler, RK2 with local time stepping, CFL $= 0.9$, MUSCL/minmod, HLLC.

**Expected behaviour.** With inflow along $+\hat{\boldsymbol{z}}$ and $\boldsymbol{\Omega} = \omega\hat{\boldsymbol{x}}$, the Coriolis acceleration $-2\,\boldsymbol{\Omega}\times\boldsymbol{w}$ deflects the throughflow toward $+\hat{\boldsymbol{y}}$ (since $\hat{\boldsymbol{x}}\times\hat{\boldsymbol{z}} = -\hat{\boldsymbol{y}}$). The centrifugal contribution $\rho\omega^2 \boldsymbol{r}_{\perp}$ adds a streamwise body force ($+\hat{\boldsymbol{z}}$) that grows linearly with $z$ and a transverse component antisymmetric about the channel mid-plane $y = 0.01$. The case verifies sign and magnitude of both source terms together and that the inflow/outflow BCs remain stable when they are active.

## Rotating annulus (`test/3D/pressure-centrifugal`)

A full annulus of quiescent air rotates about its axis. With $\boldsymbol{w}=\boldsymbol{0}$ initially the Coriolis term is inactive, and the centrifugal source $\rho\omega^2\boldsymbol{r}_{\perp}$ — together with its work contribution $\boldsymbol{w}\cdot\boldsymbol{S}_{\text{mom}}$ in the energy equation — drives the transient toward an isentropic rigid-body equilibrium that admits a closed-form pressure profile.

**Geometry and grid.** A single annular block aligned with the rotation axis $\hat{\boldsymbol{x}}$. In the $(y,z)$ plane the mesh wraps a full $2\pi$ azimuthally; the radial coordinate is $r_\perp = \sqrt{y^2 + z^2}$.

| Index | Direction | Extent | Cells |
|---|---|---|---|
| $i$ | axial ($x$) | $[0,\,0.02]$ m | 10 |
| $j$ | radial | $r_\perp \in [0.05,\,0.15]$ m | 20 |
| $k$ | azimuthal | $\theta \in [0,\,2\pi]$, $\Delta\theta = 18°$ per cell | 20 |


**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0,0,0)$ on the axis (which lies *outside* the fluid domain).

**Initial state.** Quiescent air at $T = 300$ K, $p = 1.054725\times 10^{5}$ Pa, $\boldsymbol{w}=\boldsymbol{0}$.

**Boundary conditions (from `test/3D/pressure-centrifugal/input.ini`).**

| Face | Location | Type | Role |
|---|---|---|---|
| 1, 2 | axial ends ($x = 0,\,0.02$ m) | symmetry | axially homogeneous closure |
| 3, 4 | radial ends ($r_\perp = 0.05$ and $0.15$ m) | symmetry | inviscid wall (slip) at hub and shroud |
| 5, 6 | azimuthal periodic pair ($\theta = 0$) | periodic | closes the full ring via `src/lib/numerics/fluxes/bc/Lib_BC_Fluxes_Rotational.f90` since the domain is a curved rectangle. |

**Solver.** Euler, RK3 with local time stepping, CFL $= 2.0$, MUSCL/minmod, HLLC.

**Closed-form reference.** For inviscid steady rigid-body rotation with no heat transfer, the equilibrium is **isentropic** ($p/\rho^\gamma$ uniform). Combining $\nabla p = \rho\omega^2 \boldsymbol{r}_\perp$ with the isentropic relation gives

$$\frac{p(r_\perp)}{p_{\rm ref}} = \left[1 + \frac{\gamma-1}{2}\,\frac{\omega^2\bigl(r_\perp^2 - r_{\rm ref}^2\bigr)}{\gamma R\, T_{\rm ref}}\right]^{\gamma/(\gamma-1)}, \qquad \frac{T(r_\perp)}{T_{\rm ref}} = \left(\frac{p(r_\perp)}{p_{\rm ref}}\right)^{(\gamma-1)/\gamma}.$$

The reference state $(p_{\rm ref}, T_{\rm ref})$ at $r_{\rm ref}$ is fixed by total-mass conservation against the initial uniform field.

**Expected behaviour.** The centrifugal source builds a radially stratified field with $p$ increasing monotonically with $r_\perp$ at constant entropy, while $\boldsymbol{w}$ relaxes back to zero. The case verifies that

- the centrifugal contribution is wired into both the momentum and energy residuals with the correct sign and magnitude,
- the rotational-periodic BC closes the ring with the proper velocity rotation,
- the converged pressure and temperature profiles match the isentropic closed form above (radial $L_2$ error against the analytical reference is the verification metric). The initial solution in **ic.tec** stores pressure profile that should balance the centrifugal contribution, i.e. the final solution should be similar to the initial one. (Residuals start very low)

**Sensitivity note.** At $\omega = 100$ rad/s the rim Mach number is $\omega r_{\rm out}/\sqrt{\gamma R T} \approx 0.04$, giving a pressure variation of only $\sim 0.1\%$ across the annulus. Verification at this setting requires tight residual convergence and a precise mass-conservation pinning of $p_{\rm ref}$; running at higher $\omega$ (e.g. 1000 rad/s, $\sim 10\%$ variation) gives the comparison more dynamic range.