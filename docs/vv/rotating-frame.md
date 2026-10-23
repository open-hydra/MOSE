# Rotating Frame Verification

This page collects the rotating-frame verification cases that exercise the Coriolis and centrifugal source terms in isolation. Both tests verify the source-term coupling in [the SRF theory page](../theory/rotating-frame.md).

## Coriolis channel (`test/3D/coriolis-channel`)

A straight rectangular channel carries an axial throughflow of air in a frame rotating about an axis aligned with $\hat{\boldsymbol{x}}$. The rotation axis is offset transversely so that $\boldsymbol{r}_{\perp}$ has both an antisymmetric component across the cross-section (in $y$) and a streamwise component (in $z$); both Coriolis and centrifugal source terms are therefore active, and the case exercises the full body-force coupling rather than isolating one term.

**Geometry and grid.** Single block covering $[0,0.01]\times[0,0.02]\times[0,0.1]$ m with $5\times20\times30$ cells (uniform spacing).

**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0.005,\,0.010,\,0)$ m. With $\boldsymbol{r} = (x-0.005,\,y-0.01,\,z)$, the perpendicular component is $\boldsymbol{r}_{\perp} = (0,\,y-0.01,\,z)$: antisymmetric in $y$ across the cross-section but strictly positive in $z$ along the channel.

**Boundary conditions (from `test/3D/coriolis-channel/input.ini`).**

| Face | Location | Type | Role |
|---|---|---|---|
| 1–4 | lateral sides | `symmetry` | slip walls |
| 5 | inlet, $z=0$ | `inlet` | mass flux $g = 12.25$ kg m⁻² s⁻¹ at $T = 300$ K, normal to the face |
| 6 | outlet, $z=0.1$ m | `outlet` | static pressure $p = 1.054725\times10^{5}$ Pa |

The inflow boundary lies on the rotation axis and injects the flow in the normal direction following the rotation. The outlet is what sets the pressure level of this subsonic flow: with zero-gradient `extrapolation` on every open face nothing does, and the mean pressure rises without bound instead of settling.

**Initial state.** Uniform air at $T = 300$ K, $p = 1.054725\times10^{5}$ Pa, moving along the channel at $w = 10$ m/s (`INPUT/ic.tec`, written by ATLAS ICB from `[ICB-Block1]`).

**Solver.** Euler, RK2 with local time stepping, CFL $= 0.9$, MUSCL/minmod, HLLC.

**Expected behaviour.** With inflow along $+\hat{\boldsymbol{z}}$ and $\boldsymbol{\Omega} = \omega\hat{\boldsymbol{x}}$, the Coriolis acceleration $-2\,\boldsymbol{\Omega}\times\boldsymbol{w}$ points toward $+\hat{\boldsymbol{y}}$ (since $\hat{\boldsymbol{x}}\times\hat{\boldsymbol{z}} = -\hat{\boldsymbol{y}}$). The centrifugal contribution $\rho\omega^2 \boldsymbol{r}_{\perp}$ adds a streamwise body force ($+\hat{\boldsymbol{z}}$) that grows linearly with $z$ and a transverse component antisymmetric about the channel mid-plane $y = 0.01$. Between slip walls the throughflow cannot be deflected, so the steady solution is a uniform stream $w$ along $\hat{\boldsymbol{z}}$ whose pressure field balances the two body forces:

$$\frac{\partial p}{\partial y} = \rho\left[2\,\omega\,w + \omega^2\,(y - 0.01)\right], \qquad \frac{\partial p}{\partial z} = \rho\,\omega^2\,z .$$

The Coriolis term dominates the cross-channel gradient, $2\rho\omega w \approx 2450$ Pa/m, and the centrifugal term raises the pressure by $\rho\omega^2 z^2/2 \approx 61$ Pa along the channel. The case verifies sign and magnitude of both source terms together and that the inflow/outflow BCs remain stable when they are active.

**Result.** The case has no measuring script and no pass criterion: the numbers below are measured on its solution. The run stops when its density residual has fallen by the default threshold, $10^{-10}$ of its value at the first iteration, after 122 856 iterations with the GNU 12.2 RELEASE build; with gfortran 13 TESTING it stops after 101 892, since the residual falls while oscillating. The mass flux is 12.250 kg m⁻² s⁻¹ at the inlet and 12.255 at the outlet, with $w$ between 9.78 and 10.24 m/s and no cross-flow ($|v| < 0.1$ m/s away from the outlet). At mid length the pressure difference across the channel is 45.1 Pa, against 46.6 Pa from the balance above. Away from the walls, the cross-channel gradient is within 0.3 % of the balance over the first third of the channel and within 2.6 % over the middle third; the streamwise gradient is within 2 % from the third cell layer to two thirds of the length. Both fall short towards the outlet (by up to 20 % and 35 % in the last layers), where the imposed pressure is uniform and cannot hold the cross-channel gradient. The cross-channel gradient is also 6–10 % short in the two cell rows next to each wall. Over the region that leaves out 5 layers at the inlet and 8 at the outlet, the $L_2$ deviation from the balance is 5.0 % for $\partial p/\partial y$ and 3.3 % for $\partial p/\partial z$.

## Rotating annulus (`test/3D/pressure-centrifugal`)

A full annulus of quiescent air rotates about its axis. With $\boldsymbol{w}=\boldsymbol{0}$ initially the Coriolis term is inactive, and the centrifugal source $\rho\omega^2\boldsymbol{r}_{\perp}$ — together with its work contribution $\boldsymbol{w}\cdot\boldsymbol{S}_{\text{mom}}$ in the energy equation — drives the transient toward an isentropic rigid-body equilibrium that admits a closed-form pressure profile.

**Geometry and grid.** A single annular block aligned with the rotation axis $\hat{\boldsymbol{x}}$. In the $(y,z)$ plane the mesh wraps a full $2\pi$ azimuthally; the radial coordinate is $r_\perp = \sqrt{y^2 + z^2}$.

| Index | Direction | Extent | Cells |
|---|---|---|---|
| $i$ | axial ($x$) | $[0,\,0.02]$ m | 10 |
| $j$ | radial | $r_\perp \in [0.05,\,0.15]$ m | 20 |
| $k$ | azimuthal | $\theta \in [0,\,2\pi]$, $\Delta\theta = 18°$ per cell | 20 |


**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0,0,0)$ on the axis (which lies *outside* the fluid domain).

**Initial state.** Quiescent air at $T = 300$ K, $p = 1.054725\times 10^{5}$ Pa, $\boldsymbol{w}=\boldsymbol{0}$: the uniform field of `[ICB-Block1]`, which ATLAS ICB writes to `INPUT/ic.tec`.

**Boundary conditions (from `test/3D/pressure-centrifugal/input.ini`).**

| Face | Location | Type | Role |
|---|---|---|---|
| 1, 2 | axial ends ($x = 0,\,0.02$ m) | symmetry | axially homogeneous closure |
| 3, 4 | radial ends ($r_\perp = 0.05$ and $0.15$ m) | symmetry | inviscid wall (slip) at hub and shroud |
| 5, 6 | azimuthal periodic pair ($\theta = 0$) | periodic | closes the full ring: the two faces coincide, so the pair is a connection of the block with itself (type 201, no rotation). |

**Solver.** Euler, RK3 with local time stepping, CFL $= 2.0$ with implicit residual smoothing (`irs-beta = 0.5`), MUSCL/minmod, HLLC.

**Closed-form reference.** For inviscid steady rigid-body rotation with no heat transfer, the equilibrium is **isentropic** ($p/\rho^\gamma$ uniform). Combining $\nabla p = \rho\omega^2 \boldsymbol{r}_\perp$ with the isentropic relation gives

$$\frac{p(r_\perp)}{p_{\rm ref}} = \left[1 + \frac{\gamma-1}{2}\,\frac{\omega^2\bigl(r_\perp^2 - r_{\rm ref}^2\bigr)}{\gamma R\, T_{\rm ref}}\right]^{\gamma/(\gamma-1)}, \qquad \frac{T(r_\perp)}{T_{\rm ref}} = \left(\frac{p(r_\perp)}{p_{\rm ref}}\right)^{(\gamma-1)/\gamma}.$$

The reference state $(p_{\rm ref}, T_{\rm ref})$ at $r_{\rm ref}$ is fixed by total-mass conservation against the initial uniform field.

**Expected behaviour.** The centrifugal source builds a radially stratified field with $p$ increasing monotonically with $r_\perp$ at constant entropy, while $\boldsymbol{w}$ should relax back to zero. The case verifies that

- the centrifugal contribution is wired into both the momentum and energy residuals with the correct sign and magnitude,
- the periodic connection (`201`) closes the ring,
- the converged pressure and temperature profiles match the isentropic closed form above (radial $L_2$ error against the analytical reference is the verification metric).

**Result.** The case has no measuring script and no pass criterion: the numbers below are measured on its solution, with the volume-weighted relative $L_2$ error of $p$ and $T$ as the metric named above. From the uniform start the density residual rises from $1.3\times10^{-4}$ at the first iteration to $1.7\times10^{-3}$ at the 51st, falls to $10^{-8}$ in 40 000 iterations and then only creeps ($4.4\times10^{-9}$ at 200 000). It does not reach the default threshold, $10^{-10}$ of its first value, so the run has to be stopped with `iter-threshold`. At 200 000 iterations:

| Quantity | MOSE | Closed form |
|---|---|---|
| $p$, range over the cells | 105 413 – 105 530 Pa | 105 414 – 105 528 Pa |
| $T$, range over the cells | 299.951 – 300.059 K | 299.953 – 300.045 K |
| Pressure rise hub to shroud | 112.7 Pa | 113.5 Pa |
| Volume-weighted relative $L_2$ error | $p$: $1.6\times10^{-5}$, $T$: $1.4\times10^{-5}$ | – |
| Total mass against the initial one | $+2.4\times10^{-6}$ | 0 |

The fluid does not come fully to rest. A tangential velocity remains in the cell rows next to hub and shroud, up to 1.66 m/s at the shroud against a rim speed $\omega r_{\rm out} = 15$ m/s, and about 0.05 m/s in the core; it grows slowly with the iterations and goes with the 0.7 % deficit in the pressure rise.

**Sensitivity note.** At $\omega = 100$ rad/s the rim Mach number is $\omega r_{\rm out}/\sqrt{\gamma R T} \approx 0.04$, giving a pressure variation of only $\sim 0.1\%$ across the annulus. Verification at this setting requires tight residual convergence and a precise mass-conservation pinning of $p_{\rm ref}$.