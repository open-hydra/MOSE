# Single Rotating Frame (SRF) Theory

## Overview

The Single Rotating Frame (SRF) model allows simulation of the compressible Navier–Stokes equations in a reference frame rotating rigidly at constant angular velocity $\boldsymbol{\Omega}$ [rad/s] about a fixed axis. This is particularly useful for turbomachinery applications where the computational domain rotates with an impeller or rotor.

## Velocity Definition

The velocity field stored in the solution vector is the **relative (rotating-frame) velocity**:

$$\boldsymbol{w} = \boldsymbol{u}_{\text{abs}} - \boldsymbol{\Omega} \times \boldsymbol{r}$$

where:
- $\boldsymbol{u}_{\text{abs}}$ is the absolute (inertial-frame) velocity
- $\boldsymbol{\Omega}$ is the angular velocity vector
- $\boldsymbol{r}$ is the position vector from the axis of rotation

This formulation ensures that no modifications are required to the standard convective or diffusive flux modules—they operate directly on the relative velocity without change.

## Governing Equations in the Rotating Frame

### Continuity and Species Equations

The continuity and species mass fraction equations are **unchanged** from the inertial frame:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{w}) = 0$$

$$\frac{\partial (\rho Y_k)}{\partial t} + \nabla \cdot (\rho Y_k \boldsymbol{w}) = -\nabla \cdot \boldsymbol{J}_k + \dot{m}_k$$

### Momentum Equation

The momentum equation includes additional source terms accounting for Coriolis and centrifugal accelerations (per unit volume):

$$\frac{\partial (\rho \boldsymbol{w})}{\partial t} + \nabla \cdot (\rho \boldsymbol{w} \boldsymbol{w}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \boldsymbol{S}_{\text{mom}}$$

where the source term is:

$$\boldsymbol{S}_{\text{mom}} = -2\rho (\boldsymbol{\Omega} \times \boldsymbol{w}) + \rho \omega^2 \boldsymbol{r}_{\perp}$$

**Components:**

- **Coriolis force:** $-2\rho (\boldsymbol{\Omega} \times \boldsymbol{w})$
- **Centrifugal force:** $\rho \omega^2 \boldsymbol{r}_{\perp}$

where $\boldsymbol{r}_{\perp} = \boldsymbol{r} - (\boldsymbol{r} \cdot \boldsymbol{e}) \boldsymbol{e}$ is the position vector perpendicular to the rotation axis (with $\boldsymbol{e}$ the unit vector along the axis), and $\omega = |\boldsymbol{\Omega}|$.

### Energy Equation

The energy equation (with specific total energy $E = c_v T + \frac{1}{2} |\boldsymbol{w}|^2$) receives a source term:

$$\frac{\partial (\rho E)}{\partial t} + \nabla \cdot (\rho E \boldsymbol{w}) = -\nabla \cdot p\boldsymbol{w} + \nabla \cdot (\boldsymbol{\tau} \cdot \boldsymbol{w}) - \nabla \cdot \boldsymbol{q} + \boldsymbol{w} \cdot \boldsymbol{S}_{\text{mom}}$$

The work done by body forces simplifies to:

$$\boldsymbol{S}_{E} = \boldsymbol{w} \cdot (\rho \omega^2 \boldsymbol{r}_{\perp})$$

**Note:** The Coriolis force does no work since $\boldsymbol{w} \cdot (\boldsymbol{\Omega} \times \boldsymbol{w}) = 0$ identically.

### RANS Turbulence Equations

The transport equations themselves are **unchanged** from the inertial frame, but every closure that consumes the vorticity (rather than the strain rate) must use the **absolute** vorticity. For constant $\boldsymbol{\Omega}$ the kinematic identity

$$\boldsymbol{\omega}_{\rm abs} = \boldsymbol{\omega}_{\rm rel} + 2\boldsymbol{\Omega}$$

follows from $\nabla\times(\boldsymbol{\Omega}\times\boldsymbol{r}) = 2\boldsymbol{\Omega}$, and MOSE applies it wherever vorticity enters a closure:

| Model | Vorticity-dependent term |
|---|---|
| Spalart–Allmaras (baseline) | $\tilde S = \|\boldsymbol{\omega}\|$ in the production term |
| SA + Spalart–Shur (SA-RC) | $\tilde r$ in the rotation/curvature factor $f_{r1}$ |
| Menter SST $k$–$\omega$ | $O = \sqrt{2 W_{ij}W_{ij}}$ in the SST-RC denominator and $f_{r1}$ |
| Wilcox 2006 $k$–$\omega$ | $f_\beta$ vortex-stretching limiter, $W_{ij}W_{jk}\hat S_{ki}$ |
| SSG–LRR (RSM) | $W_{ij}^{\rm abs}$ in the $C_5$ pressure-strain term |

**Strain-rate-based production is frame-invariant.** The $S_{ij}$-driven production used by SST and Wilcox 2006 needs no correction — only the vorticity-dependent auxiliary terms above receive the $2\boldsymbol{\Omega}$ shift.

**Caveat for baseline SA at high $\Omega$.** Because $|\boldsymbol{\omega}_{\rm abs}| = 2|\boldsymbol{\Omega}|$ in a region of pure rigid-body co-rotation, baseline SA produces a non-zero $\tilde S$ there and will generate eddy viscosity from the frame rotation alone. This is the correct inertial-frame answer (the absolute vorticity really is non-zero in that region), but in flows where streamline curvature does not enhance turbulence the production is spurious. SA-RC's $f_{r1}$ is designed to attenuate exactly this contribution — so the recommended pairing in a rotating frame is **SA-RC**, not baseline SA.

**References.**

- Speziale, C.G. (1989) *"Turbulence modeling in noninertial frames of reference,"* Theoretical and Computational Fluid Dynamics **1**, 3–19. — General statement of why vorticity-dependent closures require the absolute vorticity in a rotating frame.
- Spalart, P.R. and Shur, M. (1997) *"On the sensitization of turbulence models to rotation and curvature,"* Aerospace Science and Technology **1**(5), 297–302. — SA-RC; absolute vorticity in $\tilde r$.
- Smirnov, P.E. and Menter, F.R. (2009) *"Sensitization of the SST turbulence model to rotation and curvature by applying the Spalart–Shur correction term,"* ASME Journal of Turbomachinery **131**, 041010. — SST-RC.
- Wilcox, D.C. (2006) *Turbulence Modeling for CFD*, 3rd ed., DCW Industries. — Frame indifference and rotating-frame closures, §2.4.

## Viscous Fluxes

The stress tensor in the rotating frame uses the relative velocity:

$$\boldsymbol{\tau} = \mu \left( \nabla\boldsymbol{w} + (\nabla\boldsymbol{w})^T - \frac{2}{3}(\nabla \cdot \boldsymbol{w}) \mathbf{I} \right)$$

The viscous energy flux is computed as $\boldsymbol{\tau} \cdot \boldsymbol{w}$, and the heat conduction flux remains $\boldsymbol{q} = -k \nabla T$.

Since all quantities already use the relative velocity stored in the solution, **no modifications are required** to standard diffusive flux routines.

## Boundary Conditions

### Standard Boundaries
Inflow/outflow and other standard boundary conditions are applied directly to the relative velocity field $\boldsymbol{w}$ without modification. When defining the inflow properties it is necessary to specify them in the **laboratory frame of reference**. The correction in the rotating frame is performed automatically by [`RF_Convert_BC_Inflow`](../../src/lib/physics/Lib_RotatingFrame.f90).

### Wall Boundaries
On rotating walls, the no-slip condition applies to the **absolute velocity**. When the wall rotates with the frame at angular velocity $\boldsymbol{\Omega}$, the relative velocity of fluid at a wall point is:

$$\boldsymbol{w}_{\text{wall}} = \boldsymbol{u}_{\text{abs, wall}} - (\boldsymbol{\Omega} \times \boldsymbol{r}_{\text{wall}})$$

For a perfectly rotating solid body, $\boldsymbol{u}_{\text{abs, wall}} = \boldsymbol{\Omega} \times \boldsymbol{r}_{\text{wall}}$, hence:

$$\boldsymbol{w}_{\text{wall}} = \boldsymbol{0}$$

(i.e., zero relative velocity at a wall moving with the frame).

## Physical Interpretation

The rotating-frame formulation is equivalent to solving the full inertial-frame problem but with fictitious body forces representing:

1. **Coriolis acceleration:** Deflects moving fluid perpendicular to both the rotation axis and the velocity vector.
2. **Centrifugal acceleration:** Acts radially outward from the rotation axis.

For steady-state turbomachinery flows, a steady rotating-frame solution avoids the need to track time-periodic structures in the inertial frame, enabling efficient steady-state analysis.

## Implementation Notes

- The computational domain is assumed to rotate as a **rigid body** at constant $\boldsymbol{\Omega}$.
- The axis and angular velocity magnitude/direction must be specified at initialization.
- The origin of the coordinate system (where $\boldsymbol{r} = \boldsymbol{0}$) should coincide with the axis of rotation.
- All flux routines (convective and viscous) operate on the relative velocity without modification.

## Verification cases

Two 3D tests under `test/3D/` exercise the Coriolis and centrifugal source terms in isolation, so that each contribution to $\boldsymbol{S}_{\text{mom}}$ can be verified independently against a known balance.

### Coriolis channel (`test/3D/coriolis-channel`)

A straight rectangular channel carries an axial throughflow of air in a frame rotating about an axis aligned with $\hat{\boldsymbol{x}}$. The rotation axis is offset transversely so that $\boldsymbol{r}_{\perp}$ has both an antisymmetric component across the cross-section (in $y$) and a streamwise component (in $z$); both Coriolis and centrifugal source terms are therefore active, and the case exercises the full body-force coupling rather than isolating one term.

**Geometry and grid.** Single block covering $[0,0.01]\times[0,0.02]\times[0,0.1]$ m with $5\times20\times30$ cells (uniform spacing).

**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0.005,\,0.010,\,0)$ m. With $\boldsymbol{r} = (x-0.005,\,y-0.01,\,z)$, the perpendicular component is $\boldsymbol{r}_{\perp} = (0,\,y-0.01,\,z)$: antisymmetric in $y$ across the cross-section but strictly positive in $z$ along the channel.

**Boundary conditions (from [input.ini](../../test/3D/coriolis-channel/input.ini)).**

| Face | Location | type |
|---|---|---|
| 1–4 | lateral sides | `outflow` | 11 |
| 5 | inlet, $z=0$ | `inflow` | 4 (subsonic inflow with total conditions) |
| 6 | outlet, $z=0.1$ m | `outflow` | 11 |

The inlet section sets `T = 300`, `g = 12.25`, `yAir = 1.0`. The numeric type IDs above are the user-side preset codes from `input.ini`; faces 1–4 and face 6 share the same preset (11). The inflow boundary lays on the rotation axis and injects the flow in the normal direction following the rotation.

**Solver.** Euler, RK2 with local time stepping, CFL $= 0.9$, MUSCL/minmod, HLLC.

**Expected behaviour.** With inflow along $+\hat{\boldsymbol{z}}$ and $\boldsymbol{\Omega} = \omega\hat{\boldsymbol{x}}$, the Coriolis acceleration $-2\,\boldsymbol{\Omega}\times\boldsymbol{w}$ deflects the throughflow toward $+\hat{\boldsymbol{y}}$ (since $\hat{\boldsymbol{x}}\times\hat{\boldsymbol{z}} = -\hat{\boldsymbol{y}}$). The centrifugal contribution $\rho\omega^2 \boldsymbol{r}_{\perp}$ adds a streamwise body force ($+\hat{\boldsymbol{z}}$) that grows linearly with $z$ and a transverse component antisymmetric about the channel mid-plane $y = 0.01$. The case verifies sign and magnitude of both source terms together and that the inflow/outflow BCs remain stable when they are active.

### Rotating annulus (`test/3D/pressure-centrifugal`)

A full annulus of quiescent air rotates about its axis. With $\boldsymbol{w}=\boldsymbol{0}$ initially the Coriolis term is inactive, and the centrifugal source $\rho\omega^2\boldsymbol{r}_{\perp}$ — together with its work contribution $\boldsymbol{w}\cdot\boldsymbol{S}_{\text{mom}}$ in the energy equation — drives the transient toward an isentropic rigid-body equilibrium that admits a closed-form pressure profile.

**Geometry and grid.** A single annular block aligned with the rotation axis $\hat{\boldsymbol{x}}$. In the $(y,z)$ plane the mesh wraps a full $2\pi$ azimuthally; the radial coordinate is $r_\perp = \sqrt{y^2 + z^2}$.

| Index | Direction | Extent | Cells |
|---|---|---|---|
| $i$ | axial ($x$) | $[0,\,0.02]$ m | 10 |
| $j$ | radial | $r_\perp \in [0.05,\,0.15]$ m | 20 |
| $k$ | azimuthal | $\theta \in [0,\,2\pi]$, $\Delta\theta = 18°$ per cell | 20 |


**Rotating frame.** $\omega = 100$ rad/s, axis $\hat{\boldsymbol{x}}$, origin $(0,0,0)$ on the axis (which lies *outside* the fluid domain).

**Initial state.** Quiescent air at $T = 300$ K, $p = 1.054725\times 10^{5}$ Pa, $\boldsymbol{w}=\boldsymbol{0}$.

**Boundary conditions (from [input.ini](../../test/3D/pressure-centrifugal/input.ini)).**

| Face | Location | Type | Role |
|---|---|---|---|
| 1, 2 | axial ends ($x = 0,\,0.02$ m) | symmetry | axially homogeneous closure |
| 3, 4 | radial ends ($r_\perp = 0.05$ and $0.15$ m) | symmetry | inviscid wall (slip) at hub and shroud |
| 5, 6 | azimuthal periodic pair ($\theta = 0$) | periodic | closes the full ring via [`Lib_BC_Fluxes_Rotational.f90`](../../src/lib/numerics/fluxes/bc/Lib_BC_Fluxes_Rotational.f90) since the domain is a curved rectangle. |

**Solver.** Euler, RK3 with local time stepping, CFL $= 2.0$, MUSCL/minmod, HLLC.

**Closed-form reference.** For inviscid steady rigid-body rotation with no heat transfer, the equilibrium is **isentropic** ($p/\rho^\gamma$ uniform). Combining $\nabla p = \rho\omega^2 \boldsymbol{r}_\perp$ with the isentropic relation gives

$$\frac{p(r_\perp)}{p_{\rm ref}} = \left[1 + \frac{\gamma-1}{2}\,\frac{\omega^2\bigl(r_\perp^2 - r_{\rm ref}^2\bigr)}{\gamma R\, T_{\rm ref}}\right]^{\gamma/(\gamma-1)}, \qquad \frac{T(r_\perp)}{T_{\rm ref}} = \left(\frac{p(r_\perp)}{p_{\rm ref}}\right)^{(\gamma-1)/\gamma}.$$

The reference state $(p_{\rm ref}, T_{\rm ref})$ at $r_{\rm ref}$ is fixed by total-mass conservation against the initial uniform field.

**Expected behaviour.** The centrifugal source builds a radially stratified field with $p$ increasing monotonically with $r_\perp$ at constant entropy, while $\boldsymbol{w}$ relaxes back to zero. The case verifies that

- the centrifugal contribution is wired into both the momentum and energy residuals with the correct sign and magnitude,
- the rotational-periodic BC closes the ring with the proper velocity rotation
- the converged pressure and temperature profiles match the isentropic closed form above (radial $L_2$ error against the analytical reference is the verification metric). The initial solution in **ic.tec** stores pressure profile that should balance the centrifugal contribution, i.e. the final solution should be similar to the initial one. (Residuals start very low)

**Sensitivity note.** At $\omega = 100$ rad/s the rim Mach number is $\omega r_{\rm out}/\sqrt{\gamma R T} \approx 0.04$, giving a pressure variation of only $\sim 0.1\%$ across the annulus. Verification at this setting requires tight residual convergence and a precise mass-conservation pinning of $p_{\rm ref}$; running at higher $\omega$ (e.g. 1000 rad/s, $\sim 10\%$ variation) gives the comparison more dynamic range.
