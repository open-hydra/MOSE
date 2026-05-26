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
Inflow/outflow and other standard boundary conditions are applied directly to the relative velocity field $\boldsymbol{w}$ without modification. When defining the inflow properties it is necessary to specify them in the **laboratory frame of reference**. The correction in the rotating frame is then performed automatically.

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

The rotating-frame verification cases are documented in [the V&V section](../vv/rotating-frame.md).
