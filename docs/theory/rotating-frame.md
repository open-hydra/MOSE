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

The energy equation (with specific total energy $E = e + \frac{1}{2} |\boldsymbol{w}|^2$, $e$ the specific internal energy of the gas) receives a source term:

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
| SSG–LRR (RSM; in the source code, `turbulence` does not accept it) | $W_{ij}^{\rm abs}$ in the $C_5$ pressure-strain term |

In SA-RC and SST-RC the rate of change of the strain in $\tilde r$ is its steady convective part $\mathbf{u}\cdot\nabla S_{ij}$: the frame enters through the absolute vorticity, but the frame-rotation term $(\epsilon_{imn}S_{jn} + \epsilon_{jmn}S_{in})\,\Omega_m$ of the published $\tilde r$ is not included.

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

A wall at rest in the laboratory (a casing or a hub that does not turn) is listed in `[MOSE-Rotating-Frame]` as `stationary-face-<n> = <block> <face>`, numbered from 1; the walls `301` and `302` of that block face take it. The other records of that face may be of the types that never read whether their face is stationary, whose boundary condition is then the same as without the key: block connections, periodic and chimera faces (`101`, `102`, `200`, `201`), extrapolation (`400`), the ambient boundary (`407`), the mapped boundary (`410`) and unused faces (`0`); so a casing face that is a wall in part and a block connection in part can be stationary. A stationary face that covers a type that would take the velocity of the wall at rest in part only (the inviscid flux of a symmetry `300`, of the inlets and outlets `401`–`406`, `408` and `420`, and of the walls `103`, `104` and `502`–`505`), that has no wall `301` or `302`, or that names a face with no record in `bc.txt` stops the run at set-up. Such a wall moves in the rotating frame with

$$\boldsymbol{w}_{\text{wall}} = -\boldsymbol{\Omega} \times (\boldsymbol{x} - \boldsymbol{x}_o)$$

($\boldsymbol{x}_o$ the `origin` of the frame), a velocity linear in $\boldsymbol{x}$: its no-slip condition and its ghost cells take this velocity, and its wall gradient also the differences of this velocity along the wall, from its gradient $-[\boldsymbol{\Omega}\times]$ (see [Gradient Computation for Diffusive Fluxes](numerics.md#gradient-computation-for-diffusive-fluxes)). The energy flux through such a wall is the heat flux plus the work of the wall stress on the wall velocity, $(\boldsymbol{\tau}\cdot\mathbf{n})\cdot\boldsymbol{w}_{\text{wall}}$, as a face between cells carries $(\boldsymbol{\tau}\cdot\mathbf{n})\cdot\boldsymbol{w}$ ([Viscous Fluxes](#viscous-fluxes)); on a wall that turns with the frame $\boldsymbol{w}_{\text{wall}} = \boldsymbol{0}$ and the energy flux is the heat flux alone.

A stationary wall must be a surface of revolution about the axis of the frame, as a casing or a hub is: its velocity in the frame is then tangent to it, $\boldsymbol{w}_{\text{wall}}\cdot\mathbf{n} = 0$, so that no mass crosses it and the pressure does no work on it. A wall that is not, such as a strut at rest in the laboratory, moves across the mesh of the rotating frame and cannot be one of its faces.

## Physical Interpretation

The rotating-frame formulation is equivalent to solving the full inertial-frame problem but with fictitious body forces representing:

1. **Coriolis acceleration:** Deflects moving fluid perpendicular to both the rotation axis and the velocity vector.
2. **Centrifugal acceleration:** Acts radially outward from the rotation axis.

For steady-state turbomachinery flows, a steady rotating-frame solution avoids the need to track time-periodic structures in the inertial frame, enabling efficient steady-state analysis.

## Implementation Notes

- The computational domain is assumed to rotate as a **rigid body** at constant $\boldsymbol{\Omega}$.
- The axis and angular velocity magnitude/direction must be specified at initialization. When `input.ini` is read again
  during the run (every `ini-diter` iterations), the `omega`, `axis` and `origin` read again replace those of the set-up
  in every term at once (the log of the run prints them): the Coriolis and centrifugal sources, the velocity of the
  stationary walls, of the inlets and of the far field, and the power of the torque file. The frame is taken at constant
  $\boldsymbol{\Omega}$ between two readings: the force $-\rho\,\dot{\boldsymbol{\Omega}} \times \boldsymbol{r}$ of a
  frame that speeds up is not modelled. A change of `rotational-frame` or of the stationary faces, which changes the
  equations or the boundary conditions, is an error of the reading: a `stationary-face-<n>` key that names another face, or one that
  the set-up does not have. A key that the file read again leaves out keeps its face, as every parameter left out keeps
  its value. The keys are counted from 1 up to the first one missing, at the start and at every reading: a key after a
  gap of the numbering is not read. A value of `axis` or `origin` read again that is not finite (NaN or infinity) or does not
  read as numbers, a `stationary-face-<n>` value that does not read, or an axis of zero length, is an error of the reading
  too. With `ini-error = warn`, the default, a reading with an error is not taken and the frame goes on as it was; with
  `ini-error = stop` the run stops ([Input](../user/input.md)).
- The axis of rotation passes through `origin` of `[MOSE-Rotating-Frame]` (the origin of the coordinates by default), and $\boldsymbol{r}$ is measured from it.
- All flux routines (convective and viscous) operate on the relative velocity without modification.

The rotating-frame verification cases are documented in [the V&V section](../vv/rotating-frame.md).
