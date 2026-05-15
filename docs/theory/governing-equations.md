# Governing Equations

MOSE solves the compressible Navier–Stokes equations for multi-species ideal-gas
mixtures in conservative form on structured multi-block grids.
Both the inviscid (Euler) and the viscous (Navier–Stokes) limits are supported.

The governing system is written in conservative form as

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla\!\cdot\!\mathbf{F}^{c}(\mathbf{U})
= \nabla\!\cdot\!\mathbf{F}^{v}(\mathbf{U})
+ \mathbf{S}
$$

where $\mathbf{U}$ is the vector of conservative variables,
$\mathbf{F}^{c}$ the convective (inviscid) flux,
$\mathbf{F}^{v}$ the viscous (diffusive) flux,
and $\mathbf{S}$ a source term vector (chemistry, turbulence, soot).

---

## State Vectors

The flow state is described by two complementary representations: **primitive variables**
used for physical interpretation and I/O, and **conservative variables** for the temporal
integration and flux evaluation in the solver.

### Primitive variables

The primitive state vector is

$$
\mathbf{P} = [\rho_1, \rho_2, \ldots, \rho_{N_s}, u, v, w, p, f_\text{soot}, f_\text{pass}, \mathbf{q}_\text{rans}]
$$

where the components are:

| Component | Symbol | Description |
|:-----------:|:------:|----------|
| Species densities | $\rho_s$ ($s=1,\ldots,N_s$) | Partial mass density of each species |
| Velocity | $u, v, w$ | Cartesian velocity components |
| Pressure | $p$ | Thermodynamic pressure |
| Soot variable | $f_\text{soot}$ | Soot volume fraction or total soot density (optional) |
| Passive scalar | $f_\text{pass}$ | Inert tracer field for mixing studies (optional) |
| RANS variables | $\mathbf{q}_\text{rans}$ | Turbulence model variables (optional) |

**Soot Model Variables (if active):**
  
- $f_\text{soot} = \rho_\text{soot} / \rho$: volume fraction of soot particles  
- Additional soot moments or precursor scalars may be tracked

**RANS Model Variables (if active):**
  
- **Spalart–Allmaras**: $\tilde{\nu}$ (modified eddy viscosity)  
- **$k$–$\omega$ (SST)**: $k$ (turbulent kinetic energy), $\omega$ (specific dissipation rate)  

### Conservative variables

Conservation laws require the flux-divergence form, which motivates the conservative
variable vector:

$$
\mathbf{U} = [\rho_1, \rho_2, \ldots, \rho_{N_s}, \rho\,u, \rho\,v, \rho\,w, \rho\,E_0, f_\text{soot}, f_\text{pass}, \mathbf{q}_\text{rans}]
$$

where the first $N_s + 4$ components correspond to mass-momentum-energy conservation:

$$
\rho = \sum_{s=1}^{N_s} \rho_s \quad \text{(total density)}
$$

$$
E_0 = e + \frac{1}{2}(u^2 + v^2 + w^2) \quad \text{(total specific energy)}
$$

$$
e = h - \frac{p}{\rho} \quad \text{(internal specific energy)}
$$

The mixture specific enthalpy is computed as a mass-fraction-weighted sum:

$$
h = \sum_{s=1}^{N_s} Y_s\, h_s(T), \quad Y_s = \frac{\rho_s}{\rho}
$$

where species enthalpies $h_s(T)$ are obtained from NASA polynomial correlations
(see [Thermodynamic and Transport Properties](thermo.md)).

**Passive and Turbulence Variables:**

Soot, passive scalars, and RANS variables do not enter the energy equation
or equation of state. They are advected with the bulk flow and obey scalar transport equations with their own source terms.

### Transformation between conservative and primitive variables

Forward transformation ($\mathbf{P} \to \mathbf{U}$): Given primitive variables, the transformation is straightforward: species densities remain unchanged, and the conservative energy is computed by combining internal energy with kinetic energy.

Inverse transformation ($\mathbf{U} \to \mathbf{P}$): Recovering primitive variables from conservative ones requires solving for temperature implicitly, since the ideal-gas law couples pressure, density, and temperature through the equation of state. The algorithm employs Newton–Raphson iteration:

1. Extract densities and momentum:  
   $\rho_s = U_s$ for each species; $\rho = \sum_s \rho_s$  
   
2. Recover velocity:  
   $(u, v, w) = (U_u, U_v, U_w) / \rho$
   
3. Compute specific internal energy:  
   $e = \frac{U_E}{\rho} - \frac{1}{2}(u^2 + v^2 + w^2)$
   
4. Solve for temperature iteratively:  
   Newton–Raphson iterations with the implicit relation:
   $$T^{(k+1)} = T^{(k)} + \frac{e - e(T^{(k)})}{c_v(T^{(k)})}$$
   where $e(T)$ and $c_v(T)$ are from NASA polynomial fits.  
   Iterate until relative error drops below $10^{-6}$ or a specified tolerance.
   
5. Compute pressure:  
   From the ideal-gas equation of state,
   $$p = \rho\, R_\text{mix}\, T$$
   where $R_\text{mix} = \sum_s Y_s\, R_s$ is the mixture gas constant (mass-fraction-weighted).

---

## Convective Fluxes

The convective flux through a face of outward unit normal $\hat{\mathbf{n}}$ is

$$
\mathbf{F}^{c}\!\cdot\!\hat{\mathbf{n}} = [
  \rho_1\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}), 
  \ldots, 
  \rho_{N_s}\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}),
  \rho\,u\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_x,
  \rho\,v\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_y,
  \rho\,w\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}}) + p\,\hat{n}_z,
  \rho\,H_0\,(\mathbf{v}\!\cdot\!\hat{\mathbf{n}})
]
$$

where $H_0 = h + \tfrac{1}{2}\lvert\mathbf{v}\rvert^2$ is the total
specific enthalpy.  The normal mass flux is: $\dot{m} = \rho\,(\mathbf{v}\,\cdot\,\hat{\mathbf{n}})$

---

## Viscous Fluxes

### Stress tensor

The viscous stress tensor follows the Newtonian fluid model with the
Stokes hypothesis $\lambda = -\tfrac{2}{3}\mu$:

$$
\tau_{ij} =
\mu\!\left(\frac{\partial v_i}{\partial x_j}
+ \frac{\partial v_j}{\partial x_i}\right)
- \frac{2}{3}\,\mu\,(\nabla\!\cdot\!\mathbf{v})\,\delta_{ij}
$$

where $\mu = \mu_\ell + \mu_t$ is the sum of laminar and eddy viscosity.

For two-equation RANS models the isotropic part of the Reynolds stress
is included:

$$
\tau_{ij}^{\text{2eq}} = \tau_{ij} - \tfrac{2}{3}\,\rho\,k\,\delta_{ij}
$$

### Species diffusion

Each species has a diffusive mass flux

$$
\mathbf{j}_s = -\rho\,D_{m,s}\,\nabla Y_s
$$

where the mixture-averaged diffusion coefficient is

$$
D_{m,s} = \frac{\mu_\ell}{\rho\,\text{Sc}}
         + \frac{\mu_t}{\rho\,\text{Sc}_t}
$$

A mass-flux correction enforces $\sum_s \mathbf{j}_s = 0$.

### Heat conduction

The effective thermal conductivity is

$$
\kappa = k_\ell + \frac{\mu_t\,c_p}{\text{Pr}_t}
$$

### Diffusive flux vector

The complete viscous flux projected onto $\hat{\mathbf{n}}$ is

$$
\mathbf{F}^{v}\!\cdot\!\hat{\mathbf{n}} = [
  -\mathbf{j}_1 \cdot \hat{\mathbf{n}},
  \ldots,
  -\mathbf{j}_{N_s} \cdot \hat{\mathbf{n}},
  \boldsymbol{\tau}\!\cdot\!\hat{\mathbf{n}},
  (\boldsymbol{\tau}\!\cdot\!\mathbf{v})\!\cdot\!\hat{\mathbf{n}} + \kappa\,\nabla T\!\cdot\!\hat{\mathbf{n}} + \rho\!\sum_s D_{m,s}\,(\nabla Y_s\!\cdot\!\hat{\mathbf{n}})\,h_s
]
$$

The energy flux therefore comprises viscous work, heat conduction,
and species enthalpy diffusion.

---

## Source Terms

Beyond convection and diffusion, the flow may be driven by source terms in chemistry,
turbulence, and soot formation. These are summarized below:

| Physical Process | Affected Variables | Integration Strategy |
|--------|-------------------|----------|
| Chemical kinetics | Species partial densities $\rho_s$ | Strang operator splitting; finite-rate or quasi-steady-state chemistry |
| Turbulence | RANS variables $k, \omega, \tilde{\nu}$, etc. | Explicit source in spatial residual; wall-bounded corrections |
| Soot transport | Soot volume fraction $f_\text{soot}$ | Explicit source; includes nucleation, surface reactions, coagulation |

---
