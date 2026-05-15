# Turbulence Modelling

MOSE provides Reynolds-Averaged Navier–Stokes (RANS) closures for
turbulent flows.  Three families of eddy-viscosity models are
available, ranging from a one-equation model to two-equation
formulations.  All models are selected at run time via the input file.

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
\Bigr]}_{\text{Diffusion}}
- \underbrace{c_{w1}\,f_w\,\rho\!\left(\frac{\tilde\nu}{y}\right)^{\!2}}_{\text{Destruction}}
$$

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

where $\Omega = \sqrt{2\,W_{ij}\,W_{ij}}$ is the vorticity magnitude and
$y$ is the distance to the nearest wall.

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
\tilde{r} = \frac{C_D}{D_A^2}
$$

where $C_D$ involves the material derivative interaction between the
strain and vorticity tensors:

$$
C_D = W_{ij}\,S_{jk}\!\left(\frac{DS_{ki}}{Dt} - \frac{DW_{ki}}{Dt}\right)
$$

The correction factor:

$$
f_{r1} = (1 + c_{r1})\,\frac{2\,r^\ast}{1 + r^\ast}\,
  \bigl(1 - c_{r3}\arctan(c_{r2}\,\tilde{r})\bigr) - c_{r1}
$$

bounded to $[0,\;1.25]$, with constants $c_{r1} = 1$, $c_{r2} = 2$,
$c_{r3} = 1$.  Production is multiplied by $\max(1,\,f_{r1})$.

---

## Menter SST $k$–$\omega$ — Two-Equation Model

The Shear Stress Transport model blends a $k$–$\omega$ formulation
(near walls) with a $k$–$\varepsilon$-like behaviour (in the
freestream) using blending functions.

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
\omega_\text{wall} = \frac{6\,\nu}{0.075\,y^2}
$$

### Energy coupling *(optional)*

When enabled, the turbulence kinetic energy production/destruction
contributes to the mean-flow energy equation:

$$
\frac{\partial(\rho e)}{\partial t} \mathrel{+}= S_k
$$

---

## Wilcox 2006 $k$–$\omega$ — Two-Equation Model

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

The $\omega$-destruction coefficient is modified by the Baerten
correction:

$$
\beta = \beta_0\,f_\beta, \qquad
X_\omega = \frac{|W_{ij}\,W_{jk}\,S_{ki}|}{(\beta^\ast\omega)^3},
\qquad
f_\beta = \frac{1 + 85\,X_\omega}{1 + 100\,X_\omega}
$$

### Cross-diffusion *(conditional)*

Cross-diffusion is included only when $\nabla k\!\cdot\!\nabla\omega \le 0$:

$$
\text{CD} = \frac{\sigma_d\,\rho}{\omega}\,
  (\nabla k\!\cdot\!\nabla\omega)
\qquad \text{if } \nabla k\!\cdot\!\nabla\omega \le 0
$$

### Production limiter

More permissive than SST:

$$
P_k = \min\!\bigl(\mu_t\,S^2,\;20\,\beta^\ast\,\rho\,\omega\,k\bigr)
$$

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
| `RANS_Extrapolate_Wall` | Ghost-cell extrapolation for RANS variables |

This design allows switching models at run time without recompilation.

### Blowing / wall-suction correction

An optional correction (`blowing_corr`) modifies the wall boundary-layer
treatment for SST and Wilcox 2006 models in the presence of wall
injection or suction.

---

## References

1. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for
   aerodynamic flows," AIAA-92-0439, 1992.
2. F. R. Menter, "Two-equation eddy-viscosity turbulence models for
   engineering applications," *AIAA J.*, 32(8), 1994.
3. F. R. Menter, M. Kuntz, R. Langtry, "Ten years of industrial
   experience with the SST turbulence model," in *Turbulence, Heat
   and Mass Transfer 4*, Begell House, 2003.
4. D. C. Wilcox, *Turbulence Modeling for CFD*, 3rd ed., DCW Industries,
   2006.
5. P. R. Spalart, M. L. Shur, "On the sensitization of turbulence models
   to rotation and curvature," *Aerosp. Sci. Technol.*, 1(5), 1997.
