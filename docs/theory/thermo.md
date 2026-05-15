# Thermodynamic and Transport Properties

## Overview

MOSE models a mixture of $N_s$ thermally perfect gases following the thermal equation of state:

$$
p = \rho R_\text{mix} T
$$

where $p$ is pressure, $\rho$ is density, $R_\text{mix}$ is the mixture gas constant, and $T$ is temperature.

Thermodynamic and transport properties are provided by the [FLINT](https://github.com/MarcoGrossi92/FLINT) library. Two evaluation methods are available:

- **Native tabulated data** (default) — temperature-varying tabulated properties for each species: specific heat capacity $c_{p}$, enthalpy $h$, entropy $s$, dynamic viscosity $\mu$, and thermal conductivity $k$. This is the recommended path for production runs.

- **Cantera integration** (optional) — delegates property evaluation to the [Cantera](https://cantera.org) library. Useful for validation and benchmarking.

---

## Thermodynamic Properties

For tabulated properties, linear interpolation is applied to retrieve any property for each species, then mixture thermodynamic quantities are computed using mass-weighted averaging [1,2].

**Mixture density:**

$$
\rho = \sum_{s=1}^{N_s}\rho_s
$$

**Mixture gas constant:**

$$
R_\text{mix} = \sum_{s=1}^{N_s} Y_s R_s = \sum_{s=1}^{N_s} \frac{\rho_s}{\rho} R_s
$$

where $Y_s = \rho_s/\rho$ is the mass fraction of species $s$, and $R_s = R_u/M_s$ is the specific gas constant ($R_u = 8314.46$ J/(kmol·K) is the universal gas constant and $M_s$ is the molecular weight).

**Specific heat capacities:**

$$
c_{p,\text{mix}} = \sum_{s=1}^{N_s} Y_s c_{p,s}, \qquad c_{v,\text{mix}} = c_{p,\text{mix}} - R_\text{mix}
$$

**Heat capacity ratio:**

$$
\gamma = \frac{c_{p,\text{mix}}}{c_{v,\text{mix}}} = \frac{c_{p,\text{mix}}}{c_{p,\text{mix}} - R_\text{mix}}
$$

**Speed of sound:**

$$
a = \sqrt{\gamma R_\text{mix} T}
$$

**Specific enthalpy:**

$$
h(\rho_s, T) = \sum_{s=1}^{N_s} Y_s h_s(T)
$$

**Total specific enthalpy (including kinetic energy):**

$$
h_0(\rho_s, T, \mathbf{u}) = h(\rho_s, T) + \frac{1}{2} |\mathbf{u}|^2
$$

**Specific internal energy:**

$$
e(\rho_s, T) = h(\rho_s, T) - \frac{p}{\rho} = \sum_{s=1}^{N_s} Y_s h_s(T) - R_\text{mix} T
$$

**Total specific internal energy:**

$$
e_0(\rho_s, T, \mathbf{u}) = e(\rho_s, T) + \frac{1}{2}|\mathbf{u}|^2
$$

---

## Transport Properties

Mixture transport properties are computed using **Wilke's mixing rule** [3,4], which accounts for molecular interactions between different species.

**Wilke's interaction parameter:**

$$
\phi_{ij} = \frac{\left[1 + \left(\mu_i/\mu_j\right)^{1/2} \left(M_j/M_i\right)^{1/4}\right]^2}{\sqrt{8 \left(1 + M_i/M_j\right)}}
$$

where $\mu_i$ is the dynamic viscosity of species $i$ and $M_i$ is its molecular weight.

**Mixture dynamic viscosity:**

$$
\mu_\text{mix} = \sum_{i=1}^{N_s} \frac{X_i \mu_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}}
$$

**Mixture thermal conductivity:**

$$
k_\text{mix} = \sum_{i=1}^{N_s} \frac{X_i k_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}}
$$

**Mole–mass fraction conversion:**

$$
X_s = \frac{Y_s / M_s}{\sum_{j=1}^{N_s} Y_j / M_j}
$$

---

## References

[1] Poinsot, T., and Veynante, D. *Theoretical and Numerical Combustion*, 3rd edition. Published by the authors, 2012.

[2] Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*, 2nd edition. John Wiley & Sons, 2003.

[3] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." *The Journal of Chemical Physics*, vol. 18, no. 4, 1950, pp. 517–519.

[4] Bird, R. B., Stewart, W. E., and Lightfoot, E. N. *Transport Phenomena*, 2nd edition. John Wiley & Sons, 2002.

[5] Blazek, J. *Computational Fluid Dynamics: Principles and Applications*, 3rd edition. Butterworth-Heinemann, 2015.

---
