# Thermodynamic and Transport Properties

## Overview

MOSE models a mixture of $N_s$ thermally perfect gases following the thermal equation of state:

$$
p = \rho R_\text{mix} T
$$

where $p$ is pressure, $\rho$ is density, $R_\text{mix}$ is the mixture gas constant, and $T$ is temperature.

Thermodynamic and transport properties are provided by the [FLINT](https://github.com/MarcoGrossi92/FLINT) library. Two evaluation methods are available:

- **Native tabulated data** (default) — temperature-varying tabulated properties for each species: specific heat capacity $c_{p}$, enthalpy $h$, entropy $s$, dynamic viscosity $\mu$, and thermal conductivity $k$, plus an optional table of binary diffusion coefficients $\mathcal{D}_{sj}$ for the [multicomponent diffusion](#species-diffusion) closure. This is the recommended path for production runs.

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

The transport properties consumed by the [viscous
fluxes](governing-equations.md#viscous-fluxes) are the mixture viscosity
$\mu_\ell$, the mixture thermal conductivity $k_\ell$, and the per-species mass
diffusion coefficients $D_s$. The first two are mixture averages of the tabulated
single-species data; the third has two selectable closures (see below).

### Viscosity and thermal conductivity

Mixture viscosity and conductivity are computed with **Wilke's mixing rule** [3,4],
which accounts for molecular interactions between different species.

**Wilke's interaction parameter:**

$$
\phi_{ij} = \frac{\left[1 + \left(\mu_i/\mu_j\right)^{1/2} \left(M_j/M_i\right)^{1/4}\right]^2}{\sqrt{8 \left(1 + M_i/M_j\right)}}
$$

where $\mu_i$ is the dynamic viscosity of species $i$ and $M_i$ is its molecular weight.

**Mixture dynamic viscosity:**

$$
\mu_\ell = \sum_{i=1}^{N_s} \frac{X_i \mu_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}}
$$

**Mixture thermal conductivity:**

$$
k_\text{Wilke} = \sum_{i=1}^{N_s} \frac{X_i k_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}}
$$

**Mole–mass fraction conversion:**

$$
X_s = \frac{Y_s / M_s}{\sum_{j=1}^{N_s} Y_j / M_j}
$$

The laminar conductivity actually used in the energy flux is selected by the
laminar Prandtl number $Pr_\ell$:

$$
k_\ell =
\begin{cases}
k_\text{Wilke} & Pr_\ell \le 0 \quad (\text{computed mixture conductivity}) \\[4pt]
\dfrac{\mu_\ell\, c_{p,\text{mix}}}{Pr_\ell} & Pr_\ell > 0 \quad (\text{constant-Prandtl override})
\end{cases}
$$

The override $k_\ell = \mu_\ell c_{p}/Pr_\ell$ is the path used, for example, to
impose a unity Lewis number ($Pr_\ell = Sc = 1$).

### Species diffusion

Each species carries a diffusive mass flux $\mathbf{j}_s = -\rho\,D_s\,\nabla Y_s$
(closed by a mass-conservation correction, see [governing
equations](governing-equations.md#species-diffusion)). The laminar diffusion
coefficient $D_s$ has **two closures**, selected by the sign of the laminar Schmidt
number $Sc$:

=== "Constant Schmidt ($Sc > 0$)"

    A single diffusivity, tied to the mixture viscosity, is shared by all species:

    $$
    D_s^{\ell} = \frac{\mu_\ell}{\rho\,Sc}
    $$

    Robust and inexpensive; appropriate when differential diffusion between species
    is unimportant. Combined with the $Pr_\ell$ override above it reproduces a
    constant-Lewis closure ($Le_s = Sc/Pr_\ell$).

=== "Mixture-averaged multicomponent ($Sc \le 0$)"

    Each species gets its own composition- and temperature-dependent coefficient
    from the **Curtiss–Hirschfelder mixture-averaged rule** [2,6]:

    $$
    D_s^{\ell} = \frac{1 - X_s}{\displaystyle\sum_{j \ne s} X_j / \mathcal{D}_{sj}}
    $$

    where $\mathcal{D}_{sj}(T,p)$ are the binary diffusion coefficients. As a cell
    approaches a pure species ($X_s \to 1$) the ratio is indeterminate; MOSE then
    falls back to the mean binary diffusivity of species $s$, which is immaterial
    because $\nabla Y_s \to 0$ there and the flux is fixed by the mass-conservation
    correction.

In both cases the turbulent contribution is added on top with a constant turbulent
Schmidt number,

$$
D_s = D_s^{\ell} + \frac{\mu_t}{\rho\,Sc_t},
$$

with $\mu_t$ the [eddy viscosity](turbulence.md) (zero for laminar runs).

#### Binary diffusion coefficients

The binary coefficients $\mathcal{D}_{sj}$ are a kinetic-theory **pair** property:
composition-independent and, for ideal gases, functions of temperature and pressure
only, scaling exactly as $\mathcal{D}\propto 1/p$. They are supplied as an offline
table over the temperature grid at a reference pressure $p_\text{ref}$ (one column
per unique unordered pair $s\!<\!j$, since $\mathcal{D}_{sj}=\mathcal{D}_{js}$) and
rescaled to the local pressure at run time:

$$
\mathcal{D}_{sj}(T,p) = \frac{p_\text{ref}}{p}\,\mathcal{D}_{sj}(T,p_\text{ref}).
$$

The reference pressure is embedded in the table header, so the same table is valid
at atmospheric and elevated pressures. The table (`diffusion.dat`) is generated by
the preprocessor alongside the thermodynamic and transport tables; the
multicomponent closure requires it to be present.

!!! note "Mixture rule: $(1-X_s)$ vs $(1-Y_s)$"
    MOSE uses the mole-fraction numerator $(1 - X_s)$ — the original
    Curtiss–Hirschfelder form. The Chemkin/Kee convention (and Cantera's
    `mix_diff_coeffs`) instead use the mass-fraction numerator $(1 - Y_s)$, with the
    same denominator and binary data. The two differ by $\sim\!10\text{–}20\%$ per
    species in mixtures of light and heavy molecules. The $(1 - Y_s)$ variant is
    kept as commented reference code in the diffusion routine for users who prefer
    the Chemkin/Cantera convention.

---

## References

[1] Poinsot, T., and Veynante, D. *Theoretical and Numerical Combustion*, 3rd edition. Published by the authors, 2012.

[2] Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*, 2nd edition. John Wiley & Sons, 2003.

[3] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." *The Journal of Chemical Physics*, vol. 18, no. 4, 1950, pp. 517–519.

[4] Bird, R. B., Stewart, W. E., and Lightfoot, E. N. *Transport Phenomena*, 2nd edition. John Wiley & Sons, 2002.

[5] Blazek, J. *Computational Fluid Dynamics: Principles and Applications*, 3rd edition. Butterworth-Heinemann, 2015.

[6] Hirschfelder, J. O., Curtiss, C. F., and Bird, R. B. *Molecular Theory of Gases and Liquids*. John Wiley & Sons, 1954.

---
