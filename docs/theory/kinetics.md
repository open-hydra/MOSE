# Finite-Rate Kinetics

## Overview

MOSE computes species mass source terms for chemically reacting flows using finite-rate kinetics, accounting for elementary, three-body, and pressure-dependent reactions. The kinetics framework is provided by the [FLINT](https://github.com/MarcoGrossi92/FLINT) library and supports both general (runtime-loaded) and mechanism-specific (hard-coded) evaluation routines.

---

## Mathematical Formulation

The rate of change of molar concentration $c_i$ for species $i$ is governed by contributions from all $N_r$ reactions:

$$
\frac{dc_i}{dt} = \sum_{r=1}^{N_r} \nu_{i,r} \, R_r
$$

where:

- $c_i$ = molar concentration of species $i$ (mol/m³)
- $\nu_{i,r} = \nu'_{i,r} - \nu''_{i,r}$ = net stoichiometric coefficient of species $i$ in reaction $r$
- $\nu''_{i,r}$ = stoichiometric coefficient of species $i$ as a reactant
- $\nu'_{i,r}$ = stoichiometric coefficient of species $i$ as a product
- $R_r$ = rate of progress for reaction $r$ (mol/(m³·s))

For a general elementary reaction:

$$
\sum_{i=1}^{N_s} \nu''_{i,r} \, X_i \rightleftharpoons \sum_{i=1}^{N_s} \nu'_{i,r} \, X_i
$$

The mass production rate $\dot{\omega}_i$ (kg/(m³·s)) used in the species conservation equations is related to the molar rate by:

$$
\dot{\omega}_i = M_i \frac{dc_i}{dt}
$$

where $M_i$ is the molecular weight of species $i$.

---

## Reaction Rate (Arrhenius)

The rate of progress for reaction $r$ follows the modified Arrhenius form [1,2]:

$$
R_r = k_{f,r}(T) \prod_{j=1}^{N_s} c_j^{\nu''_{j,r}} - k_{b,r}(T) \prod_{j=1}^{N_s} c_j^{\nu'_{j,r}}
$$

where the forward rate coefficient is:

$$
k_{f,r}(T) = A_r \, T^{b_r} \exp\!\left(-\frac{E_{a,r}}{R_u T}\right)
$$

**Parameters:**

| Symbol | Description |
|--------|-------------|
| $A_r$ | Pre-exponential factor (units depend on reaction order) |
| $b_r$ | Temperature exponent (dimensionless) |
| $E_{a,r}$ | Activation energy (J/mol) |
| $R_u$ | Universal gas constant, 8.314 J/(mol·K) |

The backward rate coefficient is computed from thermodynamic equilibrium:

$$
k_{b,r}(T) = \frac{k_{f,r}(T)}{K_{c,r}(T)}
$$

where $K_{c,r}(T)$ is the equilibrium constant in concentration units [3].

---

## Pressure-Dependent Reactions

### Three-Body Reactions

Three-body (chaperon) reactions involve a collision partner $M$ that provides or absorbs energy without being consumed [4]:

$$
A + B + M \rightleftharpoons AB + M
$$

The effective third-body concentration is:

$$
[M] = \sum_{j=1}^{N_s} \alpha_{j,r} \, c_j
$$

where $\alpha_{j,r}$ is the third-body efficiency of species $j$ in reaction $r$ (default $\alpha_{j,r} = 1$ unless specified otherwise).

### Lindemann Falloff

Lindemann reactions describe the transition between low-pressure (termolecular) and high-pressure (bimolecular) limiting kinetics [5,6]:

$$
k_r(T, P) = k_{\infty,r}(T) \left(\frac{P_r}{1 + P_r}\right) F_r
$$

where:

- $k_{\infty,r}(T)$ = high-pressure limit rate coefficient
- $k_{0,r}(T)$ = low-pressure limit rate coefficient
- $P_r$ = reduced pressure:

$$
P_r = \frac{k_{0,r}(T) \, [M]}{k_{\infty,r}(T)}
$$

- $F_r$ = broadening factor (see Troe formulation below; $F_r = 1$ for simple Lindemann)

**Limiting behavior:**

| Regime | Condition | Effective rate |
|--------|-----------|----------------|
| Low pressure | $P_r \ll 1$ | $k_r \approx k_{0,r}[M]$ (third-order) |
| High pressure | $P_r \gg 1$ | $k_r \approx k_{\infty,r}$ (second-order) |

### Troe Formulation

The Troe formulation provides an empirical broadening factor $F_r$ that better matches experimental falloff curves [7,8]:

$$
\log_{10} F_r = \frac{\log_{10} F_{\text{cent},r}}{1 + \left[\dfrac{\log_{10} P_r + c}{n - d\,(\log_{10} P_r + c)}\right]^2}
$$

where the centering factor is:

$$
F_{\text{cent},r} = (1-a) \exp\!\left(-\frac{T}{T^{***}}\right) + a \exp\!\left(-\frac{T}{T^*}\right) + \exp\!\left(-\frac{T^{**}}{T}\right)
$$

and the auxiliary constants are:

$$
c = -0.4 - 0.67 \log_{10} F_{\text{cent},r}, \qquad n = 0.75 - 1.27 \log_{10} F_{\text{cent},r}, \qquad d = 0.14
$$

The four Troe parameters ($a$, $T^*$, $T^{**}$, $T^{***}$) are fitted to experimental data for each reaction.

!!! note "Simplified Troe form"
    When only three parameters are available ($T^{**} \to \infty$), the centering factor reduces to:
    $F_{\text{cent},r} = (1-a) \exp(-T/T^{***}) + a \exp(-T/T^*)$

---

## ODE Integration

The species production rates form a stiff system of ordinary differential equations that must be integrated at each grid cell during operator-split time advancement. MOSE uses the stiff ODE solvers provided by the [OSlo](https://github.com/MarcoGrossi92/OSlo) library:

| Solver | Type | Description |
|--------|------|-------------|
| `H-radau5` | Implicit Runge-Kutta | Radau IIA method, 5th order (default) |
| `sdirk4b` | Singly-DIRK | 4th-order L-stable SDIRK |
| `ros4` | Rosenbrock | 4th-order Rosenbrock method |

---

## References

[1] Arrhenius, S. "Über die Reaktionsgeschwindigkeit bei der Inversion von Rohrzucker durch Säuren." *Zeitschrift für Physikalische Chemie*, vol. 4, 1889, pp. 226–248.

[2] Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*, 2nd edition. John Wiley & Sons, 2003.

[3] Smith, G. P., et al. "GRI-Mech 3.0." http://www.me.berkeley.edu/gri_mech/

[4] Baulch, D. L., et al. "Evaluated Kinetic Data for Combustion Modeling: Supplement II." *Journal of Physical and Chemical Reference Data*, vol. 34, no. 3, 2005, pp. 757–1397.

[5] Lindemann, F. A. "Discussion on 'The Radiation Theory of Chemical Action'." *Transactions of the Faraday Society*, vol. 17, 1922, pp. 598–599.

[6] Gilbert, R. G., Luther, K., and Troe, J. "Theory of Thermal Unimolecular Reactions in the Fall-off Range. II." *Berichte der Bunsengesellschaft für Physikalische Chemie*, vol. 87, 1983, pp. 169–177.

[7] Troe, J. "Predictive Possibilities of Unimolecular Rate Theory." *Journal of Physical Chemistry*, vol. 83, no. 1, 1979, pp. 114–126.

[8] Gilbert, R. G., Smith, S. C., and Jordan, M. J. T. *UNIMOL: Calculation of Fall-off Curves for Unimolecular and Recombination Reactions*. Blackwell Scientific Publications, 1990.

---
