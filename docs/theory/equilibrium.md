# Chemical Equilibrium

## Overview

MOSE can compute chemical equilibrium compositions for multi-component gas mixtures as an alternative to finite-rate kinetics. The equilibrium solver handles the UV problem: constant internal energy $U$ and volume $V$ (equivalently, constant specific internal energy and density).

The implementation follows the NASA CEA (Chemical Equilibrium with Applications) methodology developed by Gordon and McBride [1,2], which is the standard approach for aerospace and combustion equilibrium calculations. The solver is provided by the [FLINT](https://github.com/MarcoGrossi92/FLINT) library.

---

## Mathematical Formulation

For a closed system at constant internal energy $U$ and volume $V$, chemical equilibrium corresponds to the state of maximum entropy $S$ [3]. The optimization problem is:

$$
\max_{\mathbf{n}, T} \; S(\mathbf{n}, T) = -\frac{1}{T} \sum_{i=1}^{N_s} n_i \left[\bar{h}_i^{\circ}(T) - T\bar{s}_i^{\circ}(T) - R_u T \ln\!\left(\frac{x_i P}{P^{\circ}}\right)\right]
$$

subject to:

**Element conservation:**

$$
\sum_{i=1}^{N_s} a_{ij} \, n_i = b_j \quad \text{for } j = 1, 2, \ldots, N_e
$$

**Internal energy constraint:**

$$
\sum_{i=1}^{N_s} n_i \, \bar{u}_i(T) = U
$$

**Equation of state:**

$$
P V = \left(\sum_{i=1}^{N_s} n_i\right) R_u T
$$

**Non-negativity:**

$$
n_i \geq 0 \quad \text{for all } i
$$

where:

| Symbol | Description |
|--------|-------------|
| $\mathbf{n} = [n_1, \ldots, n_{N_s}]^T$ | Species mole numbers |
| $T$ | Temperature (unknown) |
| $\bar{h}_i^{\circ}(T)$ | Standard-state molar enthalpy of species $i$ |
| $\bar{s}_i^{\circ}(T)$ | Standard-state molar entropy of species $i$ |
| $\bar{u}_i(T) = \bar{h}_i^{\circ}(T) - R_u T$ | Molar internal energy of species $i$ |
| $x_i = n_i / \sum_k n_k$ | Mole fraction of species $i$ |
| $a_{ij}$ | Atoms of element $j$ in one molecule of species $i$ |
| $b_j$ | Total moles of element $j$ (conserved) |
| $P^{\circ}$ | Standard pressure (1 bar) |
| $N_s$, $N_e$ | Number of species and elements |

---

## Solution Method

The constrained optimization problem is solved using the method of Lagrange multipliers combined with Newton–Raphson iteration on an extended system that includes temperature as an unknown [1,3].

### Lagrangian Formulation

Introduce Lagrange multipliers $\boldsymbol{\lambda} = [\lambda_1, \ldots, \lambda_{N_e}]^T$ for element conservation and $\lambda_U$ for the energy constraint:

$$
\mathcal{L}(\mathbf{n}, T, \boldsymbol{\lambda}, \lambda_U) = G(\mathbf{n}, T, P) - \sum_{j=1}^{N_e} \lambda_j \left(\sum_{i=1}^{N_s} a_{ij} n_i - b_j\right) - \lambda_U \left(\sum_{i=1}^{N_s} n_i \bar{u}_i(T) - U\right)
$$

where $G$ is the Gibbs free energy and $P$ is determined from the equation of state.

At equilibrium, the stationarity conditions yield:

$$
\mu_i(T, P, \mathbf{x}) = \sum_{j=1}^{N_e} \lambda_j a_{ij} + \lambda_U \bar{u}_i(T) \quad \text{for all } i
$$

$$
\sum_{i=1}^{N_s} a_{ij} n_i = b_j \quad \text{for all } j
$$

$$
\sum_{i=1}^{N_s} n_i \bar{u}_i(T) = U
$$

where $\mu_i$ is the chemical potential of species $i$.

### Newton–Raphson Iteration

The CEA algorithm uses a reduced Newton–Raphson method that iterates on temperature and composition simultaneously [1,2]:

1. Variable transformation — use $\ln n_i$ and $\ln T$ to enforce positivity
2. Reduced system — eliminate dependent variables using element constraints
3. Jacobian construction — compute derivatives with respect to $\{\ln n_i, \ln T\}$
4. Linear solve — solve for the Newton correction $\mathbf{J} \, \Delta \mathbf{z} = -\mathbf{r}$
5. Line search — apply damping factor $\alpha \in (0, 1]$ to ensure descent and feasibility
6. Convergence check — stop when residuals satisfy $\|\mathbf{r}\| < \epsilon$ (typically $\epsilon = 10^{-10}$)

---

## References

[1] Gordon, S., and McBride, B. J. "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 1: Analysis." NASA Reference Publication 1311, 1994.

[2] McBride, B. J., and Gordon, S. "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 2: Users Manual and Program Description." NASA Reference Publication 1311, 1996.

[3] Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*, 2nd edition. John Wiley & Sons, 2003.

---
