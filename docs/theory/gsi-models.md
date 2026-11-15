# Gas–Surface Interaction: Surface Models in Detail

This page gives the surface energy and mass balances of the models summarised in [Gas–Surface Interaction](gsi.md): pyrolysis,
solid propellant combustion, heterogeneous surface reactions, and the pyrolysing ablator.

## Pyrolysis

The pyrolysis model represents thermal decomposition of a solid fuel: the polymer chain is broken by heat into gaseous fragments, which leave the surface and burn in the boundary layer.

The mass flux is a function of wall temperature,

$$
\dot m=\dot m(T_w),
$$

and the surface energy balance is

$$
q_{\mathrm{cond}}+q_{\mathrm{rad}}
=
\dot m(T_w)
\left[
\Delta h_p+c_s(T_w-T_i)
\right].
$$

The wall temperature is therefore an unknown of the coupled problem. The bracket is the total energy cost of producing a unit mass of pyrolysis gas: the heat of decomposition $\Delta h_p$ plus the sensible heat needed to bring the solid from its bulk temperature $T_i$ to the reacting surface temperature $T_w$.

All three models share the same interface — given $T_w$ they return $\dot m$, the species source terms $\omega_s = \dot m Y_s$, and $q_{\mathrm{surf}}$ — and differ only in the kinetic law. The product composition $Y_s$ is **prescribed per boundary**, not derived from a decomposition mechanism: the user supplies the mass fractions of the pyrolysis gas, which then reacts in the gas phase through the ordinary finite-rate [kinetics](kinetics.md).

### HTPB

HTPB is the workhorse hybrid fuel and the binder of AP/Al composite solid propellants. Its decomposition proceeds by two distinct chemical routes — depolymerisation of the butadiene backbone at lower temperature, and scission of the resulting fragments above it — so the regression rate is fitted by a **two-branch Arrhenius law** in the surface regression velocity $\dot r$, converted to a mass flux through the solid density:

$$
\dot m = \rho_f\,\dot r,
\qquad
\dot r =
\begin{cases}
A_1
\exp\!\left(-\dfrac{E_1}{\mathcal RT_w}\right),
& T_w\leq T_s,\\[2ex]
A_2
\exp\!\left(-\dfrac{E_2}{\mathcal RT_w}\right),
& T_w>T_s.
\end{cases}
$$

| Symbol | Value | Units | Description |
|---|---|---|---|
| $\rho_f$ | 960 | kg/m³ | solid density |
| $c_s$ | 1632 | J/(kg·K) | solid specific heat |
| $\Delta h_p$ | 1.1 × 10⁶ | J/kg | heat of pyrolysis |
| $T_i$ | 298.15 | K | bulk solid temperature |
| $A_1$ | 3.965 | m/s | pre-exponential, low-temperature branch |
| $E_1$ | 55.8564 × 10⁶ | J/kmol | activation energy, low-temperature branch ($E_1/\mathcal R$ = 6718 K) |
| $A_2$ | 11.04 × 10⁻³ | m/s | pre-exponential, high-temperature branch |
| $E_2$ | 20.54344 × 10⁶ | J/kmol | activation energy, high-temperature branch ($E_2/\mathcal R$ = 2471 K) |
| $T_s$ | 722 | K | switch temperature between branches |

The two branches are matched at $T_s$: evaluated there they agree to within 0.2 %, so $\dot r(T_w)$ is continuous to the precision of the fit and the secant iteration does not see a jump when it crosses the switch. The high-temperature branch is much less sensitive to $T_w$ ($E_2 \approx E_1/2.7$), which is the physical statement that above $T_s$ the rate is limited by heat supply rather than by chemistry.

### HDPE

HDPE is described by a single Arrhenius-type law returning the mass flux directly, without an explicit solid density:

$$
\dot m
=
A\exp\!\left(
-\frac{E}{2\mathcal RT_w}
\right).
$$

| Symbol | Value | Units | Description |
|---|---|---|---|
| $c_s$ | 1255.2 | J/(kg·K) | solid specific heat |
| $\Delta h_p$ | 2.72 × 10⁶ | J/kg | heat of pyrolysis |
| $T_i$ | 298.15 | K | bulk solid temperature |
| $A$ | 4.5888 × 10⁶ | kg/(m²·s) | pre-exponential factor |
| $E$ | 251.03973 × 10⁶ | J/kmol | activation energy ($E/2\mathcal R$ = 15096 K) |

!!! note "The factor of two in the exponent"
    The $2\mathcal R T_w$ in the denominator is not a typographical quirk: it is the signature of a **thermal-wave** (in-depth) pyrolysis analysis, in which the regression rate comes out as the square root of the volumetric reaction rate, $\dot r \propto \sqrt{A'\exp(-E/\mathcal R T_w)}$, halving the effective activation energy. The HDPE law is that result with the remaining thermal factors absorbed into $A$; the PP law below keeps them explicit.

### PP

Polypropylene uses the same thermal-wave physics written out in full. The pyrolysis front is treated as a reaction zone embedded in the conduction wave that precedes it, and the regression rate follows from balancing the two:

$$
\dot m = \rho_f
\sqrt{
\frac{A\,e^{-\Xi}}{\Xi}\,
\frac{\alpha}{\Theta}
},
\qquad
\Xi=\frac{E}{\mathcal RT_w},
\qquad
\alpha = \frac{k_s}{\rho_f c_s}.
$$

The factor $\Theta$ is the dimensionless thickness of the reacting layer, defined by the depth over which the reaction rate falls by two decades:

$$
\Theta=
\sqrt{
\ln(100)\left(1-\frac{T_i}{T_w}+\beta\right)-\beta
},
\qquad
\beta = \frac{\Delta h_p}{c_s T_w}.
$$

| Symbol | Value | Units | Description |
|---|---|---|---|
| $\rho_f$ | 910 | kg/m³ | solid density |
| $c_s$ | 1700 | J/(kg·K) | solid specific heat |
| $k_s$ | 0.2 | W/(m·K) | solid thermal conductivity ($\alpha$ = 1.29 × 10⁻⁷ m²/s) |
| $\Delta h_p$ | 2.4823 × 10⁶ | J/kg | heat of pyrolysis |
| $T_i$ | 298.15 | K | bulk solid temperature |
| $A$ | 2.12 × 10¹⁵ | 1/s | pre-exponential factor |
| $E$ | 2.12 × 10⁸ | J/kmol | activation energy ($E/\mathcal R$ = 25498 K) |

$\Theta$ is a weak, slowly varying function of $T_w$, so the temperature dependence is dominated as usual by $e^{-\Xi/2}$ — consistent with the HDPE form.

---

## Solid propellant combustion

The grain combustion model stands apart from the other three: it is the only one whose mass flux does **not** come from a surface energy balance. A solid propellant carries its own oxidizer, so the surface is a deflagration front whose rate is set by the flame standing off it, and both the rate and the flame temperature are supplied as measured properties of the propellant rather than derived from the gas state.

The regression rate follows Saint-Robert's (Vieille's) law. MOSE evaluates it in **mass-flux form, normalised by a reference pressure**:

$$
\boxed{
\dot m = a\left(\frac{p_c}{p_{\mathrm{ref}}}\right)^{\!n}
}
$$

so that $a$ is the propellant mass flux at $p_c = p_{\mathrm{ref}}$. The wall is then held at the adiabatic flame temperature and the combustion products are injected with a prescribed composition:

$$
T_w = T_{\mathrm{af}},
\qquad
\omega_s = \dot m\,Y_{s,\mathrm{af}}.
$$

| Symbol | Description |
|---|---|
| $T_{\mathrm{af}}$ | adiabatic flame temperature of the propellant, imposed as the wall temperature |
| $a$ | burn-rate coefficient, as a **mass flux** in kg/(m²·s) |
| $n$ | pressure exponent |
| $p_{\mathrm{ref}}$ | reference pressure at which $a$ is quoted |
| $Y_{s,\mathrm{af}}$ | mass fractions of the combustion products |

!!! warning "The coefficient is a mass flux, not a linear rate"
    Burn-rate data is normally tabulated as a *linear* rate $\dot r = a_r\,p^n$ in mm/s. The coefficient of the model is $a = \rho_p\,a_r$, evaluated at $p_{\mathrm{ref}}$ and expressed in kg/(m²·s). The grain density is **not** applied by the model, so a linear coefficient in its place underestimates the injected mass by three orders of magnitude.

Unlike the gasifying surfaces, whose wall composition is inherited from the adjacent cell under the boundary-layer assumption, the grain surface imposes the **product composition** on the wall state, and the wall transport properties used for the conductive flux are evaluated from it. This is the physically right choice here: the gas at a burning grain surface *is* the flame product mixture.

!!! note "Enthalpy of the injected products"
    The energy flux returns to the gas the enthalpy the products carry,
    $\dot m\,h_{\mathrm{af}}$ with $h_{\mathrm{af}} = \sum_s Y_{s,\mathrm{af}}\,h_s(T_{\mathrm{af}})$,
    evaluated against the thermodynamic tables at the flame temperature: it follows
    from $T_{\mathrm{af}}$ and the product composition, and is the dominant term of the
    surface energy flux.

---

## Heterogeneous surface reactions

The surface-reaction model represents the heterogeneous attack of a **carbon surface** — graphite or carbon–carbon, as used for nozzle throat inserts — by the oxidizing species of the gas. Unlike pyrolysis, no material property is prescribed: the mass flux is an outcome of the reaction mechanism and the local gas composition, and the product composition is fixed by stoichiometry rather than by the user.

### Reaction set

Five surface reactions consume solid carbon:

| Reaction | Oxidizer | Rate |
|---|---|---|
| $\mathrm{C_{(s)}} + \mathrm{H_2O} \rightarrow \mathrm{CO} + \mathrm{H_2}$ | steam | $\dot m_{\mathrm{H_2O}}$ |
| $\mathrm{C_{(s)}} + \mathrm{CO_2} \rightarrow 2\,\mathrm{CO}$ | Boudouard | $\dot m_{\mathrm{CO_2}}$ |
| $\mathrm{C_{(s)}} + \tfrac12\mathrm{O_2} \rightarrow \mathrm{CO}$ | molecular oxygen | $\dot m_{\mathrm{O_2}}$ |
| $\mathrm{C_{(s)}} + \mathrm{OH} \rightarrow \mathrm{CO} + \mathrm{H}$ | hydroxyl | $\dot m_{\mathrm{OH}}$ |
| $\mathrm{C_{(s)}} + \mathrm{O} \rightarrow \mathrm{CO}$ | atomic oxygen | $\dot m_{\mathrm{O}}$ |

This requires the gas mixture to contain $\mathrm{H}$, $\mathrm{OH}$, $\mathrm{CO}$, $\mathrm{CO_2}$, $\mathrm{H_2}$, $\mathrm{H_2O}$, $\mathrm{O_2}$ and $\mathrm{O}$ — the species set of a hydrocarbon or hydrogen combustion product. A run with a wall `505` or `506` whose gas lacks one of them stops at set-up with a message that names the missing species; without $\mathrm{O_2}$ in the gas ($p_{\mathrm{O_2}} = 0$) the oxygen rate is 0 ($\chi = 0$). The steam and Boudouard routes usually dominate in a solid-motor throat, since $\mathrm{H_2O}$ and $\mathrm{CO_2}$ are major products while free $\mathrm{O}$, $\mathrm{OH}$ and $\mathrm{O_2}$ are minor.

### Rate laws

Rates are written in terms of **partial pressures in atmospheres**, $p_s = \rho_s R_s T_w / p_{\mathrm{ref}}$ with $p_{\mathrm{ref}}$ = 1.01325 × 10⁵ Pa, and return a carbon mass flux in kg/(m²·s):

$$
\dot m_{\mathrm{H_2O}} = k_{\mathrm{H_2O}}\sqrt{p_{\mathrm{H_2O}}},
\qquad
\dot m_{\mathrm{CO_2}} = k_{\mathrm{CO_2}}\sqrt{p_{\mathrm{CO_2}}},
$$

$$
\dot m_{\mathrm{OH}} = k_{\mathrm{OH}}\,p_{\mathrm{OH}},
\qquad
\dot m_{\mathrm{O}} = k_{\mathrm{O}}\,p_{\mathrm{O}}.
$$

The half-order dependence of the steam and Boudouard routes reflects dissociative adsorption of the oxidizer on the carbon surface; the radical routes are first order because $\mathrm{OH}$ and $\mathrm{O}$ react on impact. The radical rate coefficients carry no activation energy at all — they are collision-limited, with the $T_w^{-1/2}$ of a gas-kinetic flux:

| Coefficient | Expression | Activation temperature |
|---|---|---|
| $k_{\mathrm{H_2O}}$ | $4.80\times10^{5}\exp\!\left(-E/\mathcal RT_w\right)$, $E$ = 288.00 × 10⁶ J/kmol | 34 640 K |
| $k_{\mathrm{CO_2}}$ | $9.00\times10^{3}\exp\!\left(-E/\mathcal RT_w\right)$, $E$ = 285.00 × 10⁶ J/kmol | 34 279 K |
| $k_{\mathrm{OH}}$ | $3.61\times10^{2}\,T_w^{-1/2}$ | — |
| $k_{\mathrm{O}}$ | $6.655\times10^{2}\,T_w^{-1/2}$ | — |

The very large activation temperatures of the $\mathrm{H_2O}$ and $\mathrm{CO_2}$ routes are why throat erosion is negligible at moderate wall temperature and grows abruptly once the insert approaches its operating temperature.

### Oxidation: two-site mechanism

Molecular oxygen is handled with a two-site surface mechanism — the Nagle–Strickland-Constable form — in which the carbon surface presents reactive sites of two types whose relative population shifts with temperature and oxygen pressure:

$$
\dot m_{\mathrm{O_2}} =
\chi\,\frac{k_5\,p_{\mathrm{O_2}}}{1+k_6\,p_{\mathrm{O_2}}}
+
(1-\chi)\,k_7\,p_{\mathrm{O_2}},
\qquad
\chi=\frac{1}{1+\dfrac{k_8}{k_7\,p_{\mathrm{O_2}}}}
$$

where $\chi$ is the fraction of the surface covered by the more reactive site type. The first term saturates at high $p_{\mathrm{O_2}}$ (site-limited); the second stays linear. Bianchi and Nasuti (J. Propuls. Power 29(3), 2013, eq. 7) print $1+k_6\,\chi$ in the denominator; their own Fig. 1a, at 1 atm, follows the form above within about 0.1 % from 1000 to 4000 K and lies up to a factor 4.7 below the printed form (at 3000 K): the print is a misprint. MOSE had $1+k_6\,p_{\mathrm{O_2}}\,\chi$ in the denominator from September 2026 until this correction (`CHANGELOG.md`).

| Coefficient | Expression | Note |
|---|---|---|
| $k_5$ | $2.40\times10^{3}\exp\!\left(-125.60\times10^{6}/\mathcal RT_w\right)$ | reactive-site rate |
| $k_6$ | $2.13\times10^{1}\exp\!\left(+17.17\times10^{6}/\mathcal RT_w\right)$ | site saturation — note the **positive** exponent |
| $k_7$ | $5.35\times10^{-1}\exp\!\left(-63.64\times10^{6}/\mathcal RT_w\right)$ | less-reactive-site rate |
| $k_8$ | $1.81\times10^{7}\exp\!\left(-406.10\times10^{6}/\mathcal RT_w\right)$ | thermal site conversion |

The positive exponent of $k_6$ is physical, not a sign slip: it describes desorption-controlled coverage, which *falls* with temperature.

### Source terms and surface energy

Each route contributes to the gas-phase species according to its stoichiometry, scaled by the molecular-weight ratio to the carbon atom consumed:

$$
\omega_s = \frac{W_s}{W_{\mathrm C}}\,\nu_s\,\dot m_{(\cdot)},
\qquad
\dot m = \sum_s \omega_s,
\qquad
q_{\mathrm{surf}} = \sum_s \omega_s\,h_s(T_w).
$$

The net mass flux is the **sum of the species sources**, which is exactly the carbon removed from the surface: the oxidizers consumed carry negative $\omega_s$, the products positive, and the difference is the eroded solid. The surface energy term is the enthalpy flux of that species exchange evaluated at the wall temperature, so an exothermic reaction set heats the surface and feeds back into the wall-temperature solution.

## Pyrolysing ablator (BC 506)

A charring ablator (a carbon-phenolic insert, say) is modelled at its surface as the carbon wall of the heterogeneous
reactions above plus the gas of the pyrolysing virgin material, which percolates through the char and is blown into the
boundary layer (Rotondi, PhD thesis, §2.4.2–2.4.5). With the char mass flux $\dot m_{ch}$ of the reactions and the
pyrolysis-gas mass flux $\dot m_{pg}$, the surface mass and species balances are (eqs. 2.100, 2.101, 2.115)

$$
\dot m_w = \dot m_{ch} + \dot m_{pg},
\qquad
\dot m_{pg} = \phi\,\dot m_{ch},
\qquad
\phi = \frac{\rho_v}{\rho_{ch}} - 1,
\qquad
\omega_s = \omega_{ch,s} + \dot m_{pg}\,y_{pg,s}(p_w, T_w):
$$

in steady state the material between the virgin and the char densities leaves as pyrolysis gas, whose composition
$y_{pg}$ is the equilibrium of the elemental makeup of the virgin material at the wall pressure and temperature, read
from a table (`INPUT/<prefix>pyrolysis.dat`, written by ATLAS GPB, read by FLINT `FLINT_Lib_PyrolysisGas`: bilinear in
$(p, T)$, clamped to its range; a species of the gas absent from the table gets 0). Its enthalpy is
$h_{pg} = \sum_s y_{pg,s}\,h_s(T_w)$ (eq. 2.116).

The solid conducts in steady state, so the heat it takes from the
surface is $q_{cond} = \dot m_{ch} h_{ch} + \dot m_{pg} h_{pg} - (\dot m_{ch} + \dot m_{pg})\,h_{v,0}$ (eq. 2.95),
$h_{v,0}$ the enthalpy of the virgin material at its undisturbed temperature on the datum of the gas tables, and the
surface energy balance (eq. 2.96) is solved for $T_w$ as for the reactions, with the heat the surface absorbs

$$
q_{\mathrm{surf}} = \sum_s \omega_s\,h_s(T_w) - \dot m_w\,h_{v,0}
= q_{ch} + \dot m_{pg}\,h_{pg} - (\dot m_{ch} + \dot m_{pg})\,h_{v,0}.
$$

The table is read again at every trial temperature of the balance, at the pressure of the boundary cell (the wall
pressure of the boundary-layer assumption). With $\phi = 0$ no pyrolysis gas is blown and the table is not read, but the
virgin-material term stays, $-\dot m_{ch} h_{v,0}$; with $h_{v,0} = 0$ too, BC 506 is BC 505 bit for bit. The erosion
rate follows from the total mass flux of the wall output, $\dot r = \dot m_w/\rho_v$ (eq. 2.91).

Differences from the original implementation of this model (a bit-for-bit comparison with it is not a test): the secant
of the wall temperature takes 20 iterations with a perturbation of $10^{-4}$ (there 10 and $10^{-3}$); the wall pressure
is part of the GSI flux, with no separate reflective flux; RANS runs take the blowing-wall treatment; any face may be the
wall (there the fourth only); $q_{\mathrm{rad}}$ is prescribed, without re-radiation; the virgin-material term
$-\dot m_w h_{v,0}$ enters the balance of every record, also with $\phi = 0$ (there only when $\phi \neq 0$), and
$h_{v,0}$ is an input of the record (there fixed at $-8.57\cdot10^5$ J/kg); the interpolation of the table is
bilinear; the table is in Pa and names its species; its default treatment is gas-only equilibrium (the elements of the
makeup all leave as gas), the treatment of the thesis (graphite allowed to form and removed, the gas renormalised) is
the `strip` treatment of the generator.

Not modelled: an imposed wall temperature or conduction flux instead of steady
conduction, soot in the blown gas (the `condensed` treatment, refused by the wall), $\dot m_{pg}$ from a solid solver.

---

## References

## References

The works cited are listed in [Gas–Surface Interaction](gsi.md#references).
