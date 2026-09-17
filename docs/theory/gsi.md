# Gas–Surface Interaction

Gas–surface interaction (GSI) describes the coupling between the gas phase and a thermally active or chemically reactive surface. A GSI boundary differs from an inert impermeable wall because the surface can exchange **energy and chemical species** with the gas and can, in some physical regimes, also exchange net mass.

Two broad physical classes are possible:

* **Consumptive surfaces**, for which material is removed from the surface and transferred to the gas, as in melting, pyrolysis, or surface oxidation.
* **Catalytic surfaces**, for which heterogeneous reactions modify the gas composition and reaction rates without net consumption of the surface material.

MOSE implements the **consumptive/ablative class** through three surface models:

| Model | Surface physics | Wall temperature | Surface mass flux |
|-------------------|---------------------------------------------------------|-----------------------------------------|--------------------------------------|
| Melting | phase change and heating of the solid | prescribed at the melting temperature | obtained from the energy balance |
| Pyrolysis | thermally activated decomposition of a solid | solved from the coupled surface balance | temperature-dependent kinetic law |
| Surface reactions | heterogeneous reaction of a carbon surface with the gas | solved from the coupled surface balance | obtained from the reaction mechanism |

These models couple surface thermodynamics and chemistry to the gas-phase conservation equations through the wall mass, species, momentum, and energy fluxes.

---

## Propulsive context

The GSI models in MOSE were selected for **chemical rocket propulsion**, where the chamber or nozzle wall is not an inert boundary but an active participant in the flow. Three distinct propulsion problems are covered.

### Hybrid rocket fuel grains

In a hybrid rocket motor a solid fuel grain is burned with a liquid or gaseous oxidizer injected into the port. There is no premixed propellant: the fuel enters the flow **only** through gasification of the grain surface, and the fuel mass flow rate is therefore an *output* of the boundary condition rather than an input.

Classically this is closed with the Marxman correlation

$$
\dot r = a\,G^{\,n}
$$

where $G$ is the port mass flux. This is an empirical fit whose coefficients must be measured for every fuel, oxidizer and geometry. The GSI models replace it with the underlying physics: the surface energy balance sets the wall temperature, a material-specific kinetic law sets the gasification rate, and the injected mass feeds back on the near-wall transport. The regression rate is then *predicted* rather than prescribed.

Three grain materials are implemented, spanning the two families used in practice:

* **HTPB** — hydroxyl-terminated polybutadiene, the reference hybrid fuel and the standard binder of composite solid propellants. Represented by the two-branch pyrolysis law.
* **HDPE** and **PP** — high-density polyethylene and polypropylene: inexpensive, well-characterised thermoplastics widely used in laboratory-scale and additively manufactured grains. Represented by thermal-wave pyrolysis laws.
* **Liquefying fuels** (paraffin waxes) — represented by the melting model, which supplies the phase-change mass flux (see the caveat on entrainment below).

### Solid rocket propellant

For a solid rocket motor the propellant burning surface itself is not directly modelled due to the complex burning phenomen occurring in a very thin space. Instead, a boundary condition is used to impose mass and energy fluxes. MOSE implements the standard APN law, the mathematical shorthand for Saint Robert’s Law (also called Vieille’s Law), which dictates the linear burning rate of a solid propellant as a function of combustion chamber pressure:

$$
\dot{r} = a\,p^n
$$

where $\dot{r}$ is the linear burning rate (typically in mm/s or in/s), representing how fast the propellant surface regresses perpendicularly to itself; $a$ is the burn rate coefficient (or constant); p is the gas pressure; n is the pressure exponent (or burning index).

### Carbon oxidation

In the rocket propulsion scenario, the carbon oxidation is particularly relevant for the **nozzle throat erosion**: the graphite or carbon–carbon throat insert is attacked chemically by the oxidizing species of the combustion products of hybrid and solid propellants ($\mathrm{H_2O}$, $\mathrm{CO_2}$, $\mathrm{OH}$, $\mathrm{O}$, $\mathrm{O_2}$). Erosion enlarges the throat, lowering chamber pressure and delivered specific impulse over the burn. This is the target of the surface-reaction model.

The same model applies to carbon–carbon and graphite surfaces in any high-enthalpy oxidizing environment like atmospheric re-entry.

---

## Physical picture

The defining feature of a GSI boundary is that the surface response depends on the gas state, while the resulting surface fluxes modify the gas.

$$
\boxed{
\text{gas state}
\;\longrightarrow\;
\text{surface heat / chemical loading / gas pressure}
\;\longrightarrow\;
\text{surface response}
\;\longrightarrow\;
\text{gas-side fluxes}
}
$$

* For the melting, pyrolysis, and surface reactions models, the surface depends on the heat from conduction and radiation. These models are thus thermally coupled GSI. 

* A much simpler relation holds for solid propellants where the material is dependent on the gas pressure.

The central relation for all thermally coupled GSI models is the surface energy balance

$$
\boxed{
q_{\mathrm{cond}}+q_{\mathrm{rad}}
=
q_{\mathrm{surf}}
}
$$

where $q_{\mathrm{cond}}$ is the heat transferred from the gas to the surface, $q_{\mathrm{rad}}$ is externally supplied radiative heating, and $q_{\mathrm{surf}}$ is the energy associated with the surface process.

For melting, $q_{\mathrm{surf}}$ is the energy required to heat and melt the material.

For pyrolysis, it includes sensible heating and the heat of decomposition,

$$
q_{\mathrm{surf}}
=
\dot m
\left[
\Delta h_p+c_s(T_w-T_i)
\right].
$$

For heterogeneous surface reactions, the surface energy exchange follows from the enthalpy carried by the reacting species,

$$
q_{\mathrm{surf}}
=
\sum_s \omega_s h_s(T_w),
$$

as represented by the implemented surface-reaction model.

The wall temperature is therefore either imposed by the physical model or determined as part of the coupled gas–surface problem. In the latter case the balance is a single nonlinear equation in $T_w$ — $q_{\mathrm{cond}}$ falls with increasing $T_w$ while $q_{\mathrm{surf}}$ rises steeply through the Arrhenius factor — and MOSE solves it with a secant iteration at every boundary face and every time step.

---

# Thermally GSI physics

## Melting

The melting model represents a surface held at a prescribed melting temperature,

$$
T_w=T_m.
$$

The available surface energy is used to raise the solid from its initial temperature $T_i$ to the melting point and then supply the latent heat of melting:

$$
\boxed{
\dot m=
\frac{q_{\mathrm{cond}}+q_{\mathrm{rad}}}
{\Delta h_m+c_s(T_m-T_i)}
}
$$

The released material is injected into the gas with a prescribed species composition,

$$
\omega_s=\dot m\,Y_s.
$$

Thus the implemented melting model represents **phase-change-driven ablation**: the wall temperature is prescribed and the recession/mass flux follows from the energy balance.

**Input parameters** (BC 503 data line): $c_s$, $T_m$, $T_i$, $\Delta h_m$, $q_{\mathrm{rad}}$, $\varepsilon_w$, $Y_{1\ldots N_s}$. All material properties are supplied per boundary, so the model is not tied to a particular wax — any material characterised by a melting point, a latent heat and a solid specific heat can be used.

!!! warning "Liquefying fuels: melt flux is not the full regression rate"
    Paraffin-based hybrid fuels regress three to four times faster than classical polymers because the low-viscosity melt layer becomes unstable and sheds droplets into the gas stream — the *entrainment* mechanism. The model above supplies only the mass gasified by the surface energy balance; entrainment of liquid is **not** modelled. For a liquefying fuel the predicted regression rate is therefore a lower bound, and the model is best read as the thermal part of the problem.

---

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

| Symbol | Input | Description |
|---|---|---|
| $T_{\mathrm{af}}$ | `Taf` | adiabatic flame temperature of the propellant, imposed as the wall temperature |
| $a$ | `a` | burn-rate coefficient, as a **mass flux** in kg/(m²·s) |
| $n$ | `n` | pressure exponent |
| $p_{\mathrm{ref}}$ | `pRef` | reference pressure at which $a$ is quoted |
| $Y_{s,\mathrm{af}}$ | `massf(1:nsc)` | mass fractions of the combustion products |

!!! warning "The coefficient is a mass flux, not a linear rate"
    Burn-rate data is normally tabulated as a *linear* rate $\dot r = a_r\,p^n$ in mm/s. The input expected here is $a = \rho_p\,a_r$, evaluated at $p_{\mathrm{ref}}$ and expressed in kg/(m²·s). The grain density is **not** applied internally — `rhoGrain` is read from the data line but never used — so supplying a linear coefficient underestimates the injected mass by three orders of magnitude.

Unlike the gasifying surfaces, whose wall composition is inherited from the adjacent cell under the boundary-layer assumption, the grain surface imposes the **product composition** on the wall state, and the wall transport properties used for the conductive flux are evaluated from it. This is the physically right choice here: the gas at a burning grain surface *is* the flame product mixture.

!!! note "Enthalpy of the injected products"
    The energy flux returns to the gas the enthalpy the products carry,
    $\dot m\,h_{\mathrm{af}}$ with $h_{\mathrm{af}} = \sum_s Y_{s,\mathrm{af}}\,h_s(T_{\mathrm{af}})$,
    evaluated against the thermodynamic tables at the flame temperature. It is *not*
    read from the boundary-condition file: earlier versions took it from an ATLAS
    field that was later dropped, leaving it at zero and injecting the propellant gas
    cold — the dominant term in the surface energy flux, silently absent. Deriving it
    from $T_{\mathrm{af}}$ and the product composition removes the input entirely.

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

This requires the gas mixture to contain $\mathrm{H}$, $\mathrm{OH}$, $\mathrm{CO}$, $\mathrm{CO_2}$, $\mathrm{H_2}$, $\mathrm{H_2O}$, $\mathrm{O_2}$ and $\mathrm{O}$ — the species set of a hydrocarbon or hydrogen combustion product. The steam and Boudouard routes usually dominate in a solid-motor throat, since $\mathrm{H_2O}$ and $\mathrm{CO_2}$ are major products while free $\mathrm{O}$, $\mathrm{OH}$ and $\mathrm{O_2}$ are minor.

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

where $\chi$ is the fraction of the surface covered by the more reactive site type. The first term saturates at high $p_{\mathrm{O_2}}$ (site-limited); the second stays linear.

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

---

# Coupling to the gas phase

The surface models modify the gas through the boundary fluxes.

When there is net material release,

$$
\mathbf u_b
=
\frac{\dot m}{\rho_w}\mathbf n
$$

defines the corresponding blowing velocity.

The released material contributes to the gas species, mass, momentum, and enthalpy.

The resulting transpiration also modifies the near-wall momentum and heat transfer, so the interaction is two-way coupled:

$$
\boxed{
\text{surface chemistry}
\leftrightarrow
\text{boundary-layer transport}
}
$$

For turbulent calculations, the blowing associated with the reactive surface is correspondingly coupled to the wall treatment.

Because the wall injects mass, an ablating face is **not** closed with the reflective (symmetry) flux used for an impermeable wall: the complete boundary flux is assembled from the wall state, carrying the true wall pressure, the viscous stress and the enthalpy of the injected gas. Substituting a mirror-state flux leaves an error in the wall-adjacent cell that does not vanish under mesh refinement, identifying it as a flux inconsistency rather than discretization error; the [verification case](../vv/2D-ablating-wall.md#results) checks the injected mass flux through the wall-adjacent cell for exactly this reason.

## Inert surfaces

Nothing guarantees that a GSI boundary is ablating: if the surface is a net *loser* of heat, the melting balance returns a non-positive mass flux, and the wall would otherwise be made to swallow gas. Every GSI face is therefore tested before its flux is applied, and a non-ablating face falls back to an impermeable wall — reflective in the inviscid flux, isothermal at the surface temperature in the viscous flux.

!!! note "Pyrolysis never switches off exactly"
    An Arrhenius law returns a strictly positive rate at any finite temperature, so a pyrolysing surface always blows *something*. The guard protects against non-physical or non-converged states rather than against a cold wall; at low temperature the mass flux simply becomes negligible.

## Surface recession

The mass, species and energy leaving the solid are supplied to the gas, but **the surface itself does not move**: the mesh is fixed and the grain port does not open up over the burn. Regression rate is a reported output, not a geometric evolution. Motor-scale simulations over a significant fraction of the burn therefore require the geometry to be updated externally.

## Summary

The implemented GSI capability can therefore be summarized as:

| Physical capability                       |            Implemented            |
| ----------------------------------------- | :-------------------------------: |
| Thermally active surface                  |                 ✓                 |
| External radiative heating                |                 ✓                 |
| Solid propellant grain combustion (APN law) |               ✓                 |
| Melting / phase-change ablation           |                 ✓                 |
| Polymer pyrolysis                         |                 ✓                 |
| Temperature-dependent regression laws     |                 ✓                 |
| Heterogeneous carbon reactions            |                 ✓                 |
| Multi-species gas coupling                |                 ✓                 |
| Surface-generated mass flux               |                 ✓                 |
| Species injection from the surface        |                 ✓                 |
| Blowing / transpiration coupling          |                 ✓                 |
| Turbulence coupling with blowing          |                 ✓                 |
| Inert-surface fallback to impermeable wall |                 ✓                 |
| Liquid-layer entrainment (liquefying fuels) |          Not implemented          |
| In-depth solid conduction / charring      |          Not implemented          |
| Surface recession (moving boundary)       |          Not implemented          |
| Erosive burning (rate coupled to cross-flow) |        Not implemented          |
| Surface re-radiation $\varepsilon\sigma T_w^4$ |   Parsed, not applied         |
| Purely catalytic, non-consumptive surface | **Not documented as implemented** |

## Modelling scope and limitations

* **The solid is thermally thin in the model sense.** There is no in-depth conduction solver: the solid is characterised by a bulk temperature $T_i$ and a specific heat, and the sensible term $c_s(T_w-T_i)$ stands in for the whole thermal history. Transients in which the heat soak into the solid matters are outside the model.
* **No char layer.** Pyrolysis products are released at the surface; a porous char through which gas percolates is not represented.
* **Radiative input is prescribed.** $q_{\mathrm{rad}}$ is a constant per boundary, not coupled to a radiation solution or to the local gas state.
* **Pyrolysis constants are compiled in**, and the pyrolysis gas composition is prescribed rather than derived from a decomposition mechanism.
* **Grain combustion is a prescribed-rate model.** Flame temperature, burn-rate law and product composition are all inputs; the propellant flame is not resolved, and the burn rate responds only to pressure — erosive burning, the enhancement of the rate by cross-flow over the grain, is not represented.

## Verification

The melting and pyrolysis boundary conditions are verified against a closed-form solution on a blown-column configuration, in both face orientations, together with the inert-surface fallback: see [2D Ablating Wall](../vv/2D-ablating-wall.md).

The grain-combustion boundary condition is verified on the same configuration, against Saint-Robert's law and the flame temperature it imposes. The surface-reaction model is not covered by the verification suite.

## References

1. Marxman, G. A., and Gilbert, M., "Turbulent Boundary Layer Combustion in the Hybrid Rocket," *9th Symposium (International) on Combustion*, 1963.
2. Karabeyoglu, M. A., Altman, D., and Cantwell, B. J., "Combustion of Liquefying Hybrid Propellants: Part 1, General Theory," *Journal of Propulsion and Power*, Vol. 18, No. 3, 2002.
3. Lengellé, G., "Thermal Degradation Kinetics and Surface Pyrolysis of Vinyl Polymers," *AIAA Journal*, Vol. 8, No. 11, 1970.
4. Nagle, J., and Strickland-Constable, R. F., "Oxidation of Carbon between 1000–2000 °C," *Proceedings of the Fifth Carbon Conference*, Vol. 1, 1962.
5. Bradley, D., Dixon-Lewis, G., El-Din Habik, S., and Mushi, E. M. J., "The Oxidation of Graphite Powder in Flame Reaction Zones," *20th Symposium (International) on Combustion*, 1984.
6. Sutton, G. P., and Biblarz, O., *Rocket Propulsion Elements*, 9th ed., Wiley, 2016.
