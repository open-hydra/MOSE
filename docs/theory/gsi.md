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

For a solid rocket motor the propellant burning surface itself is not directly modelled due to the complex burning phenomena occurring in a very thin space. Instead, a boundary condition is used to impose mass and energy fluxes. MOSE implements the standard APN law, the mathematical shorthand for Saint Robert’s Law (also called Vieille’s Law), which dictates the linear burning rate of a solid propellant as a function of combustion chamber pressure:

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

## Thermally coupled surface models

### Melting

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

**Model parameters**: $c_s$, $T_m$, $T_i$, $\Delta h_m$, $q_{\mathrm{rad}}$, $\varepsilon_w$, $Y_{1\ldots N_s}$. All material properties are given per boundary, so the model is not tied to a particular wax — any material characterised by a melting point, a latent heat and a solid specific heat can be used.

!!! warning "Liquefying fuels: melt flux is not the full regression rate"
    Paraffin-based hybrid fuels regress three to four times faster than classical polymers because the low-viscosity melt layer becomes unstable and sheds droplets into the gas stream — the *entrainment* mechanism. The model above supplies only the mass gasified by the surface energy balance; entrainment of liquid is **not** modelled. For a liquefying fuel the predicted regression rate is therefore a lower bound, and the model is best read as the thermal part of the problem.

---

### Pyrolysis

The pyrolysis model represents the thermal decomposition of a solid fuel: heat breaks the polymer chain into gaseous fragments, which
leave the surface and burn in the boundary layer. The mass flux is a function of the wall temperature, which the surface energy
balance fixes. The kinetic laws and the balance are in [Gas–Surface Interaction: Surface Models in Detail](gsi-models.md#pyrolysis).

### Solid propellant combustion

The grain model stands apart from the others: its mass flux does not come from a surface energy balance. A solid propellant carries
its own oxidizer, so the surface is a deflagration front whose rate is set by the flame standing off it. The rate, in the mass-flux form
of Saint-Robert's law, and the flame temperature are supplied as properties of the propellant. The model is in [Gas–Surface Interaction: Surface Models in Detail](gsi-models.md#solid-propellant-combustion).

### Heterogeneous surface reactions

The surface-reaction model represents the heterogeneous attack of a carbon surface (graphite or carbon–carbon, as in nozzle throat
inserts) by the oxidizing species of the gas. No material property is prescribed: the mass flux follows from the reaction mechanism
and the local gas composition, and the product composition from stoichiometry. The mechanism is in [Gas–Surface Interaction: Surface Models in Detail](gsi-models.md#heterogeneous-surface-reactions).

### Pyrolysing ablator (BC 506)

A charring ablator is modelled at its surface as the carbon wall of the heterogeneous reactions plus the gas of the pyrolysing virgin
material, which percolates through the char and is blown into the boundary layer. The balances are in [Gas–Surface Interaction: Surface Models in Detail](gsi-models.md#pyrolysing-ablator-bc-506).

## Coupling to the gas phase

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

### Inert surfaces

Nothing guarantees that a GSI boundary is ablating: if the surface is a net *loser* of heat, the melting balance returns a non-positive mass flux, and the wall would otherwise be made to swallow gas. Every GSI face is therefore tested before its flux is applied, and a non-ablating face falls back to an impermeable wall — reflective in the inviscid flux, isothermal at the surface temperature in the viscous flux.

!!! note "Pyrolysis never switches off exactly"
    An Arrhenius law returns a strictly positive rate at any finite temperature, so a pyrolysing surface always blows *something*. The guard protects against non-physical or non-converged states rather than against a cold wall; at low temperature the mass flux simply becomes negligible.

### Surface recession

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
| Pyrolysing ablator: char reactions and blown pyrolysis gas, steady conduction | ✓ |
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
* **No char layer.** Pyrolysis products are released at the surface; a porous char through which gas percolates is not represented. The pyrolysing ablator (506) takes the char and the pyrolysis gas into its surface balances, with steady conduction in the solid, not an in-depth solution.
* **Reynolds-stress turbulence at a reacting wall.** With the SSG/LRR Reynolds-stress model a wall 505 or 506 is not covered: no test runs it, and the blowing-wall treatment of its turbulence wall values was not checked for that model.
* **Radiative input is prescribed.** $q_{\mathrm{rad}}$ is a constant per boundary, not coupled to a radiation solution or to the local gas state.
* **Pyrolysis constants are compiled in**, and the pyrolysis gas composition is prescribed rather than derived from a decomposition mechanism.
* **Grain combustion is a prescribed-rate model.** Flame temperature, burn-rate law and product composition are all inputs; the propellant flame is not resolved, and the burn rate responds only to pressure — erosive burning, the enhancement of the rate by cross-flow over the grain, is not represented.

## Verification

The melting and pyrolysis boundary conditions are verified against a closed-form solution on a blown-column configuration, in both face orientations, together with the inert-surface fallback: see [2D Ablating Wall](../vv/2D-ablating-wall.md).

The grain-combustion boundary condition is verified on the same configuration, against Saint-Robert's law and the flame temperature it imposes. The surface reactions (505) and the pyrolysing ablator (506) are verified on the same column by the relations of the wall model at the state the solver reached (the column has no steady state: the reactions consume its oxidisers), together with the identity of 506 without pyrolysis and 505, and by the unit checks `GSIBradley` and `GSIAblator`.

!!! info "Going deeper"
    [Gas–Surface Interaction: Surface Models in Detail](gsi-models.md): the energy and mass balances of pyrolysis, solid propellant combustion, heterogeneous surface reactions and the pyrolysing ablator.

## References

[1] Marxman, G. A., and Gilbert, M., "Turbulent Boundary Layer Combustion in the Hybrid Rocket," *9th Symposium (International) on Combustion*, 1963.

[2] Karabeyoglu, M. A., Altman, D., and Cantwell, B. J., "Combustion of Liquefying Hybrid Propellants: Part 1, General Theory," *Journal of Propulsion and Power*, Vol. 18, No. 3, 2002.

[3] Lengellé, G., "Thermal Degradation Kinetics and Surface Pyrolysis of Vinyl Polymers," *AIAA Journal*, Vol. 8, No. 11, 1970.

[4] Nagle, J., and Strickland-Constable, R. F., "Oxidation of Carbon between 1000–2000 °C," *Proceedings of the Fifth Carbon Conference*, Vol. 1, 1962.

[5] Bradley, D., Dixon-Lewis, G., El-Din Habik, S., and Mushi, E. M. J., "The Oxidation of Graphite Powder in Flame Reaction Zones," *20th Symposium (International) on Combustion*, 1984.

[6] Sutton, G. P., and Biblarz, O., *Rocket Propulsion Elements*, 9th ed., Wiley, 2016.
