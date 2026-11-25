# Viscous walls

The standard viscous wall types are:

| Type | Wall condition | Data line |
|---:|---|---|
| `301` | Adiabatic or prescribed heat flux | `q_wall, k_rough, eps_wall` |
| `302` | Isothermal wall | `T_wall, k_rough, eps_wall` |

| Field | Meaning |
|---|---|
| `q_wall` | Wall heat flux in W/m²; positive into the fluid. Use zero for adiabatic. |
| `T_wall` | Prescribed wall temperature in K. |
| `k_rough` | Sand-grain roughness height in metres. Zero is a smooth wall; the roughness model is available with Spalart–Allmaras. |
| `eps_wall` | Wall emissivity, between zero and one. |

The wall output is written to `OUTPUT/wall.*` when wall output is enabled. The reported traction is the force of the gas on the wall, with the normal directed into the gas. The reported heat flux is positive into the gas.

For melting, pyrolysis, ablation, or surface chemistry, see [Gas–surface interaction](gas-surface.md). Those are separate boundary types (`503`–`506`) with additional physical-model requirements.

For exact defaults, signs, and output variable definitions, see [Viscous-wall records in the complete reference](../boundary-conditions.md#viscous-walls).
