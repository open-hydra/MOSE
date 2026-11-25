# Gas–surface interaction boundaries

Types `502`–`506` couple the gas flow to a solid surface, propellant grain, or surface chemistry model.

| Type | Model |
|---:|---|
| `502` | Solid rocket motor (SRM) grain boundary |
| `503` | Melting surface |
| `504` | Pyrolysis surface |
| `505` | Surface reactions |
| `506` | Pyrolysing ablator: surface reactions plus gas blown through the char layer |

## Important model requirements

- Type `502` can be used with either Euler or Navier–Stokes equations. It can represent a blowing wall or, in mode 2, an inflow.
- Types `503`–`506` require `equations = navier-stokes`. With Euler equations, they fall back to an impermeable slip wall and no ablation takes place.
- Types `505` and `506` require a thermodynamic table containing O2, CO2, CO, H2O, H2, OH, O, and H.
- The gas composition and thermal state must be compatible with the thermodynamic tables. MOSE checks required properties during setup.

For a `502` grain in mode 1, injected gas carries the turbulence values of a smooth wall when a RANS model is active. In mode 2, the face behaves as an inlet (`404`) from setup onward and reads the specific turbulence values from its data line.

When the flow cannot sustain a surface process, MOSE can use the inert-surface fallback described in the [gas–surface theory](../../theory/gsi.md#inert-surfaces).

## Files and parameters

Type `506` can read the pyrolysis-gas composition from `INPUT/<prefix>pyrolysis.dat` when `phi > 0`. Its `phi` parameter is the ratio of pyrolysis-gas to char mass flux; `h_v0` is the virgin-material enthalpy on the datum of the gas tables.

The complete record table gives the remaining fields, their order, defaults, units, and model-specific checks. See [Gas–surface interaction records](../boundary-conditions.md) and the [surface-model theory](../../theory/gsi-models.md).
