# Inlets and outlets

MOSE provides several inlet types, a back-pressure outlet, an ambient boundary, a time-varying assigned state, and a choked-nozzle boundary.

## Type summary

| Type | Meaning | Data line starts with |
|---:|---|---|
| `401` | Inlet at total temperature and pressure | `T0, p0, alpha, beta, rel_fac` |
| `402` | Inlet with time-varying total pressure | `T0, p0_file, alpha, beta, rel_fac` |
| `403` | Inlet at specified mass flux and total temperature | `T0, g, alpha, beta, rel_fac` |
| `404` | Inlet at specified mass flux and static temperature | `T, g, alpha, beta, rel_fac` |
| `405` | Supersonic inlet | `mach, T, p, alpha, beta, rel_fac` |
| `406` | Back-pressure outlet | `p_back, rel_fac` |
| `407` | Ambient boundary | `T0, p0, p_back, alpha, beta, rel_fac` |
| `408` | Inlet at specified normal velocity and static temperature | `T, un, alpha, beta, rel_fac` |
| `410` | Assigned, time-varying state from file | `filename [periodic]` |
| `420` | Choked-nozzle exit | `T0, p0, p_sub, p_sup, g, rel_fac` |

Multi-species cases append one mass fraction per species where required. RANS cases append the specific turbulence variables required by the model.

## Positive direction convention

Mass flux `g` and normal velocity `un` are positive **into the domain** on every face. An inlet only admits gas; it does not become an outlet when the flow reverses. A face where gas leaves should use an outlet (`406`) or an appropriate extrapolating boundary (`400`). A zero or invalid inlet flux/velocity is rejected during setup.

## Behaviour by family

- **Total-condition inlets (`401`, `402`, `407`):** the prescribed reservoir state drives inflow when its total pressure is sufficient to turn the boundary-cell flow back into the domain. If no inflow occurs, `401` and `402` act as slip walls; `407` acts as an outlet at `p_back`.
- **Mass-flux, velocity, and supersonic inlets (`403`, `404`, `405`, `408`):** these types do not permit outflow. If the boundary-cell gas moves outward or is at rest, the face is a slip wall.
- **Back-pressure outlet (`406`):** subsonic outflow uses the relaxed back pressure. Supersonic outflow takes the state from the boundary cell. When the boundary-cell gas moves inward, the face is a slip wall.
- **Ambient (`407`):** behaves as an inlet at total conditions when the surroundings drive inflow and as an outlet at `p_back` when they do not.
- **Choked nozzle (`420`):** selects the subsonic, internal-shock, or design supersonic state from the boundary-cell pressure. The pressure limits must satisfy `0 < p_sup <= p_sub < p0`; the specified mass flux must not exceed the choked reservoir value.
- **Assigned state (`410`):** reads the state across each face cell from a Tecplot ASCII file. It is supported on the finest grid only, not with a multigrid sequence.

## Flow direction angles

The angles `alpha` and `beta` describe the direction of the incoming velocity. The token `normal` requests the inward face normal for the corresponding angle. The exact direction convention for a supersonic inlet is described in the [complete inlet/outlet reference](../boundary-conditions.md).

## Time-varying input files

Type `402` reads a text table of time and total pressure. Type `410` reads a Tecplot ASCII file containing the assigned state at each time. Times must be finite and in increasing order; a leading `periodic` line makes the type `402` pressure table repeat. MOSE interpolates the table linearly in time.

Pressure values are in Pa. MOSE does not infer or convert units from the file.

For the exact column order, required species/turbulence fields, table format, and setup checks, see [Inlet/outlet record details](../boundary-conditions.md).
