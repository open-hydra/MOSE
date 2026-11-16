# Multirate Time Stepping

The multirate integrator lets every cell advance with a time step close to its own stability limit, instead of the smallest step of the mesh. Cells are grouped in levels: level $l$ takes steps $2^l$ times smaller than level 0, and the levels are coupled at the faces between them ([theory](../../theory/multirate.md#multirate-time-stepping)).

The verification answers four questions:

1. **With a single level, is it the same solver?** Yes: the multirate step then reproduces the single-rate step bit for bit.
2. **Does it conserve mass, momentum and energy?** Yes, to round-off ($10^{-12}$), also across faces between levels, blocks and processes.
3. **Is it second order in time with several levels?** Yes on a fixed level map, and close to it when the faces between levels move with the flow: see [Temporal order](order.md).
4. **On real flows, is it as accurate as the single rate?** Yes, within 5% and usually much closer: see the [viscous shock tube](daru-tenaud.md) and the [nozzle start-up](nozzle-startup.md).

## How the accuracy is measured

No manufactured solution exercises the multirate integrator. The order in time is measured instead on a fixed mesh, against a single-rate run with a much smaller step on the same mesh. This isolates the error of the time integration and of the coupling between levels from the error of the spatial scheme. The spatial scheme needs no separate check: with one level the integrator *is* the single-rate solver, whose order is verified by the [manufactured solutions](../2D-mms.md) and the [isentropic vortex](../2D-isentropic-vortex.md).

## Accuracy

| Property | Case | Result |
|---|---|---|
| Order in time, fixed level map | 8 cases: acoustic and entropy waves, vortex across blocks, laminar plate, hot square, skewed wall, SA and SST channels | 1.91–2.19 ([details](order.md#fixed-level-map)) |
| Order in time, moving faces between levels | simple wave, mesh and step refined together | $\ge 1.9$ in $L^1$ ([details](order.md#moving-faces-between-levels)) |
| Order at boundaries | duct with an unsteady inlet (`402`, `410`) | $\ge 1.8$ in the inlet cells |
| Shocks against the exact solution | Sod tube, planar shock in 2D | error $\le 1.05\times$ the single rate's; shock within one cell; no new oscillations |
| Viscous flow with walls | Daru–Tenaud shock tube | within 0.2% of the single rate; wall density within 0.8% of the published solution ([details](daru-tenaud.md)) |
| Moving shock across levels | nozzle start-up | error $\le 1.05\times$ the single rate's ([details](nozzle-startup.md)) |
| Deep hierarchies | 9 and 13 levels | second order, stable over 2000 macro steps |
| Symmetry | pressure pulse in a symmetric box | symmetric to $10^{-10}$ |
| Turbulence positivity | SA and SST channels with stiff near-wall sources | $\tilde\nu$, $k$, $\omega$ never negative |

## Implementation

These checks prove that the integrator does exactly what it should, and they are exact (bit for bit or to round-off):

| Property | Result |
|---|---|
| One level = single rate, for every accepted boundary condition | identical bit for bit |
| Conservation with levels, blocks, axisymmetric faces, SA and SST | drift $\le 10^{-12}$ |
| Uniform flow on a curved, skewed mesh (free stream preservation) | preserved to $10^{-13}$ |
| One block against several blocks with the same levels | same solution to $10^{-11}$ |
| OpenMP threads, MPI processes and their combinations | identical bit for bit |
| Restart in the middle of a run | identical bit for bit |
| Order of the stages in one macro step | equal to the reference sequence |
| Cells that each level reads from its neighbours (halo) | equal to an independent construction |
| Fixed and adaptive macro step with the same value | identical bit for bit |
| Non-local stability condition of Gnedin et al. | detected and corrected as designed |
| Configurations the integrator does not support (35 of them) | rejected with a clear message |

The step of every cell is limited with a safety factor of 2: with a factor of 1 a shock that strengthens during a macro step can outrun the step chosen at its start ([theory](../../theory/multirate-method.md#non-local-cfl-condition)).

The cases, the checks and their exact acceptance criteria, and how to run them, are in [`test/multirate/README.md`](https://github.com/open-hydra/MOSE/blob/main/test/multirate/README.md).
