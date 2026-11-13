# Multirate Time Stepping

In a time-accurate run with the global step, every cell advances with the smallest admissible step of the mesh,
$\Delta t = \min_i \Delta t_i$. When a few cells (near walls, in a throat, in a refined band) set that minimum, most of the work
goes into advancing the other cells with steps far below their own limit.

The **multirate** integrator lets each cell advance with a step close to its own limit. On a fixed mesh it keeps the second
order in time and the exact conservation of the single-rate scheme.

## The idea

Each cell $c$ gets a **level** $\ell_c \in \{0, \dots, L\}$ and advances with sub-steps

$$
h_\ell = H/2^{\ell},
$$

where $H$ is the **macro step**. Cells with a small admissible step sit at a deep level and take many small sub-steps; the
others take few. All cells meet again at the end of the macro step.

## How the levels are chosen

The level of a cell comes from its admissible local step $\Delta t^{\text{cap}}_c$ (the
[local time step](time-integration.md#local-time-step)) and a safety factor $\sigma \ge 1$, 2 by default:

$$
\ell_c = \bigl\lceil \log_2\bigl(\sigma H/\Delta t^{\text{cap}}_c\bigr) \bigr\rceil ,
\qquad
H = \min_c \Delta t^{\text{cap}}_c\,2^{\ell_c}/\sigma ,
$$

with the level kept between 0 and the deepest level allowed. Two rules are then imposed:

- **nesting**: the levels of two cells that share a face differ by at most one;
- **equal levels across block connections**: the two sides of a connection always advance together.

The map of levels is either rebuilt from the local steps at every macro step (a *recomputed map*, for production runs) or
built once and kept (a *fixed map*, for studies of the order in time).

## The scheme

Within a macro step the integrator is a recursion of Heun's scheme (the two-stage SSP Runge–Kutta of
[Time Integration](time-integration.md)). A cell at level $\ell$ does two stages per sub-step; between the two stages of a
coarse cell, the finer cells next to it advance two sub-steps of their own. A coarse cell next to a finer one reads values
interpolated to the stage time of the fine cell, so that both sides see the same state at the same time.

## Conservation: the flux correction

At a face between a coarse and a fine cell, the two sides would integrate different fluxes. The coarse side is corrected so
that it receives exactly the flux that the fine side received over the same time. Mass, momentum and energy are then conserved
to round-off, as in the single-rate scheme.

## Keeping the steps admissible

The CFL condition is not local: within a step a wave crosses the faces of a cell and reaches neighbours that advance with their
own steps. Gnedin, Semenov and Kravtsov require that, if the admissible step of a cell at time $t_0$ is $\Delta t$, its
immediate neighbours end their step at or before $t_0 + \Delta t$ (condition (4)). MOSE keeps to it in two ways:

- the **interface margin** keeps the coarse step at a face between levels within $\Delta t^{\text{cap}}/\kappa$ of the finer
  cell, with $\kappa \ge 1.1$;
- the **redo**: a macro step that violates the condition is rejected and redone with finer levels or a smaller $H$.

## What is guaranteed

- **Conservation** to round-off, on every map.
- **Second order in time** on a fixed mesh: proven.
- **Stability**: proven for linear two-cell fluxes (first-order upwind, central diffusion), and for the limited reconstruction
  of the advection equation in one dimension (two levels, and three levels with a band of one cell). Elsewhere it rests on the
  tests of [Multirate Verification](../vv/multirate.md).
- With one level, the macro step is the single-rate step to the bit.

The input keys are in the [User Guide](../user/input.md#multirate-time-stepping), which also lists the configurations the
integrator runs and those that stop at set-up.

!!! info "Going deeper"
    - [Multirate: Method in Detail](multirate-method.md): the construction of the levels and of $H$, the checks of the local
      steps, the non-local CFL condition, the recursion, boundaries, turbulence models and MPI.
    - [Multirate Analysis](multirate-analysis.md): the proofs of the order, stability and monotonicity.
