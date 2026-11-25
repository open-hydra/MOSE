# Multigrid and FAS cycles

MOSE can use multiple nested grids to obtain an initial solution on the finest grid more efficiently. It supports **grid sequencing** and, for supported steady cases, a **Full Approximation Scheme (FAS) cycle**.

## Grid sequencing

Set `levels` in `[MOSE-Multigrid]` to the number of grid levels. The default is one level (no multigrid sequence). MOSE starts on the coarsest level and proceeds to the finest.

```ini
[MOSE-Multigrid]
levels = 3
level3-iter = 2000
level2-iter = 2000
level1-iter = 100000
```

`level<n>-iter` sets the iteration limit for level `n`. These values are read again during the run; lowering a level's limit below the iterations already completed ends that level after the next runtime read.

Each coarser level halves the cell count in each direction that has more than one cell. The finest-grid cell counts must therefore be divisible by `2^(levels-1)` in `i`, and in `j` and `k` when those directions have more than one cell. MOSE checks this during setup.

Provide a boundary-condition file for each grid level:

| File | Grid |
|---|---|
| `INPUT/bc.txt` | Level 1, finest grid |
| `INPUT/bc2.txt` | Level 2 |
| `INPUT/bc3.txt` | Level 3 |
| `INPUT/bc4.txt` | Level 4 |

## Grid sequencing or FAS?

With `cycle = none` (the default), coarse grids are used only to generate a starting solution for the next finer grid. This is grid sequencing.

With `cycle = v` or `cycle = w`, MOSE visits coarser levels repeatedly while solving the current level. A V-cycle visits each coarser level once; a W-cycle visits the coarser levels twice.

```ini
[MOSE-Multigrid]
levels = 3
cycle = w
```

FAS cycles are intended for steady runs without finite-rate or equilibrium chemistry. They are not a generic option for time-accurate simulations.

## FAS parameters

| Parameter | Default | Purpose |
|---|---:|---|
| `cycle` | `none` | `none`, `v`, or `w`. |
| `cycle-start` | `0` | Iteration of the current level at which the cycle begins. Before this, the level runs as with grid sequencing. |
| `pre-iter` | `2` | Smoothing iterations before visiting a coarser level. |
| `post-iter` | `0` | Smoothing iterations after coarse-grid correction. `pre-iter + post-iter` must be at least one. |
| `coarse-iter` | `1` | Smoothing iterations on the coarsest level in a cycle. |
| `damp-restriction` | `1.0` | Factor applied to the restricted residual in the forcing term. |
| `damp-prolongation` | `0.5` | Factor applied to the prolongated correction. |
| `cycle-turbulence` | `.false.` | Also solve supported one- or two-equation RANS equations on coarse levels; otherwise turbulence values stay frozen at their restricted state. |
| `cycle-turbulence-start` | `0` | Delay, in iterations after `cycle-start`, before turbulence equations join the cycle. |

## Restrictions checked by MOSE

The FAS cycle is rejected during setup when any of these apply:

- the run is time-accurate;
- the chemistry model is not `frozen`;
- there is only one grid level;
- `pre-iter + post-iter` is less than one;
- `cycle-turbulence = true` is combined with the SSG/LRR Reynolds-stress model.

The cycle's working arrays are allocated at startup. A run started with `cycle = none` cannot enable a V- or W-cycle later by editing `input.ini`; the runtime read is rejected. Turning an active cycle off by setting `cycle = none` is allowed.

## Practical tuning

Start with the defaults. Change one setting at a time and compare convergence with grid sequencing. If the residual stalls above the grid-sequencing result, return to the defaults before trying another change.

The project guidance recommends keeping `coarse-iter = 1`; using `damp-prolongation = 1` can accelerate a case but can also make the cycle stall or stop. With Spalart–Allmaras, test `cycle-turbulence = true`; for the k–omega models, keep it off unless a case-specific validation supports it. In general, avoid implicit residual smoothing with the cycle (`irs-beta = 0`).

For the rationale and validated settings, see FAS-cycle verification and the [multigrid theory](../theory/time-integration.md#multigrid-acceleration).

## Restart and multigrid

A normal completed run writes its final restart solution on the finest grid. A restart from that solution continues on the finest level; it does not repeat the coarse-grid sequence. Intermediate files such as `field-level<N>` or `field-prolongated` are diagnostic/output products, not the usual restart checkpoint. See [Restart a simulation](restart.md) and [Restarts in detail](restart.md#limits).
