# Changelog

All notable changes to MOSE are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- The key `diffusive-metric` of `[MOSE-Numerics]` chooses the face metric of the viscous fluxes at interior faces:
  `inverse-mean`, the default, the inverse of the mean of the inverse metrics of the two cells, whose face-normal edge is the
  vector between the cell centres; or `mean`, the mean of the two metrics, as before.
- Probes can record `mi_tilde`, the first turbulence variable, which the solution files name `rhoNuTilde`
  (Spalart–Allmaras), `rhoK` (two-equation models) or `ru'u'` (Reynolds-stress model).

### Changed

- The interior faces of the viscous fluxes take by default the inverse of the mean of the inverse metrics. The mean of the two
  metrics divided the difference across a face by the mean of 1/h instead of the distance between the cell centres: it
  overstated the face-normal gradient, and with it the viscous, conductive and diffusive fluxes, by (1 + r)²/(4r) at a cell-size
  ratio r (4.2 % at r = 1.5, 8.9 % at r = 1.8), at any resolution. Viscous runs change wherever neighbouring cells differ in
  size; `diffusive-metric = mean` keeps the previous metric.
- Set-up errors stop the run before the initial condition is read, on every MPI rank: the root rank found them only after
  reading the initial condition, and stopped by itself. Errors of the transport data are reported there, such as a mixture of
  more than one species with `schmidt` = 0, its default, without `INPUT/diffusion.dat`, which failed later, in the first step of
  the run, when the diffusive flux read the table that was not there; a single species needs no such table, as it has no species
  diffusion.
- A restart appends to the probe files, and creates those that are missing; it stopped when one was missing.
- Each rank keeps the state of the remote chimera donors of its own receivers only, and none of the remote sources of manifold
  faces (`501`), which it never reads; the chimera exchange sends each donor cell once per rank.
- The diagnostic output of a steady run with a single species needs `sol-variables` to name `thermo` (the default) or
  `transport`; otherwise the run stops at its first solution output with a message that says so, where it named only the
  size of a buffer.

### Fixed

- Faces between blocks and chimera faces take the inverse of the mean of the inverse metrics: across a large size jump at a block
  connection the mean of the two metrics made the normal gradient far too steep, and the explicit diffusion unstable on the large
  cell.
- With MPI the tangential neighbours of the source cell of a block connection, which on the edge of a face can be ghost cells
  filled from another rank, are taken after the ghost exchange: they were one stage old, and MPI runs drifted from the serial run
  at block edges and junctions.
- With MPI the output variables of a block are computed on the rank that owns it: on the root rank the eddy viscosity of a
  Spalart–Allmaras run with `sol-variables = transport` read the metrics of the blocks of the other ranks, which the root rank
  does not allocate, and the velocity gradients at block edges read ghost cells that were never filled.
- At chimera faces (`102`) the viscous flux takes the tangential gradient of the donor side, and in RANS runs the chimera ghost
  cells take the wall distance and roughness of the receiver instead of zero.
- With a turbulence model the diagnostic output writes the residuals of the turbulence variables and declares as many variables
  as it writes: it declared nres + 2 variables and named and filled seven.
- An unknown probe variable, or a probe file that cannot be opened, stops every rank with a message: the variable was left
  undefined, and the run wrote to a unit never opened.
- A degenerate cell (a zero determinant of the metric) and a rejected interface state stop the run with a message that locates
  them: the corner nodes of the cell; the block, the face, the two states of the interface and its stencil. The degenerate cell
  was announced with a fallback that never ran.
- The diagnostic output of a steady run names and fills the soot residuals (`rho*Yp`, `rho*Np`) between ρE and the turbulence
  variables: with a soot model it named fewer variables than it wrote, and with a turbulence model too it stopped on a run-time
  error in builds with bounds checking.
