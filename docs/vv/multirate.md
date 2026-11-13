# Multirate Verification

Verification of the multirate (local time stepping) integrator described in
[Time Integration](../theory/time-integration.md#multirate-time-stepping).  The checks live in
`test/multirate/mr-*`, build their own small cases with `test/multirate/mr-common/gen_mr_cases.py` (perfect-gas air
of `test/common/Air`, Euler equations, HLLC, MUSCL with the van Leer limiter; the laminar Navier–Stokes cases use the
viscous air of `test/common/Air-viscous` or a viscosity of their own, the Spalart–Allmaras and SST cases the air of
`test/common/Air-viscous`) and need only numpy.  They are registered with CTest, one test per check, all with the label `multirate`: the test is named after its check,
`Multirate.OneLevel` runs `mr-onelevel`, `Multirate.Cfl4` runs `mr-cfl4`, and so on (`Multirate.Smoke` is `mr-smoke`).
`mr-order`, `mr-deep` and `mr-daru` have a short form (`Multirate.Order`, `Multirate.Deep`, `Multirate.Daru`) and
longer forms, heavy tests that are run by hand (`Multirate.OrderFull`, `Multirate.DeepFull` with `MR_FORM=full`;
`Multirate.DaruMedium`, `./run.sh medium`; `Multirate.DaruFull`, `./run.sh full`). `Multirate.Hybrid` (`mr-hybrid`) and
`Multirate.RestartHybrid` (`mr-restart-hybrid`) need an MPI build. The tier of each test, light, medium or heavy, is
in [`test/tiers.txt`](https://github.com/open-hydra/MOSE/blob/main/test/tiers.txt) and in the
[test plan](https://github.com/open-hydra/MOSE/blob/main/test/TESTING-PLAN.md).

```bash
ctest --test-dir build -L multirate -LE heavy --output-on-failure   # the multirate tests that a pull request runs
ctest --test-dir build -L multirate --output-on-failure             # all of them, the heavy ones included
ctest --test-dir build -R Multirate.Cfl4 -V                         # one check
cd test/multirate/mr-cfl4 && ./run.sh                               # the same, from the source tree
```

## Cases

| Case | Mesh | State | Levels with safety 2 |
|---|---|---|---|
| `mr1d-acu` | 1D periodic, 284 cells in bands of 10 / 4 / 1.6 / 0.64 / 1.6 / 4 / 10 mm | right-running acoustic wave (amplitude $10^{-3}$) on $u_0 = 50$ m/s | 0, 1, 2, 3 (`max-level = 3`) |
| `mr1d-ent` | as `mr1d-acu` | entropy wave $\rho_0(1 + 0.1\sin)$, $u = 100$ m/s, uniform $p$ | 0, 1, 2, 3 |
| `mr1d-sod` | as `mr1d-acu`, closed by slip walls | Sod shock tube, discontinuity at $x = 0.25$ m | recomputed at every macro step |
| `mr1d-mix` | as `mr1d-acu` | frozen mixture of four species (CH4, H2O, N2, O2): the oxidizer mixture of `test/1D/diff1` with a slug of its fuel mixture for $0.2 < x < 0.3$ m, laminar with a Schmidt number of 0.7, 300 K, $10^5$ Pa, $u = 100$ m/s: two contact discontinuities | 0, 1, 2, 3 |
| `mr1d-blast` | 1D, 400 uniform cells of 2.5 mm, closed by slip walls | hot gas ($64\,T_0$, $5\,p_0$, $T_0$ = 100 K) for $x < 0.1$ m in cold gas at rest: a blast wave whose front enters coarser cells (Gnedin, Semenov and Kravtsov, section 2) | 3 in the hot gas, 0 in the cold gas, recomputed at every macro step |
| `mr1d-deep` | 1D periodic, bands of 4 cells with step ratio 2.1 | acoustic wave | 0 … 8 (`max-level = 8`) |
| `mr2d-box` | 2D closed box, 66 × 66 cells, refined cross | Gaussian pressure pulse, CFL 0.2 | 0, 1, 2, 3 |
| `mr2d-free` | 2D periodic, curvilinear, stretched and skewed, 48 × 48 | uniform flow ($u = 100$, $v = 50$ m/s), CFL 0.1 | 0, 1, 2 (`max-level = 2`) |
| `mr1d-deep1` | as `mr1d-deep`, bands of one cell | acoustic wave | 0 … 8; the halo of a level $l$ holds cells of levels $\le l-2$ |
| `mr1d-deep12`, `mr1d-deep12-1` | 13 levels: bands of 4 cells (of one cell), finest cell 1.36 µm | acoustic wave | 0 … 12 (`max-level = 12`) |
| `mr1d-deep12-3b` | `mr1d-deep12` in three blocks, connections on the edge of two bands and inside the finest band | acoustic wave | 0 … 12 |
| `mr1d-uni-3b` | the mesh and the blocks of `mr1d-acu-3b` | a uniform flow, 1.2 kg/m³, $10^5$ Pa, 50 m/s | as `mr1d-acu-3b` |
| `mr1d-acu-3b` | `mr1d-acu` in three blocks, connections inside a band and on the edge of two bands | as `mr1d-acu` | 0 … 3; the cell next to the edge connection is raised to the level of its partner |
| `mr1d-acu-3bi` | `mr1d-acu` in three blocks, both connections inside bands | as `mr1d-acu` | the map of `mr1d-acu` |
| `mr2d-box-4b` | `mr2d-box` in 2 × 2 blocks (vertex of four blocks in the refined cross) | as `mr2d-box` | 0 … 3 |
| `mr2d-box-4bx` | `mr2d-box` in 2 × 2 blocks cut one cell from the edge of the finest band: with the blocks on different MPI processes a face between levels is one cell from a process boundary | as `mr2d-box` | 0 … 3 |
| `mr2d-wedge-visc` | the field with dilatation and swirl of part 4 of the test `MMSResiduals` on the 2-D axisymmetric wedge (1°, faces of type 200, axis of type 300), 32 × 32 cells | laminar, `Air-viscous` | 0 |
| `mr2d-wedge-band` | the 2-D axisymmetric wedge (1°, faces of type 200, axis of type 300) with its rows in $r$ in bands of 8 rows of 0.125, 0.25, 0.5 and 1 mm from the axis, 24 columns of 1 mm, the other faces slip walls | rest air of `Air-viscous` with a solid-body swirl and a Gaussian pressure pulse; laminar, CFL 0.5, VNN 0.3 | 2, 2, 1, 0 by bands from the axis (`max-level = 2`): faces between levels on cells with records 200 |
| `mr2d-box-sym`, `mr2d-box-sym-4b` | as `mr2d-box`, one and 2 × 2 blocks | pressure pulse at the centre of symmetry of the mesh | 0 … 3 |
| `mr2d-free-4b` | `mr2d-free` in 2 × 2 blocks | uniform flow, CFL 0.1 | 0, 1, 2 |
| `mr2d-vortex-3b` | 2D periodic, three blocks in $x$, the central one 4 times finer in $x$ (1.25 mm against 5 mm) | isentropic vortex on $u = 100$ m/s, CFL 0.2 | 0, 1 (`max-level = 1`): the central block and the columns next to the connections at level 1 |
| `mr2d-shock` | as `mr2d-box`, closed by slip walls | planar Riemann problem at $x = 0.2$ m (Sod-like), CFL 0.2 | 0 … 3, fixed or recomputed |
| `mr1d-402`, `mr1d-410` | 1D duct, 0.64 mm cells at the inlet, back-pressure outlet (`406`) | inlet `402` with $p_0(t)$, or `410` with a tabulated acoustic wave | 3 at the inlet … 0 |
| `mr1d-401`, `mr1d-403`, `mr1d-405`, `mr1d-408`, `mr1d-420` | the duct of `mr1d-402` | inlet `401` (total conditions), `403` (total temperature and mass flux), `405` (supersonic, Mach 1.5), `408` (static temperature and velocity) or `420` (total conditions and mass flux), with the data of the steady state | 3 at the inlet … 0 |
| `mr1d-502` | the duct of `mr1d-402` | an SRM grain in mode 2 at the inlet (an inflow: the mass-flux inlet `404` with the mass flux of Saint-Robert's law and the static temperature of total `T_af` at the pressure of the boundary cell), `mdot_ref` the mass flux of the steady state, `n` = 0.5, `p_ref` its pressure, `T_af` its total temperature | 3 at the inlet … 0 |
| `mr1d-502c` | the duct of `mr1d-402` | the grain of `mr1d-502` with `n` = 0 (its mass flux `mdot_ref` at every pressure), the duct closed at its end by a slip wall (`300`), the gas at rest at the start | 3 at the grain … 0 |
| `mr2d-plate` | laminar flat plate of `test/2D/viscous/flat-plate-laminar`: 64 × 64 cells down to 16 µm at the wall, inlet `404`, ambient `407`, symmetry and adiabatic wall (`301`) | uniform start at Mach 0.2, air of `test/common/Air-viscous`, CFL and VNN 0.5 | 0 … 3 by rows (`max-level = 3`) |
| `mr2d-plate-2b`, `mr2d-plate-4b` | the plate in two blocks (connection across the boundary layer) and in 2 × 2 blocks (vertex of four blocks in the boundary layer) | as `mr2d-plate` | as `mr2d-plate` |
| `mr2d-plate-302` | `mr2d-plate` with its wall isothermal (`302`) at 330 K | as `mr2d-plate` | as `mr2d-plate` |
| `mr2d-plate-2br` | `mr2d-plate-2b` with block 2 turned by 180 degrees about $z$ (its $i$ and $j$ reversed): the connection across the boundary layer joins face 4 of both blocks with the orientation `-1 0 0 1` | as `mr2d-plate` | as `mr2d-plate` |
| `mr2d-hot`, `mr2d-hot-4b` | 2D periodic parallelogram, 60 × 60 skewed cells of 5 mm, one and 2 × 2 blocks | hot square ($8\,T_0$, smooth edges) carried by $u = 100$, $v = 50$ m/s, viscosity 0.3 Pa s, CFL 0.2 | 0, 1 (`max-level = 1`): a square of level 1 with sharp corners |
| `mr3d-hot`, `mr3d-hot-8b` | 3D periodic, 12 × 12 × 12 cells of 10 mm on an oblique affine mesh, one block and 2 × 2 × 2 blocks | hot cube ($8\,T_0$) carried by $u = 100$, $v = 50$, $w = 30$ m/s, viscosity 0.3 Pa s, CFL 0.2 | 0, 1: a cube of level 1 with sharp edges and vertices |
| `mr2d-chan-sa`, `mr2d-chan-sst` | wall channel of `test/checks/time-integration/turbulence-sources`, periodic along the flow (4 cells of 1 mm), adiabatic wall, symmetry plane at 0.553 mm, 16 rows from 4 µm growing by 1.25 | $u = U\tanh(y/0.2\,\text{mm})$ with no driving force; Spalart–Allmaras with $\tilde\nu$ twice the near-wall equilibrium, or SST with near-wall $k$ and $\omega$ (wall condition `asymptotic`); CFL 0.8, VNN 0.4 | 0 … 3 by rows |
| `mr2d-chanrough-sa` | as `mr2d-chan-sa`, rough wall (sand grain 10 µm) | as `mr2d-chan-sa`, SA-rough | 0 … 3 by rows |
| `mr2d-chanpos-sa`, `mr2d-chanpos-sst` | the positivity channel of `test/checks/time-integration/turbulence-sources`: 2 µm at the wall, 34 rows | near-wall turbulence far above its equilibrium ($\tilde\nu$ ×100; $\omega$ ×10 with the wall condition `practical`): the destruction acts faster than the transport | 0 … 3, recomputed at every macro step |
| `mr2d-plate-sa`, `mr2d-plate-sst`, `mr2d-plate-sa-4b`, `mr2d-plate-sst-4b` | turbulent flat plate of `test/multirate/mr-common/plate` on its coarse mesh (68 × 48 cells, 4 µm at the wall), one and 2 × 2 blocks | uniform start at Mach 0.2 with Spalart–Allmaras or SST; CFL 0.5, VNN 0.4 | by rows |

The 1D and the cartesian 2D cases are slabs one cell thick between slip walls, 10 mm thick, so that every step of the
tests stays below the slab thickness divided by the speed of sound: above it the wall mode across the slab is unstable
with Heun: the Riemann fluxes of the two walls, on a small velocity $w$ across the slab, give pressures $p \pm \rho a w$ and
$dw/dt = -2aw/\Delta z$, and Heun is stable on the negative real axis for $|\lambda h| \le 2$, that is $h \le \Delta z/a$
(the curvilinear `mr2d-free` cases keep 1 mm, where the steps at CFL 0.1 stay below the limit).

## Tests and acceptance criteria

| Check | What it checks | Criterion |
|---|---|---|
| `mr-smoke` | the light test: `mr1d-acu` with one level against the single rate (20 steps); `mr1d-sod` with four levels (`max-level` 3, `safety` 2, map rebuilt at every macro step), 20 macro steps in one go and as 10 + restart (`vtk raw`) + 10; `mr1d-acu` with cfl 1.2 and `max-level` 3 (`nonlocal-cfl = redo`, the default); `mr1d-sod` with one level and cfl 20 (a run that diverges); a 1-D tube of 8 cells with 400 ends and the gas entering through face 1 at 0.5 a (the case `count` of `test/checks/boundary-conditions/extrapolation-400`), 20 Heun steps at first order, with one level and with the single rate | fields and residual history identical byte by byte; mass and energy drift $\le 10^{-12}$ with the four levels in use; output trees identical byte by byte; non-zero exit code, the message `[MOSE-Multirate] cfl` and no solution; exit code 1, the message of the state check of the multirate (a value that is not finite or a negative pressure) and no solution; the lines of the gas entering through the 400 faces and the end line of the run with one level are those of the single rate, word for word: one line (iteration 1, faces 1) and 40 flux evaluations, one a stage |
| `mr-onelevel` | `max-level = 0`, `safety = 1` against the single-rate step (`mr1d-acu`, `mr1d-sod`, `mr2d-box`, `mr1d-acu-3b`, `mr2d-box-4b`, and the ducts `mr1d-402`, `mr1d-410` with their data at the stage times, and `mr1d-401`, `mr1d-403`, `mr1d-405`, `mr1d-408`, `mr1d-420` and the SRM grain in its inflow mode, `mr1d-502`, so that every boundary type that the integrator accepts runs); laminar: `mr2d-plate`, `mr2d-plate-2b` and `mr2d-plate-4b` with `diffusive-metric = mean`, `mr2d-plate` and `mr2d-plate-4b` also with `stencil` (the default), `mr2d-hot-4b`, `mr3d-hot`; Spalart–Allmaras and SST: `mr2d-chan-sa`, `mr2d-chan-sst`, `mr2d-chanrough-sa`, `mr2d-plate-sa`, `mr2d-plate-sst`, `mr2d-plate-sa-4b`, and the wall channels with `SA-RC` and with `SST-RC`; the axisymmetric wedge with viscous stress on its faces, `mr2d-wedge-visc`; a 1-D slab between symmetry j and k faces, planes normal to y and z, whose i faces are inclined, `mr1d-skew` (y and z momentum frozen, with the residual the pressure on the inclined faces gives them), and a 1-D duct whose top j and k faces are inclined, `mr1d-nozzle` (not planes normal to y and z: y and z momentum integrated); the plate with the records of `bc.txt` in reverse order, `mr2d-plate-rev`, and with an isothermal wall, `mr2d-plate-302`; the per-cell copies of the turbulence sources, of the local step (`Dtcap_Cell`, wedge term and bound of the stiffness of the turbulence sources included) and of the state check, the select of the ghost cells of a record (`Ghost_Record` against `Fill_Ghost_Cell`, every branch with the calls under a condition: `Ghost_Wall` at the walls and the 400 faces), the select of the boundary fluxes of a record (`BC_Flux_Record` against the branches of `BC_Fluxes` for the types that the integrator accepts), the copy of the loop of the connection fluxes (one flux per face) and the face metrics (the stored metric read in place, or the computed one) and call of the viscous flux of each direction, used with more than one level, against the single-rate routines statement by statement; the MPI tags of the multirate against those of the single-rate exchanges | fields, wall output and residual history identical byte by byte; the set-up lines that say that the y and z velocity of `mr1d-skew` are frozen and those of `mr1d-nozzle` integrated; copies identical to their originals; no MPI tag shared |
| `mr-conservation` | mass, momentum (periodic cases) and total energy after 400 macro steps (`mr1d-acu`, `mr1d-ent`, `mr1d-sod`, `mr2d-box`, `mr1d-acu-3b`, `mr2d-box-4b`; laminar: `mr2d-hot`, `mr2d-hot-4b`, `mr3d-hot` with the levels recomputed every 5 macro steps; Spalart–Allmaras and SST: `mr2d-chan-sa`, `mr2d-chan-sst`, and Spalart–Allmaras on a rough wall: `mr2d-chanrough-sa`, mass and energy); `diffusive-metric = mean` with levels (`mr2d-chan-sa`, `max-level` 3, mass and energy; `mr3d-hot`, `max-level` 1, mass, momentum and energy); levels recomputed every 5 macro steps; a `macro-step-target` below the adaptive step; positive control with `conservative-correction = .false.`; the axisymmetric wedge with levels (`mr2d-wedge-band`, `max-level` 2, mass and energy) with its own positive control; `mr1d-blast` with `max-level` 10, 200 macro steps, the depth of the map following the octaves of the local steps; `mr1d-mix`, `max-level` 3, the levels recomputed every 5 macro steps: the mass of each species, the total mass and the momentum, with its own positive control; with `max-level` 3 and the levels recomputed every 5 macro steps, `SA-R` and `sst-production = compressible` on the wall channels (mass and energy), `k-energy-coupling` on the SST channel (mass; with $k$ in the total energy the sources of $k$ enter the energy, whose drift is reported), each against the run of the same channel without the option, and `safety` 1.5 on `mr1d-acu-3b` (mass, momentum and energy); the SRM grain in its inflow mode with levels, `mr1d-502c` (`max-level` 3), the mass of the domain against $M(0) + \dot m_\mathrm{ref} A t$, $A$ the area of the grain, with its own positive control | relative drift $\le 10^{-12}$ (for `mr1d-502c` the relative difference to $M(0) + \dot m_\mathrm{ref} A t$); the runs with `SA-R`, `sst-production = compressible` and `k-energy-coupling` not identical to the runs without the option; the controls must drift by $\ge 10^{-9}$; with `mr1d-blast` at least two depths printed and cells at level 0 in every macro step |
| `mr-order` | temporal order with levels fixed from $H_0 = 2^L\,\Delta t_{\min}/\sigma$; short form: macro steps $H_0$ and $H_0/2$, reference single rate with $\Delta t_{\min}/16$; full form: macro steps $H_0 \dots H_0/8$, reference with $\Delta t_{\min}/64$; the level map (for the interface cells) is read from the log of the run with $H_0$ (`trace = .true.`); `mr1d-acu`, `mr1d-ent` and `mr2d-vortex-3b` (whose interface cells include those next to the connections); laminar: `mr2d-plate` from the uniform start ($L = 3$), `mr2d-hot` ($L = 1$) and `mr2d-wall-s`, a no-slip wall on smoothly skewed cells (skew angle at the wall 14 to 37 degrees) whose levels change along the wall ($L = 2$); Spalart–Allmaras and SST: `mr2d-chan-sa`, `mr2d-chan-sst` ($L = 3$; $\rho$, $u$, $v$, $p$ and $\tilde\nu$, $k$, $\omega$) | last pair of the form: $L^1 \ge 1.9$ ($\ge 1.8$ for the turbulent channels); $L^\infty$ within 4 cells of the faces between levels (floor 5 % of the global $L^\infty$) $\ge 1.8$; for the entropy wave $u$ and $p$ stay uniform to $10^{-12}$; for the hot square and the skewed wall also $L^\infty$ on the cells next to the faces between levels (same floor) $\ge 1.8$ |
| `mr-sod` | Sod tube up to $t = 1$ ms with the levels recomputed at every macro step, against the exact Riemann solution and the single rate on $x \ge 0.15$ m | conservation $\le 10^{-12}$; $L^1(\rho - \rho_{\text{exact}}) \le 1.05\times$ the single rate's; shock within one local cell; no more local extrema of $\rho$ than the single rate |
| `mr-deep` | 9 levels (`mr1d-deep`, `mr1d-deep1`) and 13 levels in three blocks (`mr1d-deep12-3b`): order (as `mr-order`, $T = 16 H_0$), conservation after 200 macro steps, stability over 2000 macro steps; the short form leaves the order and the stability with 13 levels to the full form | as `mr-order`; drift $\le 10^{-12}$; final amplitude $\le 1.01\times$ the initial |
| `mr-threads` | 1 against 3 and 4 OpenMP threads (`mr1d-sod` with recomputed levels, `mr2d-box`, `mr1d-acu-3b`, `mr2d-box-4b`; laminar: `mr2d-plate-4b`, `mr2d-hot-4b`, `mr3d-hot`; Spalart–Allmaras and SST: `mr2d-plate-sa-4b`, `mr2d-chan-sst`) | identical byte by byte |
| `mr-trace` | sequence of the 30 stages of one macro step with levels 0 … 3 (level, stage, abscissa, and the fraction $\eta$ of the step of every coarser level, `theta` in the trace) | equal to the reference sequence of the recursion in `bench_sequence_L3.txt` |
| `mr-halo` | halo cells and ghost records of every level (`mr1d-acu`, `mr1d-deep`, `mr2d-box`, `mr1d-deep1`, `mr1d-deep12-1`, `mr1d-acu-3b`, `mr1d-deep12-3b`, `mr2d-box-4b`; laminar, with the diagonals of the diffusive faces and the tangential neighbours of the partner cells of the connections: `mr2d-plate-4b`, `mr2d-hot`, `mr2d-hot-4b`, `mr3d-hot`; Spalart–Allmaras and SST, with the stencil of the turbulence sources: `mr2d-plate-sa-4b`, `mr2d-chan-sst`; the axisymmetric wedge with levels, records 200 on faces between levels: `mr2d-wedge-band`), printed by MOSE in the log with `trace = .true.` (lines `MR-LEVEL`, `MR-HALO`); in the viscous cases, at every stage of the first macro step, the viscous variables that the viscous fluxes of the stage read (the ten cells of every diffusive face of the level, the cell and its face neighbours of every connection and wedge record of the level) against those of their primitive variables at that point (lines `MR-VISC`) | identical to an independent construction from the stencils; with bands of one cell, levels $\le l-2$ in the halo of every level $l \ge 2$; on the hot square and the hot cube, halo cells reached only through the diagonals; viscous cases: one `MR-VISC` line per stage, every value read equal to the one of its primitive variables (none stale) |
| `mr-free` | uniform flow on the curvilinear periodic mesh, one block and 2 × 2 blocks, 100 macro steps; informative run at CFL 0.8 | deviation $\le 10^{-13}$ |
| `mr-blocks` | levels across the connections (`mr1d-acu-3b`, `mr1d-deep12-3b`, `mr2d-box-4b`); one block against three blocks with the same map (`mr1d-acu`, `mr1d-acu-3bi`), 100 macro steps | equal levels on the two sides of every connection, jumps $\le 1$; difference $\le 10^{-11}$ and $\le 10\times$ that of the single rate |
| `mr-bc` | ducts `mr1d-402`, `mr1d-410`: temporal order as `mr-order` at the four inlet cells, on the finest level; the order compares runs that share the data, so the first cell is also compared with the table at the same time (its total pressure with `402`, its pressure with `410`), at the peak of the data (4 macro steps, about a quarter of the period) and at the end of the run | last pair of the inlet $L^\infty$ (floor 5 % of the global $L^\infty$) $\ge 1.8$ for $\rho$, $u$, $p$; $u > 0$ in the four inlet cells; difference from the table at the peak at most half of the distance the table has moved from $t = 0$, at the end at most the amplitude of the table |
| `mr-moving-levels` | a simple wave of the Euler equations whose faces between levels travel with it (`mr1d-simple-<n>`: 1D periodic on $[0, 1]$, $c = 1 + 0.05\,\varphi$, $\varphi = (1 - \cos 2\pi x)/2$, $u = 2(c - 1)/(\gamma - 1)$, non-dimensional gas of `test/common/nondim-Air`, up to $T = 0.5$, before the wave breaks at $t = 1.06$), $n = 256$, 512, 1024 and 2048 cells with $\Delta x$ and $H$ halved together; multirate with `max-level` 1, `safety` 2, the map rebuilt at every macro step and $H = 0.565\,dt_0$ ($dt_0$ the single-rate step at $t = 0$), so that the faces between levels lie where $u + c$ is halfway between its extremes and move with the wave; the same with the map of the first macro step kept, and the single rate with the step $H$ and with $dt_{\min}$ | against the exact solution by characteristics: $L^1$ order of $\rho$, $u$ and $p$ on the last pair $\ge 1.9$ with the map rebuilt; both levels in the map of the first and of the last macro step of every multirate run; mass drift $\le 10^{-12}$; the $L^\infty$ orders reported, global and on the cells next to the moving faces |
| `mr-shock` | `mr2d-shock` up to $t = 1$ ms, levels fixed and recomputed, against the exact Riemann solution and the single rate on $x \ge 0.25$ m | conservation $\le 10^{-12}$; $L^1(\rho - \rho_{\text{exact}}) \le 1.05\times$ the single rate's; no more local extrema of $\rho$ on any row |
| `mr-sym` | `mr2d-box-sym`, `mr2d-box-sym-4b`, multirate and single rate, 100 macro steps | mirror symmetry in $x$, $y$ and about the diagonal: every run within $10^{-10}$ of the excursion of each variable, and the multirate within 10 times the asymmetry of the single rate on the same case (round-off floor: the single rate alone reaches $3$–$4\cdot10^{-12}$ in 100 steps) |
| `mr-reject` | 35 configurations outside the scope of the integrator, one run each (the shock detector implied by the Riemann solver or set by `shock-detector`; a boundary type without a multirate flux, the types `409` and `501` stopping every run earlier, at the reading of `bc.txt`, and the walls `505` and `506` at set-up, on the species of their surface reactions; `SAcomp`, whose compressibility correction the stages do not evaluate, `Wilcox2006` on the wall channel, and `SA-RC` and `SST-RC` with `max-level` 3; an SRM grain in mode 1; a restart with `levels = fixed` or with `fixed-macro-step`; with `max-level` > 0 and `nonlocal-cfl = redo`, `levels = recompute` and no `fixed-macro-step`, cfl above 1, with van Leer, first-order and minmod, and vnn above 0.5 in a Navier–Stokes run; `max-attempts` 0, below the rule `>= 1` of the key registry); control runs, Euler, laminar and Spalart–Allmaras, vnn 0.5 in a laminar run and 0.6 there with `max-level` 0 and in an Euler run, cfl 1.2 and, in a laminar run, vnn 0.6 with `nonlocal-cfl = off` and `count`, cfl 1.2 with `nonlocal-cfl = redo` and `levels = fixed` or a `fixed-macro-step`, and `levels = fixed` with `fixed-macro-step` above the step of the map, with `max-level` 3 and with `max-level` 10 on `mr1d-acu`, whose map has depth 3; `input.ini` read again during the run on the laminar plate (`mr2d-plate`, iteration 30 of 40): `max-level` changed with `ini-error` left out and with `stop`, `enabled = .true.` written during a single-rate run, with `nonlocal-cfl = redo` a `cfl` above 1 with `ini-error` left out and with `stop`, and with first-order, and a `vnn` above 0.5, with `nonlocal-cfl = off` a `cfl` above 1 and a `vnn` above 0.5, `time-accurate = .false.`, `cfl-rise-threshold` 1000, a `cfl` within the bound (above the floor of $\kappa$ and at it), and on the wall channel with SST (iteration 800 of 1000) `k-energy-coupling = .true.` | non-zero exit code and the message `[MOSE-Multirate] <item>` (for `max-attempts` 0 the message of the key registry); the controls complete, those with `levels = fixed` above the step of the map with the warning of the cells beyond their limit once and its run line; the reading that changes a key of `[MOSE-Multirate]` or a key fixed for every run (`time-accurate`, `k-energy-coupling`), or takes `cfl` and `vnn` (with `redo`) or `cfl-rise-threshold` out of the scope of the integrator, is not taken (exit code 0 at the last iteration, the warning that names the key, no `[ERROR]`, the residual history of the run without the change byte for byte) or, with `ini-error = stop`, stops the run at iteration 30 with the message; with `off` the readings of a `cfl` and a `vnn` above the bounds are taken (exit code 0 at the last iteration, no warning); the reading of a `cfl` within the bound is taken, and the run equals 30 macro steps and a restart with that `cfl` byte for byte (the interface margin follows the `cfl`) |
| `mr-positivity` | `mr2d-chanpos-sa`, `mr2d-chanpos-sst`, levels recomputed at every macro step, `max-level` 3 with `safety` 2 and `max-level` 0 with `safety` 1 (the stiffest cells at their limit), 200 macro steps; positive control: both cases with `max-level` 3, `safety` 2 and `fixed-macro-step` 48, whose finest cells run at 12 times their limit | the runs complete; no stage makes $\tilde\nu$, $k$ or $\omega$ negative (no count of negative values in the solver output, which prints it only when a stage made one negative); the local step was limited by the turbulence sources (smallest step below half the transport step of the initial state); the final field holds the turbulence variables (`rhoNuTilde`; `rhoK`, `rhoOmega`); the positive controls complete and print the count |
| `mr-hybrid` | `mr2d-plate-4b` (laminar plate in 2 × 2 blocks, a vertex of four blocks in the boundary layer), `mr2d-plate-sa-4b` (Spalart–Allmaras) and `mr2d-plate-sst-4b` (SST), `mr3d-hot-8b` (the hot cube in 2 × 2 × 2 blocks, a vertex of eight blocks in the cells of level 1), `mr2d-box-4b`, `mr1d-acu-3b`, `mr2d-box-4bx` and `mr2d-plate-2br` (the connection with the orientation `-1 0 0 1`), single rate 300 steps and multirate (`safety` 2, `max-level` 3) 100 macro steps, `mr2d-plate-sa-4b` with `SA-RC` and `max-level` 0 on 1 and 2 processes, and `mr1d-acu-3b` with `max-level` 10 (local steps over 3.97 octaves, with cells of the slowest octave in every block: depth 3) and `mr1d-deep12-3b` with `max-level` left out, 20 macro steps (local steps over 12.85 octaves, those of its middle block over 8.56: depth 10 from the largest local step of all the processes), on 3 processes against 1, MPI processes × OpenMP threads 2×1, 4×1, 1×2, 1×4, 2×2, 4×2 (SA, SST and `mr3d-hot-8b`: 2×1, 4×1, 2×2; `mr2d-plate-2br`: 2×1, 2×2; the others 2 and 4 or 3 processes) against 1×1; the 1×1 runs of `mr2d-plate-2br` against those of `mr2d-plate-2b`; `mr1d-blast100-3b` (the strong blast wave in three blocks) with `nonlocal-cfl = redo` and `force-redo` on 3 processes against 1; `mr1d-blast-wide-3b` (the blast tube with far coarse cells, in three blocks) with `max-level` and `nonlocal-cfl` left out (10 levels, `redo`) up to 0.6 ms, 3×1 and 3×2 against 1×1; `mr1d-blast-hot-3b` (the hot gas at 9990 K in three blocks, the contact in block 2) with `redo` on 3 processes; `mr1d-blast-wide-3b` with `max-attempts` 3, 1×1, 3×1 and 3×2; `mr1d-sod-pneg-3b` (a cell at a negative pressure, in block 2) on 3 processes; on 2×2 the probe library of `MPIThreadFunnel`, which the tests build on Linux only (elsewhere this part is left out); `mr2d-plate-4b` on 2 processes with blocks 1-3 on process 0 and block 4 on process 1 by the map, every macro step printed (label `needs-mpi`) | fields, wall output and residual history identical byte by byte, the 1×1 multirate runs with cells on more than one level; `mr2d-plate-2br` against `mr2d-plate-2b`, block 2 read with $i$ and $j$ reversed: the nodes identical and every cell value of $\rho$, $u$, $v$, $w$, $p$ within $10^{-10}$ of its scale, single rate and multirate; `mr2d-box-4bx` on 4 processes conserves mass and energy within $10^{-12}$; `mr1d-acu-3b` with `max-level` 10 has $L = 3$ in every macro step on 1 and 3 processes, and `mr1d-deep12-3b` $L = 10$ and $H = 512\,\Delta t_{min}$; the blast on 3 processes has at least one attempt with a violation; `mr1d-blast-wide-3b` identical on 1×1, 3×1 and 3×2 with at least one attempt rejected for a state that is not physical; `mr1d-blast-hot-3b` stops below $2^{-10}$ of the first attempt after 11 attempts with the lines of the face of process 1, every process stopping with the others (exit code 1, one line of the stop per process, no process ended by a signal); `mr1d-blast-wide-3b` with `max-attempts` 3 stops at the third rejection with the same lines of the limit, of the causes (1 / 0 / 2) and of the last rejection on 1×1, 3×1 and 3×2, one line of the stop per process; `mr1d-sod-pneg-3b` stops at the start of the first macro step with the line of the cell (block 2, cell 5) from process 1 and one line of the stop per process; with the probe library, no MPI call from a thread other than the main one, and the largest MPI tag at most 32767 (the least `MPI_TAG_UB` of the MPI standard) on every rank; with the uneven map, named in the log, on 1 and 2 threads, the exchange wait of the `Timing` block above zero and at least the exchange time per level (mean over the processes, summed over the levels, above zero at some level) over `Solver`, the exchange start at least zero and, with the exchange wait and the collective wait, at most 100 %, the compute spread over ranks reported, and in the end-of-run times per level over the processes minimum ≤ mean ≤ maximum with the values of process 0 between minimum and maximum; on one process none of these lines |
| `mr-fixedstep` | `mr1d-acu`, `mr2d-box`, `fixed-macro-step` 1 and 2 against the adaptive macro step with `macro-step-target` equal to the same $H$, one macro step; the depth of the map, 20 macro steps: `mr1d-acu` (local steps over 3.97 octaves) with `max-level` 10 and with `max-level` left out against `max-level` 3, `mr1d-deep12` with `max-level` left out against 10, `mr1d-simple-64` (uniform mesh) with `max-level` 10 against 0; `mr1d-acu` with `levels = fixed` and `fixed-macro-step` 1, `max-level` 10 against 3; `mr1d-acu` with `macro-step-target` $7\,\Delta t_{first}$ (within the step of the slowest cells over $\sigma$), `max-level` 10 against 4 with the same target, and $20\,\Delta t_{first}$, `max-level` 10 against 3 without a target | level map, macro step and cells per level identical; no cell beyond its limit; outputs and residual history identical byte by byte, with $L$ = 3, 3, 10, 0 ($H$ = 4, 4, 512, 0.5 $\Delta t_{min}$), 3, 4 and 3 ($H = 4\,\Delta t_{min}$) in every macro step, and cells at level 0 in every macro step of the depth runs |
| `mr-rpost` | $r_{post}$ and $r_{nb}$ (section [Non-local CFL condition](../theory/time-integration.md#non-local-cfl-condition)) on `mr1d-sod` up to 1 ms with `max-level` 0 and `safety` 1, and with `max-level` 3 and `safety` 8; on `mr1d-blast` up to 0.6 ms with `max-level` 3 and `safety` 2; on `mr1d-acu` with one level and a fixed macro step of twice the step of the first macro step, 3 macro steps | one level: the largest $r_{nb}$ equals the largest $r_{post}$ (relative $10^{-14}$), the same macro steps above 1, at least one; `safety` 8: no macro step above 1, so no run line of $r_{post}$ or $r_{nb}$ (the line of the largest of both in the macro steps shows that they are evaluated) and no line of their counts in the macro steps; no run line of $h/(\Delta t^{\text{cap}}/\sigma)$ in these two runs, and with the fixed macro step that line with its 3 macro steps above the limit; blast: $r_{nb} > 1$ in at least one macro step, and in at least one macro step a cell with $r_{nb} > 1 \ge r_{post}$ |
| `mr-cfl4` | condition (4) of Gnedin et al. (section [Non-local CFL condition](../theory/time-integration.md#non-local-cfl-condition)): `mr1d-acu-3b`, 100 macro steps, with `nonlocal-cfl = count` without and with `interface-margin`, with `redo` and with `redo` and `force-redo`; `mr1d-uni-3b` (a uniform flow on the same mesh), cfl 0.5, with `count` and with `redo`; `mr1d-acu-3b` with `redo` where $\kappa$ is at its floor (van Leer with cfl 0.5, `max-level` 10 and 3, and 0.25; minmod with 0.6; first order with 0.8) and above it (van Leer with 0.6); `mr1d-blast100` (a blast wave whose front crosses about 2.5 cells per macro step) up to 0.4 ms with `count` and the margin, with `redo` (1 and 3 threads), with `force-redo`, and with `redo` and solutions every `sol-dtime` = 4 µs, one run up to 0.4 ms and one up to 10 `sol-dtime` (40 µs); `mr1d-sod` up to 1 ms with `redo`; `mr1d-blast100` up to 0.4 ms with `redo` and `max-level` 10 (the depth of the map from the octaves of the local steps, 3 at the start); `mr1d-acu-3b`, one macro step with each limiter and with first order; `mr1d-blast-wide` (the blast tube with far coarse cells: local steps over about 11 octaves) up to 3 ms with `max-level` and `nonlocal-cfl` left out (10 levels, `redo`), 1 and 3 threads; `mr1d-sod-hot` (a cell at $2\cdot10^4$ K, outside the tables of the gas) with `redo`; `mr1d-blast-hot` (the hot gas at 9990 K, within the tables up to $10^4$ K) up to 0.4 ms with `redo`; `max-attempts` 1, 2, 3 and 4 on `mr1d-blast-wide` (its first macro step has 4 attempts), 11 and 10 on `mr1d-blast-hot`, 3 on `mr1d-acu-3b` with `count`, and 1 on `mr1d-acu-3b` with `force-redo`; `mr1d-sod-pneg` and `mr1d-sod-p0` (a cell at a negative pressure and one at zero pressure, at rest) with `redo`, and `mr1d-sod-nneg` (a cell at a negative density and pressure) and `mr1d-sod-vac` (a cell of the finest band at a density of $10^{-307}$ kg/m³, about 300 K) with `redo` and the Navier–Stokes equations | acu-3b: without the margin at least one attempt with a violation, with it none; `redo` redoes no macro step and its outputs equal those of `count` with the margin byte by byte; `force-redo` redoes all 100 macro steps once and changes no byte; `mr1d-uni-3b`: no violation, a tie with the limit at the rounding not being one; at the floor of $\kappa$ no macro step redone, the $\kappa$ printed is 1.1 (1.2 above it), and the margin raises at least one cell; blast: with the margin alone at least one attempt with a violation; with `redo` at least one macro step redone, the attempts with a violation are the rejected ones (attempts minus macro steps), mass and energy drift $\le 10^{-12}$, 1 and 3 threads identical byte by byte; `force-redo` changes no byte (violations plus forced rejections = attempts minus macro steps); with solutions every $4\,\mu$s (nearly every macro step cut to land on an output instant) at least one landed macro step whose $H$ a redo halved, the 100 numbered solutions at $k\,$`sol-dtime` and the final one at `time-threshold` bit for bit, one residual row per macro step; up to 10 `sol-dtime`, where the last macro step lands on an output instant and on the end of the run, the 10 numbered solutions, the last one written at the last iteration, at $k\,$`sol-dtime` and the final one at `time-threshold` bit for bit; Sod: the same balance, at most 2 attempts per macro step; blast with `max-level` 10: at least one macro step whose map a redo made deeper than that of its first attempt, the same balance, mass and energy drift $\le 10^{-12}$; the $\theta$ of the limiter in the interface margin at the first macro step is that of the limiter of `input.ini`, which the integrator takes once, at its first macro step ($1/2$ minmod, $\sqrt{2}/2$ van Albada, 1 van Leer, MC and superbee, 0 first order); in every run the kernel of the certificate gives the local step of the macro step back at $t_n$ (no mismatch above $10^{-12}$); a run line of the certificate or of the redo only where its count is not zero (none with the margin on acu-3b, none of the redo with `redo` on acu-3b), and no line of the certificate in the macro steps of `redo` on acu-3b (they come only with a violation or a macro step redone); `mr1d-blast-wide`: at least one attempt rejected for a state that is not physical, the attempts with a violation plus those rejected for a state that is not physical equal to the attempts minus the macro steps, mass and energy drift $\le 10^{-12}$, 1 and 3 threads identical byte by byte, and in every macro step whose rejected attempts were all rejected for a state that is not physical one halving of $H$ for each; `mr1d-sod-hot`: exit code 1 at the first stage of the first macro step with the lines of the face that is not physical, no macro step done; `mr1d-blast-hot`: at least one macro step with an attempt rejected for a face that is not physical, then the stop below $2^{-10}$ of the first attempt after 11 attempts, all rejected for a state that is not physical, with the lines of the first face of the last attempt; `max-attempts` 1, 2 and 3 on `mr1d-blast-wide`: exit code 1 in its first macro step, causes that add up to the limit and do not decrease from one limit to the next, 1 / 0 / 2 with 3, the line of the cause of the last rejection, no macro step done and no line of the floor; 4: the outputs and the run line of the redo of the run without the key, byte by byte; `mr1d-blast-hot` with 11: the lines of the attempts and of the stop of the run without the key, character by character; with 10: the stop at the 10th rejection of the second macro step, causes 0 / 0 / 10; `count` with the key: its warning once and the outputs of the run without the key; `force-redo` with 1: the outputs of the run without the key byte by byte (the attempt rejected by `force-redo` is not counted); `mr1d-sod-pneg`, `mr1d-sod-p0`, `mr1d-sod-nneg` and `mr1d-sod-vac`: exit code 1 at the start of the first macro step with the cell and its local step (not a number, infinite, negative, zero), no macro step done, no rejected attempt and no face that is not physical |
| `mr-restart` | runs in one go against runs in two parts with a restart (`newrun = false`, `vtk raw`) in between: `mr1d-sod`, `mr2d-box-4b` and `mr2d-plate-sa`, multirate and the one-level path, and `mr2d-shock`, multirate, the map recomputed at every macro step; `mr1d-sod`, multirate, the map recomputed every 4 macro steps (restart on a multiple of 4, and without a restart, with `nonlocal-cfl = off`, every macro step printed); `mr1d-sod` and `mr2d-box-4b`, multirate and the one-level path, also with `tecplot ascii`, the restart from `field.tec`; output instants $t_k = k\,$`sol-dtime` and `time-threshold`; `mr2d-shock` with three restarts in 400 macro steps; `mr1d-blast` with `max-level` 10, the depth of the map following the octaves of the local steps, 100 = 50 + restart + 50 | OUTPUT trees identical byte by byte (numbered and final solutions, residual history); with `tecplot ascii` the restart reads `field.tec` and no run writes a vtk file; TIME of the numbered solutions $= k\,$`sol-dtime` and of the final one $=$ `time-threshold` bit for bit; one residual row per macro step, the cut ones included; mass drift $\le 10^{-12}$ over the restarts; with `mr1d-blast` at least two depths in the run in one go and cells at level 0 in each of its macro steps; with `recompute-every` 4 the cells per level change, and only at the macro steps $4m + 1$ |
| `mr-daru` | the viscous shock tube of Daru and Tenaud, Re = 200, lower half of the tube (`mr2d-daru-coarse`, 50 × 25 cells; `mr2d-daru`, 100 × 50, with `./run.sh medium`; `mr2d-daru-full`, 500 × 250, with `./run.sh full`): single rate at $dt_0/4$ (reference) and $dt_0$, multirate with $H = 4\,dt_0$ and $H = dt_0$, ten numbered solutions | every cell within its limit; with $H = dt_0$ the $L^1$ density error and the probe errors $\le 1.05\times$ the single rate's at every instant; the largest density jump on the wall and symmetry rows within one local cell of the reference and no more local extrema than the single rate; mass and energy drift $\le 10^{-12}$ |
| `mr-nozzle` | the start-up of a quasi-1D convergent-divergent nozzle (`mr1d-nozzle-su`, 371 cells): inlet total pressure from $10^5$ to $2.6\cdot10^5$ Pa in 8 ms, outlet at $10^5$ Pa; single rate at $dt_0/4$ (reference), $dt_0$ and $dt_0/2$, multirate with $H = dt_0$ and with the production macro step $H$ and $H/2$ ($\sigma = 2$ and 4), $dt_0 \le 0.45\,\Delta t_{\min}(0)$ | D1: ten solutions with $\rho, p > 0$ at the same times, every cell of the fixed-step runs within its limit; in the throat band: D2 the error of $H = dt_0$ at most 1.05 times the single rate's wherever it has levels, D3 the production runs within the single-rate error plus the step differences and with an observed order $\ge 0.8$ where the single rates converge; D4 the moving shock (jump of the reference $\ge 1.5$) of $H/2$ within $s + \Delta_i + 1$ faces of the reference ($s$ the shift of the single rates, $\Delta_i$ that between $H$ and $H/2$) |
| `mr-restart-hybrid` | `mr2d-box-4b` and `mr2d-plate-4b` on 2 processes × 2 threads, 100 macro steps against 50 + restart + 50 (label `needs-mpi`) | OUTPUT trees identical byte by byte |

No manufactured solution in space and time exercises the multirate integrator.  Its order in time with levels is
measured on a fixed mesh against a single-rate reference with $\Delta t_{\min}/64$ on the same mesh (`mr-order`, `mr-deep`,
`mr-bc`; $\Delta t_{\min}/16$ in the short forms), so these tests measure the error of the time integration and of the
coupling between levels, not that of the
spatial discretisation.  The order tests run on a fixed mesh with a fixed map, where the second order is proven
([Recursive Heun scheme](../theory/time-integration.md#recursive-heun-scheme)); `mr-moving-levels` measures the order with faces between
levels that travel with a wave, where the error is bounded only for linear two-cell fluxes (see [Order with moving faces between levels](#order-with-moving-faces-between-levels)).  The spatial discretisation is the one of the single-rate solver: with one level the integrator is
the single-rate step bit for bit (`mr-onelevel`), whose order is verified by the
[manufactured solution](2D-mms.md) and the [isentropic vortex](2D-isentropic-vortex.md).  Against an exact solution the
multirate error is compared with the single-rate error on one mesh, without an order, in the Riemann problems of
`mr-sod` and `mr-shock`; `mr-moving-levels` measures its order against the exact solution of a simple wave, with the mesh and the
macro step refined together.  The axisymmetric wedge (faces of type 200) runs with one level (`mr2d-wedge-visc`,
`mr-onelevel`) and with levels (`mr2d-wedge-band`, faces between levels on cells with records 200: conservation of mass
and energy with a positive control in `mr-conservation`, halo and ghost records in `mr-halo`).

The full forms decide on the last pair of macro steps, the closest to the asymptotic range
(Appendix A, Remark 2: next to a face between levels the local error is of first order); the short forms
decide on the pair $H_0$, $H_0/2$ with the same thresholds. On the eight cases of `mr-order` the orders of that pair in
the full form are 1.93 to 2.15 ($L^1$) and 1.91 to 2.19 (interface); the smallest margin is on the entropy wave ($L^1$
1.93 against 1.9): a change of the spatial scheme or of the case that lowers the ratio of its errors with $H_0$ and
$H_0/2$ by about 2 % fails the short form, and the full form then decides on the last pair.

## Order in time with levels

The full form of `mr-order` runs every case with the levels fixed from the initial state and the macro steps $H_0$,
$H_0/2$, $H_0/4$ and $H_0/8$ up to $T = 64 H_0$, against a single-rate reference with $\Delta t_{\min}/64$ on the same
mesh. The figures are written by `test/multirate/mr-order/plot_order.py` from the results of the full form.

<figure>
  --8<-- "docs/vv/images/multirate-order.svg"
  Density error of the eight cases against the macro step $H$, relative to the error with $H_0$: left, the $L^1$ error over
  all cells; right, the $L^\infty$ error over the cells within 4 cells of a face between levels (floor 5 % of the global
  $L^\infty$). Dashed: the slope of second order. The legend gives the number of cells of each case and, in brackets,
  the observed order of the last pair; the test asks $\ge 1.9$ in $L^1$ ($\ge 1.8$ for the turbulent channels) and
  $\ge 1.8$ near the faces between levels. Full form.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-wall-s-levels.svg"
  Level map of `mr2d-wall-s` (48 × 8 cells; $\sigma = 2$, $L = 2$, the run with $H_0$, which both forms make; levels 0, 1
  and 2 in grey, blue and green): bands of cells of 8, 2.8, 1.28, 2.8 and 8 mm along
  the no-slip wall (black line, $y = 0$) take the levels 0, 1, 2, 1 and 0, so four faces between levels cut the wall, on cells
  skewed by 14 to 37 degrees at the wall; from the second row up the two columns next to the band of level 2 take level 2
  as well. A wall face of a finer cell reads the ghost cells of its coarser neighbour along
  the wall, which the stages of the finer level refresh at their own times.
</figure>

## Order with moving faces between levels

`mr-moving-levels` refines a simple wave of the Euler equations jointly in space and time ($\Delta x$ and $H$ halved together) up to
$T = 0.5$, against its exact solution by characteristics. With `max-level` 1, `safety` 2 and $H = 0.565\,dt_0$ the faces
between levels lie where the characteristic speed $u + c$ is halfway between its extremes and travel with the wave, on the
characteristics that start at $x = 1/4$ and $3/4$: the errors of order $H^2$ that a macro step leaves in the cells next to a
face between levels ([Time Integration](../theory/time-integration.md#appendix-a-order-and-stability-of-the-multirate-macro-step),
Remark 2) then stay on the same part of the wave, a case that the analysis covers only for linear two-cell fluxes, with an
error $O(H)$ in the maximum norm and $O(H^2)$ in the volume-weighted $L^1$ norm
([Time Integration](../theory/time-integration.md#corollary-every-sequence-of-maps), Corollary), not for the limited
reconstruction and HLLC of this test. The test decides on the $L^1$
order of the run with the map rebuilt at every macro step; the run with the map of the first macro step kept and the
single rate with the macro step $H$ and with $dt_{\min}$ are there for comparison. The figure is written by
`test/multirate/mr-moving-levels/plot_moving_levels.py` from the results of the test.

<figure>
  --8<-- "docs/vv/images/multirate-moving-levels-order.svg"
  Density error of the simple wave against the number of cells $n$ ($\Delta x$ and $H$ halved together, $T = 0.5$): left,
  the $L^1$ error; right, the $L^\infty$ error, for the multirate with the map rebuilt at every macro step (thick), with the
  map of the first macro step kept, and for the single rate with the macro step $H$ and with $dt_{\min}$ (dashed). Black:
  the slopes of second (dashed) and first order (dotted). The legend gives the observed orders of the last pair, $L^1$ and
  $L^\infty$; the test asks $\ge 1.9$ in $L^1$ for the multirate with the map rebuilt. In $L^\infty$ every run, the single
  rate included, stays near 1.3, with the largest error at the maximum of the wave. The test has one
  form.
</figure>

The global $L^\infty$ error of every run has its largest value at the maximum of the wave ($x \approx 0.15$ at $T$),
away from the faces between levels, so it does not isolate the cells next to the moving faces, where the analysis
([Time Integration](../theory/time-integration.md#appendix-a-order-and-stability-of-the-multirate-macro-step),
Remark 2) expects the error to decay more slowly than second order. The test reports that local error as information,
without deciding on it: on the cells within 4 cells of the two faces at the final time, placed from the exact solution,
the $L^\infty$ order of $\rho$, $u$ and $p$ is 1.80, 1.75 and 1.67 on the three pairs with the map rebuilt, while on
the same cells the single rate with the macro step $H$ gives 2.02 to 2.03 on the last pair, and the run with the map
kept, whose faces do not move, 2.10 to 2.11 next to its faces. The test is one-dimensional, inviscid, periodic, in one block and with two
levels; the moving-face case in several dimensions, with several levels, with boundaries or with the viscous terms is
not measured.

## Viscous shock tube of Daru and Tenaud

`mr-daru` runs the viscous shock tube of Daru and Tenaud (Computers & Fluids 30, 89–113, 2001, as Zhou, Xu and Liu state it;
Re = 200, the lower half of the tube) four times on the same mesh: the single rate with $dt_0/4$ (the reference) and with $dt_0$, and the multirate with the macro steps
$H = 4\,dt_0$ (levels up to 3) and $H = dt_0$ (the core cells at the step of the single rate, the cells next to the wall
smaller), each with ten numbered solutions. The short form (`Multirate.Daru`, 50 × 25 cells) runs in CI; the medium
form (`Multirate.DaruMedium`, 100 × 50 cells, `./run.sh medium`) and the full form (`Multirate.DaruFull`,
500 × 250 cells, `./run.sh full`) are heavy tests. The short and the medium forms are judged against their
own reference only; the full form is the resolution at which the wall density is compared with table 1 of Zhou, Xu
and Liu (Physics of Fluids 30, 016102,
2018; arXiv:1705.09062). The figures are written by `test/multirate/mr-daru/plot_daru.py` from the results of the test;
each names its mesh and form in its legend. `./run.sh full` makes the four runs of the full form on one process, in one
block; the figures and the numbers of the full form below come from these four runs made on 12 MPI processes, with the
mesh cut along $x$ into 12 blocks of 41 or 42 columns.

<figure>
  --8<-- "docs/vv/images/multirate-daru-medium-errors.svg"
  Medium form (100 × 50 cells). Left: the density error $e = \sum V\,|\rho - \rho_\mathrm{ref}| / \sum V\,\rho_\mathrm{ref}$
  of the single rate with $dt_0$ and of the two multirate runs at the ten output times. Right: the error of the
  multirate run with $H = dt_0$ over that of the single rate, which the test asks to stay within 1.05 at every output
  time.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-medium-levels.svg"
  Medium form (100 × 50 cells): the share of the cells on the levels 1 to $L$ at every macro step of the two multirate
  runs, stacked (level 0 holds the rest), and their work: the cell updates over those of the single rate with $dt_0$
  over the same time.
</figure>

On the full form the four runs pass the same checks, the error of the multirate run with $H = dt_0$ within 1.0016 times
that of the single rate at every output time. At the final time their wall densities, as $\rho/(\rho_R/1.2)$ with
$\rho_R$ the density of the low-pressure state (the scaling of the non-dimensional problem of the references), differ
from one another by less than $10^{-5}$ (relative) at the 20 abscissae of table 1 of Zhou, Xu and Liu, points of their
1500 × 750 solution, and agree with the table within 0.79 % on average and 4.2 % at most, at $x = 0.6123$.

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-errors.svg"
  Full form (500 × 250 cells). Left: the density error of the single rate with $dt_0$ and of the two multirate runs at
  the ten output times. Right: the error of the multirate run with $H = dt_0$ over that of the single rate, within the
  limit 1.05 of the test at every output time.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-levels.svg"
  Full form (500 × 250 cells): the share of the cells on the levels 1 to $L$ at every macro step of the two multirate
  runs, stacked (level 0 holds the rest), and their work: 0.27 for $H = 4\,dt_0$ and 1.01 for $H = dt_0$.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-wall.svg"
  Full form (500 × 250 cells), final time: the density along the bottom wall of the four runs, which lie on one another,
  and the 20 points of table 1 of Zhou, Xu and Liu.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-density.svg"
  Full form (500 × 250 cells), final time: the density field of the multirate run with $H = 4\,dt_0$ for
  $0.3 \le x \le 1$ and $y \le 0.3$.
</figure>

## Start-up of a nozzle

`mr-nozzle` starts a quasi-1D convergent-divergent nozzle (`mr1d-nozzle-su`: 371 cells in bands of 10, 4, 1.6 and
0.64 mm up to the throat band, then 1.6, 4 and 10 mm in the divergent; area three times that of the throat at both
ends) by raising the inlet total pressure from $10^5$ to $2.6\cdot10^5$ Pa in 8 ms against an outlet at $10^5$ Pa. The
flow starts and chokes at the throat, and the shock that forms in the divergent moves down it, across the face between
the band of 4 mm and that of 10 mm, a face between levels. Six runs on the same mesh: the single rate with $dt_0/4$
(the reference), $dt_0$ and $dt_0/2$; the multirate with $H = dt_0$, whose throat band takes level 1 once the flow is
fast; the multirate with its production macro step $H = 2^L\,\Delta t_{\min}/\sigma$ and with $H/2$ ($\sigma$ = 2 and
4), which land on the output instants. $dt_0 \le 0.45\,\Delta t_{\min}(0)$, below the smallest cell step of the run
(about $0.5\,\Delta t_{\min}(0)$ when the throat is choked), so that every cell of the fixed-step runs stays within its
limit, which the test checks. The criteria are those of the table above; the throat band is where the cells take the
finest levels and where the error follows the step, and the shock is judged where it is a shock, at the outputs where
the density of the reference jumps by at least 1.5 across it. A run in which the finer levels read the coarser cells at
the start of their step instead of at the stage time fails D2; D3, which asks the production runs only to converge with
the step (order at least 0.8: their maps move with the flow, where the second order is not proven, see
[Order with moving faces between levels](#order-with-moving-faces-between-levels)), does not see it. The conservative
correction at the faces between levels is the business of `mr-conservation` (a run without it passes this check). The figures are written by
`test/multirate/mr-nozzle/plot_nozzle.py` from the results of the test.

<figure>
  --8<-- "docs/vv/images/multirate-nozzle-short-shock.svg"
  371 cells, the form of the test: the position of the moving shock (the face of largest compression in the divergent) at
  the outputs where it is a shock, the density of the reference jumping by at least 1.5 across it (from the seventh
  output on); grey lines: the faces between the bands of 1.6 and 4 mm and of 4 and 10 mm, which are faces between
  levels.
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-nozzle-short-errors.svg"
  371 cells, the form of the test. Left: the density error in the throat band of every run against the reference. Right:
  the error of the multirate run with $H = dt_0$ over that of the single rate at the outputs where it has levels, which
  the test bounds by 1.05.
</figure>
