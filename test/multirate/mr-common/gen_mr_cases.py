#!/usr/bin/env python3
"""Generator of the small cases of the multirate tests (test/multirate/mr-*): 1D, 2D and 3D; Euler, laminar
Navier-Stokes, Spalart-Allmaras and SST; perfect-gas air of test/common/Air (cp = 1004.5 J/kg/K, M = 28.970418 kg/kmol)
and the gases that the cases below name; the gas constant is the solver's (mrcheck.RU, FLINT).

Cases (one directory each, files: INPUT/ic.tec, INPUT/bc.txt, input.ini; phase.txt and thermo.dat are copied from
test/common/Air of the MOSE tree given with --air):
  mr1d-acu   1D periodic single block, banded mesh h = 10 / 4 / 1.6 / 0.64 / 1.6 / 4 / 10 mm (ratio 2.5) that gives the
             levels 0|1|2|3|2|1|0 with sigma = 2, L = 3; right-running acoustic wave eps = 1e-3 on u0 = 50 m/s
  mr1d-ent   same mesh, entropy wave rho = rho0 (1 + 0.1 sin), u = 100 m/s, p = p0
  mr1d-sod   same mesh, closed by slip walls (300), Sod state (1.2 kg/m3, 1e5 Pa | 0.15 kg/m3, 1e4 Pa) with the
             discontinuity at x = 0.25 m
  mr1d-mix   same mesh, periodic, laminar with a Schmidt number of 0.7 and chemistry = frozen: the oxidizer mixture
             of test/1D/diff1 with a slug of its fuel mixture for 0.2 < x < 0.3 m, at 1e5 Pa and 300 K, u = 100 m/s: two
             contact discontinuities that cross faces between levels; four species (CH4, H2O, N2, O2), gas files copied
             from test/1D/diff1/INPUT
  mr1d-blast 1D tube of 1 m, 400 uniform cells of 2.5 mm, closed by slip walls (300); cold gas at rest (1e5 Pa, 100 K)
             and hot gas at rest for x < 0.1 m (5e5 Pa, 6400 K: sound speed 8 times higher): a blast wave into cold gas,
             as the point explosion of Gnedin, Semenov and Kravtsov (2018, section 2), whose front enters cells of
             coarser levels (shock Mach 1.907; behind it (u + a)/a1 = 2.42, from the shock-tube relation)
  mr1d-blast100 the same with the hot gas at 1e7 Pa (100 times the cold one): shock Mach 5.93, behind it
             (u + a)/a1 = 7.59 and T = 778 K; in one macro step the front crosses about 2.5 cells, beyond the face
             neighbours covered by the interface margin (test mr-cfl4: the redo of condition (4))
  mr1d-blast100-3b the same in three blocks cut after cells 133 and 266 (the front crosses the first connection at
             about 0.2 ms; test mr-hybrid: the redo on several processes)
  mr1d-blast-wide the tube of mr1d-blast for its first 200 cells of 2.5 mm, then 8 cells whose width doubles from 5 to
             640 mm and 2 cells of 640 mm (210 cells, x up to 3.055 m), with the states of mr1d-blast: the far cells
             admit a step 2^8 times that of the near cold cells, so the local steps span about 11 octaves and the map
             of the default max-level (10) is 10 levels deep; with it a macro step of the front takes a state that is
             not physical inside an attempt (test mr-cfl4: the attempt is rejected and redone)
  mr1d-blast-wide-3b the same in three blocks cut after cells 100 and 200 (the front starts in the first; test
             mr-hybrid)
  mr1d-blast-hot mr1d-blast with the hot gas at 9990 K (5e5 Pa), within the tables of the gas (1 to 10000 K): at the
             contact the reconstructed state of the face of the last hot cell leaves the tables after the first stage
             (test mr-cfl4: attempts rejected for a face that is not physical)
  mr1d-blast-hot-3b the same in three blocks cut after cells 20 and 200: the contact inside block 2 (test mr-hybrid: the
             stop below 2**-10 of the first attempt on several processes, the first face that is not physical on process 1)
  mr1d-sod-hot mr1d-sod with cell 10 (in the band of 10 mm) at 2e4 K and the pressure of its neighbours, outside the
             tables of the gas (1 to 10000 K): the faces of that cell are not physical at the start of the first macro
             step (test mr-cfl4: the run stops at the first stage, as the single-rate step does)
  mr1d-sod-pneg mr1d-sod with cell 10 at the pressure of its neighbours with the opposite sign (-1e5 Pa): its local step is
             not a number (test mr-cfl4: the run stops at the start of the first macro step)
  mr1d-sod-p0 mr1d-sod with cell 10 at zero pressure, the gas at rest: its local step is infinite (test mr-cfl4)
  mr1d-sod-nneg mr1d-sod with cell 10 at a density of -1e-6 kg/m3 and a pressure of -0.0861 Pa (about 300 K): run with the
             Navier-Stokes equations, its diffusivity mu/rho is negative and so is its local step (test mr-cfl4)
  mr1d-sod-vac mr1d-sod with cell 100 (in the band of 0.64 mm) at a density of 1e-307 kg/m3, a normal number, and a
             pressure of 8.6115e-303 Pa (about 300 K): run with the Navier-Stokes equations, its diffusion rate overflows
             and its local step is zero (test mr-cfl4)
  mr1d-sod-pneg-3b mr1d-sod-pneg in three blocks cut after cells 5 and 20: the cell is cell 5 of block 2 (test mr-hybrid:
             the stop on several processes, the cell on process 1)
  mr1d-deep  1D periodic single block, deep hierarchy L = 8 (9 levels): coarse h0 = 10 mm (93 cells), then bands of 4 cells
             with h_l = h0 / 2.1^l for l = 1..7 on each side and a band of 8 cells with h_8 in the middle; with sigma = 2
             and L = 8 the level rule gives ceil(8 - 1.0704 (8 - l)) = l; acoustic wave eps = 1e-3 on u0 = 50 m/s
  mr2d-free  2D single block periodic in x and y (201 on itself), curvilinear stretched mesh 48 x 48:
             x = Lx s(xi) + A sin(2 pi eta), y = Ly s(eta) + A sin(2 pi xi), s(q) = q + b sin(2 pi q)/(2 pi), b = 0.6
             (cell size ratio 4 -> levels), Lx = Ly = 1 m, A = 0.05 m; uniform flow u = 100, v = 50 m/s, p0, T0 (free stream)
  mr2d-box   2D closed box 0.9984 m x 0.9984 m, banded mesh h = 50 / 20 / 8 / 3.2 / 8 / 20 / 50 mm in x and in y
             (refined cross in the middle, 66 x 66 cells), Gaussian pressure pulse
             (amplitude 1e-2 p0, radius 0.08 m) at (0.35, 0.40) m on rest air
  mr3d-sector  3D sector of an annulus around the x axis, 8 x 4 x 2 cells (x in [0, 0.08] m, r in [0.05, 0.09] m,
               theta in [-2.5, 2.5] deg), rest air, k faces rotational periodic (200), other faces slip walls
Cases of several blocks and boundaries (split_blocks cuts a single block into pieces joined by
connections 101, periodic 201 across the periodic ends, records ordered by block and face as BC_Fluxes requires):
  mr1d-deep1     mr1d-deep with bands of ONE cell at levels 1..7 and two cells at level 8 (109 cells): the halo of a
                 level l then holds cells of levels m <= l - 2
  mr1d-acu-3b    mr1d-acu in three blocks cut after cells 15 (inside the 10 mm band) and 50 (on the edge of the 4 mm
                 and 1.6 mm bands, where equal levels across the connection raise cell 50 from level 1 to 2)
  mr1d-acu-3bi   mr1d-acu in three blocks cut after cells 15 and 140 (both inside bands, same level map as mr1d-acu)
  mr1d-uni-3b    the mesh and the blocks of mr1d-acu-3b with a uniform flow, 1.2 kg/m3, 1e5 Pa, 50 m/s (test mr-cfl4: the local
                 steps stay the same up to the rounding)
  mr1d-deep12    deep hierarchy with L = 12 (13 levels): as mr1d-deep with bands of 4 cells for the levels 1..11 on each
                 side and 8 cells at level 12 (189 cells; finest cell 1.36 um)
  mr1d-deep12-1  the same with bands of ONE cell at levels 1..11 and two cells at level 12 (117 cells)
  mr1d-deep12-3b mr1d-deep12 in three blocks cut after cells 105 (edge of the level 3 and 4 bands: cell 105 is raised to
                 level 4) and 141 (middle of the level 12 band)
  mr2d-box-4b    mr2d-box in 2 x 2 blocks cut after cell 33 in i and in j (vertex of 4 blocks in the refined cross)
  mr2d-box-4bx   mr2d-box in 2 x 2 blocks cut after cell 18 in i and in j: the cut is one cell from the edge between
                 the 8 mm and the 3.2 mm bands (cells 17 | 18), so with the blocks on different MPI processes a face
                 between levels is one cell from a process boundary (test mr-hybrid)
  mr2d-wedge-visc  the field with dilatation and swirl of part 4 of the test MMSResiduals on the 2-D axisymmetric wedge
                 (viscous stress on the wedge faces), 32 x 32 cells, laminar
  mr2d-box-sym   mr2d-box with the pulse at the centre of symmetry of the banded mesh, (0.4992, 0.4992) m
  mr2d-box-sym-4b  mr2d-box-sym in 2 x 2 blocks as mr2d-box-4b
  mr2d-free-4b   mr2d-free in 2 x 2 blocks cut after cell 24 in i and in j (connections 101, periodic 201)
  mr2d-vortex-3b isentropic vortex (radius 0.03 m, peak swirl 0.1 a_inf) at (0.15, 0.15) m on u = 100 m/s, periodic
                 in x and y; x: 30 cells of 5 mm | 120 cells of 1.25 mm | 30 cells of 5 mm (three blocks, the central one
                 finer by 4 in x; the vortex sits on the first connection), y: 60 cells of 5 mm. Finer by 4: with the
                 local step that adds up the directions, a ratio of 2.5 gives a step ratio of 1.84 < 2 between the blocks
                 and almost every cell at level 1; with 4 the ratio is about 2.7
  mr2d-shock     mesh of mr2d-box, planar discontinuity at x = 0.2 m (1.2 kg/m3, 1e5 Pa | 0.15 kg/m3, 1e4 Pa) at rest,
                 slip walls
  mr1d-402       1D duct, x from the inlet: 40 cells of 0.64 mm, 20 of 1.6 mm, 20 of 4 mm, 30 of 10 mm; inlet 402 at
                 i = 1 (T0 = 300 K, p0(t) = 1.2e5 Pa (1 + 0.1 sin(2 pi t / Tp)), Tp = 6.4e-5 s, table p0.dat in the
                 case directory, read from the run directory), outlet 406 at 1e5 Pa; uniform initial state = the steady
                 isentropic state of p0 = 1.2e5 Pa (as test/checks/time-integration/time-control); y faces slip
                 walls, z faces null (0)
  mr1d-410       the same duct with the inlet 410: assigned state of a right-running acoustic wave of relative amplitude
                 0.1 on the steady state, rho (1 + 0.1 s / gamma), u + 0.1 c s / gamma, p (1 + 0.1 s), s = sin(2 pi t / Tp),
                 tabulated every Tp/400 in the Q2D file q2d.dat (Tecplot ASCII read by src/lib/io/IO_BC_Q2D.f90)
  mr1d-401, mr1d-403, mr1d-405, mr1d-408, mr1d-420  the same duct with the inlet of that type on face 1 (data line of
                 docs/user/boundary-conditions.md, one species of mass fraction 1, flow normal to the face): 401 the
                 total conditions T0, p0 = 1.2e5 Pa; 403 T0 and the mass flux of the steady state; 405 a supersonic
                 inflow at Mach 1.5 with the static temperature and pressure of the steady state; 408 its static
                 temperature and velocity; 420 T0, p0 = 1.2e5 Pa, p_sub = 0.5 and p_sup = 0.25 times the outlet pressure
                 (the cell pressure stays above p_sub: the branch of the total conditions) and the mass flux of the steady
                 state; the same outlet 406 (test mr-onelevel: the boundary fluxes of these types at one level)
  mr1d-502       the same duct with an SRM grain in its inflow mode on face 1 (502, mode 2: the record becomes the
                 mass-flux inlet 404 with srm_inflow, IO_BC.f90): mdot_ref the mass flux of the steady state, n = 0.5,
                 p_ref the steady pressure, T_af = T0, one species of mass fraction 1; the same outlet 406 (test
                 mr-onelevel: the boundary flux of the grain at one level)
  mr1d-502c      the duct of mr1d-502 closed at its end by a slip wall (300) and at rest (the density and pressure of
                 the steady state, u = 0), with the burn-rate exponent n = 0 (a plateau propellant: the mass flux of
                 the grain is mdot_ref at every pressure; test mr-conservation: the mass of the domain grows by
                 mdot_ref A t)
  mr1d-nozzle    the mesh and the Sod state of mr1d-sod with the top j face (j = 2) and the top k face (k = 2) inclined:
                 the height and the thickness grow linearly from 10 mm at x = 0 to 15 mm at x = LDOM, the bottom faces
                 stay plane; every face a slip wall (300). A 1-D case with symmetry j and k faces that are not planes
                 normal to y and z: the set-up keeps the y and the z velocity (freeze_v and freeze_w off, with a
                 warning), and the step integrates the residual that the pressure on the inclined faces gives them
                 (test mr-onelevel)
  mr1d-skew      the mesh and the Sod state of mr1d-sod with the nodes of the top j row moved by SKEW_J = 1 mm along x
                 and those of the top k row by SKEW_K = 2 mm: the j and k faces stay planes normal to y and z, the i faces
                 are planes inclined to x (normal along (1, -0.1, -0.2)); every face a slip wall (300). A 1-D case with
                 symmetry j and k faces (freeze_v and freeze_w of the single-rate step): the pressure on the inclined i
                 faces gives the y and the z momentum a residual that the step does not integrate (test mr-onelevel)
  mr2d-plate-rev mr2d-plate (below) with the records of bc.txt in reverse order: in the corner cells of the wall (301) and
                 of the ambient faces (407) two boundary faces add to the same rows of the residual, which BC_Fluxes sums
                 in the order of the faces 1..6 whatever the order of bc.txt (test mr-onelevel)
Cases of the laminar Navier-Stokes equations (input.ini with equations = navier-stokes, INI_VISC):
  mr2d-plate     the laminar flat plate of test/2D/viscous/flat-plate-laminar (64 x 64 cells down to 16 um at the wall,
                 uniform start at Mach 0.2, inlet 404, ambient 407, symmetry 300 and adiabatic wall 301, k faces 200),
                 air of test/common/Air-viscous; files copied as text (plate_case), cfl and vnn 0.5
  mr2d-plate-2b  the plate in two blocks cut in the middle of j (x): connection 101 across the boundary layer
  mr2d-plate-2br mr2d-plate-2b with block 2 turned by 180 degrees about z, its i and j reversed (turn_block_text): the
                 cut is face 4 of both blocks, connected with the orientation -1 0 0 1 (the first tangential direction of
                 the face, i, reversed)
  mr2d-plate-4b  the plate in 2 x 2 blocks (vertex of four blocks in the boundary layer;
                 split_text)
  mr2d-hot       skewed periodic parallelogram 60 x 60, hot square (8 T0) advected by u = 100, v = 50 m/s, viscosity
                 0.3 Pa s: square of level 1 with sharp corners (halo cells reached only through the diagonals), cfl 0.2
  mr2d-hot-4b    mr2d-hot in 2 x 2 blocks cut through the middle of the square
  mr2d-wall-s    a no-slip adiabatic wall (301) on smoothly skewed cells whose levels change along the wall: one block
                 48 x 8, periodic in x (201) with bands of 8 / 2.8 / 1.28 / 2.8 / 8 mm (levels 0 | 1 | 2 | 1 | 0 along
                 the wall, in x, sigma = 2, L = 2), rows of 16 mm, slip wall (300) on top; x' = x + s(y) (1 + 0.5
                 sin(2 pi x/Lx)), s(y) = 0.5 lambda (1 - exp(-y/lambda)), lambda = 16 mm: wall nodes fixed, skew angle
                 at the wall 14 to 37 degrees, smooth grid lines; u = 100 (1 - exp(-y/8 mm)) m/s, T = T0 (1 + 0.25 (1 -
                 cos(2 pi x/Lx))), viscosity HOT_MU (wall_s_block)
  mr2d-plate-mu  the plate with mu = 1e-3 Pa s and Pr = 1 (the gas of the [GPB-Phase1] block of that test): steady
                 state of the explicit local step in a few 1e5 iterations
  mr1d-acu-sa    mr1d-acu with the variable of Spalart-Allmaras (rejection of the RANS runs)
  mr2d-plate-302 the plate with its adiabatic wall (301) made an isothermal wall (302) at PLATE_TW (test mr-onelevel)
  mr2d-wedge-band  the 2-D axisymmetric wedge of 1 degree (k faces of type 200, axis row at 1e-6 of the first cell height
                 on face 3, type 300) with bands of rows in r from the axis, WB_RBANDS: 8 rows each of 0.125, 0.25, 0.5 and
                 1 mm (r up to 15 mm), 24 columns of 1 mm in x; every other face a slip wall (300), so that mass and energy
                 are conserved; rest air of test/common/Air-viscous at T0 and P0 with a solid-body swirl (w = WB_SWIRL r/r_max)
                 and a Gaussian pressure pulse (1e-2 P0, radius 3 mm, isentropic) at (12, 6) mm; laminar, cfl 0.5, vnn 0.3.
                 With safety 2 and max-level 2 the levels are 2 | 2 | 1 | 0 from the axis outward: faces between levels
                 on cells with records 200 (tests mr-conservation and mr-halo)
  mr3d-hot       the 3D twin of mr2d-hot: 12 x 12 x 12 cells of 10 mm on an oblique affine mesh, periodic on the six
                 faces, a hot cube (8 T0) advected by u = 100, v = 50, w = 30 m/s, viscosity 0.3 Pa s: a cube of level 1
                 with sharp edges and vertices, cfl 0.2 (hot3d_case)
  mr3d-hot-8b    mr3d-hot in 2 x 2 x 2 blocks cut through the middle of the hot cube: a vertex of eight blocks inside the
                 cells of level 1, connections 101 and periodic 201 (test mr-hybrid)
Cases of the Spalart-Allmaras and SST models (turb_plate_case, chan_case):
  mr2d-plate-sa, mr2d-plate-sst  the turbulent flat plate of mr-common/plate on its coarse mesh
                 (68 x 48 cells, 4 um at the wall), uniform start at Mach 0.2, with Spalart-Allmaras or SST; cfl 0.5,
                 vnn 0.4
  mr2d-plate-sa-4b, mr2d-plate-sst-4b  mr2d-plate-sa and mr2d-plate-sst in 2 x 2 blocks
  mr2d-chan-sa, mr2d-chan-sst  the wall channel of the test TimeIntegration/TurbulenceSources (16 wall rows from 4 um
                 growing by 1.25, symmetry plane at 0.553 mm), periodic along the flow, with the profiles of that test
                 (Spalart-Allmaras at twice the near-wall equilibrium; SST with omega-wall-bc = asymptotic); cfl 0.8,
                 vnn 0.4
  mr2d-chanrough-sa  mr2d-chan-sa with a rough wall, sand-grain height 10 um (SA-rough)
  mr2d-chanpos-sa, mr2d-chanpos-sst  the positivity cases of that test (18 rows from 2 um growing by 1.25, then 16
                 rows of 100 um): near-wall turbulence far above its equilibrium (nu~ times 100; omega times 10 with
                 omega-wall-bc = practical), where the destruction acts faster than the transport
Viscous shock tube of Daru and Tenaud (Navier-Stokes, Re = 200; daru_case):
  mr2d-daru      the lower half [0; 1] x [0; 0.5] m of the tube, 100 x 50 cells: x uniform, y stretched towards the wall
                 y = 0 (first cell a quarter of dx), one cell of 10 mm in z; walls 301 (adiabatic, no slip) at x = 0,
                 x = 1 and y = 0, symmetry 300 at y = 0.5 and on the k faces; air of test/common/Air with a constant
                 viscosity mu = rho_ref a_ref L / Re (rho_ref = rho_R / 1.2, a_ref = a(300 K), L = 1 m) and prandtl
                 0.73; T = 300 K, 1e7 Pa for x < 0.5 and 1e5 Pa beyond, at rest; time-threshold L / a_ref; cfl 0.8,
                 vnn 0.3
  mr2d-daru-full the same tube on 500 x 250 cells, the resolution of the reference solutions
  mr2d-daru-coarse the same tube on 50 x 25 cells
Start-up of a nozzle (Euler; nozzle_case):
  mr1d-nozzle-su a quasi-1D convergent-divergent nozzle, x from 0 to 0.992 m: bands of 10 / 4 / 1.6 / 0.64 mm up to the
                 throat band (0.444-0.508 m, throat at 0.476 m), then 1.6 / 4 / 10 mm (0.508-0.700, 0.700-0.912 and
                 0.912-0.992 m: 371 cells); area A(x) / A_t = 1 + a (x - 0.476)^2 with A / A_t = 3 at both ends (two
                 coefficients), made by the top j and k faces (height and thickness 10 mm sqrt(A / A_t)) as in
                 mr1d-nozzle, every j and k face a slip wall (300); inlet 402 at i = 1 (T0 = 300 K, p0(t) from P0 to
                 2.6 P0 linearly in 8 ms, then constant: table p0.dat in the case directory), outlet 406 at P0; air at
                 rest at P0 and T0. The flow starts, chokes at the throat, and the shock that forms in the divergent
                 moves down it across the band of 4 mm into that of 10 mm (time-threshold 10 ms)
Simple wave with moving level boundaries (Euler; simple_case; test mr-moving-levels):
  mr1d-simple-<n> the 1D periodic simple wave on [0; 1] in n uniform cells: c = 1 + 0.05 phi, phi = (1 - cos 2 pi x) / 2,
                 u = 2 (c - 1) / (gamma - 1), rho = c^(2 / (gamma - 1)), p = rho c^2 / gamma, on the non-dimensional gas of
                 test/common/nondim-Air (gamma = 1.3999968 with the gas constant of the solver); exact cell averages of the
                 conserved variables; Heun, HLLC, MUSCL with van Leer, cfl 0.8, time-threshold 0.5 (the wave breaks at
                 t = 1.06)
cfl: in the cases built by make() 0.8 in 1D, 0.2 in mr2d-box, 0.1 in mr2d-free and 0.4 in mr3d-sector, 0.2 in the other 2D cases
(box, vortex, shock, hot, wall-s) and 0.1 in mr2d-free-4b; mr3d-hot has 0.2, and the cases with a function of their
own (plates, channels, wedges, Daru and Tenaud, nozzle, simple wave) set theirs.
The y and z faces of the 1D cases and the z faces of the 2D case are slip walls (300). One cell in z: 10 mm thick in the
cartesian cases (1D and 2D; the 1D cases are 10 mm in y too), 1 mm in mr2d-free and mr2d-free-4b (see DZ below). Usage: gen_mr_cases.py --air DIR_COMMON_AIR --out DIR [case ...]   (refuses existing directories)
"""
import sys
sys.dont_write_bytecode = True
import argparse
import math
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mrcheck import RU     # noqa: E402  universal gas constant of the solver (FLINT), one place for all the tests
BANDS_1D = [(0.0, 0.3, 0.01), (0.3, 0.38, 0.004), (0.38, 0.444, 0.0016), (0.444, 0.508, 0.00064),
            (0.508, 0.572, 0.0016), (0.572, 0.652, 0.004), (0.652, 0.992, 0.01)]
LDOM = 0.992
P0, T0 = 1.0e5, 300.0
# Thickness of the slabs (one cell between two slip walls): the wall mode across a slab, dw/dt = -2 a w / dz, is
# unstable with Heun for a step above dz/a. The cartesian slabs (1D: y and z; 2D: z) are
# 10 mm thick, so every step of the multirate tests stays below dz/a (about 29 us; the largest macro step is about 7 us);
# mr2d-free keeps 1 mm, where the steps at cfl 0.1 stay below dz/a (2.3 us against 2.9 us).
DZ = 0.01
DZ_FREE = 0.001
# mr1d-skew: move along x of the nodes of the top j row and of the top k row of the slab [m]
SKEW_J, SKEW_K = 0.001, 0.002


def gas(air):
    mw = float(open(os.path.join(air, 'phase.txt')).read().split()[-1])
    R = RU / mw
    th = np.loadtxt(os.path.join(air, 'thermo.dat'), skiprows=4)
    cp = float(th[0, 1])
    assert np.all(th[:, 1] == cp), 'cp not constant in the table'
    return R, cp / (cp - R)


def banded(bands):
    xs = [bands[0][0]]
    for a, b, h in bands:
        n = int(round((b - a) / h))
        assert abs(n * h - (b - a)) < 1e-12, (a, b, h)
        for k in range(1, n + 1):
            xs.append(a + k * h if k < n else b)
    return np.array(xs)


def write_tec(path, blocks, names):
    """blocks: list of dict(X, Y, Z nodal arrays (K, J, I); cell vars (K-1, J-1, I-1))"""
    with open(path, 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"  ' + ' '.join('"%s"' % n for n in names) + '\n')
        for b, B in enumerate(blocks, 1):
            K, J, I = B['X'].shape
            f.write(' ZONE  T = B%d-IG, I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED),'
                    ' SOLUTIONTIME=0.000000000000000E+000\n' % (b, I, J, K, 3 + len(names)))
            for a in (B['X'], B['Y'], B['Z']):
                f.write(''.join(' %.16E\n' % v for v in a.ravel()))
            for n in names:
                f.write(''.join(' %.16E\n' % v for v in B[n].ravel()))


def rec(b, i, j, k, f, t, prop=None):
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, f, t)
    if prop is not None:
        s += '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % prop
    return s


def bc_single_block(ni, nj, xperiodic):
    """one block (ni, nj, 1): faces 1-2 periodic on itself (201) or slip walls; faces 3-6 slip walls; order by face"""
    out = []
    for j in range(1, nj + 1):
        out.append(rec(1, 1, j, 1, 1, 201, (1, ni, j, 1, 2, 1, 0, 0, 1)) if xperiodic else rec(1, 1, j, 1, 1, 300))
    for j in range(1, nj + 1):
        out.append(rec(1, ni, j, 1, 2, 201, (1, 1, j, 1, 1, 1, 0, 0, 1)) if xperiodic else rec(1, ni, j, 1, 2, 300))
    for f in (3, 4):
        for i in range(1, ni + 1):
            out.append(rec(1, i, 1 if f == 3 else nj, 1, f, 300))
    for f in (5, 6):
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                out.append(rec(1, i, j, 1, f, 300))
    return ''.join(out)


def reverse_records(text):
    """the records of a bc.txt in reverse order: a record is a line of six integers and the property lines after it"""
    recs = []
    for ln in text.splitlines(True):
        t = ln.split()
        if len(t) == 6 and all(x.lstrip('-').isdigit() for x in t):
            recs.append([ln])
        else:
            recs[-1].append(ln)
    return ''.join(''.join(r) for r in reversed(recs))


def slab_1d(xn):
    I = xn.size
    X = np.empty((2, 2, I)); Y = np.empty((2, 2, I)); Z = np.empty((2, 2, I))
    for k in range(2):
        for j in range(2):
            X[k, j] = xn; Y[k, j] = (0.0, 0.01)[j]; Z[k, j] = (-0.5 * DZ, 0.5 * DZ)[k]
    return X, Y, Z


def split_blocks(B, cuts_i, cuts_j, per_i, per_j):
    """Pieces of the single block B (nodal X, Y, Z of shape (2, J, I); cell variables (1, J-1, I-1)) cut after the cells
    cuts_i along i and cuts_j along j, numbered with i fastest. Faces between pieces: connection 101; faces at the ends:
    periodic 201 to the piece at the other end (per_i, per_j) or slip wall 300; k faces: slip walls. Connection data
    'bs is js ks fs d11 d12 d21 d22' with the identity orientation (1 0 0 1). Records ordered by block, then face 1..6
    (BC_Fluxes takes the records of face f of block b as one contiguous range)."""
    _, J, I = B['X'].shape
    ei = [0] + list(cuts_i) + [I - 1]
    ej = [0] + list(cuts_j) + [J - 1]
    nbi, nbj = len(ei) - 1, len(ej) - 1

    def bid(bi, bj):
        return bj * nbi + bi + 1

    def ni_(bi):
        return ei[bi + 1] - ei[bi]

    def nj_(bj):
        return ej[bj + 1] - ej[bj]

    blocks, out = [], []
    for bj in range(nbj):
        for bi in range(nbi):
            P = {c: B[c][:, ej[bj]:ej[bj + 1] + 1, ei[bi]:ei[bi + 1] + 1] for c in ('X', 'Y', 'Z')}
            for v in ('rho1', 'u', 'v', 'w', 'p'):
                P[v] = B[v][:, ej[bj]:ej[bj + 1], ei[bi]:ei[bi + 1]]
            blocks.append(P)
    for bj in range(nbj):
        for bi in range(nbi):
            b, ni, nj = bid(bi, bj), ni_(bi), nj_(bj)
            for j in range(1, nj + 1):                                              # face 1, i = 1
                if bi > 0:
                    out.append(rec(b, 1, j, 1, 1, 101, (bid(bi - 1, bj), ni_(bi - 1), j, 1, 2, 1, 0, 0, 1)))
                elif per_i:
                    out.append(rec(b, 1, j, 1, 1, 201, (bid(nbi - 1, bj), ni_(nbi - 1), j, 1, 2, 1, 0, 0, 1)))
                else:
                    out.append(rec(b, 1, j, 1, 1, 300))
            for j in range(1, nj + 1):                                              # face 2, i = ni
                if bi < nbi - 1:
                    out.append(rec(b, ni, j, 1, 2, 101, (bid(bi + 1, bj), 1, j, 1, 1, 1, 0, 0, 1)))
                elif per_i:
                    out.append(rec(b, ni, j, 1, 2, 201, (bid(0, bj), 1, j, 1, 1, 1, 0, 0, 1)))
                else:
                    out.append(rec(b, ni, j, 1, 2, 300))
            for i in range(1, ni + 1):                                              # face 3, j = 1
                if bj > 0:
                    out.append(rec(b, i, 1, 1, 3, 101, (bid(bi, bj - 1), i, nj_(bj - 1), 1, 4, 1, 0, 0, 1)))
                elif per_j:
                    out.append(rec(b, i, 1, 1, 3, 201, (bid(bi, nbj - 1), i, nj_(nbj - 1), 1, 4, 1, 0, 0, 1)))
                else:
                    out.append(rec(b, i, 1, 1, 3, 300))
            for i in range(1, ni + 1):                                              # face 4, j = nj
                if bj < nbj - 1:
                    out.append(rec(b, i, nj, 1, 4, 101, (bid(bi, bj + 1), i, 1, 1, 3, 1, 0, 0, 1)))
                elif per_j:
                    out.append(rec(b, i, nj, 1, 4, 201, (bid(bi, 0), i, 1, 1, 3, 1, 0, 0, 1)))
                else:
                    out.append(rec(b, i, nj, 1, 4, 300))
            for f in (5, 6):                                                        # k faces
                for j in range(1, nj + 1):
                    for i in range(1, ni + 1):
                        out.append(rec(b, i, j, 1, f, 300))
    return blocks, ''.join(out)


def acu_block(gam, rho0, c0, eps=1.0e-3):
    """Single block of mr1d-acu: banded 1D mesh, right-running acoustic wave eps = 1e-3 on u0 = 50 m/s."""
    xn = banded(BANDS_1D)
    X, Y, Z = slab_1d(xn)
    xc = 0.5 * (xn[1:] + xn[:-1])
    u0 = 50.0
    s = np.sin(2 * np.pi * xc / LDOM)
    p = P0 * (1 + eps * s); rho = rho0 * (1 + eps * s / gam); u = u0 + eps * c0 / gam * s
    sh = (1, 1, xc.size)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
            'w': np.zeros(sh), 'p': p.reshape(sh)}


def deep_block(nband, ndeep, gam, rho0, c0, L=8):
    """Single block of the mr1d-deep family: 93 cells of h0 = 10 mm, then bands of nband cells with h_l = h0/2.1**l for
    l = 1..L-1 on each side of a band of ndeep cells with h_L; acoustic wave eps = 1e-3 on u0 = 50 m/s, periodic. With
    safety 2 and max-level L the level rule gives ceil(l - 0.0704 (L - l)) = l on every band (L <= 14)."""
    h0, r = 0.01, 2.1
    widths = [h0] * 93
    for l in range(1, L):
        widths += [h0 / r ** l] * nband
    widths += [h0 / r ** L] * ndeep
    for l in range(L - 1, 0, -1):
        widths += [h0 / r ** l] * nband
    xn = np.concatenate([[0.0], np.cumsum(widths)])
    ldom = float(xn[-1])
    X, Y, Z = slab_1d(xn)
    xc = 0.5 * (xn[1:] + xn[:-1])
    eps, u0 = 1.0e-3, 50.0
    s = np.sin(2 * np.pi * xc / ldom)
    p = P0 * (1 + eps * s); rho = rho0 * (1 + eps * s / gam); u = u0 + eps * c0 / gam * s
    sh = (1, 1, xc.size)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
            'w': np.zeros(sh), 'p': p.reshape(sh)}


BOX_BANDS = [(0.0, 0.3, 0.05), (0.3, 0.4, 0.02), (0.4, 0.448, 0.008), (0.448, 0.5504, 0.0032),
             (0.5504, 0.5984, 0.008), (0.5984, 0.6984, 0.02), (0.6984, 0.9984, 0.05)]


def box_mesh():
    xn = banded(BOX_BANDS); yn = xn.copy()
    I, J = xn.size, yn.size
    X = np.empty((2, J, I)); Y = np.empty((2, J, I)); Z = np.empty((2, J, I))
    for k in range(2):
        X[k] = xn[None, :]; Y[k] = yn[:, None]; Z[k] = (-0.5 * DZ, 0.5 * DZ)[k]
    xc = 0.5 * (xn[1:] + xn[:-1]); yc = 0.5 * (yn[1:] + yn[:-1])
    XC, YC = np.meshgrid(xc, yc)
    return X, Y, Z, XC, YC


def box_block(xp, yp, gam, rho0):
    """Single block of mr2d-box: banded 2D mesh, Gaussian pressure pulse (1e-2 p0, radius 0.08 m) at (xp, yp) on rest air."""
    X, Y, Z, XC, YC = box_mesh()
    r2 = ((XC - xp) ** 2 + (YC - yp) ** 2) / 0.08 ** 2
    p = P0 * (1 + 1e-2 * np.exp(-r2)); rho = rho0 * (p / P0) ** (1 / gam)
    sh = (1,) + XC.shape
    z = np.zeros(sh)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': z, 'v': z, 'w': z, 'p': p.reshape(sh)}


def free_block(rho0):
    """Single block of mr2d-free: curvilinear stretched and skewed 48 x 48 mesh, uniform flow u = 100, v = 50 m/s."""
    n, b, A = 48, 0.6, 0.05
    q = np.linspace(0.0, 1.0, n + 1)
    sq = q + b * np.sin(2 * np.pi * q) / (2 * np.pi)
    XI, ET = np.meshgrid(q, q)              # (J, I)
    SX, SY = np.meshgrid(sq, sq)
    Xn = SX + A * np.sin(2 * np.pi * ET)
    Yn = SY + A * np.sin(2 * np.pi * XI)
    X = np.empty((2, n + 1, n + 1)); Y = np.empty((2, n + 1, n + 1)); Z = np.empty((2, n + 1, n + 1))
    for k in range(2):
        X[k] = Xn; Y[k] = Yn; Z[k] = (-0.5 * DZ_FREE, 0.5 * DZ_FREE)[k]
    sh = (1, n, n)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.full(sh, rho0), 'u': np.full(sh, 100.0), 'v': np.full(sh, 50.0),
            'w': np.zeros(sh), 'p': np.full(sh, P0)}


def vortex_block(R, gam, rho0):
    """Single block of mr2d-vortex-3b before the cuts (see the module docstring). Isentropic vortex: swirl
    eps a_inf r exp((1 - r^2)/2) with r = distance/Rv, T/T_inf = 1 - (gamma - 1)/2 eps^2 exp(1 - r^2), rho and p
    isentropic (an exact steady solution of the Euler equations in the frame of the vortex)."""
    xn = np.concatenate([[0.0], np.cumsum([0.005] * 30 + [0.00125] * 120 + [0.005] * 30)])
    yn = np.concatenate([[0.0], np.cumsum([0.005] * 60)])
    I, J = xn.size, yn.size
    X = np.empty((2, J, I)); Y = np.empty((2, J, I)); Z = np.empty((2, J, I))
    for k in range(2):
        X[k] = xn[None, :]; Y[k] = yn[:, None]; Z[k] = (-0.5 * DZ, 0.5 * DZ)[k]
    xc = 0.5 * (xn[1:] + xn[:-1]); yc = 0.5 * (yn[1:] + yn[:-1])
    XC, YC = np.meshgrid(xc, yc)
    rv, eps = 0.03, 0.1
    a_inf = math.sqrt(gam * R * T0)
    dx, dy = (XC - 0.15) / rv, (YC - 0.15) / rv
    r2 = dx ** 2 + dy ** 2
    f = np.exp(0.5 * (1.0 - r2))
    tr = 1.0 - 0.5 * (gam - 1.0) * eps ** 2 * np.exp(1.0 - r2)
    sh = (1, J - 1, I - 1)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': (rho0 * tr ** (1.0 / (gam - 1.0))).reshape(sh),
            'u': (100.0 - eps * a_inf * dy * f).reshape(sh), 'v': (eps * a_inf * dx * f).reshape(sh),
            'w': np.zeros(sh), 'p': (P0 * tr ** (gam / (gam - 1.0))).reshape(sh)}


HOT_L, HOT_N, HOT_SKEW, HOT_TR, HOT_MU = 0.3, 60, 0.5, 8.0, 0.3


def hot_block(R, gam):
    """Single block of mr2d-hot (before the cuts of mr2d-hot-4b): 2D periodic parallelogram, 60 x 60 cells of 5 mm,
    x = xi + 0.5 eta, y = eta (skewed cells, |cos| of the angle between the edges 0.45: the face metric mixes the
    directions, so a diffusive face reads its tangential neighbours for every variable). State: uniform flow u = 100,
    v = 50 m/s at p0, temperature T0 (1 + 7 f) with f = (Sx^-8 + Sy^-8)^(-1/8), a smooth minimum of the plateaus
    Sx(xi), Sy(eta) (1 for |s - 1/2| <= 0.1, cosine to 0 within 0.05 more, s = coordinate / 0.3 m): a hot square,
    advected, where the speed of sound doubles; the level map (sigma = 2, L = 1) is a square of level 1 with sharp
    corners, where a diffusive face of the corner cell reads a coarser cell that only the diagonal reaches (the cells
    a halo along the axes misses). Viscosity HOT_MU = 0.3 Pa s (transport table of the case,
    Pr 0.72; thermal diffusivity up to 4 m2/s in the square, diffusion number about 0.16 at the step of the hot cells):
    the heat flux across the edges of the square is of the order of the convective enthalpy flux, and an error of the
    diffusive coupling at the corners shows in the cells next to the faces between levels."""
    L, n = HOT_L, HOT_N
    q = np.linspace(0.0, L, n + 1)
    XI, ET = np.meshgrid(q, q)              # (J, I)
    X = np.empty((2, n + 1, n + 1)); Y = np.empty((2, n + 1, n + 1)); Z = np.empty((2, n + 1, n + 1))
    for k in range(2):
        X[k] = XI + HOT_SKEW * ET; Y[k] = ET; Z[k] = (-0.5 * DZ, 0.5 * DZ)[k]
    qc = 0.5 * (q[1:] + q[:-1])
    XC, EC = np.meshgrid(qc, qc)

    def plateau(s):
        s = np.abs(s)
        return np.where(s <= 0.1, 1.0, np.where(s >= 0.15, 0.0, 0.5 * (1.0 + np.cos(np.pi * (s - 0.1) / 0.05))))
    sx, sy = plateau(XC / L - 0.5), plateau(EC / L - 0.5)
    inside = (sx > 0.0) & (sy > 0.0)
    f = np.zeros(sx.shape)
    f[inside] = (sx[inside] ** -8 + sy[inside] ** -8) ** (-1.0 / 8.0)
    T = T0 * (1.0 + (HOT_TR - 1.0) * f)
    sh = (1, n, n)
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': (P0 / (R * T)).reshape(sh), 'u': np.full(sh, 100.0),
            'v': np.full(sh, 50.0), 'w': np.zeros(sh), 'p': np.full(sh, P0)}


WALLS_BANDS = [(8.0e-3, 8), (2.8e-3, 8), (1.28e-3, 16), (2.8e-3, 8), (8.0e-3, 8)]
WALLS_NY, WALLS_HY, WALLS_LAM, WALLS_TW, WALLS_U0, WALLS_DELTA = 8, 16.0e-3, 16.0e-3, 0.5, 100.0, 8.0e-3


def wall_s_block(R, gam):
    """Single block of mr2d-wall-s: a no-slip adiabatic wall on smoothly skewed cells whose levels change along the wall.
    x periodic, bands of WALLS_BANDS (cell widths 8, 2.8, 1.28, 2.8, 8 mm, Lx = 193.28 mm), WALLS_NY rows of
    WALLS_HY in y (Ly = 128 mm), the wall at y = 0 (face 3). With sigma = 2 and L = 2 the step limit of the cells
    (convective and viscous) puts the bands on the levels 0 | 1 | 2 | 1 | 0 along x, so four faces
    between levels cut the wall: a wall face of a finer cell reads, for its tangential derivatives, the ghost cells of its
    neighbour along the wall on a coarser level, which the stages of the finer level must refresh at their own times.
    Smooth skew (a shear with a wave): x' = x + s(y) (1 + 0.5 sin(2 pi
    x/Lx)), s(y) = WALLS_TW WALLS_LAM (1 - exp(-y/WALLS_LAM)): the wall nodes stay in place, the angle of the cells at the
    wall goes from 14 to 37 degrees along x and decreases away from it, the grid lines have no corner, the widths stay
    within 13 % of those of the bands, and the two periodic faces move by the same s(y). State: p0, u = WALLS_U0
    (1 - exp(-y/WALLS_DELTA)) (zero at the wall), v = 0, T = T0 (1 + 0.25 (1 - cos(2 pi x/Lx))); viscosity HOT_MU, as
    mr2d-hot, so that the viscous coupling across the faces between levels is of the order of the convective one."""
    widths = []
    for h, n in WALLS_BANDS:
        widths += [h] * n
    xn = np.concatenate([[0.0], np.cumsum(widths)])
    lx = xn[-1]
    yn = WALLS_HY * np.arange(WALLS_NY + 1)
    XN, YN = np.meshgrid(xn, yn)            # (J, I)
    s = WALLS_TW * WALLS_LAM * (1.0 - np.exp(-YN / WALLS_LAM))
    XS = XN + s * (1.0 + 0.5 * np.sin(2.0 * np.pi * XN / lx))
    X = np.empty((2,) + XN.shape); Y = np.empty_like(X); Z = np.empty_like(X)
    for k in range(2):
        X[k] = XS; Y[k] = YN; Z[k] = (-0.5 * DZ, 0.5 * DZ)[k]
    xc = 0.25 * (XS[1:, 1:] + XS[1:, :-1] + XS[:-1, 1:] + XS[:-1, :-1])
    yc = 0.25 * (YN[1:, 1:] + YN[1:, :-1] + YN[:-1, 1:] + YN[:-1, :-1])
    T = T0 * (1.0 + 0.25 * (1.0 - np.cos(2.0 * np.pi * xc / lx)))
    sh = (1,) + xc.shape
    return {'X': X, 'Y': Y, 'Z': Z, 'rho1': (P0 / (R * T)).reshape(sh),
            'u': (WALLS_U0 * (1.0 - np.exp(-yc / WALLS_DELTA))).reshape(sh), 'v': np.zeros(sh), 'w': np.zeros(sh),
            'p': np.full(sh, P0)}


DUCT_BANDS = [(0.64e-3, 40), (1.6e-3, 20), (4.0e-3, 20), (10.0e-3, 30)]
DUCT_P0, DUCT_PB, DUCT_TP, DUCT_EPS = 1.2e5, 1.0e5, 6.4e-5, 0.1
DUCT_TEND = 4.0e-4          # the tables cover 0 <= t <= DUCT_TEND (the order runs end at 64 H0, about 2.6e-4 s)
DUCT_M405 = 1.5             # Mach number of the supersonic inflow of mr1d-405


def duct_state(R, gam):
    """Steady isentropic state of the duct (inlet T0, p0 = DUCT_P0; outlet DUCT_PB): rho, u, p, c."""
    mach = math.sqrt(2.0 / (gam - 1.0) * ((DUCT_P0 / DUCT_PB) ** ((gam - 1.0) / gam) - 1.0))
    t = T0 / (1.0 + 0.5 * (gam - 1.0) * mach ** 2)
    c = math.sqrt(gam * R * t)
    return DUCT_PB / (R * t), mach * c, DUCT_PB, c


def duct_case(d, inlet, R, gam, closed=False):
    """mr1d-402 (inlet = 402), mr1d-410 (inlet = 410) and the ducts with the inlets 401, 403, 405, 408, 420 and the SRM
    grain of mode 2 (502): block, boundary file and, for 402 and 410, the time table in the case directory; closed
    (mr1d-502c): the grain with n = 0, the gas at rest and a slip wall (300) instead of the outlet."""
    widths = []
    for h, n in DUCT_BANDS:
        widths += [h] * n
    xn = np.concatenate([[0.0], np.cumsum(widths)])
    X, Y, Z = slab_1d(xn)
    n = xn.size - 1
    rho, u, p, c = duct_state(R, gam)
    sh = (1, 1, n)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.full(sh, rho), 'u': np.full(sh, 0.0 if closed else u), 'v': np.zeros(sh),
         'w': np.zeros(sh), 'p': np.full(sh, p)}
    head = '%9d %7d %7d %7d %7d %7d\n'
    if inlet == 402:
        # table read by Series_Data_m (values in Pa); file name without '/': a slash ends a list-directed read
        nstep = int(round(DUCT_TEND / DUCT_TP * 2000))
        with open(os.path.join(d, 'p0.dat'), 'w') as f:
            for k in range(nstep + 1):
                t = DUCT_TEND * k / nstep
                f.write('%.15E %.15E\n' % (t, DUCT_P0 * (1.0 + DUCT_EPS * math.sin(2 * math.pi * t / DUCT_TP))))
        out = [head % (1, 1, 1, 1, 1, 402) + '   %.1f, p0.dat, normal, normal, 1.0, 1.0\n' % T0]
    elif inlet != 410:
        # the other inlets on face 1, with the data of the steady state (static temperature t, mass flux g)
        t, g = p / (rho * R), rho * u
        data = {401: '%.1f, %.6E, normal, normal, 1.0, 1.0' % (T0, DUCT_P0),
                403: '%.1f, %.6E, normal, normal, 1.0, 1.0' % (T0, g),
                405: '%.6E, %.6E, %.6E, normal, normal, 1.0, 1.0' % (DUCT_M405, t, p),
                408: '%.6E, %.6E, normal, normal, 1.0, 1.0' % (t, u),
                420: '%.1f, %.6E, %.6E, %.6E, %.6E, 1.0, 1.0' % (T0, DUCT_P0, 0.5 * DUCT_PB, 0.25 * DUCT_PB, g),
                # mdot_ref, n, p_ref, T_af, mode, massf: Saint-Robert's law gives the steady mass flux at the steady
                # pressure, and the products enter with the total temperature T0 (mode 2: an inflow)
                502: '%.6E, %.1f, %.6E, %.1f, 2.0, 1.0' % (g, 0.0 if closed else 0.5, p, T0)}[inlet]
        out = [head % (1, 1, 1, 1, 1, inlet) + '   ' + data + '\n']
    else:
        # Q2D file (IO_BC_Q2D.f90): one zone per instant, title B<block>_T<n>, I = J = 2 nodes of the face patch
        # (face 1: (j, k) = (1, 1)), node coordinates x, y, z, then one cell value per variable
        nstep = int(round(DUCT_TEND / DUCT_TP * 400))
        xs = ' '.join('%.6E' % v for v in (0.0, 0.0, 0.0, 0.0))
        ys = ' '.join('%.6E' % v for v in (0.0, 0.01, 0.0, 0.01))
        zs = ' '.join('%.6E' % v for v in (-0.5 * DZ, -0.5 * DZ, 0.5 * DZ, 0.5 * DZ))
        with open(os.path.join(d, 'q2d.dat'), 'w') as f:
            f.write('VARIABLES = "x" "y" "z" "rho(1)" "u" "v" "w" "p"\n')
            for k in range(nstep + 1):
                t = DUCT_TEND * k / nstep
                s = DUCT_EPS * math.sin(2 * math.pi * t / DUCT_TP)
                f.write('ZONE T="B1_T%d", I=2, J=2, SOLUTIONTIME=%.15E\n' % (k + 1, t))
                f.write('%s\n%s\n%s\n' % (xs, ys, zs))
                for v in (rho * (1.0 + s / gam), u + c * s / gam, 0.0, 0.0, p * (1.0 + s)):
                    f.write('%.15E\n' % v)
        out = [head % (1, 1, 1, 1, 1, 410) + 'q2d.dat\n']
    out.append(head % (1, n, 1, 1, 2, 300) if closed else head % (1, n, 1, 1, 2, 406) + '%.6E, 1.0\n' % DUCT_PB)
    for f in (3, 4):
        for i in range(1, n + 1):
            out.append(head % (1, i, 1, 1, f, 300))
    for f in (5, 6):
        for i in range(1, n + 1):
            out.append(head % (1, i, 1, 1, f, 0))
    return B, ''.join(out)


ONE_LINE = {101, 103, 201, 301, 302, 501, 502, 503, 504, 505, 506, 410, 420} | set(range(401, 409))


def split_text(ic, bc, PI, PJ):
    """Single 2D block (K = 1 cell) of a MOSE case split into PI x PJ blocks (cuts in the middle of i and of j), on the
    text of INPUT/ic.tec (one zone, DATAPACKING=BLOCK; the numbers are copied as strings, so the values are identical
    to the bit) and of INPUT/bc.txt (every record goes to the piece holding its cell, with local indices and its data
    line; the new faces become pairs of connections 101 with the identity orientation 1 0 0 1). Blocks numbered
    b = 1 + p + PI q; records ordered by block, face, k, j, i."""
    import re
    txt = ic.split('\n')
    zi = [i for i, l in enumerate(txt) if l.lstrip().upper().startswith('ZONE')]
    assert len(zi) == 1, 'one zone required'
    head, zone = txt[:zi[0]], txt[zi[0]]
    dims = {k: int(v) for k, v in re.findall(r'\b([IJK])\s*=\s*(\d+)', zone)}
    I, J, K = dims['I'], dims['J'], dims['K']
    loc = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED\)', zone)
    nnod = int(loc.group(2)) - int(loc.group(1)) + 1
    ncel = int(loc.group(4)) - int(loc.group(3)) + 1
    toks = ' '.join(txt[zi[0] + 1:]).split()
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    assert len(toks) == nnod * nn + ncel * nc, 'number of values'
    nod = [toks[v * nn:(v + 1) * nn] for v in range(nnod)]
    cel = [toks[nnod * nn + v * nc: nnod * nn + (v + 1) * nc] for v in range(ncel)]
    ni, nj, nk = I - 1, J - 1, max(K - 1, 1)
    ci = [round(ni * k / PI) for k in range(PI + 1)]
    cj = [round(nj * k / PJ) for k in range(PJ + 1)]
    out = list(head)

    def bid(p, q):
        return 1 + p + PI * q
    for q in range(PJ):
        for p in range(PI):
            i0, i1, j0, j1 = ci[p], ci[p + 1], cj[q], cj[q + 1]
            z = re.sub(r'\bI\s*=\s*\d+', 'I=%d' % (i1 - i0 + 1), zone)
            z = re.sub(r'\bJ\s*=\s*\d+', 'J=%d' % (j1 - j0 + 1), z)
            z = re.sub(r'T\s*=\s*[^,]+,', 'T = B%d-IG,' % bid(p, q), z, count=1)
            out.append(z)
            for v in range(nnod):
                for k in range(K):
                    for j in range(j0, j1 + 1):
                        out.extend(nod[v][k * I * J + j * I + i0: k * I * J + j * I + i1 + 1])
            for v in range(ncel):
                for k in range(nk):
                    for j in range(j0, j1):
                        out.extend(cel[v][k * ni * nj + j * ni + i0: k * ni * nj + j * ni + i1])
    lines = bc.split('\n')
    recs, p0 = [], 0
    while p0 < len(lines):
        s = lines[p0].split()
        if not s:
            p0 += 1
            continue
        b, i, j, k, f, t = map(int, s[:6])
        p0 += 1
        prop = None
        if t in ONE_LINE:
            prop = lines[p0]
            p0 += 1
        assert t not in (101, 102, 201), 'connections in the single block are not handled'
        recs.append((i, j, k, f, t, prop))
    assert len(recs) == 2 * (nj * nk + ni * nk + ni * nj), 'number of records'
    new = {}
    for (i, j, k, f, t, prop) in recs:
        p = next(p for p in range(PI) if ci[p] < i <= ci[p + 1])
        q = next(q for q in range(PJ) if cj[q] < j <= cj[q + 1])
        new.setdefault(bid(p, q), []).append((i - ci[p], j - cj[q], k, f, t, prop))
    for q in range(PJ):
        for p in range(PI):
            b = bid(p, q)
            nip, njq = ci[p + 1] - ci[p], cj[q + 1] - cj[q]
            for k in range(1, nk + 1):
                if p + 1 < PI:            # cut in i: face 2 of (p, q) <-> face 1 of (p + 1, q)
                    for jl in range(1, njq + 1):
                        new[b].append((nip, jl, k, 2, 101, '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (bid(p + 1, q), 1, jl, k, 1, 1, 0, 0, 1)))
                        new[bid(p + 1, q)].append((1, jl, k, 1, 101, '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (b, nip, jl, k, 2, 1, 0, 0, 1)))
                if q + 1 < PJ:            # cut in j: face 4 of (p, q) <-> face 3 of (p, q + 1)
                    for il in range(1, nip + 1):
                        new[b].append((il, njq, k, 4, 101, '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (bid(p, q + 1), il, 1, k, 3, 1, 0, 0, 1)))
                        new[bid(p, q + 1)].append((il, 1, k, 3, 101, '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (b, il, njq, k, 4, 1, 0, 0, 1)))
    bcout = []
    for b in sorted(new):
        for (i, j, k, f, t, prop) in sorted(new[b], key=lambda r: (r[3], r[2], r[1], r[0])):
            bcout.append('%9d%8d%8d%8d%8d%8d' % (b, i, j, k, f, t))
            if prop is not None:
                bcout.append(prop)
    return '\n'.join(out) + '\n', '\n'.join(bcout) + '\n'


def turn_block_text(ic, bc, b):
    """Block b of a multi-block 2D case (ic.tec and bc.txt as split_text writes them: one zone per block, DATAPACKING=BLOCK)
    turned by 180 degrees about z: its i and j reversed, k unchanged, so the block stays right-handed. The nodes and the cell
    values are moved as text (the same numbers to the bit); the records of block b take the indices (ni + 1 - i, nj + 1 - j)
    and the faces 1 <-> 2, 3 <-> 4 (5 and 6 unchanged); in every connection record (101, 201) whose cell or partner is in
    block b the partner indices and face follow the turn, and the orientation d11 d12 d21 d22 (row r: the offsets in the
    partner, along its two tangential directions, of a step along the tangential direction r of the face of the record;
    Fill_BC_Ghost_Connection in Lib_Ghost.f90) changes sign in the rows of the directions i and j of block b (the cell of the
    record in block b) and in the columns of the directions i and j of block b (the partner in block b). The tangential
    directions of a face: j and k for faces 1-2, i and k for 3-4, i and j for 5-6. Records ordered by block, face, k, j, i."""
    import re
    txt = ic.split('\n')
    zi = [n for n, l in enumerate(txt) if l.lstrip().upper().startswith('ZONE')]
    z0 = zi[b - 1]
    z1 = zi[b] if b < len(zi) else len(txt)
    zone = txt[z0]
    dims = {k: int(v) for k, v in re.findall(r'\b([IJK])\s*=\s*(\d+)', zone)}
    I, J, K = dims['I'], dims['J'], dims['K']
    loc = re.search(r'VARLOCATION=\(\[(\d+)-(\d+)\]=NODAL,\[(\d+)-(\d+)\]=CELLCENTERED\)', zone)
    nnod = int(loc.group(2)) - int(loc.group(1)) + 1
    ncel = int(loc.group(4)) - int(loc.group(3)) + 1
    toks = ' '.join(txt[z0 + 1:z1]).split()
    ni, nj, nk = I - 1, J - 1, max(K - 1, 1)
    nn, nc = I * J * K, ni * nj * nk
    assert len(toks) == nnod * nn + ncel * nc, 'number of values of block %d' % b
    out = []
    for v in range(nnod):
        a = toks[v * nn:(v + 1) * nn]
        out += [a[k * I * J + (J - 1 - j) * I + (I - 1 - i)] for k in range(K) for j in range(J) for i in range(I)]
    for v in range(ncel):
        a = toks[nnod * nn + v * nc:nnod * nn + (v + 1) * nc]
        out += [a[k * ni * nj + (nj - 1 - j) * ni + (ni - 1 - i)] for k in range(nk) for j in range(nj) for i in range(ni)]
    txt = txt[:z0 + 1] + out + txt[z1:]
    while txt and txt[-1] == '':
        txt.pop()
    fmap = {1: 2, 2: 1, 3: 4, 4: 3, 5: 5, 6: 6}
    tang = {1: 'jk', 2: 'jk', 3: 'ik', 4: 'ik', 5: 'ij', 6: 'ij'}
    lines = bc.split('\n')
    recs, p0 = [], 0
    while p0 < len(lines):
        t = lines[p0].split()
        if not t:
            p0 += 1
            continue
        rb, i, j, k, f, ty = map(int, t[:6])
        p0 += 1
        prop = None
        if ty in ONE_LINE:
            prop = lines[p0]
            p0 += 1
        if ty in (101, 201):
            bs, is_, js, ks, fs, d11, d12, d21, d22 = map(int, prop.split())
            d = [[d11, d12], [d21, d22]]
            if rb == b:
                d = [[-x for x in row] if tang[fmap[f]][r] in 'ij' else row for r, row in enumerate(d)]
            if bs == b:
                is_, js, fs = ni + 1 - is_, nj + 1 - js, fmap[fs]
                d = [[-row[c] if tang[fs][c] in 'ij' else row[c] for c in range(2)] for row in d]
            prop = '%8d%8d%8d%8d%8d%8d%8d%8d%8d' % (bs, is_, js, ks, fs, d[0][0], d[0][1], d[1][0], d[1][1])
        elif ty in (102, 104):
            raise SystemExit('turn_block_text: chimera records are not handled')
        if rb == b:
            i, j, f = ni + 1 - i, nj + 1 - j, fmap[f]
        recs.append((rb, f, k, j, i, ty, prop))
    bcout = []
    for (rb, f, k, j, i, ty, prop) in sorted(recs, key=lambda r: r[:5]):
        bcout.append('%9d%8d%8d%8d%8d%8d' % (rb, i, j, k, f, ty))
        if prop is not None:
            bcout.append(prop)
    return '\n'.join(txt) + '\n', '\n'.join(bcout) + '\n'


def plate_case(d, air, PI, PJ, mu=None):
    """mr2d-plate (PI = PJ = 1) and its splits: the laminar flat plate of test/2D/viscous/flat-plate-laminar (64 x 64
    cells, x from -0.061 to 0.305 m along j, y from 0.03 m (i = 1) to the wall (i = 64), cells from 16 um at the wall to
    2.6 mm; uniform start at Mach 0.2; inlet 404 on face 3, ambient 407 on faces 1 and 4, symmetry 300 for x < 0 and
    adiabatic wall 301 for x > 0 on face 2, k faces 200 on the planar slab 1 m thick), air of test/common/Air-viscous
    (mu = 1.1858686e-5 Pa s, k = 0.016545 W/m/K, cp = 1004.5 J/kg/K). Files copied as text (split_text for the
    splits); the input.ini of the multirate tests (INI_VISC). With mu given (mr2d-plate-mu), a constant transport table
    with that viscosity and Pr = 1 instead, over the temperature range of the thermo table (1 ... 5000 K): mu = 1e-3 Pa s
    is the viscosity of the [GPB-Phase1] block of the input.ini of that test (read by its case generator, not by MOSE),
    with which the cell Reynolds number at the wall is about 6 and the explicit steady run with the local step reaches
    the steady state in a few 1e5 iterations (with the air of Air-viscous it would take about 1e7)."""
    import re
    common = os.path.dirname(os.path.abspath(air))
    src = os.path.join(common, '..', '2D', 'viscous', 'flat-plate-laminar', 'INPUT')
    ic = open(os.path.join(src, 'ic.tec')).read()
    bc = open(os.path.join(src, 'bc.txt')).read()
    # the VARIABLES line of that file has the names without quotes, which the Tecplot reader of ORION does not count
    # (MOSE stops: "Number of variables in IOfield does not match"): names quoted, numbers unchanged
    lines = ic.split('\n')
    for n, ln in enumerate(lines):
        if ln.lstrip().upper().startswith('VARIABLES'):
            pre, _, rest = ln.partition('=')
            lines[n] = pre + '=' + ' '.join(t if t.startswith('"') else '"%s"' % t
                                            for t in re.findall(r'"[^"]*"|\S+', rest))
            break
    ic = '\n'.join(lines)
    if (PI, PJ) != (1, 1):
        ic, bc = split_text(ic, bc, PI, PJ)
    open(os.path.join(d, 'INPUT', 'ic.tec'), 'w').write(ic)
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(bc)
    for f in ('phase.txt', 'thermo.dat', 'transport.dat'):
        if f == 'transport.dat' and mu is not None:
            cp = float(np.loadtxt(os.path.join(common, 'Air-viscous', 'thermo.dat'), skiprows=4)[0, 1])
            with open(os.path.join(d, 'INPUT', f), 'w') as fo:
                fo.write('TITLE = "Transport Properties"\nVARIABLES = "Temperature", "Viscosity", "Conductivity"\n'
                         'ZONE T="Air"\nI=5000, F=POINT\n')
                fo.write(''.join('%.1f %.10e %.10e\n' % (t, mu, mu * cp) for t in range(1, 5001)))
            continue
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(common, 'Air-viscous', f), 'rb').read())
    open(os.path.join(d, 'input.ini'), 'w').write(INI_VISC.format(iters=100, cfl=0.5, vnn=0.5))


# viscous cases: Navier-Stokes with a laminar Schmidt number of 0.7 (with one species the mass diffusion flux is zero
# whatever Sc)
INI_VISC = """[MOSE-Parameters]
iter-threshold = {iters}
newrun = true

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7

[MOSE-Numerics]
cfl = {cfl}
vnn = {vnn}
cfl-rise-threshold = 0
time-scheme = RK2
time-accurate = true
irs-beta = 0.0
riemann-solver = HLLC
flux-limiter = vanleer
space-reconstruction = MUSCL
integration-variables = cons

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
res-diter = 1000000000
"""


INI = """[MOSE-Parameters]
iter-threshold = {iters}
newrun = true

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = {cfl}
cfl-rise-threshold = 0
time-scheme = RK2
time-accurate = true
irs-beta = 0.0
riemann-solver = HLLC
flux-limiter = vanleer
space-reconstruction = MUSCL
integration-variables = cons

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
res-diter = 1000000000
"""



# ---- Spalart-Allmaras and SST, and a 3D viscous case ------------------------------------------------------------------

# input.ini of the RANS cases: Navier-Stokes with the turbulence model; omega-wall-bc read by SST only
INI_RANS = """[MOSE-Parameters]
iter-threshold = {iters}
newrun = true

[MOSE-Physics]
equations = navier-stokes
turbulence = {model}
schmidt = 0.7

[MOSE-Turbulence]
omega-wall-bc = {wall}

[MOSE-Numerics]
cfl = {cfl}
vnn = {vnn}
cfl-rise-threshold = 0
time-scheme = RK2
time-accurate = true
irs-beta = 0.0
riemann-solver = HLLC
flux-limiter = vanleer
space-reconstruction = MUSCL
integration-variables = cons

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
res-diter = 1000000000
"""


def copy_viscous_air(d, air):
    """phase.txt, thermo.dat and transport.dat of test/common/Air-viscous (mu = 1.1858686e-5 Pa s, k = 0.016545 W/m/K,
    cp = 1004.5 J/kg/K), copied byte by byte"""
    common = os.path.dirname(os.path.abspath(air))
    for f in ('phase.txt', 'thermo.dat', 'transport.dat'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(common, 'Air-viscous', f), 'rb').read())


def turb_plate_case(d, air, model, PI, PJ):
    """mr2d-plate-sa, mr2d-plate-sst and the split mr2d-plate-sa-4b and mr2d-plate-sst-4b: the turbulent flat plate of
    mr-common/plate on its coarse mesh (68 x 48 cells; x from -0.333 to 2 m, y from the wall to 1 m, cells from 4 um
    at the wall), uniform start at Mach 0.2 and 555 K; inlet 404, outlet 406, symmetry 300 for x < 0 and on the top
    face, adiabatic wall 301 for x > 0, null k faces (0). SST: ic-SST.tec and bc-SST.txt. SA: the initial field
    ic-SA.tec, on the same mesh, with the records of bc-SST.txt and in the inlet records the SA value rho nu~ = 1e-6 in
    place of k and omega. The three files are those of the SST and SA cases of the turbulent flat plate with boundary
    conditions of its own (test/2D/viscous/flat-plate-turbulent, since removed), byte for byte. Text copied;
    SOLUTIONTIME of the initial field set to 0 (the files carry the time of the run that wrote them). Air of
    test/common/Air-viscous (the gas files of that test are the same bytes). cfl 0.5, vnn 0.4; SST with
    omega-wall-bc = asymptotic, as that test."""
    import re
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plate')
    bc = open(os.path.join(src, 'bc-SST.txt')).read()
    if model == 'SA':
        ic = open(os.path.join(src, 'ic-SA.tec')).read()
        old = '    0.225000E-02,    0.625000E+04,'
        n404 = len(re.findall(r'^\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+404\s*$', bc, re.M))
        assert bc.count(old) == n404 > 0, 'inlet records of the SST plate'
        bc = bc.replace(old, '    0.100000E-05,')
    else:
        ic = open(os.path.join(src, 'ic-SST.tec')).read()
    ic, n = re.subn(r'SOLUTIONTIME=\s*[-+0-9.EeDd]+', 'SOLUTIONTIME=0.000000000000000E+000', ic, count=1)
    assert n == 1, 'SOLUTIONTIME'
    if (PI, PJ) != (1, 1):
        ic, bc = split_text(ic, bc, PI, PJ)
    open(os.path.join(d, 'INPUT', 'ic.tec'), 'w').write(ic)
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(bc)
    copy_viscous_air(d, air)
    open(os.path.join(d, 'input.ini'), 'w').write(INI_RANS.format(iters=100, model=model, wall='asymptotic', cfl=0.5,
                                                                   vnn=0.4))


# wall channel of the test TimeIntegration/TurbulenceSources
# (test/checks/time-integration/turbulence-sources/make_case.py): 2-D,
# periodic along the flow (4 cells of 1 mm, connection 201 on itself), adiabatic no-slip wall (301) at y = 0, symmetry
# (300) at the top, null k faces (0), slab 1 mm; air of test/common/Air-viscous at 300 K and 1e5 Pa, u = U tanh(y/delta_u),
# no force that drives it: the wall slows the flow down and the turbulence variables relax. Constants of that file.
CH_NI, CH_DX, CH_DZ = 4, 1.0e-3, 1.0e-3
CH_MU = 1.1858686e-5
CH_T0, CH_P0, CH_U, CH_DELTA_U = 300.0, 1.0e5, 50.0, 2.0e-4
CH_KAPPA, CH_BETA1, CH_BETASTAR = 0.41, 0.075, 0.09


def chan_case(d, air, model, kind, ks=0.0):
    """mr2d-chan-sa, mr2d-chan-sst (kind 'order') and mr2d-chanpos-sa, mr2d-chanpos-sst (kind 'pos'): the wall channel
    of the test TimeIntegration/TurbulenceSources with the profiles of that test. mr2d-chanrough-sa: mr2d-chan-sa
    with a rough wall, sand-grain height ks = 10 um (SA-rough, Aupoix and Spalart 2003, as in Lib_Spalart.f90).
    order: SA nu~ = 2 kappa u_tau y exp(-y/0.2 mm), twice the near-wall equilibrium; SST k ~ y^2 (peak at 0.1 mm),
    omega = 6 nu/(beta_1 y^2) plus a log-layer part u_tau/(sqrt(beta*) kappa y), omega-wall-bc = asymptotic. Mesh: the
    16 wall rows of that test, from 4 um growing by 1.25 (symmetry plane at 0.553 mm, inside the shear layer and the
    turbulent layer). Left out: the thin cells of that test at the symmetry plane, which there keep the single-rate step
    independent of the turbulence (the order runs of the multirate use fixed macro steps), and the 8 outer rows of
    120 um, a free stream where u does not change: there the error of the single-rate reference at dt_min/64 is its own
    round-off (3e-11 m/s, the same against the single rate at dt_min/2 and dt_min/4), which the volume-weighted L1 of u
    picks up. With the finest cells at the wall the levels change every few rows inside
    the turbulent layer.
    pos: the positivity cases of that test, mesh and profiles unchanged: 18 rows from 2 um growing by 1.25, then 16 rows
    of 100 um; SA nu~ a hundred times the equilibrium within 20 um of the wall; SST omega ten times the asymptotic value
    within 20 um, omega-wall-bc = practical; the destruction acts faster than the transport, and only the stiffness in
    the local step keeps the stages positive.
    The turbulence fades out as exp(-(y/0.4 mm)^4): beyond 1 mm nu~, k and omega take uniform free-stream values
    (rho nu~ = rho k = 2e-10, twice the realizability floor; omega = 1e-2 1/s). cfl 0.8, vnn 0.4, as that test."""
    R, _ = gas(air)
    rho0 = CH_P0 / (R * CH_T0)
    nu = CH_MU / rho0
    utau = math.sqrt(nu * CH_U / CH_DELTA_U)
    if kind == 'order':
        sizes = [4.0e-6 * 1.25 ** k for k in range(16)]
    else:
        sizes = [2.0e-6 * 1.25 ** k for k in range(18)] + [1.0e-4] * 16
    y = np.concatenate(([0.0], np.cumsum(sizes)))
    nj = len(sizes)
    yc = 0.5 * (y[1:] + y[:-1])
    fade = np.exp(-(yc / 4.0e-4) ** 4)
    u = CH_U * np.tanh(yc / CH_DELTA_U)
    if model == 'SA':
        if kind == 'order':
            nut = 2.0 * CH_KAPPA * utau * yc * np.exp(-yc / 2.0e-4) * fade + 2.0e-10 / rho0
        else:
            nut = (100.0 * CH_KAPPA * utau * yc * np.exp(-yc / 2.0e-5) + 2.0 * CH_KAPPA * utau * yc * np.exp(-yc / 2.0e-4)) \
                * fade + 2.0e-10 / rho0
        rans = {'mi_tilde': rho0 * nut}
        wall = 'asymptotic'
    else:
        ka = utau ** 2 / math.sqrt(CH_BETASTAR)
        k = ka * (yc / 1.0e-4) ** 2 * np.exp(2.0 * (1.0 - yc / 1.0e-4)) * fade + 2.0e-10 / rho0
        w10 = 1.0 if kind == 'order' else 10.0
        om = (6.0 * nu / (CH_BETA1 * yc ** 2) * (1.0 + (w10 - 1.0) * np.exp(-yc / 2.0e-5))
              + utau / (math.sqrt(CH_BETASTAR) * CH_KAPPA * yc)) * fade + 1.0e-2
        rans = {'kappa': rho0 * k, 'omega': rho0 * om}
        wall = 'asymptotic' if kind == 'order' else 'practical'
    ni = CH_NI
    xs = np.arange(ni + 1) * CH_DX
    X = np.empty((2, nj + 1, ni + 1)); Y = np.empty_like(X); Z = np.empty_like(X)
    for kk in range(2):
        for j in range(nj + 1):
            X[kk, j] = xs; Y[kk, j] = y[j]; Z[kk, j] = kk * CH_DZ
    sh = (1, nj, ni)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.full(sh, rho0), 'u': np.repeat(u[:, None], ni, axis=1).reshape(sh),
         'v': np.zeros(sh), 'w': np.zeros(sh), 'p': np.full(sh, CH_P0)}
    names = ['rho1', 'u', 'v', 'w', 'p']
    for nm, val in rans.items():
        B[nm] = np.repeat(val[:, None], ni, axis=1).reshape(sh)
        names.append(nm)
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], names)
    out = []
    for j in range(1, nj + 1):
        out.append(rec(1, 1, j, 1, 1, 201, (1, ni, j, 1, 2, 1, 0, 0, 1)))
    for j in range(1, nj + 1):
        out.append(rec(1, ni, j, 1, 2, 201, (1, 1, j, 1, 1, 1, 0, 0, 1)))
    for i in range(1, ni + 1):
        out.append(rec(1, i, 1, 1, 3, 301) + '    0.000000E+00,    %.6E,    0.000000E+00,\n' % ks)
    for i in range(1, ni + 1):
        out.append(rec(1, i, nj, 1, 4, 300))
    for f in (5, 6):
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                out.append(rec(1, i, j, 1, f, 0))
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(''.join(out))
    copy_viscous_air(d, air)
    open(os.path.join(d, 'input.ini'), 'w').write(INI_RANS.format(iters=100, model=model, wall=wall, cfl=0.8, vnn=0.4))


HOT3_N, HOT3_H = 12, 0.01


def hot3d_case(d, air, cut=False):
    """mr3d-hot (cut: mr3d-hot-8b, in 2 x 2 x 2 blocks of 6 x 6 x 6 cells): the 3D twin of mr2d-hot: 12 x 12 x 12 cells of 10 mm on an oblique affine mesh x = xi + 0.5 eta + 0.3 zeta, y = eta + 0.4 zeta,
    z = zeta, periodic (201, the block connected to itself) on the six faces; uniform flow u = 100, v = 50, w = 30 m/s
    at p0, T = T0 (1 + 7 f) with f = (Sx^-8 + Sy^-8 + Sz^-8)^(-1/8), the smooth minimum of three plateaus (the plateau
    of mr2d-hot in each direction): a hot cube where the speed of sound doubles. With sigma = 2 and L = 1 the level map
    is a cube of level 1 with sharp edges and vertices, where the diffusive faces read coarser cells that only the
    diagonals reach, in k too. Viscosity HOT_MU = 0.3 Pa s, Pr 0.72; air of test/common/Air; cfl 0.2, vnn 0.5."""
    R, gam = gas(air)
    n, h = HOT3_N, HOT3_H
    q = np.arange(n + 1) * h
    ZE, ET, XI = np.meshgrid(q, q, q, indexing='ij')          # (K, J, I)
    X = XI + 0.5 * ET + 0.3 * ZE
    Y = ET + 0.4 * ZE
    Z = ZE.copy()

    def plateau(s):
        s = np.abs(s)
        return np.where(s <= 0.1, 1.0, np.where(s >= 0.15, 0.0, 0.5 * (1.0 + np.cos(np.pi * (s - 0.1) / 0.05))))
    c = (np.arange(n) + 0.5) / n - 0.5
    SK, SJ, SI = np.meshgrid(plateau(c), plateau(c), plateau(c), indexing='ij')
    inside = (SI > 0.0) & (SJ > 0.0) & (SK > 0.0)
    f = np.zeros(SI.shape)
    f[inside] = (SI[inside] ** -8 + SJ[inside] ** -8 + SK[inside] ** -8) ** (-1.0 / 8.0)
    T = T0 * (1.0 + (HOT_TR - 1.0) * f)
    sh = (n, n, n)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': P0 / (R * T), 'u': np.full(sh, 100.0), 'v': np.full(sh, 50.0),
         'w': np.full(sh, 30.0), 'p': np.full(sh, P0)}
    out = []
    if cut:
        m = n // 2
        pieces = []

        def bid(p, q, r):
            return 1 + p % 2 + 2 * (q % 2) + 4 * (r % 2)
        for r in range(2):
            for q in range(2):
                for p in range(2):
                    b = bid(p, q, r)
                    nod = (slice(r * m, r * m + m + 1), slice(q * m, q * m + m + 1), slice(p * m, p * m + m + 1))
                    cel = (slice(r * m, r * m + m), slice(q * m, q * m + m), slice(p * m, p * m + m))
                    piece = {nm: B[nm][nod] for nm in ('X', 'Y', 'Z')}
                    piece.update({nm: B[nm][cel] for nm in ('rho1', 'u', 'v', 'w', 'p')})
                    pieces.append(piece)
                    # a face inside the domain is a connection 101, a face on the end of the domain a periodic 201
                    tm = (101 if p == 1 else 201, 101 if p == 0 else 201, 101 if q == 1 else 201, 101 if q == 0 else 201,
                          101 if r == 1 else 201, 101 if r == 0 else 201)
                    idx = range(1, m + 1)
                    out += [rec(b, 1, j, k, 1, tm[0], (bid(p - 1, q, r), m, j, k, 2, 1, 0, 0, 1)) for k in idx for j in idx]
                    out += [rec(b, m, j, k, 2, tm[1], (bid(p + 1, q, r), 1, j, k, 1, 1, 0, 0, 1)) for k in idx for j in idx]
                    out += [rec(b, i, 1, k, 3, tm[2], (bid(p, q - 1, r), i, m, k, 4, 1, 0, 0, 1)) for k in idx for i in idx]
                    out += [rec(b, i, m, k, 4, tm[3], (bid(p, q + 1, r), i, 1, k, 3, 1, 0, 0, 1)) for k in idx for i in idx]
                    out += [rec(b, i, j, 1, 5, tm[4], (bid(p, q, r - 1), i, j, m, 6, 1, 0, 0, 1)) for j in idx for i in idx]
                    out += [rec(b, i, j, m, 6, tm[5], (bid(p, q, r + 1), i, j, 1, 5, 1, 0, 0, 1)) for j in idx for i in idx]
        write_tec(os.path.join(d, 'INPUT', 'ic.tec'), pieces, ['rho1', 'u', 'v', 'w', 'p'])
    else:
        write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], ['rho1', 'u', 'v', 'w', 'p'])
        out += [rec(1, 1, j, k, 1, 201, (1, n, j, k, 2, 1, 0, 0, 1)) for k in range(1, n + 1) for j in range(1, n + 1)]
        out += [rec(1, n, j, k, 2, 201, (1, 1, j, k, 1, 1, 0, 0, 1)) for k in range(1, n + 1) for j in range(1, n + 1)]
        out += [rec(1, i, 1, k, 3, 201, (1, i, n, k, 4, 1, 0, 0, 1)) for k in range(1, n + 1) for i in range(1, n + 1)]
        out += [rec(1, i, n, k, 4, 201, (1, i, 1, k, 3, 1, 0, 0, 1)) for k in range(1, n + 1) for i in range(1, n + 1)]
        out += [rec(1, i, j, 1, 5, 201, (1, i, j, n, 6, 1, 0, 0, 1)) for j in range(1, n + 1) for i in range(1, n + 1)]
        out += [rec(1, i, j, n, 6, 201, (1, i, j, 1, 5, 1, 0, 0, 1)) for j in range(1, n + 1) for i in range(1, n + 1)]
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(''.join(out))
    for fn in ('phase.txt', 'thermo.dat'):
        open(os.path.join(d, 'INPUT', fn), 'wb').write(open(os.path.join(air, fn), 'rb').read())
    cp_air = float(np.loadtxt(os.path.join(air, 'thermo.dat'), skiprows=4)[0, 1])
    with open(os.path.join(d, 'INPUT', 'transport.dat'), 'w') as fo:
        fo.write('TITLE = "Transport Properties"\nVARIABLES = "Temperature", "Viscosity", "Conductivity"\n'
                 'ZONE T="Air"\nI=10000, F=POINT\n')
        fo.write(''.join('%.1f %.10e %.10e\n' % (t, HOT_MU, HOT_MU * cp_air / 0.72) for t in range(1, 10001)))
    open(os.path.join(d, 'input.ini'), 'w').write(INI_VISC.format(iters=100, cfl=0.2, vnn=0.5))


TURB_PLATES = {'mr2d-plate-sa': ('SA', 1, 1), 'mr2d-plate-sst': ('SST', 1, 1), 'mr2d-plate-sa-4b': ('SA', 2, 2),
               'mr2d-plate-sst-4b': ('SST', 2, 2)}
CHANNELS = {'mr2d-chan-sa': ('SA', 'order'), 'mr2d-chan-sst': ('SST', 'order'),
            'mr2d-chanpos-sa': ('SA', 'pos'), 'mr2d-chanpos-sst': ('SST', 'pos'),
            'mr2d-chanrough-sa': ('SA', 'order', 1.0e-5)}


PLATES = {'mr2d-plate': (1, 1), 'mr2d-plate-2b': (1, 2), 'mr2d-plate-4b': (2, 2), 'mr2d-plate-mu': (1, 1, 1.0e-3)}


def wedge_case(d, air):
    """mr2d-wedge-visc: the polynomial field with dilatation and swirl of part 4 of the test MMSResiduals on the 2-D
    axisymmetric wedge of MOSE (1 degree, k faces of type 200, axis at j = 1 with type 300, 400 on the other faces;
    32 x 32 cells, x and r in [0, 3 mm]), written by its own generator (test/2D/viscous/mms/residuals/wedge_viscous.py,
    write_case), so that the multirate with one level is checked against the single-rate step on the viscous stress of
    the wedge faces; air of test/common/Air-viscous (nu = 1e-2 m2/s at the density of the case); time-accurate
    Heun, cfl 0.5, vnn 0.3."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'wedge_case', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2D', 'viscous', 'mms',
                                   'residuals', 'wedge_viscous.py'))
    wv = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wv)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        wv.write_case('swirl', 32)
    finally:
        os.chdir(cwd)
    copy_viscous_air(d, air)
    open(os.path.join(d, 'input.ini'), 'w').write(INI_VISC.format(iters=100, cfl=0.5, vnn=0.3))


WB_NX, WB_HX = 24, 1.0e-3                                              # mr2d-wedge-band: columns in x
WB_RBANDS = [(0.125e-3, 8), (0.25e-3, 8), (0.5e-3, 8), (1.0e-3, 8)]    # rows in r from the axis: height, number
WB_DELTA = math.radians(1.0)                                           # angle of the wedge
WB_SWIRL = 10.0                                                        # swirl velocity at r_max [m/s]
PLATE_TW = 330.0                                                       # wall temperature of mr2d-plate-302 [K]


def wedge_band_case(d, air):
    """mr2d-wedge-band: the 2-D axisymmetric wedge with bands of rows in r (see the list of the cases). Nodes
    (x, r cos(delta/2), -+ r sin(delta/2)), as test/2D/viscous/mms/residuals/wedge_viscous.py writes its wedge; records
    by face 1 ... 6, as BC_Fluxes requires."""
    R, gam = gas(air)
    rn = [0.0]
    for h, n in WB_RBANDS:
        rn += [rn[-1] + h * (q + 1) for q in range(n)]
    rn[0] = 1.0e-6 * WB_RBANDS[0][0]
    rn = np.array(rn)
    xn = WB_HX * np.arange(WB_NX + 1)
    I, J = xn.size, rn.size
    c, s = math.cos(0.5 * WB_DELTA), math.sin(0.5 * WB_DELTA)
    X = np.empty((2, J, I)); Y = np.empty((2, J, I)); Z = np.empty((2, J, I))
    for k in range(2):
        X[k] = xn[None, :]
        Y[k] = rn[:, None] * c
        Z[k] = (-1.0, 1.0)[k] * rn[:, None] * s
    xc = 0.5 * (xn[1:] + xn[:-1]); rc = 0.5 * (rn[1:] + rn[:-1])
    XC, RC = np.meshgrid(xc, rc)
    rho0 = P0 / (R * T0)
    p = P0 * (1 + 1e-2 * np.exp(-((XC - 0.012) ** 2 + (RC - 0.006) ** 2) / 0.003 ** 2))
    rho = rho0 * (p / P0) ** (1 / gam)
    sh = (1,) + XC.shape
    z = np.zeros(sh)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': z, 'v': z,
         'w': (WB_SWIRL * RC / rn[-1]).reshape(sh), 'p': p.reshape(sh)}
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], ['rho1', 'u', 'v', 'w', 'p'])
    nx, nr = I - 1, J - 1
    out = [rec(1, 1, j, 1, 1, 300) for j in range(1, nr + 1)]
    out += [rec(1, nx, j, 1, 2, 300) for j in range(1, nr + 1)]
    out += [rec(1, i, 1, 1, 3, 300) for i in range(1, nx + 1)]
    out += [rec(1, i, nr, 1, 4, 300) for i in range(1, nx + 1)]
    for f in (5, 6):
        out += [rec(1, i, j, 1, f, 200) for j in range(1, nr + 1) for i in range(1, nx + 1)]
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(''.join(out))
    copy_viscous_air(d, air)
    open(os.path.join(d, 'input.ini'), 'w').write(INI_VISC.format(iters=100, cfl=0.5, vnn=0.3))


def wall_302(bc):
    """the records 301 (adiabatic wall: q_wall, k_rough, eps_wall) of a boundary file made 302 (isothermal wall: T_wall,
    k_rough, eps_wall) at PLATE_TW, the other values of the data line kept"""
    lines = bc.split('\n')
    out, k = [], 0
    while k < len(lines):
        t = lines[k].split()
        if len(t) == 6 and t[5] == '301':
            out.append(lines[k][:lines[k].rindex('301')] + '302')
            vals = lines[k + 1].split(',')
            out.append(','.join(['    %.6E' % PLATE_TW] + vals[1:]))
            k += 2
            continue
        out.append(lines[k])
        k += 1
    return '\n'.join(out)


DARU = {'mr2d-daru-coarse': (50, 25), 'mr2d-daru': (100, 50), 'mr2d-daru-full': (500, 250)}

INI_DARU = """[MOSE-Parameters]
iter-threshold = 100000000
time-threshold = {tfin!r}
newrun = true

[MOSE-Physics]
equations = navier-stokes
prandtl = 0.73
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.8
vnn = 0.3
cfl-rise-threshold = 0
time-scheme = RK2
time-accurate = true
irs-beta = 0.0
riemann-solver = HLLC
flux-limiter = vanleer
space-reconstruction = MUSCL
integration-variables = cons

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
sol-overwrite = false
shell-diter = 1000000000
ini-diter = 1000000000
res-diter = 1000000000
"""


NOZ_BANDS = [(0.0, 0.3, 0.01), (0.3, 0.38, 0.004), (0.38, 0.444, 0.0016), (0.444, 0.508, 0.00064),
             (0.508, 0.700, 0.0016), (0.700, 0.912, 0.004), (0.912, 0.992, 0.01)]
NOZ_XT, NOZ_AIN, NOZ_AOUT = 0.476, 3.0, 3.0          # throat [m]; area at the inlet and at the outlet over the throat
NOZ_RATIO, NOZ_TRAMP, NOZ_TEND = 2.6, 8.0e-3, 1.0e-2  # p0 at the end of the ramp over P0; ramp [s]; time-threshold [s]


def nozzle_case(d, air):
    """mr1d-nozzle-su: quasi-1D convergent-divergent nozzle started by a ramp of the total pressure (module docstring)."""
    R, gam = gas(air)
    xn = banded(NOZ_BANDS)
    a = np.where(xn < NOZ_XT, (NOZ_AIN - 1.0) / NOZ_XT ** 2, (NOZ_AOUT - 1.0) / (xn[-1] - NOZ_XT) ** 2)
    h = 0.01 * np.sqrt(1.0 + a * (xn - NOZ_XT) ** 2)
    X, Y, Z = slab_1d(xn)
    Y[:, 1, :] = h
    Z[1] = -0.5 * DZ + h
    n = xn.size - 1
    sh = (1, 1, n)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.full(sh, P0 / (R * T0)), 'u': np.zeros(sh), 'v': np.zeros(sh),
         'w': np.zeros(sh), 'p': np.full(sh, P0)}
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], ['rho1', 'u', 'v', 'w', 'p'])
    # table read by Series_Data_m (values in Pa), up to 1.2 times the time-threshold
    nstep = 2000
    with open(os.path.join(d, 'p0.dat'), 'w') as f:
        for k in range(nstep + 1):
            t = 1.2 * NOZ_TEND * k / nstep
            f.write('%.15E %.15E\n' % (t, P0 * (1.0 + (NOZ_RATIO - 1.0) * min(t / NOZ_TRAMP, 1.0))))
    out = [rec(1, 1, 1, 1, 1, 402) + '   %.1f, p0.dat, normal, normal, 1.0, 1.0\n' % T0,
           rec(1, n, 1, 1, 2, 406) + '%.6E, 1.0\n' % P0]
    out += [rec(1, i, 1, 1, f, 300) for f in (3, 4, 5, 6) for i in range(1, n + 1)]
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(''.join(out))
    for f in ('phase.txt', 'thermo.dat'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(air, f), 'rb').read())
    ini = INI.format(iters=100, cfl=0.8).replace('iter-threshold = 100\n',
                                                 'iter-threshold = 100\ntime-threshold = %r\n' % NOZ_TEND)
    open(os.path.join(d, 'input.ini'), 'w').write(ini)


def daru_case(d, air, nx, ny):
    """Viscous shock tube of Daru and Tenaud, the lower half [0; 1] x [0; 0.5] m on nx x ny cells: x uniform; y_j =
    0.5 (1 - tanh(beta (1 - j/ny)) / tanh(beta)) with beta such that the first cell at the wall y = 0 is dx / 4; one
    cell of DZ in z. Walls 301 (q = 0) at x = 0, x = 1 and y = 0, symmetry 300 at y = 0.5 and on the k faces. Air of
    test/common/Air (cp constant) with mu = rho_ref a_ref L / Re, Re = 200, rho_ref = rho_R / 1.2, a_ref = a(300 K),
    L = 1 m, and the conductivity of Pr = 0.73 (the key prandtl sets k = mu cp / Pr). T = 300 K everywhere, p = 100 P0
    for x < 0.5 and P0 beyond, at rest; time-threshold = L / a_ref, one unit of the non-dimensional time."""
    R, gam = gas(air)
    a_ref = math.sqrt(gam * R * T0)
    rho_r = P0 / (R * T0)
    mu = (rho_r / 1.2) * a_ref * 1.0 / 200.0
    cp = gam * R / (gam - 1.0)
    ly, dy1 = 0.5, 0.25 / nx

    def first_cell(b):
        return ly * (1.0 - math.tanh(b * (1.0 - 1.0 / ny)) / math.tanh(b))
    lo, hi = 1e-6, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if first_cell(mid) > dy1:
            lo = mid
        else:
            hi = mid
    beta = 0.5 * (lo + hi)
    jj = np.arange(ny + 1)
    yn = ly * (1.0 - np.tanh(beta * (1.0 - jj / ny)) / math.tanh(beta))
    yn[0], yn[-1] = 0.0, ly
    xn = np.linspace(0.0, 1.0, nx + 1)
    X = np.empty((2, ny + 1, nx + 1)); Y = np.empty_like(X); Z = np.empty_like(X)
    for k in range(2):
        X[k] = xn[None, :]
        Y[k] = yn[:, None]
        Z[k] = (-0.5 * DZ, 0.5 * DZ)[k]
    xc = 0.5 * (xn[1:] + xn[:-1])
    p = np.where(np.broadcast_to(xc[None, :] < 0.5, (ny, nx)), 100.0 * P0, P0)
    sh = (1, ny, nx)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': (p / (R * T0)).reshape(sh), 'u': np.zeros(sh), 'v': np.zeros(sh),
         'w': np.zeros(sh), 'p': p.reshape(sh)}
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], ['rho1', 'u', 'v', 'w', 'p'])
    wall = '   0.0, 0.0, 0.0\n'
    out = [rec(1, 1, j, 1, 1, 301) + wall for j in range(1, ny + 1)]
    out += [rec(1, nx, j, 1, 2, 301) + wall for j in range(1, ny + 1)]
    out += [rec(1, i, 1, 1, 3, 301) + wall for i in range(1, nx + 1)]
    out += [rec(1, i, ny, 1, 4, 300) for i in range(1, nx + 1)]
    out += [rec(1, i, j, 1, f, 300) for f in (5, 6) for j in range(1, ny + 1) for i in range(1, nx + 1)]
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(''.join(out))
    for f in ('phase.txt', 'thermo.dat'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(air, f), 'rb').read())
    with open(os.path.join(d, 'INPUT', 'transport.dat'), 'w') as fo:
        fo.write('TITLE = "Transport Properties"\nVARIABLES = "Temperature", "Viscosity", "Conductivity"\n'
                 'ZONE T="Air"\nI=10000, F=POINT\n')
        fo.write(''.join('%.1f %.12f %.12f \n' % (1.0 + k, mu, mu * cp / 0.73) for k in range(10000)))
    open(os.path.join(d, 'input.ini'), 'w').write(INI_DARU.format(tfin=1.0 / a_ref))


# Simple wave with moving level boundaries (test mr-moving-levels): c = 1 + DELTA phi, phi = (1 - cos 2 pi x) / 2,
# u = 2 (c - 1) / (gamma - 1), rho = c^(2 / (gamma - 1)), p = rho c^2 / gamma: the Riemann invariant J- = u - 2 c / (gamma - 1) is
# the same everywhere, so the wave is simple and c is constant along the characteristics x = x0 + (u + c) t, for the gamma of the
# gas the solver computes (with gamma = 1.4, u = 5 DELTA phi and rho = c^5). Gas of test/common/nondim-Air: its molar mass gives
# R = 1 / 1.4 with a universal gas constant of 8314.51; with the exact SI value of FLINT and of mrcheck.RU,
# R = 0.7142820 and gamma = 1.3999968, which the case and test mr-moving-levels take from gas(). The cell values are the primitive variables
# of the exact cell averages of the conserved variables (Gauss-Legendre with 8 points), as MOSE writes them.
SIMPLE_DELTA = 0.05
SIMPLE_GL = np.polynomial.legendre.leggauss(8)
INI_SIMPLE = """[MOSE-Parameters]
iter-threshold = 100000000
time-threshold = 0.5
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.8
vnn = 0.3
cfl-rise-threshold = 0
time-scheme = RK2
time-accurate = true
irs-beta = 0.0
riemann-solver = HLLC
flux-limiter = vanleer
space-reconstruction = MUSCL
integration-variables = cons

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
res-diter = 1000000000
"""


def simple_state(x0, gam):
    """rho, u, p of the simple wave at the foot x0 of the characteristic (the initial state at x0)"""
    phi = 0.5 * (1.0 - np.cos(2.0 * np.pi * x0))
    c = 1.0 + SIMPLE_DELTA * phi
    u = 2.0 * (c - 1.0) / (gam - 1.0)
    rho = c ** (2.0 / (gam - 1.0))
    return rho, u, rho * c * c / gam


def simple_averages(xn, gam, t=0.0, x0_of=None):
    """cell averages of (rho, rho u, E) of the simple wave at the time t on the nodes xn (x0_of(x, t) inverts the
    characteristics; at t = 0 x0 = x)"""
    a, b = xn[:-1], xn[1:]
    xm, hw = 0.5 * (a + b), 0.5 * (b - a)
    X = xm[:, None] + hw[:, None] * SIMPLE_GL[0][None, :]
    rho, u, p = simple_state(X if x0_of is None else x0_of(X, t), gam)
    E = p / (gam - 1.0) + 0.5 * rho * u * u
    w = 0.5 * SIMPLE_GL[1][None, :]
    return (rho * w).sum(1), (rho * u * w).sum(1), (E * w).sum(1)


def simple_case(d, air, n):
    """mr1d-simple-<n>: the simple wave on [0; 1] in n uniform cells, one block periodic in x (201 on itself), slip y and z
    faces of the 1D slabs; gas of test/common/nondim-Air (next to the air given with --air); T = 0.5 before the breaking at
    t = 1.06. The level boundaries of the multirate with L = 1, sigma = 2 and fixed-macro-step 0.565 lie where
    lambda = u + c = 1 + 3 DELTA, phi = 1/2, and travel with the characteristics of x0 = 1/4 and 3/4"""
    nd = os.path.join(os.path.dirname(os.path.abspath(air)), 'nondim-Air')
    R, gam = gas(nd)
    xn = np.linspace(0.0, 1.0, n + 1)
    X, Y, Z = slab_1d(xn)
    m, mu, E = simple_averages(xn, gam)
    u = mu / m
    p = (gam - 1.0) * (E - 0.5 * m * u * u)
    sh = (1, 1, n)
    B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': m.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh), 'w': np.zeros(sh),
         'p': p.reshape(sh)}
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], ['rho1', 'u', 'v', 'w', 'p'])
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(bc_single_block(n, 1, True))
    for f in ('phase.txt', 'thermo.dat'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(nd, f), 'rb').read())
    open(os.path.join(d, 'input.ini'), 'w').write(INI_SIMPLE)


# mass fractions of the two sides of test/1D/diff1 (its build_ic.py, table 1 of the reference there)
MIX_OXIDIZER = {'CH4': 0.0, 'H2O': 0.1, 'N2': 0.758, 'O2': 0.142}
MIX_FUEL = {'CH4': 0.214, 'H2O': 0.0, 'N2': 0.591, 'O2': 0.195}


def mix_case(d, air):
    """mr1d-mix: see the list of the cases. In 400 macro steps of about 5 us the slug moves by about 0.2 m. The partial
    densities follow from p = rho RU T / W_mix with the molar masses of phase.txt, in its order."""
    src = os.path.join(os.path.dirname(os.path.abspath(air)), '..', '1D', 'diff1', 'INPUT')
    for f in ('phase.txt', 'thermo.dat', 'transport.dat', 'composition.txt'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(src, f), 'rb').read())
    species = [ln.split() for ln in open(os.path.join(src, 'phase.txt')).read().splitlines()[1:] if ln.strip()]
    xn = banded(BANDS_1D)
    X, Y, Z = slab_1d(xn)
    xc = 0.5 * (xn[1:] + xn[:-1])
    slug = (xc > 0.2) & (xc < 0.3)
    fractions = [np.where(slug, MIX_FUEL[s], MIX_OXIDIZER[s]) for s, _ in species]
    rho = P0 / (RU * T0) / sum(y / float(w) for y, (_, w) in zip(fractions, species))
    sh = (1, 1, xc.size)
    B = {'X': X, 'Y': Y, 'Z': Z, 'u': np.full(sh, 100.0), 'v': np.zeros(sh), 'w': np.zeros(sh), 'p': np.full(sh, P0)}
    names = []
    for k, y in enumerate(fractions, 1):
        names.append('rho%d' % k)
        B[names[-1]] = (y * rho).reshape(sh)
    write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], names + ['u', 'v', 'w', 'p'])
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(bc_single_block(xc.size, 1, True))
    ini = INI_VISC.format(iters=100, cfl=0.8, vnn=0.5).replace('schmidt = 0.7\n', 'schmidt = 0.7\nchemistry = frozen\n')
    open(os.path.join(d, 'input.ini'), 'w').write(ini)


def make(case, air, out):
    d = os.path.join(out, case)
    if os.path.exists(d):
        print('exists', d)
        return
    os.makedirs(os.path.join(d, 'INPUT'))
    if case in PLATES:
        plate_case(d, air, *PLATES[case])
        print('created', d)
        return
    if case == 'mr2d-plate-2br':
        plate_case(d, air, 1, 2)
        fi, fb = os.path.join(d, 'INPUT', 'ic.tec'), os.path.join(d, 'INPUT', 'bc.txt')
        ic, bc = turn_block_text(open(fi).read(), open(fb).read(), 2)
        open(fi, 'w').write(ic)
        open(fb, 'w').write(bc)
        print('created', d)
        return
    if case == 'mr2d-plate-rev':
        plate_case(d, air, 1, 1)
        f = os.path.join(d, 'INPUT', 'bc.txt')
        rev = reverse_records(open(f).read())
        open(f, 'w').write(rev)
        print('created', d)
        return
    if case in TURB_PLATES:
        turb_plate_case(d, air, *TURB_PLATES[case])
        print('created', d)
        return
    if case in CHANNELS:
        chan_case(d, air, *CHANNELS[case])
        print('created', d)
        return
    if case in ('mr3d-hot', 'mr3d-hot-8b'):
        hot3d_case(d, air, cut=case == 'mr3d-hot-8b')
        print('created', d)
        return
    if case == 'mr1d-nozzle-su':
        nozzle_case(d, air)
        print('created', d)
        return
    if case.startswith('mr1d-simple-') and case[len('mr1d-simple-'):].isdigit():
        simple_case(d, air, int(case[len('mr1d-simple-'):]))
        print('created', d)
        return
    if case == 'mr2d-wedge-visc':
        wedge_case(d, air)
        print('created', d)
        return
    if case == 'mr2d-wedge-band':
        wedge_band_case(d, air)
        print('created', d)
        return
    if case == 'mr2d-plate-302':
        plate_case(d, air, 1, 1)
        f = os.path.join(d, 'INPUT', 'bc.txt')
        bc = wall_302(open(f).read())
        open(f, 'w').write(bc)
        print('created', d)
        return
    if case in DARU:
        daru_case(d, air, *DARU[case])
        print('created', d)
        return
    if case == 'mr1d-mix':
        mix_case(d, air)
        print('created', d)
        return
    R, gam = gas(air)
    rho0 = P0 / (R * T0)
    c0 = math.sqrt(gam * R * T0)
    names = ['rho1', 'u', 'v', 'w', 'p']
    blocks = None                       # set by the cases of several blocks (split_blocks)
    if case == 'mr1d-deep':
        B = deep_block(4, 8, gam, rho0, c0)
        bc = bc_single_block(B['p'].shape[2], 1, True)
    elif case in ('mr1d-acu', 'mr1d-ent', 'mr1d-sod'):
        xn = banded(BANDS_1D)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        if case == 'mr1d-acu':
            eps, u0 = 1.0e-3, 50.0
            s = np.sin(2 * np.pi * xc / LDOM)
            p = P0 * (1 + eps * s); rho = rho0 * (1 + eps * s / gam); u = u0 + eps * c0 / gam * s
            per = True
        elif case == 'mr1d-ent':
            s = np.sin(2 * np.pi * xc / LDOM)
            rho = rho0 * (1 + 0.1 * s); u = np.full_like(xc, 100.0); p = np.full_like(xc, P0)
            per = True
        elif case == 'mr1d-sod':
            left = xc < 0.25
            rho = np.where(left, 1.2, 0.15); p = np.where(left, 1.0e5, 1.0e4); u = np.zeros_like(xc)
            per = False
        else:
            raise SystemExit('unknown case ' + case)
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        bc = bc_single_block(xc.size, 1, per)
    elif case in ('mr1d-nozzle', 'mr1d-skew'):
        xn = banded(BANDS_1D)
        X, Y, Z = slab_1d(xn)
        if case == 'mr1d-nozzle':
            Y[:, 1, :] = 0.01 * (1.0 + 0.5 * xn / LDOM)
            Z[1] = -0.5 * DZ + DZ * (1.0 + 0.5 * xn / LDOM)
        else:
            X[:, 1, :] += SKEW_J
            X[1] += SKEW_K
        xc = 0.5 * (xn[1:] + xn[:-1])
        left = xc < 0.25
        rho = np.where(left, 1.2, 0.15); p = np.where(left, 1.0e5, 1.0e4); u = np.zeros_like(xc)
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        bc = bc_single_block(xc.size, 1, False)
    elif case == 'mr2d-free':
        blocks, bc = split_blocks(free_block(rho0), [], [], True, True)
    elif case == 'mr2d-box':
        B = box_block(0.35, 0.40, gam, rho0)
        bc = bc_single_block(B['p'].shape[2], B['p'].shape[1], False)
    elif case == 'mr3d-sector':
        ni, nj, nk = 8, 4, 2
        xs = np.linspace(0.0, 0.08, ni + 1)
        rs = np.linspace(0.05, 0.09, nj + 1)
        ts = np.radians(np.linspace(-2.5, 2.5, nk + 1))
        X = np.empty((nk + 1, nj + 1, ni + 1)); Y = np.empty_like(X); Z = np.empty_like(X)
        for k in range(nk + 1):
            for j in range(nj + 1):
                X[k, j] = xs
                Y[k, j] = rs[j] * np.cos(ts[k])
                Z[k, j] = rs[j] * np.sin(ts[k])
        sh = (nk, nj, ni)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.full(sh, rho0), 'u': np.zeros(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': np.full(sh, P0)}
        out = []
        for f in (1, 2):
            for k in range(1, nk + 1):
                for j in range(1, nj + 1):
                    out.append(rec(1, 1 if f == 1 else ni, j, k, f, 300))
        for f in (3, 4):
            for k in range(1, nk + 1):
                for i in range(1, ni + 1):
                    out.append(rec(1, i, 1 if f == 3 else nj, k, f, 300))
        for f in (5, 6):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    out.append(rec(1, i, j, 1 if f == 5 else nk, f, 200))
        bc = ''.join(out)
    elif case == 'mr1d-deep1':
        B = deep_block(1, 2, gam, rho0, c0)
        bc = bc_single_block(B['p'].shape[2], 1, True)
    elif case in ('mr1d-acu-3b', 'mr1d-acu-3bi'):
        blocks, bc = split_blocks(acu_block(gam, rho0, c0), [15, 50 if case == 'mr1d-acu-3b' else 140], [], True, False)
    elif case == 'mr1d-uni-3b':
        blocks, bc = split_blocks(acu_block(gam, 1.2, c0, eps=0.0), [15, 50], [], True, False)
    elif case == 'mr1d-acu-sa':
        # mr1d-acu with the variable of Spalart-Allmaras (rho nu~ = 3 mu, free stream): only for the rejection of the
        # RANS runs (test mr-reject)
        B = acu_block(gam, rho0, c0)
        B['mi_t'] = np.full(B['p'].shape, 3.0 * 1.8e-5)
        names = names + ['mi_t']
        bc = bc_single_block(B['p'].shape[2], 1, True)
    elif case == 'mr1d-deep12':
        B = deep_block(4, 8, gam, rho0, c0, L=12)
        bc = bc_single_block(B['p'].shape[2], 1, True)
    elif case == 'mr1d-deep12-1':
        B = deep_block(1, 2, gam, rho0, c0, L=12)
        bc = bc_single_block(B['p'].shape[2], 1, True)
    elif case == 'mr1d-deep12-3b':
        blocks, bc = split_blocks(deep_block(4, 8, gam, rho0, c0, L=12), [105, 141], [], True, False)
    elif case in ('mr2d-box-4b', 'mr2d-box-4bx', 'mr2d-box-sym', 'mr2d-box-sym-4b'):
        xp, yp = (0.35, 0.40) if case in ('mr2d-box-4b', 'mr2d-box-4bx') else (0.4992, 0.4992)
        B = box_block(xp, yp, gam, rho0)
        if case == 'mr2d-box-sym':
            bc = bc_single_block(B['p'].shape[2], B['p'].shape[1], False)
        elif case == 'mr2d-box-4bx':
            blocks, bc = split_blocks(B, [18], [18], False, False)
        else:
            blocks, bc = split_blocks(B, [33], [33], False, False)
    elif case == 'mr2d-free-4b':
        blocks, bc = split_blocks(free_block(rho0), [24], [24], True, True)
    elif case == 'mr2d-vortex-3b':
        blocks, bc = split_blocks(vortex_block(R, gam, rho0), [30, 150], [], True, True)
    elif case == 'mr2d-hot':
        blocks, bc = split_blocks(hot_block(R, gam), [], [], True, True)
    elif case == 'mr2d-hot-4b':
        blocks, bc = split_blocks(hot_block(R, gam), [30], [30], True, True)
    elif case == 'mr2d-wall-s':
        blocks, bc = split_blocks(wall_s_block(R, gam), [], [], True, False)
        # face 3 (y = 0): adiabatic no-slip wall 301 (property line of the wall: qw and two unused values, as daru_case)
        out = []
        for ln in bc.splitlines(True):
            t = ln.split()
            if len(t) == 6 and t[4] == '3' and t[5] == '300':
                out.append(ln[:ln.rindex('300')] + '301\n' + '   0.0, 0.0, 0.0\n')
            else:
                out.append(ln)
        bc = ''.join(out)
    elif case == 'mr2d-shock':
        X, Y, Z, XC, YC = box_mesh()
        left = XC < 0.2
        sh = (1,) + XC.shape
        z = np.zeros(sh)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': np.where(left, 1.2, 0.15).reshape(sh), 'u': z, 'v': z, 'w': z,
             'p': np.where(left, 1.0e5, 1.0e4).reshape(sh)}
        bc = bc_single_block(XC.shape[1], XC.shape[0], False)
    elif case in ('mr1d-401', 'mr1d-402', 'mr1d-403', 'mr1d-405', 'mr1d-408', 'mr1d-410', 'mr1d-420', 'mr1d-502'):
        B, bc = duct_case(d, int(case[-3:]), R, gam)
    elif case == 'mr1d-502c':
        B, bc = duct_case(d, 502, R, gam, closed=True)
    elif case in ('mr1d-blast', 'mr1d-blast100', 'mr1d-blast100-3b'):
        xn = np.linspace(0.0, 1.0, 401)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        hot = xc < 0.1
        p = np.where(hot, (5.0 if case == 'mr1d-blast' else 100.0) * P0, P0)
        T = np.where(hot, 6400.0, 100.0)
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': (p / (R * T)).reshape(sh), 'u': np.zeros(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        if case == 'mr1d-blast100-3b':
            blocks, bc = split_blocks(B, [133, 266], [], False, False)
        else:
            bc = bc_single_block(xc.size, 1, False)
    elif case in ('mr1d-blast-wide', 'mr1d-blast-wide-3b'):
        # the nodes summed one cell at a time, from x = 0
        h0 = 2.5e-3
        xn = [0.0]
        for _ in range(200):
            xn.append(xn[-1] + h0)
        for q in range(1, 9):
            xn.append(xn[-1] + h0 * 2 ** q)
        for _ in range(2):
            xn.append(xn[-1] + h0 * 2 ** 8)
        xn = np.array(xn)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        hot = xc < 0.1
        p = np.where(hot, 5.0 * P0, P0)
        T = np.where(hot, 6400.0, 100.0)
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': (p / (R * T)).reshape(sh), 'u': np.zeros(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        if case == 'mr1d-blast-wide-3b':
            blocks, bc = split_blocks(B, [100, 200], [], False, False)
        else:
            bc = bc_single_block(xc.size, 1, False)
    elif case in ('mr1d-blast-hot', 'mr1d-blast-hot-3b'):
        xn = np.linspace(0.0, 1.0, 401)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        hot = xc < 0.1
        p = np.where(hot, 5.0 * P0, P0)
        T = np.where(hot, 9990.0, 100.0)
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': (p / (R * T)).reshape(sh), 'u': np.zeros(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        if case == 'mr1d-blast-hot-3b':
            blocks, bc = split_blocks(B, [20, 200], [], False, False)
        else:
            bc = bc_single_block(xc.size, 1, False)
    elif case == 'mr1d-sod-hot':
        xn = banded(BANDS_1D)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        left = xc < 0.25
        rho = np.where(left, 1.2, 0.15); p = np.where(left, 1.0e5, 1.0e4); u = np.zeros_like(xc)
        rho[9] = p[9] / (R * 2.0e4)                 # cell 10: 2e4 K, above the tables
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        bc = bc_single_block(xc.size, 1, False)
    elif case in ('mr1d-sod-pneg', 'mr1d-sod-p0', 'mr1d-sod-pneg-3b', 'mr1d-sod-nneg', 'mr1d-sod-vac'):
        xn = banded(BANDS_1D)
        X, Y, Z = slab_1d(xn)
        xc = 0.5 * (xn[1:] + xn[:-1])
        left = xc < 0.25
        rho = np.where(left, 1.2, 0.15); p = np.where(left, 1.0e5, 1.0e4); u = np.zeros_like(xc)
        if case == 'mr1d-sod-nneg':
            rho[9], p[9] = -1.0e-6, -0.0861                # cell 10: density and pressure negative, about 300 K
        elif case == 'mr1d-sod-vac':
            rho[99], p[99] = 1.0e-307, 8.6115e-303         # cell 100, in the band of 0.64 mm: about 300 K
        else:
            p[9] = 0.0 if case == 'mr1d-sod-p0' else -p[9]   # cell 10: zero pressure, or that of its neighbours negated
        sh = (1, 1, xc.size)
        B = {'X': X, 'Y': Y, 'Z': Z, 'rho1': rho.reshape(sh), 'u': u.reshape(sh), 'v': np.zeros(sh),
             'w': np.zeros(sh), 'p': p.reshape(sh)}
        if case == 'mr1d-sod-pneg-3b':
            blocks, bc = split_blocks(B, [5, 20], [], False, False)
        else:
            bc = bc_single_block(xc.size, 1, False)
    else:
        raise SystemExit('unknown case ' + case)
    if blocks is not None:
        write_tec(os.path.join(d, 'INPUT', 'ic.tec'), blocks, names)
    else:
        write_tec(os.path.join(d, 'INPUT', 'ic.tec'), [B], names)
    open(os.path.join(d, 'INPUT', 'bc.txt'), 'w').write(bc)
    for f in ('phase.txt', 'thermo.dat'):
        open(os.path.join(d, 'INPUT', f), 'wb').write(open(os.path.join(air, f), 'rb').read())
    # constant transport properties (mu = 1.8e-5 Pa s, Pr = 0.72), so that the viscous set-up of MOSE can run (the
    # Euler runs do not read them); same table format as test/2D/viscous/mms/INPUT/transport.dat
    cp_air = float(np.loadtxt(os.path.join(air, 'thermo.dat'), skiprows=4)[0, 1])
    viscous = case in ('mr2d-hot', 'mr2d-hot-4b', 'mr2d-wall-s')
    mu, pr = (HOT_MU if viscous else 1.8e-5), 0.72
    with open(os.path.join(d, 'INPUT', 'transport.dat'), 'w') as f:
        f.write('TITLE = "Transport Properties"\nVARIABLES = "Temperature", "Viscosity", "Conductivity"\n'
                'ZONE T="Air"\nI=10000, F=POINT\n')
        f.write(''.join('%.1f %.10e %.10e\n' % (t, mu, mu * cp_air / pr) for t in range(1, 10001)))
    cfl = {'mr2d-box': 0.2, 'mr2d-free': 0.1, 'mr3d-sector': 0.4, 'mr2d-box-4b': 0.2, 'mr2d-box-4bx': 0.2,
           'mr2d-box-sym': 0.2, 'mr2d-box-sym-4b': 0.2, 'mr2d-free-4b': 0.1, 'mr2d-vortex-3b': 0.2, 'mr2d-shock': 0.2,
           'mr2d-hot': 0.2, 'mr2d-hot-4b': 0.2, 'mr2d-wall-s': 0.2}.get(case, 0.8)
    if viscous:
        open(os.path.join(d, 'input.ini'), 'w').write(INI_VISC.format(iters=100, cfl=cfl, vnn=0.5))
    else:
        open(os.path.join(d, 'input.ini'), 'w').write(INI.format(iters=100, cfl=cfl))
    print('created', d)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--air', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('cases', nargs='*', default=['mr1d-acu', 'mr1d-ent', 'mr1d-sod', 'mr2d-box'])
    a = ap.parse_args()
    for c in a.cases:
        make(c, a.air, a.out)
