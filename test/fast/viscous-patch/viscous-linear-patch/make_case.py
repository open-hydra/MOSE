"""Inputs and checks of the ViscousLinearPatch and ViscousWallPatch fast tests (patch test of the viscous fluxes).

A linear velocity field has a uniform viscous stress and a linear temperature field a uniform heat flux (mu and k constant):
the viscous force and the heat conduction of every cell are zero. A face gradient that is exact for linear fields gives
zero in every cell, to round-off, on any mesh: skewed, stretched, with jumps of cell size, with folds at the faces between
blocks, at points where 3 or 5 blocks meet, next to walls whose temperature varies along the wall. The difference
approximations of the field and of the coordinates must be the same (M. Vinokur, NASA CR-177416, 1986, p. 26).
This is a test of consistency, the patch test of finite elements, not a test of the order of accuracy.

Meshes (unit coordinates, scaled by L = 1 cm and moved to x >= 1, y >= 2; 2-D meshes are one cell thick in z, with
symmetry (300) k faces and the same in-plane move for both node layers; the interior nodes are moved at random, with
fixed seeds, by up to 0.25 of the cell width on fold2, 0.15 of the cell size on box3d and sym, 0.15/n on star3 and star5
(n cells along a block side, polygon of unit circumradius), 0.2 of the cell size on rand1, not at all on jump, paraT and
paraQ):
  fold2  two blocks of 8 x 8 cells joined at x = 0.5 (block faces 2 and 1); the j lines of block 2 rise as y = t + 0.25 s t,
         so the grid lines fold at the face between the blocks;
  jump   one block of 8 x 8 cells whose sizes alternate in the ratio 1.5 in x and in y, sheared by x += 0.5 y, nodes not
         moved: the face gradients of mean are not exact there, but for a linear field their errors
         cancel in the sums of the interior cells; on this mesh the check bites in the cells next to the walls;
  jump2  as jump with the sizes 1 and 2 (units of 1/16: the block is 0.75 wide), the first cells next to the faces 1
         and 3 half as wide as the second ones, unit 2^-7 m instead of 1 cm: every coordinate is a dyadic fraction, so
         that the ghost nodes extrapolated from the first three node planes are exact with any compiler (as
         3 N0 - 3 N1 + N2 they would fall on the face, a ghost cell of zero width);
  fan    one block of 8 x 8 cells whose j lines spread away from the face 1: the j spacing grows 7/4 times across the
         first cell (from 1/16 at x = 0 to 7/64 at x = 1/8, then uniform), unit 2^-7 m (dyadic coordinates, as jump2):
         the ghost cells extrapolated behind the face 1 have edges along it of 1/4 and -1/2 of those of the face, the
         first layer narrower than the first cell, the second turned over (no flux reads its metric);
  fan4   as fan with the j spacing growing 4 times across the first cell (1/16 at x = 0, 1/4 from x = 1/8 on): the
         ghost cells behind the face 1 have edges along it of -2 and -5 times those of the face (-1/2 and -5/4 of
         those of the node plane after it), both layers turned over (det -1/5 and -7/5 of the first cell); only the
         viscous flux of a chimera face (102) would read the metric of the first, and the face 1 is a 400 face (mom)
         or a wall 302 (ene);
  star3  a triangle cut into 3 blocks of 8 x 8 cells that meet at its centre (a point shared by 3 blocks);
  star5  a pentagon cut into 5 blocks of 8 x 8 cells that meet at its centre (a point shared by 5 blocks);
  box3d  two blocks of 6 x 6 x 6 cells joined at x = 0.5, the j lines of block 2 rising as y = t + 0.2 s t;
  sym    one block of 8 x 8 cells sheared by x += 0.35 y on a symmetry plane (300) at y = 2;
  paraT  one block of 8 x 8 parallelogram cells sheared by 20 degrees, walls 302 at y = 2 and 3, faces 400 at the sides,
  paraQ  the same with walls 301: controls, run ene only, where every face metric is exact (every cell has the same
         metric, T_w uniform along the 302 walls of paraT, q = 0 on the 301 walls of paraQ): mean passes
         them too, so the test is not stricter than the defects it is written for;
  rand1  one block of 8 x 8 cells sheared by x += 0.3 y (replica only).

Runs:
  mom    uniform p = 10 Pa and T, so uniform rho; u = U0 + G (x - x0), G with shear and dilatation; outer faces 400 (no
         viscous flux there). The low pressure keeps the round-off of the pressure forces, p (face area)/volume, far
         below the viscous scale. Checked: the momentum of the cells without a 400 face.
  ene    u = 0, uniform rho, T = T0 + g . (x - x0), p = rho R T (so the extrapolated ghost cells stay on the linear field);
         outer faces of type 302 (T_w = T at the face centre, the mean of its four nodes), 301 (q into the fluid: k g . n
         with the face normal of the solver on faces 2, 4, 6, -k g . n on faces 1, 3, 5) or 400, mixed so that the faces
         between blocks end on every kind of ghost cell.
         Checked: the energy of the cells without a 400 face.
  momw   uniform p = 10 Pa and T, u = (a, b, c) (y - 2 cm), zero on the wall y = 2 cm (302, T_w = T) of rand1 and box3d,
         whose wall cells have nodes moved at random; the other outer faces 400. Checked: the momentum of the wall cells.
  momr   as momw in a frame turning at OMEGA_R about an axis through the origin (z on the 2-D meshes, (1, 2, 2)/3 on
         box3d), with the wall y = 2 cm a stationary face of the frame (stationary-face): there the wall moves with
         -Omega x x, which varies along it, and u = -Omega x x + (a, b, c) (y - 2 cm) takes that velocity on the
         wall. The sources of the frame are the same in the reference step. On rand1 and box3d with the default face
         metric; on parallelogram cells, where the metric of the wall cell is exact too, with mean on paraT and on
         paraQ, whose wall is an adiabatic 301. Checked: the momentum of the wall cells; on paraT and paraQ,
         at p = 1 Pa (P_WORK), also the energy of the cells, against the work of the uniform stress, tau : grad u (see
         work below).
  rk3    u = 0, uniform p = 10 Pa, T linear (rho = p/(R T)): steady for the Euler equations and, with a face gradient
         exact for linear fields, for the Navier-Stokes equations, at every stage (the low pressure: see P_LOW).
         Three-stage Runge-Kutta step, 12 x 12 cells per block (star3, star5 only). Checked are the cells within 4 cells
         of the point shared by the blocks, whose stencils at every stage hold the ghost cells and the cells of the
         other blocks at that point.
  enep   as ene, but with a uniform pressure p = 10 Pa, so rho = p/(R T) is not uniform: the ghost cells of the walls and
         400 faces must extrapolate T itself (and the 302 walls take T_w), not the density and the pressure, for T to
         stay on the linear field. Checked: the energy of the cells without a 400 face.
  hot    paraT only, uniform T, p = 10 Pa + 4e4 Pa/m (x - 1 cm), walls 302 at T_w = T/2: the first cell is hotter than
         1.5 T_w and the pressure more than doubles between the first two cells at the 400 faces, so the floor of the
         wall temperature and the cut of the extrapolation act. Checked: the line of the run that counts them.
  ene103, ene104, mom103, mom104
         rand1 and box3d, the wall y = 2 cm a wall coupled to another solver (103, or 104 with its donor lines), the
         other outer faces 400. MOSE gives such a wall no viscous flux (the other solver adds it), and fills its ghost
         cells as those of an adiabatic wall (Ghost_Wall). ene: as ene, with T linear along the wall, so that no heat
         crosses it: the energy of the wall cells is zero. mom: as momw: the momentum of the wall cells is minus the
         force of the uniform stress on the wall face, -tau . A/V (A the outward area vector of the face). Both read the
         ghost cells of the wall through the faces between the wall cells, at the positions of the extrapolated (T) and
         of the wall-valued (u) quantities.
  momio1, momio2, momio3, eneio1, eneio2, eneio3
         every outer face an inlet or an outlet, one type per block face (OPEN_SETS): rand1, set 1, faces 1 to 4 401,
         402, 403, 404; set 2 405, 406, 407, 420; set 3 408 on the four faces; box3d, set 1, faces 1, 3, 4, 5, 6 of
         block 1 401, 402, 403, 404, 405 and faces 2, 3, 4, 5, 6 of block 2 406, 407, 408, 408, 420. MOSE gives these
         faces an inviscid flux only, and fills their ghost cells with the cell (401-407, 420: Ghost_ZG_Extrapolate) or
         with 3 P1 - 3 P2 + P3 of the primitive variables, rho_s, u, p (408: Ghost_Extrapolate). mom: as mom; ene:
         u = 0, a uniform density, that of p = 10 Pa at T0, and T linear, so that p = rho R T is linear and the
         ghost cells of 408 hold T on the linear field (with a density that varies they would not: they extrapolate
         rho_s and p, not T). The cells with faces of these types must give minus the flux of the uniform stress (mom)
         or of the uniform heat conduction (ene) through those faces, -tau . A/V or -k g . A/V summed over them, and
         the other cells zero; they read the ghost cells through the faces between the cells of a face (and, on rand1,
         through their symmetry faces in k), at the positions that the face metric gives them. The property lines hold
         the state of the cell of the face (open_line), so that the inviscid flux of the face, the same in the two
         steps, stays at the scale of the fluxes of the cells; the low pressure keeps its round-off far below the
         viscous scale, as in mom.
  mom410, enep410, sa410, turb410
         rand1 and box3d, the plane y = 2 cm a face of type 410 (assigned state from a file), the other outer faces 400.
         The file INPUT/q2d-410.dat, written here as ATLAS writes it (one zone per block with the nodes of the face and
         one value per cell of the face: rho(1), u, v, w, p), holds the state of the field at the centre of the ghost
         cell of each cell of the face, the cell mirrored through the face, 2 x_f - x_1 (x_f the mean of the four nodes
         of the face, x_1 the centre of the cell): the cell across the face of the solver. MOSE gives that face no
         viscous flux; the faces between the cells next to it and the source terms of Spalart-Allmaras read its ghost
         cells, the state of the file at 2 x_f - x_1 and, for the turbulence variable, which the file does not give,
         the value of the cell per unit mass at x_1. mom410: as mom, the momentum of the cells of the face is minus the
         force of the uniform stress on the face, -tau . A/V (A the outward area vector of the face); enep410: as enep,
         their energy is -k g . A/V. sa410: Spalart-Allmaras, the fluid at rest at uniform p and T (rho uniform), nu~
         linear: no wall, so the wall distance is 1e8 m and production and destruction stay below 1e-20 of the scale;
         the right-hand side of rho nu~ is the flux (mu + rho nu~_f)/sigma (grad nu~ . A) of every face but the 410 one,
         with rho nu~_f the mean of the two cells (the state of the face), plus rho cb2/sigma |grad nu~|^2 (sa_replica
         below). turb410: Spalart-Allmaras, nu~ uniform, a uniform velocity along the grid lines across the face (the
         lines along which j varies), so that no gas crosses the faces 1, 2, 5, 6, into the domain through the 410 face,
         whose file has 1.2 times the density of the cells: the gradient of nu~ is zero and every flux of rho nu~ is nu~
         times the mass flux (the MUSCL slopes of rho nu~ are nu~ times those of rho, and every limiter is homogeneous
         of degree one), so r(rho nu~) = nu~ r(rho) in every cell (one species). Checked: the cells without a 400 face.

Observable: one explicit step from the same state with the Navier-Stokes equations and one reference step (mom, mom410:
the Euler equations; ene, momw, momr, momio, eneio, enep410: Navier-Stokes with mu = k = 1e-30 and q = 0 on the 301
faces, because the solver refuses walls 301 and 302 with the Euler equations and the ghost cells must stay those of the
walls; rk3, sa410: the Navier-Stokes step only, whose inviscid part is zero for this state; turb410: the Navier-Stokes
step only, whose inviscid part the check takes as nu~ r(rho)). The diagnostic output holds, per cell, the residual over
the step, r, and the step dt; r/dt of Navier-Stokes minus r/dt of the reference is the viscous right-hand side per unit
volume. It is compared with the scale mu |G| (sum of face areas)/volume (mom, momw, momr, momio) or k |g| (sum of face
areas)/volume (ene, rk3, eneio); threshold 1e-11 (for the energy of momr next to the 302 wall of paraT, plus a bound of
the round-off of the temperature of the cell: see work). The cells are grouped by the faces they have: interior (no
boundary face), connection (a face between blocks; corner: two of them, the cells at a point shared by 3 or 5 blocks),
symmetry (300), wall 301, wall 302; in the runs momio and eneio a cell with faces of inlets and outlets is in the group
of each of their types.

ViscousLinearPatch checks every group but wall 302 and coupled, plus another metric: on the rand1 mesh, run mom,
diffusive-metric = mean must give, cell by cell, the viscous force of an independent replica of its face gradients
(written here from Lib_Diffusive.f90 and Lib_Metrics.f90 at bff5a51, with the ghost values of the 400 faces of
Ghost_Wall), within 1e-10 of the scale. ViscousWallPatch checks the cells with a 302 face: walls whose temperature varies
along the wall, on skewed cells (ene), and the velocity at a no-slip wall whose cells have nodes moved at random (momw),
also at a wall that does not turn with a rotating frame, whose velocity in the frame varies along it, with every face
metric (momr; on paraQ the cells with a 301 face); and the cells next to walls 301 and 302 with T linear at uniform
pressure (enep), the count of the limits of the wall ghost cells (hot), and the cells next to a face of type 410 (mom410,
enep410, sa410, turb410). ViscousOpenPatch checks every group of the runs momio and eneio.

Usage (from a test directory):
  python3 -B make_case.py case MESH RUN [METRIC]     INPUT/ic.tec; INPUT/bc-ns.txt, input-ns.ini for the Navier-Stokes
                                                      step; INPUT/bc-ref.txt, input-ref.ini, INPUT/transport-ref.dat
                                                      for the reference step; METRIC: mean (no key
                                                      when absent: the default of the solver)
  python3 -B make_case.py check MESH RUN NS REF GROUPS prints the largest error of each group and exits with 0 when the
                                                      groups listed (comma-separated) are below the threshold
  python3 -B make_case.py replica MESH NS REF METRIC  compares the viscous force with the replica of METRIC (mom run)
  python3 -B make_case.py work MESH RUN NS REF GROUPS  run momr on paraT or paraQ: the energy of the cells of the groups
                                                      against tau : grad u, exits with 0 when below the threshold
  python3 -B make_case.py limits LOG                  the run hot: the counts of the floor and of the cut are positive
  python3 -B make_case.py open MESH RUN NS REF        runs momio and eneio: prints the largest error of each group and
                                                      exits with 0 when every group is below the threshold
"""
import math
import os
import re
import sys

import numpy as np

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous: constant viscosity, Pa s
KCOND = 0.016544513989              # test/common/Air-viscous: constant conductivity, W/(m K)
P0, T0 = 1.0e5, 300.0
RHO0 = P0 / (RGAS * T0)
P_LOW = 10.0                        # mom, rk3: a low pressure. The round-off of the pressure forces, p (face area)/volume,
                                    # stays far below the viscous scale mu |G| (face area)/volume (mom: p/(mu |G|) ~ 4e3);
                                    # the velocity it gives at the first stage carries an energy flux (rho e + p) u, of
                                    # order p, into the later stages (rk3: 3e-10 of the scale at 1e5 Pa, 1e-13 at 10 Pa)
P_WORK = 1.0                        # momr on paraT and paraQ, whose energy is checked (work): the round-off of the energy
                                    # fluxes, (rho e + p) u (face area)/volume, is 1e-13 of the scale of the work at 1 Pa,
                                    # 1e-12 at 10 Pa (gfortran; 1e-12 and 4e-12 with ifx)
LREF = 1.0e-2                       # m
LREF2 = 2.0**-7                     # m, jump2, fan and fan4: dyadic coordinates
OFFSET = np.array([1.0, 2.0, 0.0])  # unit coordinates: keep y > 0 away from the axis
TOL = 1.0e-11
TOL_REPLICA = 1.0e-10
N_ULP_T = 6                         # momr on paraT: the temperature of a cell next to the 302 wall can differ from T_w by up
                                    # to 6 units in the last place of T_w (see work)
FMT = '{: .17E}'
MESHES = ('fold2', 'jump', 'jump2', 'fan', 'fan4', 'star3', 'star5', 'box3d', 'sym', 'paraT', 'paraQ', 'rand1')
RUNS = ('mom', 'ene', 'rk3', 'momw', 'momr', 'enep', 'hot', 'ene103', 'ene104', 'mom103', 'mom104',
        'momio1', 'momio2', 'momio3', 'eneio1', 'eneio2', 'eneio3',
        'mom410', 'enep410', 'sa410', 'turb410')
BC_CONNECT, BC_SYM, BC_Q, BC_T, BC_OPEN = 101, 300, 301, 302, 400
BC_COUPLED = (103, 104)                    # walls coupled to another solver (104: with the donor cells of a chimera)
BC_INOUT = (401, 402, 403, 404, 405, 406, 407, 408, 420)   # inlets and outlets: an inviscid flux only
OPEN_SETS = {                              # runs momioN, eneioN: (mesh, N) -> {(block, face): type}
    ('rand1', 1): {(0, 1): 401, (0, 2): 402, (0, 3): 403, (0, 4): 404},
    ('rand1', 2): {(0, 1): 405, (0, 2): 406, (0, 3): 407, (0, 4): 420},
    ('rand1', 3): {(0, 1): 408, (0, 2): 408, (0, 3): 408, (0, 4): 408},
    ('box3d', 1): {(0, 1): 401, (0, 3): 402, (0, 4): 403, (0, 5): 404, (0, 6): 405,
                   (1, 2): 406, (1, 3): 407, (1, 4): 408, (1, 5): 408, (1, 6): 420},
}
P0_TABLE = 'INPUT/p0-402.dat'              # the table of p0 of the 402 faces (time, p0 in bar)
GAMMA_AIR = 1.4                            # the 405 faces: the Mach number of the speed of the cell, |u|/sqrt(gamma R T)
BC_MAPPED = 410                            # assigned state from a file
RUNS_410 = ('mom410', 'enep410', 'sa410', 'turb410')
RUNS_SA = ('sa410', 'turb410')             # Spalart-Allmaras: the turbulence variable rho nu~ after p
Q2D_FILE = 'INPUT/q2d-410.dat'             # the file of the 410 face
SIGMA, CB2 = 2.0 / 3.0, 0.622              # Spalart-Allmaras (Lib_Spalart.f90)
NU0 = 0.1                                  # sa410, turb410: nu~ at x0, m^2/s (chi = rho nu~/mu about 1 at P_LOW)
GNU2 = np.array([2.0, -1.5, 0.0])          # sa410: gradient of nu~, m/s, on rand1 ...
GNU3 = np.array([2.0, -1.5, 1.0])          # ... and on box3d
RHO_RATIO_410 = 1.2                        # turb410: density of the file over that of the cells
U_410 = 1.0                                # turb410: speed of the flow, m/s


# ----------------------------------------------------------------------------------------------------------------------
# Meshes: list of blocks, each an array of nodes X[i, j, k, :] (unit coordinates), and the type of each outer face
# ----------------------------------------------------------------------------------------------------------------------
def perturb(X, scale, seed, fixed):
    """Move the nodes not in `fixed` (a boolean mask) by up to `scale` in each in-plane direction (and in z in 3-D)."""
    rng = np.random.RandomState(seed)
    d = rng.uniform(-1.0, 1.0, X.shape) * scale
    if X.shape[2] == 2:                  # one cell thick: the same in-plane move for the two node layers
        d[..., 2] = 0.0
        d[:, :, 1, :] = d[:, :, 0, :]
    d[fixed] = 0.0
    return X + d


def boundary_mask(X):
    m = np.zeros(X.shape[:3], bool)
    m[0, :, :] = m[-1, :, :] = True
    m[:, 0, :] = m[:, -1, :] = True
    if X.shape[2] > 2:
        m[:, :, 0] = m[:, :, -1] = True
    return m


def extrude(XY, dz):
    """2-D nodes (ni+1, nj+1, 2) -> one cell thick in z."""
    ni1, nj1 = XY.shape[:2]
    X = np.zeros((ni1, nj1, 2, 3))
    X[:, :, :, 0] = XY[:, :, None, 0]
    X[:, :, :, 1] = XY[:, :, None, 1]
    X[:, :, 1, 2] = dz
    return X


def mesh_fold2(n=8):
    s = np.linspace(0.0, 1.0, n + 1)
    S, T = np.meshgrid(s, s, indexing='ij')
    b1 = np.stack([0.5 * S, T], axis=-1)
    b2 = np.stack([0.5 + 0.5 * S, T + 0.25 * S * T], axis=-1)
    h = 0.5 / n
    blocks = [extrude(b1, h), extrude(b2, h)]
    blocks = [perturb(X, 0.25 * h, 11 + b, boundary_mask(X)) for b, X in enumerate(blocks)]
    return blocks


def alternating(n, r):
    h = 1.0 / n
    sizes = np.array([2.0 / (1.0 + r) if k % 2 == 0 else 2.0 * r / (1.0 + r) for k in range(n)]) * h
    x = np.concatenate([[0.0], np.cumsum(sizes)])
    x[-1] = 1.0
    return x


def mesh_jump(n=8):
    a = alternating(n, 1.5)
    S, T = np.meshgrid(a, a, indexing='ij')
    return [extrude(np.stack([S + 0.5 * T, T], axis=-1), 1.0 / n)]


def mesh_jump2(n=8):
    a = np.concatenate([[0.0], np.cumsum([1.0 if k % 2 == 0 else 2.0 for k in range(n)])]) / 16.0
    S, T = np.meshgrid(a, a, indexing='ij')
    return [extrude(np.stack([S + 0.5 * T, T], axis=-1), 1.0 / n)]


def mesh_fan(n=8, top=7.0 / 64.0):
    """The j spacing 1/16 on the face x = 0 and top (7/64, or 1/4 for fan4) from the next node plane on: a fan across the
    first cell."""
    XY = np.zeros((n + 1, n + 1, 2))
    for i in range(n + 1):
        XY[i, :, 0] = i / n
        XY[i, :, 1] = np.arange(n + 1) * (1.0 / 16.0 if i == 0 else top)
    return [extrude(XY, 1.0 / n)]


def mesh_star(m, n):
    """A regular m-gon cut into m quadrilateral blocks that meet at its centre C. Block b: (xi, eta) = (0, 0) at C, (1, 0)
    at the middle of the side before vertex b, (1, 1) at vertex b, (0, 1) at the middle of the side after it, so that the
    eta = 0 side of block b+1 is the xi = 0 side of block b (the same nodes)."""
    ang = [2.0 * math.pi * b / m + 0.3 for b in range(m)]
    V = [np.array([math.cos(a), math.sin(a)]) for a in ang]
    E = [0.5 * (V[b] + V[(b + 1) % m]) for b in range(m)]
    C = np.zeros(2)
    t = np.linspace(0.0, 1.0, n + 1)
    spoke = [np.array([C + tt * (E[b] - C) for tt in t]) for b in range(m)]     # C -> E_b, shared by blocks b and b+1
    blocks = []
    for b in range(m):
        Eb, Em = E[b], E[(b - 1) % m]
        XY = np.zeros((n + 1, n + 1, 2))
        for i in range(n + 1):
            for j in range(n + 1):
                x, y = t[i], t[j]
                XY[i, j] = (1 - x) * (1 - y) * C + x * (1 - y) * Em + x * y * V[b] + (1 - x) * y * Eb
        XY[:, 0] = spoke[(b - 1) % m]          # eta = 0: C -> E_{b-1}
        XY[0, :] = spoke[b]                    # xi = 0: C -> E_b
        blocks.append(XY)
    h = 0.5 / n
    out = []
    for b, XY in enumerate(blocks):
        X = extrude(XY, h)
        out.append(perturb(X, 0.15 * h * 2, 31 + b, boundary_mask(X)))
    return out


def mesh_box3d(n=6):
    s = np.linspace(0.0, 1.0, n + 1)
    S, T, U = np.meshgrid(s, s, s, indexing='ij')
    b1 = np.stack([0.5 * S, T, U], axis=-1)
    b2 = np.stack([0.5 + 0.5 * S, T + 0.2 * S * T, U], axis=-1)
    h = 1.0 / n
    return [perturb(X, 0.15 * h, 51 + b, boundary_mask(X)) for b, X in enumerate([b1, b2])]


def mesh_rand1(n=8):
    s = np.linspace(0.0, 1.0, n + 1)
    S, T = np.meshgrid(s, s, indexing='ij')
    X = extrude(np.stack([S + 0.3 * T, T], axis=-1), 1.0 / n)
    return [perturb(X, 0.2 / n, 91, boundary_mask(X))]


def mesh_para(n=8):
    """Parallelogram cells, sheared by 20 degrees, nodes not moved: every cell has the same metric."""
    s = np.linspace(0.0, 1.0, n + 1)
    S, T = np.meshgrid(s, s, indexing='ij')
    return [extrude(np.stack([S + math.tan(math.radians(20.0)) * T, T], axis=-1), 1.0 / n)]


def mesh_sym(n=8):
    s = np.linspace(0.0, 1.0, n + 1)
    S, T = np.meshgrid(s, s, indexing='ij')
    X = extrude(np.stack([S + 0.35 * T, T], axis=-1), 1.0 / n)
    fixed = boundary_mask(X)
    return [perturb(X, 0.15 / n, 71, fixed)]


def build(mesh, run):
    n2 = 12 if run == 'rk3' else 8
    if mesh == 'fold2':
        blocks = mesh_fold2()
    elif mesh == 'jump':
        blocks = mesh_jump()
    elif mesh == 'jump2':
        return [LREF2 * (X + OFFSET) for X in mesh_jump2()]
    elif mesh == 'fan':
        return [LREF2 * (X + OFFSET) for X in mesh_fan()]
    elif mesh == 'fan4':
        return [LREF2 * (X + OFFSET) for X in mesh_fan(top=0.25)]
    elif mesh == 'star3':
        blocks = mesh_star(3, n2)
    elif mesh == 'star5':
        blocks = mesh_star(5, n2)
    elif mesh == 'box3d':
        blocks = mesh_box3d()
    elif mesh == 'sym':
        blocks = mesh_sym()
    elif mesh in ('paraT', 'paraQ'):
        blocks = mesh_para()
    elif mesh == 'rand1':
        blocks = mesh_rand1()
    else:
        sys.exit('unknown mesh ' + mesh)
    return [LREF * (X + OFFSET) for X in blocks]


def coupled_type(run):
    """103 or 104 for the runs on a coupled wall (ene103, mom104, ...), else None."""
    return int(run[-3:]) if run[-3:] in ('103', '104') else None


def open_run(run):
    """('mom' or 'ene', N) for the runs on inlets and outlets, momioN and eneioN, else None."""
    m = re.match(r'^(mom|ene)io([1-9])$', run)
    return (m.group(1), int(m.group(2))) if m else None


def outer_type(mesh, run, b, face):
    """Type of the outer face `face` of block b (faces 5, 6 of 2-D meshes are symmetry planes)."""
    if open_run(run):                    # inlets and outlets, one type per block face
        return OPEN_SETS[(mesh, open_run(run)[1])][(b, face)]
    if coupled_type(run):                # one wall, the plane y = 2 cm, coupled to another solver
        return coupled_type(run) if face == 3 else BC_OPEN
    if run in RUNS_410:                  # one face of type 410, the plane y = 2 cm
        return BC_MAPPED if face == 3 else BC_OPEN
    if run == 'mom':
        if mesh == 'sym' and face == 3:
            return BC_SYM
        return BC_OPEN
    if run in ('momw', 'momr'):          # one wall, the plane y = 2 cm: 302, or 301 (q = 0) on paraQ
        if face != 3:
            return BC_OPEN
        return BC_Q if mesh == 'paraQ' else BC_T
    if run == 'rk3':
        mix = {0: (BC_T, BC_T), 1: (BC_Q, BC_T), 2: (BC_OPEN, BC_Q), 3: (BC_T, BC_OPEN), 4: (BC_Q, BC_Q)}
        return mix[b][0 if face == 2 else 1]
    # ene
    if mesh == 'fold2':
        return {(0, 1): BC_T, (0, 3): BC_Q, (0, 4): BC_T, (1, 2): BC_OPEN, (1, 3): BC_T, (1, 4): BC_Q}[(b, face)]
    if mesh in ('jump', 'jump2', 'fan', 'fan4', 'rand1'):
        return {1: BC_T, 2: BC_Q, 3: BC_T, 4: BC_OPEN}[face]
    if mesh == 'paraT':
        return {1: BC_OPEN, 2: BC_OPEN, 3: BC_T, 4: BC_T}[face]
    if mesh == 'paraQ':
        return {1: BC_OPEN, 2: BC_OPEN, 3: BC_Q, 4: BC_Q}[face]
    if mesh in ('star3', 'star5'):
        mix = {0: (BC_T, BC_T), 1: (BC_Q, BC_T), 2: (BC_OPEN, BC_Q), 3: (BC_T, BC_OPEN), 4: (BC_Q, BC_Q)}
        return mix[b][0 if face == 2 else 1]
    if mesh == 'box3d':
        return {(0, 1): BC_T, (0, 3): BC_Q, (0, 4): BC_T, (0, 5): BC_T, (0, 6): BC_Q,
                (1, 2): BC_OPEN, (1, 3): BC_T, (1, 4): BC_Q, (1, 5): BC_Q, (1, 6): BC_T}[(b, face)]
    if mesh == 'sym':
        return {1: BC_T, 2: BC_OPEN, 3: BC_SYM, 4: BC_T}[face]
    sys.exit('no outer types for ' + mesh)


# ----------------------------------------------------------------------------------------------------------------------
# Geometry of the solver: cell centre = mean of the 8 nodes; faces from the nodes as in Compute_Norm_Area
# ----------------------------------------------------------------------------------------------------------------------
def centres(X):
    return 0.125 * (X[:-1, :-1, :-1] + X[1:, :-1, :-1] + X[:-1, 1:, :-1] + X[1:, 1:, :-1] +
                    X[:-1, :-1, 1:] + X[1:, :-1, 1:] + X[:-1, 1:, 1:] + X[1:, 1:, 1:])


def face_nodes(X, face, i, j, k):
    """The four nodes of face `face` of cell (i, j, k) (1-based cell indices)."""
    if face in (1, 2):
        a = i - 1 if face == 1 else i
        return [X[a, j - 1, k - 1], X[a, j, k - 1], X[a, j, k], X[a, j - 1, k]]
    if face in (3, 4):
        a = j - 1 if face == 3 else j
        return [X[i - 1, a, k - 1], X[i, a, k - 1], X[i, a, k], X[i - 1, a, k]]
    a = k - 1 if face == 5 else k
    return [X[i - 1, j - 1, a], X[i, j - 1, a], X[i, j, a], X[i - 1, j, a]]


def face_signs(X):
    """signi, signj, signk of Compute_Norm_Area: orientation of the normals towards growing index."""
    out = []
    d1 = X[1, 1, 1] - X[1, 0, 0]; d2 = X[1, 0, 1] - X[1, 1, 0]
    c = 0.25 * (X[1, 1, 1] + X[1, 0, 1] + X[1, 0, 0] + X[1, 1, 0]) - 0.25 * (X[0, 1, 1] + X[0, 0, 1] + X[0, 0, 0] + X[0, 1, 0])
    out.append(math.copysign(1.0, np.dot(np.cross(d1, d2), c)))
    d1 = X[0, 1, 1] - X[1, 1, 0]; d2 = X[1, 1, 1] - X[0, 1, 0]
    c = 0.25 * (X[1, 1, 1] + X[0, 1, 1] + X[0, 1, 0] + X[1, 1, 0]) - 0.25 * (X[1, 0, 1] + X[0, 0, 1] + X[0, 0, 0] + X[1, 0, 0])
    out.append(math.copysign(1.0, np.dot(np.cross(d1, d2), c)))
    d1 = X[0, 1, 1] - X[1, 0, 1]; d2 = X[0, 0, 1] - X[1, 1, 1]
    c = 0.25 * (X[1, 1, 1] + X[0, 1, 1] + X[0, 0, 1] + X[1, 0, 1]) - 0.25 * (X[1, 1, 0] + X[0, 1, 0] + X[0, 0, 0] + X[1, 0, 0])
    out.append(math.copysign(1.0, np.dot(np.cross(d1, d2), c)))
    return out


def face_normal_area(X, sgn, d, i, j, k):
    """Unit normal (towards growing index d) and area of the face of index (i, j, k) in direction d (node indices of the
    face as in Compute_Norm_Area: an i face (i, j, k) lies between cells (i, j, k) and (i+1, j, k), 1-based)."""
    if d == 1:
        d1 = X[i, j, k] - X[i, j - 1, k - 1]; d2 = X[i, j - 1, k] - X[i, j, k - 1]
    elif d == 2:
        d1 = X[i - 1, j, k] - X[i, j, k - 1]; d2 = X[i, j, k] - X[i - 1, j, k - 1]
    else:
        d1 = X[i - 1, j, k] - X[i, j - 1, k]; d2 = X[i - 1, j - 1, k] - X[i, j, k]
    v = 0.5 * np.cross(d1, d2)
    a = math.sqrt(float(np.dot(v, v)))
    return v / a * sgn[d - 1], a


def cell_faces(X, sgn, i, j, k):
    """Outward area vectors of the six faces of cell (i, j, k), faces 1..6."""
    out = []
    for face in range(1, 7):
        d = (face + 1) // 2
        fi, fj, fk = i, j, k
        if face == 1:
            fi = i - 1
        elif face == 3:
            fj = j - 1
        elif face == 5:
            fk = k - 1
        n, a = face_normal_area(X, sgn, d, fi, fj, fk)
        out.append((n * a) if face % 2 == 0 else (-n * a))
    return out


def cell_volume(X, i, j, k):
    """Compute_Metric_Tensor: the volume of the cell, the integral of det J of the trilinear map of its nodes N1..N8,
    1/6 sum_d e_d . (S_d- + S_d+), e_d the mean edge of direction d and S_d-, S_d+ the area vectors (diagonal x diagonal)/2
    of the two faces of direction d, along +d; the same formula as the solver (docs/theory/numerics.md, Cell volume), in
    absolute value as the solver stores it."""
    N = [X[i - 1, j - 1, k - 1], X[i - 1, j - 1, k], X[i - 1, j, k - 1], X[i - 1, j, k],
         X[i, j - 1, k - 1], X[i, j - 1, k], X[i, j, k - 1], X[i, j, k]]
    r = [p - N[0] for p in N]                                       # r[m] = N(m+1) - N1
    e1 = 0.25 * (r[4] + (r[5] - r[1]) + (r[6] - r[2]) + (r[7] - r[3]))
    e2 = 0.25 * (r[2] + (r[3] - r[1]) + (r[6] - r[4]) + (r[7] - r[5]))
    e3 = 0.25 * (r[1] + (r[3] - r[2]) + (r[5] - r[4]) + (r[7] - r[6]))
    v = (np.dot(e1, np.cross(r[3], r[1] - r[2]) + np.cross(r[7] - r[4], r[5] - r[6]))
         + np.dot(e2, np.cross(r[5], r[4] - r[1]) + np.cross(r[7] - r[2], r[6] - r[3]))
         + np.dot(e3, np.cross(r[6], r[2] - r[4]) + np.cross(r[7] - r[1], r[3] - r[5])))
    return abs(float(v)) / 12.0


# ----------------------------------------------------------------------------------------------------------------------
# Connectivity: the faces two blocks share become pairs of connection records (101)
# ----------------------------------------------------------------------------------------------------------------------
def dims(X):
    return X.shape[0] - 1, X.shape[1] - 1, X.shape[2] - 1


def face_cells(X, face):
    ni, nj, nk = dims(X)
    if face in (1, 2):
        return [((1 if face == 1 else ni), j, k) for k in range(1, nk + 1) for j in range(1, nj + 1)]
    if face in (3, 4):
        return [(i, (1 if face == 3 else nj), k) for k in range(1, nk + 1) for i in range(1, ni + 1)]
    return [(i, j, (1 if face == 5 else nk)) for j in range(1, nj + 1) for i in range(1, ni + 1)]


def key_of(nodes):
    return tuple(sorted(tuple(np.round(p / LREF, 9)) for p in nodes))


def connections(blocks):
    """{(b, face, i, j, k): (bs, fs, is, js, ks)} for every face shared by two cells of different blocks."""
    seen = {}
    for b, X in enumerate(blocks):
        twod = X.shape[2] == 2
        for face in range(1, 7):
            if twod and face >= 5:
                continue
            for (i, j, k) in face_cells(X, face):
                seen.setdefault(key_of(face_nodes(X, face, i, j, k)), []).append((b, face, i, j, k))
    conn = {}
    for v in seen.values():
        if len(v) == 2:
            conn[v[0]] = v[1]
            conn[v[1]] = v[0]
        elif len(v) > 2:
            sys.exit('face shared by more than two cells')
    return conn


def tangential(face):
    """Index directions (0-based) tangential to `face`, first then second, as the solver orders them."""
    return {1: (1, 2), 2: (1, 2), 3: (0, 2), 4: (0, 2), 5: (0, 1), 6: (0, 1)}[face]


def orientation(blocks, conn, rec):
    """d11, d12, d21, d22 of a connection record: the step of the source cell, in the two tangential directions of the
    source face, for a step of +1 of the master cell along its first and second tangential direction."""
    b, face, i, j, k = rec
    bs, fs, is_, js, ks = conn[rec]
    out = []
    for t in tangential(face):
        idx = [i, j, k]
        lim = dims(blocks[b])[t]
        if lim == 1:                     # one cell: the k direction of a 2-D mesh maps onto itself
            out.append((0, 1) if fs in (1, 2, 3, 4) else (1, 0))
            continue
        step = 1 if idx[t] < lim else -1
        idx[t] += step
        nb = conn.get((b, face, idx[0], idx[1], idx[2]))
        if nb is None or nb[0] != bs or nb[1] != fs:
            sys.exit('the face between blocks does not continue along a tangential direction')
        ds = [nb[2] - is_, nb[3] - js, nb[4] - ks]
        ts = tangential(fs)
        out.append((step * ds[ts[0]], step * ds[ts[1]]))
    return out[0][0], out[0][1], out[1][0], out[1][1]


# ----------------------------------------------------------------------------------------------------------------------
# Fields
# ----------------------------------------------------------------------------------------------------------------------
X0 = LREF * (OFFSET + np.array([0.5, 0.5, 0.25]))
G2 = np.array([[120.0, -80.0, 0.0], [60.0, 150.0, 0.0], [0.0, 0.0, 0.0]])
G3 = np.array([[120.0, -80.0, 50.0], [60.0, 150.0, -70.0], [-40.0, 90.0, 110.0]])
GSYM = np.array([[120.0, 0.0, 0.0], [0.0, 150.0, 0.0], [0.0, 0.0, 0.0]])
U0 = np.array([1.0, 0.5, 0.0])
g2 = np.array([2000.0, -1500.0, 0.0])
g3 = np.array([2000.0, -1500.0, 1000.0])
gsym = np.array([2000.0, 0.0, 0.0])
gparaT = np.array([0.0, 2000.0, 0.0])      # normal to the walls y = 2, 3 cm of paraT: T_w uniform along each wall
gparaQ = np.array([2000.0, 0.0, 0.0])      # along the walls of paraQ: q = k g . n = 0 there
gwall2 = np.array([2000.0, 0.0, 0.0])      # ene103, ene104: along the wall y = 2 cm (2-D) ...
gwall3 = np.array([2000.0, 0.0, 1000.0])   # ... and on box3d
GP_HOT = 4.0e4                             # run hot: pressure gradient along x, Pa/m (50 Pa per cell of paraT)
OMEGA_R = 100.0                            # run momr: angular speed of the frame, rad/s


def axis_r(mesh):
    """Run momr: the axis of the frame as input.ini gives it, (1, 2, 2) on box3d, z on the 2-D meshes (no velocity along
    z, which their symmetry planes k would reflect)."""
    return np.array([1.0, 2.0, 2.0]) if mesh == 'box3d' else np.array([0.0, 0.0, 1.0])


def omega_r(mesh):
    """Run momr: the angular velocity of the frame, OMEGA_R times the unit axis, with the operations of the solver."""
    a = axis_r(mesh)
    return OMEGA_R * (a / math.sqrt(float(a.dot(a))))


def field(mesh, run):
    if open_run(run):                    # momio: the field of mom; eneio: the gradient of T of ene
        run = open_run(run)[0]
    if run in ('momw', 'momr', 'mom103', 'mom104'):
        # zero on the wall y = 2 cm, with shear and dilatation: u = (a, b, c) (y - 2 cm)
        G = np.array([[0.0, 120.0, 0.0], [0.0, 60.0, 0.0], [0.0, 40.0 if mesh == 'box3d' else 0.0, 0.0]])
        x0 = X0.copy()
        x0[1] = LREF * OFFSET[1]
        if run == 'momr':
            # plus W x = -Omega x x, the velocity in the frame of a point at rest (the origin of the frame at 0)
            om = omega_r(mesh)
            W = np.array([[0.0, om[2], -om[1]], [-om[2], 0.0, om[0]], [om[1], -om[0], 0.0]])
            return {'G': G + W, 'u0': W.dot(x0), 'x0': x0}
        return {'G': G, 'u0': np.zeros(3), 'x0': x0}
    if run in ('mom', 'mom410'):
        G = GSYM if mesh == 'sym' else (G3 if mesh == 'box3d' else G2)
        u0 = {'box3d': np.array([1.0, 0.5, 0.3]), 'sym': np.array([1.0, 0.0, 0.0])}.get(mesh, U0)
        x0 = X0.copy()
        if mesh == 'sym':
            x0[1] = LREF * OFFSET[1]           # u_y odd about the symmetry plane y = 2 cm
        return {'G': G, 'u0': u0, 'x0': x0}
    if run in ('ene103', 'ene104'):      # along the wall y = 2 cm: no heat crosses it
        return {'g': gwall3 if mesh == 'box3d' else gwall2, 'x0': X0.copy()}
    if run == 'sa410':                   # nu~ linear
        return {'gn': GNU3 if mesh == 'box3d' else GNU2, 'x0': X0.copy()}
    if run == 'turb410':                 # along the grid lines across the face, the lines along which j varies (rand1:
                                         # sheared by x += 0.3 y), so that no gas crosses the faces 1, 2, 5, 6
        d = np.array([0.0, 1.0, 0.0]) if mesh == 'box3d' else np.array([0.3, 1.0, 0.0])
        return {'u': U_410 * d / math.sqrt(float(d.dot(d))), 'x0': X0.copy()}
    g = {'sym': gsym, 'box3d': g3, 'paraT': gparaT, 'paraQ': gparaQ}.get(mesh, g2)
    return {'g': g, 'x0': X0.copy()}


def state(mesh, run, x):
    f = field(mesh, run)
    if open_run(run) and open_run(run)[0] == 'ene':
        rho = P_LOW / (RGAS * T0)        # uniform: the 408 ghost cells hold T on the linear field
        return rho, np.zeros(3), rho * RGAS * (T0 + f['g'].dot(x - f['x0']))
    if open_run(run) or run in ('mom', 'momw', 'momr', 'mom103', 'mom104', 'mom410'):
        u = f['u0'] + f['G'].dot(x - f['x0'])
        p = P_WORK if (run == 'momr' and mesh in ('paraT', 'paraQ')) else P_LOW
        return p / (RGAS * T0), u, p
    if run == 'sa410':
        return P_LOW / (RGAS * T0), np.zeros(3), P_LOW
    if run == 'turb410':
        return P_LOW / (RGAS * T0), f['u'], P_LOW
    if run == 'hot':
        p = P_LOW + GP_HOT * (x[0] - LREF * OFFSET[0])
        return p / (RGAS * T0), np.zeros(3), p
    T = T0 + f['g'].dot(x - f['x0'])
    if run in ('ene', 'ene103', 'ene104'):
        return RHO0, np.zeros(3), RHO0 * RGAS * T
    return P_LOW / (RGAS * T), np.zeros(3), P_LOW


def rans(mesh, run, x, rho):
    """Spalart-Allmaras runs: rho nu~ at x, for the density rho (sa410: nu~ linear; turb410: nu~ uniform)."""
    f = field(mesh, run)
    if run == 'sa410':
        return rho * (NU0 + f['gn'].dot(x - f['x0']))
    return rho * NU0


def file_state(mesh, run, x):
    """The state that the file of the 410 face gives at x, the centre of the ghost cell: that of the field (turb410: with
    RHO_RATIO_410 times its density)."""
    rho, u, p = state(mesh, run, x)
    if run == 'turb410':
        rho = RHO_RATIO_410 * rho
    return rho, u, p


def write_q2d(blocks, mesh, run):
    """The file of the 410 face y = 2 cm (face 3 of every block) in the form that ATLAS writes (map_strip_to_block in
    src/BCB/file_unwrapped.f90): VARIABLES "x" "y" "z" and the names of the state, one zone per block, T="B<block>_T1" with
    the nodes of the face, I along i and J along k, and one value of each variable per cell of the face, i fastest, at one
    time. The value of a cell is the state of the field at the centre of its ghost cell, 2 x_f - x_1, x_f the mean of the
    four nodes of the face and x_1 that of the eight nodes of the cell, as Ghost_Position computes them."""
    names = ('rho(1)', 'u', 'v', 'w', 'p')
    with open(Q2D_FILE, 'w') as out:
        out.write(' VARIABLES ="x" "y" "z"' + ''.join(f' "{n}"' for n in names) + '\n')
        for b, X in enumerate(blocks):
            ni, _, nk = dims(X)
            out.write(f' ZONE T="B{b + 1}_T1", I={ni + 1}, J={nk + 1}, K=1, DATAPACKING=BLOCK, '
                      f'VARLOCATION=([1-3]=NODAL,[4-{3 + len(names)}]=CELLCENTERED), SOLUTIONTIME= 0.00000000E+00\n')
            for axis in range(3):
                out.write(''.join(' ' + FMT.format(X[i, 0, k, axis]) + '\n' for k in range(nk + 1) for i in range(ni + 1)))
            C = centres(X)
            cells = []
            for k in range(1, nk + 1):
                for i in range(1, ni + 1):
                    xf = 0.25 * sum(face_nodes(X, 3, i, 1, k))
                    cells.append(file_state(mesh, run, 2.0 * xf - C[i - 1, 0, k - 1]))
            cols = ([c[0] for c in cells], [c[1][0] for c in cells], [c[1][1] for c in cells],
                    [c[1][2] for c in cells], [c[2] for c in cells])
            for col in cols:
                out.write(''.join(' ' + FMT.format(v) + '\n' for v in col))


# ----------------------------------------------------------------------------------------------------------------------
# Case files
# ----------------------------------------------------------------------------------------------------------------------
def outward_area(X, sgn, face, i, j, k):
    """Unit normal toward growing index (as the solver's) and outward area vector of face `face` of cell (i, j, k)."""
    d = (face + 1) // 2
    fi, fj, fk = i - (face == 1), j - (face == 3), k - (face == 5)
    n, a = face_normal_area(X, sgn, d, fi, fj, fk)
    return n, (n * a if face % 2 == 0 else -n * a)


def open_line(mesh, run, X, sgn, face, i, j, k, t):
    """Property line of an inlet or outlet face (boundary-conditions.md, Inlet / Outlet) from the state of its cell:
    T0 or T the temperature of the cell, p0, p or p_back its pressure, the mass flux g and the velocity un its own across
    the face, rho u . n and u . n (n the normal of the solver, toward growing index), a 405 face the Mach number and the
    direction of the velocity of the cell (alpha, beta as angles: Supersonic_Inflow takes the velocity along them, so the
    face gets the velocity of the cell, where with the token normal the direction would come from the normal of the
    solver, toward growing index), the other inlets the token normal, rel_fac 1 and the mass fraction of the single
    species 1. Each face takes the branch of its type that the state of its cell selects: an inflow where the gas
    of the cell moves in, else an outflow or the slip wall of a face the gas does not enter (eneio: the gas is at rest;
    the 420 faces, with p0 = p, are slip walls in both runs). So the inviscid flux of the face stays at the scale of the
    flux of the cell. The ghost cells of these types do not read the line."""
    rho, u, p = state(mesh, run, centres(X)[i - 1, j - 1, k - 1])
    T = p / (rho * RGAS)
    n, _ = outward_area(X, sgn, face, i, j, k)
    un = float(u.dot(n))
    tail = '1.0, 1.0,'                   # rel_fac, the mass fraction of the single species
    if t == 401:
        return f'    {T:.17E}, {p:.17E}, normal, normal, {tail}'
    if t == 402:
        return f'    {T:.17E}, {P0_TABLE}, normal, normal, {tail}'
    if t in (403, 404):
        return f'    {T:.17E}, {rho * un:.17E}, normal, normal, {tail}'
    if t == 405:
        mach = float(np.linalg.norm(u)) / math.sqrt(GAMMA_AIR * RGAS * T)
        alpha, beta = math.atan2(u[1], u[0]), math.atan2(u[2], math.hypot(u[0], u[1]))
        return f'    {mach:.17E}, {T:.17E}, {p:.17E}, {alpha:.17E}, {beta:.17E}, {tail}'
    if t == 406:
        return f'    {p:.17E}, 1.0,'
    if t == 407:
        return f'    {T:.17E}, {p:.17E}, {p:.17E}, normal, normal, {tail}'
    if t == 408:
        return f'    {T:.17E}, {un:.17E}, normal, normal, {tail}'
    if t == 420:                         # T0, p0, p_sub, p_sup, g: p0 = p, the branch of a face the gas does not enter
        return f'    {T:.17E}, {p:.17E}, 0.0, 0.0, 0.0, {tail}'
    sys.exit(f'no property line for type {t}')


def write_p0_table(blocks, mesh, run, conn):
    """The table of p0 of the 402 faces (time, p0 in bar, constant): the pressure of the cells (mom: uniform), or just
    below the lowest pressure of the cells of the 402 faces (eneio: the gas at rest does not enter)."""
    ps = []
    for b, X in enumerate(blocks):
        C = centres(X)
        for face in range(1, 7):
            if X.shape[2] == 2 and face >= 5:
                continue
            for (i, j, k) in face_cells(X, face):
                if (b, face, i, j, k) not in conn and outer_type(mesh, run, b, face) == 402:
                    ps.append(state(mesh, run, C[i - 1, j - 1, k - 1])[2])
    if not ps:
        return
    p0 = min(ps) if open_run(run)[0] == 'mom' else min(ps) * (1.0 - 1.0e-6)
    with open(P0_TABLE, 'w') as out:
        out.write(f'0.0 {p0 / 1.0e5:.17E}\n1.0 {p0 / 1.0e5:.17E}\n')


def write_case(mesh, run, metric):
    blocks = build(mesh, run)
    conn = connections(blocks)
    os.makedirs('INPUT', exist_ok=True)
    for p in ('INPUT/ic.tec', 'INPUT/bc.txt', 'INPUT/bc-ns.txt', 'INPUT/bc-ref.txt', 'INPUT/transport-ref.dat',
              P0_TABLE, Q2D_FILE):
        if os.path.lexists(p):
            os.remove(p)
    for b, X in enumerate(blocks):
        ni, nj, nk = dims(X)
        vmin = min(cell_volume(X, i, j, k) for k in range(1, nk + 1) for j in range(1, nj + 1) for i in range(1, ni + 1))
        if vmin <= 0.0:
            sys.exit(f'block {b + 1}: a cell with volume {vmin}')
    sa = run in RUNS_SA
    with open('INPUT/ic.tec', 'w') as out:
        out.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"' + (' "rhoNuTilde"' if sa else '') + '\n')
        for b, X in enumerate(blocks):
            ni, nj, nk = dims(X)
            out.write(f'ZONE T="Block{b + 1}", I={ni + 1}, J={nj + 1}, K={nk + 1}, '
                      f'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-{9 if sa else 8}]=CELLCENTERED)\n')
            for axis in range(3):
                out.write(''.join(' ' + FMT.format(X[i, j, k, axis]) + '\n'
                                  for k in range(nk + 1) for j in range(nj + 1) for i in range(ni + 1)))
            C = centres(X)
            cells = [state(mesh, run, C[i, j, k]) for k in range(nk) for j in range(nj) for i in range(ni)]
            cols = ([c[0] for c in cells], [c[1][0] for c in cells], [c[1][1] for c in cells],
                    [c[1][2] for c in cells], [c[2] for c in cells])
            if sa:
                cols += ([rans(mesh, run, C[i, j, k], cells[i + ni * (j + nj * k)][0])
                          for k in range(nk) for j in range(nj) for i in range(ni)],)
            for col in cols:
                out.write(''.join(' ' + FMT.format(v) + '\n' for v in col))
    if open_run(run):
        write_p0_table(blocks, mesh, run, conn)
    if run in RUNS_410:
        write_q2d(blocks, mesh, run)
    f = field(mesh, run)
    for tag in ('ns', 'ref'):
        rows = []
        for b, X in enumerate(blocks):
            sgn = face_signs(X)
            twod = X.shape[2] == 2
            for face in range(1, 7):
                for (i, j, k) in face_cells(X, face):
                    rec = (b, face, i, j, k)
                    head = f'{b + 1:8d}{i:8d}{j:8d}{k:8d}{face:8d}'
                    if twod and face >= 5:
                        rows.append(head + f'{BC_SYM:8d}\n')
                    elif rec in conn:
                        bs, fs, is_, js, ks = conn[rec]
                        d11, d12, d21, d22 = orientation(blocks, conn, rec)
                        rows.append(head + f'{BC_CONNECT:8d}\n' +
                                    f'{bs + 1:8d}{is_:8d}{js:8d}{ks:8d}{fs:8d}{d11:8d}{d12:8d}{d21:8d}{d22:8d}\n')
                    else:
                        t = outer_type(mesh, run, b, face)
                        if t == BC_T:
                            xf = 0.25 * sum(face_nodes(X, face, i, j, k))
                            tw = T0 if run in ('momw', 'momr') else T0 + f['g'].dot(xf - f['x0'])
                            if run == 'hot':
                                tw = 0.5 * T0
                            rows.append(head + f'{t:8d}\n' + f'    {tw:.17E},    0.000000E+00,    0.000000E+00,\n')
                        elif t == BC_Q:
                            d = (face + 1) // 2
                            fi, fj, fk = i, j, k
                            if face == 1:
                                fi = i - 1
                            elif face == 3:
                                fj = j - 1
                            elif face == 5:
                                fk = k - 1
                            nrm, _ = face_normal_area(X, sgn, d, fi, fj, fk)
                            into = 1.0 if face % 2 == 0 else -1.0     # q_wall is positive into the fluid
                            q = into * KCOND * f['g'].dot(nrm) if ('g' in f and tag == 'ns') else 0.0
                            rows.append(head + f'{t:8d}\n' + f'    {q:.17E},    0.000000E+00,    0.000000E+00,\n')
                        elif t == 103:           # the cell of the other solver: here the wall cell (MOSE does not use it)
                            rows.append(head + f'{t:8d}\n' +
                                        f'{b + 1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{1:8d}{0:8d}{0:8d}{1:8d}\n')
                        elif t == 104:           # one donor cell in each of the two lists of a chimera: here the wall cell
                            rows.append(head + f'{t:8d}\n' + f'{1:8d}{1:8d}\n' +
                                        2 * f'{b + 1:8d}{i:8d}{j:8d}{k:8d}    1.0\n')
                        elif t in BC_INOUT:
                            rows.append(head + f'{t:8d}\n' + open_line(mesh, run, X, sgn, face, i, j, k, t) + '\n')
                        elif t == BC_MAPPED:     # the name of its file
                            rows.append(head + f'{t:8d}\n' + Q2D_FILE + '\n')
                        else:
                            rows.append(head + f'{t:8d}\n')
        with open(f'INPUT/bc-{tag}.txt', 'w') as out:
            out.write(''.join(rows))
    # Reference run: the Euler equations (mom), or Navier-Stokes with mu = k = 1e-30 and q = 0 on the 301 faces (ene:
    # the solver refuses walls 301 and 302 with the Euler equations, and their ghost cells must stay those of the walls)
    tiny = []
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'common', 'Air-viscous',
                                  'transport.dat')):
        p = line.split()
        if len(p) == 3 and re.match(r'^[0-9.Ee+-]+$', p[0]):
            tiny.append(f'{p[0]} 1.0E-30 1.0E-30\n')
        else:
            tiny.append(line)
    with open('INPUT/transport-ref.dat', 'w') as out:
        out.write(''.join(tiny))
    key = '' if metric is None else f'diffusive-metric = {metric}\n'
    scheme = 'RK3' if run == 'rk3' else 'euler'
    rot, frame = '', ''
    if sa:                               # Spalart-Allmaras, without the point-implicit factor of steady runs (it scales r)
        rot = 'turbulence = SA\n'
        frame = '\n[MOSE-Turbulence]\npoint-implicit = false\n'
    if run == 'momr':                    # the frame of momr; the wall y = 2 cm (face 3 of every block) does not turn
        ax = axis_r(mesh)
        rot = 'rotational-frame = rigid-body\n'
        frame = ('\n[MOSE-Rotating-Frame]\n'
                 f'omega = {OMEGA_R}\naxis = {ax[0]} {ax[1]} {ax[2]}\n'
                 'origin = 0.0 0.0 0.0\n' + ''.join(f'stationary-face-{b + 1} = {b + 1} 3\n' for b in range(len(blocks))))
    for tag, eq in (('ns', 'navier-stokes'), ('ref', 'euler' if run in ('mom', 'mom410') else 'navier-stokes')):
        with open(f'input-{tag}.ini', 'w') as out:
            out.write('[MOSE-Parameters]\niter-threshold = 1\nnewrun = true\nres-threshold = 1e-300\n\n'
                      f'[MOSE-Physics]\nequations = {eq}\nschmidt = 0.7\n' + rot + '\n'
                      '[MOSE-Numerics]\ncfl = 0.5\nvnn = 0.3\ntime-accurate = false\n'
                      f'time-scheme = {scheme}\nintegration-variables = cons\nspace-reconstruction = MUSCL\n'
                      'flux-limiter = vanleer\nriemann-solver = HLLC\n' + key + frame)
    print(f'{mesh}/{run}: {len(blocks)} blocks, {sum(np.prod(dims(X)) for X in blocks)} cells, {len(rows)} records')


# ----------------------------------------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------------------------------------
def read_zones(path):
    """Cell values of each zone of a MOSE ASCII Tecplot file (BLOCK packing): list of dicts name -> array (i fastest)."""
    with open(path) as fh:
        text = fh.read()
    head, _, _ = text.partition('ZONE')
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', head.split('=', 1)[1])]
    zones = []
    parts = re.split(r'(?m)^\s*ZONE', text)[1:]
    for part in parts:
        line, _, body = part.partition('\n')
        ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', line).group(1)) for k in ('I', 'J', 'K'))
        vals = np.array(body.split(), float)
        nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1)
        z, pos = {}, 0
        for q, name in enumerate(names):
            cnt = nn if q < 3 else nc
            z[name] = vals[pos:pos + cnt]
            pos += cnt
        z['dims'] = (ni - 1, nj - 1, max(nk - 1, 1))
        zones.append(z)
    return zones


def viscous_rhs(ns, eu, comps):
    """Per zone: array (ncell, len(comps)) of r/dt (Navier-Stokes) - r/dt (reference); eu None: Navier-Stokes alone."""
    out = []
    for q, zn in enumerate(ns):
        a = np.stack([zn[c] / zn['dt'] for c in comps], axis=1)
        if eu is not None:
            ze = eu[q]
            a = a - np.stack([ze[c] / ze['dt'] for c in comps], axis=1)
        out.append(a)
    return out


def groups(blocks, mesh, run, conn):
    """Group of every cell (list per block, i fastest), or None when not checked."""
    out = []
    for b, X in enumerate(blocks):
        ni, nj, nk = dims(X)
        twod = X.shape[2] == 2
        g = []
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    kinds = []
                    for face in range(1, 7):
                        if twod and face >= 5:
                            continue
                        on = ((face == 1 and i == 1) or (face == 2 and i == ni) or (face == 3 and j == 1) or
                              (face == 4 and j == nj) or (face == 5 and k == 1) or (face == 6 and k == nk))
                        if not on:
                            continue
                        kinds.append('conn' if (b, face, i, j, k) in conn else outer_type(mesh, run, b, face))
                    if BC_OPEN in kinds:
                        g.append(None)
                    elif any(t in kinds for t in BC_COUPLED):
                        g.append('coupled')
                    elif run == 'rk3':
                        g.append('centre' if (i <= 4 and j <= 4) else None)
                    elif BC_T in kinds:
                        g.append('wall302')
                    elif BC_Q in kinds:
                        g.append('wall301')
                    elif BC_MAPPED in kinds:
                        g.append('mapped410')
                    elif BC_SYM in kinds:
                        g.append('symmetry')
                    elif kinds.count('conn') >= 2:
                        g.append('corner')
                    elif 'conn' in kinds:
                        g.append('connection')
                    else:
                        g.append('interior')
        out.append(g)
    return out


def scales(blocks, mesh, run):
    """Per block and cell (i fastest): coef (sum of face areas)/volume, coef = mu |G| (momentum), the largest
    (mu + rho nu~)/sigma times |grad nu~| (sa410), nu~ rho |u| (turb410: nu~ times the mass flux per unit area), k |g|
    (energy)."""
    f = field(mesh, run)
    if run in ('mom', 'momw', 'momr', 'mom103', 'mom104', 'mom410'):
        coef = MU * float(np.linalg.norm(f['G']))
    elif run == 'sa410':
        rho = P_LOW / (RGAS * T0)
        rnmax = max(rans(mesh, run, c, rho) for X in blocks for c in centres(X).reshape(-1, 3))
        coef = (MU + rnmax) / SIGMA * float(np.linalg.norm(f['gn']))
    elif run == 'turb410':
        coef = NU0 * P_LOW / (RGAS * T0) * float(np.linalg.norm(f['u']))
    else:
        coef = KCOND * float(np.linalg.norm(f['g']))
    out = []
    for X in blocks:
        ni, nj, nk = dims(X)
        sgn = face_signs(X)
        s = []
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    area = sum(math.sqrt(float(np.dot(a, a))) for a in cell_faces(X, sgn, i, j, k))
                    s.append(coef * area / cell_volume(X, i, j, k))
        out.append(np.array(s))
    return out


def check(mesh, run, ns_path, eu_path, wanted):
    blocks = build(mesh, run)
    conn = connections(blocks)
    if run in ('mom', 'momw', 'momr', 'mom103', 'mom104', 'mom410'):
        comps = ['rhou', 'rhov', 'rhow']
    elif run == 'sa410':
        comps = ['rhoNuTilde']
    elif run == 'turb410':
        comps = ['rhoNuTilde', 'rho']
    else:
        comps = ['rhoe']
    rhs = viscous_rhs(read_zones(ns_path), None if run in ('rk3', 'sa410', 'turb410') else read_zones(eu_path), comps)
    if run == 'turb410':                 # r(rho nu~) - nu~ r(rho)
        rhs = [a[:, :1] - NU0 * a[:, 1:] for a in rhs]
    grp = groups(blocks, mesh, run, conn)
    sc = scales(blocks, mesh, run)
    expected, on = None, ()
    if run in ('mom103', 'mom104'):
        expected, on = coupled_wall_rhs(blocks, mesh, run), ('coupled',)
    elif run in ('mom410', 'enep410'):
        expected, on = coupled_wall_rhs(blocks, mesh, run), ('mapped410',)
    elif run == 'sa410':
        expected, on = sa_replica(blocks, mesh, run, conn), ('interior', 'connection', 'mapped410')
    worst = {}
    for b in range(len(blocks)):
        for c, gname in enumerate(grp[b]):
            if gname is None:
                continue
            want = expected[b][c] if (expected is not None and gname in on) else 0.0
            e = float(np.max(np.abs(rhs[b][c] - want))) / sc[b][c]
            worst[gname] = max(worst.get(gname, 0.0), e)
    ok = True
    for gname in sorted(worst):
        good = worst[gname] <= TOL
        if gname in wanted:
            ok = ok and good
        mark = ('ok  ' if good else 'FAIL') if gname in wanted else '    '
        print(f'    {mesh:6s} {run:4s} {gname:10s} largest |viscous rhs| / scale = {worst[gname]:.3e}  [{mark}]')
    missing = [g for g in wanted if g not in worst]
    if missing:
        print(f'    {mesh} {run}: no cells in the groups {missing}')
    return ok


def coupled_wall_rhs(blocks, mesh, run):
    """Runs mom103, mom104, mom410: per block and cell (i fastest), the momentum right-hand side per unit volume of a cell
    whose face 3 is the coupled wall or the 410 face, -tau . A/V: a linear velocity has a uniform stress tau, whose forces
    on the faces of a cell sum to zero, and MOSE gives that face no viscous flux (A, the outward area vector of the face).
    Run enep410: the energy right-hand side, -k g . A/V, from the uniform heat flux of a linear temperature."""
    f = field(mesh, run)
    if run == 'enep410':
        tau = KCOND * f['g'].reshape(1, 3)       # one row: the energy flux per unit area along A is k g . A/|A|
    else:
        G = f['G']
        tau = MU * (G + G.T) - 2.0 / 3.0 * MU * np.trace(G) * np.eye(3)
    out = []
    for X in blocks:
        ni, nj, nk = dims(X)
        sgn = face_signs(X)
        e = []
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    a3 = cell_faces(X, sgn, i, j, k)[2]
                    e.append(-tau.dot(a3) / cell_volume(X, i, j, k) if j == 1 else np.zeros(len(tau)))
        out.append(e)
    return out


def sa_replica(blocks, mesh, run, conn):
    """Run sa410: per block and cell (i fastest), the right-hand side of rho nu~ per unit volume, (1/V) times the sum over
    the faces with a viscous flux (every face but the 410 one; the k faces of a 2-D mesh, symmetry planes across which nu~
    does not vary, give zero) of (mu + rho nu~_f)/sigma (grad nu~ . A), A the outward area vector of the face and rho nu~_f
    the mean of the two cells (Compute_Diffusive_Flux takes the state of a face as the mean of its two cells; across a face
    between blocks the cell of the other block), plus the cb2 term of the source, rho cb2/sigma |grad nu~|^2
    (SA_Source_Blk: grad rho = 0; with no wall, production and destruction stay below 1e-20 of the scale). Cells with a
    400 face, which are not checked, get the sum of their other faces."""
    gn = field(mesh, run)['gn']
    rho = P_LOW / (RGAS * T0)
    cent = [centres(X) for X in blocks]
    out = []
    for b, X in enumerate(blocks):
        n = dims(X)
        sgn = face_signs(X)
        twod = X.shape[2] == 2
        e = []
        for k in range(1, n[2] + 1):
            for j in range(1, n[1] + 1):
                for i in range(1, n[0] + 1):
                    faces = cell_faces(X, sgn, i, j, k)
                    rc = rans(mesh, run, cent[b][i - 1, j - 1, k - 1], rho)
                    tot = 0.0
                    for face in range(1, 7):
                        if twod and face >= 5:
                            continue
                        d = (face + 1) // 2 - 1
                        idx = [i, j, k]
                        idx[d] += 1 if face % 2 == 0 else -1
                        if 1 <= idx[d] <= n[d]:
                            xn = cent[b][idx[0] - 1, idx[1] - 1, idx[2] - 1]
                        elif (b, face, i, j, k) in conn:
                            bs, _, is_, js, ks = conn[(b, face, i, j, k)]
                            xn = cent[bs][is_ - 1, js - 1, ks - 1]
                        else:                    # an outer face: 410 (no viscous flux) or 400
                            continue
                        rf = 0.5 * (rc + rans(mesh, run, xn, rho))
                        tot += (MU + rf) / SIGMA * float(gn.dot(faces[face - 1]))
                    e.append(np.array([tot / cell_volume(X, i, j, k) + rho * CB2 / SIGMA * float(gn.dot(gn))]))
        out.append(e)
    return out


def check_open(mesh, run, ns_path, eu_path):
    """Runs momio, eneio: the right-hand side per unit volume of every cell against minus the flux of the uniform stress
    (mom: tau . A summed over the faces of the cell that are inlets or outlets, which have no viscous flux) or of the
    uniform heat conduction (ene: k g . A), over the volume; zero in the cells without such faces. Groups: one per type of
    inlet or outlet (the cells with a face of that type), interior, connection, corner; exits with 0 when every group is
    below the threshold."""
    blocks = build(mesh, run)
    conn = connections(blocks)
    kind = open_run(run)[0]
    f = field(mesh, run)
    if kind == 'mom':
        G = f['G']
        flux = MU * (G + G.T) - 2.0 / 3.0 * MU * np.trace(G) * np.eye(3)       # tau: the force on A is tau . A
        comps = ['rhou', 'rhov', 'rhow']
    else:
        flux = (KCOND * f['g']).reshape(1, 3)                                 # the heat conducted in through A, k g . A
        comps = ['rhoe']
    rhs = viscous_rhs(read_zones(ns_path), read_zones(eu_path), comps)
    sc = scales(blocks, mesh, kind)
    worst = {}
    for b, X in enumerate(blocks):
        ni, nj, nk = dims(X)
        sgn = face_signs(X)
        twod = X.shape[2] == 2
        c = 0
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    faces = cell_faces(X, sgn, i, j, k)
                    want = np.zeros(len(comps))
                    types, nconn = set(), 0
                    for face in range(1, 7):
                        on = ((face == 1 and i == 1) or (face == 2 and i == ni) or (face == 3 and j == 1) or
                              (face == 4 and j == nj) or (face == 5 and k == 1) or (face == 6 and k == nk))
                        if not on or (twod and face >= 5):
                            continue
                        if (b, face, i, j, k) in conn:
                            nconn += 1
                            continue
                        t = outer_type(mesh, run, b, face)
                        types.add(t)
                        want -= flux.dot(faces[face - 1]) / cell_volume(X, i, j, k)
                    e = float(np.max(np.abs(rhs[b][c] - want))) / sc[b][c]
                    names = ([str(t) for t in types] if types else
                             ['interior' if nconn == 0 else ('connection' if nconn == 1 else 'corner')])
                    for gname in names:
                        worst[gname] = max(worst.get(gname, 0.0), e)
                    c += 1
    ok = True
    for gname in sorted(worst, key=lambda g: (not g.isdigit(), g)):
        good = worst[gname] <= TOL
        ok = ok and good
        print(f'    {mesh:6s} {run:6s} {gname:10s} largest |viscous rhs - expected| / scale = {worst[gname]:.3e}'
              f'  [{"ok  " if good else "FAIL"}]')
    missing = [str(t) for t in set(OPEN_SETS[(mesh, open_run(run)[1])].values()) if str(t) not in worst]
    if missing:
        print(f'    {mesh} {run}: no cells next to the types {missing}  [FAIL]')
        ok = False
    return ok


def work(mesh, run, ns_path, eu_path, wanted):
    """Run momr on parallelogram cells (paraT, paraQ): the viscous right-hand side of the energy, per unit volume, is
    tau : grad u in every cell. A linear velocity has a uniform stress tau, so the net work of the stress on a cell, the
    sum over its faces of (tau . n) . u A, is (tau : grad u) V when u is the velocity at the centre of each face and the
    faces are plane (divergence theorem); on parallelogram cells the velocity of a face between cells, the mean of the
    velocities of the two cells, is the one at its centre, and the face on a wall that moves takes the wall velocity.
    Compared with the scale mu |G| U (sum of face areas)/volume, U the largest speed of the cell centres; threshold
    1e-11.

    A cell next to the 302 wall of paraT also takes the heat flux of the wall, k (T - T_w)/d over its area A_w, with d
    the distance of the cell centre from the wall face: on these parallelogram cells that is the gradient the solver takes
    across the wall face. T = T_w = T0 here, but the solver computes the temperature of the cell as p/(rho R), with R from
    the mass fractions (f_Rtot), and T_w is the T0 of bc.txt, read exactly. From the T0 of this script, T holds 8
    roundings of relative size at most u = 2^-53: the density written here, p/(RGAS T0), 2; RGAS, 1; R of FLINT, Runiv/W,
    1; f_Rtot, R rho/rho, 2; p/(rho R) of the wall flux, 2. Fast math (the RELEASE builds) may turn each of the three
    divisions of the solver into a multiplication by a rounded reciprocal: 11 roundings. So |T - T_w| <= 11 u T0, 6.4
    units in the last place (ulp) of T0 = 300 K, and T - T_w, the difference of two numbers between 256 and 512, is a
    whole number of ulp: at most N_ULP_T = 6. Each ulp adds k ulp(T0) A_w/(d V) to the energy of the cell, 1.9e-11 of
    the scale on paraT, where k T0 is 6e5 times mu |G| U d. The threshold of these cells is 1e-11 plus N_ULP_T times
    that term over the scale."""
    blocks = build(mesh, run)
    conn = connections(blocks)
    f = field(mesh, run)
    G = f['G']
    tau = MU * (G + G.T) - 2.0 / 3.0 * MU * np.trace(G) * np.eye(3)
    phi = float(np.sum(tau * G))
    rhs = viscous_rhs(read_zones(ns_path), read_zones(eu_path), ['rhoe'])
    grp = groups(blocks, mesh, run, conn)
    speed = max(float(np.linalg.norm(f['u0'] + G.dot(c - f['x0']))) for X in blocks for c in centres(X).reshape(-1, 3))
    worst, limit, within = {}, {}, {}
    for b, X in enumerate(blocks):
        ni, nj, nk = dims(X)
        sgn = face_signs(X)
        C = centres(X)
        c = 0
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    gname = grp[b][c]
                    if gname is not None:
                        faces = cell_faces(X, sgn, i, j, k)
                        area = sum(math.sqrt(float(np.dot(a, a))) for a in faces)
                        vol = cell_volume(X, i, j, k)
                        scale = MU * float(np.linalg.norm(G)) * speed * area / vol
                        tol = TOL
                        if gname == 'wall302':       # face 3, the wall y = 2 cm: the round-off of T through its heat flux
                            aw = math.sqrt(float(np.dot(faces[2], faces[2])))
                            xw = 0.25 * sum(face_nodes(X, 3, i, j, k))
                            d = abs(float(np.dot(C[i - 1, j - 1, k - 1] - xw, faces[2]))) / aw
                            tol += N_ULP_T * float(np.spacing(T0)) * KCOND * aw / (d * vol) / scale
                        e = abs(float(rhs[b][c][0]) - phi) / scale
                        worst[gname] = max(worst.get(gname, 0.0), e)
                        limit[gname] = min(limit.get(gname, tol), tol)
                        within[gname] = within.get(gname, True) and e <= tol
                    c += 1
    ok = True
    for gname in sorted(worst):
        good = within[gname]
        if gname in wanted:
            ok = ok and good
        mark = ('ok  ' if good else 'FAIL') if gname in wanted else '    '
        bound = '' if limit[gname] == TOL else f' (threshold {limit[gname]:.3e}: 1e-11 and {N_ULP_T} ulp of T_w)'
        print(f'    {mesh:6s} {run:4s} {gname:10s} largest |energy rhs - tau:grad u| / scale = {worst[gname]:.3e}{bound}'
              f'  [{mark}]')
    missing = [g for g in wanted if g not in worst]
    if missing:
        print(f'    {mesh} {run}: no cells in the groups {missing}')
    return ok


def limits(log):
    """The run hot: the line of the run that counts the ghost fills with the cut of the extrapolation and with the floor
    of the 302 wall temperature; both must be positive."""
    m = None
    for line in open(log, errors='replace'):
        m = re.search(r'fills with the extrapolation cut:\s*(\d+), with the 302 temperature at T_w/2:\s*(\d+)', line) or m
    if m is None:
        print('    paraT  hot  no line counting the limits of the ghost cells  [FAIL]')
        return False
    cut, floor = int(m.group(1)), int(m.group(2))
    ok = cut > 0 and floor > 0
    print(f'    paraT  hot  ghost fills with the cut: {cut}, with the floor T_w/2: {floor}  [{"ok  " if ok else "FAIL"}]')
    return ok


# ----------------------------------------------------------------------------------------------------------------------
# Replica of the face gradients of mean (Lib_Diffusive.f90, Lib_Metrics.f90 at bff5a51) for the mom run; with METRIC =
# inverse-mean the replica takes the inverse-mean metric, which diffusive-metric does not give to the interior faces
# ----------------------------------------------------------------------------------------------------------------------
def inverse3(A):
    """Inverse3 of Lib_Diffusive.f90: adjugate over the determinant."""
    B = np.empty((3, 3))
    B[0, 0] = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    B[0, 1] = -(A[0, 1] * A[2, 2] - A[0, 2] * A[2, 1])
    B[0, 2] = A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]
    B[1, 0] = -(A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
    B[1, 1] = A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]
    B[1, 2] = -(A[0, 0] * A[1, 2] - A[0, 2] * A[1, 0])
    B[2, 0] = A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]
    B[2, 1] = -(A[0, 0] * A[2, 1] - A[0, 1] * A[2, 0])
    B[2, 2] = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    det = A[0, 0] * B[0, 0] + A[0, 1] * B[1, 0] + A[0, 2] * B[2, 0]
    return B / det


def cell_metric(X, i, j, k):
    """M = A^-T of Compute_Metric_Tensor: rows of A are the mean edges along i, j, k."""
    N = [X[i - 1, j - 1, k - 1], X[i - 1, j - 1, k], X[i - 1, j, k - 1], X[i - 1, j, k],
         X[i, j - 1, k - 1], X[i, j - 1, k], X[i, j, k - 1], X[i, j, k]]
    A = np.empty((3, 3))
    A[0] = 0.25 * (N[4] - N[0] + N[5] - N[1] + N[6] - N[2] + N[7] - N[3])
    A[1] = 0.25 * (N[2] - N[0] + N[3] - N[1] + N[6] - N[4] + N[7] - N[5])
    A[2] = 0.25 * (N[1] - N[0] + N[3] - N[2] + N[5] - N[4] + N[7] - N[6])
    return np.linalg.inv(A).T


def extrapolation_ratio(X, face, i, j):
    """s = l0/l1 of the MUSCL distances behind face 1-4 of cell (i, j) of a 2-D block (Extrapolation_Ratio of
    Lib_Metrics.f90): mean edges across the face of the cell, of the next one and of the ghost cell whose nodes are
    extrapolated as 2 N0 - N1 (the edge of the cell itself)."""
    ni, nj, _ = dims(X)
    planes = {1: (0, 1, 2), 2: (ni, ni - 1, ni - 2), 3: (0, 1, 2), 4: (nj, nj - 1, nj - 2)}[face]
    e1 = np.zeros(3); e2 = np.zeros(3); eg = np.zeros(3)
    for t in ((j - 1, j) if face <= 2 else (i - 1, i)):
        for k in (0, 1):
            n0, n1, n2 = ((X[a, t, k] if face <= 2 else X[t, a, k]) for a in planes)
            e1 += n1 - n0
            e2 += n2 - n1
            eg += n1 - n0
    l = [float(np.linalg.norm(0.25 * e)) for e in (eg, e1, e2)]
    return (l[0] + l[1]) / (l[1] + l[2])


def replica(mesh, ns_path, eu_path, metric):
    """mom run, one block, outer faces 400 (no viscous flux; ghost values u1 + s (u1 - u2), s = l0/l1 of the MUSCL
    distances: Ghost_Wall), k faces 300 (copies)."""
    run = 'mom'
    blocks = build(mesh, run)
    if len(blocks) != 1:
        sys.exit('the replica is written for one block')
    X = blocks[0]
    ni, nj, nk = dims(X)
    sgn = face_signs(X)
    C = centres(X)
    f = field(mesh, run)
    U = np.zeros((ni + 2, nj + 2, 3))
    for i in range(ni):
        for j in range(nj):
            U[i + 1, j + 1] = f['u0'] + f['G'].dot(C[i, j, 0] - f['x0'])
    for j in range(1, nj + 1):
        U[0, j] = U[1, j] + extrapolation_ratio(X, 1, 1, j) * (U[1, j] - U[2, j])
        U[ni + 1, j] = U[ni, j] + extrapolation_ratio(X, 2, ni, j) * (U[ni, j] - U[ni - 1, j])
    for i in range(1, ni + 1):
        U[i, 0] = U[i, 1] + extrapolation_ratio(X, 3, i, 1) * (U[i, 1] - U[i, 2])
        U[i, nj + 1] = U[i, nj] + extrapolation_ratio(X, 4, i, nj) * (U[i, nj] - U[i, nj - 1])
    M = {}
    for i in range(1, ni + 1):
        for j in range(1, nj + 1):
            M[(i, j)] = cell_metric(X, i, j, 1)
    lam = -2.0 / 3.0 * MU
    force = np.zeros((ni + 2, nj + 2, 3))
    for d in (1, 2):
        rng_i = range(1, ni) if d == 1 else range(1, ni + 1)
        rng_j = range(1, nj + 1) if d == 1 else range(1, nj)
        for i in rng_i:
            for j in rng_j:
                L, R = (i, j), ((i + 1, j) if d == 1 else (i, j + 1))
                g = np.zeros((3, 3))                     # g[v, h]: difference of u_v along index direction h
                g[:, d - 1] = U[R] - U[L]
                t = 2 - d                                # the other in-plane direction (0-based)
                if d == 1:
                    g[:, t] = 0.25 * (U[i, j + 1] - U[i, j - 1] + U[i + 1, j + 1] - U[i + 1, j - 1])
                else:
                    g[:, t] = 0.25 * (U[i + 1, j] - U[i - 1, j] + U[i + 1, j + 1] - U[i - 1, j + 1])
                if metric == 'inverse-mean':
                    Mf = inverse3(0.5 * (inverse3(M[L]) + inverse3(M[R])))
                else:
                    Mf = 0.5 * (M[L] + M[R])
                grad = g.dot(Mf)                         # grad[v, m] = sum_h g[v, h] Mf[h, m]
                div = grad[0, 0] + grad[1, 1] + grad[2, 2]
                tau = MU * (grad + grad.T) + lam * div * np.eye(3)
                n, a = face_normal_area(X, sgn, d, i, j, 1)
                flux = a * tau.dot(n)
                force[L] += flux
                force[R] -= flux
    ns, eu = read_zones(ns_path), read_zones(eu_path)
    rhs = viscous_rhs(ns, eu, ['rhou', 'rhov', 'rhow'])[0]
    sc = scales(blocks, mesh, run)[0]
    worst, c = 0.0, 0
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            exp = force[i, j] / cell_volume(X, i, j, 1)
            worst = max(worst, float(np.max(np.abs(rhs[c] - exp))) / sc[c])
            c += 1
    good = worst <= TOL_REPLICA
    print(f'    {mesh} mom diffusive-metric = {metric}: largest |solver - replica| / scale = {worst:.3e}  '
          f'[{"ok  " if good else "FAIL"}]')
    return good


def main():
    what = sys.argv[1]
    if what == 'case':
        mesh, run = sys.argv[2], sys.argv[3]
        write_case(mesh, run, sys.argv[4] if len(sys.argv) > 4 else None)
    elif what == 'check':
        mesh, run, ns, eu, wanted = sys.argv[2:7]
        sys.exit(0 if check(mesh, run, ns, eu, wanted.split(',')) else 1)
    elif what == 'replica':
        mesh, ns, eu, metric = sys.argv[2:6]
        sys.exit(0 if replica(mesh, ns, eu, metric) else 1)
    elif what == 'limits':
        sys.exit(0 if limits(sys.argv[2]) else 1)
    elif what == 'work':
        mesh, run, ns, eu, wanted = sys.argv[2:7]
        sys.exit(0 if work(mesh, run, ns, eu, wanted.split(',')) else 1)
    elif what == 'open':
        mesh, run, ns, eu = sys.argv[2:6]
        sys.exit(0 if check_open(mesh, run, ns, eu) else 1)
    else:
        sys.exit('unknown request ' + what)


if __name__ == '__main__':
    main()
