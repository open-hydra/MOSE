"""Inputs and checks of the ViscousLinearPatch and ViscousWallPatch fast tests (patch test of the viscous fluxes).

A linear velocity field has a uniform viscous stress and a linear temperature field a uniform heat flux (mu and k constant):
the viscous force and the heat conduction of every cell are zero. A face gradient that is exact for linear fields gives
zero in every cell, to round-off, on any mesh: skewed, stretched, with jumps of cell size, with folds at the faces between
blocks, at points where 3 or 5 blocks meet, next to walls whose temperature varies along the wall. The difference
approximations of the field and of the coordinates must be the same (M. Vinokur, J. Comput. Phys. 81, 1989, 1-52, p. 26).
This is a test of consistency, the patch test of finite elements, not a test of the order of accuracy.

Meshes (unit coordinates, scaled by L = 1 cm and moved to x >= 1, y >= 2; 2-D meshes are one cell thick in z, with
symmetry (300) k faces and the same in-plane move for both node layers; the interior nodes are moved at random, with
fixed seeds, by up to 0.25 of the cell width on fold2, 0.15 of the cell size on box3d and sym, 0.15/n on star3 and star5
(n cells along a block side, polygon of unit circumradius), 0.2 of the cell size on rand1, not at all on jump, paraT and
paraQ):
  fold2  two blocks of 8 x 8 cells joined at x = 0.5 (block faces 2 and 1); the j lines of block 2 rise as y = t + 0.25 s t,
         so the grid lines fold at the face between the blocks;
  jump   one block of 8 x 8 cells whose sizes alternate in the ratio 1.5 in x and in y, sheared by x += 0.5 y, nodes not
         moved: the face gradients of inverse-mean and mean are not exact there, but for a linear field their errors
         cancel in the sums of the interior cells; on this mesh the check bites in the cells next to the walls;
  jump2  as jump with the sizes 1 and 2 (units of 1/16: the block is 0.75 wide), the first cells next to the faces 1
         and 3 half as wide as the second ones, unit 2^-7 m instead of 1 cm: every coordinate is a dyadic fraction, so
         that the ghost nodes extrapolated from the first three node planes are exact with any compiler (as
         3 N0 - 3 N1 + N2 they would fall on the face, a ghost cell of zero width);
  star3  a triangle cut into 3 blocks of 8 x 8 cells that meet at its centre (a point shared by 3 blocks);
  star5  a pentagon cut into 5 blocks of 8 x 8 cells that meet at its centre (a point shared by 5 blocks);
  box3d  two blocks of 6 x 6 x 6 cells joined at x = 0.5, the j lines of block 2 rising as y = t + 0.2 s t;
  sym    one block of 8 x 8 cells sheared by x += 0.35 y on a symmetry plane (300) at y = 2;
  paraT  one block of 8 x 8 parallelogram cells sheared by 20 degrees, walls 302 at y = 2 and 3, faces 400 at the sides,
  paraQ  the same with walls 301: controls, run ene only, where every face metric is exact (every cell has the same
         metric, T_w uniform along the 302 walls of paraT, q = 0 on the 301 walls of paraQ): inverse-mean and mean pass
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

Observable: one explicit step from the same state with the Navier-Stokes equations and one reference step (mom: the Euler
equations; ene, momw: Navier-Stokes with mu = k = 1e-30 and q = 0 on the 301 faces, because the solver refuses walls 301 and
302 with the Euler equations and the ghost cells must stay those of the walls; rk3: the Navier-Stokes step only, whose
inviscid part is zero for this state). The diagnostic output holds, per cell, the residual over the step, r, and the
step dt; r/dt of Navier-Stokes minus r/dt of the reference is the viscous right-hand side per unit volume. It is
compared with the scale mu |G| (sum of face areas)/volume (mom, momw) or k |g| (sum of face areas)/volume (ene, rk3);
threshold 1e-11. The cells are grouped by the faces they have: interior (no boundary face), connection (a face between
blocks; corner: two of them, the cells at a point shared by 3 or 5 blocks), symmetry (300), wall 301, wall 302.

ViscousLinearPatch checks every group but wall 302, plus the old metrics: on the rand1 mesh, run mom, diffusive-metric =
inverse-mean and mean must give, cell by cell, the viscous force of an independent replica of their face gradients
(written here from Lib_Diffusive.f90 and Lib_Metrics.f90 at bff5a51, with the ghost values of the 400 faces of
Ghost_Wall), within 1e-10 of the scale. ViscousWallPatch checks the cells with a 302 face: walls whose temperature varies
along the wall, on skewed cells (ene), and the velocity at a no-slip wall whose cells have nodes moved at random (momw); and
the cells next to walls 301 and 302 with T linear at uniform pressure (enep), and the count of the limits of the wall ghost
cells (hot).

Usage (from a test directory):
  python3 -B make_case.py case MESH RUN [METRIC]     INPUT/ic.tec; INPUT/bc-ns.txt, input-ns.ini for the Navier-Stokes
                                                      step; INPUT/bc-ref.txt, input-ref.ini, INPUT/transport-ref.dat
                                                      for the reference step; METRIC: inverse-mean or mean (no key
                                                      when absent: the default of the solver)
  python3 -B make_case.py check MESH RUN NS REF GROUPS prints the largest error of each group and exits with 0 when the
                                                      groups listed (comma-separated) are below the threshold
  python3 -B make_case.py replica MESH NS REF METRIC  compares the viscous force with the replica of METRIC (mom run)
  python3 -B make_case.py limits LOG                  the run hot: the counts of the floor and of the cut are positive
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
LREF = 1.0e-2                       # m
LREF2 = 2.0**-7                     # m, jump2: dyadic coordinates
OFFSET = np.array([1.0, 2.0, 0.0])  # unit coordinates: keep y > 0 away from the axis
TOL = 1.0e-11
TOL_REPLICA = 1.0e-10
FMT = '{: .17E}'
MESHES = ('fold2', 'jump', 'jump2', 'star3', 'star5', 'box3d', 'sym', 'paraT', 'paraQ', 'rand1')
RUNS = ('mom', 'ene', 'rk3', 'momw', 'enep', 'hot')
BC_CONNECT, BC_SYM, BC_Q, BC_T, BC_OPEN = 101, 300, 301, 302, 400


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


def outer_type(mesh, run, b, face):
    """Type of the outer face `face` of block b (faces 5, 6 of 2-D meshes are symmetry planes)."""
    if run == 'mom':
        if mesh == 'sym' and face == 3:
            return BC_SYM
        return BC_OPEN
    if run == 'momw':                    # one wall, the plane y = 2 cm of rand1 and box3d
        return BC_T if face == 3 else BC_OPEN
    if run == 'rk3':
        mix = {0: (BC_T, BC_T), 1: (BC_Q, BC_T), 2: (BC_OPEN, BC_Q), 3: (BC_T, BC_OPEN), 4: (BC_Q, BC_Q)}
        return mix[b][0 if face == 2 else 1]
    # ene
    if mesh == 'fold2':
        return {(0, 1): BC_T, (0, 3): BC_Q, (0, 4): BC_T, (1, 2): BC_OPEN, (1, 3): BC_T, (1, 4): BC_Q}[(b, face)]
    if mesh in ('jump', 'jump2', 'rand1'):
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


def tvol(v, i1, i2, i3, i4):
    a, b, c, d = v[i1], v[i2], v[i3], v[i4]
    return abs(np.dot(b - a, np.cross(c - a, d - a))) / 6.0


def cell_volume(X, i, j, k):
    """Compute_Metric_Tensor: five tetrahedra on the nodes N1..N8."""
    N = [X[i - 1, j - 1, k - 1], X[i - 1, j - 1, k], X[i - 1, j, k - 1], X[i - 1, j, k],
         X[i, j - 1, k - 1], X[i, j - 1, k], X[i, j, k - 1], X[i, j, k]]
    v = [None, N[0], N[4], N[2], N[6], N[1], N[5], N[3], N[7]]    # vx(1..8) = N1, N5, N3, N7, N2, N6, N4, N8
    return (tvol(v, 1, 2, 3, 5) + tvol(v, 2, 4, 3, 8) + tvol(v, 5, 8, 6, 2) + tvol(v, 5, 7, 8, 3) + tvol(v, 5, 8, 2, 3))


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
GP_HOT = 4.0e4                             # run hot: pressure gradient along x, Pa/m (50 Pa per cell of paraT)


def field(mesh, run):
    if run == 'momw':
        # zero on the wall y = 2 cm, with shear and dilatation: u = (a, b, c) (y - 2 cm)
        G = np.array([[0.0, 120.0, 0.0], [0.0, 60.0, 0.0], [0.0, 40.0 if mesh == 'box3d' else 0.0, 0.0]])
        x0 = X0.copy()
        x0[1] = LREF * OFFSET[1]
        return {'G': G, 'u0': np.zeros(3), 'x0': x0}
    if run == 'mom':
        G = GSYM if mesh == 'sym' else (G3 if mesh == 'box3d' else G2)
        u0 = {'box3d': np.array([1.0, 0.5, 0.3]), 'sym': np.array([1.0, 0.0, 0.0])}.get(mesh, U0)
        x0 = X0.copy()
        if mesh == 'sym':
            x0[1] = LREF * OFFSET[1]           # u_y odd about the symmetry plane y = 2 cm
        return {'G': G, 'u0': u0, 'x0': x0}
    g = {'sym': gsym, 'box3d': g3, 'paraT': gparaT, 'paraQ': gparaQ}.get(mesh, g2)
    return {'g': g, 'x0': X0.copy()}


def state(mesh, run, x):
    f = field(mesh, run)
    if run in ('mom', 'momw'):
        u = f['u0'] + f['G'].dot(x - f['x0'])
        return P_LOW / (RGAS * T0), u, P_LOW
    if run == 'hot':
        p = P_LOW + GP_HOT * (x[0] - LREF * OFFSET[0])
        return p / (RGAS * T0), np.zeros(3), p
    T = T0 + f['g'].dot(x - f['x0'])
    if run == 'ene':
        return RHO0, np.zeros(3), RHO0 * RGAS * T
    return P_LOW / (RGAS * T), np.zeros(3), P_LOW


# ----------------------------------------------------------------------------------------------------------------------
# Case files
# ----------------------------------------------------------------------------------------------------------------------
def write_case(mesh, run, metric):
    blocks = build(mesh, run)
    conn = connections(blocks)
    os.makedirs('INPUT', exist_ok=True)
    for p in ('INPUT/ic.tec', 'INPUT/bc.txt', 'INPUT/bc-ns.txt', 'INPUT/bc-ref.txt', 'INPUT/transport-ref.dat'):
        if os.path.lexists(p):
            os.remove(p)
    for b, X in enumerate(blocks):
        ni, nj, nk = dims(X)
        vmin = min(cell_volume(X, i, j, k) for k in range(1, nk + 1) for j in range(1, nj + 1) for i in range(1, ni + 1))
        if vmin <= 0.0:
            sys.exit(f'block {b + 1}: a cell with volume {vmin}')
    with open('INPUT/ic.tec', 'w') as out:
        out.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"\n')
        for b, X in enumerate(blocks):
            ni, nj, nk = dims(X)
            out.write(f'ZONE T="Block{b + 1}", I={ni + 1}, J={nj + 1}, K={nk + 1}, '
                      f'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED)\n')
            for axis in range(3):
                out.write(''.join(' ' + FMT.format(X[i, j, k, axis]) + '\n'
                                  for k in range(nk + 1) for j in range(nj + 1) for i in range(ni + 1)))
            C = centres(X)
            cells = [state(mesh, run, C[i, j, k]) for k in range(nk) for j in range(nj) for i in range(ni)]
            cols = ([c[0] for c in cells], [c[1][0] for c in cells], [c[1][1] for c in cells],
                    [c[1][2] for c in cells], [c[2] for c in cells])
            for col in cols:
                out.write(''.join(' ' + FMT.format(v) + '\n' for v in col))
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
                            tw = T0 if run == 'momw' else T0 + f['g'].dot(xf - f['x0'])
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
                            q = into * KCOND * f['g'].dot(nrm) if (run != 'mom' and tag == 'ns') else 0.0
                            rows.append(head + f'{t:8d}\n' + f'    {q:.17E},    0.000000E+00,    0.000000E+00,\n')
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
    for tag, eq in (('ns', 'navier-stokes'), ('ref', 'euler' if run == 'mom' else 'navier-stokes')):
        with open(f'input-{tag}.ini', 'w') as out:
            out.write('[MOSE-Parameters]\niter-threshold = 1\nnewrun = true\nres-threshold = 1e-300\n\n'
                      f'[MOSE-Physics]\nequations = {eq}\nschmidt = 0.7\n\n'
                      '[MOSE-Numerics]\ncfl = 0.5\nvnn = 0.3\ntime-accurate = false\n'
                      f'time-scheme = {scheme}\nintegration-variables = cons\nspace-reconstruction = MUSCL\n'
                      'flux-limiter = vanleer\nriemann-solver = HLLC\n' + key)
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
                    elif run == 'rk3':
                        g.append('centre' if (i <= 4 and j <= 4) else None)
                    elif BC_T in kinds:
                        g.append('wall302')
                    elif BC_Q in kinds:
                        g.append('wall301')
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
    f = field(mesh, run)
    coef = MU * float(np.linalg.norm(f['G'])) if run in ('mom', 'momw') else KCOND * float(np.linalg.norm(f['g']))
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
    comps = ['rhou', 'rhov', 'rhow'] if run in ('mom', 'momw') else ['rhoe']
    rhs = viscous_rhs(read_zones(ns_path), None if run == 'rk3' else read_zones(eu_path), comps)
    grp = groups(blocks, mesh, run, conn)
    sc = scales(blocks, mesh, run)
    worst = {}
    for b in range(len(blocks)):
        for c, gname in enumerate(grp[b]):
            if gname is None:
                continue
            e = float(np.max(np.abs(rhs[b][c]))) / sc[b][c]
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
# Replica of the face gradients of inverse-mean and mean (Lib_Diffusive.f90, Lib_Metrics.f90 at bff5a51) for the mom run
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
    else:
        sys.exit('unknown request ' + what)


if __name__ == '__main__':
    main()
