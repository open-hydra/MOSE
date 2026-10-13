"""Inputs of the WallDistance fast test, and the check of its output.

The distance d to the nearest wall enters the Spalart-Allmaras model, and it must be the minimum distance from the cell centre
to the wall surface, which near a convex corner is the distance to the corner (NASA Turbulence Modeling Resource,
Spalart-Allmaras page). The solver does not write d, but it can be read back from the model: with the fluid at rest
(vorticity zero) and nu~ uniform, the right-hand side of the equation for rho nu~ in a cell away from the walls has no
convection and no diffusion, and it is (Spalart-Allmaras without ft2, as in MOSE)
    rho (nu~/d)^2 [ cb1 fv2/kappa^2 - cw1 fw(r) ],     r = min(1/fv2, 10),
since S~ = nu~ fv2/(kappa^2 d^2). One explicit Euler step with local time stepping writes, in OUTPUT/diagnostic.tec, the
increment of rho nu~ and the step of each cell; their ratio is the right-hand side, and so d. Air of
test/common/Air-viscous at 300 K and 1 kPa, nu~ = nu (chi = 1, so fv2 > 0 and r = 10): at this low density the round-off
of the fluxes of the gas at rest stays far below the source.

Three single-block cases, the wall on the j = 1 face, faces 300 elsewhere (and 200 on the k faces of the wedge):
  corner  2-D, a convex corner of 90 degrees between two walls: the image of the upper half plane under w = z^(3/2), so
          that the fluid takes three quadrants and the solid the fourth; in the second quadrant the nearest wall point is
          the corner.
  arc     2-D, a convex circular wall (radius 1, 12 straight faces over 120 degrees) under a mesh whose radial lines lean
          (theta shifts with r), so that the cell centres are not above the face centres.
  lip     the lip of a nozzle on a 1-degree axisymmetric wedge: in the meridional plane the nozzle wall, curved, meets the
          lip face at a convex corner of 90 degrees at x = 0, r = 2; in the jet below the lip the nearest wall point is the
          lip. The distance is the one in the meridional plane, from the wall curve.
Nine more cases have periodic boundaries, across which the wall goes on: the distance is the minimum over the periodic
images of the wall too.
  periodic  2-D channel periodic in x (connections 201 between the i faces, period 1), wall on the j = 1 face for x in
            [0, 0.5] only, symmetry on the rest; near x = 1 the nearest wall is the image of the wall that starts at x = 1.
  pipe      the same on a 1-degree axisymmetric wedge: pipe of radius 0.5 periodic in x, wall on the outer face for x in
            [0, 0.5] only; the distance is the meridional one. The axis nodes are at 1e-6 of the radial size of the
            first cell, as in the meshes of MOSE, so that no face has zero area.
  sector    3-D sector of 60 degrees of an annulus (radii 0.5 and 1), closed by the rotational periodicity (type 200 on
            the k faces), wall on the inner face for the first half of the sector only; near the end of the sector the
            nearest wall is the image of the wall in the next sector, turned by 60 degrees.
  noisy     the channel of periodic, with the nodes of the plane i = NI moved by up to 1e-8 in x and y (differently for
            every j): the two sides of the periodic connection match only to within the precision of the mesh file, as in
            a mesh written in single precision. The translation of the images then differs from the period by up to
            1.4e-8, and the distance by as much: the tolerance is 1e-6, for a distance of at least 0.03.
  span      3-D slab (1 x 0.5 x 0.6) closed in k by type 200 with the plane k = NK the plane k = 0 translated by 0.6 in z
            (periodic in k), wall on the j = 1 face for z in [0, 0.3] only; near z = 0.6 the nearest wall is the image of
            the wall that starts at z = 0.6.
  thin      3-D channel periodic in x and in z (connections 201 between the i faces and between the k faces), both periods
            THIN_P = 0.005, a hundredth of its height THIN_H = 0.5, wall on the whole j = 1 face: the images of the wall are
            the same plane, so the distance is y of the centre, but the periodic translations are many times shorter than
            the mesh.
  far       the channel of periodic with every node moved by FAR = 10 in x, ten periods away from the origin.
  staggered 2-D staggered row of plane blades: one block of 24 x 10 cells, x in [-0.5, 1.5] and y in [1.7 x, 1.7 x + 1],
            periodic in y (connections 201 between the j faces, translation (0, 1, 0)) ahead of and behind the blade, wall
            on the j = 1 face for x in [0, 1]: a blade from (0, 0) to (1, 1.7), 1.7 periods long along y, symmetry on the
            j = 10 face above it. The wall faces lie in different periods, and the nearest wall of a cell can be the
            blade or its images above and below.
  oblique   3-D block of 4 x 4 x 2 cells, the parallelepiped on the edges E1 = (1, 0, 0), E2 = (-1/2, sqrt(3)/2, 0)
            and E3 = (-1/2, -sqrt(3)/2, OBLIQUE_EPS = 0.02), periodic along the three (connections 201 between the i
            faces, between the k faces, and between the j faces for the second half of the cells along i), wall on the
            j = 1 face for the first half, symmetry on the j = 4 face above it. The three translations are nearly
            dependent: no pairwise step shortens them, but E1 + E2 + E3 = (0, 0, 0.02) is fifty times shorter, and the
            images of the wall are planes 0.02 apart. The reference takes the images along the reduced basis
            (0, 0, 0.02), E1, -E2.
The reference is the exact distance from the cell centre (mean of the nodes; in the wedge in the coordinates x and
r = sqrt(y^2 + z^2)) to the wall of the mesh (straight faces between the wall nodes, plane in every case) and to its
periodic images, by a scan of all faces. It is checked on the cells that do not touch the wall (in one step the wall
reaches only the cells next to it, through the diffusion across the wall face and the central difference of the cb2
term): the largest relative difference must be below 1e-9 (1e-6 in the case noisy).

Usage (from the test directory):
  python3 -B make_case.py case NAME            INPUT/ic.tec and INPUT/bc.txt of corner, arc, lip, periodic, pipe, sector,
                                               noisy, span, thin, far, staggered or oblique
  python3 -B make_case.py error NAME FILE      prints the largest relative error of d read from FILE (a diagnostic.tec)
  python3 -B make_case.py check N1=E1 ...      exit status 0 when the error of every case is below its tolerance
  python3 -B make_case.py case twist           INPUT/ic.tec and INPUT/bc.txt of the 3-D case with twisted wall faces
  python3 -B make_case.py twist LOG            exit status 0 when LOG reports the twist of the wall faces of the twist case

The twist case is a 3-D block whose wall, on the j = 1 face, is the surface y = 2 (x - 0.5)(z - 0.5): every wall face is a
bilinear surface with twist 2 dx dz = 1/32 along y, and the four triangles of a face lie within 1/512 of it. MOSE reports
at start-up the largest of that bound over the distance of the cell next to the face; it must match the value computed
here (the distance of that cell to the four triangles of the faces, as MOSE takes it) within 1e-3, the printed precision.
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P0 = 300.0, 1.0e3              # low density: the round-off of the fluxes at rest stays far below the source
RHO = P0 / (RGAS * T0)
CHI = 1.0                           # nu~/nu: fv2 > 0, r = min(1/fv2, 10) = 10
NUT = CHI * MU / RHO
CB1, CB2, SIGMA, KAPPA = 0.1355, 0.622, 2.0 / 3.0, 0.41
CW2, CW3, CV1 = 0.3, 2.0, 7.1
CW1 = CB1 / KAPPA ** 2 + (1.0 + CB2) / SIGMA
DELTA = math.radians(1.0)           # wedge angle of the lip case
LZ = 0.1                            # thickness of the 2-D cases
TOL = 1.0e-9
TOLS = {"noisy": 1.0e-6}            # the two sides of the periodic connection match to 1e-8 only
NOISE = 1.0e-8                      # case noisy: the largest displacement of the nodes of the plane i = NI
SKIP = 1                            # rows of cells next to the wall left out of the check (j = 1 .. SKIP)

SIZES = {"corner": (24, 16), "arc": (12, 12), "lip": (24, 16)}


def corner_map(xi, eta):
    """w = z^(3/2) on the closed upper half plane; the real axis goes to the two walls, exactly."""
    if eta == 0.0:
        return (xi ** 1.5, 0.0) if xi >= 0.0 else (0.0, -((-xi) ** 1.5))
    rho, phi = math.hypot(xi, eta), math.atan2(eta, xi)
    return rho ** 1.5 * math.cos(1.5 * phi), rho ** 1.5 * math.sin(1.5 * phi)


def nodes(name):
    """Node coordinates (x, y) of the 2-D cases, or (x, r) of the wedge, as a list of rows j = 0..NJ of NI + 1 points."""
    ni, nj = SIZES[name]
    pts = []
    for j in range(nj + 1):
        row = []
        for i in range(ni + 1):
            if name == "corner":
                row.append(corner_map(-1.0 + 2.0 * i / ni, (j / nj) ** 1.2))
            elif name == "arc":
                r = 1.0 + 0.6 * (j / nj) ** 1.2
                th = math.radians(150.0 - 10.0 * i) + 0.5 * (r - 1.0)
                row.append((r * math.cos(th), r * math.sin(th)))
            else:
                u, v = corner_map(-1.0 + 2.0 * i / ni, (j / nj) ** 1.2)
                big_x, big_y = -u, -v                   # solid in the second quadrant: nozzle wall Y = 0, X < 0
                if big_x < 0.0:                         # lip face X = 0, Y > 0
                    big_y += 0.3 * big_x * big_x        # curved nozzle wall
                row.append((big_x, 2.0 + big_y))        # x, r: lip at x = 0, r = 2
        pts.append(row)
    return pts


def fmt(v):
    return f"{v:.16E}"


def coords(name):
    """The node coordinates as MOSE reads them: (x, y, z) of the k = 0 and k = 1 planes, through the text of ic.tec."""
    pts = nodes(name)
    out = {}
    for k in range(2):
        for j, row in enumerate(pts):
            for i, (a, b) in enumerate(row):
                if name == "lip":
                    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
                    xyz = (a, b * c, (-1.0, 1.0)[k] * b * s)
                else:
                    xyz = (a, b, (0.0, LZ)[k])
                out[(i, j, k)] = tuple(float(fmt(v)) for v in xyz)
    return out


def sorted_bc(text):
    """The records of bc.txt face by face, block by block (the documented order), each with its data lines; the order
    within a face is kept. A record line holds six integers; the lines after it, up to the next one, are its data."""
    recs = []
    for line in text.splitlines(True):
        tok = line.split()
        if len(tok) == 6 and all(t.lstrip("-").isdigit() for t in tok):
            recs.append([(int(tok[0]), int(tok[4])), line])
        else:
            recs[-1][1] += line
    recs.sort(key=lambda r: r[0])
    return "".join(r[1] for r in recs)


def write_case(name):
    ni, nj = SIZES[name]
    xyz = coords(name)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    nc = ni * nj
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fmt(xyz[(i, j, k)][c])}\n")
        for v in (RHO, 0.0, 0.0, 0.0, P0, RHO * NUT):
            f.write(f" {fmt(v)}\n" * nc)

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    out = [rec(1, j, 1, 300) for j in range(1, nj + 1)] + [rec(ni, j, 2, 300) for j in range(1, nj + 1)]
    for i in range(1, ni + 1):
        out += [rec(i, 1, 3, 301), "    0.000000E+00,    0.000000E+00,    0.000000E+00,\n"]
    out += [rec(i, nj, 4, 300) for i in range(1, ni + 1)]
    kface = 200 if name == "lip" else 0
    for face in (5, 6):
        out += [rec(i, j, face, kface) for j in range(1, nj + 1) for i in range(1, ni + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write(sorted_bc("".join(out)))


def plane(name, p):
    """(x, y) of a node in the plane where the distance is measured: x, r in the wedge."""
    x, y, z = p
    return (x, math.hypot(y, z)) if name == "lip" else (x, y)


def segment_distance(p, a, b):
    """Distance from p to the segment ab, and whether the nearest point is an end of it."""
    ax, ay = b[0] - a[0], b[1] - a[1]
    px, py = p[0] - a[0], p[1] - a[1]
    l2 = ax * ax + ay * ay
    t = 0.0 if l2 == 0.0 else (px * ax + py * ay) / l2
    end = t <= 0.0 or t >= 1.0
    t = min(max(t, 0.0), 1.0)
    return math.hypot(px - t * ax, py - t * ay), end


def geometry(name):
    """Per cell (i, j), 1-based: the reference distance, whether its nearest wall point is a wall node (a corner), and the
    distance of the old method (cell centre to the nearest wall-face centre, in 3-D)."""
    ni, nj = SIZES[name]
    xyz = coords(name)
    wall = [plane(name, xyz[(i, 0, 0)]) for i in range(ni + 1)]
    centres3 = [[sum(xyz[(i + a, 0, k)][c] for a in (0, 1) for k in (0, 1)) / 4.0 for c in range(3)]
                for i in range(ni)]
    out = {}
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            corners = [xyz[(i - 1 + a, j - 1 + b, k)] for a in (0, 1) for b in (0, 1) for k in (0, 1)]
            pc = [plane(name, q) for q in corners]
            p = (sum(q[0] for q in pc) / 8.0, sum(q[1] for q in pc) / 8.0)
            best, at_node = None, False
            for s in range(ni):
                d, end = segment_distance(p, wall[s], wall[s + 1])
                if best is None or d < best:
                    best, at_node = d, end
            c3 = [sum(q[c] for q in corners) / 8.0 for c in range(3)]
            old = min(math.sqrt(sum((c3[c] - f[c]) ** 2 for c in range(3))) for f in centres3)
            out[(i, j)] = (best, at_node, old)
    return out


PERIODIC = {"periodic": (16, 10, 1), "pipe": (16, 10, 1), "sector": (4, 6, 12), "noisy": (16, 10, 1),
            "span": (4, 6, 12), "thin": (2, 8, 2), "far": (16, 10, 1), "staggered": (24, 10, 1), "oblique": (4, 4, 2)}
OBLIQUE_EPS = 0.02                  # case oblique: the edges E1, E2, E3 of the block, E1 + E2 + E3 = (0, 0, eps)
OBLIQUE = ((1.0, 0.0, 0.0), (-0.5, 0.5 * math.sqrt(3.0), 0.0), (-0.5, -0.5 * math.sqrt(3.0), OBLIQUE_EPS))
THIN_P, THIN_H = 0.005, 0.5         # case thin: periods in x and z, and height
FAR = 10.0                          # case far: the shift in x of the channel of periodic
STAGGER, BLADE = 1.7, (7, 18)       # case staggered: slope of the blade, and its cells on the j = 1 face (x in [0, 1])
SECTOR = math.radians(60.0)
SPAN = 0.6                          # period in z of the case span
WALL_LINE = "    0.000000E+00,    0.000000E+00,    0.000000E+00,\n"


def periodic_nodes(name):
    """Node coordinates (x, y, z) as MOSE reads them, by (i, j, k)."""
    ni, nj, nk = PERIODIC[name]
    out = {}
    for k in range(nk + 1):
        for j in range(nj + 1):
            for i in range(ni + 1):
                x = i / ni if name != "sector" else 0.1 * i
                if name == "periodic":
                    xyz = (x, 0.5 * (j / nj) ** 1.5, (0.0, LZ)[k])
                elif name == "noisy":
                    y = 0.5 * (j / nj) ** 1.5
                    if i == ni:                                  # the far side of the period, matched to 1e-8 only
                        x, y = x + NOISE * math.sin(1.0 + 3.0 * j), y + NOISE * math.cos(2.0 + 5.0 * j)
                    xyz = (x, y, (0.0, LZ)[k])
                elif name == "span":
                    xyz = (0.25 * i, 0.5 * (j / nj) ** 1.5, SPAN * k / nk)
                elif name == "thin":
                    xyz = (THIN_P * i / ni, THIN_H * (j / nj) ** 1.5, THIN_P * k / nk)
                elif name == "far":
                    xyz = (x + FAR, 0.5 * (j / nj) ** 1.5, (0.0, LZ)[k])
                elif name == "staggered":
                    xs = -0.5 + 2.0 * i / ni
                    xyz = (xs, STAGGER * xs + j / nj, (0.0, LZ)[k])
                elif name == "oblique":
                    xyz = tuple(i / ni * OBLIQUE[0][c] + j / nj * OBLIQUE[1][c] + k / nk * OBLIQUE[2][c]
                                for c in range(3))
                elif name == "pipe":
                    r = 0.5 * (j / nj) ** 0.8 if j > 0 else 1.0e-6 * 0.5 / nj ** 0.8   # axis at j = 0, wall at r = 0.5
                    xyz = (x, r * math.cos(0.5 * DELTA), (-1.0, 1.0)[k] * r * math.sin(0.5 * DELTA))
                else:
                    r = 0.5 + 0.5 * (j / nj) ** 1.5             # inner wall at j = 0
                    th = SECTOR * k / nk
                    xyz = (x, r * math.cos(th), r * math.sin(th))
                out[(i, j, k)] = tuple(float(fmt(v)) for v in xyz)
    return out


def write_periodic_case(name):
    ni, nj, nk = PERIODIC[name]
    xyz = periodic_nodes(name)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    nc = ni * nj * nk
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fmt(xyz[(i, j, k)][c])}\n")
        for v in (RHO, 0.0, 0.0, 0.0, P0, RHO * NUT):
            f.write(f" {fmt(v)}\n" * nc)

    def rec(i, j, k, face, t):
        out = '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)
        return out + WALL_LINE if t == 301 else out

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    out = []
    for k in range(1, nk + 1):
        for j in range(1, nj + 1):
            if name in ("sector", "span", "staggered"):
                out += [rec(1, j, k, 1, 300), rec(ni, j, k, 2, 300)]
            else:
                out += [rec(1, j, k, 1, 201), con(1, ni, j, k, 2, 1, 0, 0, 1)]
                out += [rec(ni, j, k, 2, 201), con(1, 1, j, k, 1, 1, 0, 0, 1)]
        for i in range(1, ni + 1):
            if name in ("periodic", "noisy", "far"):
                out += [rec(i, 1, k, 3, 301 if i <= ni // 2 else 300), rec(i, nj, k, 4, 300)]
            elif name == "thin":
                out += [rec(i, 1, k, 3, 301), rec(i, nj, k, 4, 300)]
            elif name == "staggered" and BLADE[0] <= i <= BLADE[1]:
                out += [rec(i, 1, k, 3, 301), rec(i, nj, k, 4, 300)]
            elif name == "staggered":
                out += [rec(i, 1, k, 3, 201), con(1, i, nj, k, 4, 1, 0, 0, 1)]
                out += [rec(i, nj, k, 4, 201), con(1, i, 1, k, 3, 1, 0, 0, 1)]
            elif name == "oblique" and i <= ni // 2:
                out += [rec(i, 1, k, 3, 301), rec(i, nj, k, 4, 300)]
            elif name == "oblique":
                out += [rec(i, 1, k, 3, 201), con(1, i, nj, k, 4, 1, 0, 0, 1)]
                out += [rec(i, nj, k, 4, 201), con(1, i, 1, k, 3, 1, 0, 0, 1)]
            elif name == "pipe":
                out += [rec(i, 1, k, 3, 300), rec(i, nj, k, 4, 301 if i <= ni // 2 else 300)]
            else:
                out += [rec(i, 1, k, 3, 301 if k <= nk // 2 else 300), rec(i, nj, k, 4, 300)]
    kface = {"periodic": 0, "pipe": 200, "sector": 200, "noisy": 0, "span": 200, "thin": 201, "far": 0,
             "staggered": 0, "oblique": 201}[name]
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            if kface == 201:
                out += [rec(i, j, 1, 5, 201), con(1, i, j, nk, 6, 1, 0, 0, 1)]
                out += [rec(i, j, nk, 6, 201), con(1, i, j, 1, 5, 1, 0, 0, 1)]
            else:
                out += [rec(i, j, 1, 5, kface), rec(i, j, nk, 6, kface)]
    with open("INPUT/bc.txt", "w") as f:
        f.write(sorted_bc("".join(out)))


def triangle_distance(p, a, b, c):
    """Distance from p to the triangle abc in 3-D: the plane when the projection falls inside, else the nearest side."""
    def sub(u, v):
        return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

    def dot(u, v):
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    def cross(u, v):
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]

    def seg(a, b):
        ab, ap = sub(b, a), sub(p, a)
        t = min(max(dot(ap, ab) / dot(ab, ab), 0.0), 1.0)
        return math.sqrt(sum((ap[i] - t * ab[i]) ** 2 for i in range(3)))

    n = cross(sub(b, a), sub(c, a))
    n2 = dot(n, n)
    if all(dot(cross(sub(v, u), sub(p, u)), n) >= 0.0 for u, v in ((a, b), (b, c), (c, a))):
        return abs(dot(sub(p, a), n)) / math.sqrt(n2)
    return min(seg(a, b), seg(b, c), seg(c, a))


def periodic_reference(name):
    """Exact distance of every cell (i, j, k), 1-based, that does not touch the wall, to the wall and its images."""
    ni, nj, nk = PERIODIC[name]
    xyz = periodic_nodes(name)
    out = {}
    if name == "thin":
        # the wall and its images are the plane y = 0: the distance is y of the centre
        for k in range(1, nk + 1):
            for j in range(2, nj + 1):
                for i in range(1, ni + 1):
                    corners = [xyz[(i - 1 + a, j - 1 + b, k - 1 + c)] for a in (0, 1) for b in (0, 1) for c in (0, 1)]
                    out[(i, j, k)] = sum(q[1] for q in corners) / 8.0
        return out
    if name == "staggered":
        # the blade and its images shifted by k in y, k = -3 .. 3 (the nearest wall of every cell is among them)
        wall = [(plane("corner", xyz[(i, 0, 0)]), plane("corner", xyz[(i + 1, 0, 0)]))
                for i in range(BLADE[0] - 1, BLADE[1])]
        for j in range(2, nj + 1):
            for i in range(1, ni + 1):
                corners = [xyz[(i - 1 + a, j - 1 + b, k)] for a in (0, 1) for b in (0, 1) for k in (0, 1)]
                pc = [plane("corner", q) for q in corners]
                p = (sum(q[0] for q in pc) / 8.0, sum(q[1] for q in pc) / 8.0)
                out[(i, j, 1)] = min(segment_distance((p[0], p[1] - shift), a, b)[0]
                                     for shift in range(-3, 4) for a, b in wall)
        return out
    if name == "oblique":
        # the wall faces (j = 0, the first half along i) and their images by sum n_l R_l, |n_l| <= 3, R the reduced
        # basis (0, 0, eps) = E1 + E2 + E3, E1, -E2: the coordinates of p - q (p a cell centre, q a wall point) along R
        # are at most 1, 2, 2, and those of a nearest lattice vector differ from them by at most beta = 1/2, 1, 1
        faces = [[xyz[(i, 0, k)], xyz[(i + 1, 0, k)], xyz[(i + 1, 0, k + 1)], xyz[(i, 0, k + 1)]]
                 for k in range(nk) for i in range(ni // 2)]
        e1, e2, e3 = OBLIQUE
        red = [[e1[c] + e2[c] + e3[c] for c in range(3)], list(e1), [-v for v in e2]]
        shifts = [[a * red[0][c] + b * red[1][c] + g * red[2][c] for c in range(3)]
                  for a in range(-3, 4) for b in range(-3, 4) for g in range(-3, 4)]
        for k in range(1, nk + 1):
            for j in range(2, nj + 1):
                for i in range(1, ni + 1):
                    corners = [xyz[(i - 1 + a, j - 1 + b, k - 1 + c)] for a in (0, 1) for b in (0, 1) for c in (0, 1)]
                    p = [sum(q[c] for q in corners) / 8.0 for c in range(3)]
                    out[(i, j, k)] = min(min(triangle_distance([p[c] - s[c] for c in range(3)], f[0], f[1], f[2]),
                                             triangle_distance([p[c] - s[c] for c in range(3)], f[0], f[2], f[3]))
                                         for s in shifts for f in faces)
        return out
    if name in ("periodic", "pipe", "noisy", "far"):
        jw = nj if name == "pipe" else 0
        wall = [((plane("lip" if name == "pipe" else "corner", xyz[(i, jw, 0)])),
                 (plane("lip" if name == "pipe" else "corner", xyz[(i + 1, jw, 0)]))) for i in range(ni // 2)]
        for j in range(1, nj + 1):
            if (name != "pipe" and j == 1) or (name == "pipe" and j == nj):
                continue
            for i in range(1, ni + 1):
                corners = [xyz[(i - 1 + a, j - 1 + b, k)] for a in (0, 1) for b in (0, 1) for k in (0, 1)]
                pc = [plane("lip" if name == "pipe" else "corner", q) for q in corners]
                p = (sum(q[0] for q in pc) / 8.0, sum(q[1] for q in pc) / 8.0)
                out[(i, j, 1)] = min(segment_distance((p[0] - shift, p[1]), a, b)[0]
                                     for shift in (-2.0, -1.0, 0.0, 1.0, 2.0) for a, b in wall)
        return out
    faces = []
    for k in range(nk // 2):
        for i in range(ni):
            faces.append([xyz[(i, 0, k)], xyz[(i + 1, 0, k)], xyz[(i + 1, 0, k + 1)], xyz[(i, 0, k + 1)]])
    images = []
    for m in range(-3, 4):
        cm, sm = math.cos(m * SECTOR), math.sin(m * SECTOR)
        for f in faces:
            if name == "span":
                images.append([(q[0], q[1], q[2] + m * SPAN) for q in f])
            else:
                images.append([(q[0], q[1] * cm - q[2] * sm, q[1] * sm + q[2] * cm) for q in f])
    for k in range(1, nk + 1):
        for j in range(2, nj + 1):
            for i in range(1, ni + 1):
                corners = [xyz[(i - 1 + a, j - 1 + b, k - 1 + c)] for a in (0, 1) for b in (0, 1) for c in (0, 1)]
                p = [sum(q[c] for q in corners) / 8.0 for c in range(3)]
                out[(i, j, k)] = min(min(triangle_distance(p, f[0], f[1], f[2]), triangle_distance(p, f[0], f[2], f[3]))
                                     for f in images)
    return out


def periodic_error(name, path):
    data, ni, nj = read_zone(path)
    ref = periodic_reference(name)
    f = factor()
    nip, njp, nkp = PERIODIC[name]
    assert (ni, nj) == (nip, njp), (ni, nj)
    err = 0.0
    for (i, j, k), d0 in ref.items():
        n = ((k - 1) * nj + (j - 1)) * ni + (i - 1)
        rhs = data["rhoNuTilde"][n] / data["dt"][n]
        d = math.sqrt(RHO * NUT * NUT * f / rhs) if rhs < 0.0 else float("inf")
        err = max(err, abs(d / d0 - 1.0))
    print(f"{err:.3e} {len(ref)} 0")


TWIST = (8, 6, 8)


def twist_nodes():
    """Node coordinates (x, y, z) of the twist case as MOSE reads them, by (i, j, k)."""
    ni, nj, nk = TWIST
    out = {}
    for k in range(nk + 1):
        for j in range(nj + 1):
            for i in range(ni + 1):
                x, z = i / ni, k / nk
                eta = 1.0e-3 * (1.6 ** j - 1.0) / 0.6
                out[(i, j, k)] = tuple(float(fmt(v)) for v in (x, 2.0 * (x - 0.5) * (z - 0.5) + eta, z))
    return out


def write_twist_case():
    ni, nj, nk = TWIST
    xyz = twist_nodes()
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    nc = ni * nj * nk
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fmt(xyz[(i, j, k)][c])}\n")
        for v in (RHO, 0.0, 0.0, 0.0, P0, RHO * NUT):
            f.write(f" {fmt(v)}\n" * nc)
    out = []
    for k in range(1, nk + 1):
        for j in range(1, nj + 1):
            out += ['%9d %7d %7d %7d %7d %7d\n' % (1, 1, j, k, 1, 300), '%9d %7d %7d %7d %7d %7d\n' % (1, ni, j, k, 2, 300)]
        for i in range(1, ni + 1):
            out += ['%9d %7d %7d %7d %7d %7d\n' % (1, i, 1, k, 3, 301), WALL_LINE]
            out += ['%9d %7d %7d %7d %7d %7d\n' % (1, i, nj, k, 4, 300)]
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            out += ['%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, 5, 300), '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, nk, 6, 300)]
    with open("INPUT/bc.txt", "w") as f:
        f.write(sorted_bc("".join(out)))


def twist_expected():
    """Largest min(|T|/16, |T_n|/2) of the wall faces over the distance of their cell, with the distance MOSE takes:
    the cell centre to the four triangles between the face centre and its sides, minimum over the wall faces."""
    ni, nj, nk = TWIST
    xyz = twist_nodes()

    def sub(u, v):
        return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

    def dot(u, v):
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    def cross(u, v):
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]

    faces = {}
    for k in range(nk):
        for i in range(ni):
            faces[(i + 1, k + 1)] = [xyz[(i, 0, k)], xyz[(i + 1, 0, k)], xyz[(i + 1, 0, k + 1)], xyz[(i, 0, k + 1)]]

    def fan(p, c):
        m = [sum(q[n] for q in c) / 4.0 for n in range(3)]
        return min(triangle_distance(p, m, c[n], c[(n + 1) % 4]) for n in range(4))

    best = 0.0
    for (i, k), c in faces.items():
        t = [c[0][n] - c[1][n] + c[2][n] - c[3][n] for n in range(3)]
        bound = math.sqrt(dot(t, t)) / 16.0
        nrm = cross(sub(c[2], c[0]), sub(c[3], c[1]))
        if all(dot(cross(sub(c[(n + 1) % 4], c[n]), sub(c[(n + 2) % 4], c[(n + 1) % 4])), nrm) > 0.0 for n in range(4)):
            bound = min(bound, abs(dot(t, nrm)) / (2.0 * math.sqrt(dot(nrm, nrm))))
        corners = [xyz[(i - 1 + a, b, k - 1 + cc)] for a in (0, 1) for b in (0, 1) for cc in (0, 1)]
        p = [sum(q[n] for q in corners) / 8.0 for n in range(3)]
        d = min(fan(p, f) for f in faces.values())
        best = max(best, bound / d)
    return best


def twist_check(log):
    want = twist_expected()
    got = None
    for line in open(log, errors="replace"):
        m = re.search(r"twisted wall faces, error bound up to\s+([0-9.Ee+-]+)", line)
        if m:
            got = float(m.group(1))
    if got is None:
        print(f"    no report of the twist in {log} (expected {want:.3e})")
        sys.exit(1)
    print(f"    reported {got:.3e}, expected {want:.3e}")
    sys.exit(0 if abs(got / want - 1.0) <= 1.0e-3 else 1)


def fv1(chi):
    return chi ** 3 / (chi ** 3 + CV1 ** 3)


def factor():
    """F = cb1 fv2/kappa^2 - cw1 fw: the right-hand side of rho nu~ at rest is rho nu~^2 F / d^2."""
    fv2 = 1.0 - CHI / (1.0 + CHI * fv1(CHI))
    r = min(1.0 / fv2, 10.0)
    g = r + CW2 * (r ** 6 - r)
    fw = g * ((1.0 + CW3 ** 6) / (g ** 6 + CW3 ** 6)) ** (1.0 / 6.0)
    return CB1 * fv2 / KAPPA ** 2 - CW1 * fw


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the cells per side."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, ni - 1, nj - 1


def measured(name, path):
    """Per cell (i, j): d read from the solver, reference, nearest point at a node, old method."""
    data, ni, nj = read_zone(path)
    assert (ni, nj) == SIZES[name], (ni, nj)
    geo = geometry(name)
    f = factor()
    out = {}
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            n = (j - 1) * ni + (i - 1)
            rhs = data["rhoNuTilde"][n] / data["dt"][n]
            d = math.sqrt(RHO * NUT * NUT * f / rhs) if rhs < 0.0 else float("nan")
            out[(i, j)] = (d,) + geo[(i, j)]
    return out


def error(name, path):
    m = measured(name, path)
    cells = [(i, j) for (i, j) in m if j > SKIP]
    err = max(abs(m[c][0] / m[c][1] - 1.0) for c in cells)
    at_node = sum(1 for c in cells if m[c][2])
    print(f"{err:.3e} {len(cells)} {at_node}")


def check(pairs):
    ok = True
    for pair in pairs:
        name, err = pair.split("=")
        ok = ok and float(err) <= TOLS.get(name, TOL)
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case" and sys.argv[2] == "twist":
        write_twist_case()
    elif what == "case":
        (write_periodic_case if sys.argv[2] in PERIODIC else write_case)(sys.argv[2])
    elif what == "error":
        (periodic_error if sys.argv[2] in PERIODIC else error)(sys.argv[2], sys.argv[3])
    elif what == "check":
        check(sys.argv[2:])
    elif what == "twist":
        twist_check(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
