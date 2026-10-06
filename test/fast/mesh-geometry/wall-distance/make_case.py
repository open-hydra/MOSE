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
The reference is the exact distance from the cell centre (mean of the nodes; in the wedge in the coordinates x and
r = sqrt(y^2 + z^2)) to the wall of the mesh (straight faces between the wall nodes), by a scan of all faces. It is checked
on the cells that do not touch the wall (in one step the wall reaches only the cells next to it, through the diffusion
across the wall face and the central difference of the cb2 term): the largest relative difference must be below 1e-9.

Usage (from the test directory):
  python3 -B make_case.py case NAME            INPUT/ic.tec and INPUT/bc.txt of corner, arc or lip
  python3 -B make_case.py error NAME FILE      prints the largest relative error of d read from FILE (a diagnostic.tec)
  python3 -B make_case.py check E1 E2 ...      exit status 0 when every error is below the tolerance
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
        f.write("".join(out))


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


def check(errs):
    ok = all(float(e) <= TOL for e in errs)
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "error":
        error(sys.argv[2], sys.argv[3])
    elif what == "check":
        check(sys.argv[2:])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
