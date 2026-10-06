"""Inputs of the MeshType fast test, and the checks of its outputs.

MOSE takes the type of the mesh from its nodes: 1-D, 2-D or 3-D from the number of cells in j and k, and on a 2-D mesh
an axisymmetric wedge when the k = 1 node plane is the k = 0 plane turned about the x axis (the angle delthe of the
wedge); on a 3-D mesh delthe is the angle of the sectors closed by the rotational periodicity (type 200).
The nodes are checked relative to the thickness of the block, so that their
round-off, with the digits they are written with, does not change the type: relative to the move of the nodes per cell in k
(the thickness of a 2-D block). Ten cases:
  slab      2-D plane slab, z from -0.5 to 0.5, y falling from 1 (j = 0) to 0.5 (j = NJ), wall (301) on the j = NJ face,
            symmetry (300) elsewhere. Its k = 1 plane is a translation of the k = 0 plane, not a rotation, so the mesh
            is plane. The wall distance, read back from the Spalart-Allmaras source (gas at rest, nu~ uniform, one
            explicit Euler step, as in WallDistance), must be the plane distance y_c - 0.5 within 1e-9; on a wedge it
            would be measured from the x axis.
  flipwedge 2-D wedge of 1 degree, right-handed with the axis on the j = NJ face (r = 0 exactly there) and the wall on
            the j = 0 face (r = 1): a wedge, so the distance read back must be the meridional one, 1 - r_c, within 1e-9.
  wedge5    two blocks of a wedge of 5 degrees as flipwedge, that do not touch (x from 0 to 1.2 and from 2 to 3.2), with
            radii falling as (j/NJ)^1.5 and (j/NJ)^1.12, every coordinate written with 5 significant digits: a wedge,
            whose nodes match a turn only to the round-off of those digits, and whose blocks give angles that differ
            by it, 1.7e-6 rad (their probe nodes are at different radii). The distance read back must be the
            meridional distance to the wall segments of the nodes as written within 1e-9.
  farslab   2-D plane slab far from the x axis: y falling from 1.001 (j = 0) to 1 (j = NJ, the wall), z from 0 to 1e-4,
            17 digits: plane (a turn about the x axis by 1e-4 rad nearly fits its nodes, a translation fits them
            exactly), so the distance must be the plane one, y_c - 1, within 1e-9.
  coarsewedge  the block of flipwedge with its k = 0 plane at z = 0 and its k = 1 plane turned by -1 degree, every
            coordinate written with 2 significant digits: its nodes are nearer a turn than a translation, but match it
            only to 3 % of the thickness of the block, so the run must stop at set-up with a non-zero exit code and a
            message that asks for more digits.
  tiltslab  the plane slab of slab with its k = 1 plane tilted, z from 0.5 at j = 0 to 0.52 at j = NJ: its k planes are
            not parallel, nearer a translation than a turn but only to 2e-2 of its thickness, so the run must stop at
            set-up with a non-zero exit code and a message on the mesh type.
  sectorshift100, sectorshift20  a 3-D sector of 90 degrees about x, 2 x 2 x K cells (x from 0 to 0.2, r from 0.5 to 1), type
            200 on its k faces, whose last k node plane does not match the first one turned: its nodes are turned by a
            fraction f of a k cell more, but the nodes (0,0) and (0,1); K = 100 and f = 0.5, K = 20 and f = 0.1. The nodes
            match the turn to a few thousandths of the thickness of the block, but to half a cell, or a tenth, of the cells in
            k: the run must stop at set-up with a non-zero exit code and a message on the mesh type.
  mixed     two blocks that do not touch, a wedge of 1 degree and a plane slab: a mesh that is neither a wedge nor
            plane, so the run must stop at set-up with a non-zero exit code and a message on the mesh type.
  slab3d    3-D slab, 6 x 6 x 4 cells, y falling with j, type 200 on the k faces (periodicity in k with no rotation, the
            k = NK plane being a translation of the k = 0 plane) and extrapolation (400) on the other faces; uniform
            flow (100, 50, 0) m/s. Five explicit steps must leave the flow uniform within 1e-12 of |u|; a rotation of
            the periodic images would turn v into w.
Gas: test/common/Air-viscous (slab, flipwedge, at 300 K and 1 kPa with nu~ = nu) and test/common/Air (mixed, slab3d, at
1 bar and 300 K).

Usage (from the test directory):
  python3 -B make_case.py case NAME          INPUT/ic.tec and INPUT/bc.txt of slab, flipwedge, wedge5, farslab,
                                             coarsewedge, tiltslab, sectorshift100, sectorshift20, mixed or slab3d
  python3 -B make_case.py error NAME FILE    largest relative error of the distance read from FILE (a diagnostic.tec)
  python3 -B make_case.py uniform FILE       largest change of u, v, w over |u| in FILE (a field.tec), and of rho and p
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P_SA, P_EU = 300.0, 1.0e3, 1.0e5
RHO_SA = P_SA / (RGAS * T0)
RHO_EU = P_EU / (RGAS * T0)
CHI = 1.0
NUT = CHI * MU / RHO_SA
CB1, CB2, SIGMA, KAPPA = 0.1355, 0.622, 2.0 / 3.0, 0.41
CW2, CW3, CV1 = 0.3, 2.0, 7.1
CW1 = CB1 / KAPPA ** 2 + (1.0 + CB2) / SIGMA
DELTA = math.radians(1.0)
DELTA5 = math.radians(5.0)
DIGITS = {"wedge5": 5, "coarsewedge": 2}   # significant digits of the coordinates in ic.tec; 17 for the other cases
U0 = (100.0, 50.0, 0.0)
TOL_D = 1.0e-9
TOL_U = 1.0e-12


def fmt(v, digits=17):
    return f"{v:.{digits - 1}E}"


def as_read(v, digits=17):
    """The value as MOSE reads it from the text of ic.tec."""
    return float(fmt(v, digits))


def slab_nodes():
    ni, nj = 12, 8
    xs = [i / ni for i in range(ni + 1)]
    ys = [0.5 + 0.5 * (1.0 - j / nj) ** 1.3 for j in range(nj + 1)]          # 1 at j = 0, 0.5 at j = NJ
    return [[[(x, y, (-0.5, 0.5)[k]) for x in xs] for y in ys] for k in range(2)]


def flipwedge_nodes():
    ni, nj = 12, 10
    xs = [0.1 * i for i in range(ni + 1)]
    rs = [1.0 - (j / nj) ** 1.5 for j in range(nj)] + [0.0]                     # 1 at j = 0, exactly 0 at j = NJ
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    zsign = (1.0, -1.0)                                                         # right-handed with the axis at j = NJ
    return [[[(x, r * c, zsign[k] * r * s) for x in xs] for r in rs] for k in range(2)]


def wedge5_blocks():
    out = []
    for x0, p in ((0.0, 1.5), (2.0, 1.12)):
        ni, nj = 12, 10
        xs = [x0 + 0.1 * i for i in range(ni + 1)]
        rs = [1.0 - (j / nj) ** p for j in range(nj)]                           # 1 at j = 0
        rs += [1.0e-6 * rs[-1]]                                                 # axis at j = NJ
        c, s = math.cos(0.5 * DELTA5), math.sin(0.5 * DELTA5)
        zsign = (1.0, -1.0)
        out.append([[[(x, r * c, zsign[k] * r * s) for x in xs] for r in rs] for k in range(2)])
    return out


def coarsewedge_nodes():
    ni, nj = 12, 10
    xs = [0.1 * i for i in range(ni + 1)]
    rs = [1.0 - (j / nj) ** 1.5 for j in range(nj)]
    rs += [1.0e-6 * rs[-1]]
    turn = (0.0, -DELTA)                                                        # right-handed, as flipwedge
    return [[[(x, r * math.cos(turn[k]), r * math.sin(turn[k])) for x in xs] for r in rs] for k in range(2)]


def sector_shift_nodes(nk, frac, ni=2, nj=2, span=0.5 * math.pi):
    """Nodes [k][j][i] of a sector of span about x, nk cells in k; the nodes of the last k plane turned by frac of a k cell
    more, but (0,0) and (0,1), the node the angle of the block is taken from."""
    nodes = []
    for k in range(nk + 1):
        plane = []
        for j in range(nj + 1):
            r = 0.5 + 0.5 * j / nj
            row = []
            for i in range(ni + 1):
                phi = span * k / nk
                if k == nk and (i, j) not in ((0, 0), (0, 1)):
                    phi += frac * span / nk
                row.append((0.1 * i, r * math.cos(phi), r * math.sin(phi)))
            plane.append(row)
        nodes.append(plane)
    return nodes


def farslab_nodes():
    ni, nj = 12, 8
    xs = [i / ni for i in range(ni + 1)]
    ys = [1.0 + 0.001 * (1.0 - j / nj) for j in range(nj + 1)]                  # 1.001 at j = 0, 1 at j = NJ
    return [[[(x, y, (0.0, 1.0e-4)[k]) for x in xs] for y in ys] for k in range(2)]


def mixed_nodes():
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    b1 = [[[(0.1 * i, (0.5 + 0.125 * j) * c, (-1.0, 1.0)[k] * (0.5 + 0.125 * j) * s) for i in range(5)]
           for j in range(5)] for k in range(2)]
    b2 = [[[(2.0 + 0.1 * i, 0.5 + 0.125 * j, (0.0, 0.1)[k]) for i in range(5)] for j in range(5)] for k in range(2)]
    return [b1, b2]


def tiltslab_nodes():
    nodes = slab_nodes()
    nj = len(nodes[1]) - 1
    nodes[1] = [[(x, y, 0.5 + 0.02 * j / nj) for (x, y, z) in row] for j, row in enumerate(nodes[1])]
    return nodes


def slab3d_nodes():
    ni, nj, nk = 6, 6, 4
    return [[[(0.1 * i, 1.0 - 0.5 * j / nj, 0.05 * k) for i in range(ni + 1)] for j in range(nj + 1)]
            for k in range(nk + 1)]


def write_ic(blocks, cellvalues, names, digits=17):
    """blocks: list of node arrays [k][j][i] of (x, y, z); cellvalues: per block, list of columns; the coordinates with
    the given number of significant digits."""
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  ' + " ".join(f'"{n}"' for n in names) + "\n")
        nv = len(names)
        for b, nodes in enumerate(blocks):
            nk, nj, ni = len(nodes), len(nodes[0]), len(nodes[0][0])
            f.write(f" ZONE  T = B{b + 1}-IG, I={ni}, J={nj}, K={nk}, DATAPACKING=BLOCK, "
                    f"VARLOCATION=([1-3]=NODAL,[4-{3 + nv}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for comp in range(3):
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            f.write(f" {fmt(nodes[k][j][i][comp], digits)}\n")
            ncell = (ni - 1) * (nj - 1) * max(nk - 1, 1)
            for v in cellvalues[b]:
                f.write(f" {fmt(v)}\n" * ncell)


def rec(b, i, j, k, face, t):
    return '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, face, t)


WALL_LINE = "    0.000000E+00,    0.000000E+00,    0.000000E+00,\n"


def faces_2d(b, ni, nj, types):
    """Records of a 2-D block; types = (face1, face2, face3, face4, face56); a wall (301) gets its property line."""
    out = []
    def add(i, j, face, t):
        out.append(rec(b, i, j, 1, face, t))
        if t == 301:
            out.append(WALL_LINE)
    for j in range(1, nj + 1):
        add(1, j, 1, types[0])
        add(ni, j, 2, types[1])
    for i in range(1, ni + 1):
        add(i, 1, 3, types[2])
        add(i, nj, 4, types[3])
    for face in (5, 6):
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                add(i, j, face, types[4])
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
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    sa = [RHO_SA, 0.0, 0.0, 0.0, P_SA, RHO_SA * NUT]
    eu = [RHO_EU, 0.0, 0.0, 0.0, P_EU]
    if name == "slab":
        write_ic([slab_nodes()], [sa], ["rho1", "u", "v", "w", "p", "rhoNuTilde"])
        recs = faces_2d(1, 12, 8, (300, 300, 300, 301, 300))
    elif name == "flipwedge":
        write_ic([flipwedge_nodes()], [sa], ["rho1", "u", "v", "w", "p", "rhoNuTilde"])
        recs = faces_2d(1, 12, 10, (300, 300, 301, 300, 200))
    elif name == "wedge5":
        write_ic(wedge5_blocks(), [sa, sa], ["rho1", "u", "v", "w", "p", "rhoNuTilde"], DIGITS[name])
        recs = faces_2d(1, 12, 10, (300, 300, 301, 300, 200)) + faces_2d(2, 12, 10, (300, 300, 301, 300, 200))
    elif name == "farslab":
        write_ic([farslab_nodes()], [sa], ["rho1", "u", "v", "w", "p", "rhoNuTilde"])
        recs = faces_2d(1, 12, 8, (300, 300, 300, 301, 300))
    elif name == "tiltslab":
        write_ic([tiltslab_nodes()], [eu], ["rho1", "u", "v", "w", "p"])
        recs = faces_2d(1, 12, 8, (300, 300, 300, 300, 300))
    elif name == "coarsewedge":
        write_ic([coarsewedge_nodes()], [eu], ["rho1", "u", "v", "w", "p"], DIGITS[name])
        recs = faces_2d(1, 12, 10, (300, 300, 300, 300, 200))
    elif name in ("sectorshift100", "sectorshift20"):
        nk, frac = {"sectorshift100": (100, 0.5), "sectorshift20": (20, 0.1)}[name]
        write_ic([sector_shift_nodes(nk, frac)], [eu], ["rho1", "u", "v", "w", "p"])
        recs = []
        for k in range(1, nk + 1):
            for j in range(1, 3):
                recs += [rec(1, 1, j, k, 1, 300), rec(1, 2, j, k, 2, 300)]
            for i in range(1, 3):
                recs += [rec(1, i, 1, k, 3, 300), rec(1, i, 2, k, 4, 300)]
        for j in range(1, 3):
            for i in range(1, 3):
                recs += [rec(1, i, j, 1, 5, 200), rec(1, i, j, nk, 6, 200)]
    elif name == "mixed":
        write_ic(mixed_nodes(), [eu, eu], ["rho1", "u", "v", "w", "p"])
        recs = faces_2d(1, 4, 4, (300, 300, 300, 300, 200)) + faces_2d(2, 4, 4, (300, 300, 300, 300, 300))
    elif name == "slab3d":
        write_ic([slab3d_nodes()], [[RHO_EU, U0[0], U0[1], U0[2], P_EU]], ["rho1", "u", "v", "w", "p"])
        ni, nj, nk = 6, 6, 4
        recs = []
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                recs += [rec(1, 1, j, k, 1, 400), rec(1, ni, j, k, 2, 400)]
            for i in range(1, ni + 1):
                recs += [rec(1, i, 1, k, 3, 400), rec(1, i, nj, k, 4, 400)]
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                recs += [rec(1, i, j, 1, 5, 200), rec(1, i, j, nk, 6, 200)]
    else:
        sys.exit(f"unknown case {name}")
    with open("INPUT/bc.txt", "w") as f:
        f.write(sorted_bc("".join(recs)))


def read_zones(path):
    """Cell values of every zone of a MOSE ASCII file (BLOCK packing), names quoted or not, with the cells per
    direction: a list of (values, ni, nj)."""
    with open(path) as f:
        text = f.read()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', text.splitlines()[0].split('=', 1)[1])]
    zones = []
    for part in re.split(r'(?m)^\s*ZONE', text)[1:]:
        head, _, body = part.partition('\n')
        ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', head).group(1)) for k in ('I', 'J', 'K'))
        values = [float(v) for v in " ".join(body.splitlines()).split()]
        nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
        for k, name in enumerate(names):
            cnt = nn if k < 3 else nc
            out[name] = values[pos:pos + cnt]
            pos += cnt
        zones.append((out, ni - 1, nj - 1))
    return zones


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the cells per direction."""
    return read_zones(path)[0]


def factor():
    fv1 = CHI ** 3 / (CHI ** 3 + CV1 ** 3)
    fv2 = 1.0 - CHI / (1.0 + CHI * fv1)
    r = min(1.0 / fv2, 10.0)
    g = r + CW2 * (r ** 6 - r)
    fw = g * ((1.0 + CW3 ** 6) / (g ** 6 + CW3 ** 6)) ** (1.0 / 6.0)
    return CB1 * fv2 / KAPPA ** 2 - CW1 * fw


def segment_distance(p, a, b):
    """Distance from the point p to the segment from a to b, in the plane (x, r)."""
    dx, dr = b[0] - a[0], b[1] - a[1]
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dr) / (dx * dx + dr * dr)))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dr)


def reference(name):
    """Exact distance of every cell (b, i, j), 1-based, that does not touch the wall, from the nodes as written."""
    blocks = {"slab": [slab_nodes()], "flipwedge": [flipwedge_nodes()], "wedge5": wedge5_blocks(),
              "farslab": [farslab_nodes()]}[name]
    digits = DIGITS.get(name, 17)
    read = [[[[tuple(as_read(v, digits) for v in q) for q in row] for row in plane] for plane in nodes] for nodes in blocks]
    # wedge5: the wall segments of the meridional plane, between the nodes of the j = 0 face (both k planes have the
    # same x and r there)
    walls = [((nodes[0][0][i - 1][0], math.hypot(*nodes[0][0][i - 1][1:])),
              (nodes[0][0][i][0], math.hypot(*nodes[0][0][i][1:])))
             for nodes in read for i in range(1, len(nodes[0][0]))]
    out = {}
    for b, nodes in enumerate(read, 1):
        nj, ni = len(nodes[0]) - 1, len(nodes[0][0]) - 1
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                corners = [nodes[k][jj][ii] for k in (0, 1) for jj in (j - 1, j) for ii in (i - 1, i)]
                if name == "slab":
                    if j == nj:
                        continue
                    out[(b, i, j)] = sum(q[1] for q in corners) / 8.0 - 0.5
                elif name == "farslab":
                    if j == nj:
                        continue
                    out[(b, i, j)] = sum(q[1] for q in corners) / 8.0 - 1.0
                elif name == "flipwedge":
                    if j == 1:
                        continue
                    out[(b, i, j)] = 1.0 - sum(math.hypot(q[1], q[2]) for q in corners) / 8.0
                else:
                    if j == 1:
                        continue
                    pc = (sum(q[0] for q in corners) / 8.0, sum(math.hypot(q[1], q[2]) for q in corners) / 8.0)
                    out[(b, i, j)] = min(segment_distance(pc, a, c) for a, c in walls)
    return out


def error(name, path):
    zones = read_zones(path)
    ref = reference(name)
    f = factor()
    err = 0.0
    for (b, i, j), d0 in ref.items():
        data, ni, nj = zones[b - 1]
        n = (j - 1) * ni + (i - 1)
        rhs = data["rhoNuTilde"][n] / data["dt"][n]
        d = math.sqrt(RHO_SA * NUT * NUT * f / rhs) if rhs < 0.0 else float("nan")
        err = max(err, abs(d / d0 - 1.0)) if d == d else float("inf")
    print(f"{err:.3e} {len(ref)}")


def uniform(path):
    data, ni, nj = read_zone(path)
    names = list(data)
    rho, p = data[names[3]], data["p"]
    du = max(max(abs(a - b) for a in data[n]) for n, b in zip(("u", "v", "w"), U0)) / math.sqrt(sum(x * x for x in U0))
    dr = max(abs(a / RHO_EU - 1.0) for a in rho)
    dp = max(abs(a / P_EU - 1.0) for a in p)
    print(f"{max(du, dr, dp):.3e}")


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "error":
        error(sys.argv[2], sys.argv[3])
    elif what == "uniform":
        uniform(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
