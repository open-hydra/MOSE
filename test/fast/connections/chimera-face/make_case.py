"""Inputs and checks of the ChimeraFace fast test.

The viscous flux of a chimera face (102) is the mean of two fluxes: one with the central differences along the face of
the row of the cell, one with those of the row of the ghost cells, filled from the donors. Beyond the edges of a block
face the row of the ghost cells has no cell; there the stencil takes the row extrapolated linearly, and in the k
direction of a 2-D mesh the ghost cell turned by the wedge angle, as the k neighbours of a cell are. Each case is one
evaluation of the residual: one step of a steady run with the forward Euler scheme, read back from the diagnostic output,
whose columns are the changes of the conservative variables over the step, -dt R / V, and dt; R / V is the residual per
unit volume.

Two blocks of cells of side H overlap by OVER cells in x: block 1 covers x from 0 to NI H, block 2 from (NI - OVER) H to
(2 NI - OVER) H, both on the same rows. Face 2 of block 1 and face 1 of block 2 are chimera faces whose two layers of ghost
cells coincide with cells of the other block (one donor each, volume fraction 1). The other faces extrapolate (400), which
is exact for a linear field and carries no viscous flux. Air of test/common/Air-viscous at rest in the mean, p and T
uniform:
  plane     2-D, NI x NJ square cells, one cell thick in z (null k faces): u = G y. The flux of y momentum through every
            x face is mu G H^2, and the y momentum residual is zero, but in the cells at the ends of the chimera faces
            (j = 1 and j = NJ) the chimera face takes 3/4 of it if the ghost cell is copied beyond the edge, and the
            residual of those cells is 1/4 of mu G / H per unit volume.
  plane3d   3-D, NI x NJ x NK cubic cells: u = G (y + z); the same for y momentum at j = 1 and NJ and for z momentum at
            k = 1 and NK.
  wedge     2-D axisymmetric wedge of 1 degree (type 200 on the k faces) whose first node row is at the radius R0: a
            radial velocity proportional to the radius, v = G y and w = G z at the cell centres. On the side of the cell
            the k neighbours are the cell turned by the wedge angle; with the ghost cell itself as the k neighbours of the
            ghost row, the chimera face takes half of the v/r of the divergence.
  wedgeref  the same wedge as one block of 2 NI - OVER cells in x, where the chimera faces are interior faces: the
            momentum residual of every cell of the two blocks of wedge must be that of the same cell here.
  fanchim   plane with the y spacing of block 1 growing 4 times across its last cell (H/4 on the chimera face 2,
            H from the node plane before it): the first layer of ghost cells extrapolated behind that face is turned
            over, and the viscous flux of a chimera face takes its metric: the run must stop at set-up with
            "ghost cell behind face 2: degenerate or inverted" and a non-zero exit code (a control). The centres of
            those ghost cells fall at negative y, outside block 2: harmless, as the donors are given by index and the
            run stops at set-up, before any ghost cell is filled.
The scale S = mu G / H of the residuals per unit volume is the x momentum residual of the cell (2, 1) of block 1 in plane,
whose one viscous x momentum flux is mu G H^2 through its upper face.

Usage (from the test directory):
  python3 -B make_case.py case NAME     INPUT/ic.tec and INPUT/bc.txt of plane, plane3d, wedge, wedgeref or fanchim
  python3 -B make_case.py check         the checks, on OUTPUT/diagnostic-NAME.tec of the four cases; exit code 1 on a failure
  python3 -B make_case.py stopped LOG   exit code 0 when LOG (fanchim) stopped at set-up on the ghost cell behind face 2
"""
import math
import os
import re
import sys

H = 1.0e-6                       # side of the cells, m
NI, NJ, NK, OVER = 8, 6, 4, 4    # cells of a block in x, y (and z in plane3d); overlap of the two blocks in x
G = 1.0e6                        # velocity gradient, 1/s
T0, P0 = 300.0, 1.0e4
RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
R0 = 20 * H                      # radius of the first node row of the wedge
WEDGE = math.radians(1.0)
TOL = 1.0e-9                     # largest residual of a consistent face, over S


def nodes(case, b):
    """Node coordinates [k][j][i] of block b."""
    if case == "wedgeref":
        ni, x0 = 2 * NI - OVER, 0.0
    else:
        ni, x0 = NI, (0.0, (NI - OVER) * H)[b - 1]
    nk = NK if case == "plane3d" else 1
    out = []
    for k in range(nk + 1):
        out.append([])
        for j in range(NJ + 1):
            row = []
            for i in range(ni + 1):
                if case.startswith("wedge"):
                    r, th = R0 + H * j, (k - 0.5) * WEDGE
                    row.append((x0 + H * i, r * math.cos(th), r * math.sin(th)))
                elif case == "fanchim" and b == 1 and i == ni:
                    row.append((x0 + H * i, 0.25 * H * j, H * k))     # the chimera face 2: a fan across the last cell
                else:
                    row.append((x0 + H * i, H * j, H * k))
            out[k].append(row)
    return out


def centre(nd, i, j, k):
    pts = [nd[k + c][j + b][i + a] for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    return tuple(sum(p[m] for p in pts) / 8.0 for m in range(3))


def state(case, x, y, z):
    rho = P0 / (RGAS * T0)
    if case in ("plane", "fanchim"):
        return rho, G * y, 0.0, 0.0, P0
    if case == "plane3d":
        return rho, G * (y + z), 0.0, 0.0, P0
    return rho, 0.0, G * y, G * z, P0


def rec(b, i, j, k, f, t):
    return f"{b:8d}{i:8d}{j:8d}{k:8d}{f:8d}{t:8d}\n"


def write_case(case):
    nblk = 1 if case == "wedgeref" else 2
    blocks = [nodes(case, b) for b in range(1, nblk + 1)]
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b, nd in enumerate(blocks, 1):
            nk, nj, ni = len(nd) - 1, len(nd[0]) - 1, len(nd[0][0]) - 1
            f.write(f" ZONE  T = B{b}, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                    "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for axis in range(3):
                f.write("".join(f" {nd[k][j][i][axis]:.17E}\n"
                                for k in range(nk + 1) for j in range(nj + 1) for i in range(ni + 1)))
            cells = [state(case, *centre(nd, i, j, k)) for k in range(nk) for j in range(nj) for i in range(ni)]
            for m in range(5):
                f.write("".join(f" {c[m]:.17E}\n" for c in cells))
    kfaces = {"plane": 0, "fanchim": 0, "plane3d": 400}.get(case, 200)
    rows = []
    for b, nd in enumerate(blocks, 1):
        nk, nj, ni = len(nd) - 1, len(nd[0]) - 1, len(nd[0][0]) - 1
        for face in range(1, 7):
            if face <= 2:
                cells = [(1 if face == 1 else ni, j, k) for k in range(1, nk + 1) for j in range(1, nj + 1)]
            elif face <= 4:
                cells = [(i, 1 if face == 3 else nj, k) for k in range(1, nk + 1) for i in range(1, ni + 1)]
            else:
                cells = [(i, j, 1 if face == 5 else nk) for j in range(1, nj + 1) for i in range(1, ni + 1)]
            chim = nblk == 2 and (b, face) in ((1, 2), (2, 1))
            for (i, j, k) in cells:
                if not chim:
                    rows.append(rec(b, i, j, k, face, kfaces if face >= 5 else 400))
                    continue
                # ghost layers 1 and 2 of the face, as cells of the other block
                if b == 1:
                    donors = [(2, NI + g - (NI - OVER), j, k) for g in (1, 2)]
                else:
                    donors = [(1, 1 - g + (NI - OVER), j, k) for g in (1, 2)]
                rows.append(rec(b, i, j, k, face, 102))
                rows.append(f"{1:8d}{1:8d}\n")
                rows.extend(f"{d[0]:8d}{d[1]:8d}{d[2]:8d}{d[3]:8d}    0.1000000000E+01\n" for d in donors)
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(rows))


def read_zones(path):
    """Cell values of every zone of a MOSE ASCII file (BLOCK packing): a list of (values by name, ni, nj, nk)."""
    with open(path) as f:
        text = f.read()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', text.splitlines()[0].split('=', 1)[1])]
    zones = []
    for part in re.split(r'(?m)^\s*ZONE', text)[1:]:
        head, _, body = part.partition('\n')
        ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', head).group(1)) for k in ('I', 'J', 'K'))
        values = [float(v) for v in " ".join(l for l in body.splitlines() if not l.startswith("DATASETAUXDATA")).split()]
        nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
        for m, name in enumerate(names):
            cnt = nn if m < 3 else nc
            out[name] = values[pos:pos + cnt]
            pos += cnt
        zones.append((out, ni - 1, nj - 1, max(nk - 1, 1)))
    return zones


def rates(path):
    """Per zone: {(i, j, k): residual per unit volume R / V of rho u, rho v, rho w, rho E}, cells from 1."""
    out = []
    for z, ni, nj, nk in read_zones(path):
        d = {}
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    n = i + ni * (j + nj * k)
                    d[(i + 1, j + 1, k + 1)] = [-z[v][n] / z["dt"][n] for v in ("rhou", "rhov", "rhow", "rhoe")]
        out.append(d)
    return out


def check():
    ok = True
    plane = rates("OUTPUT/diagnostic-plane.tec")
    scale = abs(plane[0][(2, 1, 1)][0])
    print(f"  scale S = mu G / H from the x momentum residual of the cell (2, 1) of block 1 in plane: {scale:.6e}")
    for case in ("plane", "plane3d"):
        z = plane if case == "plane" else rates("OUTPUT/diagnostic-plane3d.tec")
        s = scale if case == "plane" else abs(z[0][(2, 1, 2)][0])
        nk = NK if case == "plane3d" else 1
        edge, rest = 0.0, 0.0
        for b, d in enumerate(z, 1):
            face_i = NI if b == 1 else 1
            for (i, j, k), r in d.items():
                if (b == 1 and i == 1) or (b == 2 and i == NI):
                    continue          # next to an outer x face (400), which carries no viscous flux
                for comp, at_edge in ((1, j in (1, NJ)), (2, nk > 1 and k in (1, nk))):
                    v = abs(r[comp]) / s
                    if i == face_i and at_edge:
                        edge = max(edge, v)
                    else:
                        rest = max(rest, v)
        good = edge <= TOL and rest <= TOL
        ok = ok and good
        what = "y momentum at j = 1, NJ" + (", z momentum at k = 1, NK" if nk > 1 else "")
        print(f"  [{'ok  ' if good else 'FAIL'}] {case}: {what} of the cells at the ends of the chimera faces, over S: "
              f"{edge:.3e} (1/4 with the ghost cell copied beyond the edge); elsewhere {rest:.3e}; both <= {TOL:g}")
    wedge = rates("OUTPUT/diagnostic-wedge.tec")
    ref = rates("OUTPUT/diagnostic-wedgeref.tec")[0]
    chim, edge, rest = 0.0, 0.0, 0.0
    for b, d in enumerate(wedge, 1):
        off, face_i = (0, NI) if b == 1 else (NI - OVER, 1)
        for (i, j, k), r in d.items():
            v = max(abs(a - c) for a, c in zip(r[:3], ref[(i + off, j, k)][:3])) / scale
            if i == face_i:
                if j in (1, NJ):
                    edge = max(edge, v)
                else:
                    chim = max(chim, v)
            else:
                rest = max(rest, v)
    good = max(chim, edge, rest) <= TOL
    ok = ok and good
    print(f"  [{'ok  ' if good else 'FAIL'}] wedge against wedgeref, largest difference of the momentum residual over S: "
          f"cells of the chimera faces {chim:.3e}, the same at j = 1, NJ {edge:.3e}, other cells {rest:.3e}; "
          f"all <= {TOL:g}")
    return ok


def stopped(log):
    """fanchim: the run stopped at set-up on a ghost cell of the first layer behind the chimera face 2."""
    with open(log, errors="replace") as f:
        text = f.read()
    hit = re.search(r"block 1, cell +\d+ +\d+ +\d+, ghost cell behind face 2: degenerate or inverted", text) is not None
    print(f"  [{'ok  ' if hit else 'FAIL'}] fanchim: stop at set-up on the turned-over ghost cell behind the chimera face 2")
    return hit


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "case":
        write_case(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == "stopped":
        sys.exit(0 if stopped(sys.argv[2]) else 1)
    elif len(sys.argv) == 2 and sys.argv[1] == "check":
        sys.exit(0 if check() else 1)
    else:
        sys.exit(__doc__)
