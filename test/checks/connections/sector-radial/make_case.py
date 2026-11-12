"""Inputs of the SectorRadial fast test, and the check of its output.

A 3-D sector of an annulus closed in k by the rotational periodicity (type 200 on faces 5 and 6): the cells next to face 5
go on in the cells next to face 6 turned back by the sector angle, and those next to face 6 in the cells next to face 5
turned forward. The sector spans theta from -45 to 45 degrees, symmetric about the plane z = 0, and the flow is radial,
u_r = U sin(pi (r - r_in)/(r_out - r_in)), with uniform density and pressure: a flow that is its own mirror image across
z = 0 (v even, w odd in z). The Euler equations keep that symmetry, and so must the scheme. The angle of 90 degrees makes
the two mirror planes of the problem, z = 0 and the periodic face at 45 degrees, commute with a limiter that acts on each
Cartesian component (the first changes the sign of w, the second swaps v and w); a turn by another angle mixes v and w and
would not be reproduced exactly. After 20 steps of RK3 the cells of column k and of column NK + 1 - k must hold the same
rho, u, v and p, and opposite w, within 1e-10 of their scales.

Case: x in [0, 0.2] (4 cells), r in [0.5, 1] (6 cells), theta in [-45, 45] degrees (6 cells), slip walls (300) on the
i and j faces; U = 50 m/s, rho = 1.161 kg/m3, p = 1e5 Pa.

Usage (from the test directory):
  python3 -B make_case.py case          INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py check FILE    prints the largest departure from the mirror symmetry in FILE (a field.tec)
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
NI, NJ, NK = 4, 6, 6
X1, R0, R1 = 0.2, 0.5, 1.0
TH = math.radians(90.0)
U = 50.0
P0, T0 = 1.0e5, 300.0
RHO0 = P0 / (RGAS * T0)


def fmt(v):
    return f"{v:.16E}"


def nodes():
    """Node coordinates (x, y, z) by (i, j, k), as MOSE reads them; the planes k and NK - k are mirror images."""
    out = {}
    for k in range(NK + 1):
        th = TH * (k - NK // 2) / NK                   # k - NK/2 an integer: the planes k and NK - k exactly mirrored
        for j in range(NJ + 1):
            r = R0 + (R1 - R0) * j / NJ
            for i in range(NI + 1):
                out[(i, j, k)] = tuple(float(fmt(v)) for v in (X1 * i / NI, r * math.cos(th), r * math.sin(th)))
    return out


def write_case():
    xyz = nodes()
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    cells = []
    for k in range(1, NK + 1):
        for j in range(1, NJ + 1):
            for i in range(1, NI + 1):
                c = [sum(xyz[(i - 1 + a, j - 1 + b, k - 1 + d)][n] for a in (0, 1) for b in (0, 1) for d in (0, 1)) / 8.0
                     for n in range(3)]
                r = math.hypot(c[1], c[2])
                ur = U * math.sin(math.pi * (r - R0) / (R1 - R0))
                cells.append((RHO0, 0.0, ur * c[1] / r, ur * c[2] / r, P0))
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={NJ + 1}, K={NK + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for n in range(3):
            for k in range(NK + 1):
                for j in range(NJ + 1):
                    for i in range(NI + 1):
                        f.write(f" {fmt(xyz[(i, j, k)][n])}\n")
        for v in range(5):
            for c in cells:
                f.write(f" {fmt(c[v])}\n")

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    out = []                                            # face by face, as ATLAS writes them
    for face, i in ((1, 1), (2, NI)):
        out += [rec(i, j, k, face, 300) for k in range(1, NK + 1) for j in range(1, NJ + 1)]
    for face, j in ((3, 1), (4, NJ)):
        out += [rec(i, j, k, face, 300) for k in range(1, NK + 1) for i in range(1, NI + 1)]
    for face, k in ((5, 1), (6, NK)):
        out += [rec(i, j, k, face, 200) for j in range(1, NJ + 1) for i in range(1, NI + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, names


def check(path):
    data, names = read_zone(path)
    rho = data[names[3]]
    worst, where = 0.0, ""
    for k in range(NK):
        km = NK - 1 - k                                # the mirror column
        for j in range(NJ):
            for i in range(NI):
                n, m = (k * NJ + j) * NI + i, (km * NJ + j) * NI + i
                diffs = {"rho": abs(rho[n] - rho[m]) / RHO0, "u": abs(data["u"][n] - data["u"][m]) / U,
                         "v": abs(data["v"][n] - data["v"][m]) / U, "w": abs(data["w"][n] + data["w"][m]) / U,
                         "p": abs(data["p"][n] - data["p"][m]) / P0}
                for q, d in diffs.items():
                    if d > worst:
                        worst, where = d, f"{q} at i {i + 1}, j {j + 1}, k {k + 1}"
    print(f"{worst:.3e} {where}")


def main():
    what = sys.argv[1]
    if what == "case":
        write_case()
    elif what == "check":
        check(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
