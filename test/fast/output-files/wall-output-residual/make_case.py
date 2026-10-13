"""Inputs of the WallOutputResidual fast test, and the checks of its outputs.

A short hypersonic plate segment, laminar: 6 x 10 cells, k faces null. Along x, cells of 1 mm: the first one over a symmetry
plane (300) ahead of the plate, the other five over an isothermal wall at 300 K (302). Along y, cells growing by 1.25 from
20 um at the wall (10 cells, 0.67 mm). Supersonic inflow (405) at Mach 5, 68.3 K, 4000 Pa; supersonic outflow (406,
extrapolation) on the outflow face and on the top face. Air of test/common/Air-viscous with the Sutherland viscosity of
test/common/Air-sutherland.

The solver converges the steady state; then OUTPUT/diagnostic.tec holds the residual of every cell, which at convergence is
round-off small everywhere, wall cells included. The wall output is written just before the diagnostic one: it must not
leave anything in the residual. On every plate face the residual of rho u in the wall cell must stay below 1e-6 of the wall
shear force tau_x A of that face, and the residual of rho e below 1e-6 of the wall heat flow q_w A (A from the nodes of
OUTPUT/field.tec).

Usage (from the test directory):
  python3 -B make_case.py case                        INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py converged WALL_1 WALL_2     exit status 0 when the wall shear changed by less than 1e-11
  python3 -B make_case.py residual DIAG WALL FIELD    exit status 0 when the wall-cell residuals pass
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MACH, T_IN, P_IN, T_WALL = 5.0, 68.3, 4000.0, 300.0
NI, NS, NJ, DX, Y1, GROWTH = 6, 1, 10, 1.0e-3, 2.0e-5, 1.25


def write_case():
    y = [0.0]
    for j in range(NJ):
        y.append(y[-1] + Y1 * GROWTH ** j)
    x = [(i - NS) * DX for i in range(NI + 1)]
    rho = P_IN / (RGAS * T_IN)
    u = MACH * math.sqrt(1.4 * RGAS * T_IN)
    nc = NI * NJ
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={NJ + 1}, K=2, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(NJ + 1):
                    for i in range(NI + 1):
                        f.write(f" {(x[i], y[j], k * DX)[c]:.16E}\n")
        for col in ([rho] * nc, [u] * nc, [0.0] * nc, [0.0] * nc, [P_IN] * nc):
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    inflow = f"    {MACH:.6E}, {T_IN:.6E}, {P_IN:.6E}, normal, normal, 1.000000E+00, 1.000000E+00,\n"
    extrapolation = "    0.000000E+00, 1.000000E+00,\n"
    out = []
    for j in range(1, NJ + 1):
        out += [rec(1, j, 1, 405), inflow]
    for j in range(1, NJ + 1):
        out += [rec(NI, j, 2, 406), extrapolation]
    for i in range(1, NI + 1):
        if i <= NS:
            out.append(rec(i, 1, 3, 300))
        else:
            out += [rec(i, 1, 3, 302), f"    {T_WALL:.6E}, 0.000000E+00, 0.000000E+00,\n"]
    for i in range(1, NI + 1):
        out += [rec(i, NJ, 4, 406), extrapolation]
    for face in (5, 6):
        out += [rec(i, j, face, 0) for j in range(1, NJ + 1) for i in range(1, NI + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Values of a one-zone MOSE ASCII file (BLOCK packing): nodal x, y, z first, then the cell variables; the numbers of
    nodes along i, j, k."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, (ni, nj, nk)


def converged(w1, w2):
    a, _ = read_zone(w1)
    b, _ = read_zone(w2)
    plate = [k for k, t in enumerate(a["tauX"]) if t != 0.0]
    d = max(abs(b["tauX"][k] - a["tauX"][k]) / abs(a["tauX"][k]) for k in plate)
    print(f"    wall shear change over the last 400 iterations: {d:.3e}")
    sys.exit(0 if d < 1.0e-11 else 1)


def residual(diag, wall, field):
    d, (ni, nj, _) = read_zone(diag)
    w, _ = read_zone(wall)
    f, (fni, fnj, _) = read_zone(field)
    ncx = ni - 1
    x, z = f["x"], f["z"]
    plate = [i for i, t in enumerate(w["tauX"]) if t != 0.0]
    worst_u = worst_e = 0.0
    for i in plate:
        area = abs((x[i + 1] - x[i]) * (z[fni * fnj] - z[0]))     # face (i, j = 1): dx of its nodes times the k extent
        ru = abs(d["rhou"][i]) / abs(w["tauX"][i] * area)          # wall cell (i, j = 1) is cell i of the first row
        re_ = abs(d["rhoe"][i]) / abs(w["qw"][i] * area)
        worst_u, worst_e = max(worst_u, ru), max(worst_e, re_)
    interior = max(abs(v) for k, v in enumerate(d["rhou"]) if k >= ncx)
    scale = max(abs(w["tauX"][i] * abs((x[i + 1] - x[i]) * (z[fni * fnj] - z[0]))) for i in plate)
    print(f"    wall cells: |r(rho u)| / |tau_x A| up to {worst_u:.3e}, |r(rho e)| / |q_w A| up to {worst_e:.3e} "
          f"({len(plate)} plate faces); other rows: |r(rho u)| up to {interior / scale:.3e} of the largest tau_x A")
    sys.exit(0 if worst_u < 1.0e-6 and worst_e < 1.0e-6 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case()
    elif what == "converged":
        converged(sys.argv[2], sys.argv[3])
    elif what == "residual":
        residual(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
