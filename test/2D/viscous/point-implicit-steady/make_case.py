"""Inputs of the PointImplicitSteady test, and the checks of its outputs.

A short hypersonic flat-plate segment: 6 x 16 cells, k faces null. Along x, cells of 1 mm: the first one over a symmetry plane
(300) ahead of the plate, the other five over an isothermal wall at 300 K (302). Along y, cells growing by 1.25 from 4 um at
the wall (16 cells, 0.55 mm). Supersonic inflow (405) at Mach 5, 68.3 K, 4000 Pa; supersonic outflow (406, extrapolation) on
the outflow face and on the top face. Air of test/common/Air-viscous with the Sutherland viscosity of test/common/Air-sutherland.
Free-stream turbulence:
  SA   rho nu~ = 3 mu (mu = 4.605e-6 Pa s at 68.3 K, from the Sutherland table),
  SST  rho k = 0.5, rho omega = 248.1 (mu_t about 90 mu), omega-wall-bc = practical;
the eddy viscosity then matters in the wall layer, so the wall shear and heat flux feel the turbulence equations there.

The test converges the steady state with the point-implicit factor at (cfl, vnn) = (0.8, 0.4), then continues from it at
(0.4, 0.2). A steady state that does not depend on the time step stays where it is; the wall shear and the wall heat flux after
the second run must equal those of the first to 1e-9 (largest relative difference over the plate faces).

Usage (from the test directory):
  python3 -B make_case.py case MODEL                    INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py converged WALL_1 WALL_2       exit status 0 when the wall shear changed by less than 1e-11
  python3 -B make_case.py compare WALL_A WALL_B         wall shear and heat flux of two runs; exit status 0 when equal to 1e-9
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MACH, T_IN, P_IN, T_WALL = 5.0, 68.3, 4000.0, 300.0
NI, NS, NJ, DX, Y1, GROWTH = 6, 1, 16, 1.0e-3, 4.0e-6, 1.25
MU_IN = 4.582749e-6 + 0.3 * (4.658099e-6 - 4.582749e-6)   # Sutherland table of test/common/Air-sutherland, 68 K and 69 K
RANS_IN = {"SA": [3.0 * MU_IN], "SST": [0.5, 248.1]}


def write_case(model):
    y = [0.0]
    for j in range(NJ):
        y.append(y[-1] + Y1 * GROWTH ** j)
    x = [(i - NS) * DX for i in range(NI + 1)]
    rho = P_IN / (RGAS * T_IN)
    u = MACH * math.sqrt(1.4 * RGAS * T_IN)
    rans = RANS_IN[model]
    names = ' "rhoNuTilde"' if model == "SA" else ' "rhoK" "rhoOmega"'
    nvar, nc = 8 + len(rans), NI * NJ
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"' + names + "\n")
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={NJ + 1}, K=2, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-{nvar}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(NJ + 1):
                    for i in range(NI + 1):
                        f.write(f" {(x[i], y[j], k * DX)[c]:.16E}\n")
        for col in [[rho] * nc, [u] * nc, [0.0] * nc, [0.0] * nc, [P_IN] * nc] + [[r] * nc for r in rans]:
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    # bc.txt takes the specific values (nuTilde, k, omega): the inlet multiplies them by its density
    turb = "".join(f" {r / rho:.10E}," for r in rans)
    inflow = f"    {MACH:.6E}, {T_IN:.6E}, {P_IN:.6E}, normal, normal, 1.000000E+00, 1.000000E+00,{turb}\n"
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


def read_wall(path):
    """Wall shear tauX and heat flux qw on the plate faces of a MOSE wall output (one zone, BLOCK packing)."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    plate = [k for k, t in enumerate(out["tauX"]) if t != 0.0]
    return [out["tauX"][k] for k in plate], [out["qw"][k] for k in plate]


def largest_relative(a, b):
    return max(abs(x - y) / abs(x) for x, y in zip(a, b))


def converged(w1, w2):
    t1, _ = read_wall(w1)
    t2, _ = read_wall(w2)
    d = largest_relative(t1, t2)
    print(f"    wall shear change over the last 2000 iterations of the first run: {d:.3e}")
    sys.exit(0 if d < 1.0e-11 else 1)


def compare(wa, wb):
    ta, qa = read_wall(wa)
    tb, qb = read_wall(wb)
    dt, dq = largest_relative(ta, tb), largest_relative(qa, qb)
    print(f"    wall shear: largest relative difference {dt:.3e}; wall heat flux: {dq:.3e} (plate faces {len(ta)})")
    sys.exit(0 if max(dt, dq) <= 1.0e-9 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "converged":
        converged(sys.argv[2], sys.argv[3])
    elif what == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
