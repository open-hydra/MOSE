"""Inputs of the MMSResiduals fast test, and the check of its output.

A constructed field for the Spalart-Allmaras equation with variable density. Square of side L = 1 mm, periodic in x and in y
(connections 201 on itself), N x N cells, null k faces; air of test/common/Air-viscous (mu and k constant) at rest at 1 bar,
with smooth periodic temperature and nu~:
  T   = T0 [1 + 0.3 sin(2 pi x/L) (1 + 0.4 sin(2 pi y/L))],                       T0 = 300 K, rho = p/(R T);
  nu~ = nu0 [1 + 0.5 cos(2 pi x/L + 0.7) (1 + 0.4 cos(2 pi y/L + 0.3))],          nu0 = 1e-4 m2/s.
There is no wall, so the wall distance is the sentinel 1e8 m and the production and destruction of nu~ vanish. With the fluid at
rest the right-hand side of the conservative equation for rho nu~ at t = 0 is its diffusion part, which in the standard form
(the one of the NASA Turbulence Modeling Resource, Spalart-Allmaras page, conservation form) reads
  (1/sigma) [ div( rho (nu + nu~) grad nu~ ) + rho cb2 |grad nu~|^2 - (nu + nu~) grad rho . grad nu~ ]
    = (rho/sigma) [ div( (nu + nu~) grad nu~ ) + cb2 |grad nu~|^2 ],
nu = mu/rho. The last term on the first line is what makes the conservation form equal to the original equation.

The solver makes one explicit Euler step with local time stepping; the diagnostic output holds the increment of rho nu~ and the
step of each cell, so their ratio is the right-hand side the solver evaluated. It must match the exact one above (derivatives of
the analytic fields by fourth-order differences) to second order: on N = 32 and N = 64 the largest difference over the cells,
divided by the largest value of the right-hand side, must be below 0.02 on N = 64 and fall with an observed order above 1.5.

Usage (from the test directory):
  python3 -B make_case.py case N             INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py error FILE         prints the error of the right-hand side in OUTPUT/diagnostic.tec
  python3 -B make_case.py check E32 E64      exit status 0 when the errors on N = 32 and N = 64 pass
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P0, L = 300.0, 1.0e5, 1.0e-3
SIGMA, CB2 = 2.0 / 3.0, 0.622       # Spalart-Allmaras
NU0 = 1.0e-4


def temp(x, y):
    return T0 * (1.0 + 0.3 * math.sin(2 * math.pi * x / L) * (1.0 + 0.4 * math.sin(2 * math.pi * y / L)))


def rho(x, y):
    return P0 / (RGAS * temp(x, y))


def nut(x, y):
    return NU0 * (1.0 + 0.5 * math.cos(2 * math.pi * x / L + 0.7) * (1.0 + 0.4 * math.cos(2 * math.pi * y / L + 0.3)))


H = L * 1.0e-3                      # step of the fourth-order differences of the analytic fields


def ddx(f, x, y):
    return (-f(x + 2 * H, y) + 8 * f(x + H, y) - 8 * f(x - H, y) + f(x - 2 * H, y)) / (12 * H)


def ddy(f, x, y):
    return (-f(x, y + 2 * H) + 8 * f(x, y + H) - 8 * f(x, y - H) + f(x, y - 2 * H)) / (12 * H)


def rhs(x, y):
    """Exact right-hand side of the standard conservative equation for rho nu~ at t = 0, fluid at rest, no wall."""
    def diff(a, b):
        return MU / rho(a, b) + nut(a, b)

    div = ddx(lambda a, b: diff(a, b) * ddx(nut, a, b), x, y) + ddy(lambda a, b: diff(a, b) * ddy(nut, a, b), x, y)
    g2 = ddx(nut, x, y) ** 2 + ddy(nut, x, y) ** 2
    return rho(x, y) / SIGMA * (div + CB2 * g2)


def write_case(n):
    h = L / n
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    cells = [((i + 0.5) * h, (j + 0.5) * h) for j in range(n) for i in range(n)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={n + 1}, J={n + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(n + 1):
                    for i in range(n + 1):
                        f.write(f" {(i * h, j * h, k * h)[c]:.16E}\n")
        zero = [0.0] * (n * n)
        for col in ([rho(x, y) for x, y in cells], zero, zero, zero, [P0] * (n * n),
                    [rho(x, y) * nut(x, y) for x, y in cells]):
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(b, i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, face, t)

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    out = []
    for j in range(1, n + 1):
        out += [rec(1, 1, j, 1, 1, 201), con(1, n, j, 1, 2, 1, 0, 0, 1)]
    for j in range(1, n + 1):
        out += [rec(1, n, j, 1, 2, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    for i in range(1, n + 1):
        out += [rec(1, i, 1, 1, 3, 201), con(1, i, n, 1, 4, 1, 0, 0, 1)]
    for i in range(1, n + 1):
        out += [rec(1, i, n, 1, 4, 201), con(1, i, 1, 1, 3, 1, 0, 0, 1)]
    for face in (5, 6):
        out += [rec(1, i, j, 1, face, 0) for j in range(1, n + 1) for i in range(1, n + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the number of cells per side."""
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
    return out, ni - 1


def error(path):
    """Largest difference between the right-hand side of rho nu~ evaluated by the solver (increment over step) and the exact
    one, over the largest exact value."""
    d, n = read_zone(path)
    h = L / n
    cells = [((i + 0.5) * h, (j + 0.5) * h) for j in range(n) for i in range(n)]
    got = [r / t for r, t in zip(d["rhoNuTilde"], d["dt"])]
    exact = [rhs(x, y) for x, y in cells]
    e = max(abs(g - x) for g, x in zip(got, exact)) / max(abs(x) for x in exact)
    print(f"{e:.6e}")


def check(e32, e64):
    p = math.log2(e32 / e64) if e32 > 0 and e64 > 0 else float("nan")
    print(f"    error of the right-hand side of rho nu~: {e32:.3e} (N = 32), {e64:.3e} (N = 64), observed order {p:.2f}")
    sys.exit(0 if e64 < 0.02 and p > 1.5 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(int(sys.argv[2]))
    elif what == "error":
        error(sys.argv[2])
    elif what == "check":
        check(float(sys.argv[2]), float(sys.argv[3]))
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
