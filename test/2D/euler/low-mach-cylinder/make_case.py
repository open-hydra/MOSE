"""Inputs of the Preconditioning tests, and the check of their output.

Inviscid flow past a circular cylinder at free-stream Mach number 0.01 (air of test/common/Air, calorically perfect,
cp = 1004.5 J/(kg K); p_inf = 1e5 Pa, T_inf = 300 K), in a 2-D O-grid of one block:
  i  around the cylinder, NI cells, clockwise from the downstream axis (so that i, j, k are right-handed); faces 1 and
     2 of the block meet on the downstream axis and are joined by a block connection (101) on themselves;
  j  radial, NJ cells from the cylinder (r = R, face 3, slip wall 300) to r = R_OUT (face 4), the cells growing
     geometrically from a first cell as wide as it is long;
  k  one cell (faces 5 and 6 null, type 0).
Face 4 is an inlet at total conditions (401: T0, p0 of the free stream, flow along x) where its face centre has
cos(theta) < -0.5, the upstream third of the circle, and an outlet (406 at p_inf) elsewhere. The initial field is the
uniform free stream.

At this Mach number the flow is incompressible to O(M^2) = 1e-4: the pressure coefficient of potential flow,
Cp = (p - p_inf)/q_inf = 2 (R/r)^2 cos(2 theta) - (R/r)^4, with q_inf = p0 - p_inf, holds in the cells to that order,
to the truncation error of the grid and to the blockage of the far boundary, O((R/R_OUT)^2). An upwind flux whose
dissipation scales with the speed of sound gives pressure fluctuations of order M p_inf instead of M^2 p_inf, so Cp
errors of order 1/M; the preconditioned flux (HLLC-PC) scales its dissipation with Ur.

The preconditioned scheme (integration-variables = prec, HLLC-PC) has an incompressible limit: with Ur proportional to
the flow speed (preconditioning-eps-min proportional to the Mach number), its steady Cp field does not depend on the Mach
number, to O(M^2). The cell Cp of potential flow is not reached on this grid: the dissipation of the scheme separates
the flow at the back of the cylinder (a wake of total-pressure loss), at any Mach number, so the front of the cylinder
is compared with potential flow within a loose bound only.

Usage (from the test directory; MACH, in the environment, sets the free-stream Mach number, 0.01 by default):
  python3 -B make_case.py case                          writes INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py converged HISTORY FALL        exit 0 when the density residual of the last row of the
                                                        residual history is below FALL times that of its first row
  python3 -B make_case.py same F1 M1 F2 M2 TOL          exit 0 when the Cp of the solutions F1 (Mach number M1) and F2
                                                        (M2) differ by less than TOL in every cell with r < 1.5 R
  python3 -B make_case.py differs F1 M1 F2 M2 MIN       exit 0 when they differ by more than MIN in some cell there
  python3 -B make_case.py front F M TOL                 exit 0 when the Cp of the first row of cells of the front half
                                                        (cos theta < 0) is within TOL of potential flow
"""
import math
import os
import sys

NI, NJ = 96, 48
R, R_OUT = 0.5, 20.0
DZ = 0.01
CP, RGAS = 1004.5, 8314.46261815324 / 28.970418
GAMMA = CP / (CP - RGAS)
P_INF, T_INF = 1.0e5, 300.0
MACH = float(os.environ.get('MACH', '0.01'))       # free-stream Mach number (environment variable MACH)
P0 = P_INF * (1.0 + 0.5 * (GAMMA - 1.0) * MACH ** 2) ** (GAMMA / (GAMMA - 1.0))
T0 = T_INF * (1.0 + 0.5 * (GAMMA - 1.0) * MACH ** 2)
U_INF = MACH * math.sqrt(GAMMA * RGAS * T_INF)
RHO_INF = P_INF / (RGAS * T_INF)


def radii():
    """NJ + 1 node radii from R to R_OUT, the first cell as wide as it is long."""
    d0 = 2.0 * math.pi * R / NI
    lo, hi = 1.0, 2.0
    for _ in range(200):                       # growth ratio q: d0 (q^NJ - 1)/(q - 1) = R_OUT - R
        q = 0.5 * (lo + hi)
        if d0 * (q ** NJ - 1.0) / (q - 1.0) > R_OUT - R:
            hi = q
        else:
            lo = q
    r = [R]
    for j in range(NJ):
        r.append(r[-1] + d0 * q ** j)
    r[-1] = R_OUT
    return r


def nodes():
    """x, y of the nodes, i fastest; theta decreases with i (clockwise)."""
    r = radii()
    xs, ys = [], []
    for j in range(NJ + 1):
        for i in range(NI + 1):
            th = -2.0 * math.pi * i / NI
            xs.append(r[j] * math.cos(th))
            ys.append(r[j] * math.sin(th))
    return xs, ys


def fresh(path):
    if os.path.lexists(path):
        os.remove(path)
    return open(path, 'w')


def write_case():
    xs, ys = nodes()
    nc = NI * NJ
    with fresh('INPUT/ic.tec') as f:
        f.write(' VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"\n')
        f.write(' ZONE  T = B1-IG, I=%d, J=%d, K=2, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n' % (NI + 1, NJ + 1))
        for col in (xs + xs, ys + ys, [0.0] * len(xs) + [DZ] * len(xs)):
            f.write(''.join(' %.16E\n' % v for v in col))
        for v in (RHO_INF, U_INF, 0.0, 0.0, P_INF):
            f.write(' %.16E\n' % v * nc)
    rec = lambda i, j, face, t: ''.join('%8d' % v for v in (1, i, j, 1, face, t)) + '\n'
    inflow = '    %.17E, %.17E, 0.0, 0.0, 1.0, 1.0,\n' % (T0, P0)
    outflow = '    %.17E, 1.0,\n' % P_INF
    out = []
    for j in range(1, NJ + 1):                 # the seam: face 1 of cell (1, j) on face 2 of cell (NI, j)
        out += [rec(1, j, 1, 101), '%10d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (1, NI, j, 1, 2, 1, 0, 0, 1)]
        out += [rec(NI, j, 2, 101), '%10d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (1, 1, j, 1, 1, 1, 0, 0, 1)]
    for i in range(1, NI + 1):
        out.append(rec(i, 1, 3, 300))
    for i in range(1, NI + 1):
        th = -2.0 * math.pi * (i - 0.5) / NI
        out += [rec(i, NJ, 4, 401), inflow] if math.cos(th) < -0.5 else [rec(i, NJ, 4, 406), outflow]
    for face in (5, 6):
        for j in range(1, NJ + 1):
            for i in range(1, NI + 1):
                out.append(rec(i, j, face, 0))
    with fresh('INPUT/bc.txt') as f:
        f.write(''.join(out))


def read_field(path):
    """Cell centres (x, y), pressure and row j of the cells of the single zone of an ASCII Tecplot solution of MOSE:
    a VARIABLES line, a ZONE line, then the values (DATAPACKING=BLOCK, the coordinates nodal, the rest cell-centred)."""
    with open(path) as f:
        names = f.readline().split('=', 1)[1].replace('"', ' ').split()
        zone = f.readline()
        vals = [float(t) for t in f.read().split('DATASETAUXDATA')[0].split()]   # the restart records follow
    dims = {d: int(zone.split(' %s=' % d)[1].split(',')[0]) for d in ('I', 'J', 'K')}
    ni, nj, nk = dims['I'], dims['J'], dims['K']
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    x, y = vals[0:nn], vals[nn:2 * nn]
    p_at = 3 * nn + (names.index('p') - 3) * nc
    p = vals[p_at:p_at + nc]
    cells = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            n = [j * ni + i, j * ni + i + 1, (j + 1) * ni + i, (j + 1) * ni + i + 1]
            cells.append((sum(x[m] for m in n) / 4.0, sum(y[m] for m in n) / 4.0, p[j * (ni - 1) + i], j))
    return cells


def cp_field(path, mach):
    """Cp = (p - p_inf)/(p0 - p_inf) of the cells, with the free stream of Mach number mach; (cp, x, y, row j)."""
    p0 = P_INF * (1.0 + 0.5 * (GAMMA - 1.0) * mach ** 2) ** (GAMMA / (GAMMA - 1.0))
    return [((p - P_INF) / (p0 - P_INF), cx, cy, j) for cx, cy, p, j in read_field(path)]


def residual_fall(path):
    rows = [l.split() for l in open(path) if l.strip() and not l.lstrip().startswith(('#', 'V', 'Z', '"'))]
    first, last = float(rows[0][2]), float(rows[-1][2])
    return last / first


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    what, args = sys.argv[1], sys.argv[2:]
    if what == 'case':
        write_case()
        return
    if what == 'converged':
        fall = residual_fall(args[0])
        print('  density residual: last row / first row = %.3e (required below %s)' % (fall, args[1]))
        sys.exit(0 if fall < float(args[1]) else 1)
    if what in ('same', 'differs'):
        a = cp_field(args[0], float(args[1]))
        b = cp_field(args[2], float(args[3]))
        d = max(abs(x[0] - y[0]) for x, y in zip(a, b) if math.hypot(x[1], x[2]) < 1.5 * R)
        print('  largest |Cp(M = %s) - Cp(M = %s)| over r < 1.5 R: %.3e' % (args[1], args[3], d))
        if what == 'same':
            sys.exit(0 if d < float(args[4]) else 1)
        sys.exit(0 if d > float(args[4]) else 1)
    if what == 'front':
        worst = 0.0
        for cp, cx, cy, j in cp_field(args[0], float(args[1])):
            if j == 0 and cx < 0.0:
                r, th = math.hypot(cx, cy), math.atan2(cy, cx)
                worst = max(worst, abs(cp - (2.0 * (R / r) ** 2 * math.cos(2.0 * th) - (R / r) ** 4)))
        print('  largest |Cp - Cp_potential| on the first row of the front half: %.3f' % worst)
        sys.exit(0 if worst < float(args[2]) else 1)
    sys.exit(__doc__)


if __name__ == '__main__':
    main()
