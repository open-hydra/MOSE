"""Exact volumes of the cells of the unit check CellVolume (test/unit/cell_volume_check.f90).

Each cell is 8 nodes in the numbering of Compute_Metric_Tensor, N = 1 + d3 + 2 d2 + 4 d1 for the node at (d1, d2, d3) of the
unit cube, as binary64 values. Its volume is the integral of det J of the trilinear map of the nodes over the unit cube,
computed here in rational arithmetic from those binary64 values: the map, its Jacobian matrix and det J as polynomials in the
three cube coordinates with Fraction coefficients, then the integral of each monomial, 1/((a+1)(b+1)(c+1)). Standard library
only.

    python3 cell_volume_exact.py            prints the data block of cell_volume_check.f90
    python3 cell_volume_exact.py --table    prints, for each cell, the exact volume, det A of Compute_Metric_Tensor, the five
                                            tetrahedra in absolute value, whether det J keeps the sign of
                                            det A on the closed cell (all its Bernstein coefficients, degree 2 in each
                                            variable, of that sign or zero: a sufficient condition) and the twisted faces
"""
import itertools
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr


def node(N, d1, d2, d3):
    return N[1 + d3 + 2 * d2 + 4 * d1]


def comb(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            out[k] = out.get(k, 0) + va * vb
    return out


def padd(a, b, sb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sb * v
    return out


def det_j(N):
    """det J of the trilinear map of the nodes, a polynomial {(a, b, c): coefficient} in the cube coordinates."""
    q = {k: [Fr(c) for c in v] for k, v in N.items()}

    def w(var, d):
        e = [0, 0, 0]
        e[var] = 1
        return {tuple(e): Fr(1)} if d else {(0, 0, 0): Fr(1), tuple(e): Fr(-1)}

    x = [{}, {}, {}]
    for d1, d2, d3 in itertools.product((0, 1), repeat=3):
        basis = pmul(pmul(w(0, d1), w(1, d2)), w(2, d3))
        for m in range(3):
            x[m] = padd(x[m], {k: v * node(q, d1, d2, d3)[m] for k, v in basis.items()})

    def diff(p, var):
        out = {}
        for k, v in p.items():
            if k[var] > 0:
                kk = list(k)
                kk[var] -= 1
                out[tuple(kk)] = out.get(tuple(kk), 0) + v * k[var]
        return out

    J = [[diff(x[m], var) for var in range(3)] for m in range(3)]
    t1 = pmul(J[0][0], padd(pmul(J[1][1], J[2][2]), pmul(J[1][2], J[2][1]), -1))
    t2 = pmul(J[0][1], padd(pmul(J[1][0], J[2][2]), pmul(J[1][2], J[2][0]), -1))
    t3 = pmul(J[0][2], padd(pmul(J[1][0], J[2][1]), pmul(J[1][1], J[2][0]), -1))
    return padd(padd(t1, t2, -1), t3), q


def volume(N):
    dj, _ = det_j(N)
    return sum(v / ((k[0] + 1) * (k[1] + 1) * (k[2] + 1)) for k, v in dj.items())


def table_row(N):
    dj, q = det_j(N)
    V = sum(v / ((k[0] + 1) * (k[1] + 1) * (k[2] + 1)) for k, v in dj.items())
    coeffs = []
    for i, j, k in itertools.product(range(3), repeat=3):
        s = Fr(0)
        for a in range(i + 1):
            for b in range(j + 1):
                for c in range(k + 1):
                    p = dj.get((a, b, c), 0)
                    if p:
                        s += Fr(comb(i, a) * comb(j, b) * comb(k, c), comb(2, a) * comb(2, b) * comb(2, c)) * p
        coeffs.append(s)

    def row(pairs):
        return [sum(node(q, *h)[m] - node(q, *lo)[m] for h, lo in pairs) / 4 for m in range(3)]

    A1 = row([((1, b, c), (0, b, c)) for b in (0, 1) for c in (0, 1)])
    A2 = row([((a, 1, c), (a, 0, c)) for a in (0, 1) for c in (0, 1)])
    A3 = row([((a, b, 1), (a, b, 0)) for a in (0, 1) for b in (0, 1)])
    detA = (A1[0] * (A2[1] * A3[2] - A2[2] * A3[1]) - A1[1] * (A2[0] * A3[2] - A2[2] * A3[0])
            + A1[2] * (A2[0] * A3[1] - A2[1] * A3[0]))
    sgn = 1 if detA > 0 else -1
    keeps = all(sgn * cf >= 0 for cf in coeffs) and any(cf != 0 for cf in coeffs)

    def sub(a, b):
        return [a[m] - b[m] for m in range(3)]

    def cross(a, b):
        return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

    twisted = []
    for d in range(3):
        o = [(d + 1) % 3, (d + 2) % 3]
        for side in (0, 1):
            pts = []
            for a_, b_ in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                idx = [0, 0, 0]
                idx[d] = side
                idx[o[0]], idx[o[1]] = a_, b_
                pts.append(node(q, *idx))
            t = [pts[0][m] - pts[1][m] + pts[2][m] - pts[3][m] for m in range(3)]
            n = cross(sub(pts[2], pts[0]), sub(pts[3], pts[1]))
            if sum(t[m] * n[m] for m in range(3)) != 0:
                twisted.append("ijk"[d] + "-+"[side])
    # a sum of five tetrahedra, each in absolute value: exact only where the faces are plane
    # (docs/theory/numerics.md, Cell volume)
    v = {1: q[1], 2: q[5], 3: q[3], 4: q[7], 5: q[2], 6: q[6], 7: q[4], 8: q[8]}

    def tv(i1, i2, i3, i4):
        a, b, c = sub(v[i2], v[i1]), sub(v[i3], v[i1]), sub(v[i4], v[i1])
        cr = cross(b, c)
        return abs(sum(a[m] * cr[m] for m in range(3))) / 6

    v5 = tv(1, 2, 3, 5) + tv(2, 4, 3, 8) + tv(5, 8, 6, 2) + tv(5, 7, 8, 3) + tv(5, 8, 2, 3)
    return V, detA, v5, keeps, twisted


def cube(f=None):
    N = {}
    for d1, d2, d3 in itertools.product((0, 1), repeat=3):
        p = [float(d1), float(d2), float(d3)]
        N[1 + d3 + 2 * d2 + 4 * d1] = f(d1, d2, d3, p) if f else p
    return N


def moved(N, s):
    return {k: [v[0] + s[0], v[1] + s[1], v[2] + s[2]] for k, v in N.items()}


def cells():
    out = []
    out.append(("cube, one twisted face, T_n/h = 1",
                cube(lambda a, b, c, p: [p[0], p[1], p[2] + (1.0 if (a, b, c) == (1, 1, 1) else 0.0)])))
    out.append(("cube, one twisted face, T_n/h = 5",
                cube(lambda a, b, c, p: [p[0], p[1], p[2] + (5.0 if (a, b, c) == (1, 1, 1) else 0.0)])))
    pat = {(0, 0, 0): (0.05, -0.1, 0.0), (0, 0, 1): (-0.15, 0.05, 0.1), (0, 1, 0): (0.1, 0.05, -0.2),
           (0, 1, 1): (0.0, -0.1, 0.15), (1, 0, 0): (0.2, 0.1, 0.05), (1, 0, 1): (-0.05, -0.15, -0.1),
           (1, 1, 0): (-0.1, 0.2, 0.1), (1, 1, 1): (0.15, 0.0, -0.05)}
    out.append(("valid non-convex cell, six twisted faces",
                cube(lambda a, b, c, p: [p[m] + pat[(a, b, c)][m] for m in range(3)])))

    def collapsed(a, b, c, p):
        if (a, b, c) == (1, 0, 0):
            return [0.0, 0.0, 0.0]          # N5 on N1: the i edge at j = k = 0
        if (a, b, c) == (1, 1, 1):
            return [1.0, 1.0, 1.5]          # N8 raised
        return p

    out.append(("collapsed edge, two twisted faces", cube(collapsed)))
    quad = {(0, 0): (0.0, 0.0), (1, 0): (1.1, 0.1), (1, 1): (1.3, 0.9), (0, 1): (0.2, 1.2)}
    prism = cube(lambda a, b, c, p: [quad[(a, b)][0], quad[(a, b)][1], 0.5 * c])
    out.append(("2-D prism", prism))
    dl = math.radians(5.0)

    def wedge(mq):
        def f(a, b, c, p):
            x, r = mq[(a, b)]
            s = -1.0 if c == 0 else 1.0
            return [x, r * math.cos(dl / 2), s * r * math.sin(dl / 2)]
        return cube(f)

    out.append(("wedge cell, 5 degrees", wedge({(0, 0): (0.0, 0.8), (1, 0): (0.3, 0.75), (1, 1): (0.35, 1.1),
                                                (0, 1): (-0.05, 1.05)})))
    out.append(("wedge cell on the axis", wedge({(0, 0): (0.0, 0.0), (1, 0): (0.3, 0.0), (1, 1): (0.32, 0.2),
                                                 (0, 1): (-0.02, 0.21)})))
    out.append(("cube of cell 1 moved by 1e3", moved(out[0][1], (1000.0, -1000.0, 1000.0))))
    out.append(("2-D prism moved by 1e3", moved(prism, (1000.0, -1000.0, 1000.0))))
    out.append(("cube of cell 1 mirrored in x (left-handed)", {k: [-v[0], v[1], v[2]] for k, v in out[0][1].items()}))
    # det A > 0 at the centre and V < 0: an inverted cell that the centre does not show
    out.append(("inverted cell, det A > 0, V < 0",
                {1: [0.754373, 0.297540, 0.374380], 2: [-0.767401, -0.746476, 1.405035],
                 3: [0.295925, 1.736358, -0.344624], 4: [0.797913, 0.675172, 0.594738],
                 5: [0.344012, 0.293516, -0.571998], 6: [0.262633, 0.047749, 1.744357],
                 7: [0.206212, 1.467360, 0.169779], 8: [0.792210, 1.491336, 1.641459]}))
    return out


def dec(fr, digits=25):
    getcontext().prec = digits + 5
    return format(Decimal(fr.numerator) / Decimal(fr.denominator), f".{digits - 1}e")


def lit(x):
    """Fortran literal of a binary64 value, exact: repr gives the shortest decimal that rounds back to it."""
    return repr(x) + "_R8"


if __name__ == "__main__":
    if "--table" in sys.argv:
        print("cell | exact volume | det A | five tetrahedra | det J keeps the sign | twisted faces")
        for c, (name, N) in enumerate(cells(), 1):
            V, dA, v5, keeps, tw = table_row(N)
            print(f"{c} {name} | {float(V):.17g} | {float(dA):.12g} | {float(v5):.17g} | {keeps} | {' '.join(tw) or '-'}")
    else:
        for c, (name, N) in enumerate(cells(), 1):
            print(f"  names({c}) = '{name}'")
            for k in range(1, 9):
                print(f"  nodes(:,{k},{c}) = [ {lit(N[k][0])}, {lit(N[k][1])}, {lit(N[k][2])} ]")
            print(f"  vexact({c}) = {dec(volume(N))}_R8")
