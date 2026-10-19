"""Exact distance from a point to a bilinear face, the reference of the unit check BilinearDistance
(test/unit/bilinear_distance_check.f90, data written by bilinear_distance_cases.py), independent of the branch and bound of
MOSE (Bilinear_Face_Distance): the stationarity equations, in exact rational arithmetic, with Sturm sequences.

Face x(s,t) = c1 + s a + t b + s t T, (s,t) in [0,1]^2, a = c2 - c1, b = c4 - c1, T = c1 - c2 + c3 - c4 (corners in cyclic
order). f(s,t) = |x(s,t) - p|^2. The minimum over the closed square is attained at
  (i)   a corner or a point of the four edges (straight segments: closed form, exact);
  (ii)  an interior critical point: (x - p).x_s = 0 and (x - p).x_t = 0;
  (iii) a point of a pinched line t = t0 where x_s = a + t0 T = 0 (the whole line is one point of the surface): that point.
With u(t) = c1 - p + t b, v(t) = a + t T: x - p = u + s v, x_s = v, x_t = b + s T. Where v(t) != 0 the first equation
gives s = -N/D, N = u.v, D = v.v, and the second, times D^2, the polynomial of degree <= 5
  Q(t) = (u.b) D^2 - N D (v.b + u.T) + N^2 (v.T).
Every interior critical point with v(t) != 0 has Q(t) = 0. The real roots of Q in [0,1] are isolated with a Sturm
sequence in exact arithmetic and refined by bisection to a width of 2^-160; s and f there are evaluated exactly (the
error of f is second order in the width, since f is stationary). When Q vanishes identically (no isolated t), the
roles of s and t are swapped (the same construction in s). The inputs are the exact binary64 values (the Fraction of a
float is exact); the square root is taken in decimal to 60 digits. Standard library only.

ref_distance(p, c) -> (Decimal distance, info).
"""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import random

getcontext().prec = 60


# ------------------------------------------------------------------ polynomials with Fraction coefficients (low -> high)
def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def psub(p, q):
    return padd(p, [-c for c in q])


def pmul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return ptrim(r)


def pscale(p, c):
    return ptrim([c * a for a in p])


def peval(p, x):
    r = F(0)
    for a in reversed(p):
        r = r * x + a
    return r


def pderiv(p):
    return ptrim([i * p[i] for i in range(1, len(p))]) if len(p) > 1 else [F(0)]


def pdivmod(p, q):
    p = ptrim(p)
    q = ptrim(q)
    if len(q) == 1 and q[0] == 0:
        raise ZeroDivisionError
    out = [F(0)] * max(1, len(p) - len(q) + 1)
    r = list(p)
    while len(r) >= len(q) and not (len(r) == 1 and r[0] == 0):
        c = r[-1] / q[-1]
        k = len(r) - len(q)
        out[k] = c
        for i, b in enumerate(q):
            r[i + k] -= c * b
        r = ptrim(r[:-1]) if len(r) > 1 else [F(0)]
        if len(r) < len(q):
            break
    return ptrim(out), ptrim(r)


def is_zero(p):
    p = ptrim(p)
    return len(p) == 1 and p[0] == 0


def sturm_chain(p):
    p = ptrim(p)
    chain = [p, pderiv(p)]
    while not is_zero(chain[-1]):
        _, r = pdivmod(chain[-2], chain[-1])
        chain.append([-c for c in r])
    return chain[:-1]


def sign_changes(chain, x):
    vals = [peval(q, x) for q in chain]
    vals = [v for v in vals if v != 0]
    return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))


def squarefree(p):
    """p / gcd(p, p'): same roots, all simple."""
    a, b = ptrim(p), pderiv(p)
    if is_zero(b):
        return a
    while not is_zero(b):
        _, r = pdivmod(a, b)
        a, b = b, r
    g = a
    if len(g) == 1:
        return ptrim(p)
    q, r = pdivmod(p, g)
    assert is_zero(r)
    return q


def roots_in(p, lo=F(0), hi=F(1), width=F(1, 2 ** 160)):
    """Real roots of p in [lo, hi] (closed), each as a Fraction within width of a root (exact arithmetic)."""
    p = ptrim(p)
    if is_zero(p):
        raise ValueError("zero polynomial")
    if len(p) == 1:
        return []
    p = squarefree(p)
    if len(p) == 1:
        return []
    ch = sturm_chain(p)
    out = []

    def count(a, b):          # roots in (a, b]
        return sign_changes(ch, a) - sign_changes(ch, b)

    def rec(a, b):            # roots in (a, b]
        n = count(a, b)
        if n == 0:
            return
        if n == 1:
            # bisection on the sign of p at the right end (one simple root in (a, b]; p(a) may be zero, a root
            # outside (a, b], so the sign of p(b) is the reference)
            fb = peval(p, b)
            if fb == 0:
                out.append(b)
                return
            while b - a > width:
                m = (a + b) / 2
                fm = peval(p, m)
                if fm == 0:
                    out.append(m)
                    return
                if (fm > 0) == (fb > 0):
                    b, fb = m, fm
                else:
                    a = m
            out.append((a + b) / 2)
            return
        m = (a + b) / 2
        rec(a, m)
        rec(m, b)

    if peval(p, lo) == 0:
        out.append(lo)
    rec(lo, hi)
    return out


# ------------------------------------------------------------------ vectors of Fractions
def vF(x):
    return [F(float(c)) if not isinstance(c, F) else c for c in x]


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vdot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vsc(a, c):
    return [x * c for x in a]


def dsqrt(q):
    """Decimal square root of a nonnegative Fraction."""
    return (Decimal(q.numerator) / Decimal(q.denominator)).sqrt()


def seg_d2(p, a, b):
    ab = vsub(b, a)
    ap = vsub(p, a)
    l2 = vdot(ab, ab)
    t = F(0) if l2 == 0 else min(max(vdot(ap, ab) / l2, F(0)), F(1))
    w = vsub(ap, vsc(ab, t))
    return vdot(w, w)


def lin(c0, c1):
    """Vector polynomial c0 + c1 t as a list of 3 polynomials."""
    return [[c0[k], c1[k]] for k in range(3)]


def pvdot(u, v):
    r = [F(0)]
    for k in range(3):
        r = padd(r, pmul(u[k], v[k]))
    return r


def critical_t(c1, a, b, T, p):
    """Roots t in [0,1] of Q(t) (interior critical points with v(t) != 0), the pinch roots of D(t), or None if Q == 0."""
    u = lin(vsub(c1, p), b)
    v = lin(a, T)
    bb = [[b[k]] for k in range(3)]
    TT = [[T[k]] for k in range(3)]
    N = pvdot(u, v)
    D = pvdot(v, v)
    ub = pvdot(u, bb)
    vb = pvdot(v, bb)
    uT = pvdot(u, TT)
    vT = pvdot(v, TT)
    Q = psub(pmul(ub, pmul(D, D)), pmul(pmul(N, D), padd(vb, uT)))
    Q = padd(Q, pmul(pmul(N, N), vT))
    pinch = roots_in(D) if not is_zero(D) else None
    if is_zero(Q):
        return None, pinch, N, D
    return roots_in(Q), pinch, N, D


def ref_distance(p, c):
    """Exact-arithmetic minimum of |x(s,t) - p| over the bilinear face with corners c[0..3] (cyclic). Returns
    (Decimal distance, dict with the minimiser kind and parameters)."""
    p = vF(p)
    c = [vF(q) for q in c]
    c1, c2, c3, c4 = c
    a = vsub(c2, c1)
    b = vsub(c4, c1)
    T = vadd(vsub(c1, c2), vsub(c3, c4))
    best = None
    info = {}

    def x_at(s, t):
        return vadd(vadd(c1, vsc(a, s)), vadd(vsc(b, t), vsc(T, s * t)))

    def consider(d2, kind, st):
        nonlocal best, info
        if best is None or d2 < best:
            best = d2
            info = {"kind": kind, "st": st}

    # (i) the four edges (closed segments): corners included
    for k, (e0, e1) in enumerate(((c1, c2), (c2, c3), (c3, c4), (c4, c1))):
        consider(seg_d2(p, e0, e1), "edge", k)
    # (ii)-(iii) interior critical points, eliminating s (roots in t), or t when Q vanishes identically
    for swap in (False, True):
        if not swap:
            cc1, aa, bb_, TT = c1, a, b, T
        else:                 # x(s,t) with the roles of s and t exchanged: corners c1, c4, c3, c2
            cc1, aa, bb_, TT = c1, b, a, T
        roots, pinch, N, D = critical_t(cc1, aa, bb_, TT, p)
        if pinch:
            for t0 in pinch:      # a pinched line: one point of the surface
                xq = vadd(cc1, vsc(bb_, t0))
                w = vsub(xq, p)
                consider(vdot(w, w), "pinch", (swap, t0))
        if roots is None:
            continue              # Q identically zero: try the other elimination
        for t0 in roots:
            Dv = peval(D, t0)
            if Dv == 0:
                continue
            s0 = -peval(N, t0) / Dv
            if s0 < 0 or s0 > 1:
                continue
            ss, tt = (s0, t0) if not swap else (t0, s0)
            w = vsub(x_at(ss, tt), p)
            consider(vdot(w, w), "interior", (ss, tt))
        break
    else:
        info["note"] = "Q identically zero in both eliminations"
    return dsqrt(best), info
