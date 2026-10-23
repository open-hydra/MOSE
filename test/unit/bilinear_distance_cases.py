"""Cases of the unit check BilinearDistance (test/unit/bilinear_distance_check.f90), written to bilinear_distance_cases.txt.
Standard library only; every coordinate is a binary64 value written exactly (repr), every reference is exact
(bilinear_distance_exact.py for the distance, rational arithmetic for the gap of the fast path).

D lines, the certified minimum (Bilinear_Face_Distance):  D label px py pz c1x c1y c1z ... c4z d_exact
  faces twisted (twist 1, 0.1, 1e-3, 1e-6), plane convex, plane dart (a reflex corner), a collapsed edge, two collapsed edges
  (a segment), one point, a pinch (two opposite sides meeting), turned at random and moved by 0 or 1e3, with points on the
  normal at heights 1e-1, 1e-3, 1e-9, in a box around the face and beyond the sides and corners; degenerate faces (a collapsed
  edge, nearly collapsed by 1e-6, 1e-9, 1e-12, two collapsed edges, three coincident corners, a bowtie, a pinch, a dart, a
  needle of aspect 1e6) with points next to their degenerate feature; the focal point of a saddle and next to it; the
  adverse cases of the review of the method (beyond the apex of a collapsed edge and on it, beyond the tip of a needle and
  over it at 1e-9, a nearly plane face far from the origin under the threshold of the fast path, points on a twisted face
  and on a corner, a bowtie, a pinch, the focus of a saddle, beyond the corners).
F lines, the fast path (Face_Fast_Test):  F label c1x c1y c1z ... c4z fast gap_exact
  plane faces, which must take it: convex quadrilaterals in a plane of dyadic slopes (coplanar to the bit) and turned
  exactly then rounded to binary64; thin tapered trapezoids of aspect 10 to 1e5 (kappa_n up to 1e5), both ways; triangles
  of a collapsed edge, also far from the origin (moved by 1e3, sides 1e-3 to 1); parallelograms, 2-D extruded faces and
  turned parallelograms rounded to binary64; faces that must not: twisted faces and plane darts. gap_exact is the gap g* of
  the four triangles from the face in exact arithmetic, min(|T|/16, |T.n|/(2|n|) when the projection is convex). A face is
  kept only when its outcome is certain in the model of 2u per operation (margins of the proof in
  docs/theory/turbulence.md): a plane face with every corner turning well beyond the rounding of the test, a twisted face
  or a dart whose gap is well beyond the threshold.
Usage: python3 bilinear_distance_cases.py [OUT=bilinear_distance_cases.txt]"""
import math
import os
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bilinear_distance_exact as R  # noqa: E402

getcontext().prec = 60
U = 2.0 ** -53
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                         "bilinear_distance_cases.txt")


# ------------------------------------------------------------------ binary64 vectors (as the faces are built)
def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scl(a, c):
    return (a[0] * c, a[1] * c, a[2] * c)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def norm(a):
    return math.sqrt(dot(a, a))


def rot(v, Rm):
    return tuple(sum(Rm[i][j] * v[j] for j in range(3)) for i in range(3))


def random_rotation(rng):
    q = [rng.gauss(0, 1) for _ in range(4)]
    n = math.sqrt(sum(x * x for x in q))
    w, x, y, z = (v / n for v in q)
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]


def make_face(rng, kind, tw, shift):
    if kind == "twist":
        base = [(0, 0), (1, 0), (1, 1), (0, 1)]
        c = [(x + rng.uniform(-0.2, 0.2), y + rng.uniform(-0.2, 0.2), rng.uniform(-tw, tw)) for x, y in base]
    elif kind == "plane":
        base = [(0, 0), (1, 0), (1, 1), (0, 1)]
        c = [(x + rng.uniform(-0.2, 0.2), y + rng.uniform(-0.2, 0.2), 0.0) for x, y in base]
    elif kind == "dart":
        c = [(0, 0, 0), (1, 0, 0), (0.35 + rng.uniform(-0.05, 0.05), 0.3 + rng.uniform(-0.05, 0.05), 0), (0, 1, 0)]
    elif kind == "collapsed":
        c = [(0, 0, 0), (0, 0, 0), (1 + rng.uniform(-0.2, 0.2), 1, rng.uniform(-tw, tw)), (0, 1, rng.uniform(-tw, tw))]
    elif kind == "segment":
        c = [(0, 0, 0), (0, 0, 0), (1, 1, rng.uniform(-tw, tw)), (1, 1, rng.uniform(-tw, tw))]
        c[3] = c[2]
    elif kind == "point":
        c = [(0.3, 0.2, 0.1)] * 4
    elif kind == "pinch":
        a = (1.0, 0.2, 0.1)
        t0 = 0.4
        c4 = (0.1, 1.0, tw)
        c = [(0.0, 0.0, 0.0), a, add(c4, scl(a, -(1 - t0) / t0)), c4]
    else:
        raise ValueError(kind)
    c = [tuple(float(v) for v in q) for q in c]
    Rm = random_rotation(rng)
    sh = tuple(shift * rng.uniform(-1, 1) for _ in range(3))
    return [add(rot(q, Rm), sh) for q in c]


def bil_point(c, s, t):
    a = sub(c[1], c[0])
    b = sub(c[3], c[0])
    T = sub(add(c[0], c[2]), add(c[1], c[3]))
    return add(c[0], add(add(scl(a, s), scl(b, t)), scl(T, s * t)))


def bil_normal(c, s, t):
    a = sub(c[1], c[0])
    b = sub(c[3], c[0])
    T = sub(add(c[0], c[2]), add(c[1], c[3]))
    n = cross(add(a, scl(T, t)), add(b, scl(T, s)))
    nn = norm(n)
    return scl(n, 1.0 / nn) if nn > 0 else (0.0, 0.0, 1.0)


def make_point(rng, c, mode, h):
    if mode == "normal":
        s, t = rng.uniform(0.05, 0.95), rng.uniform(0.05, 0.95)
        return add(bil_point(c, s, t), scl(bil_normal(c, s, t), h * rng.choice((-1, 1))))
    if mode == "box":
        cen = scl(add(add(c[0], c[1]), add(c[2], c[3])), 0.25)
        return add(cen, tuple(rng.uniform(-1.5, 1.5) for _ in range(3)))
    if mode == "beyond":
        s, t = rng.uniform(-0.6, 1.6), rng.uniform(-0.6, 1.6)
        return add(bil_point(c, s, t), scl(bil_normal(c, min(max(s, 0), 1), min(max(t, 0), 1)), h))
    raise ValueError(mode)


def degenerate_faces():
    tw = 0.3
    out = [("collapsed", [(0, 0, 0), (0, 0, 0), (1.1, 1.0, tw), (0.0, 1.0, -tw)])]
    for eps in (1e-6, 1e-9, 1e-12):
        out.append((f"near_collapsed{eps:g}", [(0, 0, 0), (eps, 0.3 * eps, 0), (1.1, 1.0, tw), (0.0, 1.0, -tw)]))
    out.append(("two_collapsed", [(0, 0, 0), (0, 0, 0), (1.0, 1.0, tw), (1.0, 1.0, tw)]))
    out.append(("three", [(0, 0, 0), (0, 0, 0), (0, 0, 0), (0.2, 1.0, tw)]))
    out.append(("bowtie", [(0, 0, 0), (1.0, 0.1, tw), (0, 0, 0), (0.1, 1.0, -tw)]))
    a = (1.0, 0.2, 0.1)
    t0 = 0.4
    c4 = (0.1, 1.0, tw)
    out.append(("pinch", [(0, 0, 0), a, add(c4, scl(a, -(1 - t0) / t0)), c4]))
    out.append(("dart", [(0, 0, 0), (1, 0, 0), (0.35, 0.3, 0), (0, 1, 0)]))
    out.append(("needle", [(0, 0, 0), (1.0, 0, 1e-7), (1.0, 1e-6, 0), (0, 1e-6, 1e-7)]))
    return out


def feature_points(name, c):
    pts = []
    if name.startswith("collapsed") or name.startswith("near_collapsed") or name in ("two_collapsed", "three"):
        pts.append(add(c[0], (-0.05, -0.07, 0.02)))
    if name == "bowtie":
        pts.append(add(c[0], (0.01, 0.02, 0.5)))
    if name == "pinch":
        pp = add(scl(c[0], 0.6), scl(c[3], 0.4))
        pts.append(add(pp, (0.0, 0.0, 0.2)))
        pts.append(add(pp, (1e-3, -2e-3, 1e-3)))
    return pts


def distance_cases():
    rng = random.Random(2026)
    cases = []

    def put(label, p, c):
        cases.append((label, tuple(float(v) for v in p), [tuple(float(v) for v in q) for q in c]))

    for kind, tw in [("twist", 1.0), ("twist", 0.1), ("twist", 1e-3), ("twist", 1e-6), ("plane", 0.0), ("dart", 0.0),
                     ("collapsed", 0.3), ("segment", 0.3), ("point", 0.0), ("pinch", 0.3)]:
        for shift in (0.0, 1e3):
            for mode, h in [("normal", 1e-1), ("normal", 1e-3), ("normal", 1e-9), ("box", 0.0), ("beyond", 1e-3)]:
                c = make_face(rng, kind, tw, shift)
                p = make_point(rng, c, mode, h) if kind != "point" else add(c[0], (0.1, 0.2, 0.3))
                put(f"{kind}{tw:g}-{mode}{h:g}-s{shift:g}", p, c)
    for name, c0 in degenerate_faces():
        for shift in (0.0, 1e3):
            Rm = random_rotation(rng)
            sh = tuple(shift * rng.uniform(-1, 1) for _ in range(3))
            c = [add(rot(tuple(float(v) for v in q), Rm), sh) for q in c0]
            pts = [("feature", q) for q in feature_points(name, c)]
            for hh in (1e-3, 1e-9):
                pts.append((f"normal{hh:g}", add(bil_point(c, 0.37, 0.61), scl(bil_normal(c, 0.37, 0.61), hh))))
            cen = scl(add(add(c[0], c[1]), add(c[2], c[3])), 0.25)
            pts.append(("box", add(cen, tuple(rng.uniform(-1.5, 1.5) for _ in range(3)))))
            for lab, p in pts:
                put(f"{name}-{lab}-s{shift:g}", p, c)
    for k in (0.5, 2.0):
        c = [(-1.0, -1.0, k), (1.0, -1.0, -k), (1.0, 1.0, k), (-1.0, 1.0, -k)]
        for h in (0.999 / k, 1.0 / k, 1.001 / k):
            put(f"focal{k:g}-h{h:.6f}", (0.0, 0.0, h), c)
    # adverse cases of the review of the method
    rng2 = random.Random(77)
    c = [(0, 0, 0), (0, 0, 0), (1.1, 1.0, 0.3), (0.0, 1.0, -0.3)]
    put("adverse-collapsed-beyond-apex", (-0.1, -0.5, 0.05), c)
    put("adverse-collapsed-at-apex", (0.0, 0.0, 0.0), c)
    c = [(0, 0, 0), (1.0, 0, 1e-7), (1.0, 1e-6, 0), (0, 1e-6, 1e-7)]
    put("adverse-needle-beyond-tip", (1.5, 2e-7, 1e-3), c)
    put("adverse-needle-over-centre-1e-9", (0.5, 5e-7, 1e-9), c)
    for frac in (0.9, 0.5):
        sz = 1e-3
        tn = 2 * frac * (64 * 2 * U * 1e3)
        c = [(0.0, 0.0, 0.0), (sz, 0.0, 0.0), (sz, sz, tn), (0.0, sz, 0.0)]
        c = [add(q, (1e3, -1e3, 1e3)) for q in c]
        for h in (1e-6, 1e-7):
            put(f"adverse-nearly-plane-far-{frac}-h{h:g}", add(add(c[0], (sz / 2, sz / 2, 0.0)), (0.0, 0.0, h)), c)
    c = make_face(rng2, "twist", 0.5, 0.0)
    put("adverse-on-twisted-face", bil_point(c, 0.3, 0.7), c)
    put("adverse-on-corner", c[2], c)
    c = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.3), (1.0, 1.0, -0.3)]
    put("adverse-bowtie-near-crossing", (0.5, 0.5, 0.2), c)
    put("adverse-bowtie-over-1e-7", (0.5, 0.5, 1e-7), c)
    a = (1.0, 0.2, 0.1)
    c4 = (0.1, 1.0, 0.3)
    c = [(0.0, 0.0, 0.0), a, add(c4, scl(a, -(1 - 0.4) / 0.4)), c4]
    pp = add(scl(c[0], 0.6), scl(c[3], 0.4))
    put("adverse-pinch-at-pinch", pp, c)
    put("adverse-pinch-beyond", add(pp, (-0.3, 0.1, 0.05)), c)
    c = [(-1.0, -1.0, 2.0), (1.0, -1.0, -2.0), (1.0, 1.0, 2.0), (-1.0, 1.0, -2.0)]
    put("adverse-saddle-focus", (0.0, 0.0, 0.5), c)
    for it in range(4):
        c = make_face(rng2, "twist", rng2.choice((1.0, 0.1)), rng2.choice((0.0, 1e3)))
        k_ = rng2.randint(0, 3)
        put(f"adverse-beyond-corner{it}", add(c[k_], scl(sub(c[k_], scl(add(add(c[0], c[1]), add(c[2], c[3])), 0.25)),
                                                         rng2.uniform(0.2, 2.0))), c)
    return cases


# ------------------------------------------------------------------ the fast path: faces and their exact gap
def dsq(fr):
    return (Decimal(fr.numerator) / Decimal(fr.denominator)).sqrt()


def fsub(a, b):
    return [a[i] - b[i] for i in range(3)]


def fcross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def fdot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def gap_and_outcome(c):
    """Exact gap g* of the four triangles and the outcome of Face_Fast_Test when it is certain in the model of 2u per
    operation (True, False), or None. First-order rounding bounds of the test (docs/theory/turbulence.md): E1, E2, Ew, nu."""
    q = [[Fr(x) for x in p] for p in c]
    e = [fsub(q[(k + 1) % 4], q[k]) for k in range(4)]
    d1, d2 = fsub(q[2], q[0]), fsub(q[3], q[1])
    T = [-(e[0][i] + e[2][i]) for i in range(3)]
    n = fcross(d1, d2)
    u = Decimal(2) ** -53
    nT = dsq(fdot(T, T))
    nn2 = fdot(n, n)
    nn = dsq(nn2)
    ne = [dsq(fdot(x, x)) for x in e]
    M = Decimal(max(abs(x) for p in c for x in p))
    g1 = nT / 16
    E1 = u / 8 * (ne[0] + ne[2]) + 7 * u / 16 * nT
    convex = nn2 > 0 and all(fdot(fcross(e[k], e[(k + 1) % 4]), n) > 0 or fdot(e[k], e[k]) == 0
                             or fdot(e[(k + 1) % 4], e[(k + 1) % 4]) == 0 for k in range(4))
    g2 = None
    if convex:
        g2 = abs(Decimal(fdot(T, n).numerator) / Decimal(fdot(T, n).denominator)) / (2 * nn)
    g = min(g1, g2) if g2 is not None else g1
    thr = 8 * u * M
    # certain to take it: branch 1 with margin, or branch 2 with a normal known to 1/10, every corner turning beyond 3.5 Ew
    yes = g1 <= thr + Decimal("0.8") * E1
    if convex and not yes:
        nu = (4 + 4 * Decimal(2).sqrt()) * u * dsq(fdot(d1, d1)) * dsq(fdot(d2, d2))
        E2 = u * (ne[0] + ne[2]) + Decimal("7.5") * u * nT + nT * nu / nn
        ok_corners = True
        for k in range(4):
            k1 = (k + 1) % 4
            if ne[k] == 0 or ne[k1] == 0:
                continue
            w = fdot(fcross(e[k], e[k1]), n)
            Ew = ne[k] * ne[k1] * ((10 + 4 * Decimal(2).sqrt()) * u * nn + nu)
            if Decimal(w.numerator) / Decimal(w.denominator) < Decimal("3.5") * Ew:
                ok_corners = False
        yes = nu <= nn / 10 and ok_corners and g2 <= thr + Decimal("0.8") * E2
    if yes:
        return g, True
    # certain not to take it: each branch beyond the threshold by four times its rounding
    no = g1 > thr + 4 * E1
    if no and convex:
        nu = (4 + 4 * Decimal(2).sqrt()) * u * dsq(fdot(d1, d1)) * dsq(fdot(d2, d2))
        E2 = u * (ne[0] + ne[2]) + Decimal("7.5") * u * nT + nT * nu / nn
        no = g2 > thr + 4 * E2
    if no:
        return g, False
    return g, None


def rot_rational(rng):
    while True:
        a, b, c, d = (rng.randint(-9, 9) for _ in range(4))
        s = a * a + b * b + c * c + d * d
        if s > 0:
            break
    return [[Fr(a * a + b * b - c * c - d * d, s), Fr(2 * (b * c - a * d), s), Fr(2 * (b * d + a * c), s)],
            [Fr(2 * (b * c + a * d), s), Fr(a * a - b * b + c * c - d * d, s), Fr(2 * (c * d - a * b), s)],
            [Fr(2 * (b * d - a * c), s), Fr(2 * (c * d + a * b), s), Fr(a * a - b * b - c * c + d * d, s)]]


def dy(x, bits=20):
    if x == 0:
        return Fr(0)
    e = math.floor(math.log2(abs(x))) - bits
    return Fr(round(x / 2.0 ** e)) * Fr(2) ** e


def place(rng, quad2d, scale, offset, kind):
    """kind 'b': the plane quadrilateral turned exactly, scaled, moved, then rounded to binary64; kind 'a': on the plane
    z = a x + b y + c of dyadic slopes, every coordinate representable (coplanar to the bit)."""
    if kind == "b":
        Rm = rot_rational(rng)
        return [[float(Rm[i][0] * x * scale + Rm[i][1] * y * scale + offset[i]) for i in range(3)] for x, y in quad2d]
    mag = max(abs(float(o)) for o in offset) + 1.5 * float(scale)
    g = Fr(2) ** (math.floor(math.log2(mag)) - 44)

    def qz(v):
        return Fr(round(v / g)) * g
    for _ in range(200):
        a = Fr(rng.randint(-40, 40), 32)
        b = Fr(rng.randint(-40, 40), 32)
        out, ok = [], True
        for x, y in quad2d:
            xx, yy = qz(x * scale + offset[0]), qz(y * scale + offset[1])
            zz = a * xx + b * yy + qz(offset[2])
            if Fr(float(xx)) != xx or Fr(float(yy)) != yy or Fr(float(zz)) != zz:
                ok = False
                break
            out.append([float(xx), float(yy), float(zz)])
        if ok:
            return out
    return None


def convex_quad(rng):
    while True:
        ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(4))
        rad = [rng.uniform(0.4, 1.0) for _ in range(4)]
        q = [(dy(r * math.cos(a)), dy(r * math.sin(a))) for a, r in zip(ang, rad)]
        if all((q[(k + 1) % 4][0] - q[k][0]) * (q[(k + 2) % 4][1] - q[(k + 1) % 4][1])
               - (q[(k + 1) % 4][1] - q[k][1]) * (q[(k + 2) % 4][0] - q[(k + 1) % 4][0]) > 0 for k in range(4)):
            return q


def offset3(rng, big):
    return [dy(rng.uniform(-1, 1) * big, 30) for _ in range(3)]


def fast_faces():
    rng = random.Random(4242)
    out = []
    for kind in ("a", "b"):
        for it in range(14):
            s = 10.0 ** rng.uniform(-6, 0)
            out.append((f"plane-{kind}{it}", place(rng, convex_quad(rng), dy(s), offset3(rng, rng.choice((0.0, 1.0, 1e3))),
                                                   kind)))
        for it in range(14):
            asp = 10.0 ** (1 + 4 * it / 13)                       # aspect 10 to 1e5
            h = dy(1.0 / asp)
            s1, s2 = dy(rng.uniform(0.0, 0.45)), dy(rng.uniform(0.0, 0.45))
            q = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(1) - s1, h), (s2, h)]
            s = 10.0 ** rng.uniform(-3, 0)
            out.append((f"thin-{kind}{it}-aspect{asp:.0e}", place(rng, q, dy(s), offset3(rng, rng.choice((0.0, 1e3))), kind)))
        for it in range(10):
            asp = 10.0 ** rng.uniform(0, 4)
            p = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (dy(rng.uniform(0.2, 0.8)), dy(1.0 / asp))]
            s = 10.0 ** rng.uniform(-3, 0)
            c = place(rng, [p[0], p[1], p[2], p[2]], dy(s), offset3(rng, 1e3 if it % 2 == 0 else 0.0), kind)
            if c is not None:
                c[3] = list(c[2])
            out.append((f"triangle-{kind}{it}" + ("-far" if it % 2 == 0 else ""), c))
    for it in range(8):
        x1, y1, x2, y2 = (rng.uniform(-1, 1) * 10.0 ** rng.uniform(-3, 3) for _ in range(4))
        z0, z1 = sorted((rng.uniform(-1, 1), rng.uniform(-1, 1)))
        out.append((f"extruded{it}", [[x1, y1, z0], [x2, y2, z0], [x2, y2, z1], [x1, y1, z1]]))
    for it in range(8):
        a = (dy(rng.uniform(0.5, 1)), dy(rng.uniform(-0.3, 0.3)))
        b = (dy(rng.uniform(-0.3, 0.3)), dy(10.0 ** rng.uniform(-4, 0)))
        q = [(Fr(0), Fr(0)), a, (a[0] + b[0], a[1] + b[1]), b]
        out.append((f"parallelogram-b{it}", place(rng, q, dy(10.0 ** rng.uniform(-3, 0)),
                                                   offset3(rng, rng.choice((0.0, 1e3))), "b")))
    for it in range(10):
        tw = 10.0 ** rng.uniform(-3, 0)
        s = 10.0 ** rng.uniform(-3, 0)
        o = [rng.uniform(-1, 1) * rng.choice((0.0, 1.0, 1e3)) for _ in range(3)]
        out.append((f"twisted{it}", [[o[0], o[1], o[2]], [o[0] + s, o[1], o[2]], [o[0] + s, o[1] + s, o[2] + tw * s],
                                     [o[0], o[1] + s, o[2]]]))
    for it in range(6):
        q = [(Fr(0), Fr(0)), (Fr(1), dy(rng.uniform(0.5, 1))), (dy(rng.uniform(0.1, 0.4)), dy(rng.uniform(0.1, 0.4))),
             (dy(rng.uniform(0.5, 1)), Fr(1))]
        out.append((f"dart-b{it}", place(rng, q, dy(10.0 ** rng.uniform(-3, 0)), offset3(rng, rng.choice((0.0, 1e3))), "b")))
    return [(lab, c) for lab, c in out if c is not None]


if __name__ == "__main__":
    lines = ["# Cases of the unit check BilinearDistance, written by test/unit/bilinear_distance_cases.py: D label p(3) c(3,4) "
             "d_exact; F label c(3,4) fast gap_exact"]
    nd = 0
    for label, p, c in distance_cases():
        d, info = R.ref_distance(p, c)
        lines.append(f"D {label} " + " ".join(repr(v) for v in p) + " " + " ".join(repr(v) for q in c for v in q)
                     + f" {format(d, '.25e')}")
        nd += 1
    nf = nskip = 0
    for label, c in fast_faces():
        g, outcome = gap_and_outcome(c)
        if outcome is None:
            nskip += 1
            continue
        lines.append(f"F {label} " + " ".join(repr(float(v)) for q in c for v in q) + (" T " if outcome else " F ")
                     + format(g, '.25e'))
        nf += 1
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"{nd} distance cases, {nf} fast-path faces ({nskip} left out: outcome not certain) written to {OUT}")
