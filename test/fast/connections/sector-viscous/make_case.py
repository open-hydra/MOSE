"""Inputs of the SectorViscous fast test, and the checks of its outputs.

A 3-D sector of an annulus closed in k by the rotational periodicity (type 200 on faces 5 and 6): beyond face 6 lies the
cell next to face 5 turned forward by the sector angle, beyond face 5 the cell next to face 6 turned back. The two periodic
faces are interior faces of the whole annulus, so their viscous flux (stress, heat flux, the work of the stress) must be that
of an interior face. x in [0, 0.2] (2 cells), r in [0.5, 1], slip walls (300) on the i and j faces, viscous air of
test/common/Air-viscous (mu and k constant) at T = 300 K and p = 10 Pa. The pressure does not enter the viscous fluxes the
test checks, only their round-off: the viscous force is the difference of two residuals that hold the pressure forces of the
faces, p A/V, about 6e7 times the viscous scale at 1 bar, and their round-off (1e-16 of them) would reach the tolerance; at
10 Pa they are about 6e3 times the scale, and the round-off stays within 1e-11 of it (as in ViscousLinearPatch). This is
the only reason for the low pressure: one explicit step, no stability limit involved.

ring  The sector from -45 to 45 degrees (S45, 4 cells in k) against the sector from -90 to 90 degrees (R90, 8 cells in k)
      that holds it: node planes at theta = m pi/8 (m integer, the same bits in both meshes), radii
      r_j = 0.5 + 0.5 (j/4)^1.3 (4 cells), and the cells of S45 are the cells of R90 with theta in [-45, 45] degrees, where
      the faces at +-45 degrees are interior faces. The velocity, at the cell centres, is periodic for a turn by 90 degrees
      and holds shear on the meridian planes:
        u = U (0, y, z) + G (0, -z, y)/r^2 + (4 C y z (y^2 - z^2), 0, 0),
      a radial field with uniform stress (normal stress (2/3) mu U on a meridian plane), a potential vortex u_theta = G/r
      (tau_rtheta = -2 mu G/r^2) and an axial harmonic field u_x = C r^4 sin(4 theta) (tau_xtheta = -4 mu C r^3 on the
      planes at +-45 degrees); each of them has zero viscous force. With diffusive-metric = inverse-mean, and with
      stencil, whose metric of the periodic face takes the differences of the positions of the cell and of its image row
      as the interior face takes those of its neighbours, the two meshes compute the same face. The viscous force of a
      cell is r/dt of a Navier-Stokes step minus r/dt of an Euler step from the same state; in every cell of S45, for
      rho u, rho v, rho w and rho E, it must equal that of the same cell of R90 within 1e-9 of the scale
        S = mu (U + 2 G/r0^2 + 4 C r1^3)/(r0 dtheta)   (momentum),   S (U r1 + G/r0 + C r1^4)   (energy),
      r0 = 0.5, r1 = 1, dtheta = pi/8, with each metric. (Before the periodic face took inverse-mean with every value of
      the key, and with stencil S45 differed from R90 on these stretched radii.)
mirror The sector from -45 to 45 degrees in 2 x 8 x 8 cells (uniform radii, node planes at theta = m pi/16) with the
      radial field alone, k faces of type 200 against type 300: the radial field is its own mirror image across the
      meridian planes, so the two conditions set the same ghost cells, and the viscous stress on the planes is the
      normal stress (2/3) mu U. The default face metric. A Navier-Stokes step with each: in the cells not next to the j
      faces, r/dt of rho u, rho v and rho w must agree within 0.02 (2/3) mu U A/V, with A/V = 1/(r dtheta) the area of a k
      face over the volume (r the mid radius of the cell), and r/dt of rho E within the same times U r1. (A field with shear
      on the plane is not its own mirror image, and the symmetry condition would not set the same problem.)
shift  A sector with no angle: a block of 4 x 4 x 8 cubic cells of 1 mm periodic in z, closed in k by type 200 (shift200:
      its last k node plane is the first one translated by the period T = (0, 0, Lz)) or by periodic pairs 201 (shift201),
      the i and j faces of type 0 (no flux). u = U0 sin(2 pi z/Lz) along x, rho nu~ = 3 mu uniform, no wall;
      Navier-Stokes with Spalart-Allmaras, explicit sources (point-implicit = false), the default face metric (stencil).
      Both conditions put behind each k face the values of the cell of the same column at the other end, and the ghost
      cell sits where that cell lies translated by the period, -T behind face 5 and +T behind face 6: the two steps must
      be the same. In every cell r/dt of rho, rho u, rho v, rho w, rho E and rho nu~ must agree within 1e-9 of the largest
      |r/dt| of that variable over the cells of shift201. A ghost cell left at the cell at the other end, not translated,
      scales the derivative in z of the cells next to the k faces by 2/(2 - NK) = -1/3: in the cell metric of the source
      terms (the vorticity of Spalart-Allmaras, so rho nu~ in those cells) and in the k column of the face metric of the
      i faces (the stress tau_zx = mu du/dz, which the faces of type 0 do not balance in the cells next to them, so rho w).
stretch The mirror sector with radii stretched as the ring, r_j = 0.5 + 0.5 (j/8)^1.3 (stretchpot), or growing by 1.2
      from one cell to the next (stretchgeo), k faces of type 200, the default face metric (stencil), the radial field
      alone: a linear field that the rotation allows (u(R x) = R u(x)), so the viscous force of every cell is zero. In the
      cells not next to the j faces (slip walls 300, which the radial field does not satisfy) the force must vanish to
      round-off, below 1e-11 of (2/3) mu U A/V, also next to the periodic faces. (Before, with the inverse-mean metric of
      the periodic face, whose columns along the face are the mean edges of the cell and of its image instead of the
      differences of the centres, those cells had 2.6e-2 and 5.9e-3 of it.)
wall  A translated sector next to a wall: a block of 4 x 4 x 6 cells periodic along z (period 6 mm), closed in k by type 200,
      whose k node planes are tilted, z = z_k + 0.3 x, so that the normal of the periodic faces has a part along x; cells
      growing by 1.3 along x from a wall 301 on i = 1 (the plane x = 0), j lines bent as y = y_j + 50 x^2 (1/m), 400 on the
      other i face and on the j faces. Behind the wall a ghost cell holds the velocity at 2 x_f - x_1 and the temperature
      at x_1 + s (x_1 - x_2) (extrapolated, Ghost_Wall), two positions that the bent lines keep apart: the periodic faces of
      the cells next to the wall take one metric per group of quantities. Two linear fields that the periodicity allows:
        wallmom  u = x (a, b, c), zero on the wall, uniform T and p;   wallene  u = 0, T = T0 + g x, uniform p, the heat
                 flux of the field on the wall (q into the fluid, k g . n with the normal of the face);
      a Navier-Stokes step, and a reference step with mu = k = 1e-30 and no heat flux on the wall (the solver takes no wall
      with the Euler equations): the difference is the viscous force (wallmom) and the heat conduction (wallene), zero for
      these fields. In every cell with no face 400, the wall cells included, it must vanish below 1e-11 of
      mu |G| (sum of the face areas)/volume (wallmom, rho u, rho v, rho w) and k |g| (sum of the face areas)/volume
      (wallene, rho E). The generator prints the largest distance of the two positions of the wall ghost cells.

One explicit Euler step with local time stepping, first-order reconstruction. OUTPUT/diagnostic.tec holds, per cell, the
increment of each conserved variable over the step, and the step dt.

Usage (from the test directory):
  python3 -B make_case.py case ring45|ring90|mirror200|mirror300   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case stretchpot|stretchgeo               INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case shift200|shift201                   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case wallmom|wallene                     INPUT/ic.tec, INPUT/bc-ns.txt, INPUT/bc-ref.txt and
                                                                   INPUT/transport-ref.dat
  python3 -B make_case.py ring NS45 EU45 NS90 EU90                 largest difference of the viscous forces over the scale
  python3 -B make_case.py mirror F200 F300                         largest difference over (2/3) mu U A/V
  python3 -B make_case.py linear NAME NS EU                        largest viscous force over (2/3) mu U A/V (stretch)
  python3 -B make_case.py shift F200 F201                          largest difference over the largest value of F201
  python3 -B make_case.py wall NAME NS REF                         largest viscous force or conduction over its scale
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
X1, R0, R1 = 0.2, 0.5, 1.0
NI = 2
P0, T0 = 10.0, 300.0
RHO0 = P0 / (RGAS * T0)
U, G, C = 50.0, 14.0, 30.0
TOL_RING, TOL_MIRROR = 1.0e-9, 0.02
SHIFT_N, SHIFT_H, SHIFT_U = (4, 4, 8), 1.0e-3, 10.0    # shift: cells in i, j and k, cell size (m), U0 (m/s)
SHIFT_M = 3.0 * MU                                     # shift: rho nu~
TOL_SHIFT = 1.0e-9
TOL_LINEAR = 1.0e-11
KCOND = 0.016544513989              # test/common/Air-viscous
WALL_N, WALL_H, WALL_Q = (4, 4, 6), 1.0e-3, 1.3        # wall: cells in i, j, k, base cell size (m), growth along x
WALL_BEND, WALL_TILT = 50.0, 0.3                       # wall: y = y_j + WALL_BEND x^2, z = z_k + WALL_TILT x
WALL_G = (1.0e3, 0.6e3, -0.8e3)                        # wallmom: u = x WALL_G, 1/s
WALL_GT = 2.5e3                                        # wallene: T = T0 + WALL_GT x, K/m
TOL_WALL = 1.0e-11


def fmt(v):
    return f"{v:.16E}"


def mesh(name):
    """Radii of the node rows, node planes m (theta = m pi/den) and the field of the case."""
    if name in ("ring45", "ring90"):
        nj = 4
        rs = [R0 + (R1 - R0) * (j / nj) ** 1.3 for j in range(nj + 1)]
        ms = range(-2, 3) if name == "ring45" else range(-4, 5)
        return rs, list(ms), 8, True
    nj = 8
    if name == "stretchpot":
        rs = [R0 + (R1 - R0) * (j / nj) ** 1.3 for j in range(nj + 1)]
    elif name == "stretchgeo":
        w = [1.2 ** j for j in range(nj)]
        rs = [R0]
        for x in w:
            rs.append(rs[-1] + (R1 - R0) * x / sum(w))
    else:
        rs = [R0 + (R1 - R0) * j / nj for j in range(nj + 1)]
    return rs, list(range(-4, 5)), 16, False


def velocity(y, z, full):
    r2 = y * y + z * z
    if not full:
        return (0.0, U * y, U * z)
    return (4.0 * C * y * z * (y * y - z * z), U * y - G * z / r2, U * z + G * y / r2)


def write_case(name):
    rs, ms, den, full = mesh(name)
    nj, nk = len(rs) - 1, len(ms) - 1
    xyz = {}
    for k, m in enumerate(ms):
        th = m * math.pi / den
        for j, r in enumerate(rs):
            for i in range(NI + 1):
                xyz[(i, j, k)] = tuple(float(fmt(v)) for v in (X1 * i / NI, r * math.cos(th), r * math.sin(th)))
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    cells = []
    for k in range(1, nk + 1):
        for j in range(1, nj + 1):
            for i in range(1, NI + 1):
                c = [sum(xyz[(i - 1 + a, j - 1 + b, k - 1 + d)][n] for a in (0, 1) for b in (0, 1) for d in (0, 1)) / 8.0
                     for n in range(3)]
                cells.append((RHO0,) + velocity(c[1], c[2], full) + (P0,))
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for n in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(NI + 1):
                        f.write(f" {fmt(xyz[(i, j, k)][n])}\n")
        for v in range(5):
            for c in cells:
                f.write(f" {fmt(c[v])}\n")

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    tk = 300 if name == "mirror300" else 200
    out = []                                            # face by face
    for face, i in ((1, 1), (2, NI)):
        out += [rec(i, j, k, face, 300) for k in range(1, nk + 1) for j in range(1, nj + 1)]
    for face, j in ((3, 1), (4, nj)):
        out += [rec(i, j, k, face, 300) for k in range(1, nk + 1) for i in range(1, NI + 1)]
    for face, k in ((5, 1), (6, nk)):
        out += [rec(i, j, k, face, tk) for j in range(1, nj + 1) for i in range(1, NI + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def write_shift(kind):
    """INPUT/ic.tec and INPUT/bc.txt of shift200 (kind "200") or shift201 (kind "201")."""
    ni, nj, nk = SHIFT_N
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for n in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fmt((i * SHIFT_H, j * SHIFT_H, k * SHIFT_H)[n])}\n")
        nc = ni * nj * nk
        u = [SHIFT_U * math.sin(2.0 * math.pi * (k + 0.5) / nk) for k in range(nk) for _ in range(nj * ni)]
        for col in ([RHO0] * nc, u, [0.0] * nc, [0.0] * nc, [P0] * nc, [SHIFT_M] * nc):
            f.write("".join(f" {fmt(v)}\n" for v in col))

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    out = []                                            # face by face
    for face, i in ((1, 1), (2, ni)):
        out += [rec(i, j, k, face, 0) for k in range(1, nk + 1) for j in range(1, nj + 1)]
    for face, j in ((3, 1), (4, nj)):
        out += [rec(i, j, k, face, 0) for k in range(1, nk + 1) for i in range(1, ni + 1)]
    for face, k, ks, fs in ((5, 1, nk, 6), (6, nk, 1, 5)):
        for j in range(1, nj + 1):
            for i in range(1, ni + 1):
                if kind == "200":
                    out.append(rec(i, j, k, face, 200))
                else:                                   # source: the cell of the same column at the other end
                    out += [rec(i, j, k, face, 201),
                            '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, i, j, ks, fs, 1, 0, 0, 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the numbers of nodes."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, (ni - 1, nj - 1, nk - 1)


CONS = ("rhou", "rhov", "rhow", "rhoe")


def per_step(path):
    z, n = read_zone(path)
    return {q: [a / b for a, b in zip(z[q], z["dt"])] for q in CONS}, n


def ring(ns45, eu45, ns90, eu90):
    (a45, n45), (b45, _) = per_step(ns45), per_step(eu45)
    (a90, n90), (b90, _) = per_step(ns90), per_step(eu90)
    ni, nj, nk = n45
    sm = MU * (U + 2.0 * G / R0 ** 2 + 4.0 * C * R1 ** 3) / (R0 * math.pi / 8.0)
    se = sm * (U * R1 + G / R0 + C * R1 ** 4)
    worst = {"momentum": (0.0, None), "energy": (0.0, None)}
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                c45 = (k * nj + j) * ni + i
                c90 = ((k + 2) * nj + j) * ni + i          # the same cell in R90, two cells further in k
                for q in CONS:
                    d = abs((a45[q][c45] - b45[q][c45]) - (a90[q][c90] - b90[q][c90]))
                    key, s = ("energy", se) if q == "rhoe" else ("momentum", sm)
                    if d / s >= worst[key][0]:
                        worst[key] = (d / s, f"{q} of cell ({i + 1}, {j + 1}, {k + 1})")
    ok = True
    for key in ("momentum", "energy"):
        ratio, where = worst[key]
        ok &= ratio <= TOL_RING
        print(f"    {key}: largest |f(S45) - f(R90)| over the scale {ratio:.3e} ({where}); tolerance {TOL_RING:.0e}")
    sys.exit(0 if ok else 1)


def mirror(f200, f300):
    (a, n), (b, _) = per_step(f200), per_step(f300)
    ni, nj, nk = n
    worst = {"momentum": (0.0, None), "energy": (0.0, None)}
    dth = math.pi / 16.0
    for k in range(nk):
        th = (k - nk / 2 + 0.5) * dth
        for j in range(1, nj - 1):                     # not next to the j faces
            for i in range(ni):
                c = (k * nj + j) * ni + i
                r = R0 + (R1 - R0) * (j + 0.5) / nj         # mid radius of the cell
                s = (2.0 / 3.0) * MU * U / (r * dth)
                for q in CONS:
                    d = abs(a[q][c] - b[q][c])
                    key, sc = ("energy", s * U * R1) if q == "rhoe" else ("momentum", s)
                    if d / sc >= worst[key][0]:
                        worst[key] = (d / sc, f"{q} of cell ({i + 1}, {j + 1}, {k + 1}), theta {math.degrees(th):.2f} deg")
    ok = True
    for key in ("momentum", "energy"):
        ratio, where = worst[key]
        ok &= ratio <= TOL_MIRROR
        print(f"    {key}: largest |f(200) - f(300)| over (2/3) mu U A/V {ratio:.3e} ({where}); tolerance {TOL_MIRROR}")
    sys.exit(0 if ok else 1)


def shift(f200, f201):
    (a, n), (b, _) = read_zone(f200), read_zone(f201)
    ni, nj, _ = n
    ok = True
    for q in ("rho", "rhou", "rhov", "rhow", "rhoe", "rhoNuTilde"):
        fa = [x / t for x, t in zip(a[q], a["dt"])]
        fb = [x / t for x, t in zip(b[q], b["dt"])]
        scale = max(abs(x) for x in fb)
        d, c = max((abs(x - y), c) for c, (x, y) in enumerate(zip(fa, fb)))
        ratio = d / scale if scale > 0.0 else d
        ok &= ratio <= TOL_SHIFT
        where = f"cell ({c % ni + 1}, {c // ni % nj + 1}, {c // (ni * nj) + 1})"
        print(f"    {q}: largest |f(200) - f(201)| over the largest |f(201)| {ratio:.3e} ({where}); tolerance {TOL_SHIFT:.0e}")
    sys.exit(0 if ok else 1)


def linear(name, ns, eu):
    """stretch: the viscous force of the radial field (r/dt of NS minus r/dt of Euler) in the cells not next to the j faces,
    over (2/3) mu U A/V, A/V = 1/(r dtheta) with r the mid radius of the cell; next to the periodic faces too."""
    rs, _, _, _ = mesh(name)
    (a, n), (b, _) = per_step(ns), per_step(eu)
    ni, nj, nk = n
    dth = math.pi / 16.0
    worst, where = 0.0, None
    for k in range(nk):
        for j in range(1, nj - 1):
            for i in range(ni):
                c = (k * nj + j) * ni + i
                s = (2.0 / 3.0) * MU * U / (0.5 * (rs[j] + rs[j + 1]) * dth)
                for q in ("rhou", "rhov", "rhow"):
                    e = abs(a[q][c] - b[q][c]) / s
                    if e >= worst:
                        worst, where = e, f"{q} of cell ({i + 1}, {j + 1}, {k + 1})"
    print(f"    {name}: largest |viscous force| over (2/3) mu U A/V {worst:.3e} ({where}); tolerance {TOL_LINEAR:.0e}")
    sys.exit(0 if worst <= TOL_LINEAR else 1)


def wall_nodes():
    """Nodes of the wall case, X[i][j][k] = (x, y, z): x growing by WALL_Q from the wall x = 0, y lines bent, k planes
    tilted; the last k plane is the first one translated by the period (0, 0, NK h)."""
    ni, nj, nk = WALL_N
    xs = [WALL_H * (WALL_Q ** i - 1.0) / (WALL_Q - 1.0) for i in range(ni + 1)]
    return [[[tuple(float(fmt(v)) for v in (xs[i], j * WALL_H + WALL_BEND * xs[i] ** 2, k * WALL_H + WALL_TILT * xs[i]))
              for k in range(nk + 1)] for j in range(nj + 1)] for i in range(ni + 1)]


def vsub(a, b):
    return tuple(a[m] - b[m] for m in range(3))


def vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def vdot(a, b):
    return sum(a[m] * b[m] for m in range(3))


def cell_faces(X, i, j, k):
    """The six faces of cell (i, j, k) (0-based): (nodes, outward area vector), the area vector half the cross product of
    the diagonals (Compute_Norm_Area); face index 1..6 as MOSE (1, 2 along i, 3, 4 along j, 5, 6 along k)."""
    n = {(a, b, c): X[i + a][j + b][k + c] for a in (0, 1) for b in (0, 1) for c in (0, 1)}
    centre = tuple(sum(p[m] for p in n.values()) / 8.0 for m in range(3))
    quads = {1: [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)], 2: [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
             3: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)], 4: [(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
             5: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], 6: [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]}
    out = {}
    for f, q in quads.items():
        p = [n[c] for c in q]
        av = tuple(0.5 * c for c in vcross(vsub(p[2], p[0]), vsub(p[3], p[1])))
        xf = tuple(sum(c[m] for c in p) / 4.0 for m in range(3))
        if vdot(vsub(xf, centre), av) < 0.0:
            av = tuple(-c for c in av)
        out[f] = (p, av, xf)
    return out, centre


def cell_scale(X, i, j, k):
    """(sum of the face areas)/volume of cell (i, j, k), the volume (1/3) sum over the faces of centre . area vector."""
    faces, _ = cell_faces(X, i, j, k)
    total = sum(math.sqrt(vdot(av, av)) for _, av, _ in faces.values())
    vol = sum(vdot(xf, av) for _, av, xf in faces.values()) / 3.0
    return total / abs(vol)


def wall_state(run, x):
    if run == "wallmom":
        return RHO0, tuple(g * x[0] for g in WALL_G), P0
    t = T0 + WALL_GT * x[0]
    return P0 / (RGAS * t), (0.0, 0.0, 0.0), P0


def write_wall(run):
    """INPUT/ic.tec, INPUT/bc-ns.txt, INPUT/bc-ref.txt and INPUT/transport-ref.dat of the wall case (run wallmom or
    wallene); prints the largest distance of the two positions of the ghost cells of the wall, over the first cell."""
    ni, nj, nk = WALL_N
    X = wall_nodes()
    for p in ("INPUT/ic.tec", "INPUT/bc-ns.txt", "INPUT/bc-ref.txt", "INPUT/transport-ref.dat"):
        if os.path.lexists(p):
            os.remove(p)
    cells = []
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                _, c = cell_faces(X, i, j, k)
                rho, u, p = wall_state(run, c)
                cells.append((rho,) + u + (p,))
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for n in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fmt(X[i][j][k][n])}\n")
        for v in range(5):
            for c in cells:
                f.write(f" {fmt(c[v])}\n")

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    worst = 0.0
    for tag in ("ns", "ref"):
        out = []
        for k in range(1, nk + 1):
            for j in range(1, nj + 1):
                faces, x1 = cell_faces(X, 0, j - 1, k - 1)
                _, x2 = cell_faces(X, 1, j - 1, k - 1)
                p, av, xf = faces[1]
                q = 0.0
                if run == "wallene" and tag == "ns":
                    nrm = tuple(-c / math.sqrt(vdot(av, av)) for c in av)   # face 1: unit normal towards growing i
                    q = -KCOND * WALL_GT * nrm[0]                            # q into the fluid: k g . n on face 1 (odd)
                out.append(rec(1, j, k, 1, 301) + f"    {q:.17E},    0.000000E+00,    0.000000E+00,\n")
                if tag == "ns":
                    # the two positions of the ghost cell (Ghost_Position): 2 x_f - x_1 and x_1 + s (x_1 - x_2), with
                    # s = 2 l1/(l1 + l2) from the mean edges of the first two cells across the face (Extrapolation_Ratio)
                    lines = [(jj, kk) for jj in (j - 1, j) for kk in (k - 1, k)]
                    m1 = tuple(sum(X[1][jj][kk][m] - X[0][jj][kk][m] for jj, kk in lines) / 4.0 for m in range(3))
                    m2 = tuple(sum(X[2][jj][kk][m] - X[1][jj][kk][m] for jj, kk in lines) / 4.0 for m in range(3))
                    l1, l2 = math.sqrt(vdot(m1, m1)), math.sqrt(vdot(m2, m2))
                    s = 2.0 * l1 / (l1 + l2)
                    xm = tuple(2.0 * xf[m] - x1[m] for m in range(3))
                    xe = tuple(x1[m] + s * (x1[m] - x2[m]) for m in range(3))
                    worst = max(worst, math.sqrt(vdot(vsub(xm, xe), vsub(xm, xe))) / l1)
            for j in range(1, nj + 1):
                out.append(rec(ni, j, k, 2, 400))
            for i in range(1, ni + 1):
                out += [rec(i, 1, k, 3, 400), rec(i, nj, k, 4, 400)]
        for face, k in ((5, 1), (6, nk)):
            out += [rec(i, j, k, face, 200) for j in range(1, nj + 1) for i in range(1, ni + 1)]
        with open(f"INPUT/bc-{tag}.txt", "w") as f:
            f.write("".join(out))
    # the reference: mu = k = 1e-30 (the solver takes no wall with the Euler equations)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "..", "..", "common", "Air-viscous", "transport.dat")) as src, \
            open("INPUT/transport-ref.dat", "w") as dst:
        for line in src:
            p = line.split()
            dst.write(f"{p[0]} 1.0E-30 1.0E-30\n" if len(p) == 3 and re.match(r"^[0-9.Ee+-]+$", p[0]) else line)
    print(f"    {run}: largest distance of the two positions of the wall ghost cells over the first cell {worst:.3e}")


def wall(name, ns, ref):
    """wall: the viscous force (wallmom, rho u, rho v, rho w) or the conduction (wallene, rho E) of the cells with no face
    400, r/dt of NS minus r/dt of the reference, over mu |G| (sum of the face areas)/volume or k |g| (...)/volume."""
    (a, n), (b, _) = per_step(ns), per_step(ref)
    ni, nj, nk = n
    X = wall_nodes()
    gnorm = math.sqrt(sum(g * g for g in WALL_G))
    worst, where = 0.0, None
    for k in range(nk):
        for j in range(1, nj - 1):
            for i in range(ni - 1):
                c = (k * nj + j) * ni + i
                sc = cell_scale(X, i, j, k) * (MU * gnorm if name == "wallmom" else KCOND * WALL_GT)
                for q in (("rhou", "rhov", "rhow") if name == "wallmom" else ("rhoe",)):
                    e = abs(a[q][c] - b[q][c]) / sc
                    if e >= worst:
                        worst, where = e, f"{q} of cell ({i + 1}, {j + 1}, {k + 1})"
    what = "viscous force" if name == "wallmom" else "heat conduction"
    print(f"    {name}: largest |{what}| over its scale {worst:.3e} ({where}), cells with no face 400, the wall cells"
          f" included; tolerance {TOL_WALL:.0e}")
    sys.exit(0 if worst <= TOL_WALL else 1)


def main():
    what = sys.argv[1]
    if what == "case" and sys.argv[2].startswith("shift"):
        write_shift(sys.argv[2][len("shift"):])
    elif what == "case" and sys.argv[2].startswith("wall"):
        write_wall(sys.argv[2])
    elif what == "case":
        write_case(sys.argv[2])
    elif what == "ring":
        ring(*sys.argv[2:6])
    elif what == "mirror":
        mirror(*sys.argv[2:4])
    elif what == "linear":
        linear(*sys.argv[2:5])
    elif what == "shift":
        shift(*sys.argv[2:4])
    elif what == "wall":
        wall(*sys.argv[2:5])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
