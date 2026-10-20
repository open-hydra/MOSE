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
split The rule of type 200: it joins the two k faces of one block. Beyond face 5 lies the cell next to face 6 of the same
      block, in the same column (i, j), turned back by the angle, and beyond face 6 the cell next to face 5 turned forward;
      the record names no other block. So every record of type 200 must lie on face 5 or 6, inside its block, with a record
      of type 200 on the other k face of the block in its column, and any other layout stops at set-up. A sector from
      -22.5 to 22.5 degrees (x in [0, 0.2], r in [0.5, 1], 2 x 4 x 8 cells, faces i and j of type 300) of air at 1 bar and
      300 K at rest with a rotation periodic with the sector, u_theta = W (1 + cos(8 (theta + 22.5 degrees))/2), W = 50 m/s,
      one explicit Euler step:
        split-whole    one block, type 200 on faces 5 and 6;
        split-cutj     two blocks of 2 cells in j joined by 101, type 200 on the k faces of both: it runs, and in every cell
                       its step must be that of split-whole within 1e-12 of the scale (the largest |r/dt| of split-whole: of
                       rho, of rho E, and one for rho u, rho v and rho w);
        split-cutk     two blocks of 4 cells in k joined by 101, type 200 on face 5 of the first and on face 6 of the second
                       (a cut along k): in the column of record 49, on face 5, face 6 of block 1 is of type 101;
        split-face4    split-whole with type 200 on face 4 (j = 4): record 81 is not on a k face;
        split-outside  split-whole with record 97 (face 5, cell (1, 1, 1)) written at i = 3: outside the block;
        split-missing  split-whole with record 97 written at the cell (2, 1, 1): no record on face 5 in the column (1, 1)
                       of record 105, on face 6;
      a 2-D axisymmetric wedge of 1 degree (one cell in k, x in [0, 0.2], r in [0.5, 1], 4 x 4 cells, faces 1 to 3 of type
      300) of the same air with a radial velocity of 20 m/s and a swirl of 50 m/s:
        split-wedge4   type 200 on faces 4, 5 and 6: record 13 is not on a k face;
        split-wedge6   type 200 on face 5, type 300 on faces 4 and 6: in the column of record 17 face 6 is of type 300.
      The last six must stop at set-up, before the type of the mesh is read and before any step, with exit code 1 and the
      message that names the record (run.sh). Without the check, the sector cut along k runs with the rows of cells next to
      the cut and to its k faces wrong by up to 31 times the scale of the step, the wedge with type 200 on face 4 or on one
      k face only runs, with fluxes there that are neither those of the rotational periodicity nor those of a wall, and the
      sector with type 200 on face 4 stops at its first flux.
ring101 The remedy that the message of split gives for a ring of several blocks: the whole turn as four blocks of 90
      degrees, block 1 the cells and the field of ring45 and block b those of the sector turned by 90 (b - 1) degrees, the
      node plane at 315 degrees written with the bits of the plane at -45 degrees (as a mesher writes coincident planes),
      face 6 of each block joined by 101 to face 5 of the next one, the junction at theta = 2 pi included, no type 200. The
      field of ring is periodic for a turn by 90 degrees, so in every cell block 1 must give the step of ring45, which type
      200 closes: an Euler step within 1e-12 of the largest |r/dt| of ring45 of each variable, and the viscous force
      (diffusive-metric = stencil) within 1e-9 of the scale of ring.

One explicit Euler step with local time stepping, first-order reconstruction. OUTPUT/diagnostic.tec holds, per cell, the
increment of each conserved variable over the step, and the step dt.

Usage (from the test directory):
  python3 -B make_case.py case ring45|ring90|mirror200|mirror300   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case stretchpot|stretchgeo               INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case shift200|shift201                   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py case wallmom|wallene                     INPUT/ic.tec, INPUT/bc-ns.txt, INPUT/bc-ref.txt and
                                                                   INPUT/transport-ref.dat
  python3 -B make_case.py case split-NAME|ring101                  INPUT/ic.tec and INPUT/bc.txt (NAME: whole, cutj, cutk,
                                                                   face4, outside, missing, wedge4, wedge6)
  python3 -B make_case.py ring NS45 EU45 NS90 EU90                 largest difference of the viscous forces over the scale
  python3 -B make_case.py mirror F200 F300                         largest difference over (2/3) mu U A/V
  python3 -B make_case.py linear NAME NS EU                        largest viscous force over (2/3) mu U A/V (stretch)
  python3 -B make_case.py shift F200 F201                          largest difference over the largest value of F201
  python3 -B make_case.py wall NAME NS REF                         largest viscous force or conduction over its scale
  python3 -B make_case.py split WHOLE CUTJ                         largest difference of the steps over the scale
  python3 -B make_case.py ring101 NS45 EU45 NS101 EU101            the same for the steps and the viscous forces of block 1
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
SPLIT_N, SPLIT_P, SPLIT_W = (2, 4, 8), 1.0e5, 50.0     # split: cells in i, j and k, pressure (Pa), rotation (m/s)
SPLIT_TH0, SPLIT_TH = -math.pi / 8.0, math.pi / 4.0    # split: first k node plane and angle of the sector
SPLIT_RHO = SPLIT_P / (RGAS * T0)
WEDGE_N, WEDGE_DELTA, WEDGE_V = (4, 4), math.radians(1.0), 20.0   # split: wedge cells in i and j, angle, radial velocity
TOL_SPLIT = 1.0e-12


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


def split_node(i, j, k):
    """Node (i, j, k) of the sector of split: uniform radii, theta from SPLIT_TH0 over SPLIT_TH in SPLIT_N[2] cells."""
    ni, nj, nk = SPLIT_N
    th = SPLIT_TH0 + SPLIT_TH * k / nk
    r = R0 + (R1 - R0) * j / nj
    return tuple(float(fmt(v)) for v in (X1 * i / ni, r * math.cos(th), r * math.sin(th)))


def wedge_node(i, j, k):
    """Node (i, j, k) of the wedge of split: the two k node planes at theta = -delta/2 and +delta/2."""
    wi, wj = WEDGE_N
    th = (k - 0.5) * WEDGE_DELTA
    r = R0 + (R1 - R0) * j / wj
    return tuple(float(fmt(v)) for v in (X1 * i / wi, r * math.cos(th), r * math.sin(th)))


def split_state(c):
    """Sector of split: air at rest with a rotation periodic with the sector, at the cell centre c."""
    th = math.atan2(c[2], c[1])
    ut = SPLIT_W * (1.0 + 0.5 * math.cos(2.0 * math.pi * (th - SPLIT_TH0) / SPLIT_TH))
    return (SPLIT_RHO, 0.0, -ut * math.sin(th), ut * math.cos(th), SPLIT_P)


def wedge_state(c):
    """Wedge of split: a radial velocity and a swirl (the cell centres lie at theta = 0)."""
    return (SPLIT_RHO, 0.0, WEDGE_V, SPLIT_W, SPLIT_P)


def split_zone(f, b, node, state, lo, hi):
    """The zone of block b of INPUT/ic.tec: the nodes from lo to hi (node indices i, j, k of the whole mesh), and the
    states of its cells at their centres."""
    (i0, j0, k0), (i1, j1, k1) = lo, hi
    f.write(f" ZONE  T = B{b}-IG, I={i1 - i0 + 1}, J={j1 - j0 + 1}, K={k1 - k0 + 1}, DATAPACKING=BLOCK, "
            "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
    for n in range(3):
        for k in range(k0, k1 + 1):
            for j in range(j0, j1 + 1):
                for i in range(i0, i1 + 1):
                    f.write(f" {fmt(node(i, j, k)[n])}\n")
    cells = []
    for k in range(k0 + 1, k1 + 1):
        for j in range(j0 + 1, j1 + 1):
            for i in range(i0 + 1, i1 + 1):
                cells.append(state([sum(node(i - 1 + a, j - 1 + bb, k - 1 + d)[n] for a in (0, 1) for bb in (0, 1)
                                        for d in (0, 1)) / 8.0 for n in range(3)]))
    for v in range(5):
        f.write("".join(f" {fmt(c[v])}\n" for c in cells))


def block_rec(b, i, j, k, face, t, src=None):
    """A record of bc.txt, with the property line of a connection 101: its source block, cell, face and the identity."""
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, face, t)
    if src:
        s += '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (src + (1, 0, 0, 1))
    return s


def box_faces(b, ni, nj, nk, types):
    """The records of the six faces of block b of ni x nj x nk cells, face by face; types[f - 1] the type of face f (None:
    no record, a connection written apart)."""
    out = []
    for face, i in ((1, 1), (2, ni)):
        if types[face - 1] is not None:
            out += [block_rec(b, i, j, k, face, types[face - 1]) for k in range(1, nk + 1) for j in range(1, nj + 1)]
    for face, j in ((3, 1), (4, nj)):
        if types[face - 1] is not None:
            out += [block_rec(b, i, j, k, face, types[face - 1]) for k in range(1, nk + 1) for i in range(1, ni + 1)]
    for face, k in ((5, 1), (6, nk)):
        if types[face - 1] is not None:
            out += [block_rec(b, i, j, k, face, types[face - 1]) for j in range(1, nj + 1) for i in range(1, ni + 1)]
    return out


def write_split(name):
    """INPUT/ic.tec and INPUT/bc.txt of the case split-<name>."""
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    ni, nj, nk = SPLIT_N
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        if name in ("whole", "face4", "outside", "missing"):
            split_zone(f, 1, split_node, split_state, (0, 0, 0), (ni, nj, nk))
            out = box_faces(1, ni, nj, nk, [300, 300, 300, 200 if name == "face4" else 300, 200, 200])
            if name in ("outside", "missing"):
                # record 97, the first one of face 5 (cell (1, 1, 1)), written at i = 3 or at the cell (2, 1, 1)
                out[96] = block_rec(1, 3 if name == "outside" else 2, 1, 1, 5, 200)
        elif name == "cutk":
            h = nk // 2
            split_zone(f, 1, split_node, split_state, (0, 0, 0), (ni, nj, h))
            split_zone(f, 2, split_node, split_state, (0, 0, h), (ni, nj, nk))
            out = (box_faces(1, ni, nj, h, [300, 300, 300, 300, 200, None]) +
                   box_faces(2, ni, nj, h, [300, 300, 300, 300, None, 200]))
            for j in range(1, nj + 1):
                for i in range(1, ni + 1):
                    out += [block_rec(1, i, j, h, 6, 101, (2, i, j, 1, 5)), block_rec(2, i, j, 1, 5, 101, (1, i, j, h, 6))]
        elif name == "cutj":
            h = nj // 2
            split_zone(f, 1, split_node, split_state, (0, 0, 0), (ni, h, nk))
            split_zone(f, 2, split_node, split_state, (0, h, 0), (ni, nj, nk))
            out = (box_faces(1, ni, h, nk, [300, 300, 300, None, 200, 200]) +
                   box_faces(2, ni, h, nk, [300, 300, None, 300, 200, 200]))
            for k in range(1, nk + 1):
                for i in range(1, ni + 1):
                    out += [block_rec(1, i, h, k, 4, 101, (2, i, 1, k, 3)), block_rec(2, i, 1, k, 3, 101, (1, i, h, k, 4))]
        elif name in ("wedge4", "wedge6"):
            wi, wj = WEDGE_N
            split_zone(f, 1, wedge_node, wedge_state, (0, 0, 0), (wi, wj, 1))
            out = box_faces(1, wi, wj, 1, [300, 300, 300, 200, 200, 200] if name == "wedge4" else
                            [300, 300, 300, 300, 200, 300])
        else:
            sys.exit(f"unknown case split-{name}")
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def write_ring101():
    """INPUT/ic.tec and INPUT/bc.txt of ring101: four blocks of the nodes, cells and field of ring45, block b turned by
    90 (b - 1) degrees, joined by 101 face 6 to face 5 of the next block, block 4 to block 1."""
    rs, _, den, _ = mesh("ring45")
    nj, nk, nb = len(rs) - 1, 4, 4

    def node(i, j, m):
        mm = (m + 2) % 16 - 2                               # m = 14, at 315 degrees, is m = -2: the same bits
        th = mm * math.pi / den
        return tuple(float(fmt(v)) for v in (X1 * i / NI, rs[j] * math.cos(th), rs[j] * math.sin(th)))

    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b in range(1, nb + 1):
            m0 = 4 * (b - 1) - 2
            split_zone(f, b, lambda i, j, k: node(i, j, m0 + k), lambda c: (RHO0,) + velocity(c[1], c[2], True) + (P0,),
                       (0, 0, 0), (NI, nj, nk))
    out = []
    for b in range(1, nb + 1):
        nxt, prv = b % nb + 1, (b - 2) % nb + 1
        out += box_faces(b, NI, nj, nk, [300, 300, 300, 300, None, None])
        for face, k, src in ((5, 1, (prv, nk, 6)), (6, nk, (nxt, 1, 5))):
            out += [block_rec(b, i, j, k, face, 101, (src[0], i, j, src[1], src[2]))
                    for j in range(1, nj + 1) for i in range(1, NI + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zones(path):
    """Cell values of a MOSE ASCII file of one or more zones (BLOCK packing), zone after zone: a list of (values by
    name, cells in i, j and k)."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    zones = []
    for line in lines[1:]:
        if line.lstrip().startswith("ZONE"):
            zones.append([tuple(int(re.search(k + r'=(\d+)', line).group(1)) - 1 for k in ('I', 'J', 'K')), []])
        else:
            zones[-1][1] += [float(v) for v in line.split()]
    out = []
    for n, values in zones:
        nn, nc, z, pos = (n[0] + 1) * (n[1] + 1) * (n[2] + 1), n[0] * n[1] * n[2], {}, 0
        for k, name in enumerate(names):
            cnt = nn if k < 3 else nc
            z[name] = values[pos:pos + cnt]
            pos += cnt
        out.append((z, n))
    return out


SPLIT_Q = ("rho", "rhou", "rhov", "rhow", "rhoe")


def split(whole, cutj):
    """split: r/dt of every cell of split-cutj against the same cell of split-whole, over the largest |r/dt| of
    split-whole, of rho, of rho E, and one for rho u, rho v and rho w."""
    (a, (ni, nj, nk)), = read_zones(whole)
    fa = {q: [x / t for x, t in zip(a[q], a["dt"])] for q in SPLIT_Q}
    fb = {q: [None] * (ni * nj * nk) for q in SPLIT_Q}
    for z, (b, (mi, mj, mk)) in enumerate(read_zones(cutj)):
        for k in range(mk):
            for j in range(mj):
                for i in range(mi):
                    cl, cg = (k * mj + j) * mi + i, (k * nj + j + z * mj) * ni + i
                    for q in SPLIT_Q:
                        fb[q][cg] = b[q][cl] / b["dt"][cl]
    smom = max(abs(x) for q in ("rhou", "rhov", "rhow") for x in fa[q])
    ok = True
    for q in SPLIT_Q:
        scale = smom if q in ("rhou", "rhov", "rhow") else max(abs(x) for x in fa[q])
        d, c = max((abs(x - y), c) for c, (x, y) in enumerate(zip(fa[q], fb[q])))
        ok &= d / scale <= TOL_SPLIT
        where = f"cell ({c % ni + 1}, {c // ni % nj + 1}, {c // (ni * nj) + 1})"
        print(f"    {q}: largest |r/dt(cutj) - r/dt(whole)| over the scale {scale:.3e}: {d / scale:.3e} ({where}); "
              f"tolerance {TOL_SPLIT:.0e}")
    sys.exit(0 if ok else 1)


def ring101(ns45, eu45, ns101, eu101):
    """ring101: block 1 of the turn joined by 101 against ring45, closed by type 200, cell by cell: r/dt of the Euler step
    over the largest |r/dt| of ring45 of each variable, and the viscous force (r/dt of Navier-Stokes minus Euler) over the
    scales of ring."""
    (s_ns, n), = read_zones(ns45)
    (s_eu, _), = read_zones(eu45)
    r_ns, r_eu = read_zones(ns101)[0][0], read_zones(eu101)[0][0]
    ni, nj, _ = n

    def steps(z):
        return {q: [x / t for x, t in zip(z[q], z["dt"])] for q in SPLIT_Q}

    a_ns, a_eu, b_ns, b_eu = steps(s_ns), steps(s_eu), steps(r_ns), steps(r_eu)
    sm = MU * (U + 2.0 * G / R0 ** 2 + 4.0 * C * R1 ** 3) / (R0 * math.pi / 8.0)
    se = sm * (U * R1 + G / R0 + C * R1 ** 4)
    ok = True
    for q in SPLIT_Q:
        scale = max(abs(x) for x in a_eu[q])
        d, c = max((abs(x - y), c) for c, (x, y) in enumerate(zip(a_eu[q], b_eu[q])))
        ok &= d / scale <= TOL_SPLIT
        print(f"    euler {q}: largest |r/dt(ring101) - r/dt(ring45)| over the scale {scale:.3e}: {d / scale:.3e} "
              f"(cell ({c % ni + 1}, {c // ni % nj + 1}, {c // (ni * nj) + 1})); tolerance {TOL_SPLIT:.0e}")
    for q in CONS:
        s = se if q == "rhoe" else sm
        d, c = max((abs((x1 - x2) - (y1 - y2)), c)
                   for c, (x1, x2, y1, y2) in enumerate(zip(a_ns[q], a_eu[q], b_ns[q], b_eu[q])))
        ok &= d / s <= TOL_RING
        print(f"    viscous force {q}: largest |f(ring101) - f(ring45)| over the scale: {d / s:.3e} "
              f"(cell ({c % ni + 1}, {c // ni % nj + 1}, {c // (ni * nj) + 1})); tolerance {TOL_RING:.0e}")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case" and sys.argv[2].startswith("shift"):
        write_shift(sys.argv[2][len("shift"):])
    elif what == "case" and sys.argv[2].startswith("wall"):
        write_wall(sys.argv[2])
    elif what == "case" and sys.argv[2].startswith("split-"):
        write_split(sys.argv[2][len("split-"):])
    elif what == "case" and sys.argv[2] == "ring101":
        write_ring101()
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
    elif what == "split":
        split(*sys.argv[2:4])
    elif what == "ring101":
        ring101(*sys.argv[2:6])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
