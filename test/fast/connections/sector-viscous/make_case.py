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
      planes at +-45 degrees); each of them has zero viscous force. diffusive-metric = inverse-mean, the metric that the
      periodic face takes, so that the two meshes compute the same face. The viscous force of a cell is r/dt of a
      Navier-Stokes step minus r/dt of an Euler step from the same state; in every cell of S45, for rho u, rho v, rho w
      and rho E, it must equal that of the same cell of R90 within 1e-9 of the scale
        S = mu (U + 2 G/r0^2 + 4 C r1^3)/(r0 dtheta)   (momentum),   S (U r1 + G/r0 + C r1^4)   (energy),
      r0 = 0.5, r1 = 1, dtheta = pi/8.
mirror The sector from -45 to 45 degrees in 2 x 8 x 8 cells (uniform radii, node planes at theta = m pi/16) with the
      radial field alone, k faces of type 200 against type 300: the radial field is its own mirror image across the
      meridian planes, so the two conditions set the same ghost cells, and the viscous stress on the planes is the
      normal stress (2/3) mu U. The default face metric. A Navier-Stokes step with each: in the cells not next to the j
      faces, r/dt of rho u, rho v and rho w must agree within 0.02 (2/3) mu U A/V, with A/V = 1/(r dtheta) the area of a k
      face over the volume (r the mid radius of the cell), and r/dt of rho E within the same times U r1. (A field with shear
      on the plane is not its own mirror image, and the symmetry condition would not set the same problem.)

One explicit Euler step with local time stepping, first-order reconstruction. OUTPUT/diagnostic.tec holds, per cell, the
increment of each conserved variable over the step, and the step dt.

Usage (from the test directory):
  python3 -B make_case.py case ring45|ring90|mirror200|mirror300   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py ring NS45 EU45 NS90 EU90                 largest difference of the viscous forces over the scale
  python3 -B make_case.py mirror F200 F300                         largest difference over (2/3) mu U A/V
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


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "ring":
        ring(*sys.argv[2:6])
    elif what == "mirror":
        mirror(*sys.argv[2:4])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
