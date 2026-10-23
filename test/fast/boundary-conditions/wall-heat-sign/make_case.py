"""Inputs of the WallHeatSign fast test, and the checks of its outputs.

A channel along x between two walls, at y = 0 and y = LY, that mirror each other: the same wall on the two faces of one
index direction DIR, one face of each parity. Walls on the faces of
  j  (faces 3 and 4): a planar 2-D channel of 4 x 4 cells on 0 <= x <= LX, 0 <= y <= LY, one cell of LZ in z, with null (0)
     k faces;
  i  (faces 1 and 2): a 3-D channel of 4 x 4 x 4 cells on 0 <= x <= LX, 0 <= y <= LY, 0 <= z <= 4 LZ, the indices i, j, k
     along y, z, x;
  k  (faces 5 and 6): the same 3-D channel, the indices i, j, k along z, x, y;
the faces normal to x are slip walls (300), or extrapolation (400) when the gas slides along x, and in 3-D the faces normal
to z are slip walls. Viscous air (test/common/Air-viscous, Sutherland transport of test/common/Air-sutherland), p = 1 bar
and T = 300 K. One forward Euler step with first-order reconstruction. The wall output writes the same quantities on every
face (docs/user/boundary-conditions.md): qw is the heat flux into the gas, like q_wall, and tauX, tauY, tauZ the traction
of the gas on the wall, tau . n with n the normal into the gas. The two walls, mirror images, must then write the same
qw, and the mirror images of their tractions: the same tauX and tauZ, opposite tauY.
  flux   type 301 with q_wall = Q on both walls, positive into the fluid, the gas at rest: both walls heat the gas, so the
         cells next to them must warm up alike, and qw must be Q on both walls;
  temp   type 302 with T_wall = TW > T on both walls, the gas at rest: the walls heat the gas, and qw must be the same on the
         two walls and positive;
  shear  type 302 with T_wall = T on both walls, the gas sliding along x at U: the gas drags both walls along x, so tauX must
         be the same on the two walls and positive; tauY must be opposite on the two walls, tauZ the same (the ends of the
         channel give the cells next to them a small v, and the stress a small tauY there).

Usage (from the test directory; DIR is i, j or k):
  python3 -B make_case.py case flux|temp|shear DIR  INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py check-flux DIR FIELD      temperatures of the cells next to the two walls in FIELD (a field.tec)
  python3 -B make_case.py check-qw-flux DIR WALL    qw of the two walls in WALL (a wall.tec) of the case flux
  python3 -B make_case.py check-qw-temp DIR WALL    qw of the two walls in WALL (a wall.tec) of the case temp
  python3 -B make_case.py check-tau DIR WALL        tauX, tauY and tauZ of the two walls in WALL (a wall.tec) of the case shear
"""
import os
import re
import sys

NC = 4                               # cells along x and y, and along z in 3-D
LX, LY, LZ = 0.04, 0.004, 0.001      # LZ: the size of the cells along z
P, T = 1.0e5, 300.0
Q, TW = 1.0e4, 400.0
U = 10.0                             # case shear: the speed of the gas along x, m/s
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air
AXES = {"i": (1, 2, 0), "j": (0, 1, 2), "k": (2, 0, 1)}   # the axis (0 x, 1 y, 2 z) of i, j and k: a rotation each


def grid(d):
    """The axes of i, j and k, the cells along them and the lengths along x, y, z of the channel with walls on DIR d."""
    axes = AXES[d]
    cells = (NC, NC, 1 if d == "j" else NC)                  # along x, y, z
    return axes, [cells[a] for a in axes], cells, (LX, LY, LZ * cells[2])


def cells_on_face(n, m, side):
    """The cells (i, j, k), from 1, next to face 2 m + 1 (side 0) or 2 m + 2 (side 1), i running fastest."""
    out = []
    for k in range(1, n[2] + 1):
        for j in range(1, n[1] + 1):
            for i in range(1, n[0] + 1):
                if (i, j, k)[m] == (1 if side == 0 else n[m]):
                    out.append((i, j, k))
    return out


def write_case(kind, d):
    axes, n, cells, length = grid(d)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={n[0] + 1}, J={n[1] + 1}, K={n[2] + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            m = axes.index(axis)                             # the index that runs along this axis
            for k in range(n[2] + 1):
                for j in range(n[1] + 1):
                    for i in range(n[0] + 1):
                        f.write(f" {length[axis] * (i, j, k)[m] / cells[axis]:.15E}\n")
        for value in (P / (RGAS * T), U if kind == "shear" else 0.0, 0.0, 0.0, P):
            f.write(f" {value:.15E}\n" * (n[0] * n[1] * n[2]))
    t, v = {"flux": (301, Q), "temp": (302, TW), "shear": (302, T)}[kind]
    with open("INPUT/bc.txt", "w") as f:
        for m in range(3):
            if axes[m] == 1:                                 # the two walls
                btype = t
            elif axes[m] == 0:                               # the gas of shear leaves and enters through the faces normal to x
                btype = 400 if kind == "shear" else 300
            else:
                btype = 0 if d == "j" else 300
            for side in (0, 1):
                for i, j, k in cells_on_face(n, m, side):
                    f.write(f"       1     {i:3d}     {j:3d}     {k:3d}       {2 * m + 1 + side}     {btype:3d}\n")
                    if btype in (301, 302):
                        f.write(f"    {v:.15E}, 0.0, 0.0,\n")


def zones(path):
    """The zones of a MOSE ASCII file (BLOCK packing): {title: {name: values}}, cell-centred values after the nodal x, y, z."""
    with open(path) as f:
        text = "".join(l for l in f if not l.startswith("DATASETAUXDATA"))
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', text.splitlines()[0].split('=', 1)[1])]
    out = {}
    parts = re.split(r'(?m)^\s*ZONE', text)[1:]
    for part in parts:
        head, _, body = part.partition('\n')
        title = re.search(r'T\s*=\s*"?([^",]+)', head).group(1).strip()
        ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', head).group(1)) for k in ('I', 'J', 'K'))
        values = [float(v) for v in body.split()]
        nn, nc, z, pos = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
        for k, name in enumerate(names):
            cnt = nn if k < 3 else nc
            z[name] = values[pos:pos + cnt]
            pos += cnt
        out[title] = z
    return out


def walls(d, path):
    """The faces of the two walls, low and high, and their zones in WALL."""
    m = "ijk".index(d)
    z = zones(path)
    return 2 * m + 1, 2 * m + 2, z[f"B1F{2 * m + 1}"], z[f"B1F{2 * m + 2}"]


def check_flux(d, path):
    _, n, _, _ = grid(d)
    m = "ijk".index(d)
    temp = next(iter(zones(path).values()))["T"]
    lo, hi = ([temp[(i - 1) + n[0] * ((j - 1) + n[1] * (k - 1))] for i, j, k in cells_on_face(n, m, side)]
              for side in (0, 1))
    d_max = max(abs(a - b) for a, b in zip(lo, hi))
    rise = min(hi) - T
    ok = rise > 0 and min(lo) > T and d_max <= 1e-9 * abs(rise)
    print(f"    T next to face {2 * m + 1} (y = 0): {min(lo):.10f} K, next to face {2 * m + 2}: {min(hi):.10f} K "
          f"(start {T} K); largest difference {d_max:.3e} K")
    sys.exit(0 if ok else 1)


def check_qw_flux(d, path):
    f1, f2, z1, z2 = walls(d, path)
    q1, q2 = z1["qw"], z2["qw"]
    dev = max(abs(v - Q) for v in q1 + q2)
    print(f"    qw of face {f1} (y = 0): {min(q1):.6e} to {max(q1):.6e} W/m2, of face {f2}: {min(q2):.6e} to "
          f"{max(q2):.6e} W/m2 (q_wall = {Q:.1e} W/m2)")
    sys.exit(0 if dev <= 1e-12 * Q else 1)


def check_qw_temp(d, path):
    f1, f2, z1, z2 = walls(d, path)
    q1, q2 = z1["qw"], z2["qw"]
    dev = max(abs(a - b) for a, b in zip(q1, q2))
    ok = min(q1 + q2) > 0 and dev <= 1e-9 * max(abs(v) for v in q2)
    print(f"    qw of face {f1} (y = 0): {min(q1):.6e} to {max(q1):.6e} W/m2, of face {f2}: {min(q2):.6e} to "
          f"{max(q2):.6e} W/m2")
    sys.exit(0 if ok else 1)


def check_tau(d, path):
    f1, f2, z1, z2 = walls(d, path)
    t1, t2 = z1["tauX"], z2["tauX"]
    scale = max(abs(v) for v in t1 + t2)
    dev = max(abs(a - b) for a, b in zip(t1, t2))
    dev_y = max(abs(a + b) for a, b in zip(z1["tauY"], z2["tauY"]))      # mirror images: opposite tauY
    dev_z = max(abs(a - b) for a, b in zip(z1["tauZ"], z2["tauZ"]))
    big_y = max(abs(v) for v in z1["tauY"] + z2["tauY"])
    ok = min(t1 + t2) > 0 and max(dev, dev_y, dev_z) <= 1e-9 * scale
    print(f"    tauX of face {f1} (y = 0): {min(t1):.6e} to {max(t1):.6e} Pa, of face {f2}: {min(t2):.6e} to "
          f"{max(t2):.6e} Pa; largest |tauY| {big_y:.3e} Pa; largest difference of tauX {dev:.1e} Pa, "
          f"sum of tauY {dev_y:.1e} Pa, difference of tauZ {dev_z:.1e} Pa")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2], sys.argv[3])
    elif what == "check-flux":
        check_flux(sys.argv[2], sys.argv[3])
    elif what == "check-qw-flux":
        check_qw_flux(sys.argv[2], sys.argv[3])
    elif what == "check-qw-temp":
        check_qw_temp(sys.argv[2], sys.argv[3])
    elif what == "check-tau":
        check_tau(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
