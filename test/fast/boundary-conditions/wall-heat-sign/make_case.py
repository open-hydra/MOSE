"""Inputs of the WallHeatSign fast test, and the checks of its outputs.

A planar 2-D channel along x: NX x NY cells on 0 <= x <= LX, 0 <= y <= LY, one cell of LZ in z; slip walls (300) on the i
faces, null (0) k faces, and the same wall on both j faces: face 3 (y = 0) and face 4 (y = LY) mirror each other. Viscous air
(test/common/Air-viscous, Sutherland transport of test/common/Air-sutherland) at rest, p = 1 bar and T = 300 K. One forward Euler
step with first-order reconstruction:
  flux  type 301 with q_wall = Q on both walls, positive into the fluid (docs/user/boundary-conditions.md): both walls heat the
        gas, so the cells next to them must warm up alike. The heat flux written in the wall output, qw, is k grad T . n along
        the face normal n toward growing index: -Q on face 3, where n points into the fluid, and Q on face 4;
  temp  type 302 with T_wall = TW > T on both walls: the walls heat the gas, and qw must be opposite on the two walls and
        positive on face 4.

Usage (from the test directory):
  python3 -B make_case.py case flux|temp        INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py check-flux FIELD      temperatures of the cells next to the two walls in FIELD (a field.tec)
  python3 -B make_case.py check-qw-flux WALL    qw of the two walls in WALL (a wall.tec) of the case flux
  python3 -B make_case.py check-qw-temp WALL    qw of the two walls in WALL (a wall.tec) of the case temp
"""
import os
import re
import sys

NX, NY = 4, 4
LX, LY, LZ = 0.04, 0.004, 0.001
P, T = 1.0e5, 300.0
Q, TW = 1.0e4, 400.0
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air


def write_case(kind):
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J={NY + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(NY + 1):
                    for i in range(NX + 1):
                        f.write(f" {(LX * i / NX, LY * j / NY, LZ * k)[axis]:.15E}\n")
        for value in (P / (RGAS * T), 0.0, 0.0, 0.0, P):
            f.write(f" {value:.15E}\n" * (NX * NY))
    t, v = (301, Q) if kind == "flux" else (302, TW)
    with open("INPUT/bc.txt", "w") as f:
        for face, i in ((1, 1), (2, NX)):
            for j in range(1, NY + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     300\n")
        for face, j in ((3, 1), (4, NY)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     {t}\n")
                f.write(f"    {v:.15E}, 0.0, 0.0,\n")
        for face in (5, 6):
            for j in range(1, NY + 1):
                for i in range(1, NX + 1):
                    f.write(f"       1     {i:3d}     {j:3d}       1       {face}       0\n")


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


def check_flux(path):
    z = next(iter(zones(path).values()))
    temp = z["T"]
    lo = [temp[i] for i in range(NX)]                            # j = 1, next to face 3
    hi = [temp[(NY - 1) * NX + i] for i in range(NX)]             # j = NY, next to face 4
    d = max(abs(a - b) for a, b in zip(lo, hi))
    rise = min(hi) - T
    ok = rise > 0 and min(lo) > T and d <= 1e-9 * abs(rise)
    print(f"    T next to face 3 (y = 0): {min(lo):.10f} K, next to face 4: {min(hi):.10f} K (start {T} K); "
          f"largest difference {d:.3e} K")
    sys.exit(0 if ok else 1)


def check_qw_flux(path):
    z = zones(path)
    q3, q4 = z["B1F3"]["qw"], z["B1F4"]["qw"]
    d = max([abs(v + Q) for v in q3] + [abs(v - Q) for v in q4])
    print(f"    qw of face 3 (y = 0): {min(q3):.6e} to {max(q3):.6e} W/m2, of face 4: {min(q4):.6e} to {max(q4):.6e} W/m2 "
          f"(q_wall = {Q:.1e} W/m2)")
    sys.exit(0 if d <= 1e-12 * Q else 1)


def check_qw_temp(path):
    z = zones(path)
    q3, q4 = z["B1F3"]["qw"], z["B1F4"]["qw"]
    d = max(abs(a + b) for a, b in zip(q3, q4))
    ok = min(q4) > 0 and d <= 1e-9 * max(abs(v) for v in q4)
    print(f"    qw of face 3 (y = 0): {min(q3):.6e} to {max(q3):.6e} W/m2, of face 4: {min(q4):.6e} to {max(q4):.6e} W/m2")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check-flux":
        check_flux(sys.argv[2])
    elif what == "check-qw-flux":
        check_qw_flux(sys.argv[2])
    elif what == "check-qw-temp":
        check_qw_temp(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
