"""Inputs of the Quasi2DSymmetry fast test, and a reader of its outputs.

Quasi-2D box, one cell in z: NX x NY cells on 0 <= x <= LX, 0 <= y <= LY,
thickness LZ in z (much smaller than the cell size in x and y, as in a planar
case meshed with a thin layer). Calorically perfect air of test/common/Air at
uniform pressure and velocity (p = 1 bar, u = 100 m/s) carrying an entropy wave
of amplitude A in density (T = 300 K on average); extrapolation (400) on the i
and j faces. The wave keeps the run unsteady, so that the residual stop does
not end it after one step, while p, u and v stay uniform to round-off. The k
faces are symmetry planes (300) or null faces (0).

The z velocity starts at W0 in every cell: 1e-8 m/s, the size of the round-off
that a real mesh leaves in w. Between two symmetry planes one cell apart, w
obeys dw/dt = -2 a w / LZ (a the speed of sound): with the time step set by the
x and y cells at cfl 0.5, 2 a dt / LZ is about 2.2, outside the real stability
interval of Heun, [-2, 0], and inside that of RK3, [-2.51, 0] (a thinner layer
leaves that of RK3 too). W0 must not grow. A larger W0 (1 m/s: a field restarted from a run
already hit by the mode, or from a 3-D field) must not spoil p, u and v either.

1-D tube (ic1d, bc1d, check1d): NX cells along x, one cell in y and one in z, both LZ thick; the same
gas and wave, and y and z velocities V0 = 1e-8 m/s. With symmetry planes on all the j faces v obeys
dv/dt = -2 a v / LZ, and with symmetry planes on all the k faces w obeys dw/dt = -2 a w / LZ: neither
may grow. Extrapolation (400) on the j faces and null k faces are the control.

Usage (from the test directory):
  python3 -B make_case.py ic [W0]            INPUT/ic.tec (W0 = 1e-8 m/s if omitted)
  python3 -B make_case.py bc KTYPE           INPUT/bc.txt with type KTYPE (300 or 0) on the k faces
  python3 -B make_case.py check FILE [W0]    compare a solution file with the uniform state
  python3 -B make_case.py ic1d               INPUT/ic.tec of the 1-D tube
  python3 -B make_case.py bc1d JTYPE KTYPE   INPUT/bc.txt of the 1-D tube (JTYPE 300 or 400, KTYPE 300 or 0)
  python3 -B make_case.py check1d FILE       compare a solution file of the 1-D tube with the uniform state
  python3 -B make_case.py channel KIND        INPUT/ic.tec and INPUT/bc.txt of the channel of variable section
                                             (KIND layer: quasi-2-D, KIND cd1d: 1-D)
  python3 -B make_case.py checkchannel FILE KIND  compare a solution file with the quasi-1-D isentropic flow
"""
import sys

NX, NY = 8, 4        # cells along x and y
LX, LY = 0.04, 0.02  # box size [m]: cells of 5 mm
LZ = 1.0e-3          # k thickness [m]
T, P, U = 300.0, 1.0e5, 100.0
A = 1.0e-3           # relative amplitude of the entropy wave
W0 = 1.0e-8          # default initial z velocity [m/s]
RGAS = 8314.46261815324 / 28.970418   # universal gas constant of FLINT / molar mass of air


def write_ic(w0, path="INPUT/ic.tec"):
    import math
    rho0 = P / (RGAS * T)
    rho = [rho0 * (1.0 + A * math.sin(2.0 * math.pi * (i + 0.5) / NX))
           for _ in range(NY) for i in range(NX)]
    with open(path, "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J={NY + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(NY + 1):
                    for i in range(NX + 1):
                        f.write(f" {(LX * i / NX, LY * j / NY, LZ * k)[axis]:.15E}\n")
        f.write("".join(f" {r:.15E}\n" for r in rho))
        for value in (U, 0.0, w0, P):
            f.write(f" {value:.15E}\n" * (NX * NY))


def write_bc(ktype, path="INPUT/bc.txt"):
    """One line per boundary cell face, in order of face (block 1 only)."""
    with open(path, "w") as f:
        for face, i in ((1, 1), (2, NX)):
            for j in range(1, NY + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     400\n")
        for face, j in ((3, 1), (4, NY)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     400\n")
        for face in (5, 6):
            for j in range(1, NY + 1):
                for i in range(1, NX + 1):
                    f.write(f"       1     {i:3d}     {j:3d}       1       {face}     {ktype:3d}\n")


def write_ic_1d(path="INPUT/ic.tec"):
    """1-D tube: NX cells along x, one cell of LZ in y and in z; v = w = W0."""
    import math
    rho0 = P / (RGAS * T)
    rho = [rho0 * (1.0 + A * math.sin(2.0 * math.pi * (i + 0.5) / NX)) for i in range(NX)]
    with open(path, "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(NX + 1):
                        f.write(f" {(LX * i / NX, LZ * j, LZ * k)[axis]:.15E}\n")
        f.write("".join(f" {r:.15E}\n" for r in rho))
        for value in (U, W0, W0, P):
            f.write(f" {value:.15E}\n" * NX)


def write_bc_1d(jtype, ktype, path="INPUT/bc.txt"):
    """One line per boundary cell face of the 1-D tube, face by face."""
    with open(path, "w") as f:
        for face, i in ((1, 1), (2, NX)):
            f.write(f"       1     {i:3d}       1       1       {face}     400\n")
        for face, t in ((3, jtype), (4, jtype), (5, ktype), (6, ktype)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {t:3d}\n")


def check_1d(path):
    """Exit status 0 when p and u are still uniform and neither v nor w has grown."""
    field = read_field(path, NX, 1)
    worst = {
        "u": max(abs(u / U - 1.0) for u in field["u"]),
        "p": max(abs(p / P - 1.0) for p in field["p"]),
    }
    vmax = max(abs(v) for v in field["v"])
    wmax = max(abs(w) for w in field["w"])
    print("    largest relative deviation of p, u from the uniform state: " +
          ", ".join(f"{k} {v:.2e}" for k, v in worst.items()) + f"; max |v| = {vmax:.3e}, max |w| = {wmax:.3e} m/s")
    ok = all(v < 1.0e-12 for v in worst.values()) and vmax <= 2.0 * W0 and wmax <= 2.0 * W0
    sys.exit(0 if ok else 1)


def read_field(path, nx=None, ny=None):
    """Cell values of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    nx = NX if nx is None else nx
    ny = NY if ny is None else ny
    import re
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = (nx + 1) * (ny + 1) * 2, nx * ny, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out


def check(path, w0):
    """Exit status 0 when p, u and v are still uniform and w has not grown."""
    field = read_field(path)
    worst = {
        "u": max(abs(u / U - 1.0) for u in field["u"]),
        "v": max(abs(v) / U for v in field["v"]),
        "p": max(abs(p / P - 1.0) for p in field["p"]),
    }
    wmax = max(abs(w) for w in field["w"])
    print("    largest relative deviation of p, u, v from the uniform state: " +
          ", ".join(f"{k} {v:.2e}" for k, v in worst.items()) + f"; max |w| = {wmax:.3e} m/s")
    ok = all(v < 1.0e-12 for v in worst.values()) and wmax <= 2.0 * w0
    sys.exit(0 if ok else 1)


# Channel of variable section, steady subsonic flow (layer and cd1d): the section A(x) = 0.2 - 0.08 sin(pi x) [m] on
# 0 <= x <= 1, NXC cells along x; inflow at total conditions (401: T0, P0) on face 1 and outflow at PB (406) on face 2.
# The faces along the walls of the channel are symmetry planes (300), one of them inclined, so that the gas has a
# velocity across the cells next to it: the quasi-1-D isentropic flow of the same A(x) is the reference.
#   layer  quasi-2-D: NYC cells in y from y = 0 to y = A(x) (the j face on top inclined), one cell across the layer,
#          between the parallel planes z = CLY y and z = CLY y + LZL (k faces not normal to z; symmetry there asks w = CLY v)
#   cd1d   1-D: one cell in y from y = 0 to y = A(x) (the j face on top inclined), one cell in z of LZC, null k faces
NXC, NYC, LZC, CLY, LZL = 50, 2, 0.05, 0.5, 0.2
T0C, P0C, PBC, GAMMA = 300.0, 1.0e5, 0.95e5, 1.4


def area(x):
    import math
    return 0.2 - 0.08 * math.sin(math.pi * x)


def write_channel(kind, path_ic="INPUT/ic.tec", path_bc="INPUT/bc.txt"):
    """ic.tec and bc.txt of the channel (kind layer or cd1d), gas at rest at PB, T0."""
    ny = NYC if kind == "layer" else 1
    xs = [i / NXC for i in range(NXC + 1)]
    with open(path_ic, "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NXC + 1}, J={ny + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(ny + 1):
                    for i in range(NXC + 1):
                        x = xs[i]
                        if kind == "layer":
                            point = (x, area(x) * j / ny, CLY * area(x) * j / ny + LZL * k)
                        else:
                            point = (x, area(x) * j, LZC * k)
                        f.write(f" {point[axis]:.15E}\n")
        rho = PBC / (RGAS * T0C)
        for value in (rho, 10.0, 0.0, 0.0, PBC):
            f.write(f" {value:.15E}\n" * (NXC * ny))
    with open(path_bc, "w") as f:
        for j in range(1, ny + 1):
            f.write(f"       1       1     {j:3d}       1       1     401\n")
            f.write(f"    {T0C:.15E}, {P0C:.15E}, 0.0, 0.0, 1.0, 1.0,\n")
        for j in range(1, ny + 1):
            f.write(f"       1     {NXC:3d}     {j:3d}       1       2     406\n")
            f.write(f"    {PBC:.15E}, 1.0,\n")
        for face, j in ((3, 1), (4, ny)):
            for i in range(1, NXC + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     300\n")
        ktype = 300 if kind == "layer" else 0
        for face in (5, 6):
            for j in range(1, ny + 1):
                for i in range(1, NXC + 1):
                    f.write(f"       1     {i:3d}     {j:3d}       1       {face}     {ktype:3d}\n")


def read_cells(path, nx, ny):
    """Cell values of a one-zone MOSE ASCII field of nx x ny x 1 cells (BLOCK packing), names quoted or not."""
    import re
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = (nx + 1) * (ny + 1) * 2, nx * ny, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out

def mach_sub(ar):
    """Subsonic Mach number of the area ratio A/A* (bisection)."""
    lo, hi = 1e-8, 1.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        f = (1.0 / m) * ((2.0 / (GAMMA + 1.0)) * (1.0 + 0.5 * (GAMMA - 1.0) * m * m)) ** ((GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)))
        if f > ar:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def check_channel(path, kind):
    """Exit status 0 when p is within 0.03 P0 of the quasi-1-D isentropic flow of the same section in every cell."""
    import math
    ny = NYC if kind == "layer" else 1
    field = read_cells(path, NXC, ny)
    me = math.sqrt(2.0 / (GAMMA - 1.0) * ((P0C / PBC) ** ((GAMMA - 1.0) / GAMMA) - 1.0))
    ar_e = (1.0 / me) * ((2.0 / (GAMMA + 1.0)) * (1.0 + 0.5 * (GAMMA - 1.0) * me * me)) ** ((GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)))
    astar = area(1.0) / ar_e
    worst = 0.0
    for j in range(ny):
        for i in range(NXC):
            m = mach_sub(area((i + 0.5) / NXC) / astar)
            pex = P0C * (1.0 + 0.5 * (GAMMA - 1.0) * m * m) ** (-GAMMA / (GAMMA - 1.0))
            worst = max(worst, abs(field["p"][j * NXC + i] - pex) / P0C)
    across = max(abs(c) for c in field["w" if kind == "layer" else "v"])
    print(f"    largest |p - p of the quasi-1-D flow| = {worst:.3e} P0; largest velocity across the cells {across:.3e} m/s")
    sys.exit(0 if worst <= 0.03 else 1)


def main():
    what = sys.argv[1]
    if what == "ic":
        write_ic(float(sys.argv[2]) if len(sys.argv) > 2 else W0)
    elif what == "bc":
        write_bc(int(sys.argv[2]))
    elif what == "check":
        check(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else W0)
    elif what == "ic1d":
        write_ic_1d()
    elif what == "bc1d":
        write_bc_1d(int(sys.argv[2]), int(sys.argv[3]))
    elif what == "check1d":
        check_1d(sys.argv[2])
    elif what == "channel":
        write_channel(sys.argv[2])
    elif what == "checkchannel":
        check_channel(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
