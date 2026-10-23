"""Inputs of the stationary-wall case of the RotatingFrame fast test, and the check of its wall output.

A short tube along x in a rigid-body frame turning at OMEGA about the x axis: NX cells along 0 <= x <= LX, one cell across,
y from Y0 to Y0 + H and z from -H/2 to H/2. Calorically perfect air of test/common/Air at rest in the absolute frame:
p = 1 bar, T = 300 K, and in the rotating frame the relative velocity -Omega x r = OMEGA (0, z, -y) of the cell centre.
Face 1 (x = 0) is an adiabatic no-slip wall (301) that does not rotate with the frame (stationary-face-1 = 1 1): at rest in
the absolute frame, like the gas, and moving with -Omega x r in the rotating frame. Its normal is the axis, so the wall and
the cell next to it move alike, and the shear on it is that of the gas at rest against a wall at rest: none, apart from the
one time step the run takes. Taken as rotating with the frame (wall velocity 0), the wall would see the relative velocity
of the cell, OMEGA y_c, across half a cell: a shear of the order of mu OMEGA y_c / (dx/2). Face 2: extrapolation (400); j
faces: slip walls (300), across which the gas does not move (v = OMEGA z = 0 at the cell centre); k faces: null (0), so that
the swirl -OMEGA y_c along z meets no wall.

Also the boundary file of the first case of the test: the Sod tube of test/1D/Sod79 with its two ends (faces 1 and 2,
extrapolation 400 there) made adiabatic walls 301, so that they can be listed as stationary faces.

And a stationary face that is a wall in part (mixed): the same tube with two cells across, y from Y0 to Y0 + 2H, whose face 1
is a wall 302 at T for the cell j = 1 and, for the cell j = 2, a block connection (101) to a second block, NX2 cells along
-NX2 LX/NX <= x <= 0 at the height of that cell (faces 1: extrapolation 400, j faces: 300, k faces: 0). A connection never
reads the stationary flag, so the face can be listed as stationary; the wall part must then be at rest like the single wall
above. With mixed-symmetry the cell j = 2 of face 1 is a symmetry plane (300) instead, and there is no second block: a 300
takes the velocity of a wall at rest in its inviscid flux only, so that face cannot be stationary.

Usage (from the test directory):
  python3 -B make_case.py case             INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py check FILE       largest wall shear on face 1 in FILE (a wall.tec) over mu OMEGA y_c / (dx/2)
  python3 -B make_case.py sod-walls BC     INPUT/bc.txt: BC (the Sod tube) with walls 301 on faces 1 and 2
  python3 -B make_case.py mixed            INPUT/ic.tec and INPUT/bc.txt of the face that is a wall in part
  python3 -B make_case.py mixed-symmetry   the same with a symmetry plane in place of the connection, one block
"""
import re
import sys

NX = 10
LX = 0.05                     # tube length [m]: cells of 5 mm
Y0, H = 0.1, 0.01             # the tube spans y in [Y0, Y0 + H], z in [-H/2, H/2]
OMEGA = 100.0                 # [rad/s] about x
P, T = 1.0e5, 300.0
MU = 1.846e-5                 # viscosity of air at 300 K [Pa s], the scale of the check only
RGAS = 8314.46261815324 / 28.970418    # universal gas constant of FLINT / molar mass of air
NX2 = 2                       # cells of the second block of the mixed case, as long as those of the first


def write_case():
    import os
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):      # the Sod links of the first case of the test, not their targets
        if os.path.lexists(p):
            os.remove(p)
    rho = P / (RGAS * T)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(NX + 1):
                        f.write(f" {(LX * i / NX, Y0 + H * j, H * (k - 0.5))[axis]:.15E}\n")
        yc, zc = Y0 + 0.5 * H, 0.0    # cell centre: mean of the eight nodes
        for value in (rho, 0.0, OMEGA * zc, -OMEGA * yc, P):
            f.write(f" {value:.15E}\n" * NX)
    with open("INPUT/bc.txt", "w") as f:
        f.write("       1       1       1       1       1     301\n    0.000000E+00, 0.000000E+00, 0.000000E+00,\n")
        f.write(f"       1     {NX:3d}       1       1       2     400\n")
        for face, t in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {t:3d}\n")


def write_mixed(second):
    """Face 1 of block 1: wall 302 at T for j = 1, and for j = 2 a connection 101 to block 2 (second = 101) or a symmetry
    plane 300 (second = 300, no block 2). The gas at rest in the absolute frame, as in write_case."""
    import os
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    rho = P / (RGAS * T)
    dx = LX / NX
    blocks = [([dx * i for i in range(NX + 1)], [Y0, Y0 + H, Y0 + 2 * H])]
    if second == 101:
        blocks.append(([-dx * (NX2 - i) for i in range(NX2 + 1)], [Y0 + H, Y0 + 2 * H]))
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b, (xs, ys) in enumerate(blocks, start=1):
            ni, nj = len(xs) - 1, len(ys) - 1
            f.write(f" ZONE  T = B{b}-IG, I={ni + 1}, J={nj + 1}, K=2, DATAPACKING=BLOCK, "
                    "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for axis in range(3):
                for k in range(2):
                    for j in range(nj + 1):
                        for i in range(ni + 1):
                            f.write(f" {(xs[i], ys[j], H * (k - 0.5))[axis]:.15E}\n")
            cells = [(j, i) for j in range(nj) for i in range(ni)]
            yc = [0.5 * (ys[j] + ys[j + 1]) for j, _ in cells]      # cell centres, z_c = 0
            for column in ([rho] * len(cells), [0.0] * len(cells), [0.0] * len(cells), [-OMEGA * y for y in yc],
                           [P] * len(cells)):
                f.write("".join(f" {v:.15E}\n" for v in column))
    rec = "       {}     {:3d}     {:3d}       1       {}     {:3d}\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write(rec.format(1, 1, 1, 1, 302) + f"    {T:.6E}, 0.000000E+00, 0.000000E+00,\n")
        if second == 101:
            f.write(rec.format(1, 1, 2, 1, 101) + f"{2:8d}{NX2:8d}{1:8d}{1:8d}{2:8d}{1:8d}{0:8d}{0:8d}{1:8d}\n")
        else:
            f.write(rec.format(1, 1, 2, 1, 300))
        for j in (1, 2):
            f.write(rec.format(1, NX, j, 2, 400))
        for face, t in ((3, 300), (4, 300)):
            for i in range(1, NX + 1):
                f.write(rec.format(1, i, 1 if face == 3 else 2, face, t))
        for face in (5, 6):
            for j in (1, 2):
                for i in range(1, NX + 1):
                    f.write(rec.format(1, i, j, face, 0))
        if second == 101:
            f.write(rec.format(2, 1, 1, 1, 400))
            f.write(rec.format(2, NX2, 1, 2, 101) + f"{1:8d}{1:8d}{2:8d}{1:8d}{1:8d}{1:8d}{0:8d}{0:8d}{1:8d}\n")
            for face in (3, 4, 5, 6):
                for i in range(1, NX2 + 1):
                    f.write(rec.format(2, i, 1, face, 300 if face < 5 else 0))


def read_zone(path):
    """Values of a one-zone MOSE ASCII file (BLOCK packing): nodal x, y, z first, then the cell variables."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, max(ni - 1, 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out


def check(path):
    w = read_zone(path)
    scale = MU * OMEGA * (Y0 + 0.5 * H) / (0.5 * LX / NX)
    tau = max(max(abs(a), abs(b), abs(c)) for a, b, c in zip(w["tauX"], w["tauY"], w["tauZ"]))
    print(f"{tau / scale:.3e}")


def sod_walls(src):
    import os
    if os.path.lexists("INPUT/bc.txt"):
        os.remove("INPUT/bc.txt")
    with open(src) as f, open("INPUT/bc.txt", "w") as g:
        for line in f:
            tok = line.split()
            if len(tok) == 6 and tok[4] in ("1", "2") and tok[5] == "400":
                g.write(line.rsplit("400", 1)[0] + "301\n")
                g.write("    0.000000E+00, 0.000000E+00, 0.000000E+00,\n")
            else:
                g.write(line)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case()
    elif what == "sod-walls":
        sod_walls(sys.argv[2])
    elif what == "mixed":
        write_mixed(101)
    elif what == "mixed-symmetry":
        write_mixed(300)
    elif what == "check":
        check(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
