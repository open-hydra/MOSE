"""Part 1 of Inlets (run.sh): the far-field inlet with the flow direction normal, inputs and a reader of the outputs.

A planar 2-D duct along x: NX x NY cells on 0 <= x <= LX, 0 <= y <= LY, one cell of LZ in z; slip walls (300) on the j and k
faces. Calorically perfect air of test/common/Air at p = 1 bar and T = 300 K, flowing in along x at u = 50 m/s. Face 1 is an
ambient boundary (407: T0, p0, p_back, alpha, beta, rel_fac, massf) where the gas enters; face 2 an outlet at the back
pressure (406). The flow direction of the 407 records is given either as the token normal, the face normal, or as the
angles that it stands for: the normal of face 1 is (1, 0, 0), so alpha = atan(n_y / n_x) = 0 and beta = atan(n_z / n_x) = 0.
The two inputs are the same boundary condition and must give the same fields to the bit.

Usage (from the test directory):
  python3 -B farfield_normal.py case DIR      INPUT/ic.tec and INPUT/bc.txt, DIR = normal or angle
  python3 -B farfield_normal.py check FILE    largest |v| over u in FILE (a field.tec)
"""
import os
import re
import sys

NX, NY = 8, 4
LX, LY, LZ = 0.08, 0.04, 0.004
P, T, U = 1.0e5, 300.0, 50.0
GAMMA = 1.4
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air
CP = GAMMA * RGAS / (GAMMA - 1.0)
T0 = T + U * U / (2.0 * CP)
P0 = P * (T0 / T) ** (GAMMA / (GAMMA - 1.0))


def write_case(direction):
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    rho = P / (RGAS * T)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J={NY + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(NY + 1):
                    for i in range(NX + 1):
                        f.write(f" {(LX * i / NX, LY * j / NY, LZ * k)[axis]:.15E}\n")
        for value in (rho, U, 0.0, 0.0, P):
            f.write(f" {value:.15E}\n" * (NX * NY))
    ang = "normal" if direction == "normal" else "0.000000E+00"
    with open("INPUT/bc.txt", "w") as f:
        for j in range(1, NY + 1):
            f.write(f"       1       1     {j:3d}       1       1     407\n")
            f.write(f"    {T0:.15E}, {P0:.15E}, {P:.15E}, {ang}, {ang}, 1.0, 1.0,\n")
        for j in range(1, NY + 1):
            f.write(f"       1     {NX:3d}     {j:3d}       1       2     406\n")
            f.write(f"    {P:.15E}, 1.0,\n")
        for face, j in ((3, 1), (4, NY)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     300\n")
        for face in (5, 6):
            for j in range(1, NY + 1):
                for i in range(1, NX + 1):
                    f.write(f"       1     {i:3d}     {j:3d}       1       {face}     300\n")


def check(path):
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = (NX + 1) * (NY + 1) * 2, NX * NY, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    print(f"{max(abs(v) for v in out['v']) / U:.3e}")


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check":
        check(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
