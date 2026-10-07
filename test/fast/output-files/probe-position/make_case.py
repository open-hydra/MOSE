"""Inputs of the ProbePosition fast test, and the check of its probe files.

A planar 2-D duct along x: NX x NY cells on 0 <= x <= LX, 0 <= y <= LY, one cell of LZ in z; slip walls (300) on the j faces,
null (0) k faces, extrapolation (400) on the i faces. Calorically perfect air of test/common/Air at rest, at T = 300 K, with a
pressure that differs in every cell, p = P (1 + 0.01 i + 0.001 j) (i, j the cell indexes from 1), so that every cell has its own
values. The probes are given by position (position = x y z), each at the centre of a cell in the mid-plane z = LZ/2; a probe
records the cell whose centre is the nearest, here the cell itself, never a ghost cell.

Usage (from the test directory):
  python3 -B make_case.py case                 INPUT/ic.tec, INPUT/bc.txt and the probe sections of input.ini (appended)
  python3 -B make_case.py check FIELD          the last probe rows against the cells of FIELD, the solution of their step;
                                               exit status 0 when p and rho agree within 1e-12 relative
"""
import os
import re
import sys

NX, NY = 8, 4
LX, LY, LZ = 0.08, 0.04, 0.004
P, T = 1.0e5, 300.0
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air
CELLS = [(2, 2), (5, 3), (8, 4)]     # (i, j) of the cells probed, from 1


def pressure(i, j):
    return P * (1.0 + 0.01 * i + 0.001 * j)


def centre(i, j):
    return (LX * (i - 0.5) / NX, LY * (j - 0.5) / NY, 0.5 * LZ)


def write_case():
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
        cells = [(i, j) for j in range(1, NY + 1) for i in range(1, NX + 1)]
        f.write("".join(f" {pressure(i, j) / (RGAS * T):.15E}\n" for i, j in cells))
        for value in (0.0, 0.0, 0.0):
            f.write(f" {value:.15E}\n" * (NX * NY))
        f.write("".join(f" {pressure(i, j):.15E}\n" for i, j in cells))
    with open("INPUT/bc.txt", "w") as f:
        for j in range(1, NY + 1):
            f.write(f"       1       1     {j:3d}       1       1     400\n")
            f.write(f"       1     {NX:3d}     {j:3d}       1       2     400\n")
        for face, j in ((3, 1), (4, NY)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}     {j:3d}       1       {face}     300\n")
        for face in (5, 6):
            for j in range(1, NY + 1):
                for i in range(1, NX + 1):
                    f.write(f"       1     {i:3d}     {j:3d}       1       {face}       0\n")
    with open("input.ini", "a") as f:
        f.write("\n[MOSE-Probes]\n")
        for n, (i, j) in enumerate(CELLS, 1):
            f.write(f"probe{n} = c{i}{j}\n")
        for i, j in CELLS:
            x, y, z = centre(i, j)
            f.write(f"\n[c{i}{j}]\nvariables = rho(1) u p\nditer = 1\nposition = {x:.15E} {y:.15E} {z:.15E}\n")


def read_field(path):
    """Cell values of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = (NX + 1) * (NY + 1) * 2, NX * NY, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out


def check(path):
    field = read_field(path)
    rho = field.get("rho(1)", field.get("rho1"))
    fail = False
    for i, j in CELLS:
        c = (j - 1) * NX + (i - 1)
        # a row is the iteration and the three variables, written list-directed: it may take more than one line
        with open(f"OUTPUT/c{i}{j}.txt") as f:
            values = f.read().split()
        r_p, u_p, p_p = (float(v.replace('D', 'E')) for v in values[-3:])
        ok = abs(p_p - field["p"][c]) <= 1e-12 * abs(field["p"][c]) and abs(r_p - rho[c]) <= 1e-12 * abs(rho[c])
        fail |= not ok
        print(f"    probe at the centre of cell ({i}, {j}): p {p_p:.10e}, cell {field['p'][c]:.10e}; "
              f"rho {r_p:.10e}, cell {rho[c]:.10e}{'' if ok else '  <-- differs'}")
    sys.exit(1 if fail else 0)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case()
    elif what == "check":
        check(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
