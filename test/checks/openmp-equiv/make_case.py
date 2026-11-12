"""Inputs of the chimera case of the OpenMPEquiv fast test.

Two blocks of square cells of side H, one cell thick in z: block 2, NB x NB cells, is the background; block 1, NA x NA cells,
lies on it shifted by OFF cells in x and in y, so that every cell of block 1 and every ghost cell of its first two layers
coincide with a cell of block 2. The four side faces of block 1 are chimera boundaries (102): each ghost cell takes the cell of
block 2 it coincides with (one donor, volume fraction 1, for the first and the second ghost layer). Faces 1 and 3, 1 and 4,
2 and 3, 2 and 4 of block 1 meet at its four edges, where the ghost cell of the edge belongs to both faces. The side faces of
block 2 extrapolate (400), the z faces of both blocks are null (0). Air of test/common/Air-viscous, a smooth state:
u = U0 + G (y - C), v = V0 + G (x - C), T = T0 [1 + 0.1 (x - C)/C], p = P0, C the centre of block 2.
The face lists of block 1 come first, face by face, so that the entries of two faces meeting at an edge are more than 64
entries apart: in the loops over the boundary entries (schedule(dynamic, 64)) they fall in different chunks, which several
threads can run at the same time.

Usage (from the test directory):
  python3 -B make_case.py chimera     INPUT/ic.tec and INPUT/bc.txt of the chimera case
"""
import os

H, NA, NB, OFF = 1.0e-3, 40, 48, 4
U0, V0, G = 20.0, 15.0, 1.0e3
T0, P0 = 300.0, 1.0e5
RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
C = 0.5 * NB * H


def state(x, y):
    t = T0 * (1.0 + 0.1 * (x - C) / C)
    return P0 / (RGAS * t), U0 + G * (y - C), V0 + G * (x - C), 0.0, P0


def zone(f, name, n, x0):
    f.write(f" ZONE  T = {name}, I={n + 1}, J={n + 1}, K=2, DATAPACKING=BLOCK, "
            "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
    for axis in range(3):
        for k in range(2):
            for j in range(n + 1):
                for i in range(n + 1):
                    f.write(f" {(x0 + H * i, x0 + H * j, H * k)[axis]:.17E}\n")
    cells = [state(x0 + H * (i + 0.5), x0 + H * (j + 0.5)) for j in range(n) for i in range(n)]
    for m in range(5):
        f.write("".join(f" {c[m]:.17E}\n" for c in cells))


def chimera_rows():
    """The four side faces of block 1, face by face: ghost layers 1 and 2 from the coincident cells of block 2."""
    rows = []
    for face in (1, 2, 3, 4):
        for m in range(1, NA + 1):
            i, j = {1: (1, m), 2: (NA, m), 3: (m, 1), 4: (m, NA)}[face]
            gi, gj = {1: (0, j), 2: (NA + 1, j), 3: (i, 0), 4: (i, NA + 1)}[face]
            di, dj = {1: (-1, 0), 2: (1, 0), 3: (0, -1), 4: (0, 1)}[face]
            rows.append(f"       1     {i:3d}     {j:3d}       1       {face}     102\n       1       1\n")
            rows.append(f"       2     {gi + OFF:3d}     {gj + OFF:3d}       1    0.1000000000E+01\n")
            rows.append(f"       2     {gi + di + OFF:3d}     {gj + dj + OFF:3d}       1    0.1000000000E+01\n")
    return rows


def write_chimera():
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        zone(f, "B1-IG", NA, OFF * H)
        zone(f, "B2-IG", NB, 0.0)
    rows = chimera_rows()
    for b, n in ((1, NA), (2, NB)):
        if b == 2:
            for face, idx in ((1, lambda m: (1, m)), (2, lambda m: (n, m)), (3, lambda m: (m, 1)), (4, lambda m: (m, n))):
                for m in range(1, n + 1):
                    i, j = idx(m)
                    rows.append(f"       {b}     {i:3d}     {j:3d}       1       {face}     400\n")
        for face in (5, 6):
            for j in range(1, n + 1):
                for i in range(1, n + 1):
                    rows.append(f"       {b}     {i:3d}     {j:3d}       1       {face}       0\n")
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(rows))


if __name__ == "__main__":
    write_chimera()
