"""Inputs of the ConnectionOrder tests.

A 2-D square, x and y in [0, 1] m, one cell thick in z, as two blocks of 16 x 8 square cells cut at y = 0.5:
  block 1, y in [0, 0.5]: i along x, j along y (the cut is its face 4, j = 8);
  block 2, y in [0.5, 1]: i along -y, j along x (the cut is its face 2, i = 8), so the index directions of the two blocks
           are permuted across the cut, as in meshes where blocks meet with different orientations.
The connection records (101) of the cut pair face 4 of block 1 with face 2 of block 2: two faces of the same parity, so
neither cell is on the low-index side of the face for both blocks, and the solver must still pick the record that computes
the flux of the face without looking at the order of the records in bc.txt. All other faces have the null type (0): no
flux, so no cell sums two boundary contributions whose order follows bc.txt.

Gas: air of test/common/Air-viscous. State: at rest, rho = 1 kg/m3 and p = 1e5 Pa below the cut, rho = 0.125 kg/m3 and
p = 1e4 (1 + 0.5 sin(pi x)) Pa above it: a shock leaves the cut at once, varying along x, so the shock sensor switches at
the cells on both sides of the cut and its value there differs from side to side.

Usage (from the test directory):
  python3 -B make_case.py case sorted|reversed      INPUT/ic.tec and INPUT/bc.txt, with the bc.txt records sorted by block
                                                   and face, or the same records in reverse order
"""
import os
import sys

NI, NJ = 16, 8                    # cells of block 1 along x and y; block 2 has NJ along -y and NI along x
H = 1.0 / NI                      # cell size, m
RGAS = 8314.46261815324 / 28.970418        # gas constant of FLINT for air, J/(kg K)
BC_NULL, BC_CONNECT = 0, 101
FMT = '{: .17E}'                  # 17 significant digits: exact round trip
SIN = __import__('math').sin
PI = __import__('math').pi


def state(x, y):
    """rho, u, v, w, p at (x, y)."""
    if y < 0.5:
        return 1.0, 0.0, 0.0, 0.0, 1.0e5
    return 0.125, 0.0, 0.0, 0.0, 1.0e4 * (1.0 + 0.5 * SIN(PI * x))


def node(b, i, j):
    """(x, y) of node (i, j) of block b."""
    if b == 1:
        return i * H, j * H
    return j * H, 1.0 - i * H


def write_ic(path):
    dims = {1: (NI, NJ), 2: (NJ, NI)}
    with open(path, 'w') as f:
        f.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"\n')
        for b in (1, 2):
            ni, nj = dims[b]
            f.write(f'ZONE T="Block{b}", I={ni + 1}, J={nj + 1}, K=2, '
                    'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED)\n')
            for axis in range(3):
                for k in range(2):
                    for j in range(nj + 1):
                        for i in range(ni + 1):
                            x, y = node(b, i, j)
                            f.write(' ' + FMT.format((x, y, k * H)[axis]) + '\n')
            cells = []
            for j in range(nj):
                for i in range(ni):
                    xs = [node(b, i + a, j + c) for a in (0, 1) for c in (0, 1)]
                    cells.append(state(sum(p[0] for p in xs) / 4.0, sum(p[1] for p in xs) / 4.0))
            for var in range(5):
                f.write(''.join(' ' + FMT.format(c[var]) + '\n' for c in cells))


def write_bc(path, reverse):
    records = []                                   # (block, i, j, k, face, type, data)
    for b, (ni, nj) in ((1, (NI, NJ)), (2, (NJ, NI))):
        for face in range(1, 7):
            if face in (1, 2):
                cells = [((1 if face == 1 else ni), j) for j in range(1, nj + 1)]
            elif face in (3, 4):
                cells = [(i, (1 if face == 3 else nj)) for i in range(1, ni + 1)]
            else:
                cells = [(i, j) for j in range(1, nj + 1) for i in range(1, ni + 1)]
            for i, j in cells:
                if b == 1 and face == 4:           # cut: block 1 cell (i, NJ) <-> block 2 cell (NJ, i), face 2
                    data = f'{2:8d}{NJ:8d}{i:8d}{1:8d}{2:8d}{1:8d}{0:8d}{0:8d}{1:8d}'
                    records.append((b, i, j, 1, face, BC_CONNECT, data))
                elif b == 2 and face == 2:         # cut: block 2 cell (NJ, j) <-> block 1 cell (j, NJ), face 4
                    data = f'{1:8d}{j:8d}{NJ:8d}{1:8d}{4:8d}{1:8d}{0:8d}{0:8d}{1:8d}'
                    records.append((b, i, j, 1, face, BC_CONNECT, data))
                else:
                    records.append((b, i, j, 1, face, BC_NULL, None))
    if reverse:
        records.reverse()
    with open(path, 'w') as f:
        for (b, i, j, k, face, t, data) in records:
            f.write(f'{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}{t:8d}\n')
            if data is not None:
                f.write(data + '\n')
    return len(records)


def main():
    if len(sys.argv) != 3 or sys.argv[1] != 'case' or sys.argv[2] not in ('sorted', 'reversed'):
        sys.exit('usage: make_case.py case sorted|reversed')
    os.makedirs('INPUT', exist_ok=True)
    write_ic(os.path.join('INPUT', 'ic.tec'))
    n = write_bc(os.path.join('INPUT', 'bc.txt'), sys.argv[2] == 'reversed')
    print(f'two blocks with permuted indices across the cut, {n} boundary records, {sys.argv[2]}')


if __name__ == '__main__':
    main()
