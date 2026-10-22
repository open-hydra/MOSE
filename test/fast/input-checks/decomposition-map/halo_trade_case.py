#!/usr/bin/env python3
"""Write a case whose blocks match halo-trade.decomposition.map, a decomposition map written by the ATLAS mesh decomposer
(MDB, test/MDB/halo-trade/reference/decomposition.map of ATLAS, kept byte for byte): its 6 blocks are 14, 13, 13, 12, 12
and 8 cells along i, 16 along j and 8 along k (three parent blocks of 40, 24 and 8 cells cut by MDB).

The blocks are put in a row along x, each connected to the next one (101, identity orientation), cubic cells of 25 mm;
all the outer faces are symmetry planes (300). Air at rest, 1.16 kg/m3, with a pressure wave of 1 % along x, so that the
flow changes from the first step on. Only the block sizes have to match the map: MOSE checks ni nj nk and cells.

Writes INPUT/ic.tec and INPUT/bc.txt.  Usage: python3 -B halo_trade_case.py
"""
import math
import os

DIMS = [(14, 16, 8), (13, 16, 8), (13, 16, 8), (12, 16, 8), (12, 16, 8), (8, 16, 8)]
H = 0.025                                   # cell edge [m]
RHO, P = 1.16, 1.0e5
LENGTH = H * sum(d[0] for d in DIMS)
FMT = '{: .17E}'
BC_CONNECT, BC_SYMMETRY = 101, 300


def write_ic(path):
    x0 = 0.0
    with open(path, 'w') as f:
        f.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"\n')
        for b, (ni, nj, nk) in enumerate(DIMS, start=1):
            f.write(f'ZONE T="Block{b}", I={ni + 1}, J={nj + 1}, K={nk + 1}, '
                    'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED)\n')
            for axis in range(3):
                vals = []
                for k in range(nk + 1):
                    for j in range(nj + 1):
                        for i in range(ni + 1):
                            vals.append((x0 + i * H, j * H, k * H)[axis])
                f.write(''.join(' ' + FMT.format(v) + '\n' for v in vals))
            cells = [(RHO, 0.0, 0.0, 0.0, P * (1.0 + 0.01 * math.sin(2.0 * math.pi * (x0 + (i + 0.5) * H) / LENGTH)))
                     for k in range(nk) for j in range(nj) for i in range(ni)]
            for var in range(5):
                f.write(''.join(' ' + FMT.format(c[var]) + '\n' for c in cells))
            x0 += ni * H


def write_bc(path):
    nb = len(DIMS)
    with open(path, 'w') as f:
        for b, (ni, nj, nk) in enumerate(DIMS, start=1):
            for face in range(1, 7):
                if face in (1, 2):
                    cells = [((1 if face == 1 else ni), j, k) for k in range(1, nk + 1) for j in range(1, nj + 1)]
                elif face in (3, 4):
                    cells = [(i, (1 if face == 3 else nj), k) for k in range(1, nk + 1) for i in range(1, ni + 1)]
                else:
                    cells = [(i, j, (1 if face == 5 else nk)) for j in range(1, nj + 1) for i in range(1, ni + 1)]
                for (i, j, k) in cells:
                    if face == 2 and b < nb:
                        f.write(f'{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}{BC_CONNECT:8d}\n')
                        f.write(f'{b + 1:8d}{1:8d}{j:8d}{k:8d}{1:8d}{1:8d}{0:8d}{0:8d}{1:8d}\n')
                    elif face == 1 and b > 1:
                        pni = DIMS[b - 2][0]
                        f.write(f'{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}{BC_CONNECT:8d}\n')
                        f.write(f'{b - 1:8d}{pni:8d}{j:8d}{k:8d}{2:8d}{1:8d}{0:8d}{0:8d}{1:8d}\n')
                    else:
                        f.write(f'{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}{BC_SYMMETRY:8d}\n')


def main():
    os.makedirs('INPUT', exist_ok=True)
    write_ic(os.path.join('INPUT', 'ic.tec'))
    write_bc(os.path.join('INPUT', 'bc.txt'))
    print(f'{len(DIMS)} blocks in a row, {sum(a * b * c for a, b, c in DIMS)} cells')


if __name__ == '__main__':
    main()
