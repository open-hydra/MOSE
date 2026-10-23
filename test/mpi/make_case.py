#!/usr/bin/env python3
"""Generate a small multi-block 2-D case for the MPI tests (test/mpi/).

The domain is a rectangle split into PI x PJ blocks of NI x NJ cells each,
one cell thick in z. Block (p, q) -- p along x, q along y -- is block number
b = 1 + p + PI*q, so a 2x2 layout puts four blocks around one common vertex
and a 4x1 layout is a strip of blocks with no such vertex.

    q=1  +-------+-------+
         |   3   |   4   |
    q=0  +-------+-------+        (2x2 layout)
         |   1   |   2   |
         +-------+-------+
            p=0     p=1

Outer boundaries (--bc):
    channel : inlet (401, T0 p0) on x-min, outlet (406) on x-max,
              adiabatic no-slip wall (301) on y-min, symmetry (300) on y-max
    box     : adiabatic no-slip walls (301) all around
The two z faces of every block carry --zfaces (300 symmetry or 0 unused).
Block cuts become pairs of connection records (101) with the identity
orientation (d11 d12 d21 d22 = 1 0 0 1), as ATLAS writes them.

Initial state: uniform p, T and a velocity U0 along x, plus a smooth vortex
centred on --centre (default: the middle of the domain, i.e. the common
vertex of a 2x2 layout), so the flow near the block cuts changes from the
first stage on.

bc.txt is written grouped by block and, within a block, by face (the order
the solver documents). --reverse writes the very same records in the reverse
order: a permutation of the input that must change neither the solution nor
the wall output.

--turbulence SA adds the Spalart-Allmaras variable rhoNuTilde = rho nu~
(3 times the molecular viscosity) to the initial state, and the specific
value nu~ = rhoNuTilde / rho to the data line of the inlet.

Usage:
    make_case.py --layout 4x1 --cells 16 16 --bc channel [--zfaces 300]
                 [--reverse] [--turbulence SA] [--centre X Y] [--outdir DIR]
Writes DIR/INPUT/ic.tec and DIR/INPUT/bc.txt.
"""

import argparse
import math
import os

# --- gas: single-species air of test/common/Air-viscous (cp = 1004.5) -------
R_GAS = 8314.46261815324 / 28.970418 # [J/(kg K)]
GAMMA = 1004.5 / (1004.5 - R_GAS)

# --- reference state and geometry -------------------------------------------
P_REF = 1.0e5                        # [Pa]
T_REF = 300.0                        # [K]
BLOCK_SIZE = 1.0e-3                  # block edge [m] (square blocks)
THICKNESS = 1.0e-4                   # z extent of the single cell layer [m]

# Boundary condition type codes (src/lib/io/IO_BC.f90)
BC_UNUSED, BC_CONNECT, BC_SYMMETRY, BC_WALL_Q = 0, 101, 300, 301
BC_INLET_T0P0, BC_OUTLET = 401, 406

FMT = '{: .17E}'                     # 17 significant digits: exact round trip
MU_TILDE = 3.0 * 1.8467e-5           # rho nu~ of --turbulence SA: 3 mu of air at 300 K [Pa s]


def parse():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--layout', default='2x2', help='PIxPJ blocks, e.g. 4x1')
    ap.add_argument('--cells', nargs=2, type=int, default=[16, 16],
                    metavar=('NI', 'NJ'), help='cells per block')
    ap.add_argument('--bc', choices=['channel', 'box'], default='box')
    ap.add_argument('--zfaces', type=int, choices=[BC_UNUSED, BC_SYMMETRY],
                    default=BC_SYMMETRY)
    ap.add_argument('--u0', type=float, default=None,
                    help='uniform x velocity [m/s] (channel: 100, box: 0)')
    ap.add_argument('--vortex', type=float, default=60.0,
                    help='peak vortex velocity [m/s]')
    ap.add_argument('--centre', nargs=2, type=float, default=None,
                    metavar=('X', 'Y'), help='vortex centre [m]')
    ap.add_argument('--reverse', action='store_true',
                    help='write the bc.txt records in reverse order')
    ap.add_argument('--turbulence', choices=['none', 'SA'], default='none',
                    help='add the Spalart-Allmaras variable')
    ap.add_argument('--outdir', default='.')
    a = ap.parse_args()
    a.PI, a.PJ = (int(s) for s in a.layout.lower().split('x'))
    if a.u0 is None:
        a.u0 = 100.0 if a.bc == 'channel' else 0.0
    return a


def nodes(n, x0):
    return [x0 + BLOCK_SIZE * k / n for k in range(n + 1)]


def state(x, y, a, xc, yc, radius):
    """Primitive state (rho, u, v, w, p) at (x, y)."""
    dx, dy = x - xc, y - yc
    g = a.vortex / radius * math.exp(0.5 * (1.0 - (dx * dx + dy * dy) / radius**2))
    u = a.u0 - g * dy
    v = g * dx
    rho = P_REF / (R_GAS * T_REF)
    return rho, u, v, 0.0, P_REF


def write_ic(path, a):
    ni, nj = a.cells
    xc, yc = a.centre if a.centre else (0.5 * a.PI * BLOCK_SIZE, 0.5 * a.PJ * BLOCK_SIZE)
    radius = 0.25 * BLOCK_SIZE
    with open(path, 'w') as f:
        sa = a.turbulence == 'SA'
        f.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"' + (' "rhoNuTilde"' if sa else '') + '\n')
        for q in range(a.PJ):
            for p in range(a.PI):
                b = 1 + p + a.PI * q
                xs, ys = nodes(ni, p * BLOCK_SIZE), nodes(nj, q * BLOCK_SIZE)
                zs = [0.0, THICKNESS]
                f.write(f'ZONE T="Block{b}", I={ni + 1}, J={nj + 1}, K=2, '
                        f'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-{9 if sa else 8}]=CELLCENTERED)\n')
                for axis in range(3):
                    for k in range(2):
                        for j in range(nj + 1):
                            for i in range(ni + 1):
                                v = (xs[i], ys[j], zs[k])[axis]
                                f.write(' ' + FMT.format(v) + '\n')
                cells = [state(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]),
                               a, xc, yc, radius)
                         for j in range(nj) for i in range(ni)]
                for var in range(5):
                    for c in cells:
                        f.write(' ' + FMT.format(c[var]) + '\n')
                if sa:
                    for c in cells:
                        f.write(' ' + FMT.format(MU_TILDE) + '\n')


def outer_bc(a, face):
    """(type, data line) of an outer boundary face."""
    wall = (BC_WALL_Q, '    0.000000E+00,    0.000000E+00,    0.000000E+00,')
    if a.bc == 'box':
        return wall
    m2 = (a.u0 / math.sqrt(GAMMA * R_GAS * T_REF))**2
    t0 = T_REF * (1.0 + 0.5 * (GAMMA - 1.0) * m2)
    p0 = P_REF * (t0 / T_REF)**(GAMMA / (GAMMA - 1.0))
    turb = f'    {MU_TILDE * R_GAS * T_REF / P_REF:.17E},' if a.turbulence == 'SA' else ''   # nu~ = rho nu~ / rho
    return {
        1: (BC_INLET_T0P0, f'    {t0:.17E},    {p0:.17E},         normal,         normal,'
                           f'    1.000000E+00,    1.000000E+00,' + turb),
        2: (BC_OUTLET, f'    {P_REF:.17E},    1.000000E+00,'),
        3: wall,
        4: (BC_SYMMETRY, None),
    }[face]


def write_bc(path, a):
    ni, nj, k = a.cells[0], a.cells[1], 1
    records = []                                   # (block, i, j, k, face, type, data)
    for q in range(a.PJ):
        for p in range(a.PI):
            b = 1 + p + a.PI * q
            per_face = []
            # face 1 (i = 1) and face 2 (i = ni): block cut or outer boundary
            for face, i, pn in ((1, 1, p - 1), (2, ni, p + 1)):
                for j in range(1, nj + 1):
                    if 0 <= pn < a.PI:
                        bs, is_, fs = 1 + pn + a.PI * q, (ni if face == 1 else 1), 3 - face
                        data = f'{bs:8d}{is_:8d}{j:8d}{k:8d}{fs:8d}{1:8d}{0:8d}{0:8d}{1:8d}'
                        per_face.append((b, i, j, k, face, BC_CONNECT, data))
                    else:
                        t, data = outer_bc(a, face)
                        per_face.append((b, i, j, k, face, t, data))
            # face 3 (j = 1) and face 4 (j = nj)
            for face, j, qn in ((3, 1, q - 1), (4, nj, q + 1)):
                for i in range(1, ni + 1):
                    if 0 <= qn < a.PJ:
                        bs, js, fs = 1 + p + a.PI * qn, (nj if face == 3 else 1), 7 - face
                        data = f'{bs:8d}{i:8d}{js:8d}{k:8d}{fs:8d}{1:8d}{0:8d}{0:8d}{1:8d}'
                        per_face.append((b, i, j, k, face, BC_CONNECT, data))
                    else:
                        t, data = outer_bc(a, face)
                        per_face.append((b, i, j, k, face, t, data))
            # faces 5 and 6: the inactive third direction
            for face in (5, 6):
                for j in range(1, nj + 1):
                    for i in range(1, ni + 1):
                        per_face.append((b, i, j, k, face, a.zfaces, None))
            records.extend(per_face)
    expected = a.PI * a.PJ * 2 * (nj + ni + ni * nj)
    assert len(records) == expected, (len(records), expected)
    if a.reverse:
        records.reverse()
    with open(path, 'w') as f:
        for (b, i, j, kk, face, t, data) in records:
            f.write(f'{b:8d}{i:8d}{j:8d}{kk:8d}{face:8d}{t:8d}\n')
            if data is not None:
                f.write(data + '\n')
    return len(records)


def main():
    a = parse()
    os.makedirs(os.path.join(a.outdir, 'INPUT'), exist_ok=True)
    write_ic(os.path.join(a.outdir, 'INPUT', 'ic.tec'), a)
    n = write_bc(os.path.join(a.outdir, 'INPUT', 'bc.txt'), a)
    print(f'{a.PI}x{a.PJ} blocks of {a.cells[0]}x{a.cells[1]} cells, bc={a.bc}, '
          f'zfaces={a.zfaces}, turbulence={a.turbulence}, '
          f'{n} boundary records{" (reversed)" if a.reverse else ""}')


if __name__ == '__main__':
    main()
