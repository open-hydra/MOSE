#!/usr/bin/env python3
"""
Build one case of the manufactured-solution verification on skewed and
non-smooth grids: the mesh (two blocks with different index orientations),
INPUT/ic.tec (the exact field at the cell centres) and INPUT/bc.txt.

Grid families on the unit square (2D, one cell across in z) or cube (3D),
N cells per side:
  zigzag  the nodes of every other i line moved by +-0.3 h in y (3D: also in z,
          alternating with j, and in x, alternating with k), the displacement
          weighted by 4 s (1 - s) so that the boundary stays straight: every
          cell is skewed and neighbouring cells are skewed the opposite way;
  iface   a cell-size jump of 1.5 at y = 1/2 (3N/5 cells below, 2N/5 above)
          and grid lines kinked there: x is moved by c(y) = 0.3 y below and
          0.3 (1 - y) above (3D: z by c/2), zero on the boundary;
  cart    the Cartesian grid.
The two blocks meet with permuted indices: zigzag, the second half turned by
180 degrees about k (3D: i and k exchanged); iface, the upper part turned so
that its j face meets an i face (3D: a k face, with the index directions
permuted and reversed).

Fields (the manufactured fields of the MOSE_mms driver, src/app/mms.f90):
  periodic2d, periodic3d  periodic on every side (types 201);
  wall2d, wall3d          periodic in x (and z), walls at y = 0 (type 302,
                          T = T_w) and y = 1 (type 301, the exact wall heat
                          flux).
In 2D the k faces are symmetry planes (type 300).

The gas constant must match the thermodynamics of the solver: Runiv of FLINT
over the molar mass of test/common/Air/phase.txt.

Usage: python3 make_case.py --fam zigzag|iface|cart --dim 2|3 --n N
                            --field periodic2d|periodic3d|wall2d|wall3d [--out DIR]
"""
import argparse
import math
import os

import numpy as np

RUNIV = 8314.46261815324        # J/(kmol K), FLINT (exact SI value)
WAIR = 28.970418                # kg/kmol, test/common/Air/phase.txt
RGAS = RUNIV / WAIR
CP = 1004.5
MU = 10.0
PR = 0.72
KCOND = MU * CP / PR
OM = 2.0 * math.pi
FMT = '{: .17E}'
BC_CONN, BC_PER = 101, 201


# --------------------------------------------------------------------------------------------------
# Exact fields (the same as src/app/mms.f90 and mms_gen_fields.py)
# --------------------------------------------------------------------------------------------------
def exact(field, x, y, z):
    """rho, u, v, w, p, T at the given points."""
    if field == 'periodic2d':
        rho = 1.0 + 0.1 * np.sin(OM * x) * np.cos(OM * y)
        u = 40.0 + 8.0 * np.sin(OM * x) * np.sin(OM * y)
        v = 30.0 + 8.0 * np.cos(OM * x) * np.cos(OM * y)
        w = np.zeros_like(x)
        p = 8.0e3 + 8.0e2 * np.cos(OM * x) * np.sin(OM * y)
    elif field == 'periodic3d':
        rho = 1.0 + 0.1 * np.sin(OM * x) * np.cos(OM * y) * np.cos(OM * z)
        u = 40.0 + 8.0 * np.sin(OM * x) * np.sin(OM * y) * np.cos(OM * z)
        v = 30.0 + 8.0 * np.cos(OM * x) * np.cos(OM * y) * np.cos(OM * z)
        w = 20.0 + 8.0 * np.cos(OM * x) * np.sin(OM * y) * np.sin(OM * z)
        p = 8000.0 + 800.0 * np.cos(OM * x) * np.sin(OM * y) * np.cos(OM * z)
    elif field in ('wall2d', 'wall3d'):
        cz = np.cos(OM * z) if field == 'wall3d' else 1.0
        sy = np.sin(math.pi * y)
        u = sy * (40.0 + 8.0 * np.sin(OM * x) * cz)
        v = 8.0 * sy * np.cos(OM * x) * cz
        w = sy * (20.0 + 8.0 * np.cos(OM * x) * np.sin(OM * z)) if field == 'wall3d' else np.zeros_like(x)
        T = 28.0 + 1.4 * np.cos(OM * x) * cz + 5.6 * sy
        p = 8000.0 + 800.0 * np.cos(OM * x) * np.cos(math.pi * y) * cz
        rho = p / (RGAS * T)
        return rho, u, v, w, p, T
    else:
        raise SystemExit('unknown field ' + field)
    return rho, u, v, w, p, p / (rho * RGAS)


def dTdy(field, x, y, z):
    """Exact dT/dy of the wall fields, T = T_w(x, z) + 5.6 sin(pi y)."""
    return 5.6 * math.pi * np.cos(math.pi * y) + 0.0 * x


# --------------------------------------------------------------------------------------------------
# Grid families: nodes (ni+1, nj+1, nk+1, 3) on [0,1]^d
# --------------------------------------------------------------------------------------------------
def nodes(fam, dim, n):
    s = np.linspace(0.0, 1.0, n + 1)
    sk = s if dim == 3 else np.array([0.0, 1.0 / n])
    XI, ETA, ZE = np.meshgrid(s, s, sk, indexing='ij')
    X, Y, Z = XI.copy(), ETA.copy(), ZE.copy()
    h = 1.0 / n
    if fam == 'zigzag':
        si = np.where(np.arange(n + 1) % 2 == 0, 1.0, -1.0)
        Y = ETA + 0.3 * h * si[:, None, None] * 4.0 * ETA * (1.0 - ETA)
        if dim == 3:
            Z = ZE + 0.3 * h * si[None, :, None] * 4.0 * ZE * (1.0 - ZE)
            X = XI + 0.3 * h * si[None, None, :] * 4.0 * XI * (1.0 - XI)
    elif fam == 'iface':
        if n % 5:
            raise SystemExit('iface needs N divisible by 5')
        nl = 3 * n // 5
        yv = np.concatenate([np.linspace(0.0, 0.5, nl + 1), np.linspace(0.5, 1.0, n - nl + 1)[1:]])
        Y = np.broadcast_to(yv[None, :, None], X.shape).copy()
        c = np.where(Y <= 0.5, 0.3 * Y, 0.3 * (1.0 - Y))
        X = XI + c
        if dim == 3:
            Z = ZE + 0.5 * c
    elif fam != 'cart':
        raise SystemExit('unknown family ' + fam)
    N = np.stack([X, Y, Z], axis=-1)
    if dim == 2:
        N[:, :, 1, 2] = 1.0 / n
        N[:, :, 0, 2] = 0.0
    return N


def rotate_ij(X):
    """(i, j) -> (j, n_i - i): the block stays right-handed."""
    return np.ascontiguousarray(np.transpose(X, (1, 0, 2, 3))[:, ::-1])


def rotate_ik(X):
    return np.ascontiguousarray(np.transpose(X, (2, 1, 0, 3))[:, :, ::-1])


def flip_ij(X):
    """i and j reversed together (a turn of 180 degrees about k)."""
    return np.ascontiguousarray(X[::-1, ::-1])


def blocks(fam, dim, N):
    ni, nj = N.shape[0] - 1, N.shape[1] - 1
    if fam == 'iface':
        jl = 3 * nj // 5
        low, up = N[:, :jl + 1], N[:, jl:]
        up = rotate_ij(up) if dim == 2 else rotate_ik(rotate_ij(up))
        return [np.ascontiguousarray(low), np.ascontiguousarray(up)]
    h = ni // 2
    a, b = N[:h + 1], N[h:]
    if fam == 'zigzag':
        b = flip_ij(b) if dim == 2 else rotate_ik(b)
    return [np.ascontiguousarray(a), np.ascontiguousarray(b)]


# --------------------------------------------------------------------------------------------------
# Block faces, connections (101) and periodic translations (201)
# --------------------------------------------------------------------------------------------------
def dims(X):
    return X.shape[0] - 1, X.shape[1] - 1, X.shape[2] - 1


def centres(X):
    """Cell centre = average of the 8 corner nodes, as in the driver."""
    return 0.125 * (X[:-1, :-1, :-1] + X[1:, :-1, :-1] + X[:-1, 1:, :-1] + X[1:, 1:, :-1] +
                    X[:-1, :-1, 1:] + X[1:, :-1, 1:] + X[:-1, 1:, 1:] + X[1:, 1:, 1:])


def face_nodes(X, face):
    """The four nodes of every cell face of block face `face` (1..6), shape (n1, n2, 4, 3)."""
    if face in (1, 2):
        P = X[0 if face == 1 else -1]
    elif face in (3, 4):
        P = X[:, 0 if face == 3 else -1]
    else:
        P = X[:, :, 0 if face == 5 else -1]
    return np.stack([P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]], axis=2)


def face_cell(X, face, t1, t2):
    """1-based indices of the boundary cell for the 0-based tangential indices (t1, t2) of the face."""
    ni, nj, nk = dims(X)
    if face in (1, 2):
        return (1 if face == 1 else ni), t1 + 1, t2 + 1
    if face in (3, 4):
        return t1 + 1, (1 if face == 3 else nj), t2 + 1
    return t1 + 1, t2 + 1, (1 if face == 5 else nk)


def face_order(X, face):
    """Record order on a face: k outermost, then j or i."""
    ni, nj, nk = dims(X)
    if face in (1, 2):
        return [(j, k) for k in range(nk) for j in range(nj)]
    if face in (3, 4):
        return [(i, k) for k in range(nk) for i in range(ni)]
    return [(i, j) for j in range(nj) for i in range(ni)]


def tangential(face):
    return {1: (1, 2), 2: (1, 2), 3: (0, 2), 4: (0, 2), 5: (0, 1), 6: (0, 1)}[face]


def _keys(F, scale, shift=None):
    Q = F if shift is None else F + shift
    R = np.round(Q / scale, 7)
    out = {}
    for a in range(Q.shape[0]):
        for b in range(Q.shape[1]):
            out[tuple(sorted(tuple(p) for p in R[a, b]))] = (a, b)
    return out


def connections(blks, scale, translations=(), skip=None):
    """{(b, face, i, j, k): (type, bs, fs, is, js, ks)} for every cell face shared by two cells (101) or tied by one of
    the translations (201; the source cell is the one at x + T)."""
    table, faces = {}, []
    for b, X in enumerate(blks):
        for face in range(1, 7):
            if skip is not None and skip(b, face):
                continue
            F = face_nodes(X, face)
            faces.append((b, face, F))
            for key, (a, c) in _keys(F, scale).items():
                table.setdefault(key, []).append((b, face, a, c))
    conn = {}
    for v in table.values():
        if len(v) == 2:
            (b1, f1, a1, c1), (b2, f2, a2, c2) = v
            r1 = (b1, f1) + face_cell(blks[b1], f1, a1, c1)
            r2 = (b2, f2) + face_cell(blks[b2], f2, a2, c2)
            conn[r1] = (BC_CONN,) + r2
            conn[r2] = (BC_CONN,) + r1
        elif len(v) > 2:
            raise ValueError('a cell face shared by more than two cells')
    for T in translations:
        T = np.asarray(T, float)
        for b, face, F in faces:
            for key, (a, c) in _keys(F, scale, shift=T).items():
                hit = table.get(key)
                if not hit or len(hit) != 1:
                    continue
                bs, fs, a2, c2 = hit[0]
                r1 = (b, face) + face_cell(blks[b], face, a, c)
                r2 = (bs, fs) + face_cell(blks[bs], fs, a2, c2)
                if r1 in conn or r2 == r1:
                    continue
                conn[r1] = (BC_PER,) + r2
                conn[r2] = (BC_PER,) + r1
    return conn


def orientation(blks, conn, rec):
    """d11, d12, d21, d22: step of the source cell along the two tangential directions of its face for a +1 step of the
    record cell along its first and second tangential direction."""
    b, face, i, j, k = rec
    _, bs, fs, is_, js, ks = conn[rec]
    out = []
    for t in tangential(face):
        idx = [i, j, k]
        lim = dims(blks[b])[t]
        if lim == 1:                      # one cell: the k direction of a 2D mesh maps onto itself
            out.append((0, 1) if fs in (1, 2, 3, 4) else (1, 0))
            continue
        step = 1 if idx[t] < lim else -1
        idx[t] += step
        nb = conn.get((b, face, idx[0], idx[1], idx[2]))
        if nb is None or nb[1] != bs or nb[2] != fs:
            raise ValueError('the interface is not continuous along a tangential direction: %s' % (rec,))
        ds = [nb[3] - is_, nb[4] - js, nb[5] - ks]
        ts = tangential(fs)
        out.append((step * ds[ts[0]], step * ds[ts[1]]))
    return out[0][0], out[0][1], out[1][0], out[1][1]


def face_centre(X, face, i, j, k):
    if face in (1, 2):
        a = i - 1 if face == 1 else i
        P = [X[a, j - 1, k - 1], X[a, j, k - 1], X[a, j, k], X[a, j - 1, k]]
    elif face in (3, 4):
        a = j - 1 if face == 3 else j
        P = [X[i - 1, a, k - 1], X[i, a, k - 1], X[i, a, k], X[i - 1, a, k]]
    else:
        a = k - 1 if face == 5 else k
        P = [X[i - 1, j - 1, a], X[i, j - 1, a], X[i, j, a], X[i - 1, j, a]]
    return 0.25 * (P[0] + P[1] + P[2] + P[3]), P


# --------------------------------------------------------------------------------------------------
# Files
# --------------------------------------------------------------------------------------------------
def write_ic(path, blks, names, values):
    with open(path, 'w') as out:
        out.write('VARIABLES = "x" "y" "z" ' + ' '.join('"%s"' % n for n in names) + '\n')
        nv = len(names)
        for b, X in enumerate(blks):
            ni, nj, nk = dims(X)
            out.write('ZONE T="Block%d", I=%d, J=%d, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED)\n'
                      % (b + 1, ni + 1, nj + 1, nk + 1, 3 + nv))
            for axis in range(3):
                col = np.transpose(X[..., axis], (2, 1, 0)).ravel()
                out.write(''.join(' ' + FMT.format(v) + '\n' for v in col))
            for arr in values[b]:
                col = np.transpose(np.asarray(arr, float), (2, 1, 0)).ravel()
                out.write(''.join(' ' + FMT.format(v) + '\n' for v in col))


def write_bc(path, blks, conn, outer):
    recs = []
    for b, X in enumerate(blks):
        for face in range(1, 7):
            for (t1, t2) in face_order(X, face):
                i, j, k = face_cell(X, face, t1, t2)
                rec = (b, face, i, j, k)
                head = '%8d%8d%8d%8d%8d' % (b + 1, i, j, k, face)
                if rec in conn:
                    typ, bs, fs, is_, js, ks = conn[rec]
                    d11, d12, d21, d22 = orientation(blks, conn, rec)
                    recs.append(head + '%8d\n' % typ +
                                '%8d%8d%8d%8d%8d%8d%8d%8d%8d\n' % (bs + 1, is_, js, ks, fs, d11, d12, d21, d22))
                else:
                    typ, data = outer(b, face, i, j, k)
                    recs.append(head + '%8d\n' % typ + (data + '\n' if data else ''))
    with open(path, 'w') as out:
        out.write(''.join(recs))
    return len(recs)


def main():
    ap = argparse.ArgumentParser(description='manufactured-solution case on a skewed grid family')
    ap.add_argument('--fam', required=True, choices=('zigzag', 'iface', 'cart'))
    ap.add_argument('--dim', type=int, required=True, choices=(2, 3))
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--field', required=True, choices=('periodic2d', 'periodic3d', 'wall2d', 'wall3d'))
    ap.add_argument('--out', default='.')
    a = ap.parse_args()
    walls = a.field.startswith('wall')
    B = blocks(a.fam, a.dim, nodes(a.fam, a.dim, a.n))
    T = [(1.0, 0.0, 0.0)]
    if not walls:
        T.append((0.0, 1.0, 0.0))
    if a.dim == 3:
        T.append((0.0, 0.0, 1.0))
    skip = (lambda b, face: face >= 5) if a.dim == 2 else None
    conn = connections(B, 1e-4, translations=T, skip=skip)

    def outer(b, face, i, j, k):
        if a.dim == 2 and face >= 5:
            return 300, None
        c, P = face_centre(B[b], face, i, j, k)
        P = np.array(P)
        if walls and np.all(np.abs(P[:, 1]) < 1e-12):
            T_w = exact(a.field, *c)[5]
            return 302, '    %.17E,    0.000000E+00,    0.000000E+00,' % float(T_w)
        if walls and np.all(np.abs(P[:, 1] - 1.0) < 1e-12):
            # heat flux of the 301 wall as MOSE applies it: k grad T . n_g, with n_g the face normal oriented towards
            # increasing index; on the plane wall y = 1 it is +e_y on the faces 2, 4, 6 and -e_y on 1, 3, 5
            ng = 1.0 if face % 2 == 0 else -1.0
            q = KCOND * float(dTdy(a.field, *c)) * ng
            return 301, '    %.17E,    0.000000E+00,    0.000000E+00,' % q
        raise ValueError('boundary face not classified: %s' % ((b, face, i, j, k),))

    os.makedirs(os.path.join(a.out, 'INPUT'), exist_ok=True)
    n = write_bc(os.path.join(a.out, 'INPUT', 'bc.txt'), B, conn, outer)
    vals = []
    for X in B:
        C = centres(X)
        rho, u, v, w, p, _ = exact(a.field, C[..., 0], C[..., 1], C[..., 2])
        vals.append([rho, u, v, w, p])
    write_ic(os.path.join(a.out, 'INPUT', 'ic.tec'), B, ['rho1', 'u', 'v', 'w', 'p'], vals)
    cells = int(sum(np.prod(dims(X)) for X in B))
    print('%s %dD N=%d %s: %d blocks, %d cells, %d boundary records' % (a.fam, a.dim, a.n, a.field, len(B), cells, n))


if __name__ == '__main__':
    main()
