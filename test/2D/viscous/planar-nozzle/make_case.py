#!/usr/bin/env python3
"""
Build one case of the planar-nozzle grid-family study: the mesh, INPUT/ic.tec
and INPUT/bc.txt (numpy only).

Geometry (units of the throat half-height R_t = 0.01 m): the half-height is
    h(x) = sqrt(1 + (A_i - 1) sin^2(pi x / (2 L_c))),   -L_c <= x <= 0,
    h(x) = sqrt(1 + (A_e - 1) sin^2(pi x / (2 L_d))),    0 <= x <= L_d,
with A_i = 2.25, A_e = 1.8, L_d = 4 and L_c = L_d sqrt((A_i - 1)/(A_e - 1)) = 5
(continuous second derivative at the throat): inlet and exit half-heights
1.5 and 1.342, zero wall slope at inlet and exit, isentropic exit Mach number
1.70; depth 0.1 (one cell).

Grid families, F = y/h from -1 (lower wall) to +1 (upper wall); every level
doubles the cells of the coarse one in both directions:
  O    grid lines normal to the walls and to the centre line,
       x = xi - h'(xi) h(xi) (F^2 - 1)/2;
  H    vertical grid lines, x = xi;
  ZZ   the O grid with every other row moved in x and every other column in F
       (+-0.15 of the local spacing) and a period-6 triangle in F (+-0.10):
       deterministic, non-smooth at every level, not symmetric;
  SB2  the ZZ grid plus random terms (grid lines meeting the walls at angles
       up to +-50 degrees, kinks, noise in x and F), drawn anew on every level
       from the seed: no two levels belong to one family.
Wall-normal spacing: the first two layers of the coarse grid 1e-4 h, then a
growth of 1.2 per coarse cell up to a cap; O and H have 48 cells from each
wall to the centre line on the coarse grid, ZZ and SB2 44 below and 56 above.
The grid is split into four blocks (convergent and divergent, lower and upper
half) with different index orientations.

Boundaries: inlet x = -L_c, ambient inlet at total conditions (type 407:
p0 = 1e5 Pa, T0 = 600 K, flow direction normal to the face, freestream
turbulence of the model); outlet x = L_d, supersonic, extrapolated (type 400);
walls, isothermal at T_w = 300 K (type 302); k faces, symmetry (type 300).
Initial field: the quasi-1D isentropic flow (subsonic branch in the
convergent, supersonic in the divergent), velocity along the grid lines of the
O family.

The gas constant must match the thermodynamics of the solver: Runiv of FLINT
over the molar mass of test/common/Air-viscous/phase.txt.

Usage: python3 make_case.py --fam O|H|ZZ|SB2 --level 0|1|2 --model SA|SST [--seed S] [--out DIR]
"""
import argparse
import math
import os

import numpy as np

RT = 0.01                       # throat half-height [m]
AI, AE = 2.25, 1.8              # inlet and exit area ratios
LD = 4.0
LC = LD * math.sqrt((AI - 1.0) / (AE - 1.0))
T0, P0, TW = 600.0, 1.0e5, 300.0
RUNIV = 8314.46261815324        # J/(kmol K), FLINT (exact SI value)
WAIR = 28.970418                # kg/kmol, test/common/Air-viscous/phase.txt
RGAS = RUNIV / WAIR
GAMMA = 1004.5 / (1004.5 - RGAS)
TMR_K, TMR_W = 9.0e-9, 1.0e-6   # freestream k / a0^2 and omega mu0 / (rho0 a0^2) of SST
DZ = 0.1                        # depth [R_t]
LEVELS = [(20, 16), (40, 32), (80, 64)]                    # cells along the convergent and the divergent
CELLS_F = {'O': (48, 48), 'H': (48, 48), 'ZZ': (44, 56), 'SB2': (44, 56)}   # cells below and above the centre line
A1, R1 = 1.0e-4, 1.2            # wall-normal map: first spacing / h, growth per coarse cell
SEED = 20261007
TH_MAX = math.radians(50.0)
L_FR, A_ZZ, A_RX, A_RF, A_KF, A_SL, A_ZF = 0.2, 0.15, 0.20, 0.20, 0.10, 0.20, 0.15
D_BL, C_OUT, K_OUT = 0.1, 3.0, 0.8
FMT = '{: .17E}'
BC_CONN = 101


def fmt(v):
    return '%.17E' % v


# --------------------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------------------
def half_height(x):
    x = np.asarray(x, float)
    a = np.where(x <= 0.0, 1.0 + (AI - 1.0) * np.sin(np.pi * x / (2 * LC)) ** 2,
                 1.0 + (AE - 1.0) * np.sin(np.pi * x / (2 * LD)) ** 2)
    return np.sqrt(a)


def dhalf_height(x):
    x = np.asarray(x, float)
    da = np.where(x <= 0.0, (AI - 1.0) * np.sin(np.pi * x / LC) * np.pi / (2 * LC),
                  (AE - 1.0) * np.sin(np.pi * x / LD) * np.pi / (2 * LD))
    return da / (2.0 * half_height(x))


def xi_axis(nc, nd):
    """Abscissae of the grid lines: uniform step, the same in the convergent and the divergent (Nc:Nd = Lc:Ld)."""
    return np.concatenate([np.linspace(-LC, 0.0, nc + 1), np.linspace(0.0, LD, nd + 1)[1:]])


def wall_map(n0, a=A1, r=R1):
    """Distance from the wall F_w(s) (units of h) for s in [0, 1] of a half with n0 cells on the coarse grid: spacing a for
    the first two layers, then growth r per coarse cell up to a cap, then constant; F_w(1) = 1."""
    lr = math.log(r)

    def F(t, k):
        t = np.asarray(t, float)
        grow = 2.0 * a + a * (r ** (np.minimum(t, 2.0 + k) - 2.0) - 1.0) / lr
        return np.where(t <= 2.0, a * t, grow + a * r ** k * np.maximum(t - 2.0 - k, 0.0))
    lo, hi = 0.0, n0 - 2.0
    if float(F(n0, hi)) < 1.0:
        raise ValueError('%d cells are not enough with a = %g and r = %g' % (n0, a, r))
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if float(F(n0, m)) < 1.0:
            lo = m
        else:
            hi = m
    k = 0.5 * (lo + hi)
    return lambda s: F(np.asarray(s, float) * n0, k) / float(F(n0, k))


def coordinates_F(fam, level):
    """F = y/h of the grid rows from -1 to +1, and the row of the centre line."""
    nl0, nu0 = CELLS_F[fam]
    nl, nu = nl0 * 2 ** level, nu0 * 2 ** level
    fl, fu = wall_map(nl0), wall_map(nu0)
    below = -(1.0 - fl(1.0 - np.linspace(0.0, 1.0, nl + 1)))[::-1]
    above = 1.0 - fu(1.0 - np.linspace(0.0, 1.0, nu + 1))
    below[0], above[-1], below[-1], above[0] = -1.0, 1.0, 0.0, 0.0
    return np.concatenate([below, above[1:]]), nl


def nodes_zz(fam, level, seed):
    """Nodes (X, Y) [R_t] of the families ZZ and SB2 and the centre-line row: the O grid with deterministic alternations (ZZ),
    plus random terms drawn from the seed (SB2)."""
    nc, nd = LEVELS[level]
    xs = xi_axis(nc, nd)
    F, nl = coordinates_F(fam, level)
    XS, FF = np.meshgrid(xs, F, indexing='ij')
    H = half_height(XS)
    X = XS - dhalf_height(XS) * H * (FF ** 2 - 1.0) / 2.0
    ni, nj = X.shape
    rng = np.random.default_rng(seed + 1000 * level)
    dx = np.gradient(X, axis=0)
    dF = np.minimum(np.abs(np.diff(FF, axis=1, prepend=np.nan)), np.abs(np.diff(FF, axis=1, append=np.nan)))
    dF = np.where(np.isnan(dF), np.nanmin(np.stack([np.abs(np.diff(FF, axis=1, prepend=np.nan)),
                                                    np.abs(np.diff(FF, axis=1, append=np.nan))]), axis=0), dF)
    upper = FF > 0
    dw = (1.0 - np.abs(FF)) * H
    ell = L_FR * np.abs(dx)
    # random draws, in a fixed order, made for both families so that SB2 does not depend on ZZ
    th = np.where(upper, rng.uniform(-TH_MAX, TH_MAX, ni)[:, None], rng.uniform(-TH_MAX, TH_MAX, ni)[:, None])
    sl = np.where(upper, rng.uniform(-A_SL, A_SL, ni)[:, None], rng.uniform(-A_SL, A_SL, ni)[:, None]) * dx
    tilt = np.tan(th) * np.where(dw < ell, dw, np.maximum(2.0 * ell - dw, 0.0))
    slide = sl * np.maximum(0.0, 1.0 - dw / (2.0 * ell))
    ramp = np.clip((dw - 2.0 * ell) / ell, 0.0, 1.0)
    dy = dF * H
    ax_ = np.minimum(np.abs(dx), 1.5 * dy)
    out = dw > D_BL * H
    ax_ = np.where(out, K_OUT * np.minimum(np.abs(dx), C_OUT * dy), ax_)
    kF = np.where(out, K_OUT, 1.0)
    zz = A_ZZ * ax_ * np.where(np.arange(nj)[None, :] % 2 == 0, 1.0, -1.0) * ramp
    rx = A_RX * ax_ * rng.uniform(-1.0, 1.0, X.shape) * ramp
    tri = np.abs(((np.arange(ni) % 6) - 3.0)) / 3.0 * 2.0 - 1.0
    kf = A_KF * dF * tri[:, None] * np.where(upper, 1.0, -1.0)
    kf = (kf + A_ZF * dF * np.where(np.arange(ni)[:, None] % 2 == 0, 1.0, -1.0)) * kF
    rf = A_RF * dF * rng.uniform(-1.0, 1.0, X.shape) * kF
    if fam == 'ZZ':
        tilt, slide, rx, rf = 0.0 * tilt, 0.0 * slide, 0.0 * rx, 0.0 * rf
    inner = np.ones_like(X, bool)
    inner[[0, -1], :] = False
    X = X + (tilt + slide + zz + rx) * inner
    FFn = FF + (kf + rf) * inner
    FFn[:, [0, -1]] = FF[:, [0, -1]]
    FFn[:, nl] = FF[:, nl] + (rf[:, nl]) * inner[:, nl]
    Y = half_height(X) * FFn
    return X, Y, nl


def rotate_ij(X):
    return np.ascontiguousarray(np.transpose(X, (1, 0, 2, 3))[:, ::-1])


def flip_ij(X):
    return np.ascontiguousarray(X[::-1, ::-1])


def grid(fam, level, seed=SEED):
    nc, nd = LEVELS[level]
    if fam in ('ZZ', 'SB2'):
        X, Y, nl = nodes_zz(fam, level, seed)
    else:
        xs = xi_axis(nc, nd)
        F, nl = coordinates_F(fam, level)
        XS, FF = np.meshgrid(xs, F, indexing='ij')
        H = half_height(XS)
        X = XS - dhalf_height(XS) * H * (FF ** 2 - 1.0) / 2.0 if fam == 'O' else XS.copy()
        Y = H * FF
    N = np.zeros(X.shape + (2, 3))
    for k in range(2):
        N[:, :, k, 0] = X * RT
        N[:, :, k, 1] = Y * RT
        N[:, :, k, 2] = k * DZ * RT
    blocks = [flip_ij(N[:nc + 1, :nl + 1]), N[:nc + 1, nl:], N[nc:, :nl + 1], rotate_ij(N[nc:, nl:])]
    return [np.ascontiguousarray(b) for b in blocks]


# --------------------------------------------------------------------------------------------------
# Flow: inflow turbulence and the quasi-1D initial field
# --------------------------------------------------------------------------------------------------
def mu_sutherland(T):
    return 1.458e-6 * T ** 1.5 / (T + 110.4)        # test/common/Air-sutherland/transport.dat


def inflow_turbulence(model):
    rho0 = P0 / (RGAS * T0)
    a0 = math.sqrt(GAMMA * RGAS * T0)
    mu0 = mu_sutherland(T0)
    if model == 'SA':
        return [3.0 * mu0 / rho0]
    return [TMR_K * a0 * a0, TMR_W * rho0 * a0 * a0 / mu0]


def mach_area(A, supersonic):
    """Mach number for the area ratio A/A* >= 1, supersonic branch where `supersonic` is true (bisection)."""
    g = GAMMA
    lo = np.where(supersonic, 1.0, 1e-9)
    hi = np.where(supersonic, 20.0, 1.0)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        a = (1 / m) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * m * m)) ** ((g + 1) / (2 * (g - 1)))
        high = (a > A) == supersonic
        hi = np.where(high, m, hi)
        lo = np.where(high, lo, m)
    return 0.5 * (lo + hi)


def quasi1d(C):
    """Quasi-1D isentropic state at the points C [m]: rho, u, v, w, p."""
    x = C[..., 0] / RT
    rr = C[..., 1] / RT
    R = half_height(x)
    M = mach_area(np.maximum(R, 1.0), x > 0.0)
    g = GAMMA
    T = T0 / (1.0 + 0.5 * (g - 1.0) * M * M)
    p = P0 * (T / T0) ** (g / (g - 1.0))
    rho = p / (RGAS * T)
    vel = M * np.sqrt(g * RGAS * T)
    vr = vel * (rr / R) * dhalf_height(x)
    return rho, vel, vr, np.zeros_like(vr), p


# --------------------------------------------------------------------------------------------------
# Connections and files (as test/2D/viscous/mms/skewed/make_case.py)
# --------------------------------------------------------------------------------------------------
def dims(X):
    return X.shape[0] - 1, X.shape[1] - 1, X.shape[2] - 1


def centres(X):
    return 0.125 * (X[:-1, :-1, :-1] + X[1:, :-1, :-1] + X[:-1, 1:, :-1] + X[1:, 1:, :-1] +
                    X[:-1, :-1, 1:] + X[1:, :-1, 1:] + X[:-1, 1:, 1:] + X[1:, 1:, 1:])


def face_nodes(X, face):
    if face in (1, 2):
        P = X[0 if face == 1 else -1]
    elif face in (3, 4):
        P = X[:, 0 if face == 3 else -1]
    else:
        P = X[:, :, 0 if face == 5 else -1]
    return np.stack([P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]], axis=2)


def face_cell(X, face, t1, t2):
    ni, nj, nk = dims(X)
    if face in (1, 2):
        return (1 if face == 1 else ni), t1 + 1, t2 + 1
    if face in (3, 4):
        return t1 + 1, (1 if face == 3 else nj), t2 + 1
    return t1 + 1, t2 + 1, (1 if face == 5 else nk)


def face_order(X, face):
    ni, nj, nk = dims(X)
    if face in (1, 2):
        return [(j, k) for k in range(nk) for j in range(nj)]
    if face in (3, 4):
        return [(i, k) for k in range(nk) for i in range(ni)]
    return [(i, j) for j in range(nj) for i in range(ni)]


def tangential(face):
    return {1: (1, 2), 2: (1, 2), 3: (0, 2), 4: (0, 2), 5: (0, 1), 6: (0, 1)}[face]


def connections(blks, scale):
    table = {}
    for b, X in enumerate(blks):
        for face in range(1, 5):
            F = face_nodes(X, face)
            R = np.round(F / scale, 7)
            for a in range(F.shape[0]):
                for c in range(F.shape[1]):
                    table.setdefault(tuple(sorted(tuple(p) for p in R[a, c])), []).append((b, face, a, c))
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
    return conn


def orientation(blks, conn, rec):
    b, face, i, j, k = rec
    _, bs, fs, is_, js, ks = conn[rec]
    out = []
    for t in tangential(face):
        idx = [i, j, k]
        lim = dims(blks[b])[t]
        if lim == 1:
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


def face_points(X, face, i, j, k):
    if face in (1, 2):
        a = i - 1 if face == 1 else i
        return [X[a, j - 1, k - 1], X[a, j, k - 1], X[a, j, k], X[a, j - 1, k]]
    if face in (3, 4):
        a = j - 1 if face == 3 else j
        return [X[i - 1, a, k - 1], X[i, a, k - 1], X[i, a, k], X[i - 1, a, k]]
    a = k - 1 if face == 5 else k
    return [X[i - 1, j - 1, a], X[i, j - 1, a], X[i, j, a], X[i - 1, j, a]]


def boundary(blks, model):
    tol = 1e-9 * RT
    turb = inflow_turbulence(model)

    def outer(b, face, i, j, k):
        if face >= 5:
            return 300, None
        P = np.array(face_points(blks[b], face, i, j, k))
        if np.all(np.abs(P[:, 0] + LC * RT) < tol):
            data = ', '.join([fmt(T0), fmt(P0), fmt(P0), 'normal', 'normal', fmt(1.0), fmt(1.0)] + [fmt(v) for v in turb]) + ','
            return 407, '    ' + data
        if np.all(np.abs(P[:, 0] - LD * RT) < tol):
            return 400, None
        hw = half_height(P[:, 0] / RT) * RT
        if np.all(np.abs(np.abs(P[:, 1]) - hw) < tol):
            return 302, '    ' + ', '.join([fmt(TW), fmt(0.0), fmt(0.0)]) + ','
        raise ValueError('boundary face not classified: block %d face %d cell %s' % (b + 1, face, (i, j, k)))
    return outer


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
    ap = argparse.ArgumentParser(description='planar nozzle on a grid family')
    ap.add_argument('--fam', required=True, choices=('O', 'H', 'ZZ', 'SB2'))
    ap.add_argument('--level', type=int, required=True, choices=(0, 1, 2))
    ap.add_argument('--model', required=True, choices=('SA', 'SST'))
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--out', default='.')
    a = ap.parse_args()
    B = grid(a.fam, a.level, a.seed)
    os.makedirs(os.path.join(a.out, 'INPUT'), exist_ok=True)
    conn = connections(B, RT * 1e-3)
    n = write_bc(os.path.join(a.out, 'INPUT', 'bc.txt'), B, conn, boundary(B, a.model))
    turb = inflow_turbulence(a.model)
    names = ['rho1', 'u', 'v', 'w', 'p'] + (['mi_t'] if a.model == 'SA' else ['kappa', 'omega'])
    vals = []
    for X in B:
        st = list(quasi1d(centres(X)))
        vals.append(st + [np.full(dims(X), v) for v in turb])
    write_ic(os.path.join(a.out, 'INPUT', 'ic.tec'), B, names, vals)
    cells = int(sum(np.prod(dims(X)) for X in B))
    print('%s level %d %s: %d blocks, %d cells, %d boundary records' % (a.fam, a.level, a.model, len(B), cells, n))


if __name__ == '__main__':
    main()
