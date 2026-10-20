#!/usr/bin/env python3
"""Helpers for the multirate tests of MOSE (test/multirate/mr-*): read the VTK output, compute conserved totals, error
norms (whole domain and cells next to a level interface), observed orders, the exact Riemann solution of a perfect gas
and simple shape checks (local extrema, mirror symmetry). numpy only.

Conventions
- Output read: OUTPUT/field.vtm and OUTPUT/vtk/fieldBlock<b>.vts as written by MOSE with `sol-format = vtk raw`:
  raw appended data (`<AppendedData encoding="raw">`, each DataArray at its offset after the '_' marker, a UInt32 byte
  count then the little-endian values), the form that ORION writes; inline base64 DataArrays (UInt32 byte count, then
  the values) are read too.
- Gas: perfect gas of test/common/Air (cp = 1004.5 J/(kg K), M = 28.970418 kg/kmol), gamma = cp/(cp - R).
- Level map: lines "MR-LEVEL b i j k lev" printed by MOSE with [MOSE-Multirate] trace = .true. (first macro step) and
  read from the log of the run (solve.log); MOSE writes no file of its own for the multirate. A cell is an interface cell when one of its face neighbours inside the same
  block has a different level (interface_mask; width 4 for the interface error of the order tests, verify_order.py).
The exact Riemann solver iterates on the pressure of the star region (Newton, on the shock and rarefaction functions of
the two sides).
"""
import base64
import os
import re

import numpy as np

# Universal gas constant: the SOLVER's constant, read from this one place by every check of the tests. FLINT uses
# Runiv = 8314.46261815324 J/(kmol K), the exact SI value N_A k_B (lib/FLINT/src/lib/Lib_ThermoTransport.f90:10;
# R_i = Runiv/W_i in Load_ThermoTransport.f90:100). Another value, as 8314.51, gives a gamma
# different by 3.2e-6 and an internal energy different by a factor 1 + 8.0e-6 from the solver's: an apparent energy
# "drift" of 8e-6 * dKE/E in the checks of conservation.
RU = 8314.46261815324
CP_AIR = 1004.5
M_AIR = 28.970418
R_AIR = RU / M_AIR
GAMMA_AIR = CP_AIR / (CP_AIR - R_AIR)

_APPENDED = b'<AppendedData encoding="raw">'
_DATA_ARRAY = re.compile(rb'<DataArray([^>]*?)/?>')
_DTYPES = {'Float64': '<f8', 'Float32': '<f4', 'Int32': '<i4', 'Int64': '<i8', 'UInt8': 'u1'}


def _attributes(tag):
    """Attributes of a DataArray tag; names are quoted once (Name="rho(1)") or twice."""
    a = dict(re.findall(r'(\w+)="([^"]*)"', tag))
    m = re.search(r'Name=("+)([^"]*)\1', tag)
    if m:
        a['Name'] = m.group(2)
    return a


def read_vts(path):
    """One block of a MOSE vtk output: {'cell': {name: (nk, nj, ni) array}, 'points': (nk+1, nj+1, ni+1, 3),
    'time': TIME, 'shape': (nk, nj, ni)}. Raw appended and inline base64 encodings (module docstring)."""
    txt = open(path, 'rb').read()
    p = txt.find(_APPENDED)
    head = txt[:p] if p >= 0 else txt
    start = txt.index(b'_', p + len(_APPENDED)) + 1 if p >= 0 else None
    mh = re.search(rb'header_type="(\w+)"', head)
    hdt = '<u8' if (mh and mh.group(1) == b'UInt64') else '<u4'
    hsz = np.dtype(hdt).itemsize
    e = [int(x) for x in re.search(rb'WholeExtent="([^"]*)"', head).group(1).split()]
    ni, nj, nk = e[1] - e[0], e[3] - e[2], e[5] - e[4]
    out = {'cell': {}, 'shape': (nk, nj, ni)}
    cpos = head.find(b'<CellData>')
    for m in _DATA_ARRAY.finditer(head):
        a = _attributes(m.group(1).decode('latin-1'))
        dt = np.dtype(_DTYPES.get(a.get('type'), '<f8'))
        if a.get('format') == 'appended':
            off = start + int(a['offset'])
            nbytes = int(np.frombuffer(txt, dtype=hdt, count=1, offset=off)[0])
            arr = np.frombuffer(txt, dtype=dt, count=nbytes // dt.itemsize, offset=off + hsz)
        elif a.get('format') == 'binary':
            end = head.index(b'</DataArray>', m.end())
            raw = base64.b64decode(b''.join(head[m.end():end].split()))
            nbytes = int(np.frombuffer(raw[:hsz], dtype=hdt)[0])
            arr = np.frombuffer(raw[hsz:hsz + nbytes], dtype=dt)
        else:
            continue
        name = a.get('Name', '')
        if name == 'TIME':
            out['time'] = float(arr[0])
        elif name == 'Points':
            out['points'] = arr.astype(np.float64).reshape(nk + 1, nj + 1, ni + 1, 3)
        elif m.start() > cpos >= 0:
            out['cell'][name] = arr.reshape(nk, nj, ni)
    return out


def read_run(run_dir):
    """Blocks of the final field of a run (list of dicts, block order of field.vtm)."""
    vtm = os.path.join(run_dir, 'OUTPUT', 'field.vtm')
    files = re.findall(r'file="([^"]*)"', open(vtm).read())
    return [read_vts(os.path.join(run_dir, 'OUTPUT', 'vtk', os.path.basename(f))) for f in files]


def volumes(points):
    """Cell volumes of a structured block from its nodes (nk+1, nj+1, ni+1, 3): each hexahedron split into 6
    tetrahedra along the 000-111 diagonal (the volume of the hexahedron when its faces are planar; otherwise of the
    polyhedron whose faces are these triangles; the conservation checks run on planar faces)."""
    P = points
    c = [P[:-1, :-1, :-1], P[:-1, :-1, 1:], P[:-1, 1:, :-1], P[:-1, 1:, 1:],
         P[1:, :-1, :-1], P[1:, :-1, 1:], P[1:, 1:, :-1], P[1:, 1:, 1:]]
    # vertex index v = 4 k + 2 j + i ; diagonal from v0 (000) to v7 (111)
    tets = [(0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7)]
    vol = np.zeros(c[0].shape[:3])
    for a, b, cc, d in tets:
        A, B, C, D = c[a], c[b], c[cc], c[d]
        vol += np.einsum('...i,...i->...', B - A, np.cross(C - A, D - A)) / 6.0
    return np.abs(vol)


def partial_densities(cell):
    """The partial densities of a block, rho(1), rho(2), ..., in the order of phase.txt."""
    keys = [k for k in cell if re.fullmatch(r'rho\(\d+\)', k)]
    return [cell[k] for k in sorted(keys, key=lambda k: int(k[4:-1]))]


def species_mass(blocks, n):
    """Total mass of the species n (1, 2, ...), the n-th partial density in the order of phase.txt."""
    return sum(np.sum(partial_densities(b['cell'])[n - 1] * volumes(b['points'])) for b in blocks)


def totals(blocks, gamma=GAMMA_AIR):
    """Conserved totals: mass, x/y/z momentum, total energy of a perfect gas (with the tables of several species the
    energy of the last entry does not hold)."""
    tot = np.zeros(5)
    for b in blocks:
        V = volumes(b['points'])
        c = b['cell']
        rho, u, v, w, p = sum(partial_densities(c)), c['u'], c['v'], c['w'], c['p']
        tot += [np.sum(rho * V), np.sum(rho * u * V), np.sum(rho * v * V), np.sum(rho * w * V),
                np.sum((p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w)) * V)]
    return tot


def read_levels(run_dir, shapes, log='solve.log'):
    """Level map per block from the "MR-LEVEL b i j k lev" lines of the log of a run with trace = .true.;
    shapes = list of (nk, nj, ni)."""
    lev = [np.full(s, -1, dtype=int) for s in shapes]
    n = 0
    for line in open(os.path.join(run_dir, log), errors='replace'):
        s = line.split()
        if len(s) == 6 and s[0] == 'MR-LEVEL':
            b, i, j, k, l = (int(x) for x in s[1:6])
            lev[b - 1][k - 1, j - 1, i - 1] = l
            n += 1
    if n == 0:
        raise ValueError('no MR-LEVEL lines in %s (run with [MOSE-Multirate] trace = .true.)' % os.path.join(run_dir, log))
    return lev


def interface_mask(lev, width=1):
    """True on the `width` cells on each side (along the axes, same block) of every face between two levels:
    width = 1 gives the cells next to a level interface, width = 4 the mask of the interface error of the order tests."""
    m = np.zeros(lev.shape, dtype=bool)
    for ax in range(3):
        n = lev.shape[ax]
        if n < 2:
            continue
        a = np.moveaxis(lev, ax, 0)
        mm = np.moveaxis(m, ax, 0)
        d = a[1:] != a[:-1]                    # face p between cells p and p + 1
        for s in range(min(width, n - 1)):
            mm[:n - 1 - s] |= d[s:]            # cell p - s on the left of face p
            mm[1 + s:] |= d[:n - 1 - s]        # cell p + 1 + s on the right of face p
    return m


def errors(blocks, ref_blocks, var, mask_blocks=None):
    """L1 (volume weighted, normalised by the total volume) and Linf of blocks[var] - ref[var]; with masks, on the
    masked cells only."""
    num, vol, inf = 0.0, 0.0, 0.0
    for k, (b, r) in enumerate(zip(blocks, ref_blocks)):
        V = volumes(b['points'])
        e = np.abs(b['cell'][var] - r['cell'][var])
        if mask_blocks is not None:
            mk = mask_blocks[k]
            V, e = V[mk], e[mk]
        if e.size:
            num += np.sum(e * V)
            vol += np.sum(V)
            inf = max(inf, float(e.max()))
    return (num / vol if vol else 0.0), inf


def orders(errs):
    """Observed orders between successive halvings: log2(e_k / e_{k+1})."""
    e = np.asarray(errs, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log2(e[:-1] / e[1:])


def riemann_exact(xi, WL, WR, g=GAMMA_AIR, tol=1e-14, maxit=100):
    """Exact solution of the Riemann problem with primitive states WL = (rho, u, p), WR at similarity
    coordinate xi = (x - x_d)/t; returns (rho, u, p)."""
    rL, uL, pL = [float(v) for v in WL]
    rR, uR, pR = [float(v) for v in WR]
    cL = np.sqrt(g * pL / rL)
    cR = np.sqrt(g * pR / rR)

    def fK(p, rK, pK, cK):
        if p > pK:
            A = 2.0 / ((g + 1.0) * rK)
            B = (g - 1.0) / (g + 1.0) * pK
            sq = np.sqrt(A / (p + B))
            return (p - pK) * sq, sq * (1.0 - 0.5 * (p - pK) / (B + p))
        pr = p / pK
        return (2.0 * cK / (g - 1.0) * (pr ** ((g - 1.0) / (2.0 * g)) - 1.0),
                1.0 / (rK * cK) * pr ** (-(g + 1.0) / (2.0 * g)))

    z = (g - 1.0) / (2.0 * g)
    p = ((cL + cR - 0.5 * (g - 1.0) * (uR - uL)) / (cL / pL ** z + cR / pR ** z)) ** (1.0 / z)
    p = max(p, 1e-12)
    # measured: at most 15 iterations for pressure ratios 1e-3 to 1e3, density ratios 1e-2 to 1e2 and velocity jumps
    # within 2 c; outside that range a run without convergence stops here
    for _ in range(maxit):
        fl, dfl = fK(p, rL, pL, cL)
        fr, dfr = fK(p, rR, pR, cR)
        pn = max(p - (fl + fr + uR - uL) / (dfl + dfr), 1e-14)
        if abs(pn - p) <= tol * 0.5 * (pn + p):
            p = pn
            break
        p = pn
    else:
        raise RuntimeError('riemann_exact: Newton on p did not converge in %d iterations' % maxit)
    fl, _ = fK(p, rL, pL, cL)
    fr, _ = fK(p, rR, pR, cR)
    us, ps = 0.5 * (uL + uR) + 0.5 * (fr - fl), p
    xi = np.asarray(xi, dtype=float)
    rho, u, pp = np.empty_like(xi), np.empty_like(xi), np.empty_like(xi)
    gm = (g - 1.0) / (g + 1.0)
    for side in ('L', 'R'):
        if side == 'L':
            mask, (rK, uK, pK, cK, sgn) = xi <= us, (rL, uL, pL, cL, 1.0)
        else:
            mask, (rK, uK, pK, cK, sgn) = xi > us, (rR, uR, pR, cR, -1.0)
        x = xi[mask]
        r_, u_, p_ = np.empty_like(x), np.empty_like(x), np.empty_like(x)
        if ps > pK:
            S = uK - sgn * cK * np.sqrt((g + 1.0) / (2.0 * g) * ps / pK + (g - 1.0) / (2.0 * g))
            rs = rK * (ps / pK + gm) / (gm * ps / pK + 1.0)
            out = (x <= S) if sgn > 0 else (x >= S)
            r_[out], u_[out], p_[out] = rK, uK, pK
            r_[~out], u_[~out], p_[~out] = rs, us, ps
        else:
            cs = cK * (ps / pK) ** ((g - 1.0) / (2.0 * g))
            rs = rK * (ps / pK) ** (1.0 / g)
            SH, ST = uK - sgn * cK, us - sgn * cs
            head = (x <= SH) if sgn > 0 else (x >= SH)
            tail = (x >= ST) if sgn > 0 else (x <= ST)
            fan = ~(head | tail)
            r_[head], u_[head], p_[head] = rK, uK, pK
            r_[tail], u_[tail], p_[tail] = rs, us, ps
            xf = x[fan]
            uf = 2.0 / (g + 1.0) * (sgn * cK + 0.5 * (g - 1.0) * uK + xf)
            cf = 2.0 / (g + 1.0) * (cK + sgn * 0.5 * (g - 1.0) * (uK - xf))
            r_[fan] = rK * (cf / cK) ** (2.0 / (g - 1.0))
            u_[fan] = uf
            p_[fan] = pK * (cf / cK) ** (2.0 * g / (g - 1.0))
        rho[mask], u[mask], pp[mask] = r_, u_, p_
    return rho, u, pp


def local_extrema(q, rel_tol=1e-4):
    """Number of strict interior local extrema of a 1D array, ignoring oscillations below rel_tol of the range."""
    q = np.asarray(q, dtype=float)
    tol = rel_tol * (q.max() - q.min() + 1e-300)
    d = np.diff(q)
    d[np.abs(d) <= tol] = 0.0
    s = np.sign(d[d != 0.0])
    return int(np.sum(s[1:] != s[:-1]))


def cell_centres(points):
    P = points
    return 0.125 * (P[:-1, :-1, :-1] + P[:-1, :-1, 1:] + P[:-1, 1:, :-1] + P[:-1, 1:, 1:]
                    + P[1:, :-1, :-1] + P[1:, :-1, 1:] + P[1:, 1:, :-1] + P[1:, 1:, 1:])


def read_ic_tec(path):
    """Initial condition written by gen_mr_cases.py (Tecplot ASCII, BLOCK packing, nodal x y z then cell-centred
    variables, one value per line): list of blocks with 'points' and 'cell' like read_run (names as in the VTK output:
    rho(1), rho(2), ..., u, v, w, p)."""
    lines = open(path).read().split('\n')
    names = re.findall(r'"([^"]+)"', lines[0])
    blocks = []
    k = 1
    while k < len(lines):
        ln = lines[k]
        if ln.strip().startswith('ZONE'):
            I = int(re.search(r'I=(\d+)', ln).group(1))
            J = int(re.search(r'J=(\d+)', ln).group(1))
            K = int(re.search(r'K=(\d+)', ln).group(1))
            nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
            vals = np.array([float(x) for x in lines[k + 1:k + 1 + 3 * nn + (len(names) - 3) * nc]])
            P = np.stack([vals[0:nn], vals[nn:2 * nn], vals[2 * nn:3 * nn]], axis=-1).reshape(K, J, I, 3)
            cell = {}
            off = 3 * nn
            for nm in names[3:]:
                key = re.sub(r'^rho(\d+)$', r'rho(\1)', nm)
                cell[key] = vals[off:off + nc].reshape(K - 1, J - 1, I - 1)
                off += nc
            blocks.append({'points': P, 'cell': cell})
            k += 1 + 3 * nn + (len(names) - 3) * nc
        else:
            k += 1
    return blocks


# ---- several blocks ------------------------------------------------------------------------------------------------
GUIDE = {1: (1, 0, 0), 2: (-1, 0, 0), 3: (0, 1, 0), 4: (0, -1, 0), 5: (0, 0, 1), 6: (0, 0, -1)}  # inward normal of face f
_WITH_DATA = {101, 201, 103, 401, 402, 403, 404, 405, 406, 407, 408, 410, 420, 501, 502, 503, 504, 505, 506}


def read_bc(case_dir):
    """Records of INPUT/bc.txt of a case written by gen_mr_cases.py: list of dicts b, i, j, k, f, t (type) and, for the
    connections 101/201, bs, is_, js, ks, fs (partner cell and face)."""
    recs = []
    lines = open(os.path.join(case_dir, 'INPUT', 'bc.txt')).read().split('\n')
    k = 0
    while k < len(lines):
        s = lines[k].split()
        if len(s) != 6:
            k += 1
            continue
        b, i, j, kk, f, t = map(int, s)
        r = {'b': b, 'i': i, 'j': j, 'k': kk, 'f': f, 't': t}
        if t in (101, 201):
            q = list(map(int, lines[k + 1].split()))
            r.update(bs=q[0], is_=q[1], js=q[2], ks=q[3], fs=q[4])
        recs.append(r)
        k += 2 if t in _WITH_DATA else 1
    return recs


def connection_mask(recs, shapes, width):
    """True on the `width` cells inward from every connection face (101, 201), on both sides (each side has its own
    record); shapes = list of (nk, nj, ni). Added to the interface mask of the order tests of cases with several
    blocks (the cells next to a connection)."""
    m = [np.zeros(s, dtype=bool) for s in shapes]
    for r in recs:
        if r['t'] not in (101, 201):
            continue
        g = GUIDE[r['f']]
        nk, nj, ni = shapes[r['b'] - 1]
        for d in range(width):
            i, j, k = r['i'] + g[0] * d, r['j'] + g[1] * d, r['k'] + g[2] * d
            if 1 <= i <= ni and 1 <= j <= nj and 1 <= k <= nk:
                m[r['b'] - 1][k - 1, j - 1, i - 1] = True
    return m


def cell_key(points, ndig=9):
    """Cell centres of a block rounded to ndig decimals (metres): keys to match the cells of the same mesh cut into
    different blocks."""
    C = cell_centres(points)
    return np.round(C, ndig)


def by_centre(blocks, var):
    """dict cell centre (rounded) -> value of var, over all blocks of a field."""
    out = {}
    for b in blocks:
        C = cell_key(b['points'])
        q = b['cell'][var]
        for idx in np.ndindex(q.shape):
            out[tuple(C[idx])] = q[idx]
    return out
