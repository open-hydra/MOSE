#!/usr/bin/env python3
"""Halo of the multirate integrator built from the true stencils:
independent construction in Python, from the level map printed by MOSE (lines MR-LEVEL of the log, trace = .true.) and
from the boundary file of the case, compared with the lists printed by MOSE (lines MR-HALO of the log).
Definition, for every level l:
  S_l  = positions read by the faces of the level: c + s e_d, s = -2..2, d = 1..3 (convective faces) and, in viscous
         runs (equations = navier-stokes in the input.ini of the run), c + s e_d + t e_e, s, t = -1, 1, d < e
         (diagonals: tangential neighbours of the diffusive faces), for every cell c of level l (positions of the array
         with ghost cells, gc = 2); in viscous runs also the four tangential neighbours of the partner cell of every
         connection record of P_l (the positions copied into Pg(:,3:6) by Fill_BC_Ghost_Connection, Lib_Ghost.f90);
         in RANS runs (turbulence = SA or SST) also c + s e_d, s = -1, 1 (O_src: the central gradients of the
         turbulence sources, SA_Source_Blk and SST_Blk; inside the convective offsets, so it adds no position);
  P_l  = viscous runs: connection records (101/201) of the cells of level l;
  G_l  = boundary records whose ghost positions (depth 1 and 2 along the face normal) are in S_l;
  H_l  = cells of S_l of a coarser level, plus the source cells of the records of G_l of a coarser level (connections
         101/201: the partner cells at depth 1 and 2; other types: the boundary cell and the next 2 cells inward).
Viscous runs, at every stage of the first macro step (lines MR-VISC, one per line MR-TRACE): after the stage refreshed the
viscous variables of the cells it reads (blk%Visc), MOSE compares those that the viscous fluxes of the stage read with the
viscous variables of their primitive variables at that point; every stage must read values, and none stale.
Usage: verify_halo.py RUN_DIR CASE_DIR [--require-deep] [--require-diag]   (exit 0 if every H_l(m), G_l and P_l is
identical, and in viscous runs no stage reads a stale viscous variable)
--require-deep (case mr1d-deep1): for every level l >= 2 the halo must also hold cells of some
level m <= l - 2 (bands of one cell), so that the test exercises the stencil rule beyond the next coarser level.
--require-diag (viscous cases with convex corners in the level map): some halo cell must be reached only through
the diagonals of O_diff (the cells that a halo along the axes would miss)."""
import os
import sys
from collections import defaultdict

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import mrcheck as mc      # noqa: E402

GUIDE = {1: (1, 0, 0), 2: (-1, 0, 0), 3: (0, 1, 0), 4: (0, -1, 0), 5: (0, 0, 1), 6: (0, 0, -1)}
WITH_DATA = {101, 201, 103, 401, 402, 403, 404, 405, 406, 407, 408, 410, 420, 501, 502, 503, 504, 505, 506}


def read_bc(path):
    recs = []
    lines = open(path).read().split('\n')
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
            r.update(bs=q[0], is_=q[1], js=q[2], ks=q[3], fs=q[4], d11=q[5], d12=q[6], d21=q[7], d22=q[8])
        recs.append(r)
        k += 2 if t in WITH_DATA else 1
    return recs


def pg_sources(r):
    """positions of the partner block read by Fill_BC_Ghost_Connection for the connection record r"""
    i, j, k, a, b, c, e = r['is_'], r['js'], r['ks'], r['d11'], r['d12'], r['d21'], r['d22']
    if r['fs'] in (1, 2):
        pos = [(i, j - a, k - b), (i, j + a, k + b), (i, j - c, k - e), (i, j + c, k + e)]
    elif r['fs'] in (3, 4):
        pos = [(i - a, j, k - b), (i + a, j, k + b), (i - c, j, k - e), (i + c, j, k + e)]
    else:
        pos = [(i - a, j - b, k), (i + a, j + b, k), (i - c, j - e, k), (i + c, j + e, k)]
    return [(r['bs'],) + p for p in pos]


def physics_key(run, key):
    sec = None
    for line in open(os.path.join(run, 'input.ini')):
        s = line.strip()
        if s.startswith('['):
            sec = s
        elif sec == '[MOSE-Physics]' and '=' in s and s.split('=')[0].strip() == key:
            return s.split('=', 1)[1].split(';')[0].strip()
    return ''


def viscous_run(run):
    return physics_key(run, 'equations') == 'navier-stokes'


def rans_run(run):
    return viscous_run(run) and physics_key(run, 'turbulence') not in ('', 'none')


def main():
    run, case = sys.argv[1], sys.argv[2]
    require_deep = '--require-deep' in sys.argv[3:]
    require_diag = '--require-diag' in sys.argv[3:]
    lev = {}
    dims = {}
    for line in open(os.path.join(run, 'solve.log'), errors='replace'):
        s = line.split()
        if len(s) != 6 or s[0] != 'MR-LEVEL':
            continue
        b, i, j, k, l = map(int, s[1:6])
        lev[(b, i, j, k)] = l
        d = dims.setdefault(b, [0, 0, 0])
        d[0], d[1], d[2] = max(d[0], i), max(d[1], j), max(d[2], k)
    recs = read_bc(os.path.join(case, 'INPUT', 'bc.txt'))
    ghost = {}
    for n, r in enumerate(recs, 1):
        g0 = GUIDE[r['f']]
        for g in (1, 2):
            ghost[(r['b'], r['i'] - g0[0] * g, r['j'] - g0[1] * g, r['k'] - g0[2] * g)] = n
    L = max(lev.values())
    visc = viscous_run(run)
    rans = rans_run(run)
    print('  viscous run: %s, RANS run: %s' % (visc, rans))
    mine_H, mine_G, mine_P = defaultdict(set), defaultdict(set), defaultdict(set)
    diag_only = 0

    def halo(l, diagonals):
        S = set()
        for (b, i, j, k), lc in lev.items():
            if lc != l:
                continue
            for s in range(-2, 3):
                S.add((b, i + s, j, k)); S.add((b, i, j + s, k)); S.add((b, i, j, k + s))
            if diagonals:
                for s in (-1, 1):
                    for t in (-1, 1):
                        S.add((b, i + s, j + t, k)); S.add((b, i + s, j, k + t)); S.add((b, i, j + s, k + t))
            if rans:
                for s in (-1, 1):
                    S.add((b, i + s, j, k)); S.add((b, i, j + s, k)); S.add((b, i, j, k + s))
        P = set()
        if visc:
            P = set(n for n, r in enumerate(recs, 1) if r['t'] in (101, 201) and lev[(r['b'], r['i'], r['j'], r['k'])] == l)
            for n in P:
                S |= set(pg_sources(recs[n - 1]))
        G = set(ghost[p] for p in S if p not in lev and p in ghost)
        H = set(p for p in S if p in lev and lev[p] < l)
        for n in G:
            r = recs[n - 1]
            if r['t'] in (101, 201):
                g0 = GUIDE[r['fs']]
                src = [(r['bs'], r['is_'] + g0[0] * g, r['js'] + g0[1] * g, r['ks'] + g0[2] * g) for g in (0, 1)]
            else:
                g0 = GUIDE[r['f']]
                src = [(r['b'], r['i'] + g0[0] * g, r['j'] + g0[1] * g, r['k'] + g0[2] * g) for g in (0, 1, 2)]
            H |= set(p for p in src if p in lev and lev[p] < l)
        return H, G, P

    for l in range(L + 1):
        H, G, P = halo(l, visc)
        mine_G[l], mine_P[l] = G, P
        for p in H:
            mine_H[(l, lev[p])].add(p)
        if visc:
            Ha, Ga, _ = halo(l, False)
            print('  level %d: halo cells reached only through the diagonals %d, ghost records %d' % (
                l, len(H - Ha), len(G - Ga)))
            diag_only += len(H - Ha)
    mose_H, mose_G, mose_P = defaultdict(set), defaultdict(set), defaultdict(set)
    for line in open(os.path.join(run, 'solve.log'), errors='replace'):
        s = line.split()
        if len(s) < 2 or s[0] != 'MR-HALO' or s[1] == '#':
            continue
        s = s[1:]
        if s[0] == 'H':
            l, m, b, i, j, k = map(int, s[1:7])
            mose_H[(l, m)].add((b, i, j, k))
        elif s[0] in ('G', 'P'):
            l, n = int(s[1]), int(s[2])
            (mose_G if s[0] == 'G' else mose_P)[l].add(n)
            r = recs[n - 1]
            if (r['b'], r['i'], r['j'], r['k'], r['f'], r['t']) != tuple(map(int, s[3:9])):
                print('  record %d: MOSE %s, bc.txt %s' % (n, s[3:9], r))
                return 1
    ok = True
    for key in sorted(set(mine_H) | set(mose_H)):
        a, b = mine_H.get(key, set()), mose_H.get(key, set())
        print('  level %d, coarser level %d: halo %d cells (Python %d)%s' % (key[0], key[1], len(b), len(a),
                                                                          '' if a == b else '  DIFFERENT'))
        ok = ok and a == b
    for l in range(L + 1):
        a, b = mine_G.get(l, set()), mose_G.get(l, set())
        print('  level %d: ghost records %d (Python %d)%s' % (l, len(b), len(a), '' if a == b else '  DIFFERENT'))
        ok = ok and a == b
        a, b = mine_P.get(l, set()), mose_P.get(l, set())
        if visc or b:
            print('  level %d: connection records with Pg(:,3:6) %d (Python %d)%s' % (l, len(b), len(a),
                                                                                  '' if a == b else '  DIFFERENT'))
            ok = ok and a == b
    if require_diag:
        print('  halo cells reached only through the diagonals: %d (at least one required)' % diag_only)
        ok = ok and diag_only > 0
    if require_deep:
        for l in range(2, L + 1):
            n = sum(len(mose_H.get((l, m), set())) for m in range(l - 1))
            print('  level %d: halo cells of levels m <= %d: %d (at least one required)' % (l, l - 2, n))
            ok = ok and n > 0
    if visc:
        ok = Visc_Check(run) and ok
    return 0 if ok else 1


def Visc_Check(run):
    """viscous runs: one line MR-VISC per stage of the first macro step (as many as the lines MR-TRACE), each with values
    read and none stale"""
    stages, lines, bad = 0, [], []
    for line in open(os.path.join(run, 'solve.log'), errors='replace'):
        s = line.split()
        if not s:
            continue
        if s[0] == 'MR-TRACE':
            stages += 1
        elif s[0] == 'MR-VISC':
            f = dict(x.split('=', 1) for x in s[1:])
            lines.append(line.strip())
            if int(f['read']) == 0 or int(f['stale']) != 0:
                bad.append(line.strip())
    read = sum(int(dict(x.split('=', 1) for x in s.split()[1:])['read']) for s in lines)
    print('  viscous variables read by the fluxes: %d stages, %d checked (%d values read), %d with a stale value or none '
          'read' % (stages, len(lines), read, len(bad)))
    for x in bad[:5]:
        print('    ' + x)
    return stages > 0 and len(lines) == stages and not bad


if __name__ == '__main__':
    sys.exit(main())
