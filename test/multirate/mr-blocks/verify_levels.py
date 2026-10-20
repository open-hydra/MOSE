#!/usr/bin/env python3
"""Levels across block connections: from the level map printed by MOSE (lines
MR-LEVEL of solve.log, trace = .true.) and the boundary file of the case:
  - every connection record (101, 201) has the same level as its partner cell (a connection is never a face
    between levels);
  - the levels of two cells sharing a face differ by at most one, inside every block and across every connection;
  - with --single RUN_1B (the same mesh in one block): the multi-block map is >= the single-block map everywhere (the equal
    levels across the connections and the nesting only raise levels) and the raised cells are listed (information).
Usage: verify_levels.py RUN CASE [--single RUN_1B]"""
import argparse
import os
import sys
from collections import Counter

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def levels(run):
    lev = {}
    for line in open(os.path.join(run, 'solve.log'), errors='replace'):
        s = line.split()
        if len(s) == 6 and s[0] == 'MR-LEVEL':
            b, i, j, k, l = map(int, s[1:6])
            lev[(b, i, j, k)] = l
    if not lev:
        raise SystemExit('no MR-LEVEL lines in %s' % run)
    return lev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('case')
    ap.add_argument('--single', default=None)
    a = ap.parse_args()
    lev = levels(a.run)
    recs = mc.read_bc(a.case)
    ok = True
    ncon, bad = 0, 0
    for r in recs:
        if r['t'] not in (101, 201):
            continue
        ncon += 1
        l0, l1 = lev[(r['b'], r['i'], r['j'], r['k'])], lev[(r['bs'], r['is_'], r['js'], r['ks'])]
        if l0 != l1:
            bad += 1
            if bad <= 5:
                print('  connection record b %d (%d,%d,%d) f %d: level %d, partner b %d (%d,%d,%d): level %d' % (
                    r['b'], r['i'], r['j'], r['k'], r['f'], l0, r['bs'], r['is_'], r['js'], r['ks'], l1))
    print('  connection records: %d, with different levels on the two sides: %d (0 required)' % (ncon, bad))
    ok = ok and bad == 0
    jumps = Counter()
    for (b, i, j, k), l in lev.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            q = (b, i + d[0], j + d[1], k + d[2])
            if q in lev:
                jumps[abs(lev[q] - l)] += 1
    print('  level jumps across the internal faces of the blocks: %s (<= 1 required)' % dict(sorted(jumps.items())))
    ok = ok and max(jumps) <= 1
    print('  cells per level: %s' % dict(sorted(Counter(lev.values()).items())))
    if a.single:
        # match the cells of the two runs by their centres (VTK output of both runs)
        mb = mc.by_centre(mc.read_run(a.run), 'p')
        sb = mc.by_centre(mc.read_run(a.single), 'p')
        if set(mb) != set(sb):
            print('  the two meshes differ')
            return 1
        cen = {}
        for b, blk in enumerate(mc.read_run(a.run), 1):
            C = mc.cell_key(blk['points'])
            for idx in np.ndindex(blk['cell']['p'].shape):
                cen[tuple(C[idx])] = (b, idx[2] + 1, idx[1] + 1, idx[0] + 1)
        cen1 = {}
        for blk in mc.read_run(a.single):
            C = mc.cell_key(blk['points'])
            for idx in np.ndindex(blk['cell']['p'].shape):
                cen1[tuple(C[idx])] = (1, idx[2] + 1, idx[1] + 1, idx[0] + 1)
        lev1 = levels(a.single)
        raised, lowered = [], 0
        for x, key in cen.items():
            lm, ls = lev[key], lev1[cen1[x]]
            if lm > ls:
                raised.append((key, ls, lm))
            elif lm < ls:
                lowered += 1
        print('  against the single block: %d cells raised, %d lowered (0 required)' % (len(raised), lowered))
        for key, ls, lm in sorted(raised)[:20]:
            print('    raised: block %d cell (%d,%d,%d): level %d -> %d' % (key + (ls, lm)))
        ok = ok and lowered == 0
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
