#!/usr/bin/env python3
"""Observed temporal order of the multirate scheme.
Runs with macro steps H0, H0/2, H0/4, H0/8 against a reference run (single rate with a tiny fixed step, same mesh, same
final time). For every variable:
  e1_k   volume-weighted L1 error over all cells;
  einf_k Linf error over all cells;
  eint_k "interface" error with a floor: max(Linf over the cells within 4 cells of a face between two levels,
         0.05 x einf_k), so that an interface error far below the global one does not decide on round-off;
  the three ratios log2(e_k / e_k+1) are reported; the LAST pair decides: order_L1 >= min-l1, order_int >= min-int.
The 1-cell interface Linf is reported as information; with --min-int1 X it decides too, with the same floor:
max(Linf over the cells next to a face between two levels, 0.05 x einf_k), last pair >= X (the cells where the faces
between levels read their coarser neighbours, diagonals included).
Usage: verify_order.py --runs R0 R1 R2 R3 --ref RREF --levels RUN_WITH_TRACE [--vars rho(1),u,p]
(RUN_WITH_TRACE: a run with [MOSE-Multirate] trace = .true., whose log solve.log holds the MR-LEVEL lines of the map)
                       [--min-l1 1.9] [--min-int 1.8]
Cases of several blocks: --case CASE_DIR adds to the interface mask the cells next to the connections (4 cells
inward from every 101/201 face, both sides; 1 cell for the 1-cell mask). Boundary data: --inlet N
replaces the interface mask by the first N cells of block 1 along i (the inlet cells), with the same floor
(Linlet; Linlet0 is the inlet Linf without the floor, reported as information);
--l1-informative leaves the L1 order out of the decision (criteria that judge only the masked Linf).
--json FILE writes the errors and the orders of every variable to FILE, for the figures of the order test (plot_order.py).
"""
import argparse
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def fmt(a, f):
    return ' '.join(f % x for x in a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--levels', required=True)
    ap.add_argument('--vars', default='rho(1),u,p')
    ap.add_argument('--min-l1', type=float, default=1.9)
    ap.add_argument('--min-int', type=float, default=1.8)
    ap.add_argument('--case', default=None)
    ap.add_argument('--inlet', type=int, default=0)
    ap.add_argument('--l1-informative', action='store_true')
    ap.add_argument('--min-int1', type=float, default=None)
    ap.add_argument('--json', default=None)
    a = ap.parse_args()
    ref = mc.read_run(a.ref)
    runs = [mc.read_run(r) for r in a.runs]
    print('  final times: ' + ' '.join('%.15e' % b[0].get('time') for b in runs) + ' | ref %.15e' % ref[0].get('time'))
    shapes = [b['cell']['p'].shape for b in ref]
    lev = mc.read_levels(a.levels, shapes)
    mask4 = [mc.interface_mask(x, 4) for x in lev]
    mask1 = [mc.interface_mask(x, 1) for x in lev]
    what = 'interface mask'
    if a.case:
        recs = mc.read_bc(a.case)
        mask4 = [m | c for m, c in zip(mask4, mc.connection_mask(recs, shapes, 4))]
        mask1 = [m | c for m, c in zip(mask1, mc.connection_mask(recs, shapes, 1))]
        what = 'interface mask with the cells next to the connections'
    if a.inlet > 0:
        mask4 = [np.zeros(s, dtype=bool) for s in shapes]
        mask4[0][:, :, :a.inlet] = True
        mask1 = [m.copy() for m in mask4]
        what = 'mask of the %d inlet cells' % a.inlet
    print('  levels: %s; cells in the 4-cell %s: %d' % (
        np.bincount(np.concatenate([x.ravel() for x in lev])).tolist(), what, int(sum(m.sum() for m in mask4))))
    ok = True
    out = {'runs': a.runs, 'ref': a.ref, 'vars': {}}
    for v in a.vars.split(','):
        e1, einf, e4, e1c = [], [], [], []
        for r in runs:
            l1, li = mc.errors(r, ref, v)
            e1.append(l1)
            einf.append(li)
            e4.append(mc.errors(r, ref, v, mask4)[1])
            e1c.append(mc.errors(r, ref, v, mask1)[1])
        eint = [max(x, 0.05 * y) for x, y in zip(e4, einf)]
        eint1 = [max(x, 0.05 * y) for x, y in zip(e1c, einf)]
        o1, oi, oint, o1c = mc.orders(e1), mc.orders(einf), mc.orders(eint), mc.orders(eint1 if a.min_int1 else e1c)
        good = (a.l1_informative or o1[-1] >= a.min_l1) and oint[-1] >= a.min_int
        if a.min_int1 is not None:
            good = good and o1c[-1] >= a.min_int1
        ok = ok and good
        out['vars'][v] = {'L1': e1, 'Linf': einf, 'Lint4': eint, 'Lint1': eint1 if a.min_int1 else e1c,
                          'order_L1': o1[-1], 'order_int': oint[-1], 'ok': bool(good)}
        print('  %-7s L1      %s | orders %s' % (v, fmt(e1, '%.3e'), fmt(o1, '%.3f')))
        print('  %-7s Linf    %s | orders %s' % ('', fmt(einf, '%.3e'), fmt(oi, '%.3f')))
        print('  %-7s %s %s | orders %s  (max(mask, 5%% Linf): %s)' % (
            '', 'Linlet ' if a.inlet > 0 else 'Lint4  ', fmt(eint, '%.3e'), fmt(oint, '%.3f'),
            'mask' if all(x >= 0.05 * y for x, y in zip(e4, einf)) else 'floor used'))
        if a.min_int1 is not None:
            print('  %-7s Lint1   %s | orders %s  (max(mask, 5%% Linf): %s; last pair >= %.2f)' % (
                '', fmt(eint1, '%.3e'), fmt(o1c, '%.3f'),
                'mask' if all(x >= 0.05 * y for x, y in zip(e1c, einf)) else 'floor used', a.min_int1))
        else:
            print('  %-7s %s %s | orders %s  (information)' % (
                '', 'Linlet0' if a.inlet > 0 else 'Lint1  ', fmt(e1c, '%.3e'), fmt(o1c, '%.3f')))
        print('  %-7s last pair: L1 %.3f (%s), %s %.3f (>= %.2f)%s -> %s' % (
            '', o1[-1], 'information' if a.l1_informative else '>= %.2f' % a.min_l1,
            'inlet' if a.inlet > 0 else 'interface', oint[-1], a.min_int,
            '' if a.min_int1 is None else ', 1-cell interface %.3f (>= %.2f)' % (o1c[-1], a.min_int1),
            'ok' if good else 'FAIL'))
    if a.json:
        out['ok'] = bool(ok)
        with open(a.json, 'w') as f:
            json.dump(out, f, indent=1)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
