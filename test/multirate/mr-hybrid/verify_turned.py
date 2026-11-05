#!/usr/bin/env python3
"""The run of a case with one block turned by 180 degrees about z (its i and j reversed, gen_mr_cases.py turn_block_text)
against the run of the same case without the turn: the nodes of the turned block, read with i and j reversed, equal to the
bit, and every cell value of rho, u, v, w, p, the turned block read with i and j reversed, within TOL times the largest
absolute value over the field of the variable (rho, p) or of the three velocity components (u, v, w: a component that the
case keeps at zero, as w on a 2D mesh, has only round-off as its own scale). The two runs differ only in the order of the floating-point operations of the
turned block (the faces of a cell summed in another order, the metrics of a cell from its nodes in another order), so the
values agree to round-off; a connection whose orientation the solver read wrong would take a tangential neighbour from the
other side of the cell and change the viscous fluxes at the connection by far more.
Usage: verify_turned.py RUN RUN_TURNED BLOCK [--tol 1e-10]   (exit 0 when every check passes)"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
from mrcheck import read_run  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('turned')
    ap.add_argument('block', type=int)
    ap.add_argument('--tol', type=float, default=1e-10)
    a = ap.parse_args()
    A, B = read_run(a.run), read_run(a.turned)
    ok = len(A) == len(B)
    if not ok:
        print('  blocks: %d and %d' % (len(A), len(B)))
        return 1
    nodes_ok = True
    for n, (x, y) in enumerate(zip(A, B), 1):
        px = x['points'][:, ::-1, ::-1, :] if n == a.block else x['points']
        nodes_ok &= px.shape == y['points'].shape and bool(np.array_equal(px, y['points']))
    print('  nodes of the turned block, read with i and j reversed: %s' % ('equal to the bit' if nodes_ok else 'DIFFERENT'))
    ok &= nodes_ok
    vel = max(float(np.max(np.abs(x['cell'][c]))) for x in A for c in ('u', 'v', 'w'))
    for var in ('rho(1)', 'u', 'v', 'w', 'p'):
        scale = vel if var in ('u', 'v', 'w') else max(float(np.max(np.abs(x['cell'][var]))) for x in A)
        diff = 0.0
        for n, (x, y) in enumerate(zip(A, B), 1):
            cx = x['cell'][var][:, ::-1, ::-1] if n == a.block else x['cell'][var]
            diff = max(diff, float(np.max(np.abs(cx - y['cell'][var]))))
        good = diff <= a.tol * scale
        print('  %-6s largest difference %.3e, %.3e of the scale %.6e: %s'
              % (var, diff, diff / scale if scale > 0 else 0.0, scale, 'ok' if good else 'ABOVE %.0e' % a.tol))
        ok &= good
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
