#!/usr/bin/env python3
"""Mirror symmetry (metamorphic test): the mesh of mr2d-box is symmetric about its centre
lines and its diagonal and the pulse of mr2d-box-sym sits at the centre, so the exact discrete solution is symmetric:
  mirror in x: rho, p, v even, u odd;  mirror in y: rho, p, u even, v odd;  diagonal (i <-> j): rho, p even, u <-> v.
For every run (blocks glued back by their cell centres): largest |q(i,j) - q(mirror)| / dq over the cells, the three
mirrors and the variables rho, u, v, p, with dq the excursion of q in the final field (max - min; for u and v the
larger of the two). The runs come in pairs, multirate and single rate on the same case. Criterion: every run <= 1e-10,
and the multirate <= 10 x the single rate of its pair. (Round-off sets the floor: the single-rate step of MOSE alone
reaches 3-4e-12 of the excursion in 100 steps, because the background state p0 = 1e5 Pa is about 400 times the
excursion of the pulse; an asymmetric defect of the scheme shows at the level of its truncation error, 1e-5 and above.)
Usage: verify_sym.py MR_RUN SR_RUN [MR_RUN SR_RUN ...] [--max 1e-10] [--ratio 10]"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def glue(run):
    """Final field on the whole (nj, ni) grid: dict of 2D arrays rho, u, v, p (cells ordered by their centres)."""
    blocks = mc.read_run(run)
    xs, ys = set(), set()
    for b in blocks:
        C = mc.cell_key(b['points'])
        xs.update(C[0, 0, :, 0].tolist()); ys.update(C[0, :, 0, 1].tolist())
    xs, ys = sorted(xs), sorted(ys)
    ix, iy = {x: k for k, x in enumerate(xs)}, {y: k for k, y in enumerate(ys)}
    out = {v: np.full((len(ys), len(xs)), np.nan) for v in ('rho(1)', 'u', 'v', 'p')}
    for b in blocks:
        C = mc.cell_key(b['points'])
        for j in range(C.shape[1]):
            for i in range(C.shape[2]):
                for v in out:
                    out[v][iy[C[0, j, i, 1]], ix[C[0, j, i, 0]]] = b['cell'][v][0, j, i]
    assert not any(np.isnan(a).any() for a in out.values()), 'cells missing in %s' % run
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--max', type=float, default=1e-10)
    ap.add_argument('--ratio', type=float, default=10.0)
    a = ap.parse_args()
    if len(a.runs) % 2:
        raise SystemExit('give the runs in pairs: multirate, single rate')
    ok = True
    worst_of = {}
    for run in a.runs:
        q = glue(run)
        rho, u, v, p = q['rho(1)'], q['u'], q['v'], q['p']
        dq = {'rho': np.ptp(rho), 'p': np.ptp(p), 'vel': max(np.ptp(u), np.ptp(v))}
        res = {
            'x': max(np.abs(rho - rho[:, ::-1]).max() / dq['rho'], np.abs(p - p[:, ::-1]).max() / dq['p'],
                     np.abs(u + u[:, ::-1]).max() / dq['vel'], np.abs(v - v[:, ::-1]).max() / dq['vel']),
            'y': max(np.abs(rho - rho[::-1, :]).max() / dq['rho'], np.abs(p - p[::-1, :]).max() / dq['p'],
                     np.abs(u - u[::-1, :]).max() / dq['vel'], np.abs(v + v[::-1, :]).max() / dq['vel']),
            'diagonal': max(np.abs(rho - rho.T).max() / dq['rho'], np.abs(p - p.T).max() / dq['p'],
                            np.abs(u - v.T).max() / dq['vel'], np.abs(v - u.T).max() / dq['vel'])}
        worst = max(res.values())
        worst_of[run] = worst
        ok = ok and worst <= a.max
        print('  %-12s %dx%d cells; excursions rho %.3e, p %.3e, velocity %.3e; largest asymmetry: x %.2e, y %.2e, '
              'diagonal %.2e (<= %.0e) -> %s' % (os.path.basename(run.rstrip('/')), rho.shape[1], rho.shape[0], dq['rho'],
                                                dq['p'], dq['vel'], res['x'], res['y'], res['diagonal'], a.max,
                                                'ok' if worst <= a.max else 'FAIL'))
    for mr, sr in zip(a.runs[0::2], a.runs[1::2]):
        good = worst_of[mr] <= a.ratio * worst_of[sr]
        ok = ok and good
        print('  %s against %s: asymmetry ratio %.2f (<= %.0f) -> %s' % (
            os.path.basename(mr.rstrip('/')), os.path.basename(sr.rstrip('/')), worst_of[mr] / worst_of[sr], a.ratio,
            'ok' if good else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
