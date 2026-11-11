#!/usr/bin/env python3
"""Planar shock through the faces between levels: 2D box with the banded mesh of
mr2d-box, Riemann problem at x0 = 0.2 m (1.2 kg/m3, 0, 1e5 Pa | 0.15 kg/m3, 0, 1e4 Pa), slip walls. Each run is
compared with the exact 1D Riemann solution (gamma of the solver, mrcheck) at its own final time, on the window
x >= xmin = 0.25 m (the rarefaction reflected by the left wall stays below x = 0.15 m at 1 ms; the shock, at about
0.71 m at 1 ms, has not reached the right wall at 0.9984 m).
Checks of each multirate run: rho > 0 and p > 0; L1(rho - rho_exact) on the window (volume weighted) <= 1.05 x that of
the single rate at dt_min; on every row of cells along x, the number of local extrema of rho on the window (tolerance
1e-4 of the range of the row) not above the single rate's on the same row. Information: L1(rho_MR - rho_SR).
Usage: verify_shock.py RUN_SR RUN_MR... [--xmin 0.25]"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402

X0, WL, WR = 0.2, (1.2, 0.0, 1.0e5), (0.15, 0.0, 1.0e4)


def analyse(run, xmin):
    b = mc.read_run(run)[0]
    t = b['time']
    C = mc.cell_centres(b['points'])
    xc = C[0, :, :, 0]                                  # (nj, ni)
    V = mc.volumes(b['points'])[0]
    rho, p = b['cell']['rho(1)'][0], b['cell']['p'][0]
    re, _, _ = mc.riemann_exact((xc - X0) / t, WL, WR)
    m = xc >= xmin
    l1 = float(np.sum(np.abs(rho - re)[m] * V[m]) / np.sum(V[m]))
    ext = [mc.local_extrema(rho[j][m[j]]) for j in range(rho.shape[0])]
    return dict(t=t, l1=l1, ext=ext, pos=bool((rho > 0).all() and (p > 0).all()), rho=rho, V=V, m=m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sr')
    ap.add_argument('mr', nargs='+')
    ap.add_argument('--xmin', type=float, default=0.25)
    a = ap.parse_args()
    s = analyse(a.sr, a.xmin)
    print('  %-14s t %.6e  L1(rho - exact) %.4e  extrema per row: max %d, total %d' % (
        os.path.basename(a.sr.rstrip('/')), s['t'], s['l1'], max(s['ext']), sum(s['ext'])))
    ok = True
    for run in a.mr:
        r = analyse(run, a.xmin)
        rows = sum(1 for x, y in zip(r['ext'], s['ext']) if x > y)
        d = float(np.sum(np.abs(r['rho'] - s['rho'])[r['m']] * r['V'][r['m']]) / np.sum(r['V'][r['m']]))
        good = r['pos'] and r['l1'] <= 1.05 * s['l1'] and rows == 0
        ok = ok and good
        print('  %-14s t %.6e  L1(rho - exact) %.4e (ratio %.4f <= 1.05)  rho, p > 0: %s  rows with more extrema than '
              'the single rate: %d (extrema per row: max %d, total %d)  information: L1(rho_MR - rho_SR) %.3e -> %s' % (
                  os.path.basename(run.rstrip('/')), r['t'], r['l1'], r['l1'] / s['l1'], 'yes' if r['pos'] else 'NO',
                  rows, max(r['ext']), sum(r['ext']), d, 'ok' if good else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
