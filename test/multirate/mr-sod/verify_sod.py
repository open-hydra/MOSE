#!/usr/bin/env python3
"""Sod tube with levels: multirate run against the single-rate run and the exact
Riemann solution (gamma of the solver, mrcheck) at the final time of each run, on the window xmin <= x <= xmax that the
reflection of the rarefaction head on the left wall has not reached (default 0.15 m: the head reaches the wall at
0.73 ms and at 1 ms is still below x = 0.09 m).
Checks: L1(rho - rho_exact) of the multirate <= 1.05 x that of the single rate; shock position within one local cell
of the exact one; number of local extrema of rho (tolerance 1e-4 of the range) not above the single rate's.
Usage: verify_sod.py RUN_MR RUN_SR [--x0 0.25] [--left 1.2,0,1e5] [--right 0.15,0,1e4] [--xmin 0.15] [--xmax 1e9]
"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def line_data(run):
    b = mc.read_run(run)[0]
    xc = mc.cell_centres(b['points'])[0, 0, :, 0]
    dx = np.diff(b['points'][0, 0, :, 0])
    return b['time'], xc, dx, b['cell']['rho(1)'][0, 0, :], mc.volumes(b['points'])[0, 0, :]


def analyse(run, x0, WL, WR, xmin, xmax):
    t, xc, dx, rho, V = line_data(run)
    re, ue, pe = mc.riemann_exact((xc - x0) / t, WL, WR)
    m = (xc >= xmin) & (xc <= xmax)
    l1 = np.sum(np.abs(rho - re)[m] * V[m]) / np.sum(V[m])
    # exact shock: the rightmost jump of the exact density, located on a fine grid of the similarity variable
    # (resolution 1e-3 m/s on the shock speed, i.e. 1e-6 m at 1 ms)
    xi = np.linspace(0.0, 2000.0, 2000001)
    rf, _, _ = mc.riemann_exact(xi, WL, WR)
    jf = np.nonzero(np.abs(np.diff(rf)) > 1e-9 * rf.max())[0][-1]
    speed = 0.5 * (xi[jf] + xi[jf + 1])
    x_shock = x0 + speed * t
    rho_post, rho_pre = rf[jf], rf[jf + 1]
    mid = 0.5 * (rho_post + rho_pre)
    # numerical: rightmost crossing of the mid value
    k = np.nonzero((rho[:-1] - mid) * (rho[1:] - mid) <= 0.0)[0]
    kk = k[-1]
    xs_num = xc[kk] + (mid - rho[kk]) * (xc[kk + 1] - xc[kk]) / (rho[kk + 1] - rho[kk])
    h = dx[kk]
    return dict(t=t, l1=l1, x_shock=x_shock, xs_num=xs_num, h=h, ext=mc.local_extrema(rho[m]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mr')
    ap.add_argument('sr')
    ap.add_argument('--x0', type=float, default=0.25)
    ap.add_argument('--left', default='1.2,0,1e5')
    ap.add_argument('--right', default='0.15,0,1e4')
    ap.add_argument('--xmin', type=float, default=0.15)
    ap.add_argument('--xmax', type=float, default=1e9)
    a = ap.parse_args()
    WL = [float(x) for x in a.left.split(',')]
    WR = [float(x) for x in a.right.split(',')]
    r = analyse(a.mr, a.x0, WL, WR, a.xmin, a.xmax)
    s = analyse(a.sr, a.x0, WL, WR, a.xmin, a.xmax)
    for nm, d in (('multirate', r), ('single rate', s)):
        print('  %-11s t %.6e  L1(rho - exact) %.4e  shock %.5f (exact %.5f, local cell %.2e)  extrema %d' %
              (nm, d['t'], d['l1'], d['xs_num'], d['x_shock'], d['h'], d['ext']))
    ok = r['l1'] <= 1.05 * s['l1'] and abs(r['xs_num'] - r['x_shock']) <= r['h'] and r['ext'] <= s['ext']
    print('  L1 ratio %.4f (<= 1.05), shock offset %.2e (<= %.2e), extrema %d <= %d -> %s' %
          (r['l1'] / s['l1'], abs(r['xs_num'] - r['x_shock']), r['h'], r['ext'], s['ext'], 'ok' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
