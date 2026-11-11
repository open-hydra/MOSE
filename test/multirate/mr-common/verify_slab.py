#!/usr/bin/env python3
"""Thin slabs (one cell between two slip walls: z in the 2D cases, y and z in the 1D cases): the wall mode across a slab,
dw/dt = -2 a w / dz, is unstable with Heun when a step exceeds dz/a. The multirate tests keep every step below it (the generator makes the cartesian slabs 10 mm thick):
  --h-max   the step of the coarsest populated level of the first macro step, H/2^lmin (log lines 'Multirate: H = ...'
            and 'cells per level' of solve.log, the first of each; lmin = coarsest populated level), must be <=
            t_min/a_max, with t_min the smallest slab thickness of the mesh and a_max = max sqrt(gamma p/rho) over the
            initial and final fields; with a fixed macro step and a fixed map it is the largest step of the run;
  --w-zero  information on cartesian meshes: w (and v in 1D) stays exactly 0.
Usage: verify_slab.py CASE_DIR RUN... (--h-max | --w-zero)"""
import argparse
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def thickness(blocks):
    """Smallest thickness of the slab directions: z for every block, y too for the blocks one cell wide in y."""
    t = []
    for b in blocks:
        P = b['points']
        t.append(float(np.abs(np.diff(P[..., 2], axis=0)).min()))
        if P.shape[1] == 2:
            t.append(float(np.abs(np.diff(P[..., 1], axis=1)).min()))
    return min(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('runs', nargs='+')
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--w-zero', action='store_true')
    g.add_argument('--h-max', action='store_true')
    a = ap.parse_args()
    ok = True
    ic = mc.read_ic_tec(os.path.join(a.case, 'INPUT', 'ic.tec'))
    tmin = thickness(ic)
    for r in a.runs:
        b = mc.read_run(r)
        name = os.path.basename(r.rstrip('/'))
        if a.w_zero:
            wmax = max(float(np.abs(x['cell']['w']).max()) for x in b)
            print('  %s: max |w| = %.3e' % (name, wmax))
            ok = ok and wmax == 0.0
            continue
        log = open(os.path.join(r, 'solve.log'), errors='replace').read()
        H = float(re.search(r'Multirate: H =\s*([0-9.Ee+-]+)', log).group(1))
        cnt = [int(x) for x in re.search(r'Multirate: cells per level\s+([0-9 ]+)', log).group(1).split()]
        lmin = next(l for l, c in enumerate(cnt) if c > 0)
        hmax = H / 2.0 ** lmin
        amax = max(float(np.sqrt(mc.GAMMA_AIR * x['cell']['p'] / x['cell']['rho(1)']).max()) for x in ic + b)
        good = hmax <= tmin / amax
        print('  %s: first macro step H/2^%d = %.4e s <= slab/a_max = %.4e s (slab %.3e m, a_max %.2f m/s): %s' % (
            name, lmin, hmax, tmin / amax, tmin, amax, 'ok' if good else 'ABOVE'))
        ok = ok and good
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
