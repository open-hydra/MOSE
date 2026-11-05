#!/usr/bin/env python3
"""Variables that must stay uniform (entropy wave: u and p): largest relative deviation, max over the cells of
|q - q_initial|/max|q_initial|, of each run from the initial condition of the case.
Usage: verify_uniform.py CASE_DIR RUN... [--vars u,p] [--max 1e-12]"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--vars', default='u,p')
    ap.add_argument('--max', type=float, default=1e-12)
    a = ap.parse_args()
    ic = mc.read_ic_tec(os.path.join(a.case, 'INPUT', 'ic.tec'))
    ok = True
    for r in a.runs:
        b = mc.read_run(r)
        parts = []
        for v in a.vars.split(','):
            d = max(float(np.abs(x['cell'][v] - y['cell'][v]).max() / np.abs(y['cell'][v]).max()) for x, y in zip(b, ic))
            parts.append('%s %.2e' % (v, d))
            ok = ok and d <= a.max
        print('  %s: uniform to %s' % (os.path.basename(r.rstrip('/')), ', '.join(parts)))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
