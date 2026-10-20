#!/usr/bin/env python3
"""Free-stream preservation: largest relative deviation of the conserved variables of a run from the
uniform initial state, max over the cells of |q - q0|/|q0| for rho, rho u, rho v, rho E, and |rho w|/(rho |u|).
Usage: verify_free.py CASE_DIR RUN_DIR [--max 1e-13] [--report]"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def cons(c, g):
    rho, u, v, w, p = c['rho(1)'], c['u'], c['v'], c['w'], c['p']
    return rho, rho * u, rho * v, p / (g - 1.0) + 0.5 * rho * (u * u + v * v + w * w), rho * w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('run')
    ap.add_argument('--max', type=float, default=1e-13)
    ap.add_argument('--report', action='store_true')
    a = ap.parse_args()
    g = mc.GAMMA_AIR
    # all blocks (cases of several blocks); the blocks of the output follow those of the initial condition
    b0 = mc.read_ic_tec(os.path.join(a.case, 'INPUT', 'ic.tec'))
    b1 = mc.read_run(a.run)
    q0 = [np.concatenate([x.ravel() for x in q]) for q in zip(*[cons(b['cell'], g) for b in b0])]
    q1 = [np.concatenate([x.ravel() for x in q]) for q in zip(*[cons(b['cell'], g) for b in b1])]
    dev = {}
    for name, x0, x1 in zip(('rho', 'rho u', 'rho v', 'rho E'), q0[:4], q1[:4]):
        dev[name] = float(np.abs(x1 - x0).max() / np.abs(x0).max())
    dev['rho w'] = float(np.abs(q1[4]).max() / np.hypot(q0[1], q0[2]).max())
    worst = max(dev.values())
    print('  %s: ' % os.path.basename(a.run.rstrip('/')) + ', '.join('%s %.2e' % kv for kv in dev.items()))
    if a.report:
        return 0
    return 0 if worst <= a.max else 1


if __name__ == '__main__':
    sys.exit(main())
