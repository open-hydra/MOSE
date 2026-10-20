#!/usr/bin/env python3
"""Conservation check of the multirate tests: relative drift of the conserved totals between the
initial condition (INPUT/ic.tec of the case) and the final field of a run.
Usage: verify_conservation.py CASE_DIR RUN_DIR --quantities mass,momx,energy --max 1e-12      (criterion)
       verify_conservation.py CASE_DIR RUN_DIR --quantities mass,momx,energy --min 1e-9 --any  (positive control)
       verify_conservation.py CASE_DIR RUN_DIR --quantities energy --report                    (informative)
Quantities: mass, momx, momy, momz, energy (of a perfect gas) and mass1, mass2, ... (the mass of each species).
"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402

IDX = {'mass': 0, 'momx': 1, 'momy': 2, 'momz': 3, 'energy': 4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('run')
    ap.add_argument('--quantities', default='mass,momx,energy')
    ap.add_argument('--max', type=float, default=None)
    ap.add_argument('--min', type=float, default=None)
    ap.add_argument('--any', action='store_true')
    ap.add_argument('--report', action='store_true')
    a = ap.parse_args()
    b0, b1 = mc.read_ic_tec(os.path.join(a.case, 'INPUT', 'ic.tec')), mc.read_run(a.run)
    t0, t1 = mc.totals(b0), mc.totals(b1)
    res = []
    for q in a.quantities.split(','):
        if q in IDX:
            v0, v1 = t0[IDX[q]], t1[IDX[q]]
        else:
            n = int(q[len('mass'):])
            v0, v1 = mc.species_mass(b0, n), mc.species_mass(b1, n)
        res.append((q, abs(v1 - v0) / max(abs(v0), 1e-300)))
    print('  %s: ' % os.path.basename(a.run.rstrip('/')) + ', '.join('%s drift %.3e' % r for r in res))
    if a.report:
        return 0
    if a.max is not None:
        return 0 if all(d <= a.max for _, d in res) else 1
    if a.min is not None:
        test = any if a.any else all
        return 0 if test(d >= a.min for _, d in res) else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
