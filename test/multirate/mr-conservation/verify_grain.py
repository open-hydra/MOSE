#!/usr/bin/env python3
"""Mass balance of an SRM grain in its inflow mode (502, mode 2) at the closed end of a duct (mr1d-502c of
gen_mr_cases.py): the grain on face 1 of block 1 has the burn-rate exponent n = 0, so its mass flux is mdot_ref at
every pressure, and the other faces are slip walls (300) and null faces (0), so the mass of the domain at the time t of
the final field is M(0) + mdot_ref A t, A the area of face 1 of block 1. mdot_ref and n are read from the line of the
grain in INPUT/bc.txt of the case.
Usage: verify_grain.py CASE_DIR RUN_DIR --max 1e-12      (relative difference to M(0) + mdot_ref A t: criterion)
       verify_grain.py CASE_DIR RUN_DIR --min 1e-9       (positive control)
"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def grain_line(case_dir):
    """mdot_ref and n of the first record of type 502 in INPUT/bc.txt (its property line follows the record line)"""
    lines = open(os.path.join(case_dir, 'INPUT', 'bc.txt')).read().split('\n')
    for k, ln in enumerate(lines):
        s = ln.split()
        if len(s) == 6 and s[5] == '502':
            v = [float(x) for x in lines[k + 1].replace(',', ' ').split()]
            return v[0], v[1]
    raise SystemExit('verify_grain: no record of type 502 in %s' % case_dir)


def face1_area(points):
    """area of the face i = 1 of a block, (nk+1, nj+1, ni+1, 3) nodes: the sum over its quadrilaterals of half the norm
    of the cross product of their diagonals (the area of a planar quadrilateral)"""
    P = points[:, :, 0, :]
    diag_a = P[1:, 1:] - P[:-1, :-1]
    diag_b = P[1:, :-1] - P[:-1, 1:]
    return float(np.sum(0.5 * np.linalg.norm(np.cross(diag_a, diag_b), axis=-1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('run')
    ap.add_argument('--max', type=float, default=None)
    ap.add_argument('--min', type=float, default=None)
    a = ap.parse_args()
    mdot, n = grain_line(a.case)
    if n != 0.0:
        print('verify_grain: n = %g in the line of the grain: the balance M(0) + mdot_ref A t needs n = 0' % n)
        return 1
    b0, b1 = mc.read_ic_tec(os.path.join(a.case, 'INPUT', 'ic.tec')), mc.read_run(a.run)
    if 'time' not in b1[0]:
        print('verify_grain: no TIME in the final field of %s' % a.run)
        return 1
    m0, m1, t = mc.totals(b0)[0], mc.totals(b1)[0], b1[0]['time']
    expect = m0 + mdot * face1_area(b1[0]['points']) * t
    d = abs(m1 - expect) / abs(expect)
    print('  %s: t = %.6e s, mass %.15e, M(0) + mdot_ref A t = %.15e (M(0) = %.6e): relative difference %.3e'
          % (os.path.basename(a.run.rstrip('/')), t, m1, expect, m0, d))
    if a.max is not None:
        return 0 if d <= a.max else 1
    if a.min is not None:
        return 0 if d >= a.min else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
