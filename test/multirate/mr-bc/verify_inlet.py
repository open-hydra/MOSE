#!/usr/bin/env python3
"""The data of the time-dependent inlets reach the cells (mr-bc, cases mr1d-402 and mr1d-410). The order of the test
compares runs that share the data, so it holds also when the data do not act: a table 402 whose p0 lies far below the
pressure of the duct makes a wall of the inlet, and data that stay at their first value make a steady inlet. Two runs,
one that ends near the peak of the data (t_p, about Tp/4) and one that ends at the end of the order runs (T):
  - the flow enters in the four inlet cells, u > 0, in both runs;
  - q is the total pressure of the first cell (402) or its pressure (410), q_d(t) the value of the table at the time of
    the run (linear in time, as MOSE reads it) and A half of the range of the table:
      |q(t_p) - q_d(t_p)| <= |q_d(t_p) - q_d(0)| / 2: a first cell that stays at the first value of the data is at the
      whole distance |q_d(t_p) - q_d(0)|, while one that follows the data lags by half a cell over u + a, which at a
      peak of the data changes their value only to second order; half of that distance separates the two;
      |q(T) - q_d(T)| <= A: with an error above the signal the run does not measure the data.
Usage: verify_inlet.py CASE_DIR RUN_PEAK RUN_END"""
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def table(case):
    """Times and values of the table of the inlet, and the name of the value: p0.dat of a 402, q2d.dat of a 410."""
    path = os.path.join(case, 'p0.dat')
    if os.path.exists(path):
        t, v = np.loadtxt(path, unpack=True)
        return t, v, 'total pressure'
    lines = open(os.path.join(case, 'q2d.dat')).read().split('\n')
    t, v = [], []
    for i, line in enumerate(lines):
        m = re.search(r'SOLUTIONTIME=(\S+)', line)
        if m:
            t.append(float(m.group(1)))
            v.append(float(lines[i + 8]))     # three lines of node coordinates, then rho(1), u, v, w, p
    return np.array(t), np.array(v), 'pressure'


def first_cells(run, total, g=mc.GAMMA_AIR):
    """Time of the field of a run, u of its four inlet cells, and the total pressure or the pressure of the first cell."""
    b = mc.read_run(run)[0]
    u = b['cell']['u'].ravel()
    p, rho = b['cell']['p'].ravel()[0], b['cell']['rho(1)'].ravel()[0]
    q = p * (1.0 + 0.5 * (g - 1.0) * u[0] ** 2 * rho / (g * p)) ** (g / (g - 1.0)) if total else p
    return b['time'], u[:4], q


def main():
    case, run_peak, run_end = sys.argv[1:4]
    tt, vv, name = table(case)
    amp = 0.5 * (vv.max() - vv.min())
    ok = True
    for lab, run in (('peak of the data', run_peak), ('end of the run', run_end)):
        time, u, q = first_cells(run, name == 'total pressure')
        qd = float(np.interp(time, tt, vv))
        lim = 0.5 * abs(qd - vv[0]) if run == run_peak else amp
        good = bool(np.all(u > 0.0))
        print('  %s, t = %.4e s: u of the inlet cells %s (> 0): %s' % (lab, time, ' '.join('%.3e' % x for x in u),
                                                                       'ok' if good else 'FAIL'))
        ok = ok and good
        good = abs(q - qd) <= lim
        print('    %s of the first cell %.6e, table %.6e, |difference| %.3e <= %.3e: %s' % (name, q, qd, abs(q - qd), lim,
                                                                                        'ok' if good else 'FAIL'))
        ok = ok and good
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
