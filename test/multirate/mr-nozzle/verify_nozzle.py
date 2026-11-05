#!/usr/bin/env python3
"""Checks of the start-up of a nozzle run by run.sh (mr1d-nozzle-su, one block, ten solutions at t_k = k t_fin/10): the
reference ref = single rate at dt0/4, sr and sr2 = single rate at dt0 and dt0/2, mrb = multirate with H = dt0, mra and
mra2 = multirate with the production macro step H and H/2 (safety 2 and 4), all on the multirate path with the level map
rebuilt at every macro step. rho = rho(1), V the cell volumes; the throat region is the band of 0.64 mm cells around the
throat (0.444 m < x < 0.508 m), where the cells take the finest levels and where mrb has its levels.
  e_X(t_k) = sum_throat V |rho_X - rho_ref| / sum_throat V rho_ref; delta = the same with sr and sr2 in place of X and
  ref, Delta = the same with mra and mra2 (normalised by ref).
  D1  the six runs wrote their ten solutions with positive density and pressure, the solutions of the same k at the same
      time in every run within 40 M eps (M steps of sr per t_fin/10), and no cell of the fixed-step runs beyond its limit
      (no run line of h/(dtcap/sigma), which MOSE writes only when a macro step started with a cell above 1 + 1e-12
      times its limit; otherwise the check is not decisive). Conservation is the business of
      the test mr-conservation: here the boundaries are open.
  D2  mrb in the throat: e_mrb <= 1.05 e_sr (or <= 1e-10 when e_sr < 1e-10) at every t_k at which mrb has cells on a
      level >= 1 (the last macro step before t_k); the levels add no error to the step of the single rate.
  D3  mra in the throat: (i) e_mra2 <= 1.05 max(e_sr, 1e-10) + delta + Delta at every t_k; (ii) the observed order
      p = log2(e_mra / e_mra2) >= 0.8 at every t_k at which the single rates converge in the throat, p_sr =
      log2(delta / e_sr2) >= 0.8, and e_mra >= 1e-10. If (ii) judges no t_k the check is not decisive.
  D4  the moving shock: the face of largest compression rho_{i+1}/rho_i in the divergent beyond the throat band
      (x > 0.508 m), judged at every t_k at which the jump of ref across it (largest rho of the 4 cells downstream over
      the smallest of the 4 upstream) is at least 1.5. With s = max(|i_sr - i_ref|, |i_sr2 - i_ref|) and Delta_i =
      |i_mra - i_mra2| (faces): |i_mra2 - i_ref| <= s + Delta_i + 1, hence |i_mra - i_ref| <= s + 2 Delta_i + 1. If no
      t_k is judged the check is not decisive.
Informative: the errors of every run in the throat and on the whole nozzle, the shock positions, the levels of mra on
the cells of the shock face.
With --json FILE the times, errors, shock positions and level counts go to FILE, for plot_nozzle.py.
Usage: verify_nozzle.py RUNS_DIR CASE_DIR M [--json FILE]"""
import json
import math
import os
import re
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import mrcheck as mc             # noqa: E402

RUNS = ('ref', 'sr', 'sr2', 'mrb', 'mra', 'mra2')
FIXED = ('ref', 'sr', 'sr2', 'mrb')
THROAT = (0.444, 0.508)          # the band of 0.64 mm cells around the throat (0.476 m)
TINY = 1e-10


def solution(run, k):
    vtm = os.path.join(run, 'OUTPUT', 'field%s.vtm' % ('' if k == 10 else str(k)))
    files = re.findall(r'file="([^"]*)"', open(vtm).read())
    return mc.read_vts(os.path.join(run, 'OUTPUT', 'vtk', os.path.basename(files[0])))


def largest_ratio(run):
    """the largest h/(dtcap/sigma) of the run line, or None: MOSE prints that line only when a macro step started with a
    cell above its limit (more than 1 + 1e-12)"""
    log = open(os.path.join(run, 'solve.log'), errors='replace').read()
    m = re.search(r'Multirate run: max h/\(dtcap/sigma\) at the start of the macro steps +([0-9.Ee+-]+)', log)
    return float(m.group(1)) if m else None


def level_counts(run):
    log = open(os.path.join(run, 'solve.log'), errors='replace').read()
    return [list(map(int, m.split())) for m in re.findall(r'Multirate: cells per level\s+([0-9 ]+)', log)]


def norm(a, b, ref, V, m):
    return float((V[m] * np.abs(a[m] - b[m])).sum() / (V[m] * ref[m]).sum())


def shock(x, rho):
    """face of largest compression rho_{i+1}/rho_i beyond the throat band, and the jump across it"""
    q = rho[1:] / rho[:-1]
    q[0.5 * (x[1:] + x[:-1]) <= THROAT[1]] = -np.inf
    i = int(np.argmax(q))                                   # face between cells i and i + 1
    up = rho[max(i - 3, 0):i + 1].min()
    down = rho[i + 1:i + 5].max()
    return i, float(down / up)


def log2r(a, b):
    return math.log2(a / b) if a > 0 and b > 0 else float('nan')


def main():
    d, case, M = sys.argv[1], sys.argv[2], int(sys.argv[3])
    jpath = sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else None
    ok, decisive = True, True
    F = {}
    for r in RUNS:
        try:
            F[r] = [solution(os.path.join(d, r), k) for k in range(1, 11)]
        except (OSError, IndexError, KeyError) as e:
            print('  D1 %s: solutions missing (%s): FAIL' % (r, e))
            return 1
        pos = all(float(B['cell']['rho(1)'].min()) > 0 and float(B['cell']['p'].min()) > 0 for B in F[r])
        print('  D1 %s: ten solutions, rho and p > 0: %s' % (r, 'yes' if pos else 'NO'))
        ok &= pos
    t = [F['ref'][k]['time'] for k in range(10)]
    tol = 40 * M * 2.0 ** -52
    dt = max(abs(F[r][k]['time'] - t[k]) / t[k] for r in RUNS for k in range(10))
    print('  D1 times of the solutions equal within %.2e (largest %.2e): %s' % (tol, dt, 'yes' if dt <= tol else 'NO'))
    ok &= dt <= tol
    for r in FIXED:
        g = largest_ratio(os.path.join(d, r))
        good = g is None
        print('  D1 %s: %s' % (r, 'no macro step started with a cell beyond its limit (no run line of h/(dtcap/sigma))'
                                  if good else 'largest h/(dtcap/sigma) %.15e: BEYOND, the check is not decisive' % g))
        decisive &= good
    x = mc.cell_centres(F['ref'][0]['points'])[0][0, :, 0]
    V = mc.volumes(F['ref'][0]['points'])[0, 0]
    m = (x > THROAT[0]) & (x < THROAT[1])
    allc = np.ones_like(m)
    rho = {r: [F[r][k]['cell']['rho(1)'][0, 0] for k in range(10)] for r in RUNS}
    e = {r: [norm(rho[r][k], rho['ref'][k], rho['ref'][k], V, m) for k in range(10)] for r in RUNS if r != 'ref'}
    eg = {r: [norm(rho[r][k], rho['ref'][k], rho['ref'][k], V, allc) for k in range(10)] for r in RUNS if r != 'ref'}
    delta = [norm(rho['sr'][k], rho['sr2'][k], rho['ref'][k], V, m) for k in range(10)]
    Delta = [norm(rho['mra'][k], rho['mra2'][k], rho['ref'][k], V, m) for k in range(10)]
    print('  throat errors (0.444 m < x < 0.508 m, %d cells) and whole nozzle, at t_1 ... t_10:' % int(m.sum()))
    for r in ('sr', 'sr2', 'mrb', 'mra', 'mra2'):
        print('    e_%-4s %s | whole %s' % (r, ' '.join('%.2e' % v for v in e[r]), ' '.join('%.2e' % v for v in eg[r])))
    print('    delta  %s' % ' '.join('%.2e' % v for v in delta))
    print('    Delta  %s' % ' '.join('%.2e' % v for v in Delta))
    lv = level_counts(os.path.join(d, 'mrb'))
    fine = [sum(lv[min((k + 1) * M, len(lv)) - 1][1:]) if lv else 0 for k in range(10)]
    judged2 = [k for k in range(10) if fine[k] > 0]
    d2 = all(e['mrb'][k] <= (1.05 * e['sr'][k] if e['sr'][k] >= TINY else TINY) for k in judged2)
    print('  D2 mrb cells on levels >= 1 at t_k: %s; judged at %d instants; e_mrb / e_sr: %s: %s'
          % (' '.join(str(v) for v in fine), len(judged2),
             ' '.join('%.4f' % (e['mrb'][k] / e['sr'][k]) if e['sr'][k] > 0 else 'n/a' for k in judged2),
             'ok' if d2 else 'NO'))
    if not judged2:
        print('  D2: mrb has no cell on a level >= 1 at any t_k: not decisive')
        decisive = False
    ok &= d2
    d3i = all(e['mra2'][k] <= 1.05 * max(e['sr'][k], TINY) + delta[k] + Delta[k] for k in range(10))
    print('  D3 (i) e_mra2 <= 1.05 max(e_sr, 1e-10) + delta + Delta at every t_k: %s' % ('ok' if d3i else 'NO'))
    psr = [log2r(delta[k], e['sr2'][k]) for k in range(10)]
    pmr = [log2r(e['mra'][k], e['mra2'][k]) for k in range(10)]
    judged3 = [k for k in range(10) if psr[k] >= 0.8 and e['mra'][k] >= TINY]
    d3ii = all(pmr[k] >= 0.8 for k in judged3)
    print('  D3 (ii) observed order of the single rates p_sr: %s' % ' '.join('%.2f' % v for v in psr))
    print('  D3 (ii) observed order of mra p: %s; judged at %d instants: %s'
          % (' '.join('%.2f' % v for v in pmr), len(judged3), 'ok' if d3ii else 'NO'))
    if not judged3:
        print('  D3 (ii) judges no t_k: not decisive')
        decisive = False
    ok &= d3i and d3ii
    pos = {}
    jumps = []
    for k in range(10):
        pos_k = {r: shock(x, rho[r][k]) for r in RUNS}
        pos[k] = {r: pos_k[r][0] for r in RUNS}
        jumps.append(pos_k['ref'][1])
    judged4 = [k for k in range(10) if jumps[k] >= 1.5]
    d4 = True
    lmra = level_counts(os.path.join(d, 'mra'))
    for k in range(10):
        i = pos[k]
        s = max(abs(i['sr'] - i['ref']), abs(i['sr2'] - i['ref']))
        di = abs(i['mra'] - i['mra2'])
        a2, a1 = abs(i['mra2'] - i['ref']), abs(i['mra'] - i['ref'])
        if k in judged4:
            good = a2 <= s + di + 1
            d4 &= good
        print('  D4 t_%-2d shock face of ref %d (x = %.4f m), jump %.2f%s; sr %+d, sr2 %+d, mra %+d, mra2 %+d faces: '
              '|mra2 - ref| %d <= s + Delta_i + 1 = %d%s, |mra - ref| %d <= %d'
              % (k + 1, i['ref'], 0.5 * (x[i['ref']] + x[i['ref'] + 1]), jumps[k], '' if k in judged4 else ' (not judged)',
                 i['sr'] - i['ref'], i['sr2'] - i['ref'], i['mra'] - i['ref'], i['mra2'] - i['ref'], a2, s + di + 1,
                 (': ok' if a2 <= s + di + 1 else ': NO') if k in judged4 else '', a1, s + 2 * di + 1))
    if not judged4:
        print('  D4 judges no t_k: not decisive')
        decisive = False
    ok &= d4
    if jpath:
        with open(jpath, 'w') as f:
            json.dump({'case': os.path.basename(os.path.normpath(case)), 'cells': int(x.size), 't': t,
                       'e': e, 'whole': eg, 'delta': delta, 'Delta': Delta, 'limit': 1.05,
                       'shock_x': {r: [float(0.5 * (x[pos[k][r]] + x[pos[k][r] + 1])) for k in range(10)] for r in RUNS},
                       'jump': jumps, 'judged_d4': judged4, 'mrb_fine': fine,
                       'throat': list(THROAT), 'ok': bool(ok), 'decisive': bool(decisive)}, f, indent=1)
    if not decisive:
        print('nozzle: NOT DECISIVE')
        return 1
    print('nozzle: %s' % ('ok' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
