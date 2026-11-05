#!/usr/bin/env python3
"""Checks of the viscous shock tube of Daru and Tenaud run by run.sh (one block, ten numbered solutions at
t_k = k t_fin/10): the reference ref = single rate at dt0/4, sr = single rate at dt0, mra = multirate with H = 4 dt0 and
mrb = multirate with H = dt0, all on the multirate path with the level map rebuilt at every macro step.
  U1  every run wrote its ten solutions, with positive density and pressure;
      the solutions of the same k are at the same time in every run, within 40 M eps (M steps of sr per t_fin/10: the
      times are sums of steps, so runs with different numbers of steps differ by rounding);
      no cell beyond its limit: no run prints the line of h/(dtcap/sigma), which MOSE writes only when a macro step
      started with a cell above 1 + 1e-12 times its limit;
  U2  e_X(t_k) = sum V |rho_X - rho_ref| / sum V rho_ref; mrb, whose core cells take the step of sr: e_mrb <= 1.05 e_sr
      at every t_k (mra, whose core cells take 4 dt0, reported);
  U3  on the two control rows (the first row of cells, next to the wall y = 0, and the last, on the symmetry plane
      y = 0.5) the face with the largest |d rho / ds| of mra and mrb within one local cell of that of ref at every t_k;
  U4  the local extrema of rho on the same rows (changes of sign of the differences above 1e-4 of the range) of mra
      and mrb at most those of sr at every t_k;
  U5  at ten probes, the cells nearest to (x; y) = (0.6, 0.7, 0.8, 0.9, 0.95) x (0.0025, 0.25), at the t_k where no
      jump of ref (neighbours with max/min - 1 > 0.1) lies within 3 cells of the probe: the largest |X - ref| of rho,
      u, v, p and T of mrb at most 1.05 times that of sr for every variable (mra reported against 0.2 % of |ref|,
      0.2 % of the speed of sound for u and v);
  U6  the box is closed (walls and symmetry): mass and total energy of every run drift by at most 1e-12 from the
      initial condition.
  Informative, full form only (500 x 250 cells, the resolution of the published solutions): the density along the wall
  at t_fin divided by rho_R/1.2, at the 20 abscissae of table 1 of Zhou, Xu and Liu (arXiv:1705.09062), next to the
  values of the table. The short form (50 x 25 cells, run in CI) and the medium form (100 x 50 cells) are judged against
  their own reference ref only.
With --json FILE the form and mesh, the output times, the errors of U2 and, for the full form, the wall densities at
the abscissae of the table go to FILE, for plot_daru.py; the checks do not change.
Usage: verify_daru.py RUNS_DIR CASE_DIR M [--json FILE]"""
import json
import math
import os
import re
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import mrcheck as mc             # noqa: E402
from gen_mr_cases import gas    # noqa: E402

RUNS = ('ref', 'sr', 'mra', 'mrb')
VARS = ('rho(1)', 'u', 'v', 'p', 'T')
FULL_MESH = (500, 250)            # cells along x and y of the full form (mr2d-daru-full)
FORMS = {(50, 25): 'short', (100, 50): 'medium', FULL_MESH: 'full'}   # mr2d-daru-coarse, mr2d-daru, mr2d-daru-full
ZHOU = [(0.3030, 39.8418), (0.4490, 37.0662), (0.5230, 52.6465), (0.5730, 42.4400), (0.5930, 40.5506),
        (0.6123, 47.3367), (0.6370, 39.3203), (0.6577, 36.9558), (0.6830, 49.6513), (0.7070, 77.9810),
        (0.7317, 108.3916), (0.7543, 92.5760), (0.7790, 64.5319), (0.7957, 59.4386), (0.8183, 95.3607),
        (0.8617, 117.6452), (0.9437, 96.4287), (0.9670, 98.2689), (0.9883, 81.8465), (0.9943, 82.7077)]


def solution(run, k):
    vtm = os.path.join(run, 'OUTPUT', 'field%s.vtm' % ('' if k == 10 else str(k)))
    files = re.findall(r'file="([^"]*)"', open(vtm).read())
    return mc.read_vts(os.path.join(run, 'OUTPUT', 'vtk', os.path.basename(files[0])))


def row(B, which):
    """(s, rho) along the first (wall) or the last (symmetry) row of cells; s = distance between the cell centres"""
    C = mc.cell_centres(B['points'])[0]
    jj = 0 if which == 'wall' else C.shape[0] - 1
    xy = C[jj, :, :2]
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))])
    return s, B['cell']['rho(1)'][0, jj, :]


def jump(s, q):
    """the face with the largest |dq/ds|, and the longer of its two cells"""
    g = np.abs(np.diff(q)) / np.diff(s)
    i = int(np.argmax(g))
    ln = np.empty_like(s)
    ln[1:-1] = 0.5 * (s[2:] - s[:-2]); ln[0] = s[1] - s[0]; ln[-1] = s[-1] - s[-2]
    return 0.5 * (s[i] + s[i + 1]), max(ln[i], ln[i + 1])


def largest_ratio(run):
    """the largest h/(dtcap/sigma) of the run line, or None: MOSE prints that line only when a macro step started with a
    cell above its limit (more than 1 + 1e-12)"""
    log = open(os.path.join(run, 'solve.log'), errors='replace').read()
    m = re.search(r'Multirate run: max h/\(dtcap/sigma\) at the start of the macro steps +([0-9.Ee+-]+)', log)
    return float(m.group(1)) if m else None


def main():
    d, case, M = sys.argv[1], sys.argv[2], int(sys.argv[3])
    jpath = sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else None
    R, gam = gas(os.path.join(case, 'INPUT'))
    ok = True
    F = {}
    for r in RUNS:
        try:
            F[r] = [solution(os.path.join(d, r), k) for k in range(1, 11)]
        except (OSError, IndexError, KeyError) as e:
            print('  U1 %s: solutions missing (%s): FAIL' % (r, e))
            return 1
        pos = all(float(B['cell']['rho(1)'].min()) > 0 and float(B['cell']['p'].min()) > 0 for B in F[r])
        print('  U1 %s: ten solutions, rho and p > 0: %s' % (r, 'yes' if pos else 'NO'))
        ok &= pos
    t = [F['ref'][k]['time'] for k in range(10)]
    tol = 40 * M * 2.0 ** -52
    dt = max(abs(F[r][k]['time'] - t[k]) / t[k] for r in RUNS for k in range(10))
    print('  U1 times of the solutions equal within %.2e (largest %.2e): %s' % (tol, dt, 'yes' if dt <= tol else 'NO'))
    ok &= dt <= tol
    for r in RUNS:
        x = largest_ratio(os.path.join(d, r))
        good = x is None
        print('  U1 %s: %s' % (r, 'no macro step started with a cell beyond its limit (no run line of h/(dtcap/sigma))'
                                  if good else 'largest h/(dtcap/sigma) %.15e: BEYOND' % x))
        ok &= good
    V = mc.volumes(F['ref'][0]['points'])
    e = {}
    for r in ('sr', 'mra', 'mrb'):
        e[r] = []
        for k in range(10):
            num = float((V * np.abs(F[r][k]['cell']['rho(1)'] - F['ref'][k]['cell']['rho(1)'])).sum())
            e[r].append(num / float((V * F['ref'][k]['cell']['rho(1)']).sum()))
        print('  U2 e_%s: %s' % (r, ' '.join('%.3e' % v for v in e[r])))
    worst = max(e['mrb'][k] / e['sr'][k] if e['sr'][k] > 0 else (0.0 if e['mrb'][k] == 0 else math.inf) for k in range(10))
    print('  U2 e_mrb <= 1.05 e_sr at every t_k: %s (largest ratio %.4f)' % ('yes' if worst <= 1.05 else 'NO', worst))
    ok &= worst <= 1.05
    for which in ('wall', 'symmetry'):
        for k in range(10):
            sR, qR = row(F['ref'][k], which)
            pR, lR = jump(sR, qR)
            nS = mc.local_extrema(row(F['sr'][k], which)[1], rel_tol=1e-4)
            msg = []
            for r in ('mra', 'mrb'):
                s, q = row(F[r][k], which)
                pX, _ = jump(s, q)
                n = mc.local_extrema(q, rel_tol=1e-4)
                c3, c4 = abs(pX - pR) <= lR, n <= nS
                ok &= c3 and c4
                msg.append('%s jump %+.2e (%s), extrema %d (%s)' % (r, pX - pR, 'ok' if c3 else 'NO', n, 'ok' if c4 else 'NO'))
            print('  U3/U4 %-8s t_%-2d: ref jump at s = %.4f (cell %.2e), sr extrema %d; %s'
                  % (which, k + 1, pR, lR, nS, '; '.join(msg)))
    C = mc.cell_centres(F['ref'][0]['points'])[0]
    probes = []
    for y in (0.0025, 0.25):
        for x in (0.6, 0.7, 0.8, 0.9, 0.95):
            dd = (C[..., 0] - x) ** 2 + (C[..., 1] - y) ** 2
            probes.append(np.unravel_index(int(np.argmin(dd)), dd.shape))
    eb = {r: {v: 0.0 for v in VARS} for r in ('sr', 'mra', 'mrb')}
    ea = True
    judged = left_out = 0
    for k in range(10):
        rho = F['ref'][k]['cell']['rho(1)'][0]
        for j, i in probes:
            box = rho[max(j - 3, 0):j + 4, max(i - 3, 0):i + 4]
            big = 0.0
            if box.shape[1] > 1:
                big = max(big, float((np.maximum(box[:, 1:], box[:, :-1]) / np.minimum(box[:, 1:], box[:, :-1]) - 1).max()))
            if box.shape[0] > 1:
                big = max(big, float((np.maximum(box[1:], box[:-1]) / np.minimum(box[1:], box[:-1]) - 1).max()))
            if big > 0.1:
                left_out += 1
                continue
            judged += 1
            ref = {v: float(F['ref'][k]['cell'][v][0, j, i]) for v in VARS}
            a = math.sqrt(gam * R * ref['T'])
            for r in eb:
                for v in VARS:
                    dv = abs(float(F[r][k]['cell'][v][0, j, i]) - ref[v])
                    eb[r][v] = max(eb[r][v], dv)
                    if r == 'mra' and dv > 0.002 * (a if v in ('u', 'v') else abs(ref[v])):
                        ea = False
    print('  U5 probe instants judged %d, left out (a jump of ref within 3 cells) %d' % (judged, left_out))
    for r in eb:
        print('  U5 %s: %s' % (r, ' '.join('%s %.3e' % (v, eb[r][v]) for v in VARS)))
    u5 = judged > 0 and all(eb['mrb'][v] <= 1.05 * eb['sr'][v] for v in VARS)
    print('  U5 mrb within 1.05 times sr for every variable: %s' % ('yes' if u5 else 'NO'))
    ok &= u5
    print('  U5 mra within 0.2 %% of ref (reported): %s' % ('yes' if ea else 'no'))
    ic = mc.read_ic_tec(os.path.join(case, 'INPUT', 'ic.tec'))[0]

    def totals(B):
        c = B['cell']
        m = float((V * c['rho(1)']).sum())
        en = float((V * (c['p'] / (gam - 1.0) + 0.5 * c['rho(1)'] * (c['u'] ** 2 + c['v'] ** 2 + c['w'] ** 2))).sum())
        return m, en
    m0, e0 = totals(ic)
    for r in RUNS:
        m1, e1 = totals(F[r][9])
        dm, de = abs(m1 - m0) / m0, abs(e1 - e0) / e0
        good = max(dm, de) <= 1e-12
        print('  U6 %s: mass drift %.2e, energy drift %.2e: %s' % (r, dm, de, 'ok' if good else 'NO'))
        ok &= good
    rho_ref = 1.0e5 / (R * 300.0) / 1.2
    nj, ni = F['ref'][9]['cell']['rho(1)'].shape[1:]
    full = (ni, nj) == FULL_MESH
    form = FORMS.get((ni, nj), '%d x %d' % (ni, nj))
    cols = {}
    if full:
        print('  informative (full form, %d x %d cells): wall density / (rho_R/1.2) at t_fin, against table 1 of Zhou, Xu'
              ' and Liu' % (ni, nj))
        for r in RUNS:
            s, q = row(F[r][9], 'wall')
            xc = mc.cell_centres(F[r][9]['points'])[0][0, :, 0]
            cols[r] = np.interp([x for x, _ in ZHOU], xc, q / rho_ref)
        for n, (x, z) in enumerate(ZHOU):
            print('    x %.4f  table %8.3f  %s' % (x, z, '  '.join('%s %8.3f' % (r, cols[r][n]) for r in RUNS)))
    else:
        print('  %s form (%d x %d cells): the comparison with table 1 of Zhou, Xu and Liu is made on the full form only'
              ' (%d x %d cells)' % (form, ni, nj, FULL_MESH[0], FULL_MESH[1]))
    if jpath:
        with open(jpath, 'w') as f:
            json.dump({'case': os.path.basename(os.path.normpath(case)), 'form': form,
                       'mesh': [int(ni), int(nj)], 't': t, 'e': e, 'limit': 1.05, 'worst': worst,
                       'rho_scale': rho_ref, 'zhou': ZHOU if full else [],
                       'wall': {r: [float(x) for x in cols[r]] for r in cols}, 'ok': bool(ok)}, f, indent=1)
    print('daru-tenaud: %s' % ('ok' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
