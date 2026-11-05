#!/usr/bin/env python3
"""Errors and orders of the simple wave with moving level boundaries (test mr-moving-levels) against its exact solution.

The exact solution (derivation): for the 1D Euler equations of a perfect gas with constant gamma and constant entropy (the
initial data have p = rho c^2 / gamma with rho = c^(2 / (gamma - 1)), so p / rho^gamma = 1 / gamma everywhere), the Riemann
invariants J+- = u +- 2 c / (gamma - 1) are constant along the characteristics dx/dt = u +- c. The initial data have
u = 2 (c - 1) / (gamma - 1), so J- = -2 / (gamma - 1) is the same on every C- characteristic, hence everywhere while the solution
stays smooth; on each C+ characteristic J+ is constant too, so u and c are constant along it and it is a straight line,
x = x0 + lambda(x0) t with lambda = u + c of the initial data at x0 (a simple wave). It stays smooth while the C+ lines do not
cross, 1 + t lambda'(x0) > 0 for every x0, that is t < 1 / max(-lambda') = 1 / (pi delta (gamma + 1) / (gamma - 1)) = 1.06.

For every run RUNS/<label>-n<n> of the case CASES/mr1d-simple-<n> (gen_mr_cases.py):
- exact solution before the breaking (t < 1.06): a simple wave, c and u constant along the characteristics
  x = x0 + lambda(x0) t, lambda = u + c = 1 + delta phi(x0) (gamma + 1) / (gamma - 1), with the gamma of the gas that the solver
  computes (gen_mr_cases.gas); x0(x) by Newton (F' = 1 + t lambda' > 0 before the breaking); cell averages of the conserved
  variables with Gauss-Legendre at 8 points, then the primitive variables of those averages (as MOSE writes them), at the final
  time of the run (TIME of its field);
- L1 (mean over the cells: the cells are equal) and Linf errors of rho, u and p, and the orders between successive n;
- near the faces (information): Linf of rho, u and p over the cells within 4 cells, on each side, of the two faces between
  levels, and its orders. The faces are placed from the exact solution: the cells of level 1 are the N1 cells whose foot x0 has
  the largest phi (the local step falls as phi grows, since u and c grow with it), N1 taken from the map of the run: the last
  map printed, at the final time, when the map is rebuilt (levels = recompute), the first, at t = 0, when it is kept (levels =
  fixed, read from the input.ini of the run); a single-rate run takes the faces of the run --faces LABEL3 (default mr) at the
  same n;
- validity: exit 0 is checked by run.sh; the final time equal to T = 0.5 within 1e-12 (relative); for the multirate runs (mr,
  mrf) at least two lines "Multirate: cells per level" in the log (MOSE prints the map at the first macro step and at the
  multiples of shell-diter, which run.sh sets to the number of macro steps: the first and the last map) and cells at both levels
  in every line: a run whose cells are all at one level does not test the multirate;
- conservation: relative drift of the totals of mass, momentum and energy between INPUT/ic.tec and the final field.
- with --ref LABEL2: L1 and Linf of the difference between the run and that of LABEL2 at the same n (same final time), its Linf
  on the cells near the faces, their orders, and the distance in cells of the largest difference from the nearest face
  (reported).
The exit code is 1 when a run is not valid; with --decide also when, on the last pair, the L1 order of rho, u or p is below 1.9
or a mass drift is above 1e-12 (the decided thresholds); the Linf orders, global and near the faces, the drifts of momentum and
energy and the differences are reported. --json FILE writes the errors for plot_moving_levels.py.
Usage: verify_moving_levels.py CASES RUNS LABEL n1 n2 ... [--ref LABEL2] [--faces LABEL3] [--multirate] [--decide] [--json FILE]"""
import argparse
import json
import math
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np            # noqa: E402
import mrcheck as mc          # noqa: E402
import gen_mr_cases as gm     # noqa: E402

T_END = 0.5
L1_MIN = 1.9
DRIFT_MAX = 1e-12
NEAR = 4                       # cells on each side of a face between levels in the Linf near the faces


def gas_gamma(case_dir):
    return gm.gas(os.path.join(case_dir, 'INPUT'))[1]


def lam_coef(gam):
    """lambda = 1 + coef phi(x0): u + c with u = 2 (c - 1) / (gamma - 1), c = 1 + delta phi"""
    return gm.SIMPLE_DELTA * (gam + 1.0) / (gam - 1.0)


def x0_of_factory(gam):
    k = lam_coef(gam)

    def x0_of(x, t):
        y = np.array(x, dtype=float) - t * (1.0 + 0.5 * k)      # first guess: the mean speed
        # F' = 1 + k phi'(y) t >= 1 - pi k t: 0.53 with k = 0.3 and t <= 0.5, at most 4 iterations on x in [-1, 2]; F' tends
        # to zero at the breaking time, t = 1.06, where the iteration diverges
        for _ in range(60):
            phi = 0.5 * (1.0 - np.cos(2.0 * np.pi * y))
            dphi = np.pi * np.sin(2.0 * np.pi * y)
            F = y + (1.0 + k * phi) * t - x
            y = y - F / (1.0 + k * dphi * t)
            if np.max(np.abs(F)) < 1e-15:
                break
        if np.max(np.abs(F)) >= 1e-13:
            raise SystemExit('x0: Newton did not converge (max |F| %.1e)' % np.max(np.abs(F)))
        return y
    return x0_of


def exact_prim(xn, gam, t):
    m, mu, E = gm.simple_averages(xn, gam, t, x0_of_factory(gam))
    u = mu / m
    return m, u, (gam - 1.0) * (E - 0.5 * m * u * u)


def map_rows(run):
    """every line "Multirate: cells per level" of the log, in order (the first map, then the multiples of shell-diter)"""
    rows = re.findall(r'Multirate: cells per level((?:\s+\d+)+)', open(os.path.join(run, 'solve.log'), errors='replace').read())
    return [[int(v) for v in r.split()] for r in rows]


def map_kept(run):
    """levels = fixed in the [MOSE-Multirate] section of the input.ini of the run"""
    sec = None
    for ln in open(os.path.join(run, 'input.ini'), errors='replace'):
        s = ln.split('#')[0].strip()
        if s.startswith('['):
            sec = s
        elif sec == '[MOSE-Multirate]' and '=' in s and s.split('=')[0].strip() == 'levels':
            return s.split('=', 1)[1].strip() == 'fixed'
    return False


def faces_of(run, gam, t_end, n):
    """faces between levels placed from the exact solution (cells p | p + 1, periodic) and the mask of the cells within NEAR
    cells of them: level 1 = the N1 cells whose foot x0 has the largest phi, N1 and the time from the map of the run"""
    rows = map_rows(run)
    if not rows or len(rows[0]) < 2:
        return None, None
    kept = map_kept(run)
    n1, t = (rows[0][1], 0.0) if kept else (rows[-1][1], t_end)
    xc = (np.arange(n) + 0.5) / n
    x0 = x0_of_factory(gam)(xc, t)
    phi = 0.5 * (1.0 - np.cos(2.0 * np.pi * x0))
    lev = np.zeros(n, dtype=int)
    lev[np.argsort(-phi, kind='stable')[:n1]] = 1
    faces = np.nonzero(lev != np.roll(lev, -1))[0]          # face between the cells p and p + 1 (mod n)
    mask = np.zeros(n, dtype=bool)
    for p in faces:
        for s in range(NEAR):
            mask[(p - s) % n] = True
            mask[(p + 1 + s) % n] = True
    return faces, mask


def face_distance(i, faces, n):
    """cells between the cell i and the nearest face (0: the cell next to a face), periodic"""
    d = [min((i - int(p) - 1) % n, (int(p) - i) % n) for p in faces]
    return int(min(d)) if d else None


def where(x, faces, n):
    """centre of the cell of the largest |x| and its distance in cells from the nearest face"""
    i = int(np.argmax(x))
    return (i + 0.5) / n, face_distance(i, faces, n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cases')
    ap.add_argument('runs')
    ap.add_argument('label')
    ap.add_argument('n', type=int, nargs='+')
    ap.add_argument('--ref', default=None)
    ap.add_argument('--faces', default='mr')
    ap.add_argument('--multirate', action='store_true')
    ap.add_argument('--decide', action='store_true')
    ap.add_argument('--json', default=None)
    a = ap.parse_args()
    res, valid = [], True
    for n in a.n:
        case = os.path.join(a.cases, 'mr1d-simple-%d' % n)
        run = os.path.join(a.runs, '%s-n%d' % (a.label, n))
        gam = gas_gamma(case)
        B = mc.read_run(run)[0]
        t = B['time']
        xn = np.linspace(0.0, 1.0, n + 1)
        rho, u, p = (B['cell'][v].ravel() for v in ('rho(1)', 'u', 'p'))
        er, eu, ep = exact_prim(xn, gam, t)
        e1 = [float(np.mean(np.abs(q - e))) for q, e in ((rho, er), (u, eu), (p, ep))]
        ei = [float(np.max(np.abs(q - e))) for q, e in ((rho, er), (u, eu), (p, ep))]
        ic = mc.read_ic_tec(os.path.join(case, 'INPUT', 'ic.tec'))[0]['cell']

        def totals(r, uu, pp):
            return np.array([r.sum(), (r * uu).sum(), (pp / (gam - 1.0) + 0.5 * r * uu * uu).sum()])
        T0 = totals(ic['rho(1)'].ravel(), ic['u'].ravel(), ic['p'].ravel())
        drift = [float(x) for x in np.abs(totals(rho, u, p) - T0) / np.abs(T0)]
        rec = {'n': n, 't': t, 'e1': e1, 'ei': ei, 'drift': drift, 'gamma': gam}
        why = []
        if abs(t - T_END) > 1e-12 * T_END:
            why.append('final time %.17g, not T = 0.5' % t)
        if a.multirate:
            rows = map_rows(run)
            rec['levels'] = rows
            if len(rows) < 2:
                why.append('%d line(s) "cells per level": the first and the last map are needed' % len(rows))
            for row in rows:
                if len(row) < 2 or min(row) == 0:
                    why.append('a map with the cells at one level: %s' % row)
        faces, mask = faces_of(run if a.multirate else os.path.join(a.runs, '%s-n%d' % (a.faces, n)), gam, t, n)
        if faces is not None:
            rec['faces'] = [int(f) for f in faces]
            rec['eloc'] = [float(np.max(np.abs(q - e)[mask])) for q, e in ((rho, er), (u, eu), (p, ep))]
            rec['emax_at'] = [where(np.abs(q - e), faces, n) for q, e in ((rho, er), (u, eu), (p, ep))]
        rec['invalid'] = why
        valid &= not why
        if a.ref:
            Br = mc.read_run(os.path.join(a.runs, '%s-n%d' % (a.ref, n)))[0]
            d = [np.abs(q - Br['cell'][v].ravel()) for q, v in ((rho, 'rho(1)'), (u, 'u'), (p, 'p'))]
            rec['ref_dt'] = Br['time'] - t
            rec['d1'] = [float(np.mean(x)) for x in d]
            rec['di'] = [float(np.max(x)) for x in d]
            if faces is not None:
                rec['dloc'] = [float(np.max(x[mask])) for x in d]
                rec['dmax_at'] = [where(x, faces, n) for x in d]
        res.append(rec)
    print('mr-moving-levels %s: final time, L1 and Linf errors of rho, u, p, drift of mass, momentum, energy (gamma %.9f)'
          % (a.label, res[0]['gamma']))
    for i, r in enumerate(res):
        s = '  n %5d  t %.15f  L1 %.3e %.3e %.3e  Linf %.3e %.3e %.3e  drift %.1e %.1e %.1e' % (
            r['n'], r['t'], *r['e1'], *r['ei'], *r['drift'])
        if 'levels' in r:
            rows = r['levels']
            s += '  cells per level (first and last map) %s %s' % ((rows[0], rows[-1]) if rows else (None, None))
        print(s)
        for w in r['invalid']:
            print('    NOT VALID: ' + w)
        if 'eloc' in r:
            print('    near the faces (between the cells %s, %d cells on each side): Linf %.3e %.3e %.3e'
                  % (' and '.join('%d|%d' % (f + 1, (f + 1) % r['n'] + 1) for f in r['faces']), NEAR, *r['eloc']))
            print('    largest error at x = %s' % ', '.join('%.4f (%d cells from a face)' % w for w in r['emax_at']))
        if 'd1' in r:
            print('    difference from %s (final time %+.1e): L1 %.3e %.3e %.3e  Linf %.3e %.3e %.3e'
                  % (a.ref, r['ref_dt'], *r['d1'], *r['di']))
            if 'dloc' in r:
                print('    difference from %s near the faces: Linf %.3e %.3e %.3e; largest difference at x = %s'
                      % (a.ref, *r['dloc'], ', '.join('%.4f (%d cells from a face)' % w for w in r['dmax_at'])))
        if i:
            q = res[i - 1]
            r['o1'] = [math.log2(q['e1'][k] / r['e1'][k]) for k in range(3)]
            r['oi'] = [math.log2(q['ei'][k] / r['ei'][k]) for k in range(3)]
            s = '    order %d -> %d: L1 %.3f %.3f %.3f  Linf %.3f %.3f %.3f' % (q['n'], r['n'], *r['o1'], *r['oi'])
            if 'eloc' in r and 'eloc' in q:
                r['oloc'] = [math.log2(q['eloc'][k] / r['eloc'][k]) for k in range(3)]
                s += '  Linf near the faces %.3f %.3f %.3f' % tuple(r['oloc'])
            print(s)
            if 'd1' in r:
                od = [math.log2(q['d1'][k] / r['d1'][k]) for k in range(3)]
                odi = [math.log2(q['di'][k] / r['di'][k]) for k in range(3)]
                s = '    order of the difference from %s %d -> %d: L1 %.3f %.3f %.3f  Linf %.3f %.3f %.3f' % (
                    a.ref, q['n'], r['n'], *od, *odi)
                if 'dloc' in r and 'dloc' in q:
                    s += '  Linf near the faces %.3f %.3f %.3f' % tuple(
                        math.log2(q['dloc'][k] / r['dloc'][k]) for k in range(3))
                print(s)
    ok = valid
    if a.decide:
        o1 = res[-1].get('o1', [0.0, 0.0, 0.0])
        mass = max(r['drift'][0] for r in res)
        ok = valid and min(o1) >= L1_MIN and mass <= DRIFT_MAX
        print('mr-moving-levels %s: L1 order on the last pair %.3f %.3f %.3f (at least %.1f); runs valid: %s; largest mass drift %.1e '
              '(at most %.0e): %s' % (a.label, *o1, L1_MIN, 'yes' if valid else 'NO', mass, DRIFT_MAX, 'ok' if ok else 'FAIL'))
    if a.json:
        json.dump(res, open(a.json, 'w'), indent=1)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
