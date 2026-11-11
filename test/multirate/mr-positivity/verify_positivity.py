#!/usr/bin/env python3
"""Positivity of the turbulence variables with the multirate integrator and explicit turbulence sources. For one run:
  - the run ended with a final field;
  - the solver counted no cell update in which a stage made nu~, k or omega negative before the realizability floor:
    the run has no line 'Turbulence variables below zero after a stage (set to the floor): N', which the solver prints
    at the end of a RANS run only when N is above zero;
  - the local step was limited by the turbulence sources, so that the run tests the stiffness in the step: dt_min of the
    first macro step (log line 'Multirate: H = ... s = ... dt_min') below half the transport step of the initial state,
    computed from the initial field of the case with the rule of Compute_dt without the sources:
    1/dt = sum_d (|u_d| + c)/h_d / cfl + D_max sum_d 1/h_d^2 / vnn, D_max = max(4/3 (mu + mu_t)/rho,
    (k + mu_t cp/Pr_t)/(rho cv), SA: (mu + rho nu~)/(sigma rho)), on the cartesian channel of the case (h_x uniform, h_y per
    row; the k direction does not count in 2D); mu_t of SA from chi, of SST as rho k/omega (an upper bound of the SST
    eddy viscosity, so the transport step is not overestimated), as the test TimeIntegration/TurbulenceSources;
  - the final field holds the turbulence variables (rhoNuTilde; rhoK and rhoOmega), whose smallest values are printed.
Gas: test/common/Air-viscous (mu, k constant, cp 1004.5 J/kg/K, gas constant of FLINT).
Usage: verify_positivity.py CASE_DIR RUN_DIR"""
import math
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402

MU, KL, CP, PRT = 1.1858686e-5, 0.016544513989, 1004.5, 0.85
CV1, SIGMA = 7.1, 2.0 / 3.0
RGAS = mc.RU / 28.970418


def ini_value(path, section, key):
    sec = None
    for line in open(path):
        s = line.strip()
        if s.startswith('['):
            sec = s.strip('[]')
        elif sec == section and '=' in s and s.split('=')[0].strip() == key:
            return s.split('=', 1)[1].split(';')[0].strip()
    return None


def transport_step(case, cfl, vnn, model):
    b = mc.read_ic_tec(os.path.join(case, 'INPUT', 'ic.tec'))[0]
    P = b['points']
    hx = np.diff(P[0, 0, :, 0])            # cell widths along i (first row)
    hy = np.diff(P[0, :, 0, 1])            # cell heights along j (first column)
    c = b['cell']
    rho, u, v, p = c['rho(1)'][0], c['u'][0], c['v'][0], c['p'][0]
    cv = CP - RGAS
    snd = np.sqrt(CP / cv * p / rho)
    if model == 'SA':
        rnt = c['mi_tilde'][0]
        chi = rnt / MU
        mut = rnt * chi ** 3 / (chi ** 3 + CV1 ** 3)
        extra = (MU + rnt) / (SIGMA * rho)
    else:
        mut = c['kappa'][0] / (c['omega'][0] / rho)
        extra = np.zeros_like(rho)
    dmax = np.maximum(np.maximum(4.0 / 3.0 * (MU + mut) / rho, (KL + mut * CP / PRT) / (rho * cv)), extra)
    HX = hx[None, :]
    HY = hy[:, None]
    rate = ((np.abs(u) + snd) / HX + (np.abs(v) + snd) / HY) / cfl + dmax * (1.0 / HX ** 2 + 1.0 / HY ** 2) / vnn
    return float((1.0 / rate).min())


def main():
    case, run = sys.argv[1], sys.argv[2]
    name = os.path.basename(run.rstrip('/'))
    log = open(os.path.join(run, 'solve.log'), errors='replace').read()
    ok = True
    if not os.path.exists(os.path.join(run, 'OUTPUT', 'field.vtm')) or 'NaN or p<0' in log:
        print('  %s: [FAIL] the run did not end with a final field' % name)
        return 1
    counts = [int(x) for x in re.findall(r'[Tt]urbulence variables below zero after a stage[^:\n]*:\s*(\d+)', log)]
    if not counts:
        print('  %s: negative nu~, k or omega after a stage: none (no count in the solver output)' % name)
    else:
        print('  %s: [FAIL] negative nu~, k or omega after a stage: %d' % (name, max(counts)))
        ok = False
    m = re.search(r'Multirate: H =\s*([0-9.Ee+-]+) s =\s*([0-9.Ee+-]+) dt_min', log)
    if m is None:
        print('  %s: [FAIL] no macro step in the solver output' % name)
        return 1
    dtmin = float(m.group(1)) / float(m.group(2))
    ini = os.path.join(run, 'input.ini')
    model = ini_value(ini, 'MOSE-Physics', 'turbulence')
    tr = transport_step(case, float(ini_value(ini, 'MOSE-Numerics', 'cfl')), float(ini_value(ini, 'MOSE-Numerics', 'vnn')),
                        model)
    ratio = dtmin / tr
    print('  %s: dt_min of the first macro step %.4e s, transport step of the initial state %.4e s, ratio %.3f'
          % (name, dtmin, tr, ratio))
    if not ratio < 0.5:
        print('  %s: [FAIL] the step was not limited by the turbulence sources: the case does not test them' % name)
        ok = False
    fin = mc.read_run(run)
    turb = ['rhoNuTilde'] if model == 'SA' else ['rhoK', 'rhoOmega']
    for t in turb:
        keys = [k for k in fin[0]['cell'] if k.lower() == t.lower()]
        if keys:
            print('  %s: min %s in the final field %.6e' % (name, keys[0], min(float(x['cell'][keys[0]].min()) for x in fin)))
        else:
            print('  %s: [FAIL] no %s in the final field' % (name, t))
            ok = False
    print('  %s: %s' % (name, 'ok' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
