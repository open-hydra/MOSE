#!/usr/bin/env python3
"""Verify the Reynolds-stress decay case against its analytic solution.

The case (see build_case.py) is a uniform stream, far from any wall, carrying
isotropic turbulence through an axisymmetric duct sector. Without mean strain
the SSG/LRR-omega model reduces to

    U dk/dx     = - C_mu k omega
    U domega/dx = - beta omega^2

(the pressure-strain model vanishes with the anisotropy; far from walls F1 = 0
and beta is that of the SSG set), whose solution from the inlet values is

    omega(x) = omega_0 / (1 + beta omega_0 x/U)
    k(x)     = k_0 (1 + beta omega_0 x/U)^(-C_mu/beta)

Four properties of the solution are checked, none of them fitted:

  1. the decay of k and omega along x follows the closed form;
  2. the turbulence stays isotropic (R11 = R22 = R33, no shear stress);
  3. the flow stays uniform across the duct, with no radial velocity: on the
     side faces of the wedge the azimuthal Reynolds stress must balance the
     radial one, as the pressure does;
  4. energy is conserved: along x the total enthalpy of the mean flow rises by
     what the turbulence loses. The Reynolds stresses do work on the mean flow
     like a pressure, so the conserved sum is H + k + R11, and
     dH = -(5/3) dk for isotropic turbulence.

Exit status: 0 = pass, 1 = fail.
"""

import re
import sys

import numpy as np

C_MU, BETA = 0.09, 0.0828      # SSG (epsilon) set of the model

failures = []
checks = 0


def check(what, ok, detail=''):
    global checks
    checks += 1
    status = 'ok  ' if ok else 'FAIL'
    print(f'   [{status}] {what}{("  " + detail) if detail else ""}')
    if not ok:
        failures.append(what)


def below(what, got, tol):
    check(what, got <= tol, f'{got:.2e} (tol {tol:g})')


def read_tec(path):
    """Minimal reader for the Tecplot BLOCK files MOSE writes."""
    with open(path) as fh:
        lines = fh.read().split('\n')
    header = lines[0].split('=', 1)[1]
    names = re.findall(r'"([^"]+)"', header) or re.split(r'[\s,]+', header.strip())
    zone = lines[1]
    I, J, K = (int(re.search(rf'{c}=\s*(\d+)', zone).group(1)) for c in 'IJK')
    # the numbers of the zone; the auxiliary data that MOSE writes after it (DATASETAUXDATA, the counters of a restart)
    # are not numbers of the zone
    vals = []
    for line in lines[2:]:
        if line.lstrip().upper().startswith(('DATASETAUXDATA', 'AUXDATA')):
            continue
        vals.extend(float(v) for v in line.split())
    n_node = I * J * K
    n_cell = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)
    out, pos = {}, 0
    for idx, name in enumerate(names):
        n = n_node if idx < 3 else n_cell
        out[name] = np.array(vals[pos:pos + n])
        pos += n
    return out, (I, J, K)


def read_case():
    cfg = {}
    for line in open('INPUT/case.txt'):
        key, val = line.split(None, 1)
        cfg[key] = val.strip()
    return cfg


def main():
    cfg = read_case()
    Rii0, omega0 = float(cfg['Rii']), float(cfg['omega'])
    k0 = 1.5 * Rii0
    print(f'RSM decay [{cfg["model"]}]')

    field, (I, J, K) = read_tec('OUTPUT/field.tec')
    nx, ny = I - 1, J - 1
    cell = lambda name: field[name].reshape(ny, nx)        # [j, i]
    xn = field['x'].reshape(K, J, I)[0, 0, :]
    x = 0.5 * (xn[:-1] + xn[1:])

    rho = cell('rho(1)')
    u, v, w = cell('u'), cell('v'), cell('w')
    R = {c: cell("r" + c) / rho for c in ("u'u'", "v'v'", "w'w'", "u'v'", "u'w'", "v'w'")}
    omega = cell('rhoOmega') / rho
    k = 0.5 * (R["u'u'"] + R["v'v'"] + R["w'w'"])
    U = float(u.mean())
    inner = slice(4, nx - 2)

    # 1. decay along x (mean over the duct height). The cells next to the inlet and to the
    #    outlet are first-order (zero-gradient ghost cells), hence the looser bound on them.
    f = 1.0 + BETA * omega0 * x / U
    k_ref, omega_ref = k0 * f ** (-C_MU / BETA), omega0 / f
    err_k, err_w = np.abs(k.mean(axis=0) / k_ref - 1.0), np.abs(omega.mean(axis=0) / omega_ref - 1.0)
    below('decay of k follows the closed form (max rel. error, inner cells)', float(err_k[inner].max()), 2.0e-3)
    below('decay of omega follows the closed form (max rel. error, inner cells)', float(err_w[inner].max()), 2.0e-3)
    below('decay of k and omega, boundary cells included (max rel. error)',
          float(max(err_k.max(), err_w.max())), 2.0e-2)
    check('the turbulence does decay', k.mean(axis=0)[-1] < 0.7 * k0,
          f'k_out/k_0 = {k.mean(axis=0)[-1] / k0:.3f}')

    # 2. isotropy
    below("normal stresses stay equal (max |Rii - 2k/3| / k)",
          float(max(np.abs(R[c] - 2.0 * k / 3.0).max() for c in ("u'u'", "v'v'", "w'w'")) / k.min()), 1.0e-3)
    below('shear stresses stay zero (max |Rij| / k)',
          float(max(np.abs(R[c]).max() for c in ("u'v'", "u'w'", "v'w'")) / k.min()), 1.0e-6)

    # 3. uniformity across the duct
    below('no radial velocity (max |v| / U)', float(np.abs(v).max() / U), 1.0e-6)
    below('no azimuthal velocity (max |w| / U)', float(np.abs(w).max() / U), 1.0e-6)
    below('k uniform across the duct (max rel. spread)',
          float(np.abs(k / k.mean(axis=0) - 1.0).max()), 1.0e-6)

    # 4. energy: H + k + R11 is constant along x
    cp = cell('g') * cell('R') / (cell('g') - 1.0)
    H = (cp * cell('T') + 0.5 * (u ** 2 + v ** 2 + w ** 2)).mean(axis=0)
    km, R11 = k.mean(axis=0), R["u'u'"].mean(axis=0)
    dH, dk = H[-1] - H[0], km[-1] - km[0]
    total = H + km + R11
    below('total enthalpy + k + R11 conserved along x (max change / |dk|, inner cells)',
          float(np.abs(total[inner] - total[inner][0]).max() / abs(dk)), 1.0e-3)
    below('total enthalpy + k + R11 conserved along x, boundary cells included',
          float(np.abs(total - total[inner][0]).max() / abs(dk)), 3.0e-2)
    check('the mean flow gains what the turbulence loses (dH/dk = -5/3)',
          abs(dH / dk + 5.0 / 3.0) <= 0.05, f'dH/dk = {dH / dk:.4f}')

    print(f'   {checks - len(failures)}/{checks} checks passed')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
