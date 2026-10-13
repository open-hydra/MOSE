#!/usr/bin/env python3
"""Derive the manufactured-solution source terms of the fields periodic3d, wall2d
and wall3d of src/app/mms.f90 (chosen there with MOSE_MMS_FIELD) and emit the
Fortran expressions of their source subroutines.

Same gas and transport as mms_gen.py (the field periodic2d): calorically perfect
air (R = 8314.46261815324/28.970418, cp = 1004.5), constant mu = 10, Pr = 0.72, Stokes
hypothesis. Manufactured PRIMITIVE fields (om = 2 pi; R = Rgas, T = p/(rho R)):

  periodic3d  periodic on [0,1]^3, the field of mms_gen.py extended in z:
      rho = 1 + 0.1 sin(om x) cos(om y) cos(om z)   u = 40 + 8 sin(om x) sin(om y) cos(om z)
      v   = 30 + 8 cos(om x) cos(om y) cos(om z)    w = 20 + 8 cos(om x) sin(om y) sin(om z)
      p   = 8000 + 800 cos(om x) sin(om y) cos(om z)
  wall2d, wall3d  periodic in x (and z), walls at y = 0 and y = 1:
      u = sin(pi y) (40 + 8 sin(om x) c_z),  v = 8 sin(pi y) cos(om x) c_z,
      w = sin(pi y) (20 + 8 cos(om x) sin(om z))  (wall3d only),
      T = T_w(x, z) + 5.6 sin(pi y),  T_w = 28 + 1.4 cos(om x) c_z,
      p = 8000 + 800 cos(om x) cos(pi y) c_z,  rho = p/(R T),
    with c_z = cos(om z) in wall3d and 1 in wall2d: zero velocity and T = T_w on
    both walls, and d2u/dy2 = d2T/dy2 = 0 there.

Each field becomes a steady exact solution once the analytic residual
S = div(Fc - Fv) is added as a source (MOSE: R = CONV - DIFF - SOURCE, so the
driver does R -= S*vol). As in mms_gen.py, the analytic divergence is checked
against a centred finite-difference divergence of the fluxes, here at 20 random
points in 30-digit arithmetic; the script stops if a relative difference exceeds
1e-6.

Usage: python3 mms_gen_fields.py FIELD   (FIELD: periodic3d, wall2d or wall3d)
writes mms_<FIELD>.f90, the body of the source subroutine of that field.
"""
import random
import sys

import sympy as sp

x, y, z = sp.symbols('x y z', real=True)
Runiv = sp.Rational(831446261815324, 10**11)  # 8314.46261815324 J/(kmol K) (FLINT, exact SI N_A*k_B), as mms_gen.py
W = sp.Float('28.970418')
R = Runiv / W
cp = sp.Float('1004.5')
gamma = cp / (cp - R)
mu = sp.Float('10.0')
Pr = sp.Float('0.72')
kcond = mu * cp / Pr
lam = sp.Rational(-2, 3) * mu
om = 2 * sp.pi


def field(name):
    if name == 'periodic3d':
        rho = 1 + sp.Float('0.1') * sp.sin(om * x) * sp.cos(om * y) * sp.cos(om * z)
        u = 40 + 8 * sp.sin(om * x) * sp.sin(om * y) * sp.cos(om * z)
        v = 30 + 8 * sp.cos(om * x) * sp.cos(om * y) * sp.cos(om * z)
        w = 20 + 8 * sp.cos(om * x) * sp.sin(om * y) * sp.sin(om * z)
        p = 8000 + 800 * sp.cos(om * x) * sp.sin(om * y) * sp.cos(om * z)
        return rho, u, v, w, p
    if name in ('wall2d', 'wall3d'):
        three = name == 'wall3d'
        cz = sp.cos(om * z) if three else sp.Integer(1)
        sy = sp.sin(sp.pi * y)
        u = sy * (40 + 8 * sp.sin(om * x) * cz)
        v = 8 * sy * sp.cos(om * x) * cz
        w = sy * (20 + 8 * sp.cos(om * x) * sp.sin(om * z)) if three else sp.Integer(0)
        Tw = sp.Float('28.0') + sp.Float('1.4') * sp.cos(om * x) * cz
        T = Tw + sp.Float('5.6') * sy
        p = 8000 + 800 * sp.cos(om * x) * sp.cos(sp.pi * y) * cz
        rho = p / (R * T)
        return rho, u, v, w, p
    raise SystemExit('unknown field: ' + name)


def fluxes(rho, u, v, w, p):
    """Total fluxes F = Fc - Fv of the five conservation equations, F[eq][direction]."""
    T = p / (rho * R)
    E = p / (gamma - 1) + sp.Rational(1, 2) * rho * (u ** 2 + v ** 2 + w ** 2)
    U = [u, v, w]
    X = [x, y, z]
    g = [[sp.diff(U[a], X[b]) for b in range(3)] for a in range(3)]     # g[a][b] = d u_a / d x_b
    div = g[0][0] + g[1][1] + g[2][2]
    tau = [[mu * (g[a][b] + g[b][a]) + (lam * div if a == b else 0) for b in range(3)] for a in range(3)]
    gT = [sp.diff(T, X[b]) for b in range(3)]
    F = {'rho': [rho * U[b] for b in range(3)]}
    for a, nm in enumerate(('mx', 'my', 'mz')):
        F[nm] = [rho * U[a] * U[b] + (p if a == b else 0) - tau[a][b] for b in range(3)]
    F['E'] = [U[b] * (E + p) - sum(U[a] * tau[a][b] for a in range(3)) - kcond * gT[b] for b in range(3)]
    return F


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'periodic3d'
    rho, u, v, w, p = field(name)
    F = fluxes(rho, u, v, w, p)
    X = [x, y, z]
    names = ['rho', 'mx', 'my', 'mz', 'E']
    S = {n: sum(sp.diff(F[n][b], X[b]) for b in range(3)) for n in names}
    print('! R     =', sp.N(R, 15))
    print('! gamma =', sp.N(gamma, 15))
    print('! kcond =', sp.N(kcond, 15))
    # check: analytic divergence against centred finite differences (step 1e-7, 30 digits)
    random.seed(20261007)
    ok = True
    h = sp.Float('1e-7', 30)
    worst = 0.0
    for _ in range(20):
        pt = {x: sp.Float(random.uniform(0.05, 0.95), 30), y: sp.Float(random.uniform(0.05, 0.95), 30),
              z: sp.Float(random.uniform(0.05, 0.95), 30)}
        for n in names:
            a = sp.N(S[n].subs(pt), 30)
            fd = 0
            for b in range(3):
                pp, pm = dict(pt), dict(pt)
                pp[X[b]] += h
                pm[X[b]] -= h
                fd += (sp.N(F[n][b].subs(pp), 30) - sp.N(F[n][b].subs(pm), 30)) / (2 * h)
            scale = max(abs(a), sp.Float(1e-30))
            rel = float(abs(a - fd) / scale) if abs(a) > 1e-12 else float(abs(a - fd))
            worst = max(worst, rel)
            if rel > 1e-6:
                ok = False
    print(f'--- finite-difference check: largest relative difference {worst:.2e} over 20 points x 5 equations')
    print('CHECK', 'PASSED' if ok else 'FAILED')
    if not ok:
        raise SystemExit(1)
    # Fortran emission (pi as a declared local, as in mms_gen.py)
    pi_s = sp.Symbol('pi')
    exprs = [(f'S_{n}', S[n].subs(sp.pi, pi_s)) for n in names]
    subs, reduced = sp.cse([e for _, e in exprs], optimizations='basic', symbols=sp.numbered_symbols('t'))
    with open(f'mms_{name}.f90', 'w') as fh:
        fh.write(f'      ! {name}: {len(subs)} temporaries t0..t{len(subs) - 1}\n')
        for s, e in subs:
            fh.write(f"      {sp.fcode(e, assign_to=str(s), source_format='free', standard=95, contract=False)}\n")
        for (n, _), red in zip(exprs, reduced):
            fh.write(f"      {sp.fcode(red, assign_to=n, source_format='free', standard=95, contract=False)}\n")
    print(f'wrote mms_{name}.f90 with {len(subs)} temporaries')


if __name__ == '__main__':
    main()
