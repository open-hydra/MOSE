#!/usr/bin/env python3
"""Derive the manufactured-solution source terms for the 2D compressible
Navier-Stokes equations (calorically-perfect ideal gas, constant mu/k, Stokes
hypothesis) and emit Fortran expressions for src/app/mms.f90 and Python for
build_ic.py.

Manufactured PRIMITIVE field (periodic on [0,L]^2):
    rho = rho0 + arho*sin(om x)*cos(om y)
    u   = u0   + au  *sin(om x)*sin(om y)
    v   = v0   + av  *cos(om x)*cos(om y)
    p   = p0   + ap  *cos(om x)*sin(om y)
It is a *steady* exact solution once the analytic residual S = div(Fc - Fv)
is added as a source.  Governing sign in MOSE: R = CONV - DIFF - SOURCE, so the
driver does  R -= S*vol.
"""
import sympy as sp

x, y = sp.symbols('x y', real=True)

# ---- baked constants (MUST match the solver's thermo/transport) -------------
Runiv = sp.Rational(831451, 100)          # 8314.51 J/(kmol K)  (FLINT)
W     = sp.Float('28.970418')             # air molecular weight (phase.txt)
R     = Runiv / W                         # specific gas constant
cp    = sp.Float('1004.5')                # constant cp (thermo.dat)
gamma = cp / (cp - R)                     # -> 1.4
mu    = sp.Float('10.0')                  # constant dynamic viscosity (flat transport.dat); Re~5
Pr    = sp.Float('0.72')                  # laminar Prandtl (input.ini Prl)
kcond = mu * cp / Pr                       # constant thermal conductivity
lam   = sp.Rational(-2, 3) * mu            # Stokes bulk viscosity

L   = sp.Float('1.0')
om  = 2*sp.pi / L

# amplitudes / means
rho0, arho = sp.Float('1.0'),   sp.Float('0.1')
u0,   au   = sp.Float('40.0'),  sp.Float('8.0')
v0,   av   = sp.Float('30.0'),  sp.Float('8.0')
p0,   ap   = sp.Float('8.0e3'), sp.Float('8.0e2')   # low sound speed -> M~0.47 (avoid low-Mach degradation)

rho = rho0 + arho*sp.sin(om*x)*sp.cos(om*y)
u   = u0   + au  *sp.sin(om*x)*sp.sin(om*y)
v   = v0   + av  *sp.cos(om*x)*sp.cos(om*y)
p   = p0   + ap  *sp.cos(om*x)*sp.sin(om*y)

T = p/(rho*R)
E = p/(gamma-1) + sp.Rational(1,2)*rho*(u**2+v**2)

ux, uy = sp.diff(u,x), sp.diff(u,y)
vx, vy = sp.diff(v,x), sp.diff(v,y)
div    = ux + vy

tau_xx = lam*div + 2*mu*ux
tau_yy = lam*div + 2*mu*vy
tau_xy = mu*(uy+vx)

Tx, Ty = sp.diff(T,x), sp.diff(T,y)

# total fluxes  F = Fc - Fv   (viscous work + heat conduction in energy)
Fx_rho = rho*u
Fy_rho = rho*v
Fx_mx  = rho*u*u + p - tau_xx
Fy_mx  = rho*u*v     - tau_xy
Fx_my  = rho*u*v     - tau_xy
Fy_my  = rho*v*v + p - tau_yy
Fx_E   = u*(E+p) - (u*tau_xx + v*tau_xy) - kcond*Tx
Fy_E   = v*(E+p) - (u*tau_xy + v*tau_yy) - kcond*Ty

S_rho = sp.diff(Fx_rho,x) + sp.diff(Fy_rho,y)
S_mx  = sp.diff(Fx_mx, x) + sp.diff(Fy_mx, y)
S_my  = sp.diff(Fx_my, x) + sp.diff(Fy_my, y)
S_E   = sp.diff(Fx_E,  x) + sp.diff(Fy_E,  y)

print("! R     =", sp.N(R, 12))
print("! gamma =", sp.N(gamma, 12))
print("! kcond =", sp.N(kcond, 12))
print()

exprs = [('S_rho', S_rho), ('S_mx', S_mx), ('S_my', S_my), ('S_E', S_E)]

# --- independent finite-difference check of the flux divergence --------------
def num(expr, xv, yv):
    return float(expr.subs({x: xv, y: yv}))
def fd_div(Fx, Fy, xv, yv, h=1e-6):
    return (num(Fx, xv+h, yv) - num(Fx, xv-h, yv))/(2*h) + \
           (num(Fy, xv, yv+h) - num(Fy, xv, yv-h))/(2*h)
xv, yv = 0.123, 0.456
fd = {
    'S_rho': fd_div(Fx_rho, Fy_rho, xv, yv),
    'S_mx' : fd_div(Fx_mx,  Fy_mx,  xv, yv),
    'S_my' : fd_div(Fx_my,  Fy_my,  xv, yv),
    'S_E'  : fd_div(Fx_E,   Fy_E,   xv, yv),
}
print("--- analytic vs finite-difference divergence check ---")
ok = True
for name, e in exprs:
    a = num(e, xv, yv); f = fd[name]
    rel = abs(a-f)/(abs(a)+1e-30)
    flag = 'OK' if rel < 1e-5 else 'FAIL'
    if rel >= 1e-5: ok = False
    print(f"  {name:6s} analytic={a: .6e}  fd={f: .6e}  rel={rel:.2e}  {flag}")
print("CHECK", "PASSED" if ok else "FAILED")
print()

# --- Fortran emission (pi as a declared local, not sympy.pi) -----------------
pi_s = sp.Symbol('pi')
exprs_pi = [(n, e.subs(sp.pi, pi_s)) for n, e in exprs]
subs, reduced = sp.cse([e for _, e in exprs_pi],
                       optimizations='basic', symbols=sp.numbered_symbols('t'))
allt = sorted({str(s) for s, _ in subs})
with open('mms_fortran.f90', 'w') as fh:
    for s, expr in subs:
        fh.write(f"      {sp.fcode(expr, assign_to=str(s), source_format='free', standard=95, contract=False)}\n")
    for (name, _), red in zip(exprs_pi, reduced):
        fh.write(f"      {sp.fcode(red, assign_to=name, source_format='free', standard=95, contract=False)}\n")
print("temporaries:", ', '.join(allt))
print(open('mms_fortran.f90').read())
