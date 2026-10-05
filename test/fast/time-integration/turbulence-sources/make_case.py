"""Inputs of the TurbulenceSources fast test, and the readers of its outputs.

One geometry: a 2-D channel, periodic in x (4 cells of 1 mm, connection 201 on itself), an adiabatic no-slip wall at
y = 0 (301) and a symmetry plane at y = H (300), null k faces (0). Air of test/common/Air-viscous (mu and k constant,
Pr = 0.72), 300 K, 1 bar, flowing along x with u(y) = U tanh(y/delta_u), U = 50 m/s, delta_u = 0.2 mm, and no force
that drives it: the wall slows the flow down and the turbulence variables relax. The state depends on y only.

Cases, written into INPUT/ (one at a time) by `python3 -B make_case.py CASE`:
  order-SA   Spalart-Allmaras, nu~ = 2 kappa u_tau y exp(-y/delta), twice the near-wall equilibrium, delta = 0.2 mm;
  order-SST  SST with omega-wall-bc = asymptotic: k ~ y^2 near the wall, omega = 6 nu/(beta_1 y^2) near the wall plus a
             log-layer part u_tau/(sqrt(beta*) kappa y);
             both on a mesh with cells of 4 um at the wall and of 1.5 um at the symmetry plane. Turbulence and shear
             vanish at the symmetry plane, so its thin cells keep their state, and they set the global time step: the
             step does not depend on the turbulence, and N steps at cfl c, 2N at c/2 and 4N at c/4 end at the same time.
  pos-SA     Spalart-Allmaras with nu~ a hundred times the near-wall equilibrium within 20 um of the wall;
  pos-SST    SST with omega-wall-bc = practical, omega ten times the asymptotic value within 20 um of the wall;
             both on a mesh with cells of 2 um at the wall and no thin cells elsewhere: the turbulence decays near the
             wall faster than the transport acts, so without a step limit an explicit stage would make nu~, k or omega
             negative.
The turbulence profiles fade out as exp(-(y/0.4 mm)^4): beyond 1 mm from the wall nu~, k and omega take uniform
free-stream values (rho nu~ = rho k = 2e-10, twice the realizability floor; omega = 1e-2 1/s), so the turbulence and its
stiffness vanish there and nothing diffuses into the thin cells at the symmetry plane.

Usage (from the test directory):
  python3 -B make_case.py CASE                         INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py order MODEL SCHEME F1 F2 F4  time order of the turbulence variables from the final fields
                                                       of the runs at cfl c, c/2 and c/4 (exit status 0 when it is at
                                                       least 1.8 with RK2, 2.7 with RK3)
  python3 -B make_case.py pos CASE CFL VNN LOG         positivity check of a run (exit status 0 when it holds)
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU, KL, CP = 1.1858686e-5, 0.016544513989, 1004.5   # test/common/Air-viscous
PRT = 0.85                          # default turbulent Prandtl number
T0, P0, U = 300.0, 1.0e5, 50.0
RHO0 = P0 / (RGAS * T0)
NU = MU / RHO0
NI, DX, DZ = 4, 1.0e-3, 1.0e-3
DELTA_U = 2.0e-4
UTAU = math.sqrt(NU * U / DELTA_U)  # from the wall shear of the tanh profile
KAPPA, BETA1, BETASTAR = 0.41, 0.075, 0.09
CV1, SIGMA = 7.1, 2.0 / 3.0         # Spalart-Allmaras
NUT_INF = K_INF = 2.0e-10 / RHO0
OMEGA_INF = 1.0e-2


def fade(y):
    return math.exp(-(y / 4.0e-4) ** 4)


def geometric(first, n, ratio):
    return [first * ratio ** k for k in range(n)]


def nodes_y(kind):
    """Wall-normal nodes: geometric growth (ratio 1.25) from the wall, a uniform core, and for the order cases a
    geometric growth from thin cells at the symmetry plane."""
    if kind == "order":
        sizes = geometric(4.0e-6, 16, 1.25) + [1.2e-4] * 8 + geometric(1.5e-6, 20, 1.25)[::-1]
    else:
        sizes = geometric(2.0e-6, 18, 1.25) + [1.0e-4] * 16
    y = [0.0]
    for s in sizes:
        y.append(y[-1] + s)
    return y


def profiles(case, yc):
    """Cell values (rho, u, rans...) of a case at the cell centres yc."""
    u = [U * math.tanh(y / DELTA_U) for y in yc]
    if case.endswith("SA"):
        if case == "order-SA":
            nut = [2.0 * KAPPA * UTAU * y * math.exp(-y / 2.0e-4) * fade(y) + NUT_INF for y in yc]
        else:
            nut = [(100.0 * KAPPA * UTAU * y * math.exp(-y / 2.0e-5) + 2.0 * KAPPA * UTAU * y * math.exp(-y / 2.0e-4))
                   * fade(y) + NUT_INF for y in yc]
        return u, [[RHO0 * v for v in nut]]
    # SST: k ~ y^2 near the wall, peak at 0.1 mm; omega: viscous sublayer plus a log-layer part
    ka = UTAU ** 2 / math.sqrt(BETASTAR)
    k = [ka * (y / 1.0e-4) ** 2 * math.exp(2.0 * (1.0 - y / 1.0e-4)) * fade(y) + K_INF for y in yc]
    wall = 1.0 if case == "order-SST" else 10.0
    om = [(6.0 * NU / (BETA1 * y ** 2) * (1.0 + (wall - 1.0) * math.exp(-y / 2.0e-5))
           + UTAU / (math.sqrt(BETASTAR) * KAPPA * y)) * fade(y) + OMEGA_INF for y in yc]
    return u, [[RHO0 * v for v in k], [RHO0 * v for v in om]]


def write_case(case):
    kind = case.split("-")[0]
    y = nodes_y(kind)
    nj = len(y) - 1
    yc = [0.5 * (y[j] + y[j + 1]) for j in range(nj)]
    u, rans = profiles(case, yc)
    names = ' "rhoNuTilde"' if len(rans) == 1 else ' "rhoK" "rhoOmega"'
    nvar = 8 + len(rans)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"' + names + "\n")
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={nj + 1}, K=2, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-{nvar}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for fun in (lambda i, j, k: i * DX, lambda i, j, k: y[j], lambda i, j, k: k * DZ):
            for k in range(2):
                for j in range(nj + 1):
                    for i in range(NI + 1):
                        f.write(f" {fun(i, j, k):.16E}\n")
        cols = [[RHO0] * nj, u, [0.0] * nj, [0.0] * nj, [P0] * nj] + rans
        for col in cols:
            for j in range(nj):
                f.write(f" {col[j]:.16E}\n" * NI)
    out = []
    for j in range(1, nj + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, 1, j, 1, 1, 201))
        out.append('%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, NI, j, 1, 2, 1, 0, 0, 1))
    for j in range(1, nj + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, NI, j, 1, 2, 201))
        out.append('%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, 1, j, 1, 1, 1, 0, 0, 1))
    for i in range(1, NI + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, 1, 1, 3, 301))
        out.append('    0.000000E+00,    0.000000E+00,    0.000000E+00,\n')
    for i in range(1, NI + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, nj, 1, 4, 300))
    for face in (5, 6):
        for j in range(1, nj + 1):
            for i in range(1, NI + 1):
                out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, 0))
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def transport_step(case, cfl, vnn):
    """Smallest local time step of the initial state from convection and diffusion only (the rule of Compute_dt
    without the turbulence sources): 1/dt = sum_d (|u_d| + c)/h_d / cfl + D_max sum_d 1/h_d^2 / vnn."""
    y = nodes_y(case.split("-")[0])
    yc = [0.5 * (y[j] + y[j + 1]) for j in range(len(y) - 1)]
    u, rans = profiles(case, yc)
    c = math.sqrt(CP / (CP - RGAS) * P0 / RHO0)
    cv = CP - RGAS
    best = math.inf
    for j in range(len(yc)):
        h = y[j + 1] - y[j]
        if len(rans) == 1:
            chi = rans[0][j] / MU
            mut = rans[0][j] * chi ** 3 / (chi ** 3 + CV1 ** 3)
            extra = (MU + rans[0][j]) / (SIGMA * RHO0)
        else:
            mut = rans[0][j] / (rans[1][j] / RHO0)
            extra = 0.0
        dmax = max(4.0 / 3.0 * (MU + mut) / RHO0, (KL + mut * CP / PRT) / (RHO0 * cv), extra)
        rate = ((abs(u[j]) + c) / DX + c / h) / cfl + dmax * (1.0 / DX ** 2 + 1.0 / h ** 2) / vnn
        best = min(best, 1.0 / rate)
    return best


def read_field(path):
    """Cell values and solution time of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    time = float(re.search(r'SOLUTIONTIME=\s*([-+0-9.EeDd]+)', lines[1]).group(1).replace('D', 'E'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, time


def order(model, scheme, f1, f2, f4):
    """Observed order p = log2(e1/e2), e1 = max|q(c) - q(c/2)|, e2 = max|q(c/2) - q(c/4)|, both over max|q(c/4)|."""
    least = {"RK2": 1.8, "RK3": 2.7}[scheme]
    (a, t1), (b, t2), (c, t4) = read_field(f1), read_field(f2), read_field(f4)
    dt = max(abs(t1 - t4), abs(t2 - t4)) / t4
    print(f"    final times {t1:.15e} {t2:.15e} {t4:.15e} s (largest relative difference {dt:.1e})")
    ok = dt < 1.0e-9
    if not ok:
        print("    [FAIL] the runs do not end at the same time: the order cannot be read")
    turb = ["rhoNuTilde"] if model == "SA" else ["rhoK", "rhoOmega"]
    for name in turb + ["rho(1)", "u", "p"]:
        scale = max(abs(v) for v in c[name])
        e1 = max(abs(x - y) for x, y in zip(a[name], b[name])) / scale
        e2 = max(abs(x - y) for x, y in zip(b[name], c[name])) / scale
        p = math.log2(e1 / e2) if e1 > 0 and e2 > 0 else float("nan")
        tag = "turbulence" if name in turb else "mean flow, informative"
        print(f"    {name:9s} e(c, c/2) {e1:.3e}  e(c/2, c/4) {e2:.3e}  order {p:5.2f}   ({tag})")
        if name in turb and not p >= least:
            ok = False
    sys.exit(0 if ok else 1)


def pos(case, cfl, vnn, log):
    """The run ends, the solver reports no negative turbulence value after a stage, and its step was limited by the
    turbulence sources at some point (smallest step below half the transport step of the initial state)."""
    text = open(log, errors="replace").read()
    ok = True
    m = re.search(r'Turbulence variables below zero after a stage[^\n]*?(\d+)\s*$', text, re.M)
    if m is None:
        print("    [FAIL] no count of negative turbulence values in the solver output")
        ok = False
    else:
        n = int(m.group(1))
        print(f"    negative nu~, k or omega after a stage: {n}")
        ok = ok and n == 0
    steps = [float(s) for s in re.findall(r'Delta t =\s*([-+0-9.Ee]+)', text)]
    if steps:
        tr = transport_step(case, cfl, vnn)
        ratio = min(steps) / tr
        print(f"    smallest step {min(steps):.4e} s, transport step of the initial state {tr:.4e} s, ratio {ratio:.3f}")
        if not ratio < 0.5:
            print("    [FAIL] the step was never limited by the turbulence sources: the case does not test them")
            ok = False
    else:
        print("    [FAIL] no step in the solver output")
        ok = False
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "order":
        order(*sys.argv[2:7])
    elif what == "pos":
        pos(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif what in ("order-SA", "order-SST", "pos-SA", "pos-SST"):
        write_case(what)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
