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
             negative. Their first step must be the rule of Compute_dt on the initial state with the bound of the
             stiffness of the destruction that the state of each cell gives alone (first_step): SA 2 c_w1
             (1 + c_w3^6)^(1/6) nu~/y^2, SST max(beta*, 2 max(beta_1, beta_2)) omega.
The turbulence profiles fade out as exp(-(y/0.4 mm)^4): beyond 1 mm from the wall nu~, k and omega take uniform
free-stream values (rho nu~ = rho k = 2e-10, twice the realizability floor; omega = 1e-2 1/s), so the turbulence and its
stiffness vanish there and nothing diffuses into the thin cells at the symmetry plane.

A second geometry, for the vortex-stretching function f_beta of Wilcox 2006 (case wilcox): a cube of NW^3 cells of DXW (dyadic:
cubes to the bit), one block, extrapolation (400) on all faces, the same air at 300 K and 1 bar, the velocity u_i = GW_ij x_j
about the centre of the cube with GW zero on the diagonal (each component varies along the other two coordinates only:
divergence-free, and the normal velocity of a face is the same on its two sides, so with uniform density and pressure the HLLC
flux of rho omega through the faces of a cell adds up to zero), uniform k = KW and omega = OMW, above the limiter of omega
(omega-hat = omega). One explicit step: the change of rho omega is dt (P_omega - beta_0 f_beta rho omega^2), with
P_omega = (gamma omega/k) tau_ij du_i/dx_j = 2 gamma rho S_ij S_ij, chi_w = |W_ij W_jk S^_ki|/(beta* omega)^3,
S^_ki = S_ki - (1/2) du_m/dx_m delta_ki = S_ki here (NASA Turbulence Modeling Resource, Wilcox 2006 page; the contraction of
chi_w is that of MOSE and of Wilcox's model, the page prints W_ij W_ij S^_ki), f_beta = (1 + 85 chi_w)/(1 + 100 chi_w).

A third geometry, for the asymptotic omega wall value of the k-omega models (cases wall-Wilcox2006 and wall-SST,
omega-wall-bc = asymptotic): the channel above with NJW uniform rows of DYW (dyadic) and square cells, the air at rest at
300 K and 1 bar, uniform k = KWALL and, everywhere, omega = 6 nu/(beta d^2), the asymptotic wall value of the model (NASA
Turbulence Modeling Resource: beta_0 = 0.0708 on the Wilcox 2006 page, beta_1 = 0.075 on the SST page), d = DYW/2 the
distance of the centre of a wall cell from the wall. At rest there is no production, and with omega uniform and equal to its
wall value no omega diffuses through the wall and the cross-diffusion is zero: one explicit step changes rho omega of the
wall cells by the destruction alone, -dt beta rho omega^2 (f_beta = 1 without vorticity; the blending of SST is 1 in those
cells, where 500 nu/(d^2 omega) = 6.25, so beta = beta_1).

Usage (from the test directory):
  python3 -B make_case.py CASE                         INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py order MODEL SCHEME F1 F2 F4  time order of the turbulence variables from the final fields
                                                       of the runs at cfl c, c/2 and c/4 (exit status 0 when it is at
                                                       least 1.8 with RK2, 2.7 with RK3)
  python3 -B make_case.py pos CASE CFL VNN LOG         positivity check of a run (exit status 0 when it holds)
  python3 -B make_case.py count LOG some|none          the count of the cell updates whose stage made nu~, k or omega
                                                       negative at the end of LOG: one line with a count above zero
                                                       (some), or no line (none); exit status 0 when it is so
  python3 -B make_case.py wilcox-check FIELD           the change of rho omega in the one step of the case wilcox
                                                       (exit status 0 when it is the one of f_beta within RTOL_W)
  python3 -B make_case.py wall-check MODEL FIELD       the change of rho omega in the wall cells in the one step of the
                                                       case wall-MODEL (exit status 0 when it is the destruction alone
                                                       within RTOL_W)
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
CB1, CB2, CW3 = 0.1355, 0.622, 2.0
CW1 = CB1 / KAPPA ** 2 + (1.0 + CB2) / SIGMA
BETA2 = 0.0828                      # SST, outer beta
NUT_INF = K_INF = 2.0e-10 / RHO0
OMEGA_INF = 1.0e-2
NW, DXW = 6, 2.0 ** -10             # case wilcox: cube of NW^3 cells of DXW m
GW = ((0.0, 1000.0, 600.0), (-400.0, 0.0, 800.0), (-500.0, 300.0, 0.0))   # du_i/dx_j, 1/s, zero on the diagonal
KW, OMW = 1.0, 5000.0               # m2/s2 and 1/s, uniform; omega-hat limit C_lim sqrt(2 S:S/beta*) = 3666 1/s
GAMMA_W, BETA0_W, CLIM_W = 13.0 / 25.0, 0.0708, 7.0 / 8.0                  # Wilcox 2006 (NASA TMR)
RTOL_W = 1.0e-6
NJW, DYW, KWALL = 4, 2.0 ** -17, 1.0e-2   # cases wall-MODEL: rows, cell size (m), k (m2/s2)
BETA_WALL = {"Wilcox2006": BETA0_W, "SST": BETA1}   # beta of the omega equation next to the wall


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


def write_wilcox():
    c = 0.5 * NW * DXW
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    order = [(i, j, k) for k in range(NW) for j in range(NW) for i in range(NW)]
    vel = {}
    for (i, j, k) in order:
        x = ((i + 0.5) * DXW - c, (j + 0.5) * DXW - c, (k + 0.5) * DXW - c)
        vel[i, j, k] = [sum(GW[m][n] * x[n] for n in range(3)) for m in range(3)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoK" "rhoOmega"\n')
        f.write(f" ZONE  T = B1-IG, I={NW + 1}, J={NW + 1}, K={NW + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-10]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(NW + 1):
                for j in range(NW + 1):
                    for i in range(NW + 1):
                        f.write(f" {(DXW * i, DXW * j, DXW * k)[axis]:.17E}\n")
        f.write(f" {RHO0:.17E}\n" * NW ** 3)
        for axis in range(3):
            f.write("".join(f" {vel[q][axis]:.17E}\n" for q in order))
        f.write(f" {P0:.17E}\n" * NW ** 3)
        f.write(f" {RHO0 * KW:.17E}\n" * NW ** 3)
        f.write(f" {RHO0 * OMW:.17E}\n" * NW ** 3)
    rows = []
    for face in range(1, 7):
        for b in range(1, NW + 1):
            for m in range(1, NW + 1):
                i, j, k = {1: (1, m, b), 2: (NW, m, b), 3: (m, 1, b), 4: (m, NW, b), 5: (m, b, 1), 6: (m, b, NW)}[face]
                rows.append(f"       1     {i:3d}     {j:3d}     {k:3d}       {face}     400\n")
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(rows))


def wilcox_source():
    """Source of rho omega of the case wilcox (NASA TMR, Wilcox 2006), f_beta and chi_w."""
    S = [[0.5 * (GW[i][j] + GW[j][i]) for j in range(3)] for i in range(3)]
    W = [[0.5 * (GW[i][j] - GW[j][i]) for j in range(3)] for i in range(3)]
    div = GW[0][0] + GW[1][1] + GW[2][2]
    Sh = [[S[i][j] - (0.5 * div if i == j else 0.0) for j in range(3)] for i in range(3)]
    ss = sum(S[i][j] ** 2 for i in range(3) for j in range(3))
    sbar = [[S[i][j] - (div / 3.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
    if OMW <= CLIM_W * math.sqrt(2.0 * sum(v ** 2 for row in sbar for v in row) / BETASTAR):
        sys.exit("the case puts omega under the limiter of omega-hat")
    num = sum(W[i][j] * W[j][k] * Sh[k][i] for i in range(3) for j in range(3) for k in range(3))
    chi = abs(num / (BETASTAR * OMW) ** 3)
    fbeta = (1.0 + 85.0 * chi) / (1.0 + 100.0 * chi)
    return 2.0 * GAMMA_W * RHO0 * ss - BETA0_W * fbeta * RHO0 * OMW ** 2, fbeta, chi


def wilcox_check(path):
    try:
        z, t = read_field(path)
    except (OSError, ValueError, KeyError, IndexError, AttributeError) as err:
        print(f"    no field of the case wilcox ({err})")
        sys.exit(1)
    src, fbeta, chi = wilcox_source()
    want = t * src
    got = [v - RHO0 * OMW for v in z["rhoOmega"]]
    if len(got) != NW ** 3 or any(math.isnan(v) for v in got):
        print("    the field is not the case wilcox, or rho omega is not a number")
        sys.exit(1)
    err = max(abs(g - want) for g in got)
    print(f"    dt = {t:.6e} s, chi_w = {chi:.4f}, f_beta = {fbeta:.5f}; change of rho omega expected {want:.6e}, computed "
          f"{min(got):.6e} to {max(got):.6e} kg/(m3 s), largest difference {err / abs(want):.3e} of it")
    sys.exit(0 if err <= RTOL_W * abs(want) else 1)


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


def first_step(case, cfl, vnn):
    """First global step of a time-accurate run of the case: the smallest local step of the initial state, the rule of
    transport_step plus, as Compute_dt adds it, the bound of the stiffness of the destruction from the state of the cell
    alone (RANS_Stiffness_Bound): SA 2 c_w1 (1 + c_w3^6)^(1/6) nu~/y^2, SST max(beta*, 2 max(beta_1, beta_2)) omega,
    y the distance of the cell centre from the wall."""
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
            stiff = 2.0 * CW1 * (1.0 + CW3 ** 6) ** (1.0 / 6.0) * (rans[0][j] / RHO0) / yc[j] ** 2
        else:
            mut = rans[0][j] / (rans[1][j] / RHO0)
            extra = 0.0
            stiff = max(BETASTAR, 2.0 * max(BETA1, BETA2)) * rans[1][j] / RHO0
        dmax = max(4.0 / 3.0 * (MU + mut) / RHO0, (KL + mut * CP / PRT) / (RHO0 * cv), extra)
        rate = ((abs(u[j]) + c) / DX + c / h) / cfl + dmax * (1.0 / DX ** 2 + 1.0 / h ** 2) / vnn + stiff
        best = min(best, 1.0 / rate)
    return best


def omega_wall(model):
    """The asymptotic wall value of omega of a model at the centre of a wall cell, 6 nu/(beta d^2), d = DYW/2."""
    return 6.0 * NU / (BETA_WALL[model] * (0.5 * DYW) ** 2)


def write_wall(model):
    om = omega_wall(model)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoK" "rhoOmega"\n')
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={NJW + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-10]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for fun in (lambda i, j, k: i * DYW, lambda i, j, k: j * DYW, lambda i, j, k: k * DYW):
            for k in range(2):
                for j in range(NJW + 1):
                    for i in range(NI + 1):
                        f.write(f" {fun(i, j, k):.17E}\n")
        for v in (RHO0, 0.0, 0.0, 0.0, P0, RHO0 * KWALL, RHO0 * om):
            f.write(f" {v:.17E}\n" * (NI * NJW))
    out = []
    for j in range(1, NJW + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, 1, j, 1, 1, 201))
        out.append('%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, NI, j, 1, 2, 1, 0, 0, 1))
    for j in range(1, NJW + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, NI, j, 1, 2, 201))
        out.append('%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, 1, j, 1, 1, 1, 0, 0, 1))
    for i in range(1, NI + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, 1, 1, 3, 301))
        out.append('    0.000000E+00,    0.000000E+00,    0.000000E+00,\n')
    for i in range(1, NI + 1):
        out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, NJW, 1, 4, 300))
    for face in (5, 6):
        for j in range(1, NJW + 1):
            for i in range(1, NI + 1):
                out.append('%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, 0))
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def wall_check(model, path):
    try:
        z, t = read_field(path)
    except (OSError, ValueError, KeyError, IndexError, AttributeError) as err:
        print(f"    no field of the case wall-{model} ({err})")
        sys.exit(1)
    om = omega_wall(model)
    want = -t * BETA_WALL[model] * RHO0 * om ** 2
    got = [v - RHO0 * om for v in z["rhoOmega"][:NI]]      # the wall row, j = 1
    if len(z["rhoOmega"]) != NI * NJW or any(math.isnan(v) for v in got):
        print(f"    the field is not the case wall-{model}, or rho omega is not a number")
        sys.exit(1)
    err = max(abs(g - want) for g in got)
    print(f"    dt = {t:.6e} s, omega = {om:.6e} 1/s, beta = {BETA_WALL[model]}; change of rho omega in the wall cells expected "
          f"{want:.6e}, computed {min(got):.6e} to {max(got):.6e} kg/(m3 s), largest difference {err / abs(want):.3e} of it")
    sys.exit(0 if err <= RTOL_W * abs(want) else 1)


def read_field(path):
    """Cell values and solution time of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    time = float(re.search(r'SOLUTIONTIME=\s*([-+0-9.EeDd]+)', lines[1]).group(1).replace('D', 'E'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
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


def count_lines(text):
    """The lines of the solver output that count the cell updates whose stage made nu~, k or omega negative; the solver
    prints that line at the end of a run only when the count is not zero."""
    return [l.strip() for l in text.splitlines() if 'Turbulence variables below zero after a stage' in l]


def pos(case, cfl, vnn, log):
    """The run ends (run.sh), the solver prints no count of negative turbulence values after a stage (no stage made
    one), its step was limited by the turbulence sources at some point (smallest step below half the transport step of
    the initial state), and its first step is the rule with the bound of the stiffness of each cell (first_step), within
    1e-6."""
    text = open(log, errors="replace").read()
    ok = True
    lines = count_lines(text)
    if lines:
        print(f"    [FAIL] negative nu~, k or omega after a stage: {'; '.join(lines)}")
        ok = False
    else:
        print("    negative nu~, k or omega after a stage: none (no count in the solver output)")
    steps = [float(s) for s in re.findall(r'Delta t =\s*([-+0-9.Ee]+)', text)]
    if steps:
        tr = transport_step(case, cfl, vnn)
        ratio = min(steps) / tr
        print(f"    smallest step {min(steps):.4e} s, transport step of the initial state {tr:.4e} s, ratio {ratio:.3f}")
        if not ratio < 0.5:
            print("    [FAIL] the step was never limited by the turbulence sources: the case does not test them")
            ok = False
        want = first_step(case, cfl, vnn)
        err = abs(steps[0] / want - 1.0)
        print(f"    first step {steps[0]:.6e} s, rule with the bound of the stiffness of each cell {want:.6e} s, relative "
              f"difference {err:.1e}")
        if not err < 1.0e-6:
            print("    [FAIL] the first step is not the rule with the bound of the stiffness of the cells")
            ok = False
    else:
        print("    [FAIL] no step in the solver output")
        ok = False
    sys.exit(0 if ok else 1)


def count(log, expect):
    """The count of the cell updates whose stage made nu~, k or omega negative, at the end of a run of the count part:
    with expect = some, one line with a count above zero; with expect = none, no line."""
    lines = count_lines(open(log, errors="replace").read())
    if expect == "some":
        m = re.search(r'(\d+)\s*$', lines[0]) if len(lines) == 1 else None
        ok = m is not None and int(m.group(1)) > 0
        want = "one line, with a count above zero"
    else:
        ok = not lines
        want = "no line"
    print(f"    count of negative nu~, k or omega after a stage: {'; '.join(lines) or 'no line'} (expected {want})")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "order":
        order(*sys.argv[2:7])
    elif what == "pos":
        pos(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    elif what == "count":
        count(sys.argv[2], sys.argv[3])
    elif what in ("order-SA", "order-SST", "pos-SA", "pos-SST"):
        write_case(what)
    elif what == "wilcox":
        write_wilcox()
    elif what == "wilcox-check":
        wilcox_check(sys.argv[2])
    elif what in ("wall-Wilcox2006", "wall-SST"):
        write_wall(what[5:])
    elif what == "wall-check":
        wall_check(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
