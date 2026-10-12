"""Inputs of the MMSResiduals fast test, and the checks of its outputs.

Part 1. A constructed field for the Spalart-Allmaras equation with variable density. Square of side L = 1 mm, periodic in x and in y
(connections 201 on itself), N x N cells, null k faces; air of test/common/Air-viscous (mu and k constant) at rest at 1 bar,
with smooth periodic temperature and nu~:
  T   = T0 [1 + 0.3 sin(2 pi x/L) (1 + 0.4 sin(2 pi y/L))],                       T0 = 300 K, rho = p/(R T);
  nu~ = nu0 [1 + 0.5 cos(2 pi x/L + 0.7) (1 + 0.4 cos(2 pi y/L + 0.3))],          nu0 = 1e-4 m2/s.
There is no wall, so the wall distance is the sentinel 1e8 m and the production and destruction of nu~ vanish. With the fluid at
rest the right-hand side of the conservative equation for rho nu~ at t = 0 is its diffusion part, which in the standard form
(the one of the NASA Turbulence Modeling Resource, Spalart-Allmaras page, conservation form) reads
  (1/sigma) [ div( rho (nu + nu~) grad nu~ ) + rho cb2 |grad nu~|^2 - (nu + nu~) grad rho . grad nu~ ]
    = (rho/sigma) [ div( (nu + nu~) grad nu~ ) + cb2 |grad nu~|^2 ],
nu = mu/rho. The last term on the first line is what makes the conservation form equal to the original equation.

The solver makes one explicit Euler step with local time stepping; the diagnostic output holds the increment of rho nu~ and the
step of each cell, so their ratio is the right-hand side the solver evaluated. It must match the exact one above (derivatives of
the analytic fields by fourth-order differences) to second order: on N = 32 and N = 64 the largest difference over the cells,
divided by the largest value of the right-hand side, must be below 0.02 on N = 64 and fall with an observed order above 1.5.

Part 2. The rotation tensor W_ij = (du_i/dx_j - du_j/dx_i)/2 of a flow is a tensor: SA-RC and SST-RC, whose term r~ uses its
components, must give the same scalars for a flow and for the same flow turned in space. A cube of NR^3 cells of DXR, one block,
extrapolation (400) on all faces; the same air at p = 1 bar and T = 300 K, a velocity u = U0 + G x + H(x, x)/2 about the centre
of the cube, with vorticity along x, y and z and a strain that varies in space, so that the material derivative of the strain,
and with it r~, is not zero; rho nu~ uniform, 3 mu (SA-RC), or k and omega uniform, KR and OMR (SST-RC: in the interior
cells f_r1 lies between its bounds 0 and 1.25 and the production of k below its limiter). One explicit step, time-accurate,
first order:
  a  the field above;
  b  the same field turned by the rotation R that takes x to y, y to z and z to x: u'(x') = R u(R^T x'), on the same cube.
The cell (i, j, k) of b sits where the cell (j, k, i) of a sits, turned, and holds the same numbers, permuted. The change of
rho nu~ (rho k and rho omega with SST-RC) in the step must be the same in the two cells within RTOL_R of its largest value: the
differences come from the order of floating-point sums only.

Part 3. The rotation-curvature correction of Spalart-Allmaras (SA-RC) multiplies the production c_b1 S~ nu~ by
  f_r1 = (1 + c_r1) 2 r*/(1 + r*) [1 - c_r3 atan(c_r2 r~)] - c_r1,   r* = S/omega,
  r~ = 2 omega_ik S_jk [DS_ij/Dt + (eps_imn S_jn + eps_jmn S_in) Omega'_m] / D^4,   D^2 = (S^2 + omega^2)/2,
with omega_ij the rotation tensor of the absolute vorticity and Omega' the rotation rate of the frame (NASA Turbulence Modeling
Resource, SA-RC; c_r1 = 1, c_r2 = 12, c_r3 = 1). A planar block of N3 x N3 cells of DX3 (cubic cells, dyadic sizes),
extrapolation (400) on the i and j faces, null (0) k faces, no wall (wall distance 1e8: no destruction); the same air at
p = 1 bar, T = 300 K, rho nu~ uniform, 3 mu. One explicit step, time-accurate, first order:
  shear    a frame turning at OMEGA3 about z (rigid body), and the uniform shear u = GAMMA3 (y - y_c) in it: DS_ij/Dt = 0, so
           that r~ comes from the rotation of the frame alone. With rho and nu~ uniform and no mass flux across the cells, the
           change of rho nu~ in the step is dt rho c_b1 f_r1 S~ nu~ in every cell, with S~ = |omega_abs| (+ nu~ f_v2/(kappa d)^2,
           1e-20 here); checked within RTOL_3.

Usage (from the test directory):
  python3 -B make_case.py case N                 INPUT/ic.tec and INPUT/bc.txt of part 1
  python3 -B make_case.py error FILE             prints the error of the right-hand side in OUTPUT/diagnostic.tec
  python3 -B make_case.py check E32 E64          exit status 0 when the errors on N = 32 and N = 64 pass
  python3 -B make_case.py rotation a|b [MODEL]   INPUT/ic.tec and INPUT/bc.txt of part 2 (MODEL SA-RC, the default, or SST-RC)
  python3 -B make_case.py rotation-check FA FB   exit status 0 when the field.tec of a and b agree
  python3 -B make_case.py rc shear               INPUT/ic.tec and INPUT/bc.txt of part 3
  python3 -B make_case.py rc-check KIND FIELD    exit status 0 when the field.tec of the case passes
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P0, L = 300.0, 1.0e5, 1.0e-3
SIGMA, CB2 = 2.0 / 3.0, 0.622       # Spalart-Allmaras
NU0 = 1.0e-4
NR, DXR = 8, 1.0e-3                 # part 2: cube of NR^3 cells of DXR
U0R = (10.0, 5.0, 3.0)                                                  # m/s
GR = ((0.0, 300.0, 200.0), (-100.0, 0.0, 400.0), (250.0, -150.0, 0.0))  # du_i/dx_j, 1/s
HR = {(0, 1, 2): 2.0e4, (0, 2, 1): 2.0e4, (1, 0, 0): 1.5e4, (1, 2, 2): -1.0e4,
      (2, 0, 1): -2.5e4, (2, 1, 0): -2.5e4, (2, 1, 1): 1.2e4}          # d2u_i/dx_j dx_k, 1/(m s)
KR, OMR = 1.0, 3.0e3                # part 2 with SST-RC: k (m2/s2) and omega (1/s), uniform
RTOL_R = 1.0e-9
N3, DX3 = 8, 2.0 ** -13             # part 3: N3 x N3 cells of DX3 m, dyadic, so that the cells are cubes to the bit
GAMMA3, OMEGA3 = 1000.0, 200.0      # shear rate, 1/s, and rotation rate of the frame, rad/s
CB1, CR1, CR2, CR3 = 0.1355, 1.0, 12.0, 1.0
RTOL_3 = 1.0e-6


def temp(x, y):
    return T0 * (1.0 + 0.3 * math.sin(2 * math.pi * x / L) * (1.0 + 0.4 * math.sin(2 * math.pi * y / L)))


def rho(x, y):
    return P0 / (RGAS * temp(x, y))


def nut(x, y):
    return NU0 * (1.0 + 0.5 * math.cos(2 * math.pi * x / L + 0.7) * (1.0 + 0.4 * math.cos(2 * math.pi * y / L + 0.3)))


H = L * 1.0e-3                      # step of the fourth-order differences of the analytic fields


def ddx(f, x, y):
    return (-f(x + 2 * H, y) + 8 * f(x + H, y) - 8 * f(x - H, y) + f(x - 2 * H, y)) / (12 * H)


def ddy(f, x, y):
    return (-f(x, y + 2 * H) + 8 * f(x, y + H) - 8 * f(x, y - H) + f(x, y - 2 * H)) / (12 * H)


def rhs(x, y):
    """Exact right-hand side of the standard conservative equation for rho nu~ at t = 0, fluid at rest, no wall."""
    def diff(a, b):
        return MU / rho(a, b) + nut(a, b)

    div = ddx(lambda a, b: diff(a, b) * ddx(nut, a, b), x, y) + ddy(lambda a, b: diff(a, b) * ddy(nut, a, b), x, y)
    g2 = ddx(nut, x, y) ** 2 + ddy(nut, x, y) ** 2
    return rho(x, y) / SIGMA * (div + CB2 * g2)


def write_case(n):
    h = L / n
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    cells = [((i + 0.5) * h, (j + 0.5) * h) for j in range(n) for i in range(n)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={n + 1}, J={n + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(n + 1):
                    for i in range(n + 1):
                        f.write(f" {(i * h, j * h, k * h)[c]:.16E}\n")
        zero = [0.0] * (n * n)
        for col in ([rho(x, y) for x, y in cells], zero, zero, zero, [P0] * (n * n),
                    [rho(x, y) * nut(x, y) for x, y in cells]):
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(b, i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, face, t)

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    out = []
    for j in range(1, n + 1):
        out += [rec(1, 1, j, 1, 1, 201), con(1, n, j, 1, 2, 1, 0, 0, 1)]
    for j in range(1, n + 1):
        out += [rec(1, n, j, 1, 2, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    for i in range(1, n + 1):
        out += [rec(1, i, 1, 1, 3, 201), con(1, i, n, 1, 4, 1, 0, 0, 1)]
    for i in range(1, n + 1):
        out += [rec(1, i, n, 1, 4, 201), con(1, i, 1, 1, 3, 1, 0, 0, 1)]
    for face in (5, 6):
        out += [rec(1, i, j, 1, face, 0) for j in range(1, n + 1) for i in range(1, n + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def write_rotation(name, model="SA-RC"):
    rho0 = P0 / (RGAS * T0)
    names, turb = ('"rhoNuTilde"', (3.0 * MU,)) if model == "SA-RC" else ('"rhoK" "rhoOmega"', (rho0 * KR, rho0 * OMR))
    c = 0.5 * NR * DXR
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    va = {}                                                         # (i, j, k) -> velocity of the field a
    for k in range(NR):
        for j in range(NR):
            for i in range(NR):
                x = ((i + 0.5) * DXR - c, (j + 0.5) * DXR - c, (k + 0.5) * DXR - c)
                u = [U0R[m] + sum(GR[m][n] * x[n] for n in range(3)) for m in range(3)]
                for (m, n, l), h in HR.items():
                    u[m] += 0.5 * h * x[n] * x[l]
                va[i, j, k] = u
    if name == "b":
        # cell (i, j, k) of b takes R u of the cell (j, k, i) of a: (R u)_x = u_z, (R u)_y = u_x, (R u)_z = u_y
        vel = {(i, j, k): (va[j, k, i][2], va[j, k, i][0], va[j, k, i][1]) for (i, j, k) in va}
    else:
        vel = va
    order = [(i, j, k) for k in range(NR) for j in range(NR) for i in range(NR)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(f' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" {names}\n')
        f.write(f" ZONE  T = B1-IG, I={NR + 1}, J={NR + 1}, K={NR + 1}, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-{8 + len(turb)}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(NR + 1):
                for j in range(NR + 1):
                    for i in range(NR + 1):
                        f.write(f" {(DXR * i, DXR * j, DXR * k)[axis]:.17E}\n")
        f.write(f" {rho0:.17E}\n" * NR ** 3)
        for axis in range(3):
            f.write("".join(f" {vel[q][axis]:.17E}\n" for q in order))
        f.write(f" {P0:.17E}\n" * NR ** 3)
        for v in turb:
            f.write(f" {v:.17E}\n" * NR ** 3)
    rows = []
    for face in range(1, 7):
        for b in range(1, NR + 1):
            for a in range(1, NR + 1):
                i, j, k = {1: (1, a, b), 2: (NR, a, b), 3: (a, 1, b), 4: (a, NR, b), 5: (a, b, 1), 6: (a, b, NR)}[face]
                rows.append(f"       1     {i:3d}     {j:3d}     {k:3d}       {face}     400\n")
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(rows))


def rotation_check(fa, fb):
    """The change of rho nu~ (or of rho k and rho omega) in cell (i, j, k) of b against cell (j, k, i) of a, over the largest
    change."""
    rho0 = P0 / (RGAS * T0)
    try:
        za, na = read_zone(fa)
        zb, nb = read_zone(fb)
    except (OSError, ValueError, KeyError) as err:
        print(f"    no field to compare ({err})")
        sys.exit(1)
    if na != NR or nb != NR:
        print(f"    unexpected zones: {na} and {nb} cells per side")
        sys.exit(1)

    ok = True
    for name, start in (("rhoNuTilde", 3.0 * MU), ("rhoK", rho0 * KR), ("rhoOmega", rho0 * OMR)):
        if name not in za:
            continue

        def at(z, i, j, k):
            return z[name][i + NR * (j + NR * k)] - start

        scale = max(abs(at(za, i, j, k)) for k in range(NR) for j in range(NR) for i in range(NR))
        diff = max(abs(at(zb, i, j, k) - at(za, j, k, i)) for k in range(NR) for j in range(NR) for i in range(NR))
        print(f"    change of {name} in the step: largest {scale:.6e}; largest difference between the flow and the "
              f"flow turned {diff:.3e} ({diff / scale if scale > 0 else float('nan'):.3e} of it)")
        ok = ok and scale > 0.0 and diff <= RTOL_R * scale
    sys.exit(0 if ok else 1)


def write_rc(kind):
    rho0 = P0 / (RGAS * T0)
    yc = 0.5 * N3 * DX3
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    u = [GAMMA3 * ((j + 0.5) * DX3 - yc) for j in range(N3) for _ in range(N3)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={N3 + 1}, J={N3 + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(N3 + 1):
                    for i in range(N3 + 1):
                        f.write(f" {(DX3 * i, DX3 * j, DX3 * k)[axis]:.17E}\n")
        f.write(f" {rho0:.17E}\n" * N3 * N3)
        f.write("".join(f" {v:.17E}\n" for v in u))
        f.write(f" {0.0:.17E}\n" * N3 * N3 * 2)
        f.write(f" {P0:.17E}\n" * N3 * N3)
        f.write(f" {3.0 * MU:.17E}\n" * N3 * N3)
    rows = []
    for face, idx in ((1, lambda m: (1, m)), (2, lambda m: (N3, m)), (3, lambda m: (m, 1)), (4, lambda m: (m, N3))):
        for m in range(1, N3 + 1):
            i, j = idx(m)
            rows.append(f"       1     {i:3d}     {j:3d}       1       {face}     400\n")
    for face in (5, 6):
        for j in range(1, N3 + 1):
            for i in range(1, N3 + 1):
                rows.append(f"       1     {i:3d}     {j:3d}       1       {face}       0\n")
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(rows))


def f_r1_shear():
    """f_r1 of the uniform shear u = GAMMA3 y in a frame turning at OMEGA3 about z (NASA TMR, SA-RC), |omega_abs| and r~."""
    s = 0.5 * GAMMA3                                 # S_12 = S_21
    w = 0.5 * (GAMMA3 - 2.0 * OMEGA3)                # omega_12 of the absolute vorticity: -(omega_3 + 2 Omega)/2
    S = GAMMA3                                       # sqrt(2 S_ij S_ij)
    O = abs(GAMMA3 - 2.0 * OMEGA3)                   # |omega_abs|
    # (eps_imn S_jn + eps_jmn S_in) Omega_m = diag(-OMEGA3 GAMMA3, OMEGA3 GAMMA3, 0); omega_ik S_jk = diag(w s, -w s, 0)
    num = 2.0 * ((-OMEGA3 * GAMMA3) * (w * s) + (OMEGA3 * GAMMA3) * (-w * s))
    D2 = 0.5 * (S ** 2 + O ** 2)
    rt = num / D2 ** 2
    return (1.0 + CR1) * (2.0 * S / (S + O)) * (1.0 - CR3 * math.atan(CR2 * rt)) - CR1, O, rt


def rc_check(kind, path):
    try:
        z, n = read_zone(path)
        with open(path) as f:
            f.readline()
            t = float(re.search(r"SOLUTIONTIME\s*=\s*([-+0-9.eEdD]+)", f.readline()).group(1).replace("D", "E"))
    except (OSError, ValueError, KeyError, AttributeError) as err:
        print(f"    no field of the case {kind} ({err})")
        sys.exit(1)
    rho0 = P0 / (RGAS * T0)
    start = 3.0 * MU
    new = z["rhoNuTilde"]
    if n != N3 or any(math.isnan(v) or math.isinf(v) for v in new):
        print(f"    rho nu~ is not a number after the step, or the zone is not the case ({n} cells per side)")
        sys.exit(1)
    if kind == "shear":
        fr1, om, rt = f_r1_shear()
        want = t * rho0 * CB1 * fr1 * om * (start / rho0)
        got = [v - start for v in new]
        err = max(abs(g - want) for g in got)
        print(f"    dt = {t:.6e} s, r~ = {rt:.6e}, f_r1 = {fr1:.6f}; change of rho nu~ expected {want:.6e}, "
              f"computed {min(got):.6e} to {max(got):.6e} kg/(m s), largest difference {err / abs(want):.3e} of it")
        sys.exit(0 if err <= RTOL_3 * abs(want) else 1)
    sys.exit(f"unknown case {kind}")


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the number of cells per side."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, ni - 1


def error(path):
    """Largest difference between the right-hand side of rho nu~ evaluated by the solver (increment over step) and the exact
    one, over the largest exact value."""
    d, n = read_zone(path)
    h = L / n
    cells = [((i + 0.5) * h, (j + 0.5) * h) for j in range(n) for i in range(n)]
    got = [r / t for r, t in zip(d["rhoNuTilde"], d["dt"])]
    exact = [rhs(x, y) for x, y in cells]
    e = max(abs(g - x) for g, x in zip(got, exact)) / max(abs(x) for x in exact)
    print(f"{e:.6e}")


def check(e32, e64):
    p = math.log2(e32 / e64) if e32 > 0 and e64 > 0 else float("nan")
    print(f"    error of the right-hand side of rho nu~: {e32:.3e} (N = 32), {e64:.3e} (N = 64), observed order {p:.2f}")
    sys.exit(0 if e64 < 0.02 and p > 1.5 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(int(sys.argv[2]))
    elif what == "error":
        error(sys.argv[2])
    elif what == "check":
        check(float(sys.argv[2]), float(sys.argv[3]))
    elif what == "rotation":
        write_rotation(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "SA-RC")
    elif what == "rotation-check":
        rotation_check(sys.argv[2], sys.argv[3])
    elif what == "rc":
        write_rc(sys.argv[2])
    elif what == "rc-check":
        rc_check(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
