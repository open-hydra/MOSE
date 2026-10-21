"""Part 4 of MMSResiduals (run.sh): the viscous stress on the 2-D axisymmetric wedge, inputs and check.

Two velocity fields on the 2-D axisymmetric wedge of MOSE: one cell in k, k faces of type 200 at theta = -delta/2 and
+delta/2 (delta = 1 degree), the axis at j = 0. Cylindrical coordinates (x, r, theta); lambda = -2 mu/3.

homann  Homann's axisymmetric stagnation-point flow, an exact solution of the Navier-Stokes equations (F. Homann, ZAMM 16,
        1936, 153-164, not read; the equation as in P. Rajamanickam, arXiv:2508.17982v4, Eq. (12)): with the strain rate a,
        nu = mu/rho, L = sqrt(nu/a) and eta = x/L,
          u_x = -2 sqrt(a nu) f(eta),   u_r = a r f'(eta),   f''' + 2 f f'' - f'^2 + 1 = 0,   f(0) = f'(0) = 0,   f'(inf) = 1,
        with f''(0) = 1.3119377 (by shooting; 1.311938 in H. Guo et al., arXiv:2010.05682v1, Table 8). The velocity is
        divergence free, so with constant
        mu the viscous force per unit volume is mu lap(u):
          radial    rho a^2 r f'''(eta),
          axial     -2 rho a sqrt(a nu) f''(eta),
          azimuthal 0.
        In the radial momentum it is (1/r) d(r tau_rr)/dr + d(tau_xr)/dx - tau_thetatheta/r, with tau_thetatheta =
        2 mu u_r/r + lambda div(u) the normal viscous stress on the two faces of the wedge; here tau_rr = tau_thetatheta =
        2 mu a f'(eta), so tau_rr/r and tau_thetatheta/r cancel and each of them is far larger than the force near the axis.

swirl   A polynomial field with dilatation and swirl, the two things Homann's flow does not have (there div u = 0 and
        tau_rr = tau_thetatheta), with L = 1 mm:
          u_x = b x^2,   u_r = c r^3,   u_theta = s r^3 (1 + x/L),   div u = 2 b x + 4 c r^2.
        With constant mu the stresses are
          tau_xx = 4 mu b x + lambda div u,   tau_rr = 6 mu c r^2 + lambda div u,   tau_thetatheta = 2 mu c r^2 + lambda div u,
          tau_xr = 0,   tau_rtheta = mu (du_theta/dr - u_theta/r) = 2 mu s r^2 (1 + x/L),   tau_xtheta = mu s r^3/L,
        and the viscous force per unit volume, div(tau), is
          axial      d tau_xx/dx                                        = (4 mu + 2 lambda) b       = (8/3) mu b,
          radial     d tau_rr/dr + (tau_rr - tau_thetatheta)/r           = (12 mu + 8 lambda + 4 mu) c r = (32/3) mu c r,
          azimuthal  d tau_rtheta/dr + 2 tau_rtheta/r + d tau_xtheta/dx  = 8 mu s r (1 + x/L).
        The field is not a solution of the Navier-Stokes equations; only the viscous right-hand side is checked. The faces of
        the wedge carry -tau_thetatheta/r of the radial force, here with its lambda div(u) part, and tau_rtheta/r of the
        azimuthal one (a quarter of it); their tau_xtheta cancel.

Case (python3 -B wedge_viscous.py case FIELD N): x in [0, 3 mm], r in [0, 3 mm] with the axis nodes at r = 1e-6 dr (like the
wedge meshes of MOSE, so that no face has zero area), N x N square cells; air of test/common/Air-viscous (mu and k constant)
at uniform density and pressure (rho = 1.1858686e-3 kg/m3, T = 300 K, so nu = 1e-2 m2/s; for homann a = 1e4 1/s, L = 1 mm),
with the velocity of FIELD at the cell centres. Faces: 400 on i = 1, i = N and j = N; 300 on the axis (j = 1); 200 on k.

Observable: the solver makes one explicit Euler step with local time stepping with the Navier-Stokes equations and one with
the Euler equations, from the same state. The diagnostic output holds the increment of each conserved variable and the step
of each cell, so the difference of the two ratios increment/step is the viscous right-hand side the solver evaluates.

Checks, over the cells with 2 <= i <= N-1 and 1 <= j <= N-1 (the i = 1, i = N and j = N cells have an outer face with no
viscous flux), with errors over the largest exact force of the same component:
  homann  radial  error below 0.02 on N = 64; observed order above 1.5 between N = 32 and N = 64 on the cells with
                  r >= L/2; on the first row of cells off the axis, the error over the hoop term 2 mu a f'/r below 1e-4 on
                  N = 64. Near the axis the azimuthal derivatives of the faces i and j, taken across the rotated ghost cells,
                  carry an error of order delta^2 that the mesh does not reduce (delta^2/12 of the hoop term, 2.5e-5 for
                  delta = 1 degree), hence the separate check there;
          axial   error below 0.02 on N = 64 and observed order above 1.5;
          azimuthal  below 1e-9 of the largest radial force.
  swirl   radial and azimuthal  error below 0.02 on N = 64 and observed order above 1.5 on the cells with r >= L/2 (next to
                  the axis the error of these cubic fields is of first order in the mesh);
          axial   below 1e-9 (exact for a quadratic u_x; the tau_xtheta of the two faces cancel).

Usage (from the test directory), FIELD = homann or swirl:
  python3 -B wedge_viscous.py case FIELD N                    INPUT/ic.tec and INPUT/bc.txt
  python3 -B wedge_viscous.py error FIELD NS_FILE EU_FILE      prints the errors of one mesh, comma-separated; homann: radial,
                                                           radial on r >= L/2, radial over the hoop term on the first row,
                                                           axial, azimuthal; swirl: radial, radial on r >= L/2, azimuthal,
                                                           azimuthal on r >= L/2, axial
  python3 -B wedge_viscous.py check FIELD E32 E64              exit status 0 when the errors on N = 32 and N = 64 pass
                                                           (E32, E64: the five numbers printed by 'error')
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
RHO0, T0 = 1.1858686e-3, 300.0      # nu = MU/RHO0 = 1e-2 m2/s
A = 1.0e4                           # strain rate of homann, 1/s
NU = MU / RHO0
L = math.sqrt(NU / A)               # 1 mm
P0 = RHO0 * RGAS * T0
DELTA = math.radians(1.0)           # wedge angle
SIZE = 3.0 * L                      # x and r extent
FPP0 = 1.311937693880               # f''(0) of Homann's equation, by shooting (RK4, eta_max 6 to 10)
B_SW = 1.1e6                        # swirl: u_x = b x^2, 9.9 m/s at x = 3 mm
C_SW = 3.7e8                        # swirl: u_r = c r^3, 10 m/s at r = 3 mm
S_SW = 1.85e8                       # swirl: u_theta = s r^3 (1 + x/L), 5 m/s at r = 3 mm, x = 0
FIELDS = ("homann", "swirl")


def homann(n):
    """f, f', f'', f''' at the cell centres eta = (i + 1/2) 3/n, i = 0 .. n-1 (RK4 from the wall)."""
    deta = SIZE / L / n
    sub = 40                        # RK4 steps per half cell
    h = 0.5 * deta / sub

    def rhs(y):
        f, fp, fpp = y
        return (fp, fpp, -2.0 * f * fpp + fp * fp - 1.0)

    y = (0.0, 0.0, FPP0)
    out = []
    for step in range(1, 2 * n * sub + 1):
        k1 = rhs(y)
        k2 = rhs(tuple(a + 0.5 * h * b for a, b in zip(y, k1)))
        k3 = rhs(tuple(a + 0.5 * h * b for a, b in zip(y, k2)))
        k4 = rhs(tuple(a + h * b for a, b in zip(y, k3)))
        y = tuple(a + h / 6.0 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4))
        if step % (2 * sub) == sub:
            out.append(y + (-2.0 * y[0] * y[2] + y[1] ** 2 - 1.0,))
    return out


def radii(n):
    """Node radii: the axis row at 1e-6 dr, then j dr."""
    dr = SIZE / n
    return [1.0e-6 * dr] + [j * dr for j in range(1, n + 1)]


def swirl_velocity(xc, rc):
    """u_x, u_r, u_theta of the swirl field."""
    return B_SW * xc * xc, C_SW * rc ** 3, S_SW * rc ** 3 * (1.0 + xc / L)


def swirl_force(xc, rc):
    """Exact viscous force per unit volume of the swirl field: axial, radial, azimuthal."""
    return (8.0 / 3.0) * MU * B_SW, (32.0 / 3.0) * MU * C_SW * rc, 8.0 * MU * S_SW * rc * (1.0 + xc / L)


def write_case(field, n):
    dx = SIZE / n
    rn = radii(n)
    f = homann(n) if field == "homann" else None
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    u, v, w = [], [], []
    for j in range(n):
        rc = 0.5 * (rn[j] + rn[j + 1])
        for i in range(n):
            if field == "homann":
                u.append(-2.0 * math.sqrt(A * NU) * f[i][0])
                v.append(A * rc * f[i][1])
                w.append(0.0)
            else:
                ux, ur, ut = swirl_velocity((i + 0.5) * dx, rc)
                u.append(ux)
                v.append(ur)
                w.append(ut)
    with open("INPUT/ic.tec", "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        out.write(f" ZONE  T = B1-IG, I={n + 1}, J={n + 1}, K=2, DATAPACKING=BLOCK, "
                  "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for comp in range(3):
            for k in range(2):
                for j in range(n + 1):
                    for i in range(n + 1):
                        val = (i * dx, rn[j] * c, (-1.0, 1.0)[k] * rn[j] * s)[comp]
                        out.write(f" {val:.16E}\n")
        for col in ([RHO0] * (n * n), u, v, w, [P0] * (n * n)):
            out.write("".join(f" {x:.16E}\n" for x in col))

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    rows = []
    rows += [rec(1, j, 1, 400) for j in range(1, n + 1)]
    rows += [rec(n, j, 2, 400) for j in range(1, n + 1)]
    rows += [rec(i, 1, 3, 300) for i in range(1, n + 1)]
    rows += [rec(i, n, 4, 400) for i in range(1, n + 1)]
    for face in (5, 6):
        rows += [rec(i, j, face, 200) for j in range(1, n + 1) for i in range(1, n + 1)]
    with open("INPUT/bc.txt", "w") as out:
        out.write("".join(rows))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not, and the number of cells per side."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out, ni - 1


def viscous_rhs(ns, eu, q):
    """Viscous right-hand side of cell q (Navier-Stokes minus Euler): axial, radial, azimuthal."""
    return tuple(ns[v][q] / ns["dt"][q] - eu[v][q] / eu["dt"][q] for v in ("rhou", "rhov", "rhow"))


def errors_homann(ns_path, eu_path):
    """Errors of the viscous right-hand side of one mesh: radial, radial on r >= L/2 (both over the largest exact radial
    force), radial on the first row over the hoop term 2 mu a f'/r, axial (over the largest exact axial force), azimuthal
    (over the largest exact radial force)."""
    ns, n = read_zone(ns_path)
    eu, n2 = read_zone(eu_path)
    assert n == n2
    rn = radii(n)
    f = homann(n)
    er = er_far = er_hoop = ex = ez = 0.0
    fr_max = fx_max = 0.0
    for j in range(n):
        rc = 0.5 * (rn[j] + rn[j + 1])
        for i in range(n):
            q = j * n + i
            fr_ex = RHO0 * A * A * rc * f[i][3]
            fx_ex = -2.0 * RHO0 * A * math.sqrt(A * NU) * f[i][2]
            fr_max = max(fr_max, abs(fr_ex))
            fx_max = max(fx_max, abs(fx_ex))
            if not (1 <= i <= n - 2 and j <= n - 2):
                continue
            fx, fr, fz = viscous_rhs(ns, eu, q)
            er = max(er, abs(fr - fr_ex))
            if rc >= 0.5 * L:
                er_far = max(er_far, abs(fr - fr_ex))
            if j == 0:
                er_hoop = max(er_hoop, abs(fr - fr_ex) / (2.0 * MU * A * f[i][1] / rc))
            ex = max(ex, abs(fx - fx_ex))
            ez = max(ez, abs(fz))
    return er / fr_max, er_far / fr_max, er_hoop, ex / fx_max, ez / fr_max


def errors_swirl(ns_path, eu_path):
    """Errors of the viscous right-hand side of one mesh, over the largest exact force of the same component: radial, radial
    on r >= L/2, azimuthal, azimuthal on r >= L/2, axial."""
    ns, n = read_zone(ns_path)
    eu, n2 = read_zone(eu_path)
    assert n == n2
    rn = radii(n)
    dx = SIZE / n
    fmax = [0.0, 0.0, 0.0]
    err = [0.0, 0.0, 0.0]
    err_far = [0.0, 0.0, 0.0]
    for j in range(n):
        rc = 0.5 * (rn[j] + rn[j + 1])
        for i in range(n):
            xc = (i + 0.5) * dx
            fex = swirl_force(xc, rc)
            for k in range(3):
                fmax[k] = max(fmax[k], abs(fex[k]))
            if not (1 <= i <= n - 2 and j <= n - 2):
                continue
            fnum = viscous_rhs(ns, eu, j * n + i)
            for k in range(3):
                e = abs(fnum[k] - fex[k])
                err[k] = max(err[k], e)
                if rc >= 0.5 * L:
                    err_far[k] = max(err_far[k], e)
    return err[1] / fmax[1], err_far[1] / fmax[1], err[2] / fmax[2], err_far[2] / fmax[2], err[0] / fmax[0]


def order(e32, e64):
    return math.log2(e32 / e64) if e32 > 0 and e64 > 0 else float("nan")


def check_homann(a32, a64):
    ok = True
    good = a64[0] < 0.02
    ok = ok and good
    print(f"    radial viscous force: error {a32[0]:.3e} (N = 32), {a64[0]:.3e} (N = 64)  [{'ok' if good else 'FAIL'}]")
    good = order(a32[1], a64[1]) > 1.5
    ok = ok and good
    print(f"      on r >= L/2: {a32[1]:.3e} (N = 32), {a64[1]:.3e} (N = 64), observed order {order(a32[1], a64[1]):.2f}"
          f"  [{'ok' if good else 'FAIL'}]")
    good = a64[2] < 1.0e-4
    ok = ok and good
    print(f"      first row off the axis, over the hoop term: {a32[2]:.3e} (N = 32), {a64[2]:.3e} (N = 64)"
          f"  [{'ok' if good else 'FAIL'}]")
    good = a64[3] < 0.02 and order(a32[3], a64[3]) > 1.5
    ok = ok and good
    print(f"    axial viscous force: error {a32[3]:.3e} (N = 32), {a64[3]:.3e} (N = 64), observed order "
          f"{order(a32[3], a64[3]):.2f}  [{'ok' if good else 'FAIL'}]")
    good = max(a32[4], a64[4]) < 1.0e-9
    ok = ok and good
    print(f"    azimuthal viscous force: {a32[4]:.1e} (N = 32), {a64[4]:.1e} (N = 64)  [{'ok' if good else 'FAIL'}]")
    return ok


def check_swirl(a32, a64):
    ok = True
    for name, k in (("radial", 0), ("azimuthal", 2)):
        good = a64[k] < 0.02 and order(a32[k + 1], a64[k + 1]) > 1.5
        ok = ok and good
        print(f"    {name} viscous force: error {a32[k]:.3e} (N = 32), {a64[k]:.3e} (N = 64); on r >= L/2 {a32[k + 1]:.3e},"
              f" {a64[k + 1]:.3e}, observed order {order(a32[k + 1], a64[k + 1]):.2f}  [{'ok' if good else 'FAIL'}]")
    good = max(a32[4], a64[4]) < 1.0e-9
    ok = ok and good
    print(f"    axial viscous force: error {a32[4]:.1e} (N = 32), {a64[4]:.1e} (N = 64)  [{'ok' if good else 'FAIL'}]")
    return ok


def main():
    what = sys.argv[1]
    field = sys.argv[2] if len(sys.argv) > 2 else ""
    if field not in FIELDS:
        sys.exit(f"unknown field {field}: use one of {', '.join(FIELDS)}")
    if what == "case":
        write_case(field, int(sys.argv[3]))
    elif what == "error":
        errs = (errors_homann if field == "homann" else errors_swirl)(sys.argv[3], sys.argv[4])
        print(",".join(f"{e:.6e}" for e in errs))
    elif what == "check":
        a32 = [float(v) for v in sys.argv[3].split(",")]
        a64 = [float(v) for v in sys.argv[4].split(",")]
        ok = (check_homann if field == "homann" else check_swirl)(a32, a64)
        sys.exit(0 if ok else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
