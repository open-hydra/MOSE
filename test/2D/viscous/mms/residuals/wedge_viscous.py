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
wedge meshes of MOSE, so that no face has zero area), N x N square cells; with the option axis0, the axis nodes at r = 0,
so that the faces on the axis have zero area; with the option geo, radii growing geometrically from the axis along one
smooth map for every N, r_j = 3 mm (exp(beta j/N) - 1)/(exp(beta) - 1), beta = 32 ln 1.014 (cells growing by 1.014 on
N = 32 and by 1.007 on N = 64, from 0.79 to 1.24 times the uniform ones), the axis nodes at 1e-6 of the first cell; air of
test/common/Air-viscous (mu and k constant)
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
The same two fields on the radii of the option geo (stretched along r, as the wedges of nozzles): the same checks, and every
error above that is not round-off (homann: radial, radial on r >= L/2, first row; swirl: radial and azimuthal, also on
r >= L/2) within 2 times the same error of the uniform radii (printed next to it). The axial force of homann is printed
only: it has no part from the faces of type 200 (tau_xtheta = 0 without swirl), and next to the axis the faces j of
stretched radii take the derivative of u_r along x at the mean radius of the two rows, not at the face, an error of first
order there, with stencil as with inverse-mean (1.7e-3 and 1.2e-3 of the largest axial force on N = 64, 6.3e-5 on uniform
radii). The faces of type 200 take the derivatives
along r from the differences of the centres of the row of the cell and of its image row (diffusive-metric = stencil): before,
the mean edges of the cell and of its image, which gave (q - 1)/4 of the hoop term on the first row, 1.7e-3 here on N = 64.
Linear fields (python3 -B wedge_viscous.py radial N LAW, axial N), whose viscous force is zero:
  radial  u_r = U r (U = 5e3 1/s), N = 16 rectangular cells, radii uniform, or growing by 1.1 or 1.2 (LAW uniform, geo1.1,
          geo1.2), 400 on i = 1, i = N and j = N: one Navier-Stokes and one Euler step. The faces of type 200 are exact for
          this field (u(R x) = R u(x)); the faces i and j take the derivative in azimuth across the rotated ghost cells over
          the mean edge of the cells, cos^2(delta/2) times too long (docs/theory/numerics.md), which leaves in every cell a
          radial force of sin^2(delta/2) of the hoop term (2/3) mu U/r_c (r_c the radius of the cell centre), whatever the
          radii, and on the first row also the round-off of the mirror ghost cell of the axis at 1e-6 dr. Checked, on the cells
          with no face 400: |force/hoop term - sin^2(delta/2)| below tol on the rows j >= 2 and below 3e-6 on the first row,
          tol = max(1e-9, 10 times the largest deviation on the uniform radii, measured in the same run and printed, itself
          at most 1e-11, so that the tolerance of the stretched radii cannot float). (Before,
          with the inverse-mean metric of the faces of type 200, 8.4e-3 on the rows j >= 2 and 5.0e-2 on the first row at 1.2.)
  axial   u = (a x, 0, 0) (a = 1e3 1/s), N = 16, cells growing by 1.15 along x from a wall 301 on i = 1 (x = 0, adiabatic) and
          along r, lines j bent inside the box, so that the faces of type 200 of the cells next to the wall take the ghost
          cells of the wall in their differences along the face: a Navier-Stokes step against one with mu = k = 1e-30 (the
          solver takes no wall with the Euler equations). The field is zero on the wall and holds no velocity in azimuth, so
          the derivatives in azimuth carry no error: in every cell with no face 400, the wall and axis rows included, the
          viscous force must vanish below 1e-11 of mu |a| (sum of the face areas)/volume.
Axis at r = 0 on the wedge (homann, N = 32, option axis0): the faces on the axis have zero area and no normal, and their
ghost cells copy the cells. With diffusive-metric = inverse-mean or mean the faces of the wedge (BC_Rotational_Periodic_Visc)
take the derivative in r across the axis with the ghost cell one cell away: there they would see half of d u_r/dr, and the
hoop stress of the first row would be wrong (lambda div(u) in tau_thetatheta by mu a f'/3, 1/6 of the hoop term 2 mu a f'/r).
With stencil the copy sits at the cell and the faces take the difference on one side of it. A viscous run stops at setup
with the message on the axis at r = 0, with diffusive-metric = stencil and inverse-mean, without NaN before it.
Axis at r = 0 on a 3-D sector (python3 -B wedge_viscous.py sector N [axis0]): Homann's flow on the same box in x and r, N x N
cells, and SECTOR_K cells of DELTA around the axis, from -SECTOR_K DELTA/2 to +SECTOR_K DELTA/2, closed by type 200 on the k
faces (u_r turned to the azimuth of each cell), 300 on the axis, 400 on i = 1, i = N and j = N. With axis0 the faces on the axis
have zero area; with inverse-mean the faces of the first row see half of d u_r/dr, as on the wedge (1/6 of the hoop term):
a viscous run stops at setup with the message on the axis at r = 0 of the sector, with diffusive-metric = stencil and
inverse-mean, without NaN before it.
Polar line of a 3-D mesh (python3 -B wedge_viscous.py polar AXIS): the same box in x and r, 360 degrees around the axis in 8
cells closed by a face between the block and itself (101), 5 x 4 cells in x and r, the axis nodes at r = 0, a linear
velocity U + G x with every entry of G (G_POLAR) at the cell centres, faces 400 on i = 1, i = N and j = N and of type
AXIS on the axis (j = 1):
  300  symmetry: the faces on the axis have zero area, and their ghost cells copy the cells, which the stencil metric
       places at the cell: the viscous force of the cells with no face 400, the first row next to the axis included, must
       vanish to round-off, below 1e-11 of mu |G| (sum of the face areas)/volume, as for any linear field;
  300  with diffusive-metric = inverse-mean or mean, whose viscous force of the first row grows as the mesh is refined:
       the viscous run must stop at setup with the message that names the metric and stencil, without NaN before it;
  302  wall at 300 K: a wall on faces of zero area, whose gradient is not defined. The run must stop at setup with the
       message that names the block, the cell and the face and says to use the symmetry condition (300) on the axis,
       without NaN before it.

Usage (from the test directory), FIELD = homann or swirl:
  python3 -B wedge_viscous.py case FIELD N [axis0|geo]        INPUT/ic.tec and INPUT/bc.txt (axis0: axis nodes at r = 0; geo:
                                                           the geometric radii)
  python3 -B wedge_viscous.py error FIELD NS_FILE EU_FILE [geo]  prints the errors of one mesh, comma-separated; homann: radial,
                                                           radial on r >= L/2, radial over the hoop term on the first row,
                                                           axial, azimuthal; swirl: radial, radial on r >= L/2, azimuthal,
                                                           azimuthal on r >= L/2, axial
  python3 -B wedge_viscous.py check FIELD E32 E64              exit status 0 when the errors on N = 32 and N = 64 pass
                                                           (E32, E64: the five numbers printed by 'error')
  python3 -B wedge_viscous.py checkgeo FIELD E32 E64 G32 G64   exit status 0 when the errors G32, G64 of the geometric radii pass
                                                           the same checks and stay within 2 times those of the uniform radii
  python3 -B wedge_viscous.py radial N LAW                     INPUT/ic.tec and INPUT/bc.txt of the radial linear field
  python3 -B wedge_viscous.py radialcheck N LAW NS EU [TOL]    prints the largest deviations; with TOL, exit status 0 when they pass
  python3 -B wedge_viscous.py axial N                          INPUT/ic.tec, INPUT/bc.txt and INPUT/transport-ref.dat of the axial
                                                           linear field with the wall
  python3 -B wedge_viscous.py axialcheck N NS REF              exit status 0 when the viscous force of the axial field vanishes
  python3 -B wedge_viscous.py axisstop LOG                     exit status 0 when the run of LOG stopped at setup with the message
                                                           on the axis at r = 0 of the wedge, without NaN or trap
  python3 -B wedge_viscous.py sector N [axis0]                 INPUT/ic.tec and INPUT/bc.txt of Homann's flow on the 3-D sector
  python3 -B wedge_viscous.py sectorstop LOG                   exit status 0 when the run of LOG stopped at setup with the message
                                                           on the axis at r = 0 of the sector, without NaN or trap
  python3 -B wedge_viscous.py polar AXIS                       INPUT/ic.tec and INPUT/bc.txt of the polar mesh (AXIS: 300 or 302)
  python3 -B wedge_viscous.py linear NS_FILE EU_FILE           exit status 0 when the viscous force of the linear field on the
                                                           polar mesh vanishes (300 on the axis)
  python3 -B wedge_viscous.py polarstop LOG METRIC            exit status 0 when the run of LOG (300 on the axis, METRIC inverse-mean
                                                           or mean) stopped at setup with the message that names stencil
  python3 -B wedge_viscous.py stopped LOG                      exit status 0 when the run of LOG stopped at setup with the message
                                                           on the wall faces of zero area on the axis, without NaN or trap
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


GEO_BETA = 32.0 * math.log(1.014)   # geo: r_j = SIZE (exp(beta j/N) - 1)/(exp(beta) - 1), growth 1.014 on N = 32


def radii(n, axis0=False, geo=False):
    """Node radii: the axis row at 1e-6 dr (at 0 with axis0), then j dr; with geo the radii of one smooth geometric map
    for every N (GEO_BETA), the axis row at 1e-6 of the first cell."""
    if geo:
        rs = [SIZE * math.expm1(GEO_BETA * j / n) / math.expm1(GEO_BETA) for j in range(n + 1)]
        rs[0] = 1.0e-6 * rs[1]
        return rs
    dr = SIZE / n
    return [0.0 if axis0 else 1.0e-6 * dr] + [j * dr for j in range(1, n + 1)]


def swirl_velocity(xc, rc):
    """u_x, u_r, u_theta of the swirl field."""
    return B_SW * xc * xc, C_SW * rc ** 3, S_SW * rc ** 3 * (1.0 + xc / L)


def swirl_force(xc, rc):
    """Exact viscous force per unit volume of the swirl field: axial, radial, azimuthal."""
    return (8.0 / 3.0) * MU * B_SW, (32.0 / 3.0) * MU * C_SW * rc, 8.0 * MU * S_SW * rc * (1.0 + xc / L)


def write_case(field, n, axis0=False, geo=False):
    dx = SIZE / n
    rn = radii(n, axis0, geo)
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


SECTOR_K = 4                        # 3-D sector: cells of DELTA around the axis


def write_sector(n, kk, axis0=False):
    """Homann's flow on a 3-D sector of kk cells of DELTA around the axis, type 200 on the k faces, N x N cells in x and r."""
    dx = SIZE / n
    rn = radii(n, axis0)
    f = homann(n)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)

    def theta(k):
        return (k - 0.5 * kk) * DELTA

    u, v, w = [], [], []
    for k in range(kk):
        c, s = math.cos(theta(k + 0.5)), math.sin(theta(k + 0.5))
        for j in range(n):
            ur = A * 0.5 * (rn[j] + rn[j + 1])
            for i in range(n):
                u.append(-2.0 * math.sqrt(A * NU) * f[i][0])
                v.append(ur * f[i][1] * c)
                w.append(ur * f[i][1] * s)
    nc = n * n * kk
    with open("INPUT/ic.tec", "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        out.write(f" ZONE  T = B1-IG, I={n + 1}, J={n + 1}, K={kk + 1}, DATAPACKING=BLOCK, "
                  "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for comp in range(3):
            for k in range(kk + 1):
                c, s = math.cos(theta(k)), math.sin(theta(k))
                for j in range(n + 1):
                    for i in range(n + 1):
                        out.write(f" {(i * dx, rn[j] * c, rn[j] * s)[comp]:.16E}\n")
        for col in ([RHO0] * nc, u, v, w, [P0] * nc):
            out.write("".join(f" {x:.16E}\n" for x in col))

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    rows = []
    for k in range(1, kk + 1):
        rows += [rec(1, j, k, 1, 400) for j in range(1, n + 1)]
        rows += [rec(n, j, k, 2, 400) for j in range(1, n + 1)]
        rows += [rec(i, 1, k, 3, 300) for i in range(1, n + 1)]
        rows += [rec(i, n, k, 4, 400) for i in range(1, n + 1)]
    for face, k in ((5, 1), (6, kk)):
        rows += [rec(i, j, k, face, 200) for j in range(1, n + 1) for i in range(1, n + 1)]
    with open("INPUT/bc.txt", "w") as out:
        out.write("".join(rows))


POLAR_N = (5, 4, 8)                 # polar mesh: cells along x, r and around the axis
U_POLAR = (10.0, -3.0, 2.0)          # polar mesh: u = U_POLAR + G_POLAR x, a linear field with every gradient entry, m/s
G_POLAR = ((1.0e3, 1.2e3, -0.7e3), (0.5e3, 2.0e3, 0.9e3), (-1.1e3, 0.4e3, -1.5e3))   # 1/s


def polar_nodes():
    """Nodes of the polar mesh, X[k][j][i] = (x, y, z): the axis row at r = 0, the last plane around the axis the first one."""
    ni, nj, nk = POLAR_N
    dx, dr = SIZE / ni, SIZE / nj
    nodes = []
    for k in range(nk + 1):
        theta = 2.0 * math.pi * (k % nk) / nk
        nodes.append([[(i * dx, j * dr * math.cos(theta), j * dr * math.sin(theta)) for i in range(ni + 1)]
                      for j in range(nj + 1)])
    return nodes


def polar_cell(nodes, i, j, k):
    """The eight nodes of cell (i, j, k) (0-based), in the order of their offsets (di, dj, dk)."""
    return [nodes[k + c][j + b][i + a] for c in (0, 1) for b in (0, 1) for a in (0, 1)]


def polar_velocity(xc):
    return tuple(U_POLAR[m] + sum(G_POLAR[m][n] * xc[n] for n in range(3)) for m in range(3))


def write_polar(axis_type):
    """The 3-D mesh around the axis (see the module docstring): the linear velocity U_POLAR + G_POLAR x at the cell
    centres (the mean of the eight nodes, as MOSE takes them), axis faces of type axis_type."""
    ni, nj, nk = POLAR_N
    nodes = polar_nodes()
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    vel = []
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                cell = polar_cell(nodes, i, j, k)
                vel.append(polar_velocity([sum(p[n] for p in cell) / 8.0 for n in range(3)]))
    nc = ni * nj * nk
    with open("INPUT/ic.tec", "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        out.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                  "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for comp in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        out.write(f" {nodes[k][j][i][comp]:.16E}\n")
        for col in ([RHO0] * nc, [v[0] for v in vel], [v[1] for v in vel], [v[2] for v in vel], [P0] * nc):
            out.write("".join(f" {x:.16E}\n" for x in col))

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    def src(i, j, k, face):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, 1, 0, 0, 1)

    rows = []
    for k in range(1, nk + 1):
        for j in range(1, nj + 1):
            rows += [rec(1, j, k, 1, 400), rec(ni, j, k, 2, 400)]
    for k in range(1, nk + 1):
        for i in range(1, ni + 1):
            rows.append(rec(i, 1, k, 3, axis_type))
            if axis_type == 302:
                rows.append(f"    {T0:.6E}, 0.000000E+00, 0.000000E+00,\n")
            rows.append(rec(i, nj, k, 4, 400))
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            rows += [rec(i, j, 1, 5, 101), src(i, j, nk, 6), rec(i, j, nk, 6, 101), src(i, j, 1, 5)]
    with open("INPUT/bc.txt", "w") as out:
        out.write("".join(rows))


def polar_scale(nodes, i, j, k):
    """mu |G| (sum of the face areas)/volume of cell (i, j, k): the size of the viscous force of one face. Faces: the
    area vector is half the cross product of the diagonals, the centre the mean of the four nodes; the volume is
    (1/3) sum over the faces of centre . area vector, with the area vectors outwards."""
    n = polar_cell(nodes, i, j, k)

    def sub(a, b):
        return tuple(a[m] - b[m] for m in range(3))

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

    faces = [(0, 2, 6, 4), (1, 5, 7, 3), (0, 4, 5, 1), (2, 3, 7, 6), (0, 1, 3, 2), (4, 6, 7, 5)]   # outward, right hand
    total, vol = 0.0, 0.0
    for f in faces:
        p = [n[q] for q in f]
        av = tuple(0.5 * c for c in cross(sub(p[2], p[0]), sub(p[3], p[1])))
        xc = [sum(q[m] for q in p) / 4.0 for m in range(3)]
        total += math.sqrt(sum(c * c for c in av))
        vol += sum(xc[m] * av[m] for m in range(3)) / 3.0
    gnorm = math.sqrt(sum(G_POLAR[m][n] ** 2 for m in range(3) for n in range(3)))
    return MU * gnorm * total / abs(vol)


def polar_linear(ns_path, eu_path):
    """The linear field on the polar mesh with symmetry on the axis: the viscous force (momentum) of the cells with no
    face 400 must vanish, below 1e-11 of polar_scale, the first row next to the axis included."""
    ns, n = read_zone(ns_path)
    eu, _ = read_zone(eu_path)
    ni, nj, nk = POLAR_N
    nodes = polar_nodes()
    worst, worst_row1 = 0.0, 0.0
    for k in range(nk):
        for j in range(nj - 1):
            for i in range(1, ni - 1):
                q = (k * nj + j) * ni + i
                e = max(abs(c) for c in viscous_rhs(ns, eu, q)) / polar_scale(nodes, i, j, k)
                worst = max(worst, e)
                if j == 0:
                    worst_row1 = max(worst_row1, e)
    print(f"    linear velocity on the polar mesh, symmetry on the axis: largest |viscous force| / scale {worst:.3e}, "
          f"on the first row next to the axis {worst_row1:.3e}")
    sys.exit(0 if worst <= 1.0e-11 else 1)


def clean_stop(text):
    """The run aborted (mpi_abort_all) with no NaN and no floating-point trap before it."""
    return "ABORT:" in text and not re.search(r"NaN|Floating-point exception|SIGFPE|Program received signal", text)


def stopped(log):
    """The run stopped at setup with the message on a wall face of zero area on the axis (block 1, a cell of the first
    row, face 3) and the remedy, with no NaN and no floating-point trap before it."""
    text = open(log, errors="replace").read()
    head = re.search(r"\[ERROR\] diffusive-metric = stencil: block 1, cell +(\d+) +1 +(\d+), face 3 \(wall\): a wall on a "
                     r"face of zero area", text)
    remedy = "use the symmetry condition (300) on the axis" in text
    clean = clean_stop(text)
    print(f"    message on a wall face of the axis: {'yes' if head else 'no'}; remedy named: {'yes' if remedy else 'no'}; "
          f"clean stop without NaN or trap: {'yes' if clean else 'no'}")
    sys.exit(0 if head and remedy and clean else 1)


def axis_stop(log):
    """The viscous run on the wedge with the axis at r = 0 stopped at setup with the message on the axis (block 1, a cell
    of the first row, face 3), with no NaN and no floating-point trap before it."""
    text = open(log, errors="replace").read()
    head = re.search(r"\[ERROR\] block 1, cell +(\d+) +1 +1, face 3: a face of zero area on the 2-D axisymmetric wedge", text)
    remedy = "the axis is exactly at r = 0" in text and "Put the axis nodes at a small radius" in text
    clean = clean_stop(text)
    print(f"    message on the axis at r = 0: {'yes' if head else 'no'}; remedy named: {'yes' if remedy else 'no'}; "
          f"clean stop without NaN or trap: {'yes' if clean else 'no'}")
    sys.exit(0 if head and remedy and clean else 1)


def polar_stop(log, metric):
    """The viscous run on the polar mesh (300 on the axis) with diffusive-metric = metric (inverse-mean or mean) stopped at
    setup with the message on a face of zero area on the axis (block 1, a cell of the first row, face 3) that names stencil,
    with no NaN and no floating-point trap before it."""
    text = open(log, errors="replace").read()
    head = re.search(r"\[ERROR\] block 1, cell +(\d+) +1 +(\d+), face 3: a face of zero area \(as on the axis of a 3-D "
                     r"mesh\), "
                     r"with diffusive-metric = " + re.escape(metric), text)
    remedy = "Use diffusive-metric = stencil" in text
    clean = clean_stop(text)
    print(f"    message on the axis with {metric}: {'yes' if head else 'no'}; stencil named: {'yes' if remedy else 'no'};"
          f" clean stop without NaN or trap: {'yes' if clean else 'no'}")
    sys.exit(0 if head and remedy and clean else 1)


def sector_stop(log):
    """The viscous run on the 3-D sector with the axis at r = 0 stopped at setup with the message on the axis (block 1, a
    cell of the first row, face 3), with no NaN and no floating-point trap before it."""
    text = open(log, errors="replace").read()
    head = re.search(r"\[ERROR\] block 1, cell +(\d+) +1 +(\d+), face 3: a face of zero area on the 3-D sector", text)
    remedy = "the axis is exactly at r = 0" in text and "Put the axis nodes at a small radius" in text
    clean = clean_stop(text)
    print(f"    message on the axis at r = 0 of the sector: {'yes' if head else 'no'}; remedy named: {'yes' if remedy else 'no'};"
          f" clean stop without NaN or trap: {'yes' if clean else 'no'}")
    sys.exit(0 if head and remedy and clean else 1)


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


def errors_homann(ns_path, eu_path, geo=False):
    """Errors of the viscous right-hand side of one mesh: radial, radial on r >= L/2 (both over the largest exact radial
    force), radial on the first row over the hoop term 2 mu a f'/r, axial (over the largest exact axial force), azimuthal
    (over the largest exact radial force)."""
    ns, n = read_zone(ns_path)
    eu, n2 = read_zone(eu_path)
    assert n == n2
    rn = radii(n, geo=geo)
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


def errors_swirl(ns_path, eu_path, geo=False):
    """Errors of the viscous right-hand side of one mesh, over the largest exact force of the same component: radial, radial
    on r >= L/2, azimuthal, azimuthal on r >= L/2, axial."""
    ns, n = read_zone(ns_path)
    eu, n2 = read_zone(eu_path)
    assert n == n2
    rn = radii(n, geo=geo)
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


def check_geo(field, a32, a64, g32, g64):
    """The geometric radii: the checks of the uniform radii on G32, G64 and every error that is not round-off within 2 times
    the same error of the uniform radii on N = 64, for the components that the faces of type 200 carry: homann radial (also
    on r >= L/2 and on the first row) and azimuthal; swirl radial, azimuthal and axial. The axial force of homann is printed
    only: an axisymmetric field without swirl has no axial stress on the faces of type 200 (tau_xtheta = 0), and next to the
    axis the faces j of stretched radii take the derivative of u_r along x at the mean radius of the two rows of cells, not
    at the face, an error of first order there, with stencil as with inverse-mean (1.7e-3 and 1.2e-3 of the largest axial
    force on N = 64)."""
    if field == "swirl":
        ok = check_swirl(g32, g64)
        names = (("radial", 0), ("radial on r >= L/2", 1), ("azimuthal", 2), ("azimuthal on r >= L/2", 3))
    else:
        ok = True
        good = g64[0] < 0.02
        ok = ok and good
        print(f"    radial viscous force: error {g32[0]:.3e} (N = 32), {g64[0]:.3e} (N = 64)  [{'ok' if good else 'FAIL'}]")
        good = order(g32[1], g64[1]) > 1.5
        ok = ok and good
        print(f"      on r >= L/2: {g32[1]:.3e} (N = 32), {g64[1]:.3e} (N = 64), observed order {order(g32[1], g64[1]):.2f}"
              f"  [{'ok' if good else 'FAIL'}]")
        good = g64[2] < 1.0e-4
        ok = ok and good
        print(f"      first row off the axis, over the hoop term: {g32[2]:.3e} (N = 32), {g64[2]:.3e} (N = 64)"
              f"  [{'ok' if good else 'FAIL'}]")
        print(f"    axial viscous force (printed only, see wedge_viscous.py): error {g32[3]:.3e} (N = 32), {g64[3]:.3e}"
              f" (N = 64), observed order {order(g32[3], g64[3]):.2f}")
        good = max(g32[4], g64[4]) < 1.0e-9
        ok = ok and good
        print(f"    azimuthal viscous force: {g32[4]:.1e} (N = 32), {g64[4]:.1e} (N = 64)  [{'ok' if good else 'FAIL'}]")
        names = (("radial", 0), ("radial on r >= L/2", 1), ("first row over the hoop term", 2))
    for name, k in names:
        good = g64[k] <= 2.0 * a64[k]
        ok = ok and good
        print(f"    {name} on N = 64: {g64[k]:.3e} against {a64[k]:.3e} of the uniform radii, ratio {g64[k] / a64[k]:.2f}"
              f" (at most 2)  [{'ok' if good else 'FAIL'}]")
    return ok


RADIAL_U = 5.0e3                    # radial: u_r = RADIAL_U r, 1/s
RADIAL_DEV_UNIFORM = 1.0e-11        # radial, uniform radii: the deviation that sets the tolerance of the stretched radii
                                    # must itself be round-off (below 1e-11, as the field U was chosen for)
AXIAL_A = 1.0e3                     # axial: u_x = AXIAL_A x, 1/s
AXIAL_Q, AXIAL_BEND = 1.15, 0.1     # axial: growth of the cells along x and r, bending of the lines j inside the box


def radial_radii(n, law):
    """Node radii of the radial case: uniform or growing by the ratio of LAW (geo1.1, geo1.2), the axis row at 1e-6 dr."""
    w = [1.0] * n if law == "uniform" else [float(law[3:]) ** j for j in range(n)]
    out = [0.0]
    for x in w:
        out.append(out[-1] + SIZE * x / sum(w))
    out[0] = 1.0e-6 * (out[1] - out[0])
    return out


def write_mesh(nodes_xr, cells_u, bc_rows, path_ic="INPUT/ic.tec", path_bc="INPUT/bc.txt"):
    """INPUT/ic.tec of a wedge whose meridional nodes are nodes_xr[i][j] = (x, r) (at theta = -delta/2 and +delta/2) and
    whose cell values are cells_u (u, v, w, j-major), at the uniform density and pressure of the module; INPUT/bc.txt."""
    ni, nj = len(nodes_xr) - 1, len(nodes_xr[0]) - 1
    for p in (path_ic, path_bc):
        if os.path.lexists(p):
            os.remove(p)
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    with open(path_ic, "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        out.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K=2, DATAPACKING=BLOCK, "
                  "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for comp in range(3):
            for k in range(2):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        x, r = nodes_xr[i][j]
                        out.write(f" {(x, r * c, (-1.0, 1.0)[k] * r * s)[comp]:.16E}\n")
        nc = ni * nj
        for col in ([RHO0] * nc, [u[0] for u in cells_u], [u[1] for u in cells_u], [u[2] for u in cells_u], [P0] * nc):
            out.write("".join(f" {x:.16E}\n" for x in col))
    with open(path_bc, "w") as out:
        out.write("".join(bc_rows))


def wedge_rec(i, j, face, t):
    return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)


def wedge_cell_centre(nodes_xr, i, j):
    """Centre of the wedge cell (i, j) (0-based), the mean of its eight nodes: (x, y, 0), y = r cos(delta/2)."""
    q = [nodes_xr[i + a][j + b] for a in (0, 1) for b in (0, 1)]
    return (sum(p[0] for p in q) / 4.0, sum(p[1] for p in q) / 4.0 * math.cos(0.5 * DELTA), 0.0)


def write_radial(n, law):
    """The radial linear field u_r = RADIAL_U r on rectangular cells (uniform x, radii of LAW); faces as in write_case."""
    rn = radial_radii(n, law)
    dx = SIZE / n
    nodes = [[(i * dx, rn[j]) for j in range(n + 1)] for i in range(n + 1)]
    cells = []
    for j in range(n):
        for i in range(n):
            xc = wedge_cell_centre(nodes, i, j)
            cells.append((0.0, RADIAL_U * xc[1], 0.0))
    rows = [wedge_rec(1, j, 1, 400) for j in range(1, n + 1)] + [wedge_rec(n, j, 2, 400) for j in range(1, n + 1)]
    rows += [wedge_rec(i, 1, 3, 300) for i in range(1, n + 1)] + [wedge_rec(i, n, 4, 400) for i in range(1, n + 1)]
    for face in (5, 6):
        rows += [wedge_rec(i, j, face, 200) for j in range(1, n + 1) for i in range(1, n + 1)]
    write_mesh(nodes, cells, rows)


def radial_deviation(n, law, ns_path, eu_path):
    """Largest |force/hoop term - sin^2(delta/2)| on the rows j >= 2 and on the first row, cells with no face 400 (the
    radial force over (2/3) mu U/r_c, r_c the radius of the cell centre)."""
    ns, _ = read_zone(ns_path)
    eu, _ = read_zone(eu_path)
    rn = radial_radii(n, law)
    s2 = math.sin(0.5 * DELTA) ** 2
    dev, dev1 = 0.0, 0.0
    for j in range(n - 1):
        rc = 0.5 * (rn[j] + rn[j + 1]) * math.cos(0.5 * DELTA)
        for i in range(1, n - 1):
            fr = viscous_rhs(ns, eu, j * n + i)[1]
            e = abs(abs(fr) / ((2.0 / 3.0) * MU * RADIAL_U / rc) - s2)
            if j == 0:
                dev1 = max(dev1, e)
            else:
                dev = max(dev, e)
    return dev, dev1


def radial_check(n, law, ns_path, eu_path, tol):
    dev, dev1 = radial_deviation(n, law, ns_path, eu_path)
    ok = dev <= tol and dev1 <= 3.0e-6
    bound = f"at most {tol:.1e}"
    if law == "uniform":
        # tol is at least 10 times this deviation: only its own bound can fail
        ok = ok and dev <= RADIAL_DEV_UNIFORM
        bound = f"at most {RADIAL_DEV_UNIFORM:.0e}, so that the tolerance of the stretched radii cannot float"
    print(f"    radii {law}: |force/hoop term - sin^2(delta/2)| {dev:.3e} on the rows j >= 2 ({bound}), "
          f"{dev1:.3e} on the first row (at most 3e-6)  [{'ok' if ok else 'FAIL'}]")
    return ok


def axial_nodes(n):
    """Meridional nodes of the axial case: x growing by AXIAL_Q from the wall x = 0, r growing by AXIAL_Q from the axis
    (axis row at 1e-6 of the first cell), the lines j bent inside the box (x = 0 and x = SIZE kept)."""
    wx = [AXIAL_Q ** a for a in range(n)]
    xs = [0.0]
    for w in wx:
        xs.append(xs[-1] + SIZE * w / sum(wx))
    rs = [0.0]
    for w in wx:
        rs.append(rs[-1] + SIZE * w / sum(wx))
    rs[0] = 1.0e-6 * rs[1]
    return [[(xs[i] + AXIAL_BEND * rs[j] * math.sin(math.pi * xs[i] / xs[-1]), rs[j]) for j in range(n + 1)]
            for i in range(n + 1)]


def write_axial(n):
    """The axial linear field u = (AXIAL_A x, 0, 0) next to a wall 301 on i = 1 (adiabatic), 400 on i = N and j = N, 300 on
    the axis, 200 on k; INPUT/transport-ref.dat with mu = k = 1e-30 for the reference step."""
    nodes = axial_nodes(n)
    cells = []
    for j in range(n):
        for i in range(n):
            xc = wedge_cell_centre(nodes, i, j)
            cells.append((AXIAL_A * xc[0], 0.0, 0.0))
    rows = []
    for j in range(1, n + 1):
        rows.append(wedge_rec(1, j, 1, 301) + "    0.000000E+00,    0.000000E+00,    0.000000E+00,\n")
    rows += [wedge_rec(n, j, 2, 400) for j in range(1, n + 1)]
    rows += [wedge_rec(i, 1, 3, 300) for i in range(1, n + 1)] + [wedge_rec(i, n, 4, 400) for i in range(1, n + 1)]
    for face in (5, 6):
        rows += [wedge_rec(i, j, face, 200) for j in range(1, n + 1) for i in range(1, n + 1)]
    write_mesh(nodes, cells, rows)
    here = os.path.dirname(os.path.abspath(__file__))
    if os.path.lexists("INPUT/transport-ref.dat"):
        os.remove("INPUT/transport-ref.dat")
    with open(os.path.join(here, "..", "..", "..", "..", "common", "Air-viscous", "transport.dat")) as src, \
            open("INPUT/transport-ref.dat", "w") as dst:
        for line in src:
            p = line.split()
            dst.write(f"{p[0]} 1.0E-30 1.0E-30\n" if len(p) == 3 and re.match(r"^[0-9.Ee+-]+$", p[0]) else line)


def wedge_cell_scale(nodes_xr, i, j):
    """(sum of the face areas)/volume of the wedge cell (i, j) (0-based): area vectors as half the cross product of the
    diagonals, volume (1/3) sum over the faces of centre . area vector (outwards)."""
    c, s = math.cos(0.5 * DELTA), math.sin(0.5 * DELTA)
    node = {}
    for a in (0, 1):
        for b in (0, 1):
            x, r = nodes_xr[i + a][j + b]
            for k in (0, 1):
                node[(a, b, k)] = (x, r * c, (-1.0, 1.0)[k] * r * s)
    nodes = [[[node[(a, b, k)] for k in (0, 1)] for b in (0, 1)] for a in (0, 1)]
    faces = [[(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)], [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
             [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)], [(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
             [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]]
    centre = [sum(node[key][m] for key in node) / 8.0 for m in range(3)]
    total, vol = 0.0, 0.0
    for f in faces:
        p = [nodes[a][b][k] for a, b, k in f]
        d1 = [p[2][m] - p[0][m] for m in range(3)]
        d2 = [p[3][m] - p[1][m] for m in range(3)]
        av = [0.5 * (d1[1] * d2[2] - d1[2] * d2[1]), 0.5 * (d1[2] * d2[0] - d1[0] * d2[2]), 0.5 * (d1[0] * d2[1] - d1[1] * d2[0])]
        xf = [sum(q[m] for q in p) / 4.0 for m in range(3)]
        if sum((xf[m] - centre[m]) * av[m] for m in range(3)) < 0.0:
            av = [-v for v in av]
        total += math.sqrt(sum(v * v for v in av))
        vol += sum(xf[m] * av[m] for m in range(3)) / 3.0
    return total / abs(vol)


def axial_check(n, ns_path, ref_path):
    """The viscous force (momentum, r/dt of NS minus r/dt of the reference) of the cells with no face 400, the wall and axis
    rows included, over mu |a| (sum of the face areas)/volume: below 1e-11."""
    ns, _ = read_zone(ns_path)
    ref, _ = read_zone(ref_path)
    nodes = axial_nodes(n)
    worst, where = 0.0, None
    for j in range(n - 1):
        for i in range(n - 1):
            sc = MU * AXIAL_A * wedge_cell_scale(nodes, i, j)
            for k, f in enumerate(viscous_rhs(ns, ref, j * n + i)):
                if abs(f) / sc >= worst:
                    worst, where = abs(f) / sc, f"component {k + 1} of cell ({i + 1}, {j + 1})"
    ok = worst <= 1.0e-11
    print(f"    axial field next to the wall: largest |viscous force| over its scale {worst:.3e} ({where}), at most 1e-11"
          f"  [{'ok' if ok else 'FAIL'}]")
    return ok


def main():
    what = sys.argv[1]
    if what == "polar":
        write_polar(int(sys.argv[2]))
        return
    if what == "stopped":
        stopped(sys.argv[2])
    if what == "axisstop":
        axis_stop(sys.argv[2])
    if what == "sector":
        write_sector(int(sys.argv[2]), SECTOR_K, len(sys.argv) > 3 and sys.argv[3] == "axis0")
        return
    if what == "sectorstop":
        sector_stop(sys.argv[2])
    if what == "polarstop":
        polar_stop(sys.argv[2], sys.argv[3])
    if what == "linear":
        polar_linear(sys.argv[2], sys.argv[3])
    if what == "radial":
        write_radial(int(sys.argv[2]), sys.argv[3])
        return
    if what == "radialtol":
        dev, _ = radial_deviation(int(sys.argv[2]), "uniform", sys.argv[3], sys.argv[4])
        print(f"{max(1.0e-9, 10.0 * dev):.3e}")
        return
    if what == "radialcheck":
        sys.exit(0 if radial_check(int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], float(sys.argv[6])) else 1)
    if what == "axial":
        write_axial(int(sys.argv[2]))
        return
    if what == "axialcheck":
        sys.exit(0 if axial_check(int(sys.argv[2]), sys.argv[3], sys.argv[4]) else 1)
    field = sys.argv[2] if len(sys.argv) > 2 else ""
    if field not in FIELDS:
        sys.exit(f"unknown field {field}: use one of {', '.join(FIELDS)}")
    opt = sys.argv[4] if len(sys.argv) > 4 else ""
    if what == "case":
        write_case(field, int(sys.argv[3]), opt == "axis0", opt == "geo")
    elif what == "error":
        geo = len(sys.argv) > 5 and sys.argv[5] == "geo"
        errs = (errors_homann if field == "homann" else errors_swirl)(sys.argv[3], sys.argv[4], geo)
        print(",".join(f"{e:.6e}" for e in errs))
    elif what == "check":
        a32 = [float(v) for v in sys.argv[3].split(",")]
        a64 = [float(v) for v in sys.argv[4].split(",")]
        ok = (check_homann if field == "homann" else check_swirl)(a32, a64)
        sys.exit(0 if ok else 1)
    elif what == "checkgeo":
        a32, a64, g32, g64 = ([float(v) for v in sys.argv[k].split(",")] for k in (3, 4, 5, 6))
        sys.exit(0 if check_geo(field, a32, a64, g32, g64) else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
