"""Inputs of the BCInput fast test, and a comparison of its outputs.

The property lines of INPUT/bc.txt must be read as written, or refused with a message: never read wrong in silence.

Inlet 402. One-dimensional duct, 100 cells on 0 <= x <= 1 m, subsonic inlet 402 (total pressure from a table in time, total
temperature) at i = 1, back-pressure outlet (406) at i = 100, slip walls (300) on the j faces, null (0) k faces,
calorically perfect air of test/common/Air. The uniform initial state is the isentropic expansion from (T0, P0) to PB;
the table raises the total pressure by 5 % over the run, so the inlet drives a wave into the duct.

The property line of the inlet names the table file. The same table under four names must give the same run, bit for
bit: p0.dat (the reference), INPUT/tables/p0.dat (a directory), a name longer than 32 characters, and a quoted name with
a blank in it. A table that does not exist, or whose times decrease (decreasing.dat: the rows of p0.dat in reverse
order), must stop the run at set-up. Before the first time of a table the inlet takes its first value: a table that
starts at t = 1e-4 (late.dat, the rows of p0.dat from the second on) must give the run of the same table with its first
value repeated at t = 0 (early.dat), bit for bit.

Inlets with a mass flux or a normal velocity. The mass flux g of the inlets 403, 404 and 420 and the velocity un of the
inlet 408 are positive into the domain on every face: zero, or a negative value on a face 1, 3 or 5, must stop the run at
set-up (an inlet only lets gas in). The case of p0.dat with its inlet at i = 1 (face 1, record 1) of the type and value
that SIGN[KIND] gives: neg404, neg408, neg420 a negative value, zero403, zero404, zero408, zero420 zero. The record of the
nozzle 420 is that of an exit-to-throat area ratio of 2 fed from (T0, P0), calorically perfect (nozzle_record below), with
-g or 0 as its g.

Nozzles 420 that cannot be the record of a nozzle. The thresholds p_sub and p_sup of a nozzle are required and must
satisfy 0 < p_sup <= p_sub < p0, its g must not be above the choked mass flux G* of the reservoir, and its T0 must not be
above one kelvin below the last row of the tables of the gas, where their interpolation ends: otherwise the run stops at
set-up. The case of p0.dat with its inlet (face 1, record 1) the nozzle of area ratio 2, changed as NOZZLE[KIND] says
(G* = 2 g of that nozzle). Fed from the last row of the tables, T_TOP, the same nozzle has the same thresholds and its g
times sqrt(T0/T_TOP), as G* of a calorically perfect gas.

Types that MOSE does not apply. A record of type 409 (a forced outlet), 501 (a manifold) or 506 (surface reactions and
pyrolysis) must stop the run at set-up, with the line of bc.txt, the record and its type: MOSE does not apply these types,
and the face would get no flux. The case of p0.dat with its outlet (record 2, face 2, line 3 of bc.txt) of type 409 or 501,
and the duct of the grain (srm, below) with its outlet of type 506: 409 without a property line, as the reader counts it
(the boundary-condition builder of ATLAS does not write 409); 501 and 506 as that builder writes them, 501 with the source
block and face as two integers, 506 with the pyrolysis and surface-reaction ids as integers, then q_rad, eps_wall and the
mass fraction of the single species.

Usage (from the test directory):
  python3 -B make_case.py case NAME      INPUT/ic.tec, INPUT/bc.txt naming the table NAME, and the table under that name
  python3 -B make_case.py same A B       exit status 0 when the fields A and B hold the same values (cell data, bit for bit)
  python3 -B make_case.py bad KIND       the case of p0.dat with the line of bc.txt that BAD[KIND] gives
  python3 -B make_case.py sign KIND      the case of p0.dat with its inlet of the type and value that SIGN[KIND] gives
  python3 -B make_case.py nozzle KIND    the case of p0.dat with its inlet the nozzle of area ratio 2 changed as
                                         NOZZLE[KIND] says
  python3 -B make_case.py srm [rest] [sst]   the duct with a viscous gas and nu~, its inlet an SRM grain boundary (502);
                                         rest: the gas at rest at T0 and PB; sst: k and omega instead of nu~
  python3 -B make_case.py unapplied TYPE the case of p0.dat (TYPE 409 or 501) or the duct of the grain (506) with its outlet
                                         (record 2) of TYPE
  python3 -B make_case.py srmcheck DIAG   exit status 0 when the change of rho nu~ over one step is that of the wall value
                                         nu~ = 0 in the cell next to the grain, and that of the distance to the grain in
                                         the other cells (below)

SRM grain with Spalart-Allmaras. Its turbulence values are those of a smooth wall, nu~_w = 0, and the ghost cells behind it
hold -nu~ of the cell (2 nu~_w - nu~_1) at the reflection of the cell centre, as behind a wall; the grain is a wall for the
wall distance d, the only one of the duct. With the gas at rest and rho nu~ = m uniform, nu~ = nu (chi = 1: fv2 > 0), the
right-hand side of rho nu~ in a cell is its source rho nu~^2 F / d^2 (F = cb1 fv2/kappa^2 - cw1 fw(r), r = min(1/fv2, 10):
Spalart-Allmaras at rest, as in the WallDistance test), d = (i - 1/2) h in cell i, and in the cell next to the grain also
the flux through the grain, -(mu/sigma) (2 nu~/h)/h per unit volume, and the cb2 term of the gradient nu~/h that the ghost
cell gives, (cb2/sigma) rho (nu~/h)^2. Times sigma rho h^2 / m: cb2 m - 2 mu + sigma m F/(i - 1/2)^2 in the first cell,
sigma m F/(i - 1/2)^2 in the others. A ghost cell that extrapolates nu~ (3 q1 - 3 q2 + q3, uniform) gives no gradient there
(-2 mu instead of cb2 m - 2 mu), and a grain left out of the wall distance no source (d is then 1e8, with no wall).
"""
import math
import os
import sys

N = 100
LX, LY, LZ = 1.0, 0.01, 1.0e-3
T0, P0, PB = 300.0, 1.2e5, 1.0e5
GAMMA = 1.4
RGAS = 8314.46261815324 / 28.970418


def uniform_state():
    mach = math.sqrt(2.0 / (GAMMA - 1.0) * ((P0 / PB) ** ((GAMMA - 1.0) / GAMMA) - 1.0))
    t = T0 / (1.0 + 0.5 * (GAMMA - 1.0) * mach ** 2)
    c = math.sqrt(GAMMA * RGAS * t)
    return PB / (RGAS * t), mach * c, PB


def write_case(name, state=None):
    rho, u, p = state or uniform_state()
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={N + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(N + 1):
                        f.write(f" {(LX * i / N, LY * j, LZ * (k - 0.5))[axis]:.15E}\n")
        for value in (rho, u, 0.0, 0.0, p):
            f.write(f" {value:.15E}\n" * N)
    written = f"'{name}'" if " " in name else name
    with open("INPUT/bc.txt", "w") as f:
        f.write(f"       1       1       1       1       1     402\n{T0}, {written}, normal, normal, 1.0, 1.0\n")
        f.write(f"       1     {N:3d}       1       1       2     406\n{PB:.6E}, 1.0\n")
        for face, btype in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, N + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {btype:3d}\n")
    if name == "missing.dat":
        return
    if os.path.dirname(name):
        os.makedirs(os.path.dirname(name), exist_ok=True)
    rows = [(1.0e-4 * n, P0 * (1.0 + 0.05 * n / 10) / 1.0e5) for n in range(11)]   # bar: 1.2 bar rising by 5 % in 1 ms
    if name == "decreasing.dat":
        rows.reverse()
    elif name == "late.dat":
        rows = rows[1:]
    elif name == "early.dat":
        rows = [(0.0, rows[1][1])] + rows[1:]
    with open(name, "w") as f:
        for t, p in rows:
            f.write(f"{t:.15E} {p:.15E}\n")


BAD = {                                              # line of INPUT/bc.txt (0-based) and what it becomes
    "empty": (3, f"{PB:.6E}, ,\n"),                 # line 4, the outlet (record 2): rel_fac empty
    "short": (3, f"{PB:.6E}\n"),                    # line 4, the outlet: rel_fac missing
    "word": (3, f"{PB:.6E}, x\n"),                  # line 4, the outlet: a word for rel_fac
    "angle": (1, f"{T0}, p0.dat, nromal, normal, 1.0, 1.0\n"),   # line 2, the inlet (record 1): alpha misspelt
}


def write_bad(kind):
    write_case("p0.dat")
    n, text = BAD[kind]
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[n] = text
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


CP = 1004.5                                          # test/common/Air: cp, J/(kg K)
T_TOP = 10000.0                                      # test/common/Air: the last row of the tables, K
G_IN, UN_IN = 120.0, 100.0                           # a mass flux, kg/(m2 s), and a velocity, m/s
SIGN = {                                             # record 1 (face 1): type, value (None: -g of the nozzle record)
    "neg404": (404, -G_IN), "neg408": (408, -UN_IN), "neg420": (420, None),
    "zero403": (403, 0.0), "zero404": (404, 0.0), "zero408": (408, 0.0), "zero420": (420, 0.0),
}


def nozzle_record(ae_at):
    """p_sub, p_sup and g of a nozzle of exit-to-throat area ratio ae_at fed from (T0, P0), calorically perfect (cp = CP):
    the exit pressures of the subsonic and of the supersonic isentropic flows, and g = G*/ae_at, G* the choked mass flux.
    The exit Mach numbers by bisection on the two branches of the area ratio A/A*(M), decreasing for M < 1 and increasing
    for M > 1, from brackets that hold the root for 1 < ae_at < 1e5 (200 halvings: the bracket shrinks to the spacing of
    the doubles)."""
    gam = CP / (CP - RGAS)
    e = (gam + 1.0) / (2.0 * (gam - 1.0))

    def area(m):
        return ((2.0 / (gam + 1.0)) * (1.0 + 0.5 * (gam - 1.0) * m * m)) ** e / m

    def root(lo, hi, increasing):
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (area(mid) < ae_at) == increasing:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def pressure(m):
        return P0 * (1.0 + 0.5 * (gam - 1.0) * m * m) ** (-gam / (gam - 1.0))

    gstar = P0 * math.sqrt(gam / (RGAS * T0)) * (2.0 / (gam + 1.0)) ** e
    return pressure(root(1.0e-6, 1.0, False)), pressure(root(1.0, 50.0, True)), gstar / ae_at


def write_sign(kind):
    write_case("p0.dat")
    btype, value = SIGN[kind]
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[0] = f"       1       1       1       1       1     {btype}\n"
    if btype == 420:                                 # T0, p0, p_sub, p_sup, g, rel_fac, mass fraction
        p_sub, p_sup, g = nozzle_record(2.0)
        lines[1] = f"{T0}, {P0}, {p_sub:.17E}, {p_sup:.17E}, {-g if value is None else value:.17E}, 1.0, 1.0\n"
    else:                                            # T0 (403) or T (404, 408), g or un, alpha, beta, rel_fac, mass fraction
        lines[1] = f"{T0}, {value}, normal, normal, 1.0, 1.0\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


NOZZLE = {                                           # record 1 (face 1): the nozzle of area ratio 2, then what changes
    "psub-empty": "p_sub empty",
    "swapped": "p_sub and p_sup swapped",
    "psub-p0": "p_sub at p0",
    "psup-zero": "p_sup zero",
    "choked": "g at 1.01 G*",
    "t0-above": "T0 at the last row of the tables, T_TOP",
}


def write_nozzle(kind):
    if kind not in NOZZLE:
        sys.exit(f"unknown nozzle {kind}")
    write_case("p0.dat")
    p_sub, p_sup, g = nozzle_record(2.0)
    t0 = T0
    if kind == "swapped":
        p_sub, p_sup = p_sup, p_sub
    elif kind == "psub-p0":
        p_sub = P0
    elif kind == "psup-zero":
        p_sup = 0.0
    elif kind == "choked":                           # G* = 2 g, the choked mass flux of the reservoir
        g = 1.01 * 2.0 * g
    elif kind == "t0-above":                         # the same nozzle fed from T_TOP: G* scales as 1/sqrt(T0)
        t0 = T_TOP
        g = g * math.sqrt(T0 / T_TOP)
    psub = " " if kind == "psub-empty" else f"{p_sub:.17E}"
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[0] = "       1       1       1       1       1     420\n"
    lines[1] = f"{t0}, {P0}, {psub}, {p_sup:.17E}, {g:.17E}, 1.0, 1.0\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


MU = 1.8467e-5                                       # a viscosity of air, Pa s: nu~ = 3 MU / rho needs only to be positive
MU_GAS = 1.1858686e-5                                # the viscosity of test/common/Air-viscous (constant), Pa s
SIGMA, CB2 = 2.0 / 3.0, 0.622                        # Spalart-Allmaras
CB1, KAPPA, CW2, CW3, CV1 = 0.1355, 0.41, 0.3, 2.0, 7.1
CW1 = CB1 / KAPPA ** 2 + (1.0 + CB2) / SIGMA


def sa_factor(chi):
    """F = cb1 fv2/kappa^2 - cw1 fw: the right-hand side of rho nu~ at rest is rho nu~^2 F / d^2 (fv2 > 0)."""
    fv1 = chi ** 3 / (chi ** 3 + CV1 ** 3)
    fv2 = 1.0 - chi / (1.0 + chi * fv1)
    r = min(1.0 / fv2, 10.0)
    g = r + CW2 * (r ** 6 - r)
    fw = g * ((1.0 + CW3 ** 6) / (g ** 6 + CW3 ** 6)) ** (1.0 / 6.0)
    return CB1 * fv2 / KAPPA ** 2 - CW1 * fw


def write_srm(rest=False, sst=False):
    """The duct of p0.dat with rho nu~ = 3 MU (mi_tilde) for the Spalart-Allmaras model, and record 1 (the inlet) an SRM
    grain boundary (502): T_af, a, n, p_ref, rho_grain, SF_geo and the mass fraction of the single species. rest: the gas
    at rest at T0 and PB, with rho nu~ = MU_GAS (nu~ = nu); sst: k = 1 m2/s2 and omega = 1e3 1/s (kappa = rho k, omega =
    rho omega) for SST instead."""
    rho, u, p = (PB / (RGAS * T0), 0.0, PB) if rest else uniform_state()
    write_case("p0.dat", (rho, u, p))
    with open("INPUT/ic.tec") as f:
        lines = f.read().splitlines(True)
    if sst:
        lines[0] = lines[0].replace('"p"\n', '"p" "rhoK" "rhoOmega"\n')
        lines[1] = lines[1].replace("[4-8]=CELLCENTERED", "[4-10]=CELLCENTERED")
        turb = f" {rho * 1.0:.15E}\n" * N + f" {rho * 1.0e3:.15E}\n" * N
    else:
        lines[0] = lines[0].replace('"p"\n', '"p" "rhoNuTilde"\n')
        lines[1] = lines[1].replace("[4-8]=CELLCENTERED", "[4-9]=CELLCENTERED")
        turb = f" {(MU_GAS if rest else 3.0 * MU):.15E}\n" * N
    with open("INPUT/ic.tec", "w") as f:
        f.write("".join(lines) + turb)
    with open("INPUT/bc.txt") as f:
        bc = f.readlines()
    bc[0] = bc[0].replace("     402\n", "     502\n")
    bc[1] = "3.0E+03, 1.0E-05, 0.4, 1.0E+06, 1.7E+03, 1.0, 1.0\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(bc))


UNAPPLIED = ("409", "501", "506")                    # types that MOSE does not apply


def write_unapplied(kind):
    """The case of p0.dat (409, 501) or the duct of the grain (506, write_srm) with record 2, the outlet on face 2 at line
    3 of INPUT/bc.txt, of type kind: 409 without a property line, as the reader counts it (the boundary-condition builder
    of ATLAS does not write 409); 501 and 506 as that builder writes them: 501 with the source block and face, 1 and 1, as
    two integers; 506 with the pyrolysis and surface-reaction ids as integers, then q_rad, eps_wall and the mass fraction of
    the single species, each value followed by a comma."""
    if kind not in UNAPPLIED:
        sys.exit(f"unknown type {kind}")
    if kind == "506":
        write_srm()
    else:
        write_case("p0.dat")
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[2] = f"       1     {N:3d}       1       1       2     {kind}\n"
    if kind == "409":
        del lines[3]
    elif kind == "501":
        lines[3] = "       1,       1,\n"
    else:
        lines[3] = f"       1,       1,{1.0e6:24.16E},{0.9:24.16E},{1.0:24.16E},\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), by variable name."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = lines[0].split('=', 1)[1].replace('"', ' ').split()
    head = lines[1].replace(',', ' ').split()
    ni, nj, nk = (int(next(w for w in head if w.startswith(k + '=')).split('=')[1]) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out


def srm_check(diag):
    """The right-hand side of rho nu~ of every cell, its increment in the diagnostic output over the step dt, times
    sigma rho h^2 / m: in the cell next to the grain cb2 m - 2 mu (the wall value nu~_w = 0 in the ghost cell) plus
    the source with the distance to the grain, d = h/2; in the other cells the source alone, d = (i - 1/2) h, but in the
    last one, next to the outlet."""
    d = read_zone(diag)
    rho, m, h = PB / (RGAS * T0), MU_GAS, LX / N
    r = [inc / t * SIGMA * rho * h * h / m for inc, t in zip(d["rhoNuTilde"], d["dt"])]
    src = [SIGMA * m * sa_factor(m / MU_GAS) / (i + 0.5) ** 2 for i in range(N)]
    expect = CB2 * m - 2.0 * MU_GAS + src[0]
    err = abs(r[0] - expect) / (2.0 * MU_GAS)
    errd = max(abs(r[i] / src[i] - 1.0) for i in range(1, N - 1))
    print(f"    rho nu~ next to the grain, times sigma rho h^2/(m dt): {r[0]:.6e}, from the wall value nu~ = 0 and the "
          f"distance h/2 {expect:.6e} (an extrapolated ghost cell: {expect - CB2 * m:.6e}; no distance to the grain: "
          f"{CB2 * m - 2.0 * MU_GAS:.6e}); difference {err:.1e} of 2 mu")
    print(f"    the other cells, the source with the distance to the grain: largest relative difference {errd:.1e} "
          f"(second cell {r[1]:.6e} against {src[1]:.6e})")
    sys.exit(0 if err < 1e-6 and errd < 1e-9 else 1)


def cell_values(path):
    with open(path) as f:
        lines = f.read().splitlines()
    values = " ".join(lines[2:]).split()
    return values[3 * (N + 1) * 4:]                  # after the node coordinates: the cell values, as written


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "bad":
        write_bad(sys.argv[2])
    elif what == "sign":
        write_sign(sys.argv[2])
    elif what == "nozzle":
        write_nozzle(sys.argv[2])
    elif what == "srm":
        write_srm("rest" in sys.argv[2:], "sst" in sys.argv[2:])
    elif what == "unapplied":
        write_unapplied(sys.argv[2])
    elif what == "srmcheck":
        srm_check(sys.argv[2])
    elif what == "same":
        a, b = cell_values(sys.argv[2]), cell_values(sys.argv[3])
        sys.exit(0 if a == b and len(a) > 0 else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
