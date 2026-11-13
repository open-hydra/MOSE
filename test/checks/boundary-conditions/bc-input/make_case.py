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

Types that MOSE does not apply. A record of type 409 (a forced outlet) or 501 (a manifold) must stop the run at set-up,
with the line of bc.txt, the record and its type: MOSE does not apply these types, and the face would get no flux. The
case of p0.dat with its outlet (record 2, face 2, line 3 of bc.txt) of type 409 or 501: 409 without a property line, as
the reader counts it (the boundary-condition builder of ATLAS does not write 409); 501 as that builder writes it, with the
source block and face as two integers.

Usage (from the test directory):
  python3 -B make_case.py case NAME      INPUT/ic.tec, INPUT/bc.txt naming the table NAME, and the table under that name
  python3 -B make_case.py same A B       exit status 0 when the fields A and B hold the same values (cell data, bit for bit)
  python3 -B make_case.py bad KIND       the case of p0.dat with the line of bc.txt that BAD[KIND] gives
  python3 -B make_case.py sign KIND      the case of p0.dat with its inlet of the type and value that SIGN[KIND] gives
  python3 -B make_case.py nozzle KIND    the case of p0.dat with its inlet the nozzle of area ratio 2 changed as
                                         NOZZLE[KIND] says
  python3 -B make_case.py srm [rest] [sst] [inflow]   the duct with a viscous gas and nu~, its inlet an SRM grain
                                         boundary (502); rest: the gas at rest at T0 and PB; sst: k and omega instead of
                                         nu~; inflow: the grain in mode 2 (an inflow) instead of mode 1 (a blowing wall)
  python3 -B make_case.py srmbad KIND    the duct of the grain with the record that SRM_BAD[KIND] gives
  python3 -B make_case.py srmtwogas [inflow]   the inviscid duct full of N2 at rest, its inlet a grain blowing in CH4
                                         (mode 1; inflow: mode 2)
  python3 -B make_case.py srmtwogascheck DIAG [DIAG1]   exit status 0 when one step gives the momentum of the blown CH4
                                         with the gas constant of CH4 at the static temperature of total T_af, and the
                                         energy of the total enthalpy h(T_af); [the changes of the grain of mode 1]
  python3 -B make_case.py srminflowcheck DIAG DIAGWALL   exit status 0 when one step from rest lets in the mass flux of
                                         the law through the inflow grain, as through the grain of mode 1 (DIAGWALL), and
                                         the grain is not a wall for the distance (below)
  python3 -B make_case.py unapplied TYPE the case of p0.dat with its outlet (record 2) of TYPE (409 or 501)
  python3 -B make_case.py refused KIND   the duct of the grain with its outlet (record 2) of the type and line of
                                         REFUSED[KIND], refused at set-up
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
    rows = [(1.0e-4 * n, P0 * (1.0 + 0.05 * n / 10)) for n in range(11)]   # Pa: P0 rising by 5 % in 1 ms
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


# The record of the grain (502): mdot_ref [kg/(m2 s)], n, p_ref [Pa], T_af [K], mode (1 a wall that blows the gas in, 2 an
# inflow), the mass fraction of the single species, and in mode 2 the turbulence values (specific: nu~, or k and omega).
SRM_RECORD = dict(mdot_ref=1.0e-5, n=0.4, p_ref=1.0e6, T_af=3.0e3, mode=1)
NU_IN, K_IN, OMEGA_IN = 1.5e-5, 1.0, 1.0e3                 # the turbulence values of the inflow grain
# Records that the reader refuses (make_case.py srmbad KIND): the change of each to SRM_RECORD, and short, the inflow
# grain with Spalart-Allmaras without its turbulence value
SRM_BAD = {
    "mode3":   dict(mode=3),
    "mode15":  dict(mode=1.5),
    "pref0":   dict(p_ref=0.0),
    "mdotneg": dict(mdot_ref=-1.0e-5),
    "inflow0": dict(mode=2, mdot_ref=0.0),
    "taf0":    dict(T_af=0.0),
    "tafhigh": dict(T_af=5.0e3),                     # the last row of test/common/Air-viscous: above 4999 K
    "short":   dict(mode=2),
}


def write_srm(rest=False, sst=False, inflow=False, bad=None):
    """The duct of p0.dat with rho nu~ = 3 MU (mi_tilde) for the Spalart-Allmaras model, and record 1 (the inlet) an SRM
    grain boundary (502): SRM_RECORD, the mass fraction of the single species and, in mode 2, the turbulence values. rest:
    the gas at rest at T0 and PB, with rho nu~ = MU_GAS (nu~ = nu); sst: k = 1 m2/s2 and omega = 1e3 1/s (kappa = rho k,
    omega = rho omega) for SST instead; inflow: the grain in mode 2, with NU_IN (or K_IN, OMEGA_IN); bad: the record of
    SRM_BAD[bad]."""
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
    rec = dict(SRM_RECORD, mode=2 if inflow else 1)
    if bad is not None:
        rec.update(SRM_BAD[bad])
    line = f"{rec['mdot_ref']:.1E}, {rec['n']}, {rec['p_ref']:.1E}, {rec['T_af']:.1E}, {rec['mode']:g}, 1.0"
    if rec["mode"] == 2 and bad != "short":
        line += f", {K_IN:.1E}, {OMEGA_IN:.1E}" if sst else f", {NU_IN:.1E}"
    with open("INPUT/bc.txt") as f:
        bc = f.readlines()
    bc[0] = bc[0].replace("     402\n", "     502\n")
    bc[1] = line + "\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(bc))


# The grain of another gas than that of the duct: the four species of test/1D/diff1 (CH4, H2O, N2, O2), the duct full of N2
# at rest at T0 and PB, the grain (mode 1, inviscid) blowing in CH4 at T_af with mdot_ref = 1 kg/(m2 s) at p_ref = PB.
WM_TWOGAS = (16.043, 18.015, 28.014, 31.998)
TAF_TWOGAS = 3.0e3


def write_srm_twogas(inflow=False):
    """The duct of write_case at rest, N2 at T0 and PB (the other species 1e-20), and record 1 (the inlet) a grain of mode
    1 (inflow: mode 2) that blows in CH4 of total temperature TAF_TWOGAS with mdot_ref = 1 kg/(m2 s) at p_ref = PB."""
    write_case("p0.dat", (PB / (8314.46261815324 / WM_TWOGAS[2] * T0), 0.0, PB))
    rho_n2 = PB / (8314.46261815324 / WM_TWOGAS[2] * T0)
    with open("INPUT/ic.tec") as f:
        lines = f.read().splitlines(True)
    nodes = lines[2:2 + 3 * 4 * (N + 1)]
    cells = lines[2 + 3 * 4 * (N + 1):]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "rho2" "rho3" "rho4" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={N + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-11]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        f.write("".join(nodes))
        for value in (1.0e-20, 1.0e-20, rho_n2, 1.0e-20):
            f.write(f" {value:.15E}\n" * N)
        f.write("".join(cells[N:]))                      # u, v, w, p of write_case
    with open("INPUT/bc.txt") as f:
        bc = f.readlines()
    bc[0] = bc[0].replace("     402\n", "     502\n")
    bc[1] = f"1.0, 0.4, {PB:.6E}, {TAF_TWOGAS:.1E}, {2 if inflow else 1}, 1.0, 0.0, 0.0, 0.0\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(bc))


def ch4_table():
    """Temperature, cp and h of CH4 from test/1D/diff1/INPUT/thermo.dat (rows T = 1, 2, ... K), as MOSE loads them."""
    rows, zone = {}, None
    with open("../../../1D/diff1/INPUT/thermo.dat") as f:
        for ln in f:
            if ln.startswith("ZONE"):
                zone = ln.split('"')[1]
                continue
            t = ln.split()
            if zone == "CH4" and len(t) == 5:
                rows[int(float(t[0]))] = (float(t[1]), float(t[2]))
    return rows


def h_ch4(rows, T):
    """h of CH4 at T, linear between the rows of the table (f_tabT)."""
    i = int(T)
    return rows[i][1] + (rows[i + 1][1] - rows[i][1]) * (T - i)


def srm_twogas_check(diag, diag_mode1=None):
    """One explicit step of the inviscid duct of write_srm_twogas from rest. The products cross the face at the static
    temperature T_w of h(T_w) + u^2/2 = h(T_af), u = mdot R T_w / p (T_af the total temperature, mdot = 1 kg/(m2 s), h of
    the table of CH4 as MOSE interpolates it, solved here by bisection). The pressure of the face, p of the cell, and that
    of the first face between cells cancel, and the cell next to the grain gains the momentum mdot u = mdot^2 R_w T_w / p
    per unit area and time: the gas constant of the face state, R_w = d(rho u) h p / (dt mdot^2 T_w), is that of the
    composition the grain imposes, CH4, not that of the gas of the cell, N2; the energy mdot (h(T_w) + u^2/2) =
    mdot h(T_af), the total enthalpy, not mdot (h(T_af) + u^2/2); the mass flux mdot. diag_mode1: the grain of mode 1 of
    the same step, whose changes those of mode 2 (an inflow, the inlet 404) must equal: rho and rho e within 1e-12
    relative; rho u within 3 ulp of the face pressure. The face pressure, p of the cell (PB), and that of the first face
    between cells cancel in the change of rho u, mdot u, PB/(mdot u) = 6.4e3 times smaller than either face flux: the two
    paths (the wall flux of mode 1, the inlet 404 of mode 2) round p + mdot u once per stage each, and RK3 sums its three
    stages with two more roundings at the size of PB, so their changes of rho u may differ by up to 3 ulp(PB) of flux,
    2.8e-12 of the change (measured with gfortran 13 on aarch64: 1.5 ulp with -O2 -g -fcheck=bounds, the TESTING flags
    of the CI, 0 with -O3 -march=native)."""
    d = read_zone(diag)
    h, mdot = LX / N, 1.0
    r_ch4, r_n2 = 8314.46261815324 / WM_TWOGAS[0], 8314.46261815324 / WM_TWOGAS[2]
    rows = ch4_table()
    h0 = h_ch4(rows, TAF_TWOGAS)
    lo, hi = TAF_TWOGAS - 50.0, TAF_TWOGAS
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if h_ch4(rows, mid) + 0.5 * (mdot * r_ch4 * mid / PB) ** 2 < h0:
            lo = mid
        else:
            hi = mid
    t_w = 0.5 * (lo + hi)
    u = mdot * r_ch4 * t_w / PB
    r_w = abs(d["rhou"][0]) / d["dt"][0] * h * PB / (mdot ** 2 * t_w)
    err = abs(r_w / r_ch4 - 1.0)
    e_in = d["rhoe"][0] / d["dt"][0] * h
    erre = abs(e_in / (mdot * h0) - 1.0)
    static = abs((mdot * (h0 + 0.5 * u * u)) / (mdot * h0) - 1.0)
    dm = d["rho"][0] / d["dt"][0] * h
    print(f"    static temperature of the face {t_w:.9f} K (total {TAF_TWOGAS:g} K, u = {u:.4f} m/s); gas constant of the "
          f"face state from the momentum: {r_w:.9f} J/(kg K), CH4 {r_ch4:.9f} (difference {err:.1e}), N2 {r_n2:.9f}")
    print(f"    energy in: {e_in:.9e} W/m2, mdot h(T_af) {mdot * h0:.9e} (difference {erre:.1e}; with T_af static: "
          f"{static:.1e}); mass flux {dm:.9f} kg/(m2 s)")
    ok = err < 1e-9 and erre < 1e-9 and abs(dm - mdot) < 1e-6
    if diag_mode1 is not None:
        w = read_zone(diag_mode1)
        dev = max(abs(d[k][0] / d["dt"][0] / (w[k][0] / w["dt"][0]) - 1.0) for k in ("rho", "rhoe"))
        dev_u = abs(d["rhou"][0] / d["dt"][0] - w["rhou"][0] / w["dt"][0]) * h / math.ulp(PB)
        print(f"    the changes of mode 2 against those of mode 1: rho and rho e within {dev:.1e} (relative), rho u within "
              f"{dev_u:.2f} ulp of the face pressure ({dev_u * math.ulp(PB) / (mdot * u):.1e} of the change)")
        ok = ok and dev < 1e-12 and dev_u <= 3.0
    sys.exit(0 if ok else 1)


UNAPPLIED = ("409", "501")                           # types that MOSE does not apply


# Records refused at set-up on the duct of the grain (make_case.py refused KIND): record 2, the outlet on face 2 at line 3
# of INPUT/bc.txt, of the type and property line of REFUSED[KIND], as the boundary-condition builder of ATLAS writes it,
# each value followed by a comma; the gas of the duct of the grain, test/common/Air-viscous, has one species.
REFUSED = {
    # reac_id, q_rad, eps_wall: the surface reactions on a gas without their eight species
    "gsi505-species": ("505", f"       1,{0.0:24.16E},{0.0:24.16E},"),
    # pyro_id, reac_id, phi, h_v0, q_rad, eps_wall (BC 506): with phi = 0 no table, the same refusal as 505
    "gsi506-species": ("506", f"       0,       1,{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},"),
    # the line of 506 of an ATLAS before the pyrolysing ablator, one species: pyro_id, reac_id, q_rad, eps_wall, massf
    "gsi506-old": ("506", f"       1,       1,{1.0e6:24.16E},{0.9:24.16E},{1.0:24.16E},"),
    # a seventh value (the line of an older ATLAS with three species holds seven numbers; with two species it holds six
    # and is read as a valid line, DEC-1)
    "gsi506-extra": ("506", f"       1,       1,{0.25:24.16E},{0.0:24.16E},{1.0e4:24.16E},{0.0:24.16E},{1.0:24.16E},"),
    "gsi506-reac": ("506", f"       0,       2,{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},"),
    "gsi506-phi": ("506", f"       0,       1,{-0.5:24.16E},{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},"),
    "gsi506-pyro": ("506", f"       0,       1,{0.25:24.16E},{0.0:24.16E},{0.0:24.16E},{0.0:24.16E},"),
}


def write_refused(kind):
    """The duct of the grain (write_srm) with record 2 of the type and property line of REFUSED[kind]."""
    if kind not in REFUSED:
        sys.exit(f"unknown kind {kind}")
    btype, line = REFUSED[kind]
    write_srm()
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[2] = f"       1     {N:3d}       1       1       2     {btype}\n"
    lines[3] = line + "\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


def write_unapplied(kind):
    """The case of p0.dat with record 2, the outlet on face 2 at line 3 of INPUT/bc.txt, of type kind: 409 without a
    property line, as the reader counts it (the boundary-condition builder of ATLAS does not write 409); 501 as that
    builder writes it, with the source block and face, 1 and 1, as two integers, each followed by a comma."""
    if kind not in UNAPPLIED:
        sys.exit(f"unknown type {kind}")
    write_case("p0.dat")
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[2] = f"       1     {N:3d}       1       1       2     {kind}\n"
    if kind == "409":
        del lines[3]
    else:
        lines[3] = "       1,       1,\n"
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


def srm_inflow_check(diag, diag_wall):
    """The grain in mode 2 (an inflow) over one step from rest, as srm_check: the face is not closed by the gas of the
    cell at rest, the products enter with the mass flux of the law at the pressure of the cell, so the density of the
    cell next to the grain grows by g dt / h, g = mdot_ref (PB/p_ref)^n (to 1e-8: the grain of mode 1, diag_wall, the
    same step from the same state, gives 3.7e-9 from it), the change of the grain of mode 1 to 1e-14, and that of the
    other cells, but the last one, not at all; the grain is not a wall for the wall distance (an inflow is applied as the
    inlet 404): with no wall in the duct the source of rho nu~ of the other cells is zero (d = 1e8), where the grain of
    mode 1 gives that of d = (i - 1/2) h."""
    d = read_zone(diag)
    w = read_zone(diag_wall)
    rho, m, h = PB / (RGAS * T0), MU_GAS, LX / N
    g = SRM_RECORD["mdot_ref"] * (PB / SRM_RECORD["p_ref"]) ** SRM_RECORD["n"]
    drho = d["rho"]
    err = abs(drho[0] / d["dt"][0] * h / g - 1.0)
    errw = abs(drho[0] / d["dt"][0] / (w["rho"][0] / w["dt"][0]) - 1.0)
    rest = max(abs(x) for x in drho[1:N - 1]) / abs(drho[0])
    r = [inc / t * SIGMA * rho * h * h / m for inc, t in zip(d["rhoNuTilde"], d["dt"])]
    src = [SIGMA * m * sa_factor(m / MU_GAS) / (i + 0.5) ** 2 for i in range(N)]
    rel = max(abs(r[i] / src[i]) for i in range(1, N - 1))
    print(f"    density next to the grain over the step, times h/dt: {drho[0] / d['dt'][0] * h:.9e} kg/(m2 s) against "
          f"g = {g:.9e}, difference {err:.1e}, from the grain of mode 1 {errw:.1e}; the other cells, but the last: "
          f"largest change {rest:.1e} of it")
    print(f"    rho nu~ of the other cells, against the source with a wall at the grain: largest ratio {rel:.1e}")
    sys.exit(0 if err < 1e-8 and errw < 1e-14 and rest == 0.0 and rel < 1e-9 else 1)


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
        write_srm("rest" in sys.argv[2:], "sst" in sys.argv[2:], "inflow" in sys.argv[2:])
    elif what == "srmbad":
        write_srm(bad=sys.argv[2])
    elif what == "srmtwogas":
        write_srm_twogas("inflow" in sys.argv[2:])
    elif what == "srmtwogascheck":
        srm_twogas_check(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif what == "srminflowcheck":
        srm_inflow_check(sys.argv[2], sys.argv[3])
    elif what == "unapplied":
        write_unapplied(sys.argv[2])
    elif what == "refused":
        write_refused(sys.argv[2])
    elif what == "srmcheck":
        srm_check(sys.argv[2])
    elif what == "same":
        a, b = cell_values(sys.argv[2]), cell_values(sys.argv[3])
        sys.exit(0 if a == b and len(a) > 0 else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
