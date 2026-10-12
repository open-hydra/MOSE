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

Usage (from the test directory):
  python3 -B make_case.py case NAME      INPUT/ic.tec, INPUT/bc.txt naming the table NAME, and the table under that name
  python3 -B make_case.py same A B       exit status 0 when the fields A and B hold the same values (cell data, bit for bit)
  python3 -B make_case.py bad KIND       the case of p0.dat with the line of bc.txt that BAD[KIND] gives
  python3 -B make_case.py srm [rest] [sst]   the duct with a viscous gas and nu~, its inlet an SRM grain boundary (502);
                                         rest: the gas at rest at T0 and PB; sst: k and omega instead of nu~
  python3 -B make_case.py srmcheck FIELD DIAG   exit status 0 when the change of rho nu~ over one step in the cell next to
                                         the grain is that of the wall value nu~ = 0 (below)

SRM grain with Spalart-Allmaras. Its turbulence values are those of a smooth wall, nu~_w = 0, and the ghost cells behind it
hold -nu~ of the cell (2 nu~_w - nu~_1) at the reflection of the cell centre, as behind a wall. With the gas at rest, rho nu~
= m uniform, and no wall for the wall distance (production and destruction vanish with it), one explicit step changes rho
nu~ in the cell next to the grain by the flux through the grain, -(mu/sigma) (2 nu~/h)/h per unit volume, and the cb2 term
of the gradient nu~/h that the ghost cell gives, (cb2/sigma) rho (nu~/h)^2: (m/(sigma rho h^2)) (cb2 m - 2 mu). A ghost
cell that extrapolates nu~ (3 q1 - 3 q2 + q3, uniform) gives no gradient there: (m/(sigma rho h^2)) (-2 mu).
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


MU = 1.8467e-5                                       # a viscosity of air, Pa s: nu~ = 3 MU / rho needs only to be positive
MU_GAS = 1.1858686e-5                                # the viscosity of test/common/Air-viscous (constant), Pa s
SIGMA, CB2 = 2.0 / 3.0, 0.622                        # Spalart-Allmaras


def write_srm(rest=False, sst=False):
    """The duct of p0.dat with rho nu~ = 3 MU (mi_tilde) for the Spalart-Allmaras model, and record 1 (the inlet) an SRM
    grain boundary (502): T_af, a, n, p_ref, rho_grain, SF_geo and the mass fraction of the single species. rest: the gas
    at rest at T0 and PB; sst: k = 1 m2/s2 and omega = 1e3 1/s (kappa = rho k, omega = rho omega) for SST instead."""
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
        turb = f" {3.0 * MU:.15E}\n" * N
    with open("INPUT/ic.tec", "w") as f:
        f.write("".join(lines) + turb)
    with open("INPUT/bc.txt") as f:
        bc = f.readlines()
    bc[0] = bc[0].replace("     402\n", "     502\n")
    bc[1] = "3.0E+03, 1.0E-05, 0.4, 1.0E+06, 1.7E+03, 1.0, 1.0\n"
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(bc))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), by variable name."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = lines[0].split('=', 1)[1].replace('"', ' ').split()
    head = lines[1].replace(',', ' ').split()
    ni, nj, nk = (int(next(w for w in head if w.startswith(k + '=')).split('=')[1]) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1), {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return out


def srm_check(field, diag):
    """The change of rho nu~ over the step in the cell next to the grain, times sigma rho h^2 / m over dt, against
    cb2 m - 2 mu (the wall value nu~_w = 0 in the ghost cell); the cell next to it, where everything is uniform, changes
    by round-off."""
    f, d = read_zone(field), read_zone(diag)
    rho, m, h = PB / (RGAS * T0), 3.0 * MU, LX / N
    r = [(new - m) / t * SIGMA * rho * h * h / m for new, t in zip(f["rhoNuTilde"], d["dt"])]
    expect = CB2 * m - 2.0 * MU_GAS
    err = abs(r[0] - expect) / (2.0 * MU_GAS)
    print(f"    change of rho nu~ next to the grain, times sigma rho h^2/(m dt): {r[0]:.6e}, from the wall value nu~ = 0 "
          f"{expect:.6e} (an extrapolated ghost cell: {-2.0 * MU_GAS:.6e}); difference {err:.1e} of 2 mu; second cell "
          f"{r[1]:.1e}")
    sys.exit(0 if err < 1e-6 and abs(r[1]) < 1e-6 * 2.0 * MU_GAS else 1)


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
    elif what == "srm":
        write_srm("rest" in sys.argv[2:], "sst" in sys.argv[2:])
    elif what == "srmcheck":
        srm_check(sys.argv[2], sys.argv[3])
    elif what == "same":
        a, b = cell_values(sys.argv[2]), cell_values(sys.argv[3])
        sys.exit(0 if a == b and len(a) > 0 else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
