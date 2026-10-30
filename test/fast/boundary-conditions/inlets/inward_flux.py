"""Part 3 of Inlets (run.sh): inlets with a prescribed mass flux or normal velocity, inputs and a reader of the outputs.

The mass flux g of the inlets 403 (with the total temperature T0), 404 (with the static temperature T) and 420 (the choked
nozzle) and the velocity un of the inlet 408 (with T) are positive into the domain on every face
(docs/user/boundary-conditions.md, common fields): on the faces 2, 4, 6, where the face normal of the solver points out of
the domain, a negative value is taken as its magnitude.

A 1-D tube along x: NX cells of DX, one cell of DX in y and in z; slip walls (300) on the j faces, null (0) k faces, the
inlet at one end and an extrapolation face (400) at the other. Calorically perfect air of test/common/Air (cp = 1004.5
J/(kg K), molar mass 28.970418) at p = 1 bar and T = 300 K, its gas moving in through the inlet face at 0.1 a (u = 0.1 a
with the inlet at i = 1, -0.1 a with the inlet at i = NX), so that the inlet lets the gas in. One forward Euler step with
first-order reconstruction: the change of the density of the boundary cell is the mass flux of its boundary face times
DT / DX, plus what the interior face takes away.

The face of an inlet 403 or 404 carries the mass flux of its record, |g|; that of an inlet 408 carries rho3 |un|, with
rho3 = p / (R T) at the pressure p of the cell and the T of the record (here the temperature of the cell); that of a nozzle
420 whose cell pressure lies between p_shock and p_sub carries |g| too: within 1e-9 of rho a, the rounding of the step.
The nozzle is fed from T0 = T and p0 = P0_420 = 2 p, and p_sub and p_sup are those of the nozzle that this reservoir
chokes at g (nozzle_thresholds below), so that the cell pressure p lies between them, and above p_shock, the pressure
behind a normal shock at the exit of that nozzle (0.61 p).

Cases (CASES below), g = 0.3 rho a and un = 0.1 a:
  imax403, imax404, imax408           the inlet on face 2 (i = NX), the value positive
  imaxneg403, imaxneg404, imaxneg408  the same with the value negative: the same field as the positive one, bit for bit
                                      (run.sh)
  imin403, imin404, imin408           the inlet on face 1 (i = 1), the value positive
  imax420                             a nozzle 420 on face 2 with g positive, the cell in the band of g (p_shock <= p < p_sub)

Usage (from the test directory):
  python3 -B inward_flux.py case NAME          INPUT/ic.tec and INPUT/bc.txt of case NAME
  python3 -B inward_flux.py check NAME FILE    mass flux into the domain through the inlet face in FILE (a field.tec)
                                               against that of the record; exit status 0 when they agree within 1e-9
                                               of rho a
"""
import math
import os
import re
import sys

NX = 8
DX = 0.01
P, T = 1.0e5, 300.0
CP = 1004.5
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air
GAMMA = CP / (CP - RGAS)             # as the thermodynamics of the solver gives it for this gas
RHO = P / (RGAS * T)
A = math.sqrt(GAMMA * RGAS * T)
G_IN = 0.3 * RHO * A                 # the mass flux of the records 403, 404 and 420
UN_IN = 0.1 * A                      # the velocity of the records 408
W = 0.1 * A                          # the speed of the gas of the cell into the domain
P0_420 = 2.0 * P                     # the total pressure of the nozzle 420 (total temperature T)
TOL = 1.0e-9

# name: (type, face of the inlet, sign of the value)
CASES = {
    "imax403": (403, 2, 1.0), "imax404": (404, 2, 1.0), "imax408": (408, 2, 1.0),
    "imaxneg403": (403, 2, -1.0), "imaxneg404": (404, 2, -1.0), "imaxneg408": (408, 2, -1.0),
    "imin403": (403, 1, 1.0), "imin404": (404, 1, 1.0), "imin408": (408, 1, 1.0),
    "imax420": (420, 2, 1.0),
}


def nozzle_thresholds(g):
    """p_sub and p_sup of the nozzle that the reservoir (T, P0_420) chokes at the mass flux g, calorically perfect (cp = CP):
    exit-to-throat area ratio G*/g, G* the choked mass flux, and the exit pressures of the subsonic and of the supersonic
    isentropic flows. The exit Mach numbers by bisection on the two branches of the area ratio A/A*(M), decreasing for M < 1
    and increasing for M > 1, from brackets that hold the root for 1 < G*/g < 1e5 (200 halvings: the bracket shrinks to the
    spacing of the doubles)."""
    e = (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0))
    ae_at = P0_420 * math.sqrt(GAMMA / (RGAS * T)) * (2.0 / (GAMMA + 1.0)) ** e / g

    def area(m):
        return ((2.0 / (GAMMA + 1.0)) * (1.0 + 0.5 * (GAMMA - 1.0) * m * m)) ** e / m

    def root(lo, hi, increasing):
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (area(mid) < ae_at) == increasing:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def pressure(m):
        return P0_420 * (1.0 + 0.5 * (GAMMA - 1.0) * m * m) ** (-GAMMA / (GAMMA - 1.0))

    return pressure(root(1.0e-6, 1.0, False)), pressure(root(1.0, 50.0, True))


def cell_velocity(name):
    """The velocity along x of the gas of the tube: into the domain through the face of the inlet."""
    return W if CASES[name][1] == 1 else -W


def write_case(name):
    btype, face, sign = CASES[name]
    value = sign * (UN_IN if btype == 408 else G_IN)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(NX + 1):
                        f.write(f" {(DX * i, DX * j, DX * k)[axis]:.15E}\n")
        for v in (RHO, cell_velocity(name), 0.0, 0.0, P):
            f.write(f" {v:.15E}\n" * NX)
    types = (btype, 400) if face == 1 else (400, btype)
    with open("INPUT/bc.txt", "w") as f:
        for fc, i, t in ((1, 1, types[0]), (2, NX, types[1])):
            f.write(f"       1     {i:3d}       1       1       {fc}     {t}\n")
            if t == 420:                                 # T0, p0, p_sub, p_sup, g, rel_fac, mass fraction
                p_sub, p_sup = nozzle_thresholds(G_IN)
                f.write(f"    {T:.15E}, {P0_420:.15E}, {p_sub:.15E}, {p_sup:.15E}, {value:.15E}, 1.0, 1.0,\n")
            elif t != 400:                               # T0 (403) or T (404, 408), then g or un
                f.write(f"    {T:.15E}, {value:.15E}, normal, normal, 1.0, 1.0,\n")
        for fc, t in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}       1       1       {fc}     {t:3d}\n")


def read_field(path):
    """Solution time and cell values of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    time = float(re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)', lines[1]).group(1).replace('D', 'E'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, out, pos = (NX + 1) * 2 * 2, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else NX
        out[name] = values[pos:pos + cnt]
        pos += cnt
    return time, out


def check(name, path):
    btype, face, sign = CASES[name]
    dt, field = read_field(path)
    rho = field.get("rho(1)", field.get("rho1"))
    u0 = cell_velocity(name)
    if face == 1:              # the interior face takes RHO u0 out of cell 1
        got = (rho[0] - RHO) * DX / dt + RHO * u0
    else:                      # the interior face brings RHO u0 into cell NX
        got = (rho[-1] - RHO) * DX / dt - RHO * u0
    want = RHO * UN_IN if btype == 408 else G_IN
    unit = RHO * A
    ok = abs(got - want) <= TOL * unit
    print(f"    face {face}: mass flux in {got / unit:.10f} rho a, record {want / unit:.10f} rho a"
          f"{'' if ok else '  <-- differs'}")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check":
        check(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
