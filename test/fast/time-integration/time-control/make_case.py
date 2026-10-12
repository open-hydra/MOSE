"""Inputs of the TimeControl fast test, and a reader of its outputs.

One-dimensional duct, 100 cells on 0 <= x <= 1 m (the mesh of test/1D/Sod79),
subsonic stagnation inlet at i = 1, back-pressure outlet (406) at i = 100,
slip walls (300) on the j faces, null (0) k faces, calorically perfect air of
test/common/Air (cp = 1004.5 J/kg/K, gamma = 1.4). The uniform initial state
is the isentropic expansion from (T0, P0) to PB, i.e. the steady state of the
duct for constant inlet data.

Part 4 runs a tube along x instead: NT cells of DXT, one cell of DXT in y and in
z; extrapolation (400) on the i faces, slip walls (300) on the j faces, null (0)
k faces. Air in a uniform flow at Mach 2 along x, at uniform pressure, with an
entropy wave: the temperature is 20 % higher in a bell around XT0, so the
density is lower there. The pressure is PT, or PT times 2^-24: with an even
power of 2 the scaling of density and pressure is exact in floating point,
products and quotients by it and the square roots of the density in the Roe
averages of the HLLC wave speeds included (sqrt(c rho) = sqrt(c) sqrt(rho),
sqrt(c) a power of 2 as well; an odd power of 2 would round them). The values
are written with 18 significant digits, so the solver reads the same doubles.

Usage (from the test directory):
  python3 -B make_case.py uniform              INPUT/ic.tec
  python3 -B make_case.py bc401                INPUT/bc.txt, inlet 401 with p0 = P0
  python3 -B make_case.py bc402 TP EPS [T1]    INPUT/bc.txt, inlet 402 with
          p0(t) = P0 (1 + EPS sin(2 pi (t - T1) / TP)), tabulated in p0.dat (bar)
          for T1 <= t <= T1 + 4 TP (T1 = 0 if omitted)
  python3 -B make_case.py dt0 CFL              time step of the uniform state
  python3 -B make_case.py tube bar|tiny        INPUT/ic.tec and INPUT/bc.txt of the
          tube of part 4, at p = PT or PT x 2^-24
"""
import math
import os
import re
import sys

N = 100              # cells along x
LX = 1.0             # duct length [m]
LY = 0.01            # j thickness [m]
LZ = 1.0e-3          # k thickness [m]
T0 = 300.0           # inlet total temperature [K]
P0 = 1.2e5           # inlet total pressure [Pa]
PB = 1.0e5           # outlet back pressure [Pa]
GAMMA = 1.4
RGAS = 8314.46261815324 / 28.970418   # universal gas constant of FLINT / molar mass of air
NT, DXT = 40, 0.01           # tube of part 4: cells along x, cell size [m]
XT0, WT = 0.1, 0.03          # centre and width of its entropy wave [m]
TT, PT = 300.0, 1.0e5        # temperature away from the wave [K], pressure [Pa]
MT = 2.0                     # Mach number of the tube
TUBE_SCALE = {"bar": 1.0, "tiny": 2.0 ** -24}


def uniform_state():
    """Static density, velocity, pressure and sound speed of the uniform state."""
    mach = math.sqrt(2.0 / (GAMMA - 1.0) * ((P0 / PB) ** ((GAMMA - 1.0) / GAMMA) - 1.0))
    t = T0 / (1.0 + 0.5 * (GAMMA - 1.0) * mach ** 2)
    c = math.sqrt(GAMMA * RGAS * t)
    return PB / (RGAS * t), mach * c, PB, c


def write_ic(path="INPUT/ic.tec"):
    rho, u, p, _ = uniform_state()
    with open(path, "w") as f:
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


def write_bc(inlet_type, inlet_line):
    with open("INPUT/bc.txt", "w") as f:
        f.write(f"       1       1       1       1       1     {inlet_type}\n{inlet_line}\n")
        f.write(f"       1     {N:3d}       1       1       2     406\n{PB:.6E}, 1.0\n")
        for face, btype in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, N + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {btype:3d}\n")


def write_p0_table(period, eps, t1):
    # piecewise-linear table read by MOSE: 2000 points per period, values in bar
    nstep = 8000
    with open("p0.dat", "w") as f:
        for n in range(nstep + 1):
            t = t1 + 4.0 * period * n / nstep
            p0 = P0 * (1.0 + eps * math.sin(2.0 * math.pi * (t - t1) / period))
            f.write(f"{t:.15E} {p0 / 1.0e5:.15E}\n")


def write_tube(name):
    c = TUBE_SCALE[name]
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    u = MT * math.sqrt(GAMMA * RGAS * TT)
    rho = []
    for i in range(NT):
        x = DXT * (i + 0.5)
        t = TT * (1.0 + 0.2 * math.exp(-((x - XT0) / WT) ** 2))
        rho.append(PT / (RGAS * t) * c)          # exact scaling: c is a power of 2
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NT + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(NT + 1):
                        f.write(f" {(DXT * i, DXT * j, DXT * k)[axis]:.17E}\n")
        f.write("".join(f" {r:.17E}\n" for r in rho))
        for value in (u, 0.0, 0.0, PT * c):
            f.write(f" {value:.17E}\n" * NT)
    with open("INPUT/bc.txt", "w") as f:
        for face, i in ((1, 1), (2, NT)):
            f.write(f"       1     {i:3d}       1       1       {face}     400\n")
        for face, t in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, NT + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {t:3d}\n")


def read_field(path):
    """Solution time and cell values of a one-zone MOSE ASCII field (BLOCK packing)."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]  # quoted or not
    time = float(lines[1].split("SOLUTIONTIME=")[1].split(",")[0])
    values = [float(v) for v in " ".join(lines[2:]).split()]
    nn, out, pos = (N + 1) * 4, {}, 0
    for k, name in enumerate(names):
        cnt = nn if k < 3 else N
        out[name] = values[pos:pos + cnt]
        pos += cnt
    out["time"] = time
    return out


def main():
    what = sys.argv[1]
    if what == "uniform":
        write_ic()
    elif what == "bc401":
        write_bc(401, f"{T0}, {P0:.6E}, normal, normal, 1.0, 1.0")
    elif what == "bc402":
        t1 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
        write_p0_table(float(sys.argv[2]), float(sys.argv[3]), t1)
        # file name without '/': a slash ends a list-directed read
        write_bc(402, f"{T0}, p0.dat, normal, normal, 1.0, 1.0")
    elif what == "tube":
        write_tube(sys.argv[2])
    elif what == "dt0":
        _, u, _, c = uniform_state()
        print(f"{float(sys.argv[2]) * (LX / N) / (u + c):.15e}")
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
