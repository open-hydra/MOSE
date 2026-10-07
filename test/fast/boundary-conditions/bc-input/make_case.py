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


def write_case(name):
    rho, u, p = uniform_state()
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
}


def write_bad(kind):
    write_case("p0.dat")
    n, text = BAD[kind]
    with open("INPUT/bc.txt") as f:
        lines = f.readlines()
    lines[n] = text
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(lines))


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
    elif what == "same":
        a, b = cell_values(sys.argv[2]), cell_values(sys.argv[3])
        sys.exit(0 if a == b and len(a) > 0 else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
