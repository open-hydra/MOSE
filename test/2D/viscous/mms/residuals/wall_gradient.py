"""Part 5 of MMSResiduals (run.sh): Spalart-Allmaras next to a wall, inputs and check.

Spalart-Allmaras next to an isothermal wall, fluid at rest. A column of 4 x 16 cells, uniform, dx = dz = 0.1 mm and
dy = 10 um; periodic in x (connections 201 on itself), isothermal wall at 300 K (302) at y = 0, symmetry (300) at the top,
null k faces. Air of test/common/Air-viscous (mu and k constant) at rest at 300 K and 1 bar, with nu~ growing linearly from
the wall:
  nu~ = a y,   a = 6.6e-5 1/s,   so that chi = nu~/nu stays below 1e-3.
With the fluid at rest S~ = nu~ fv2/(kappa^2 d^2) = a fv2/(kappa^2 y) and r = 1/fv2, and for chi -> 0 (fv2 -> 1, fw -> 1) the
terms of the Spalart-Allmaras equation per unit volume are
  production cb1 a^2/kappa^2,  destruction cw1 a^2,  diffusion (1 + cb2) a^2 / sigma  (times rho),
whose sum is zero because cw1 = cb1/kappa^2 + (1 + cb2)/sigma: the linear profile is an equilibrium of the equation, up to
terms of order chi. That holds in the cell next to the wall too, where the wall value nu~ = 0 and the linear profile give
grad nu~ = a: the ghost cell of the wall is -nu~ of the cell (NASA CFL3D user's manual, 1998, eq. H-186), and the central
difference (nu~_2 - nu~_ghost)/(2 dy) is a. A ghost that copies the cell gives a/2 there instead, and leaves a right-hand side
of -(3/4) (cb2/sigma) rho a^2 in the wall cell.

The solver makes one explicit Euler step with local time stepping: the change of rho nu~ in OUTPUT/field.tec over the step of
each cell (dt of OUTPUT/diagnostic.tec) is the right-hand side the solver evaluated. The residual column of the diagnostic is
not used: in the wall cells it also holds the wall flux that the wall output adds after the step. On every cell of the 14
rows next to the wall (the top rows feel the symmetry plane) the right-hand side must stay below 0.05 of
(cb2/sigma) rho a^2.

Usage (from the test directory):
  python3 -B wall_gradient.py case                 INPUT/ic.tec and INPUT/bc.txt
  python3 -B wall_gradient.py check FIELD DIAG     exit status 0 when the right-hand side passes
"""
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P0 = 300.0, 1.0e5
SIGMA, CB2 = 2.0 / 3.0, 0.622       # Spalart-Allmaras
NI, NJ, DX, DY = 4, 16, 1.0e-4, 1.0e-5
A = 6.6e-5
ROWS = 14


def write_case():
    rho = P0 / (RGAS * T0)
    x = [i * DX for i in range(NI + 1)]
    y = [j * DY for j in range(NJ + 1)]
    cells = [((i + 0.5) * DX, (j + 0.5) * DY) for j in range(NJ) for i in range(NI)]
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={NI + 1}, J={NJ + 1}, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(2):
                for j in range(NJ + 1):
                    for i in range(NI + 1):
                        f.write(f" {(x[i], y[j], k * DX)[c]:.16E}\n")
        n = NI * NJ
        zero = [0.0] * n
        for col in ([rho] * n, zero, zero, zero, [P0] * n, [rho * A * yc for _, yc in cells]):
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, 1, face, t)

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    out = []
    for j in range(1, NJ + 1):
        out += [rec(1, j, 1, 201), con(1, NI, j, 1, 2, 1, 0, 0, 1)]
    for j in range(1, NJ + 1):
        out += [rec(NI, j, 2, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    for i in range(1, NI + 1):
        out += [rec(i, 1, 3, 302), f"    {T0:.6E}, 0.000000E+00, 0.000000E+00,\n"]
    for i in range(1, NI + 1):
        out.append(rec(i, NJ, 4, 300))
    for face in (5, 6):
        out += [rec(i, j, face, 0) for j in range(1, NJ + 1) for i in range(1, NI + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not."""
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
    return out


def check(field, diag):
    f, d = read_zone(field), read_zone(diag)
    rho = P0 / (RGAS * T0)
    scale = CB2 / SIGMA * rho * A * A
    start = [rho * A * (j + 0.5) * DY for j in range(NJ) for i in range(NI)]
    rhs = [(new - old) / t / scale for new, old, t in zip(f["rhoNuTilde"], start, d["dt"])]
    wall = max(abs(v) for v in rhs[:NI])
    rest = max(abs(v) for v in rhs[NI:ROWS * NI])
    print(f"    right-hand side of rho nu~ over (cb2/sigma) rho a^2: wall cells {wall:.3e}, rows 2 to {ROWS} {rest:.3e} "
          f"(a ghost cell that copies the wall cell gives 0.75 there)")
    sys.exit(0 if max(wall, rest) < 0.05 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case()
    elif what == "check":
        check(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
