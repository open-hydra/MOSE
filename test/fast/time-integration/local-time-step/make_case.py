"""Inputs of the LocalTimeStep fast test, and a reader of its outputs.

Small cases, written into INPUT/ (one at a time) by `python3 -B make_case.py CASE`:
  box    closed 2-D box, 32 x 32 square cells of 10 mm, slip walls (300) on the i and j faces, null (0) k faces;
         air at rest with a Gaussian pressure pulse (1e-2 p0, radius 50 mm). Sound waves cross the square cells in
         both directions at once.
  skew   2-D single block periodic in x and y (201 on itself), 24 x 24 cells on a stretched and skewed mesh
         x = s(xi) + A sin(2 pi eta), y = s(eta) + A sin(2 pi xi), s(q) = q + b sin(2 pi q)/(2 pi), b = 0.6,
         A = 0.05 (1 m x 1 m); flow at uniform pressure and velocity (u = 100, v = 50 m/s) carrying an entropy
         wave of amplitude 1e-3 in density along xi, which keeps the run from stopping on a zero residual after
         one step. p, u and v must stay uniform.
  visc   1-D periodic Navier-Stokes, 40 cells of 1 nm, air of test/common/Air-viscous (mu, k constant, Pr = 0.72)
         at rest with an entropy wave of amplitude 1e-3 in temperature: the cell Peclet number (|u|+c) dx / nu is
         0.034, so the step is set by diffusion, and heat conduction (diffusivity k/(rho cv) = gamma nu / Pr) is
         the fastest process. The wave must decay.
  sod    1-D shock tube of test/1D/Sod79 (100 cells), run at cfl 0.99.
  rhombo 3-D single block periodic in its three directions (201 on itself), 8 x 8 x 8 rhombohedra of edge 1 nm, every
         edge angle 115 degrees (edges e1 = (1, 0, 0), e2 = (cos a, sin a, 0), e3 = (cos a, (cos a - cos^2 a)/sin a, .)),
         air of test/common/Air-viscous at rest with a temperature wave of amplitude 1e-3 and wavenumber 135 degrees per
         cell along the three indices: the mode whose decay rate exceeds 4 D sum_d |grad xi_d|^2 by 1.0365 on these cells.
         Heat conduction sets the step; with the step of 4 D sum_d |grad xi_d|^2 alone, Heun at vnn 0.5 multiplies the
         mode by 1.076 per step. The wave must decay.
  box3d  closed 3-D box, 12 x 12 x 12 cubic cells of 10 mm, slip walls (300) on the six faces; air at rest with a
         Gaussian pressure pulse (1e-2 p0, radius 30 mm): sound waves cross the cells in the three directions at once.
  species  1-D periodic Navier-Stokes, 40 cells of 1 nm, the four-species gas of test/1D/diff1 (CH4, H2O, N2, O2) as
         air (N2 and O2) at rest, 300 K, 1 bar, with a wave of amplitude 1e-3 in the O2 mass fraction; run with
         schmidt = 0.2, so the species diffusivity mu/(rho Sc) = 5 nu is the largest one and sets the step.
Calorically perfect air of test/common/Air (box, skew, sod, box3d) or test/common/Air-viscous (visc): cp = 1004.5
J/kg/K, M = 28.970418 kg/kmol, gamma = 1.4; the gas constant of FLINT.

Usage (from the test directory):
  python3 -B make_case.py box|skew|visc|box3d|species|rhombo   INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py check CASE FILE               compare a solution file with what the case must keep
"""
import math
import os
import re
import sys

RGAS = 8314.46261815324 / 28.970418
P0, T0 = 1.0e5, 300.0


def rec(b, i, j, k, f, t, prop=None):
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, k, f, t)
    if prop is not None:
        s += '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % prop
    return s


def fresh(path):
    """Remove path first, so that a link left by another case is never written through."""
    if os.path.lexists(path):
        os.remove(path)
    return open(path, "w")


def write_ic(ni, nj, x, y, z, cells, path="INPUT/ic.tec", nk=1, species=("rho1",)):
    """x, y, z(i, j, k) node functions; cells: dict of name -> list of ni*nj*nk cell values (i fastest)."""
    names = list(species) + ["u", "v", "w", "p"]
    with fresh(path) as f:
        f.write(' VARIABLES ="x" "y" "z"  ' + " ".join(f'"{n}"' for n in names) + "\n")
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-{3 + len(names)}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for fun in (x, y, z):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {fun(i, j, k):.16E}\n")
        for name in names:
            f.write("".join(f" {v:.16E}\n" for v in cells[name]))


def box():
    n, h = 32, 0.01
    rho0 = P0 / (RGAS * T0)
    gam = 1.4
    rho, p = [], []
    for j in range(n):
        for i in range(n):
            r2 = (((i + 0.5) * h - 0.12) ** 2 + ((j + 0.5) * h - 0.14) ** 2) / 0.05 ** 2
            pc = P0 * (1.0 + 1.0e-2 * math.exp(-r2))
            p.append(pc)
            rho.append(rho0 * (pc / P0) ** (1.0 / gam))
    zero = [0.0] * (n * n)
    write_ic(n, n, lambda i, j, k: i * h, lambda i, j, k: j * h, lambda i, j, k: k * h,
             {"rho1": rho, "u": zero, "v": zero, "w": zero, "p": p})
    out = []
    for face, i in ((1, 1), (2, n)):
        out += [rec(1, i, j, 1, face, 300) for j in range(1, n + 1)]
    for face, j in ((3, 1), (4, n)):
        out += [rec(1, i, j, 1, face, 300) for i in range(1, n + 1)]
    for face in (5, 6):
        out += [rec(1, i, j, 1, face, 0) for j in range(1, n + 1) for i in range(1, n + 1)]
    fresh("INPUT/bc.txt").write("".join(out))


def skew():
    n, b, A = 24, 0.6, 0.05

    def s(q):
        return q + b * math.sin(2 * math.pi * q) / (2 * math.pi)
    rho0 = P0 / (RGAS * T0)
    rho = [rho0 * (1.0 + 1.0e-3 * math.sin(2 * math.pi * (i + 0.5) / n)) for j in range(n) for i in range(n)]
    write_ic(n, n, lambda i, j, k: s(i / n) + A * math.sin(2 * math.pi * j / n),
             lambda i, j, k: s(j / n) + A * math.sin(2 * math.pi * i / n),
             lambda i, j, k: (-0.0005, 0.0005)[k],
             {"rho1": rho, "u": [100.0] * (n * n), "v": [50.0] * (n * n), "w": [0.0] * (n * n),
              "p": [P0] * (n * n)})
    out = []
    out += [rec(1, 1, j, 1, 1, 201, (1, n, j, 1, 2, 1, 0, 0, 1)) for j in range(1, n + 1)]
    out += [rec(1, n, j, 1, 2, 201, (1, 1, j, 1, 1, 1, 0, 0, 1)) for j in range(1, n + 1)]
    out += [rec(1, i, 1, 1, 3, 201, (1, i, n, 1, 4, 1, 0, 0, 1)) for i in range(1, n + 1)]
    out += [rec(1, i, n, 1, 4, 201, (1, i, 1, 1, 3, 1, 0, 0, 1)) for i in range(1, n + 1)]
    for face in (5, 6):
        out += [rec(1, i, j, 1, face, 0) for j in range(1, n + 1) for i in range(1, n + 1)]
    fresh("INPUT/bc.txt").write("".join(out))


def visc():
    n, h, amp = 40, 1.0e-9, 1.0e-3
    T = [T0 * (1.0 + amp * math.sin(2 * math.pi * (i + 0.5) / n)) for i in range(n)]
    zero = [0.0] * n
    write_ic(n, 1, lambda i, j, k: i * h, lambda i, j, k: j * h, lambda i, j, k: k * h,
             {"rho1": [P0 / (RGAS * t) for t in T], "u": zero, "v": zero, "w": zero, "p": [P0] * n})
    out = [rec(1, 1, 1, 1, 1, 201, (1, n, 1, 1, 2, 1, 0, 0, 1)), rec(1, n, 1, 1, 2, 201, (1, 1, 1, 1, 1, 1, 0, 0, 1))]
    for face in (3, 4, 5, 6):
        out += [rec(1, i, 1, 1, face, 0) for i in range(1, n + 1)]
    fresh("INPUT/bc.txt").write("".join(out))


def box3d():
    n, h = 12, 0.01
    rho0 = P0 / (RGAS * T0)
    gam = 1.4
    rho, p = [], []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                r2 = (((i + 0.5) * h - 0.05) ** 2 + ((j + 0.5) * h - 0.06) ** 2 + ((k + 0.5) * h - 0.055) ** 2) / 0.03 ** 2
                pc = P0 * (1.0 + 1.0e-2 * math.exp(-r2))
                p.append(pc)
                rho.append(rho0 * (pc / P0) ** (1.0 / gam))
    zero = [0.0] * n ** 3
    write_ic(n, n, lambda i, j, k: i * h, lambda i, j, k: j * h, lambda i, j, k: k * h,
             {"rho1": rho, "u": zero, "v": zero, "w": zero, "p": p}, nk=n)
    out = []
    for face, i in ((1, 1), (2, n)):
        out += [rec(1, i, j, k, face, 300) for k in range(1, n + 1) for j in range(1, n + 1)]
    for face, j in ((3, 1), (4, n)):
        out += [rec(1, i, j, k, face, 300) for k in range(1, n + 1) for i in range(1, n + 1)]
    for face, k in ((5, 1), (6, n)):
        out += [rec(1, i, j, k, face, 300) for j in range(1, n + 1) for i in range(1, n + 1)]
    fresh("INPUT/bc.txt").write("".join(out))


# Species of test/1D/diff1/INPUT/phase.txt, in this order: CH4, H2O, N2 (28.014 kg/kmol), O2 (31.998 kg/kmol)
Y_O2 = 0.233


def species():
    n, h, amp = 40, 1.0e-9, 1.0e-3
    yo = [Y_O2 * (1.0 + amp * math.sin(2 * math.pi * (i + 0.5) / n)) for i in range(n)]
    rho = [P0 / (8314.46261815324 * (y / 31.998 + (1.0 - y) / 28.014) * T0) for y in yo]
    zero = [0.0] * n
    write_ic(n, 1, lambda i, j, k: i * h, lambda i, j, k: j * h, lambda i, j, k: k * h,
             {"rho1": zero, "rho2": zero, "rho3": [r * (1.0 - y) for r, y in zip(rho, yo)],
              "rho4": [r * y for r, y in zip(rho, yo)], "u": zero, "v": zero, "w": zero, "p": [P0] * n},
             species=("rho1", "rho2", "rho3", "rho4"))
    out = [rec(1, 1, 1, 1, 1, 201, (1, n, 1, 1, 2, 1, 0, 0, 1)), rec(1, n, 1, 1, 2, 201, (1, 1, 1, 1, 1, 1, 0, 0, 1))]
    for face in (3, 4, 5, 6):
        out += [rec(1, i, 1, 1, face, 0) for i in range(1, n + 1)]
    fresh("INPUT/bc.txt").write("".join(out))


def rhombo():
    n, h, amp, alpha = 8, 1.0e-9, 1.0e-3, math.radians(115.0)
    ca, sa = math.cos(alpha), math.sin(alpha)
    e3y = (ca - ca * ca) / sa
    edges = ((1.0, 0.0, 0.0), (ca, sa, 0.0), (ca, e3y, math.sqrt(1.0 - ca * ca - e3y * e3y)))

    def node(c):
        return lambda i, j, k: h * (i * edges[0][c] + j * edges[1][c] + k * edges[2][c])
    theta = 2.0 * math.pi * 3.0 / 8.0
    T = [T0 * (1.0 + amp * math.cos(theta * (i + j + k))) for k in range(n) for j in range(n) for i in range(n)]
    zero = [0.0] * n ** 3
    write_ic(n, n, node(0), node(1), node(2),
             {"rho1": [P0 / (RGAS * t) for t in T], "u": zero, "v": zero, "w": zero, "p": [P0] * n ** 3}, nk=n)
    out = []
    for k in range(1, n + 1):
        for j in range(1, n + 1):
            out += [rec(1, 1, j, k, 1, 201, (1, n, j, k, 2, 1, 0, 0, 1)), rec(1, n, j, k, 2, 201, (1, 1, j, k, 1, 1, 0, 0, 1))]
    for k in range(1, n + 1):
        for i in range(1, n + 1):
            out += [rec(1, i, 1, k, 3, 201, (1, i, n, k, 4, 1, 0, 0, 1)), rec(1, i, n, k, 4, 201, (1, i, 1, k, 3, 1, 0, 0, 1))]
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            out += [rec(1, i, j, 1, 5, 201, (1, i, j, n, 6, 1, 0, 0, 1)), rec(1, i, j, n, 6, 201, (1, i, j, 1, 5, 1, 0, 0, 1))]
    fresh("INPUT/bc.txt").write("".join(out))


def read_field(path):
    """Cell values of a one-zone MOSE ASCII field (BLOCK packing), variable names quoted or not."""
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


def check(case, path):
    """Exit status 0 when the field is what the case must keep."""
    f = read_field(path)
    rhos = [name for name in f if name.startswith("rho(")]
    ok = all(math.isfinite(v) for name in rhos + ["u", "v", "p"] for v in f[name])
    if case in ("box", "box3d"):
        lo, hi = min(f["p"]) / P0, max(f["p"]) / P0
        print(f"    p/p0 in [{lo:.6f}, {hi:.6f}]")
        ok = ok and 0.98 < lo and hi < 1.02
    elif case == "skew":
        du = max(abs(u / 100.0 - 1) for u in f["u"])
        dv = max(abs(v / 50.0 - 1) for v in f["v"])
        dp = max(abs(p / P0 - 1) for p in f["p"])
        print(f"    largest relative deviation from the uniform state: u {du:.1e}, v {dv:.1e}, p {dp:.1e}")
        ok = ok and max(du, dv, dp) < 1.0e-10
    elif case in ("visc", "rhombo"):
        T = [p / (r * RGAS) for p, r in zip(f["p"], f["rho(1)"])]
        tm = sum(T) / len(T)
        amp = max(abs(t / tm - 1) for t in T)
        print(f"    amplitude of the temperature wave {amp:.3e} (initial 1.000e-03)")
        ok = ok and amp < 1.0e-3
    elif case == "species":
        rho = [sum(c) for c in zip(*(f[name] for name in rhos))]
        y = [r4 / r for r4, r in zip(f["rho(4)"], rho)]
        ym = sum(y) / len(y)
        amp = max(abs(v / ym - 1) for v in y)
        print(f"    amplitude of the O2 mass fraction wave {amp:.3e} (initial 1.000e-03)")
        ok = ok and amp < 1.0e-3
    elif case == "sod":
        lo = min(min(f["p"]), min(f["rho(1)"]))
        print(f"    smallest pressure or density {lo:.4e}")
        ok = ok and lo > 0.0
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "check":
        check(sys.argv[2], sys.argv[3])
    elif what in ("box", "skew", "visc", "box3d", "species", "rhombo"):
        globals()[what]()
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
