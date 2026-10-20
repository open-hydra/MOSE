"""Part Extrapolation400 of BoundaryConditions (run.sh): the extrapolation face (400) and the lines of the gas entering through it;
inputs, model, checks.

A 400 face takes the flux of the state of its boundary cell, F(q1), whatever the direction of the gas. With w1 the velocity of the
gas of the boundary cell along the normal out of the domain and a1 its speed of sound, where w1 < 0 the gas enters with the state
of the cell, and the run says so. The first time the gas enters through the faces of a 400 boundary (the 400 faces of one face of
a block, on each grid level) it prints the line
   400 boundary of block <b>, face <f>, gas entering at iteration <n>, faces: <m>, largest |w|/a: <x>
with the faces of the boundary where the gas enters at that stage and their largest -w1/a1. Each time the largest -w1/a1 of the
boundary at a stage is more than 10 times that of its last line, one more line follows, which counts the faces above that:
   400 boundary of block <b>, face <f>, gas entering at iteration <n>, faces with |w|/a above 10 times its last line: <m>,
   largest |w|/a: <x>
(one line in the output). At its end a run whose case has 400 faces prints the flux evaluations where the gas entered, over the
stages, and their largest -w1/a1:
   400 faces with gas entering, flux evaluations: <n>, largest |w|/a: <x>

Calorically perfect air of test/common/Air (cp = 1004.5 J/(kg K), molar mass 28.970418), p = 1 bar, T = 300 K; the total energy
per volume is rho E = p/(gamma-1) + rho |q|^2/2 (h = cp T in its table).

Part 1, one step (ONE_STEP): a 1-D tube along x of NX cells of DX, one cell of DX in y and in z, slip walls (300) on the j faces,
null (0) k faces, 400 on the i faces; a uniform state at the velocity (u0, 0, w0); one forward Euler step at first order. The j
faces carry the pressure only (v = 0) and cancel, the k faces carry nothing, the interior faces carry F(q0): with the flux of the
cell state on the 400 faces, for either sign of u0, the boundary cells do not change (every conserved quantity within TOL_CELL of
rho, rho a, rho a^2) and the flux of each 400 face out of the domain, from the change of its cell, is F(q0) . S within TOL of
rho a, rho a^2, rho a^3 (S its outward area vector). Lines: one for the face where the gas enters (face 1 for u0 > 0, face 2 for
u0 < 0), iteration 1, faces 1, x = |u0|/a; none at rest; the end line with 1 evaluation (0 at rest) and the same x.
Part 2, the end count (COUNT): the tube of part 1 with the gas entering through face 1 at 0.5 a, 20 RK3 steps at first order and
cfl 0.5: the state stays uniform, one line only (x does not grow) and the end line with 60 evaluations (3 stages a step) and x 0.5.
Part 3, the lines at each factor of 10 (TENFOLD): the tube with the gas at 1e-3 a and the pressure 10 % lower in its right half,
40 forward Euler steps at first order and cfl 0.5: the expansion speeds up the gas that enters through face 1. At least two lines
of face 1, the first at iteration 1 with x = 1e-3, each later one with the label of the later lines, at a later iteration and with
x more than ten times that of the line before; the largest x of the end line at least that of the last line and less than ten
times it; no line of another face.
Part 4, a run that stops (STOP): the gas at u = 6 a along x moves away from face 1 of type 0 (no flux), the j faces 400 with the gas
entering through face 3 at v = 1e-3 a across, one forward Euler step at first order and cfl 1.13. The boundary cell keeps only the
flux F(q0) of its face i = 2 (the j faces carry opposite fluxes of the same state, which cancel); with lambda = u dt/DX it holds
   rho' = rho (1 - lambda), (rho u)' = rho u (1 - lambda) - p lambda/u, (rho v)' = rho v (1 - lambda), (rho E)' = rho E (1 - lambda)
   - lambda p,
so v' = v, and rho' v'^2/2 = (1 - lambda) rho v^2/2, the part of (rho E)' that holds the energy of the velocity across: the two
cancel in p' = (gamma-1) ((rho E)' - (rho u)'^2/(2 rho') - rho' v'^2/2) = p (1 - lambda) - (gamma-1) p^2 lambda^2/(2 rho u^2
(1 - lambda)), which is negative for lambda > 1/(1 + sqrt((gamma-1)/(2 gamma))/M), M = u/a. The time step of a 1-D case takes the
direction along x only, dt = cfl DX/(|u| + a), so lambda = cfl M/(M + 1): at M = 6 the pressure is negative for cfl > 1.0975 and the
density positive for cfl < 7/6. At cfl 1.13 the run stops at its first stage (Integration failed); the line of face 3 (faces NX,
x 1e-3) comes before the line of the stop, and there is no end line.
Part 5, threads (THREADS, run with 1 and with 4 threads): a tube of NT = 130 cells, all j faces 400, the gas at 0.5 a along x and
1e-3 a across, one step: 130 records on face 3, more than two pieces of the dynamic schedule of 64 of the loops of the boundary
fluxes. Both runs give the same lines and the same end line: face 1 (faces 1, x 0.5) and face 3 (faces 130, x 1e-3), 131
evaluations.
Part 6, no 400 face (NO400): the tube of part 1 with slip walls (300) on its i faces too, the gas at 0.5 a along x, one forward
Euler step: no line of the gas entering and no end line.

Usage (from the test directory):
  python3 -B entering_gas.py case NAME              INPUT/ic.tec and INPUT/bc.txt of case NAME
  python3 -B entering_gas.py check NAME FIELD LOG   the checks of case NAME on FIELD (a field.tec) and LOG (the output of the run)
  python3 -B entering_gas.py stop NAME LOG RC       the checks of the run of case NAME that stops, on LOG and its exit status RC
  python3 -B entering_gas.py same LOG1 LOG4         part 5: the lines of the two runs, the same and as expected
Exit status 0 when every check passes.
"""
import math
import os
import re
import sys

NX = 8
NT = 130
DX = 0.01
P, T = 1.0e5, 300.0
CP = 1004.5
RGAS = 8314.46261815324 / 28.970418           # universal gas constant of FLINT / molar mass of air
GAMMA = CP / (CP - RGAS)                      # as the thermodynamics of the solver gives it for this gas
RHO = P / (RGAS * T)
A = math.sqrt(GAMMA * RGAS * T)
TOL = 1.0e-9
TOL_CELL = 1.0e-12
TOL_LINE = 1.0e-5                             # the lines print six digits of -w1/a1

ONE_STEP = {                                  # name: (u0 / a, w0 / a)
    "rest": (0.0, 0.0),
    "sub": (0.5, 0.0),
    "subm": (-0.5, 0.0),
    "sup": (1.327, 2.0),
    "supm": (-1.327, 2.0),
    "vac": (6.0, 0.0),
    "vacm": (-6.0, 0.0),
}
COUNT = {"count": 0.5}                        # name: u0 / a, 20 RK3 steps
NSTEP_COUNT, NSTAGE_COUNT = 20, 3
TENFOLD = {"tenfold": (1.0e-3, 0.9)}          # name: (u0 / a, pressure of the right half over P), 40 steps
STOP = {"stop": (6.0, 1.0e-3)}                # name: (u0 / a, v0 / a), the run stops at its first stage
THREADS = {"threads": (0.5, 1.0e-3)}          # name: (u0 / a, v0 / a), NT cells
NO400 = {"no400": 0.5}                        # name: u0 / a, no 400 face

FIRST = "faces"                               # the label of the count in the first line of a boundary, and in the later ones
LATER = "faces with |w|/a above 10 times its last line"
END = re.compile(r"400 faces with gas entering, flux evaluations:\s*(\d+), largest \|w\|/a:\s*(\S+)")
IN = re.compile(r"^ 400 boundary of block (\d+), face (\d+)(?:, grid level (\d+))?, gas entering at iteration (\d+), "
                r"(" + re.escape(FIRST) + "|" + re.escape(LATER) + r"): (\d+), largest \|w\|/a:\s*(\S+)$")


def total_energy(rho, vel, p):
    return p / (GAMMA - 1.0) + 0.5 * rho * sum(v * v for v in vel)


def euler_flux(rho, vel, p, s):
    """Flux of the state (rho, vel, p) through a face of area vector s, along s: (mass, momentum (3), energy)."""
    un = sum(v * c for v, c in zip(vel, s))
    h = (total_energy(rho, vel, p) + p) / rho
    return [rho * un] + [rho * un * v + p * c for v, c in zip(vel, s)] + [rho * un * h]


def ncells(name):
    return NT if name in THREADS else NX


def types(name):
    """Types of the six faces of the block of case NAME."""
    if name in STOP:
        return {1: 0, 2: 400, 3: 400, 4: 400, 5: 0, 6: 0}
    if name in THREADS:
        return {1: 400, 2: 400, 3: 400, 4: 400, 5: 0, 6: 0}
    if name in NO400:
        return {1: 300, 2: 300, 3: 300, 4: 300, 5: 0, 6: 0}
    return {1: 400, 2: 400, 3: 300, 4: 300, 5: 0, 6: 0}


def states(name):
    """rho, (u, v, w), p of every cell of case NAME at the start."""
    n = ncells(name)
    if name in ONE_STEP:
        u, w = ONE_STEP[name]
        return [(RHO, (u * A, 0.0, w * A), P)] * n
    if name in COUNT:
        return [(RHO, (COUNT[name] * A, 0.0, 0.0), P)] * n
    if name in TENFOLD:
        u, ratio = TENFOLD[name]
        return [(RHO, (u * A, 0.0, 0.0), P if i < n // 2 else ratio * P) for i in range(n)]
    if name in STOP:
        u, v = STOP[name]
        return [(RHO, (u * A, v * A, 0.0), P)] * n
    if name in NO400:
        return [(RHO, (NO400[name] * A, 0.0, 0.0), P)] * n
    u, v = THREADS[name]
    return [(RHO, (u * A, v * A, 0.0), P)] * n


def write_case(name):
    n = ncells(name)
    cells = states(name)
    for path in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(path):
            os.remove(path)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={n + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(n + 1):
                        f.write(f" {DX * (i, j, k)[axis]:.15E}\n")
        for q in range(5):
            for rho, vel, p in cells:
                f.write(f" {(rho, vel[0], vel[1], vel[2], p)[q]:.15E}\n")
    t = types(name)
    with open("INPUT/bc.txt", "w") as f:
        for face in range(1, 7):
            if face <= 2:
                f.write(f"       1     {1 if face == 1 else n:3d}       1       1       {face}     {t[face]:3d}\n")
            else:
                for i in range(1, n + 1):
                    f.write(f"       1     {i:3d}       1       1       {face}     {t[face]:3d}\n")


def read_field(path):
    """Solution time and values (nodal x, y, z, then cell values, i fastest) of a one-zone MOSE ASCII field."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    head = lines[1]
    ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', head).group(1)) for k in ('I', 'J', 'K'))
    time = float(re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)', head).group(1).replace('D', 'E'))
    values = [float(v.replace('D', 'E')) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, nc, out, pos = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1), {}, 0
    for q, name in enumerate(names):
        cnt = nn if q < 3 else nc
        out[name] = values[pos:pos + cnt]
        pos += cnt
    if pos != len(values):
        sys.exit(f"{path}: {len(values)} values, {pos} expected")
    return time, out


def cell(field, c):
    rho = field.get("rho(1)", field.get("rho1"))[c]
    return rho, (field["u"][c], field["v"][c], field["w"][c]), field["p"][c]


def conserved(rho, vel, p):
    return [rho] + [rho * v for v in vel] + [total_energy(rho, vel, p)]


def lines_in(log):
    """(block, face, grid level or 0, iteration, first, faces, x) of each line of the gas entering in LOG, in order; first is True
    when the line has the label of the first line of a boundary."""
    with open(log, errors="replace") as f:
        found = [IN.match(l.rstrip("\n")) for l in f]
    return [(int(m.group(1)), int(m.group(2)), int(m.group(3) or 0), int(m.group(4)), m.group(5) == FIRST, int(m.group(6)),
             float(m.group(7))) for m in found if m]


def end_line(log):
    """n and x of the end line of LOG, or None when it is missing or not there once."""
    with open(log, errors="replace") as f:
        found = [m for m in (END.search(l) for l in f) if m]
    return (int(found[0].group(1)), float(found[0].group(2))) if len(found) == 1 else None


def close(x, want):
    return x == 0.0 if want == 0.0 else abs(x - want) <= TOL_LINE * abs(want)


def show(g):
    return f"block {g[0]} face {g[1]} iteration {g[3]}, {FIRST if g[4] else LATER}: {g[5]}, {g[6]:.5e}"


def check_end(log, n_want, x_want):
    got = end_line(log)
    ok = got is not None and got[0] == n_want and close(got[1], x_want)
    shown = "missing" if got is None else f"{got[0]} evaluations, largest -w1/a1 {got[1]:.5e}"
    print(f"    end line: {shown} (expected {n_want}, {x_want:.5e}){'' if ok else '  <-- differs'}")
    return ok


def check_lines(log, want):
    """The lines of the gas entering in LOG: exactly one for each (face, faces, x) of WANT, in that order, all first lines of block
    1 at iteration 1, with no grid level, x within TOL_LINE."""
    got = lines_in(log)
    ok = len(got) == len(want) and all(g[0] == 1 and g[1] == f and g[2] == 0 and g[3] == 1 and g[4] and g[5] == n and close(g[6], x)
                                       for g, (f, n, x) in zip(got, want))
    shown = "; ".join(show(g) for g in got) or "none"
    wanted = "; ".join(f"block 1 face {f} iteration 1, {FIRST}: {n}, {x:.5e}" for f, n, x in want) or "none"
    print(f"    lines: {shown} (expected {wanted}){'' if ok else '  <-- differs'}")
    return ok


def check_one_step(name, path, log):
    dt, field = read_field(path)
    rho0, vel0, p0 = states(name)[0]
    u_init = conserved(rho0, vel0, p0)
    fx = euler_flux(rho0, vel0, p0, (1.0, 0.0, 0.0))
    scale_cell = [RHO, RHO * A, RHO * A, RHO * A, RHO * A * A]
    units = [RHO * A] + [RHO * A * A] * 3 + [RHO * A ** 3]
    ok = True
    for label, c, s_out, sign in (("i = 1", 0, (-1.0, 0.0, 0.0), -1.0), ("i = NX", NX - 1, (1.0, 0.0, 0.0), 1.0)):
        du = [a - b for a, b in zip(conserved(*cell(field, c)), u_init)]
        dev_cell = max(abs(d) / s for d, s in zip(du, scale_cell))
        got = [-d * DX / dt + sign * f for d, f in zip(du, fx)]
        want = euler_flux(rho0, vel0, p0, s_out)
        dev = max(abs(g - w) / un for g, w, un in zip(got, want, units))
        good = dev_cell <= TOL_CELL and dev <= TOL
        ok &= good
        print(f"    {label}: change of the cell {dev_cell:.3e} of its scale (at most {TOL_CELL:.0e}); flux out minus F(q0) . S "
              f"{dev:.3e} (at most {TOL:.0e}){'' if good else '  <-- differs'}")
    u0 = vel0[0] / A
    want_lines = [] if u0 == 0.0 else [(1 if u0 > 0.0 else 2, 1, abs(u0))]
    ok &= check_lines(log, want_lines)
    ok &= check_end(log, len(want_lines), abs(u0))
    return ok


def check_count(name, path, log):
    _, field = read_field(path)
    rho0, vel0, p0 = states(name)[0]
    u_init = conserved(rho0, vel0, p0)
    scale_cell = [RHO, RHO * A, RHO * A, RHO * A, RHO * A * A]
    dev = max(max(abs(a - b) / s for a, b, s in zip(conserved(*cell(field, c)), u_init, scale_cell)) for c in range(NX))
    good = dev <= TOL_CELL
    print(f"    after {NSTEP_COUNT} RK3 steps: largest change of a cell {dev:.3e} of its scale (at most {TOL_CELL:.0e})"
          f"{'' if good else '  <-- differs'}")
    x = COUNT[name]
    return check_lines(log, [(1, 1, x)]) & check_end(log, NSTEP_COUNT * NSTAGE_COUNT, x) & good


def check_tenfold(name, log):
    got = lines_in(log)
    face1 = [g for g in got if g[0] == 1 and g[1] == 1 and g[2] == 0]
    other = [g for g in got if not (g[0] == 1 and g[1] == 1 and g[2] == 0)]
    end = end_line(log)
    ok = len(face1) >= 2 and not other and end is not None
    ok = ok and face1[0][3] == 1 and face1[0][4] and face1[0][5] == 1 and close(face1[0][6], TENFOLD[name][0])
    ok = ok and all(not b[4] and b[3] > a[3] and b[6] > 10.0 * a[6] * (1.0 - TOL_LINE) and b[5] <= face1[0][5]
                    for a, b in zip(face1, face1[1:]))
    ok = ok and end[1] >= face1[-1][6] * (1.0 - TOL_LINE) and end[1] < 10.0 * face1[-1][6] * (1.0 + TOL_LINE)
    print("    lines of face 1: " + "; ".join(show(g) for g in face1) + f"; lines of other faces {len(other)}; end line {end} "
          f"(expected: at least two lines, the first at iteration 1 with {TENFOLD[name][0]:.5e}, each later one more than ten "
          "times the one before, the end at least the last and below ten times it, no other line)"
          + ('' if ok else '  <-- differs'))
    return ok


def check_stop(name, log, rc):
    with open(log, errors="replace") as f:
        lines = f.read().splitlines()
    first = [k for k, l in enumerate(lines) if IN.match(l)]
    stop = [k for k, l in enumerate(lines) if "Integration failed" in l]
    ends = [l for l in lines if END.search(l)]
    ok = rc != 0 and len(first) == 1 and len(stop) >= 1 and first[0] < stop[0] and not ends
    print(f"    exit status {rc} (expected not 0); line of the gas entering at line {first[0] + 1 if first else None} of the log, "
          f"stop at line {stop[0] + 1 if stop else None} (expected after it); end lines {len(ends)} (expected 0)"
          f"{'' if ok else '  <-- differs'}")
    return check_lines(log, [(3, NX, STOP[name][1])]) and ok


def check_no400(name, log):
    """No line that speaks of gas entering, neither a line of a boundary nor the end line: the case has no 400 face."""
    with open(log, errors="replace") as f:
        found = [l for l in f.read().splitlines() if "gas entering" in l]
    ok = not found
    print(f"    lines about gas entering: {len(found)} (expected none: the case has no 400 face){'' if ok else '  <-- differs'}")
    return ok


def check_same(log1, log4):
    name = "threads"
    u, v = THREADS[name]
    want = [(1, 1, u), (3, NT, v)]
    ok = check_lines(log1, want) & check_end(log1, NT + 1, u)
    same = lines_in(log1) == lines_in(log4) and end_line(log1) == end_line(log4)
    print(f"    1 and 4 threads: the same lines and the same end line: {same}{'' if same else '  <-- differs'}")
    return ok and same


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check":
        name, path, log = sys.argv[2:5]
        if name in ONE_STEP:
            ok = check_one_step(name, path, log)
        elif name in COUNT:
            ok = check_count(name, path, log)
        elif name in TENFOLD:
            ok = check_tenfold(name, log)
        elif name in NO400:
            ok = check_no400(name, log)
        else:
            sys.exit(f"unknown case {name}")
        sys.exit(0 if ok else 1)
    elif what == "stop":
        name, log, rc = sys.argv[2], sys.argv[3], int(sys.argv[4])
        if name not in STOP:
            sys.exit(f"unknown case {name}")
        sys.exit(0 if check_stop(name, log, rc) else 1)
    elif what == "same":
        sys.exit(0 if check_same(sys.argv[2], sys.argv[3]) else 1)
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
