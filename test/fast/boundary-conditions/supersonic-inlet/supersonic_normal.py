"""SupersonicInlet (run.sh): the supersonic inlet (405) with the flow direction normal, inputs and a reader of the outputs.

The 405 imposes the velocity M a d, d = (cos alpha cos beta, sin alpha cos beta, sin beta). With alpha and beta both normal, d
is the unit normal of the face into the domain; with one of them normal, that angle is the one of the normal into the domain,
alpha its azimuth and beta its elevation (alpha is 0 on a face whose normal is along z, which has no azimuth).

A 1-D tube: NX cells of DX along its axis, one cell of DX across; slip walls (300) on the j faces, null (0) k faces, the 405
at one end and an extrapolation face (400) at the other. Calorically perfect air of test/common/Air (cp = 1004.5 J/(kg K),
molar mass 28.970418) at p = 1 bar and T = 300 K, its gas moving in through the inlet face at W = 0.1 a, so that the face
takes the inflow (a 405 face is a slip wall where the gas of the cell moves out or is at rest). The 405 imposes Mach 2 at the
T and p of the gas. One forward Euler step with first-order reconstruction and HLLC: the change of the density of the inlet
cell is (F - RHO W) DT / DX, F the mass flux of the inlet face and RHO W the one that the interior face takes away.

The tube is that along x turned by an orthogonal matrix Q, nodes and velocity (TUBES, Q by columns: the images of the
directions of i, j and k):
  x  Q = I;
  m  the mirror image x -> -x (det Q = -1): i grows along -x;
  r  the turn that takes x to (1, 1, 1)/sqrt(3): a face normal with three components;
  z  the cyclic permutation that takes x to z: i grows along z, the faces 1 and 2 are normal to z.

Cases (CASES below): tube, face of the 405, alpha, beta, and what the case is checked against.
  ref                 x, face 1, 0, 0: the reference, the mass flux F_ref and the velocity w_ref of the inlet cell along the
                      normal into the domain.
  x1n, x2n, m1n, m2n  normal, normal on the faces 1 and 2 of x and m;
  r1n                 normal, normal on face 1 of r;
  x2a                 normal, 0 on face 2 of x;
  z1b, z2b            0.3, normal on the faces 1 and 2 of z:
                      the gas enters along the normal into the domain: F_ref within 1e-9 of rho a, the velocity of the inlet
                      cell along that normal w_ref within 1e-9 of a, and its component across that normal below 1e-9 of a.
  x1c, x2c            0, normal on the faces 1 and 2 of x, against ref and x2e (0, 0 on face 2);
  x1d, x2d            0.2, normal on the faces 1 and 2 of x, against x1f and x2f (0.2, 0);
  z1e                 normal, 0.3 on face 1 of z, against z1g (0, 0.3):
                      the angles the token stands for on that face: the mass flux and the velocity of the inlet cell of the
                      case with the angles as numbers, within 1e-9 of rho a and of a.

Usage (from the test directory):
  python3 -B supersonic_normal.py case NAME   INPUT/ic.tec and INPUT/bc.txt of case NAME
  python3 -B supersonic_normal.py check       every case from OUTPUT/<case>.tec against its reference; exit status 0 when all
                                              pass
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
W = 0.1 * A                          # the speed of the gas of the cell into the domain
MACH = 2.0                           # the Mach number of the 405
TOL = 1.0e-9

R3, R2, R6 = math.sqrt(3.0), math.sqrt(2.0), math.sqrt(6.0)
TUBES = {
    "x": ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "m": ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "r": ((1.0 / R3, 1.0 / R3, 1.0 / R3), (-1.0 / R2, 1.0 / R2, 0.0), (-1.0 / R6, -1.0 / R6, 2.0 / R6)),
    "z": ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
}

# name: (tube, face of the 405, alpha, beta, reference: None for ref, "normal" for the gas entering along the normal into the
# domain, or the case with the angles as numbers)
CASES = {
    "ref": ("x", 1, "0.0", "0.0", None),
    "x1n": ("x", 1, "normal", "normal", "normal"),
    "x2n": ("x", 2, "normal", "normal", "normal"),
    "m1n": ("m", 1, "normal", "normal", "normal"),
    "m2n": ("m", 2, "normal", "normal", "normal"),
    "r1n": ("r", 1, "normal", "normal", "normal"),
    "x2a": ("x", 2, "normal", "0.0", "normal"),
    "z1b": ("z", 1, "0.3", "normal", "normal"),
    "z2b": ("z", 2, "0.3", "normal", "normal"),
    "x1c": ("x", 1, "0.0", "normal", "ref"),
    "x2c": ("x", 2, "0.0", "normal", "x2e"),
    "x2e": ("x", 2, "0.0", "0.0", "-"),
    "x1d": ("x", 1, "0.2", "normal", "x1f"),
    "x2d": ("x", 2, "0.2", "normal", "x2f"),
    "x1f": ("x", 1, "0.2", "0.0", "-"),
    "x2f": ("x", 2, "0.2", "0.0", "-"),
    "z1e": ("z", 1, "normal", "0.3", "z1g"),
    "z1g": ("z", 1, "0.0", "0.3", "-"),
}


def inward(name):
    """Unit normal into the domain of the inlet face: the direction of i on face 1, the opposite one on face 2."""
    tube, face = CASES[name][0], CASES[name][1]
    s = 1.0 if face == 1 else -1.0
    return [s * c for c in TUBES[tube][0]]


def write_case(name):
    tube, face, alpha, beta, _ = CASES[name]
    q = TUBES[tube]
    vel = [W * c for c in inward(name)]
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
                        f.write(f" {q[0][axis] * (DX * i) + q[1][axis] * (DX * j) + q[2][axis] * (DX * k):.16E}\n")
        for value in (RHO, vel[0], vel[1], vel[2], P):
            f.write(f" {value:.16E}\n" * NX)
    with open("INPUT/bc.txt", "w") as f:
        for fc, i in ((1, 1), (2, NX)):
            if fc == face:
                f.write(f"       1     {i:3d}       1       1       {fc}     405\n")
                f.write(f"    {MACH:.15E}, {T:.15E}, {P:.15E}, {alpha}, {beta}, 1.0, 1.0,\n")
            else:
                f.write(f"       1     {i:3d}       1       1       {fc}     400\n")
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
    return time, out, values


def inlet_cell(name):
    """Mass flux into the domain through the inlet face and velocity of the inlet cell, from OUTPUT/<name>.tec; None
    when the file is missing or holds a value that is not finite."""
    path = os.path.join("OUTPUT", name + ".tec")
    if not os.path.isfile(path):
        return None
    dt, field, values = read_field(path)
    if not all(math.isfinite(v) for v in values):
        return None
    c = 0 if CASES[name][1] == 1 else NX - 1
    rho = field.get("rho(1)", field.get("rho1"))
    flux = (rho[c] - RHO) * DX / dt + RHO * W
    return flux, [field[v][c] for v in ("u", "v", "w")]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check():
    res = {name: inlet_cell(name) for name in CASES}
    unit_f, unit_v = RHO * A, A
    if res["ref"] is None:
        print("    [FAIL] ref: no field, or a value that is not finite")
        sys.exit(1)
    f_ref, v_ref = res["ref"]
    w_ref = dot(v_ref, inward("ref"))
    fail = f_ref <= 0.0
    print(f"    [{'FAIL' if fail else 'ok  '}] ref: mass flux in {f_ref / unit_f:.10f} rho a, velocity of the inlet cell "
          f"along the normal into the domain {w_ref / unit_v:.10f} a")
    for name, (tube, face, alpha, beta, against) in CASES.items():
        if against in (None, "-"):
            continue
        head = f"{name} ({tube}, face {face}, {alpha}, {beta})"
        if res[name] is None or (against != "normal" and res[against] is None):
            print(f"    [FAIL] {head}: no field, or a value that is not finite")
            fail = True
            continue
        flux, vel = res[name]
        if against == "normal":
            n = inward(name)
            w = dot(vel, n)
            across = math.sqrt(sum((v - w * c) ** 2 for v, c in zip(vel, n)))
            ok = abs(flux - f_ref) <= TOL * unit_f and abs(w - w_ref) <= TOL * unit_v and across <= TOL * unit_v
            text = (f"mass flux in {flux / unit_f:.10f} rho a, velocity along the normal into the domain "
                    f"{w / unit_v:.10f} a, across {across / unit_v:.2e} a")
        else:
            flux2, vel2 = res[against]
            dv = max(abs(a - b) for a, b in zip(vel, vel2))
            ok = abs(flux - flux2) <= TOL * unit_f and dv <= TOL * unit_v
            text = (f"mass flux in {flux / unit_f:.10f} rho a against {flux2 / unit_f:.10f} rho a of {against}, velocity "
                    f"within {dv / unit_v:.2e} a")
        fail |= not ok
        print(f"    [{'ok  ' if ok else 'FAIL'}] {head}: {text}")
    sys.exit(1 if fail else 0)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check":
        check()
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
