"""Part 2 of Inlets (run.sh): inlets at total conditions, inputs and a reader of the outputs.

A 1-D tube along x: NX cells of DX, one cell of DX in y and in z; slip walls (300) on the j faces, null (0) k faces.
Calorically perfect air of test/common/Air (cp = 1004.5 J/(kg K), molar mass 28.970418) at p = 1 bar and T = 300 K.
The inlets take their total conditions from a reservoir at T0 = T and p0 = RATIO p, with the flow direction normal to the
face. One forward Euler step with first-order reconstruction: the change of the density of a boundary cell is then the
mass flux of its boundary face times DT / DX, plus what the interior face takes away.

The face state of an inlet at total conditions (T0, p0, direction) is set by the total enthalpy and the total pressure of the
reservoir and by the invariant a + (gamma-1)/2 w that comes from inside along the outgoing characteristic, w the velocity
along the outward normal. The gas enters when p0 is above p_stop = p (1 + (gamma-1)/2 w_cell / a)^(2 gamma/(gamma-1)), the
pressure the gas of the cell reaches when the outgoing wave brings it to rest; the face then takes the smaller root of the
quadratic, the inflow, and, if that root has w + a3 <= 0, the choked state of normal Mach 1, a3 = a0 / sqrt(1 + (gamma-1)/2).
When the gas does not enter, an inlet (401, 402) is a slip wall, and an ambient face (407) whose cell gas moves in is an
outlet that lets no gas in: the face takes the state of the outgoing wave, a_f = a (p_f/p)^((gamma-1)/(2 gamma)) and
w_f = (Riem - a_f)/g with Riem = a + g w_cell, g = (gamma-1)/2, at p_f = min(p_back, p_stop), so that w_f >= 0 (at p_stop the gas
is at rest at the face).

Cases (CASES below):
  farfield36  ambient (407) on both ends, p0 = p_back = 36 p, gas at rest: choked inflow on both faces
  farfield2   ambient (407) on both ends, p0 = p_back = 2 p, gas at rest: subsonic inflow on both faces
  inlet2      inlet (401) on both ends, p0 = 2 p, gas at rest: subsonic inflow on both faces
  inlet36     inlet (401) on both ends, p0 = 36 p, gas at rest: choked inflow on both faces
  nofeed      inlet (401) at i = 1 with p0 = 0.2 p, extrapolation (400) at i = NX, gas moving in at 0.1 a: p0 is below
              p_stop, the reservoir cannot feed the face, which is a slip wall (no mass through it)
  outlet05    ambient (407) at i = 1 with p0 = p_back = 0.4 p, extrapolation (400) at i = NX, gas moving in at 0.5 a: p0 is below
              p_stop = 0.478 p, the face is an outlet at p_back and lets gas out
  outlet12    the same with the gas moving in at 1.2 a and p0 = p_back = 0.1 p, below p_stop = 0.147 p
  outletstop  as outlet05 with p_back = 0.6 p, above p_stop: the face is at rest at p_stop, no mass through it

Usage (from the test directory):
  python3 -B total_conditions.py case NAME           INPUT/ic.tec and INPUT/bc.txt of case NAME
  python3 -B total_conditions.py check NAME FILE     mass fluxes of the boundary faces in FILE (a field.tec) against the face
                                              state above; exit status 0 when they agree within 1e-9 of rho a
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
G = 0.5 * (GAMMA - 1.0)
RHO = P / (RGAS * T)
A = math.sqrt(GAMMA * RGAS * T)
TOL = 1.0e-9

# name: (type at i = 1, type at i = NX, p0 / p, u / a, p_back / p of a 407)
CASES = {
    "farfield36": (407, 407, 36.0, 0.0, 36.0),
    "farfield2": (407, 407, 2.0, 0.0, 2.0),
    "inlet2": (401, 401, 2.0, 0.0, None),
    "inlet36": (401, 401, 36.0, 0.0, None),
    "nofeed": (401, 400, 0.2, 0.1, None),
    "outlet05": (407, 400, 0.4, 0.5, 0.4),
    "outlet12": (407, 400, 0.1, 1.2, 0.1),
    "outletstop": (407, 400, 0.4, 0.5, 0.6),
}


def outlet_flux(p_back, w_cell):
    """Mass flux into the domain through a 407 face whose cell gas moves in (w_cell < 0) and that the ambient does not feed:
    the state of the outgoing wave at p_f = min(p_back, p_stop); 0 when Riem <= 0."""
    riem = A + G * w_cell
    if riem <= 0.0 or p_back <= 0.0:
        return 0.0
    p_stop = P * (riem / A) ** (2.0 * GAMMA / (GAMMA - 1.0))
    if p_back < p_stop:
        pf, af = p_back, A * (p_back / P) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
    else:
        pf, af = p_stop, riem                                # the gas at rest at the face
    wf = (riem - af) / G
    return -(GAMMA * pf / af ** 2) * wf


def face_flux(btype, p0, p_back, w_cell):
    """Mass flux into the domain through a face of type 401 or 407 (inlet, or ambient with back pressure p_back), cell state
    (RHO, A, P) with velocity w_cell along the outward normal."""
    flux = inflow_flux(p0, w_cell)
    if flux is not None:
        return flux
    if btype == 401:
        return 0.0                                           # slip wall
    if w_cell < 0.0:
        return outlet_flux(p_back, w_cell)
    sys.exit("a 407 face that the gas of the cell leaves is not a case of this test")


def inflow_flux(p0, w_cell):
    """Mass flux into the domain through an inlet face at total conditions (T0 = T, p0), cell state (RHO, A, P) with
    velocity w_cell along the outward normal; None when the gas does not enter."""
    a0 = math.sqrt(GAMMA * RGAS * T)
    riem = A + G * w_cell
    c = (p0 / P) ** ((GAMMA - 1.0) / GAMMA) * (A / a0) ** 2
    d = G                                                    # normal direction: 1 + E1^2 + E2^2 = 1
    ya, yc = G * G + c * d, riem * riem - c * a0 * a0
    if not (riem > 0.0 and yc < 0.0):
        return None
    w = (G * riem - math.sqrt((G * riem) ** 2 - ya * yc)) / ya
    a3 = math.sqrt(a0 * a0 - d * w * w)
    if w + a3 <= 0.0:                                        # choked: normal Mach 1
        a3 = a0 / math.sqrt(1.0 + d)
        w = -a3
    p3 = p0 * (a3 / a0) ** (2.0 * GAMMA / (GAMMA - 1.0))
    return -(GAMMA * p3 / a3 ** 2) * w


def write_case(name):
    t1, t2, ratio, mach, back = CASES[name]
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
        for value in (RHO, mach * A, 0.0, 0.0, P):
            f.write(f" {value:.15E}\n" * NX)
    with open("INPUT/bc.txt", "w") as f:
        for face, i, t in ((1, 1, t1), (2, NX, t2)):
            f.write(f"       1     {i:3d}       1       1       {face}     {t}\n")
            if t == 401:
                f.write(f"    {T:.15E}, {ratio * P:.15E}, normal, normal, 1.0, 1.0,\n")
            elif t == 407:
                f.write(f"    {T:.15E}, {ratio * P:.15E}, {back * P:.15E}, normal, normal, 1.0, 1.0,\n")
        for face, t in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {t:3d}\n")


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
    t1, t2, ratio, mach, back = CASES[name]
    dt, field = read_field(path)
    rho = field.get("rho(1)", field.get("rho1"))
    u0 = mach * A
    unit = RHO * A
    fail = False
    p_back = back * P if back is not None else None
    # face i = 1: outward normal -x, the cell velocity along it is -u0; the interior face takes RHO u0 away
    faces = [("i = 1", (rho[0] - RHO) * DX / dt + RHO * u0, face_flux(t1, ratio * P, p_back, -u0) if t1 in (401, 407) else None)]
    if t2 in (401, 407):
        faces.append(("i = NX", (rho[-1] - RHO) * DX / dt - RHO * u0, face_flux(t2, ratio * P, p_back, u0)))
    for label, got, want in faces:
        ok = abs(got - want) <= TOL * unit
        fail |= not ok
        print(f"    {label}: mass flux in {got / unit:.10f} rho a, face state {want / unit:.10f} rho a"
              f"{'' if ok else '  <-- differs'}")
    sys.exit(1 if fail else 0)


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
