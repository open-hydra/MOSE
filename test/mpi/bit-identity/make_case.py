#!/usr/bin/env python3
"""Inputs of the MPIBitIdentity test (needs-mpi tier): a small 3-D case of 2 x 2 blocks that exercises the main functions of
MOSE, to be solved on 1x1, 1x4, 2x1, 4x1 and 2x2 (MPI ranks x OpenMP threads) with outputs equal byte for byte.

Blocks of NI x NJ x NK cells of DX, the 2 x 2 layout in x and y around a vertex shared by the four blocks (block b = 1 + p +
2 q, p along x, q along y), NK cells in z; block cuts are connections (101), identity orientation, so that with 2 and 4 ranks
they join blocks of different ranks. Outer boundaries:
  x-min  404 (static temperature and mass flux) on block 1, 407 (ambient: total conditions and back pressure) on block 3;
  x-max  406 (outlet, back pressure) on block 2, 400 (extrapolation) on block 4;
  y-min  301 (adiabatic wall, q_wall = 0), y-max 302 (wall at T_WALL);
  z      300 (symmetry planes).
Air of test/common/Air-viscous at P_REF, T_REF, a velocity U0 along x with a vortex about z on the shared vertex and a weaker
swirl about x, so that every direction has gradients from the first stage on. Turbulence: Spalart-Allmaras (rho nu~ = 3 mu)
or SST (rho k, rho omega), the inlets carrying the same specific values.

Usage: make_case.py SA|SST [--outdir DIR]       writes DIR/INPUT/ic.tec and DIR/INPUT/bc.txt
"""
import argparse
import math
import os

NI, NJ, NK = 16, 16, 4
DX = 1.0 / 16.0 * 1.0e-3                # m: blocks of 1 mm x 1 mm x 0.25 mm
R_GAS = 8314.46261815324 / 28.970418            # J/(kg K): air of test/common/Air-viscous (cp = 1004.5)
GAMMA = 1004.5 / (1004.5 - R_GAS)
MU = 1.8467e-5                         # Pa s: air at 300 K
P_REF, T_REF, T_WALL = 1.0e5, 300.0, 350.0
U0, VORTEX, SWIRL = 50.0, 30.0, 10.0   # m/s
RHO = P_REF / (R_GAS * T_REF)
TURB = {"SA": [3.0 * MU / RHO], "SST": [1.0e-3 * U0 ** 2, 1.0e5]}   # specific values: nu~; k, omega
FMT = "{: .17E}"


def state(x, y, z, model):
    """rho, u, v, w, p and the conservative turbulence variables at a cell centre."""
    xc, yc, zc = NI * DX, NJ * DX, 0.5 * NK * DX
    radius = 0.25 * NI * DX
    dx, dy, dz = x - xc, y - yc, z - zc
    g = VORTEX / radius * math.exp(0.5 * (1.0 - (dx * dx + dy * dy) / radius ** 2))
    s = SWIRL / radius * math.exp(-(dy * dy + dz * dz) / radius ** 2)
    u = U0 - g * dy
    v = g * dx - s * dz
    w = s * dy
    return [RHO, u, v, w, P_REF] + [RHO * t for t in TURB[model]]


def write_ic(path, model):
    names = ' "rhoNuTilde"' if model == "SA" else ' "rhoK" "rhoOmega"'
    nvar = 8 + len(TURB[model])
    with open(path, "w") as f:
        f.write(' VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"' + names + "\n")
        for q in range(2):
            for p in range(2):
                b = 1 + p + 2 * q
                xs = [(p * NI + i) * DX for i in range(NI + 1)]
                ys = [(q * NJ + j) * DX for j in range(NJ + 1)]
                zs = [k * DX for k in range(NK + 1)]
                f.write(f' ZONE T="Block{b}", I={NI + 1}, J={NJ + 1}, K={NK + 1}, DATAPACKING=BLOCK, '
                        f'VARLOCATION=([1-3]=NODAL,[4-{nvar}]=CELLCENTERED)\n')
                for axis in range(3):
                    for k in range(NK + 1):
                        for j in range(NJ + 1):
                            for i in range(NI + 1):
                                f.write(" " + FMT.format((xs[i], ys[j], zs[k])[axis]) + "\n")
                cells = [state(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]), 0.5 * (zs[k] + zs[k + 1]), model)
                         for k in range(NK) for j in range(NJ) for i in range(NI)]
                for var in range(nvar - 3):
                    f.write("".join(" " + FMT.format(c[var]) + "\n" for c in cells))


def outer(model, b, face):
    """(type, data line) of an outer face of block b."""
    turb = "".join(f" {t:.17E}," for t in TURB[model])
    m2 = (U0 / math.sqrt(GAMMA * R_GAS * T_REF)) ** 2
    t0 = T_REF * (1.0 + 0.5 * (GAMMA - 1.0) * m2)
    p0 = P_REF * (t0 / T_REF) ** (GAMMA / (GAMMA - 1.0))
    if face == 1:
        if b == 1:
            return 404, f"    {T_REF:.17E}, {RHO * U0:.17E}, normal, normal, 1.0, 1.0,{turb}"
        return 407, f"    {t0:.17E}, {p0:.17E}, {P_REF:.17E}, normal, normal, 1.0, 1.0,{turb}"
    if face == 2:
        return (406, f"    {P_REF:.17E}, 1.0,") if b == 2 else (400, None)
    if face == 3:
        return 301, "    0.0, 0.0, 0.0,"
    if face == 4:
        return 302, f"    {T_WALL:.17E}, 0.0, 0.0,"
    return 300, None


def write_bc(path, model):
    rows = []
    for q in range(2):
        for p in range(2):
            b = 1 + p + 2 * q
            for face in range(1, 7):
                d = (face + 1) // 2
                if d == 1:
                    span = [(1 if face == 1 else NI, j, k) for k in range(1, NK + 1) for j in range(1, NJ + 1)]
                    nb = 1 + (1 - p) + 2 * q if (face == 2 and p == 0) or (face == 1 and p == 1) else None
                elif d == 2:
                    span = [(i, 1 if face == 3 else NJ, k) for k in range(1, NK + 1) for i in range(1, NI + 1)]
                    nb = 1 + p + 2 * (1 - q) if (face == 4 and q == 0) or (face == 3 and q == 1) else None
                else:
                    span = [(i, j, 1 if face == 5 else NK) for j in range(1, NJ + 1) for i in range(1, NI + 1)]
                    nb = None
                for (i, j, k) in span:
                    rows.append(f"{b:8d}{i:8d}{j:8d}{k:8d}{face:8d}")
                    if nb is not None:
                        fs = face + 1 if face % 2 == 1 else face - 1
                        si = {1: NI, 2: 1}.get(face, i)
                        sj = {3: NJ, 4: 1}.get(face, j)
                        rows[-1] += f"{101:8d}\n{nb:8d}{si:8d}{sj:8d}{k:8d}{fs:8d}{1:8d}{0:8d}{0:8d}{1:8d}"
                    else:
                        t, data = outer(model, b, face)
                        rows[-1] += f"{t:8d}" + ("" if data is None else "\n" + data)
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("model", choices=["SA", "SST"])
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.outdir, "INPUT"), exist_ok=True)
    write_ic(os.path.join(a.outdir, "INPUT", "ic.tec"), a.model)
    n = write_bc(os.path.join(a.outdir, "INPUT", "bc.txt"), a.model)
    print(f"2x2 blocks of {NI}x{NJ}x{NK} cells, turbulence {a.model}, {n} boundary records")


if __name__ == "__main__":
    main()
