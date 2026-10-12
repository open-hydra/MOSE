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

The same on cells that are not uniform, where the cell gradients of the source terms must take the differences of the cell
centres, and the ghost cells of the wall the position at which their value is right:
  stretched  4 x 12 cells growing by r = 1.2 from the wall (first cell 6 um; the rows next to the symmetry plane, the last
             two, are not checked): with the metric of the cell the central difference across the wall cell gives
             (1 + (r - 1)/4) a, and the cb2 term about 0.1 of the scale. The fluxes take the diffusion coefficient
             (mu + rho nu~)/sigma of a face from the mean of its two cells, which is the value at the midpoint of their
             centres, not at the face: on these cells the flux divergence differs from rho a^2/sigma by (rho a^2/sigma)
             ((c_top - y_top) - (c_bot - y_bot))/h, c the midpoints, y the faces (the wall face is exact), that is
             (r - 1)/4 in the wall cell and (r - 1)^2/(4 r) in the others, over cb2 in units of the scale. The check is
             on the right-hand side less this term;
  skewed     6 x 16 cells skewed by 20 degrees (x += y tan 20), the sides x = const of 400 type (their columns are not
             checked), with T and so rho linear along the wall (T = T0 + 2000 K/m (x - x_mid), wall at T(x)). The ghost
             cell of the wall takes T = 2 T_w - T_c at the reflection 2 x_f - x_c of the cell centre and the pressure of the
             cell, so its density is that of the field at that point (Ghost_Wall), and the cell metric of the density takes
             that point; a density at a point where it is not right gives grad rho a normal component, which the term
             -(nu + nu~) grad rho . grad nu~ / sigma turns into about 0.3 of the scale (a copy at the reflection, as before
             the metrics of the sources took the cell itself for it). The wall distance of MOSE is that to the nearest wall face centre, not y, on skewed cells:
             the equilibrium does not hold, so the same column at uniform T is solved too and the difference of
             (right-hand side)/rho between the two runs, over (cb2/sigma) a^2, is checked; it holds the terms in grad rho;
  skewed3d   the same in 3-D: 5 x 16 x 5 cells, x += y tan 20, z += y tan 10, T linear in x and z, the sides of 400 type.

SA-RC, the rotation-curvature correction, takes the gradient of the strain rate from the velocity gradients of the cell and
of its neighbours; next to a face between blocks the neighbour across is the ghost cell, which must hold the gradient of the
cell across (and the metric its position), so that the result does not depend on how the mesh is divided into blocks. The
uniform column with a shear flow whose strain rate varies along x, u = U0 (y/H) (1 + AMP sin(2 pi x/Lx)), v = w = 0
(DS/Dt = u dS/dx is not zero and varies across the cut), one explicit step with turbulence = SA-RC, the change of rho nu~
over the step (the only quantity the correction acts on in one step) against that of one block, periodic on itself (201):
  two       two blocks of 2 x 16 cells along x (a face between the blocks, 101, and the periodic pair, 201);
  permuted  the same two blocks with the indices of the second permuted (its i along y, j along x: the face between the
            blocks is its face 3, the periodic pair its face 4, the wall its face 1).
Both to 1e-9 of its largest value (one ulp of rho nu~ is about 1e-12 of it). A ghost cell that copies the gradient of the
cell, with the cell itself as its position (a one-sided difference), changes it by about 3e-2 in the two columns next to
the face between the blocks; the gradient of the cell across with the cell itself as its position, by about 3e-2 too.
The other variables differ next to that face in the last bits, with or without SA-RC: the fluxes of a face between blocks
are not computed as those of an interior face (ConnectionFlux).

Usage (from the test directory):
  python3 -B wall_gradient.py case [CASE [uni]]   INPUT/ic.tec and INPUT/bc.txt; CASE uniform (default), stretched, skewed,
                                              skewed3d; uni: at uniform temperature (skewed cases); shear LAYOUT: the shear
                                              column of SA-RC, LAYOUT one, two or permuted
  python3 -B wall_gradient.py check FIELD DIAG [CASE]                    exit status 0 when the right-hand side passes
  python3 -B wall_gradient.py checkdiff CASE FIELD DIAG FIELD_UNI DIAG_UNI   the same for the skewed cases
  python3 -B wall_gradient.py sarcsame DIR_ONE DIR LAYOUT   exit status 0 when the SA-RC run of DIR (field.vtm and vtk/, vtk
                                                       raw) gives the change of rho nu~ of one block in DIR_ONE
"""
import base64
import math
import os
import re
import struct
import sys

RGAS = 8314.46261815324 / 28.970418          # gas constant of FLINT for air, J/(kg K)
MU = 1.1858686e-5                   # test/common/Air-viscous
T0, P0 = 300.0, 1.0e5
SIGMA, CB2 = 2.0 / 3.0, 0.622       # Spalart-Allmaras
NI, NJ, DX, DY = 4, 16, 1.0e-4, 1.0e-5
A = 6.6e-5
ROWS = 14
R_STRETCH, DY1_STRETCH, NJ_STRETCH = 1.2, 6.0e-6, 12   # first cell 6 um: rho nu~ = 2.3e-10 there, above the floor of
                                                         # 1e-10 that MOSE keeps; 12 rows: chi below 1.6e-3 at the top
TX, TZ = 2000.0, 1000.0             # K/m: temperature gradient along the wall of the skewed cases
U0, AMP = 1.0, 0.3                  # SA-RC: the shear flow u = U0 (y/H) (1 + AMP sin(2 pi x/Lx)) of the uniform column


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


def shear_blocks(layout):
    """Blocks of the shear column of SA-RC: (ni, nj, g), g mapping node (i, j) of the block, and cell (i, j) counted from 0,
    to the node and the cell of the column."""
    half = NI // 2
    if layout == 'one':
        return [(NI, NJ, lambda i, j: (i, j))]
    second = (half, NJ, lambda i, j: (half + i, j)) if layout == 'two' else (NJ, half, lambda i, j: (half + j, i))
    return [(half, NJ, lambda i, j: (i, j)), second]


def write_shear(layout):
    """The uniform column with the shear flow of the SA-RC checks, as one block (one), as two blocks of 2 x 16 cells along x
    (two: block 1 at x in [0, 2 dx], block 2 at [2 dx, 4 dx]; the face between them, 101, and the periodic pair, 201, from
    block 1, face 1, to block 2, face 2), or as the same two blocks with block 2 written with i along y and j along x
    (permuted). Every node and cell value is that of the column, bit for bit."""
    rho = P0 / (RGAS * T0)
    half = NI // 2
    lx, hc = NI * DX, NJ * DY
    blocks = shear_blocks(layout)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        for b, (ni, nj, g) in enumerate(blocks, 1):
            f.write(f" ZONE  T = B{b}-IG, I={ni + 1}, J={nj + 1}, K=2, DATAPACKING=BLOCK, "
                    "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for c in range(3):
                for k in range(2):
                    for j in range(nj + 1):
                        for i in range(ni + 1):
                            ig, jg = g(i, j)
                            f.write(f" {(ig * DX, jg * DY, k * DX)[c]:.16E}\n")
            cells = [g(i, j) for j in range(nj) for i in range(ni)]
            n = ni * nj
            zero = [0.0] * n
            u = [U0 * ((jg + 0.5) * DY / hc) * (1.0 + AMP * math.sin(2.0 * math.pi * (ig + 0.5) * DX / lx))
                 for ig, jg in cells]
            for col in ([rho] * n, u, zero, zero, [P0] * n, [rho * A * (jg + 0.5) * DY for _, jg in cells]):
                f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(b, i, j, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (b, i, j, 1, face, t)

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    wall = f"    {T0:.6E}, 0.000000E+00, 0.000000E+00,\n"
    out = []
    if layout == 'one':
        for j in range(1, NJ + 1):
            out += [rec(1, 1, j, 1, 201), con(1, NI, j, 1, 2, 1, 0, 0, 1)]
            out += [rec(1, NI, j, 2, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    elif layout == 'two':
        for j in range(1, NJ + 1):
            out += [rec(1, 1, j, 1, 201), con(2, half, j, 1, 2, 1, 0, 0, 1)]
            out += [rec(1, half, j, 2, 101), con(2, 1, j, 1, 1, 1, 0, 0, 1)]
            out += [rec(2, 1, j, 1, 101), con(1, half, j, 1, 2, 1, 0, 0, 1)]
            out += [rec(2, half, j, 2, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    else:
        # block 2 cell (i, j) is the column cell (half + j, i): its face 3 is the face between the blocks, its face 4 the
        # periodic side; +j of block 1 is +i of block 2 (d11 = 1), +k is +k
        for j in range(1, NJ + 1):
            out += [rec(1, 1, j, 1, 201), con(2, j, half, 1, 4, 1, 0, 0, 1)]
            out += [rec(1, half, j, 2, 101), con(2, j, 1, 1, 3, 1, 0, 0, 1)]
            out += [rec(2, j, 1, 3, 101), con(1, half, j, 1, 2, 1, 0, 0, 1)]
            out += [rec(2, j, half, 4, 201), con(1, 1, j, 1, 1, 1, 0, 0, 1)]
    for b, (ni, nj, g) in enumerate(shear_blocks(layout), 1):
        if b == 2 and layout == 'permuted':                     # wall on face 1 (y = 0), symmetry on face 2
            for j in range(1, nj + 1):
                out += [rec(b, 1, j, 1, 302), wall]
                out.append(rec(b, ni, j, 2, 300))
        else:
            for i in range(1, ni + 1):
                out += [rec(b, i, 1, 3, 302), wall]
                out.append(rec(b, i, nj, 4, 300))
        for face in (5, 6):
            out += [rec(b, i, j, face, 0) for j in range(1, nj + 1) for i in range(1, ni + 1)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def vtk_cells(path):
    """Float64 cell arrays of one block of a vtk output of MOSE (DataArray "binary": base64 of an Int32 byte count and the
    little-endian data; "appended": raw AppendedData, each array an Int32 byte count and the data, at its offset)."""
    raw = open(path, "rb").read()
    head, _, appended = raw.partition(b'<AppendedData encoding="raw">')
    text = head.decode("latin-1")
    if appended:
        appended = appended[appended.index(b"_") + 1:]
    cell = re.search(r"<CellData[^>]*>(.*?)</CellData>", text, re.S)
    if not cell:
        sys.exit(f"{path}: no CellData")
    out = {}
    for m in re.finditer(r'<DataArray([^>]*?)(/>|>(.*?)</DataArray>)', cell.group(1), re.S):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        name = re.search(r'Name="+([^"]*)"+', m.group(1)).group(1)
        if attrs.get("type") != "Float64":
            continue
        if attrs.get("format") == "binary":
            blob = base64.b64decode("".join(m.group(3).split()))
            data = blob[4:4 + struct.unpack("<i", blob[:4])[0]]
        elif attrs.get("format") == "appended":
            off = int(attrs["offset"])
            data = appended[off + 4:off + 4 + struct.unpack("<i", appended[off:off + 4])[0]]
        else:
            sys.exit(f"{path}: format {attrs.get('format')} not handled")
        out[name] = struct.unpack("<%dd" % (len(data) // 8), data)
    return out


def shear_column(directory, layout):
    """Cell arrays of the SA-RC run of DIR (field.vtm and vtk/), each a dict (ig, jg) -> value over the column."""
    files = re.findall(r'file="([^"]+)"', open(os.path.join(directory, "field.vtm")).read())
    blocks = shear_blocks(layout)
    if len(files) != len(blocks):
        sys.exit(f"{directory}: {len(files)} blocks in field.vtm, {len(blocks)} expected")
    out = {}
    for name, (ni, nj, g) in zip(files, blocks):
        path = os.path.normpath(os.path.join(directory, name.replace("../", "").replace("//", "/")))
        for var, values in vtk_cells(path).items():
            cells = [g(i, j) for j in range(nj) for i in range(ni)]
            out.setdefault(var, {}).update(zip(cells, values))
    return out


def sarc_same(dir_one, dir_other, layout):
    """The change of rho nu~ over the step, which SA-RC acts on, that of one block to 1e-9 of its largest value (one ulp
    of rho nu~ is about 1e-12 of it; a ghost cell that copies the gradient of the cell changes it by about 3e-2 in the two
    columns next to the face between the blocks). The largest difference of every variable is printed: the inviscid and
    viscous fluxes of a face between blocks differ from those of an interior face in the last bits (ConnectionFlux)."""
    a, b = shear_column(dir_one, 'one'), shear_column(dir_other, layout)
    if set(a) != set(b) or not a:
        print(f"    variables differ: {sorted(a)} / {sorted(b)}")
        sys.exit(1)
    vel = max(abs(v) for var in ("u", "v", "w") for v in a[var].values())
    worst = {}
    for var in sorted(a):
        scale = vel if var in ("u", "v", "w") else max(abs(v) for v in a[var].values())
        worst[var] = max(abs(b[var][c] - x) for c, x in a[var].items()) / scale
    rho = P0 / (RGAS * T0)
    change = {lab: {c: v - rho * A * (c[1] + 0.5) * DY for c, v in col["rhoNuTilde"].items()} for lab, col in (("a", a), ("b", b))}
    scale = max(abs(v) for v in change["a"].values())
    err = max(abs(change["b"][c] - v) for c, v in change["a"].items()) / scale
    print(f"    {layout}: change of rho nu~ over the step against one block, largest difference {err:.3e} of its largest "
          f"value; largest differences of the variables (velocity over the largest velocity): "
          + ", ".join(f"{k} {v:.1e}" for k, v in worst.items()))
    sys.exit(0 if err < 1e-9 else 1)


def layout(case):
    """Nodes X[k][j][i] (x, y, z), cells per direction, and the type of the i and k sides."""
    if case == 'stretched':
        ni, nj, nk = NI, NJ_STRETCH, 1
        ys = [0.0]
        for j in range(nj):
            ys.append(ys[-1] + DY1_STRETCH * R_STRETCH ** j)
        nodes = [[[(i * DX, ys[j], k * DX) for i in range(ni + 1)] for j in range(nj + 1)] for k in range(nk + 1)]
        return nodes, (ni, nj, nk), 201, 0
    if case == 'skewed':
        ni, nj, nk = 6, NJ, 1
        tx, tz = math.tan(math.radians(20.0)), 0.0
    else:
        ni, nj, nk = 5, NJ, 5
        tx, tz = math.tan(math.radians(20.0)), math.tan(math.radians(10.0))
    nodes = [[[(i * DX + j * DY * tx, j * DY, k * DX + j * DY * tz) for i in range(ni + 1)] for j in range(nj + 1)]
             for k in range(nk + 1)]
    return nodes, (ni, nj, nk), 400, (400 if nk > 1 else 0)


def centre(nodes, i, j, k):
    """Mean of the eight nodes of cell (i, j, k), 0-based."""
    return tuple(sum(nodes[k + c][j + b][i + a][d] for a in (0, 1) for b in (0, 1) for c in (0, 1)) / 8.0 for d in range(3))


def temperature(case, x, uni):
    if case in ('uniform', 'stretched') or uni:
        return T0
    nodes, (ni, nj, nk), _, _ = layout(case)
    xm, zm = 0.5 * ni * DX, 0.5 * nk * DX
    return T0 + TX * (x[0] - xm) + (TZ * (x[2] - zm) if nk > 1 else 0.0)


def write_case_general(case, uni):
    nodes, (ni, nj, nk), side_i, side_k = layout(case)
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    cells = [centre(nodes, i, j, k) for k in range(nk) for j in range(nj) for i in range(ni)]
    with open("INPUT/ic.tec", "w") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p" "rhoNuTilde"\n')
        f.write(f" ZONE  T = B1-IG, I={ni + 1}, J={nj + 1}, K={nk + 1}, DATAPACKING=BLOCK, "
                "VARLOCATION=([1-3]=NODAL,[4-9]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for c in range(3):
            for k in range(nk + 1):
                for j in range(nj + 1):
                    for i in range(ni + 1):
                        f.write(f" {nodes[k][j][i][c]:.16E}\n")
        rho = [P0 / (RGAS * temperature(case, x, uni)) for x in cells]
        n = ni * nj * nk
        zero = [0.0] * n
        for col in (rho, zero, zero, zero, [P0] * n, [r * A * x[1] for r, x in zip(rho, cells)]):
            f.write("".join(f" {v:.16E}\n" for v in col))

    def rec(i, j, k, face, t):
        return '%9d %7d %7d %7d %7d %7d\n' % (1, i, j, k, face, t)

    def con(*a):
        return '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % a

    out = []
    for k in range(1, nk + 1):
        for j in range(1, nj + 1):
            if side_i == 201:
                out += [rec(1, j, k, 1, 201), con(1, ni, j, k, 2, 1, 0, 0, 1)]
                out += [rec(ni, j, k, 2, 201), con(1, 1, j, k, 1, 1, 0, 0, 1)]
            else:
                out += [rec(1, j, k, 1, side_i), rec(ni, j, k, 2, side_i)]
    for k in range(1, nk + 1):
        for i in range(1, ni + 1):
            xf = [0.25 * sum(nodes[k - 1 + c][0][i - 1 + a][d] for a in (0, 1) for c in (0, 1)) for d in range(3)]
            out += [rec(i, 1, k, 3, 302), f"    {temperature(case, xf, uni):.16E}, 0.000000E+00, 0.000000E+00,\n"]
            out.append(rec(i, nj, k, 4, 300))
    for j in range(1, nj + 1):
        for i in range(1, ni + 1):
            out += [rec(i, j, 1, 5, side_k), rec(i, j, nk, 6, side_k)]
    with open("INPUT/bc.txt", "w") as f:
        f.write("".join(out))


def read_zone(path):
    """Cell values of a one-zone MOSE ASCII file (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    ni, nj, nk = (int(re.search(k + r'=(\d+)', lines[1]).group(1)) for k in ('I', 'J', 'K'))
    values = [float(v) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
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


def check_general(case, field, diag):
    """stretched: the right-hand side of rho nu~ over (cb2/sigma) rho a^2, as for the uniform column."""
    nodes, (ni, nj, nk), _, _ = layout(case)
    f, d = read_zone(field), read_zone(diag)
    cells = [centre(nodes, i, j, k) for k in range(nk) for j in range(nj) for i in range(ni)]
    rho = P0 / (RGAS * T0)
    scale = CB2 / SIGMA * rho * A * A
    rhs = [(new - rho * A * x[1]) / t / scale for new, x, t in zip(f["rhoNuTilde"], cells, d["dt"])]
    # less the flux divergence error of the mean of the two cells as the face value of the diffusion coefficient
    yf = [nodes[0][j][0][1] for j in range(nj + 1)]
    yc = [cells[j * ni][1] for j in range(nj)]
    mid = [0.0] + [0.5 * (yc[j - 1] + yc[j]) for j in range(1, nj)]          # mid[j]: face j, 0 = wall (exact there)
    for j in range(nj - 1):
        top = mid[j + 1] - yf[j + 1]
        bot = (mid[j] - yf[j]) if j > 0 else 0.0
        pred = (top - bot) / (yf[j + 1] - yf[j]) / CB2
        for i in range(ni):
            rhs[j * ni + i] -= pred
    wall = max(abs(rhs[j * ni + i]) for j in range(1) for i in range(ni))
    rows = min(ROWS, nj - 2)
    rest = max(abs(rhs[j * ni + i]) for j in range(1, rows) for i in range(ni))
    print(f"    {case}: right-hand side of rho nu~ over (cb2/sigma) rho a^2, less the term of the face coefficient: wall "
          f"cells {wall:.3e}, rows 2 to {rows} {rest:.3e}")
    sys.exit(0 if max(wall, rest) < 0.05 else 1)


def check_diff(case, field, diag, field_uni, diag_uni):
    """skewed cases: (right-hand side)/rho with T linear along the wall minus the same at uniform T, over (cb2/sigma) a^2,
    in the cells not next to a 400 side."""
    nodes, (ni, nj, nk), _, side_k = layout(case)
    cells = [centre(nodes, i, j, k) for k in range(nk) for j in range(nj) for i in range(ni)]
    out = []
    for fld, dg, uni in ((field, diag, False), (field_uni, diag_uni, True)):
        f, d = read_zone(fld), read_zone(dg)
        r = []
        for new, x, t in zip(f["rhoNuTilde"], cells, d["dt"]):
            rho = P0 / (RGAS * temperature(case, x, uni))
            r.append((new - rho * A * x[1]) / t / rho)
        out.append(r)
    scale = CB2 / SIGMA * A * A
    wall, rest = 0.0, 0.0
    for k in range(nk):
        if nk > 1 and k in (0, nk - 1):
            continue
        for j in range(ROWS):
            for i in range(1, ni - 1):
                c = (k * nj + j) * ni + i
                e = abs(out[0][c] - out[1][c]) / scale
                if j == 0:
                    wall = max(wall, e)
                else:
                    rest = max(rest, e)
    print(f"    {case}: (right-hand side)/rho, T along the wall minus uniform T, over (cb2/sigma) a^2: wall cells "
          f"{wall:.3e}, rows 2 to {ROWS} {rest:.3e}")
    sys.exit(0 if max(wall, rest) < 0.05 else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        case = sys.argv[2] if len(sys.argv) > 2 else 'uniform'
        if case == 'uniform':
            write_case()
        elif case == 'shear':
            write_shear(sys.argv[3])
        else:
            write_case_general(case, len(sys.argv) > 3 and sys.argv[3] == 'uni')
    elif what == "check":
        if len(sys.argv) > 4 and sys.argv[4] != 'uniform':
            check_general(sys.argv[4], sys.argv[2], sys.argv[3])
        check(sys.argv[2], sys.argv[3])
    elif what == "checkdiff":
        check_diff(*sys.argv[2:7])
    elif what == "sarcsame":
        sarc_same(*sys.argv[2:5])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
