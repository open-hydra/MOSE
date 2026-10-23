"""Inputs of the ExitCodes fast test.

Every case starts from the 100-cell Sod tube of test/1D/Sod79 (INPUT/ic_x1.tec, INPUT/bc_x1.txt); run.sh breaks it in one way,
mostly through input.ini, so the inputs here are the Sod ones, or the Sod ones changed in one value.

Usage (from test/fast/input-checks/exit-codes):
  python3 -B make_case.py sod SOD_INPUT_DIR
      writes INPUT/ic.tec and INPUT/bc.txt, copies of the Sod ones.
  python3 -B make_case.py nan SOD_INPUT_DIR
      the same, with the velocity u of cell 50 written as NaN in INPUT/ic.tec.
  python3 -B make_case.py inf SOD_INPUT_DIR
      the same, with the coordinate y of node (i, j, k) = (20, 1, 1) written as Infinity in INPUT/ic.tec (the nodes
      in the order of the file, i first, from 0).
  python3 -B make_case.py steam SOD_INPUT_DIR PHASE
      the same, the tube holding instead the gas of the phase file PHASE (the 13 C-H-O species of
      test/2D/euler/rocket-nozzle/equilibrium): steam with 0.2 % of hydrogen, 1.6 % of oxygen and 0.2 % of
      carbon dioxide by mass, at rest, at 1e5 Pa and 1500 K in the left half and 1e4 Pa and 1200 K in the right
      half. It is far from chemical equilibrium (the equilibrium burns the hydrogen), so chemistry = equilibrium
      changes it.
  python3 -B make_case.py mixture SOD_IC OUT
      writes OUT: the 100-cell Sod tube of SOD_IC (one species, rho1) as a
      mixture of the four species of test/1D/diff1 (order of its
      INPUT/phase.txt: CH4, H2O, N2, O2) with the mass fractions of air
      (0, 0, 0.767, 0.233); density, velocity and pressure are those of the tube.
  python3 -B make_case.py soot PREMIXED_INPUT SOD_IC DEST sa|soot|plain|sa-plain
      writes DEST/phase.txt, thermo.dat, transport.dat and ic.tec: six species,
      with the names the soot model looks for (C2H2, O2, N2, H2, CO and C(gr),
      Mod_Soot.f90), whose tables are taken from test/1D/premixed-flame (C2H4
      stands for C2H2 and C for C(gr): their properties do not matter here), in
      a uniform state (1.16 kg/m3, 10 m/s, 1 bar, mass fractions 0.05, 0.233,
      0.716999, 0, 0, 1e-6) on the mesh of SOD_IC; with soot or sa the number
      density of the soot (1e14 /m3, the variable rho*Np of the soot model with
      C(gr)), with sa and sa-plain rho nu~ = 3 mu (sa-plain: no soot).
  python3 -B make_case.py diagnostic-check DIAG REF NAME...
      exit status 0 when OUTPUT/diagnostic.tec DIAG names as many variables as
      it declares and writes, the names are NAME..., and its column dt is that
      of the diagnostic REF (the same run without soot: plain, or sa-plain
      with Spalart-Allmaras) within 1e-6.
  python3 -B make_case.py duct410 FORM
      writes INPUT/ic.tec, INPUT/bc.txt, q2d-ends.dat, q2d-jfaces.dat and
      q2d-side.dat: the duct of BCInput (100 cells on 0 <= x <= 1 m, one cell
      of 0.01 m in y, two cells of 5e-4 m in z) in two blocks of 50 cells
      joined at x = 0.5 m (101), air at rest at 1e5 Pa and 300 K, with faces
      of type 410 on the six kinds of faces, read from three files in the form
      that the boundary-condition builder of ATLAS writes (one zone per block
      and time, the nodes of the face of the block, values with 17 digits, two
      times, 0 and 1e-3 s): q2d-ends.dat for face 1 of block 1 and face 2 of
      block 2, q2d-jfaces.dat for face 4 of block 1 and face 3 of block 2,
      q2d-side.dat for face 5 of block 1 and face 6 of block 2. The files hold
      the air at rest at 300 K, at 1e5 Pa for block 1 and 1.02e5 Pa for block
      2. The blocks are joined by type 101; the faces in y and z without a
      410 have the null type 0. FORM changes q2d-side.dat: valid (no change),
      swapped (the zones of the two blocks under each other's title),
      no-block-2 (no zone of block 2), no-time-2 (no zone of block 2 at the
      second time), short-zone (the zones of block 2 with I=26, its first 25
      cells), time-2-other-block (the zone B1_T2 holds the nodes and the
      values of block 2), moved-5e-5 and moved-2e-4 (every node moved by that
      along z, 0.5e-2 and 2e-2 of the shortest edge of the faces 5 and 6,
      0.01 m); digits-15 writes the coordinates of the three files with 15
      significant digits.
"""
import os
import re
import shutil
import sys

Y_AIR = (0.0, 0.0, 0.767, 0.233)
SOOT_SPECIES = (("C2H2", 26.038, "C2H4"), ("O2", 31.998, "O2"), ("N2", 28.014, "N2"), ("H2", 2.016, "H2"),
                ("CO", 28.010, "CO"), ("C(gr)", 12.011, "C"))
Y_SOOT = (0.05, 0.233, 0.716999, 0.0, 0.0, 1.0e-6)
RHO_S, U_S, P_S, NP_S, MU_S = 1.16, 10.0, 1.0e5, 1.0e14, 1.8e-5
NAN_CELL = 50
INF_NODE = (20, 1, 1)                                          # (i, j, k) from 0, as the solver counts the nodes
RU = 8314.46261815324                                          # J/(kmol K), universal gas constant of FLINT
STEAM = {"H2O": 0.98, "H2": 0.002, "O2": 0.016, "CO2": 0.002}   # mass fractions
LEFT, RIGHT = (1e5, 1500.0), (1e4, 1200.0)                     # (p [Pa], T [K]) of the two halves


def read_tube(src):
    """Header, zone line and the values of a Sod tube in block format: nodes x, y, z, then cells rho1, u, v, w, p."""
    with open(src) as f:
        lines = f.read().split('\n')
    header, zone = lines[0], lines[1]
    names = re.findall(r'"([^"]+)"', header)
    assert names == ['x', 'y', 'z', 'rho1', 'u', 'v', 'w', 'p'], names
    ni, nj, nk = (int(re.search(r'%s=\s*(\d+)' % c, zone).group(1)) for c in 'IJK')
    nodes, cells = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    vals = ' '.join(lines[2:]).split()
    assert len(vals) == 3 * nodes + 5 * cells, (len(vals), nodes, cells)
    return header, zone, nodes, cells, vals


def sod(sod_dir):
    for name, src in (("ic.tec", "ic_x1.tec"), ("bc.txt", "bc_x1.txt")):
        dst = os.path.join("INPUT", name)
        if os.path.lexists(dst):
            os.remove(dst)
        shutil.copyfile(os.path.join(sod_dir, src), dst)


def nan(sod_dir):
    sod(sod_dir)
    header, zone, nodes, cells, vals = read_tube(os.path.join(sod_dir, "ic_x1.tec"))
    vals[3 * nodes + cells + NAN_CELL - 1] = 'NaN'          # u of cell NAN_CELL
    os.remove(os.path.join("INPUT", "ic.tec"))
    with open(os.path.join("INPUT", "ic.tec"), 'w') as f:
        f.write(header + '\n' + zone + '\n')
        for v in vals:
            f.write(' %s\n' % v)


def inf(sod_dir):
    sod(sod_dir)
    header, zone, nodes, cells, vals = read_tube(os.path.join(sod_dir, "ic_x1.tec"))
    ni, nj = (int(re.search(r'%s=\s*(\d+)' % c, zone).group(1)) for c in 'IJ')
    i, j, k = INF_NODE
    vals[nodes + i + ni * (j + nj * k)] = 'Infinity'        # y of node INF_NODE: the y values follow the nodes' x values
    os.remove(os.path.join("INPUT", "ic.tec"))
    with open(os.path.join("INPUT", "ic.tec"), 'w') as f:
        f.write(header + '\n' + zone + '\n')
        for v in vals:
            f.write(' %s\n' % v)


def steam(sod_dir, phase):
    """INPUT/ic.tec: the mesh of the Sod tube (its nodal coordinates, as written), the partial densities of the
    species of PHASE, the velocity (zero) and the pressure of the steam."""
    sod(sod_dir)
    rows = [line.split() for line in open(phase).read().splitlines()[1:] if line.strip()]
    species = [row[0] for row in rows]
    weight = {row[0]: float(row[1]) for row in rows}
    unknown = sorted(set(STEAM) - set(species))
    if unknown:
        sys.exit("make_case.py: species not in %s: %s" % (phase, " ".join(unknown)))
    head, zone, body = open(os.path.join(sod_dir, "ic_x1.tec")).read().split("\n", 2)
    ni, nj, nk = (int(re.search(r"\b%s=(\d+)" % key, zone).group(1)) for key in ("I", "J", "K"))
    nodes = ni * nj * nk
    cells = (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    coords = body.split()[:3 * nodes]
    rgas = RU * sum(y / weight[s] for s, y in STEAM.items())
    state = [LEFT if c < cells // 2 else RIGHT for c in range(cells)]
    rho = [p / (rgas * t) for p, t in state]
    nvar = 3 + len(species) + 4
    names = " ".join('"rho%d"' % (s + 1) for s in range(len(species)))
    zone = re.sub(r"\[4-\d+\]=CELLCENTERED", "[4-%d]=CELLCENTERED" % nvar, zone)
    with open(os.path.join("INPUT", "ic.tec"), "w") as out:
        out.write(' VARIABLES ="x" "y" "z"  %s "u" "v" "w" "p"\n' % names)
        out.write(zone + "\n")
        out.write("\n".join(" " + v for v in coords) + "\n")
        for s in species:
            out.write("".join(" %.15E\n" % (STEAM.get(s, 0.0) * r) for r in rho))
        for _ in range(3):
            out.write("".join(" %.15E\n" % 0.0 for _ in range(cells)))
        out.write("".join(" %.15E\n" % p for p, _ in state))


def mixture(src, out):
    header, zone, nodes, cells, vals = read_tube(src)
    xyz, rho, uvwp = vals[:3 * nodes], vals[3 * nodes:3 * nodes + cells], vals[3 * nodes + cells:]
    species = [['%.15E' % (float(r) * y) for r in rho] for y in Y_AIR]
    zone = zone.replace('[4-8]=CELLCENTERED', '[4-11]=CELLCENTERED')
    assert '[4-11]=CELLCENTERED' in zone, zone
    with open(out, 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" "rho1" "rho2" "rho3" "rho4" "u" "v" "w" "p"\n')
        f.write(zone + '\n')
        for v in xyz + [s for block in species for s in block] + uvwp:
            f.write(' %s\n' % v)


def tables(text):
    """Zones of a table of test/1D/premixed-flame: name -> lines from its ZONE line on; and the lines before them."""
    lines = text.split('\n')
    idx = [i for i, l in enumerate(lines) if l.startswith('ZONE')] + [len(lines)]
    zones = {}
    for a, b in zip(idx[:-1], idx[1:]):
        zones[re.search(r'T="([^"]+)"', lines[a]).group(1)] = [l for l in lines[a:b] if l]
    return lines[:idx[0]], zones


def soot(src, sod_ic, dest, kind):
    with open(dest + '/phase.txt', 'w') as f:
        f.write('ideal-gas phase\n' + ''.join('%s %.6f\n' % (n, w) for n, w, _ in SOOT_SPECIES))
    for name in ('thermo.dat', 'transport.dat'):
        with open(src + '/' + name) as f:
            head, zones = tables(f.read())
        body = []
        for n, _, s in SOOT_SPECIES:
            rows = list(zones[s])
            rows[0] = rows[0].replace('T="%s"' % s, 'T="%s"' % n)
            body += rows
        with open(dest + '/' + name, 'w') as f:
            f.write('\n'.join(head + body) + '\n')
    with open(sod_ic) as f:
        lines = f.read().split('\n')
    zone = lines[1]
    ni, nj, nk = (int(re.search(r'%s=\s*(\d+)' % c, zone).group(1)) for c in 'IJK')
    nodes, cells = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    xyz = ' '.join(lines[2:]).split()[:3 * nodes]
    cols = [RHO_S * y for y in Y_SOOT] + [U_S, 0.0, 0.0, P_S]
    names = ' '.join('"rho%d"' % (i + 1) for i in range(len(SOOT_SPECIES))) + ' "u" "v" "w" "p"'
    if kind in ('soot', 'sa'):
        cols.append(NP_S)
        names += ' "rho*Np"'
    if kind in ('sa', 'sa-plain'):
        cols.append(3.0 * MU_S)
        names += ' "rhoNuTilde"'
    zone = zone.replace('[4-8]=CELLCENTERED', '[4-%d]=CELLCENTERED' % (3 + len(cols)))
    with open(dest + '/ic.tec', 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" %s\n' % names)
        f.write(zone + '\n')
        for v in xyz:
            f.write(' %s\n' % v)
        for c in cols:
            f.write((' %.15E\n' % c) * cells)


def read_diagnostic(path):
    """Names, declared variables, the number of variables written and the cell values of a one-zone ASCII file."""
    with open(path) as f:
        lines = f.read().split('\n')
    iv = [i for i, l in enumerate(lines) if 'VARIABLES' in l][0]
    iz = [i for i, l in enumerate(lines) if 'ZONE' in l][0]
    names = [a or b for a, b in re.findall(r'"([^"]*)"|(\S+)', ' '.join(lines[iv:iz]).split('=', 1)[1])]
    zone = lines[iz]
    ni, nj, nk = (int(re.search(r'%s=\s*(\d+)' % c, zone).group(1)) for c in 'IJK')
    nodes, cells = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    declared = int(re.search(r'\[4-(\d+)\]=CELLCENTERED', zone).group(1))
    data = [l for l in lines[iz + 1:] if l.strip() and l.lstrip()[0] in '+-.0123456789']
    vals = [float(v.replace('D', 'E')) for v in ' '.join(data).split()]
    written = 3 + (len(vals) - 3 * nodes) // cells
    cols = {n: vals[3 * nodes + (k - 3) * cells: 3 * nodes + (k - 2) * cells] for k, n in enumerate(names) if k >= 3}
    return names, declared, written, cols


def diagnostic_check(diag, ref, expected):
    try:
        names, declared, written, cols = read_diagnostic(diag)
        _, _, _, ref_cols = read_diagnostic(ref)
    except (OSError, ValueError, IndexError, AttributeError) as err:
        print('    no diagnostic to read (%s)' % err)
        return 1
    want = ['x', 'y', 'z'] + expected
    dt, dt_ref = cols.get('dt', []), ref_cols.get('dt', [])
    err = max(abs(a - b) / abs(b) for a, b in zip(dt, dt_ref)) if dt and len(dt) == len(dt_ref) else float('inf')
    print('    %d names, %d variables declared, %d written: %s; column dt against the run without soot: %.1e'
          % (len(names), declared, written, ' '.join(names[3:]), err))
    return 0 if names == want and declared == written == len(names) and err <= 1.0e-6 else 1


D410_N, D410_NK = 50, 2                   # cells of a block along x, cells along z
D410_DX, D410_LY, D410_LZ = 0.01, 0.01, 1.0e-3
D410_T = 300.0
D410_P = {1: 1.0e5, 2: 1.02e5}            # pressure of the files for each block; the cells are at 1e5 Pa
D410_TIMES = (0.0, 1.0e-3)
D410_FORMS = ('valid', 'swapped', 'no-block-2', 'no-time-2', 'short-zone', 'time-2-other-block', 'moved-5e-5', 'moved-2e-4',
              'digits-15')


def fr_p(v, digits=17):
    """A real as the format (E25.17E3) writes it, 0.ddd...E+eee; with digits = 15 the same form with 15 significant digits."""
    if v == 0.0:
        m, e = '0' * digits, 0
    else:
        mant, exp = ('%.*e' % (digits - 1, abs(v))).split('e')
        m, e = mant.replace('.', ''), int(exp) + 1
    return ('%s0.%sE%s%03d' % ('-' if v < 0.0 else '', m, '+' if e >= 0 else '-', abs(e))).rjust(25)


def duct410_nodes(b, face, dz=0.0):
    """(I, J) and the nodes of face `face` of block b, the first index fastest: (j, k) on faces 1 and 2, (i, k) on 3 and 4,
    (i, j) on 5 and 6."""
    x0 = 0.5 * (b - 1)
    xs = [x0 + D410_DX * i for i in range(D410_N + 1)]
    zs = [D410_LZ * (k / D410_NK - 0.5) for k in range(D410_NK + 1)]
    if face in (1, 2):
        xf = x0 + (0.0 if face == 1 else D410_N * D410_DX)
        dims, nodes = (2, D410_NK + 1), [(xf, D410_LY * j, z) for z in zs for j in range(2)]
    elif face in (3, 4):
        yf = 0.0 if face == 3 else D410_LY
        dims, nodes = (D410_N + 1, D410_NK + 1), [(x, yf, z) for z in zs for x in xs]
    else:
        zf = zs[0] if face == 5 else zs[-1]
        dims, nodes = (D410_N + 1, 2), [(x, D410_LY * j, zf) for j in range(2) for x in xs]
    return dims, [(x, y, z + dz) for (x, y, z) in nodes]


def duct410_zone(title, time, dims, nodes, p, digits, ncells=None):
    """A zone as the boundary-condition builder of ATLAS writes it; ncells keeps the first ncells cells along I."""
    ni, nj = dims
    if ncells is not None:
        nodes = [nd for n, nd in enumerate(nodes) if n % ni <= ncells]
        ni = ncells + 1
    out = [' ZONE T="%s", I=%d, J=%d, K=1, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), '
           'SOLUTIONTIME=%15.8E\n' % (title, ni, nj, time)]
    for c in range(3):
        out += [fr_p(nd[c], digits) + '\n' for nd in nodes]
    ncell = (ni - 1) * (nj - 1)
    for v in (p / (RU / 28.970418 * D410_T), 0.0, 0.0, 0.0, p):
        out += [fr_p(v) + '\n'] * ncell
    return ''.join(out)


def duct410(form):
    if form not in D410_FORMS:
        sys.exit('unknown form %s: %s' % (form, ', '.join(D410_FORMS)))
    files = {'q2d-ends.dat': {1: 1, 2: 2}, 'q2d-jfaces.dat': {1: 4, 2: 3}, 'q2d-side.dat': {1: 5, 2: 6}}
    digits = 15 if form == 'digits-15' else 17
    for name, faces in files.items():
        side = name == 'q2d-side.dat'
        dz = {'moved-5e-5': 5.0e-5, 'moved-2e-4': 2.0e-4}.get(form, 0.0) if side else 0.0
        zones = []
        for it, t in enumerate(D410_TIMES, start=1):
            for b in (1, 2):
                src, ncells = b, None
                if side:
                    if form == 'swapped':
                        src = 3 - b
                    elif (form == 'no-block-2' and b == 2) or (form == 'no-time-2' and b == 2 and it == 2):
                        continue
                    elif form == 'short-zone' and b == 2:
                        ncells = 25
                    elif form == 'time-2-other-block' and b == 1 and it == 2:
                        src = 2
                dims, nodes = duct410_nodes(src, faces[src], dz)
                zones.append(duct410_zone('B%d_T%d' % (b, it), t, dims, nodes, D410_P[src], digits, ncells))
        with open(name, 'w') as f:
            f.write(' VARIABLES ="x" "y" "z" "rho(1)" "u" "v" "w" "p"\n' + ''.join(zones))
    rho = 1.0e5 / (RU / 28.970418 * D410_T)
    with open(os.path.join('INPUT', 'ic.tec'), 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b in (1, 2):
            f.write(' ZONE  T = B%d-IG, I=%d, J=2, K=%d, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), '
                    'SOLUTIONTIME=0.0\n' % (b, D410_N + 1, D410_NK + 1))
            for c in range(3):
                for k in range(D410_NK + 1):
                    for j in range(2):
                        for i in range(D410_N + 1):
                            node = (0.5 * (b - 1) + D410_DX * i, D410_LY * j, D410_LZ * (k / D410_NK - 0.5))
                            f.write(' %.16E\n' % node[c])
            for v in (rho, 0.0, 0.0, 0.0, 1.0e5):
                f.write((' %.16E\n' % v) * (D410_N * D410_NK))
    q2d = {(1, 1): 'q2d-ends.dat', (2, 2): 'q2d-ends.dat', (1, 4): 'q2d-jfaces.dat', (2, 3): 'q2d-jfaces.dat',
           (1, 5): 'q2d-side.dat', (2, 6): 'q2d-side.dat'}
    rows = []
    for b in (1, 2):
        for face in range(1, 7):
            ii = [1] if face == 1 else [D410_N] if face == 2 else range(1, D410_N + 1)
            kk = [1] if face == 5 else [D410_NK] if face == 6 else range(1, D410_NK + 1)
            for k in kk:
                for i in ii:
                    head = ''.join('%8d' % v for v in (b, i, 1, k, face))
                    if (b, face) in q2d:
                        rows.append(head + '%8d\n%s\n' % (410, q2d[(b, face)]))
                    elif (b, face) in ((1, 2), (2, 1)):
                        other = (2, 1, 1, k, 1) if b == 1 else (1, D410_N, 1, k, 2)
                        rows.append(head + '%8d\n' % 101 + ''.join('%8d' % v for v in other + (1, 0, 0, 1)) + '\n')
                    else:
                        rows.append(head + '%8d\n' % 0)
    with open(os.path.join('INPUT', 'bc.txt'), 'w') as f:
        f.write(''.join(rows))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'sod':
        sod(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == 'duct410':
        duct410(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == 'nan':
        nan(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == 'inf':
        inf(sys.argv[2])
    elif len(sys.argv) == 4 and sys.argv[1] == 'steam':
        steam(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == 'mixture':
        mixture(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 6 and sys.argv[1] == 'soot':
        soot(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif len(sys.argv) > 4 and sys.argv[1] == 'diagnostic-check':
        sys.exit(diagnostic_check(sys.argv[2], sys.argv[3], sys.argv[4:]))
    else:
        sys.exit(__doc__)
