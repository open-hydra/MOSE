"""Inputs of the ExitCodes fast test.

Every case starts from the 100-cell Sod tube of test/1D/Sod79 (INPUT/ic_x1.tec, INPUT/bc_x1.txt); run.sh breaks it in one way,
mostly through input.ini, so the inputs here are the Sod ones, or the Sod ones changed in one value.

Usage (from test/fast/input-checks/exit-codes):
  python3 -B make_case.py sod SOD_INPUT_DIR
      writes INPUT/ic.tec and INPUT/bc.txt, copies of the Sod ones.
  python3 -B make_case.py nan SOD_INPUT_DIR
      the same, with the velocity u of cell 50 written as NaN in INPUT/ic.tec.
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


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'sod':
        sod(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == 'nan':
        nan(sys.argv[2])
    elif len(sys.argv) == 4 and sys.argv[1] == 'mixture':
        mixture(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 6 and sys.argv[1] == 'soot':
        soot(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif len(sys.argv) > 4 and sys.argv[1] == 'diagnostic-check':
        sys.exit(diagnostic_check(sys.argv[2], sys.argv[3], sys.argv[4:]))
    else:
        sys.exit(__doc__)
