"""The Sod tube moved by a uniform velocity, and its mirror image, for the two outer branches of HLLC (RiemannSmoke).

  python3 -B mirror_sod.py case SIGN V   INPUT/ic.tec (a file, not a link) from ../../1D/Sod79/INPUT/ic_x1.tec:
                                         SIGN = plus, the Sod tube with u = +V in every cell: every wave runs to the right
                                         and the HLLC flux of every face is that of its left state (S_L >= 0);
                                         SIGN = minus, its mirror image x -> 1 - x with u = -V: every wave runs to the left
                                         and the flux is that of the right state (S_R <= 0)
  python3 -B mirror_sod.py compare A B   the field of B must be the mirror image of the field of A: density and pressure
                                         in the reverse order of the cells, u reversed and negated, within 1e-12 of their
                                         largest values (exit status 0 when it is)

With V = 3 the fastest wave of the Sod tube, about 1.75, cannot turn any face back: S_L > 0 (plus) and S_R < 0 (minus)
everywhere for the whole run. The scheme has no preferred direction, so the two runs are mirror images to round-off.
"""
import os
import re
import sys

SOD = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '1D', 'Sod79', 'INPUT', 'ic_x1.tec')
TOL = 1.0e-12


def read(path):
    """Variable names, the number of nodes and of cells, and the values of each variable, of a one-zone ASCII Tecplot
    file with BLOCK packing (nodal coordinates, cell-centred variables)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]   # quoted or not
    zone = lines[1]
    ni, nj, nk = (int(re.search(k + r'\s*=\s*(\d+)', zone).group(1)) for k in ('I', 'J', 'K'))
    vals = ' '.join(lines[2:]).split()
    nn, nc = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    out, pos = {}, 0
    for q, name in enumerate(names):
        cnt = nn if q < 3 else nc
        out[name] = [float(v) for v in vals[pos:pos + cnt]]
        pos += cnt
    return names, lines[:2], out


def case(sign, speed):
    names, head, val = read(SOD)
    rho, p = val['rho1'], val['p']
    if sign == 'minus':
        rho, p, u = rho[::-1], p[::-1], [-speed] * len(rho)
    else:
        u = [speed] * len(rho)
    val['rho1'], val['u'], val['p'] = rho, u, p
    os.makedirs('INPUT', exist_ok=True)
    if os.path.lexists('INPUT/ic.tec'):
        os.remove('INPUT/ic.tec')                 # a link to the Sod tube: never write through it
    with open('INPUT/ic.tec', 'w') as out:
        out.write('\n'.join(head) + '\n')
        for name in names:
            out.write(''.join(f' {v:.15E}\n' for v in val[name]))


def compare(a_path, b_path):
    names, _, a = read(a_path)
    _, _, b = read(b_path)
    ok = True
    for name, flip in ((names[3], 1.0), ('p', 1.0), ('u', -1.0)):     # names[3]: the density, rho1 or rho(1)
        va, vb = a[name], b[name][::-1]
        scale = max(abs(v) for v in va)
        err = max(abs(x - flip * y) for x, y in zip(va, vb)) / scale
        good = err <= TOL
        ok = ok and good
        print(f'    {name:6s} of the mirror run against the reversed first run: {err:.3e} of the largest  '
              f'[{"ok  " if good else "FAIL"}]')
    return ok


if __name__ == '__main__':
    if sys.argv[1] == 'case':
        case(sys.argv[2], float(sys.argv[3]))
    elif sys.argv[1] == 'compare':
        sys.exit(0 if compare(sys.argv[2], sys.argv[3]) else 1)
    else:
        sys.exit('unknown request ' + sys.argv[1])
