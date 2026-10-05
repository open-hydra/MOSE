"""Inputs of the ExitCodes fast test.

Usage (from test/fast/input-checks/exit-codes):
  python3 -B make_case.py mixture SOD_IC OUT
      writes OUT: the 100-cell Sod tube of SOD_IC (one species, rho1) as a
      mixture of the four species of test/1D/diff1 (order of its
      INPUT/phase.txt: CH4, H2O, N2, O2) with the mass fractions of air
      (0, 0, 0.767, 0.233); density, velocity and pressure are those of the tube.
"""
import re
import sys

Y_AIR = (0.0, 0.0, 0.767, 0.233)


def mixture(src, out):
    with open(src) as f:
        lines = f.read().split('\n')
    header, zone = lines[0], lines[1]
    names = re.findall(r'"([^"]+)"', header)
    assert names == ['x', 'y', 'z', 'rho1', 'u', 'v', 'w', 'p'], names
    ni, nj, nk = (int(re.search(r'%s=\s*(\d+)' % c, zone).group(1)) for c in 'IJK')
    nodes, cells = ni * nj * nk, (ni - 1) * max(nj - 1, 1) * max(nk - 1, 1)
    vals = ' '.join(lines[2:]).split()
    assert len(vals) == 3 * nodes + 5 * cells, (len(vals), nodes, cells)
    xyz, rho, uvwp = vals[:3 * nodes], vals[3 * nodes:3 * nodes + cells], vals[3 * nodes + cells:]
    species = [['%.15E' % (float(r) * y) for r in rho] for y in Y_AIR]
    zone = zone.replace('[4-8]=CELLCENTERED', '[4-11]=CELLCENTERED')
    assert '[4-11]=CELLCENTERED' in zone, zone
    with open(out, 'w') as f:
        f.write(' VARIABLES ="x" "y" "z" "rho1" "rho2" "rho3" "rho4" "u" "v" "w" "p"\n')
        f.write(zone + '\n')
        for v in xyz + [s for block in species for s in block] + uvwp:
            f.write(' %s\n' % v)


if __name__ == '__main__':
    if len(sys.argv) != 4 or sys.argv[1] != 'mixture':
        sys.exit(__doc__)
    mixture(sys.argv[2], sys.argv[3])
