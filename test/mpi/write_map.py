#!/usr/bin/env python3
"""Write INPUT/decomposition.map, the block-to-rank assignment, in the format of the ATLAS mesh decomposer (MDB), for the
cases of make_case.py (PI x PJ blocks of NI x NJ x 1 cells, block b = 1 + p + PI*q).

MDB writes four comment lines, one line per block with the 13 integers
    new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank
(Fortran format 2I6, 3(2I6,2X), 3X, 3I5, I13, I6), then a '#' line and the face origin of every block (format 2I6, 6A9),
which MOSE does not read. Here no block is split: every block is its own parent and covers its whole index range.

--owners gives the rank of every block:
    lpt:N           the assignment MOSE makes without a map on N ranks: blocks by decreasing cell count (ties in block
                    order), each to the least loaded rank (the lowest one on ties);
    cyclic:N        block b on rank (b - 1) mod N;
    reverse:N       block b on rank N - 1 - ((b - 1) mod N);
    list:R1,R2,...  one rank per block.
--bad writes a wrong map, to check the errors of the reader:
    missing-row  the last block line left out;      extra-row   one block line more;
    order        the first two block lines swapped;  dims        ni of block 2 one larger (cells as in the mesh);
    cells        cells of block 2 one larger;        no-rank     12 integers per line, as MDB writes without ranks.

Usage: write_map.py --layout PIxPJ --cells NI NJ --owners SPEC [--bad KIND] [--out PATH]
"""

import argparse
import os


def parse():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--layout', required=True, help='PIxPJ blocks, e.g. 4x1')
    ap.add_argument('--cells', nargs=2, type=int, required=True, metavar=('NI', 'NJ'), help='cells per block')
    ap.add_argument('--owners', required=True, help='lpt:N, cyclic:N, reverse:N or list:R1,R2,...')
    ap.add_argument('--bad', choices=['missing-row', 'extra-row', 'order', 'dims', 'cells', 'no-rank'])
    ap.add_argument('--out', default=os.path.join('INPUT', 'decomposition.map'))
    a = ap.parse_args()
    a.PI, a.PJ = (int(s) for s in a.layout.lower().split('x'))
    return a


def owners(spec, sizes):
    kind, _, arg = spec.partition(':')
    nb = len(sizes)
    if kind == 'list':
        ranks = [int(r) for r in arg.split(',')]
        assert len(ranks) == nb, f'{len(ranks)} ranks for {nb} blocks'
        return ranks
    n = int(arg)
    if kind == 'cyclic':
        return [b % n for b in range(nb)]
    if kind == 'reverse':
        return [n - 1 - b % n for b in range(nb)]
    assert kind == 'lpt', spec
    order = sorted(range(nb), key=lambda b: -sizes[b])       # stable: ties keep block order
    load, out = [0] * n, [0] * nb
    for b in order:
        r = load.index(min(load))                             # lowest rank among the least loaded
        out[b] = r
        load[r] += sizes[b]
    return out


def block_line(new, parent, dims, rank):
    """One block line in MDB's format; rank None gives the 12-integer line MDB writes without ranks."""
    ni, nj, nk = dims
    s = f'{new:6d}{parent:6d}'
    for n in dims:
        s += f'{1:6d}{n:6d}  '
    s += f'   {ni:5d}{nj:5d}{nk:5d}{ni * nj * nk:13d}'
    return s + ('' if rank is None else f'{rank:6d}')


def main():
    a = parse()
    nb = a.PI * a.PJ
    dims = [(a.cells[0], a.cells[1], 1)] * nb
    rank = owners(a.owners, [d[0] * d[1] * d[2] for d in dims])
    lines = [block_line(b + 1, b + 1, dims[b], None if a.bad == 'no-rank' else rank[b]) for b in range(nb)]
    if a.bad == 'missing-row':
        lines.pop()
    elif a.bad == 'extra-row':
        lines.append(block_line(nb + 1, nb, dims[-1], rank[-1]))
    elif a.bad == 'order':
        lines[0], lines[1] = lines[1], lines[0]
    elif a.bad == 'dims':
        ni, nj, nk = dims[1]
        lines[1] = block_line(2, 2, (ni + 1, nj, nk), rank[1])
        lines[1] = lines[1][:-19] + f'{ni * nj * nk:13d}{rank[1]:6d}'
    elif a.bad == 'cells':
        lines[1] = lines[1][:-19] + f'{int(lines[1][-19:-6]) + 1:13d}{rank[1]:6d}'
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as f:
        f.write('# ATLAS MDB decomposition map\n')
        f.write(f'# parent blocks: {nb}  ->  new blocks: {nb}\n')
        f.write('#\n')
        f.write('# new  parent   i0   i1    j0   j1    k0   k1      ni   nj   nk        cells  rank\n')
        f.write(''.join(line + '\n' for line in lines))
        f.write('#\n')
        f.write('# face origin of every new block (parent face number, or cut)\n')
        f.write('# new  parent   face1    face2    face3    face4    face5    face6\n')
        for b in range(nb):
            f.write(f'{b + 1:6d}{b + 1:6d}' + ''.join(f'{f:>9d}' for f in range(1, 7)) + '\n')
    print(f'{a.out}: {len(lines)} block lines, ranks {",".join(str(r) for r in rank)}'
          f'{" (" + a.bad + ")" if a.bad else ""}')


if __name__ == '__main__':
    main()
