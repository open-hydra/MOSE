#!/usr/bin/env python3
"""One block against three blocks: the same mesh and state in one block (RUN_1B) and
in three blocks whose connections lie inside bands of equal level (RUN_3B), so the level maps are the same. Largest
difference over the cells, max |q_3b - q_1b| / max |q_1b| for rho, u, p (cells matched by their centres), for the
multirate pair and for the single-rate pair. Criterion: multirate <= 1e-11 and <= 10 x the single-rate difference (the
single-rate flux of a connection takes its operations in a different order than an internal face, so the single-rate
pair is not identical).
Usage: verify_one_block.py MR_1B MR_3B SR_1B SR_3B"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402


def diff(r1, r3):
    b1, b3 = mc.read_run(r1), mc.read_run(r3)
    out = {}
    for v in ('rho(1)', 'u', 'p'):
        x1, x3 = mc.by_centre(b1, v), mc.by_centre(b3, v)
        if set(x1) != set(x3):
            raise SystemExit('the meshes of %s and %s differ' % (r1, r3))
        d = max(abs(x3[c] - x1[c]) for c in x1)
        out[v] = d / max(abs(x) for x in x1.values())
    return out


def main():
    mr1, mr3, sr1, sr3 = sys.argv[1:5]
    dm, ds = diff(mr1, mr3), diff(sr1, sr3)
    print('  multirate,   one vs three blocks: ' + ', '.join('%s %.3e' % kv for kv in dm.items()))
    print('  single rate, one vs three blocks: ' + ', '.join('%s %.3e' % kv for kv in ds.items()))
    m, s = max(dm.values()), max(ds.values())
    ok = m <= 1e-11 and m <= 10.0 * s
    print('  largest: multirate %.3e (<= 1e-11 and <= 10 x %.3e) -> %s' % (m, s, 'ok' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
