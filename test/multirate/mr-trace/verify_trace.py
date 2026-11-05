#!/usr/bin/env python3
"""Stage sequence of the first macro step of MOSE's multirate integrator (lines 'MR-TRACE l=.. k=.. t=a/b theta=m:c/d
...' of the solver output, [MOSE-Multirate] trace = .true.) against the reference sequence of the recursion
(bench_sequence_L3.txt): equal element by element.
Usage: verify_trace.py SOLVE_LOG BENCH_FILE"""
import re
import sys
from fractions import Fraction


def frac(s):
    if '/' in s:
        a, b = s.split('/')
        return Fraction(int(a), int(b))
    return Fraction(int(s))


def parse(lines, prefix):
    rx = re.compile(r'^\s*' + prefix + r'l=(\d+)\s+k=(\d+)\s+t=(\S+)\s+theta=(.*)$')
    out = []
    for ln in lines:
        m = rx.match(ln.rstrip('\n'))
        if not m:
            continue
        th = tuple((int(tok.split(':')[0]), frac(tok.split(':')[1])) for tok in m.group(4).split())
        out.append((int(m.group(1)), int(m.group(2)), frac(m.group(3)), th))
    return out


def main():
    got = parse(open(sys.argv[1], errors='replace'), r'MR-TRACE\s+')
    ref = parse(open(sys.argv[2]), '')
    ok = got == ref and len(ref) > 0
    print('  bench: %d stages; MOSE: %d stages; %s' % (len(ref), len(got), 'EQUAL' if ok else 'DIFFERENT'))
    if not ok:
        for i in range(max(len(ref), len(got))):
            r = ref[i] if i < len(ref) else None
            g = got[i] if i < len(got) else None
            if r != g:
                print('  %3d bench %s | MOSE %s' % (i, r, g))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
