#!/usr/bin/env python3
"""Replace the first boundary record of a given face in bc.txt with a record of another type (set-up rejection tests
of the multirate integrator). The run directory gets its own copy of bc.txt instead of the link to the case.
Usage: set_bc.py RUN_DIR FACE TYPE 'data line' ['data line' ...]
Records in bc.txt: one line 'b i j k face type', followed by the data lines of the types that have them (101, 201, 103:
one line; 102, 401-408, 410, 420, 501-506: see src/lib/io/IO_BC.f90)."""
import os
import sys

WITH_DATA = {101, 201, 103, 401, 402, 403, 404, 405, 406, 407, 408, 410, 420, 501, 502, 503, 504, 505, 506}


def main():
    run, face, typ, data = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4:]
    path = os.path.join(run, 'INPUT', 'bc.txt')
    lines = open(path).read().splitlines()
    k = 0
    while k < len(lines):
        s = lines[k].split()
        if len(s) == 6:
            t = int(s[5])
            if int(s[4]) == face:
                s[5] = str(typ)
                old_data = 1 if t in WITH_DATA else 0
                new = [' '.join('%9s' % x for x in s)] + data
                lines[k:k + 1 + old_data] = new
                break
            k += 2 if t in WITH_DATA else 1
        else:
            k += 1
    else:
        sys.exit('no record on face %d' % face)
    os.remove(path)                 # the link to the case
    open(path, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
