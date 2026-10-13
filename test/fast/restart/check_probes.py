"""Probe rows at a restart (part C of the Restart test).

Compares the probe files of a run split by a restart with those of the same
run in one go. A restarted run numbers its iterations from 1 again, so the
rows written after the restart (those after the first PRE rows of a split
file) stand for the iterations OFFSET + n of the run in one go. Every row
must be found there, with the same values to within 1e-10 of the largest
value of each column (solution-file round-trip precision), and the probes
must read the flow: no column may be zero in every row.

Usage: python check_probes.py OFFSET <full> <split> PRE [<full> <split> PRE ...]
Exit 0 = pass, 1 = fail.
"""
import sys

TOL = 1.0e-10


def rows(path):
    out = []
    with open(path) as f:
        for line in f:
            v = [float(x) for x in line.split()]
            if v:
                out.append((int(round(v[0])), v[1:]))
    return out


offset = int(sys.argv[1])
args = sys.argv[2:]
worst, status = 0.0, 0
for i in range(0, len(args), 3):
    full_path, split_path, pre = args[i], args[i + 1], int(args[i + 2])
    try:
        full = dict(rows(full_path))
        split = rows(split_path)
    except OSError as err:
        print("    no probe file to read (%s)" % err)
        status = 1
        continue
    keyed = [(n if k < pre else n + offset, v) for k, (n, v) in enumerate(split)]
    missing = [n for n, _ in keyed if n not in full]
    if not split or missing:
        print("    %s: %d rows, %d not in the run in one go" % (split_path, len(split), len(missing)))
        status = 1
        continue
    ncol = len(next(iter(full.values())))
    scale = [max(abs(r[c]) for r in full.values()) for c in range(ncol)]
    if min(scale) == 0.0:
        print("    %s: a column is zero in every row: the probe does not read the flow" % full_path)
        status = 1
        continue
    dev = max(abs(v[c] - full[n][c]) / scale[c] for n, v in keyed for c in range(ncol))
    worst = max(worst, dev)
    if dev > TOL:
        status = 1
print("    probe rows against the run in one go: largest deviation %.1e of the largest value" % worst)
sys.exit(status)
