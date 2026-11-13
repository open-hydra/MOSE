"""Exact restart check for the vtk raw format (I3/F6).

Compares, file by file and byte by byte, the OUTPUT/ tree of an uninterrupted
run with the OUTPUT/ tree of a run that stopped half way and restarted
(newrun=false). vtk raw stores every value with all its bits, so the two trees
must be identical: the solutions written along the way (numbered files), the
final solution with its TIME, and residual-history.dat, whose first column is
the global iteration counter. Logs are not compared.

Usage: python check_restart_vtk.py <OUTPUT of the full run> <OUTPUT of the restarted run>
Exit 0 = pass, 1 = fail.
"""
import filecmp
import sys
from pathlib import Path

full, rest = Path(sys.argv[1]), Path(sys.argv[2])


def tree(root):
    return {p.relative_to(root).as_posix() for p in root.rglob("*")
            if p.is_file() and not p.name.endswith(".log")}


a, b = tree(full), tree(rest)
if not a:
    print(f"Restart (vtk raw): FAIL — no output in {full}")
    sys.exit(1)
fail = False
for name in sorted(a - b):
    print(f"  missing after the restart: {name}")
    fail = True
for name in sorted(b - a):
    print(f"  only after the restart:    {name}")
    fail = True
for name in sorted(a & b):
    if not filecmp.cmp(full / name, rest / name, shallow=False):
        print(f"  differs:                   {name}")
        fail = True
if fail:
    print("Restart (vtk raw): FAIL — the restarted run does not reproduce the uninterrupted one bit for bit")
    sys.exit(1)
print(f"Restart (vtk raw): ok — {len(a)} output files identical bit for bit")
sys.exit(0)
