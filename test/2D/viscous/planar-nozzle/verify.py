"""
Planar nozzle on grid families (O, H, ZZ, SB2), SA and SST: short form.

The O and H grids are mirror images of themselves about the centre line, but
their four blocks have different index orientations (the lower convergent
block is reversed in i and j, the upper divergent one turned by 90 degrees).
A discretisation that does not depend on the orientation of the indices gives
a solution that is mirror-symmetric at every iteration: rho, u, p, T and the
turbulence variables even in y, v odd.  The short form (./MOSE.sh testfast)
runs 300 iterations on the coarse level from the quasi-1D field and checks
  1. on O and H, with SA and SST: the largest mirror asymmetry of every
     variable, relative to its largest magnitude, is below SYM_TOL;
  2. on ZZ, with SA: the field is finite, with positive density and pressure.
The full form (./MOSE.sh test) runs every family on three levels to
convergence; its analysis (orders, limits of the families, quasi-1D
comparison) and figures belong to the V&V page.

Usage:
    python verify.py            # report, exit 0/1

Exit 0 = pass, 1 = fail.
"""
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

SYM_TOL = 1.0e-10               # largest relative mirror asymmetry allowed
RT = 0.01                       # throat half-height [m]

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError("ORION package not found. Run 'conda activate base' or install ORION.")


def centres(a):
    return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] + a[:-1, 1:, :-1] + a[1:, 1:, :-1] +
                    a[:-1, :-1, 1:] + a[1:, :-1, 1:] + a[:-1, 1:, 1:] + a[1:, 1:, 1:])


def read_cells(path):
    """Cell centres (x, y) of all blocks and the cell variables, flattened in the same order."""
    xb, yb, zb, vb, names = read_TEC(path)
    C, V = [], {n: [] for n in names[3:]}
    for x, y, var in zip(xb, yb, vb):
        C.append(np.stack([centres(x).ravel(), centres(y).ravel()], axis=-1))
        for k, n in enumerate(names[3:]):
            V[n].append(np.asarray(var[k]).ravel())
    return np.concatenate(C), {n: np.concatenate(v) for n, v in V.items()}


def mirror_asymmetry(path):
    C, V = read_cells(path)
    d, j = cKDTree(C).query(np.stack([C[:, 0], -C[:, 1]], axis=-1))
    if d.max() > 1e-9 * RT:
        return None
    out = {}
    for n, a in V.items():
        scale = float(np.max(np.abs(a)))
        if scale == 0.0 or n in ("g", "R"):
            continue
        b = -a[j] if n == "v" else a[j]
        out[n] = float(np.max(np.abs(a - b))) / scale
    return out


form = open("OUTPUT/form.txt").read().strip() if os.path.exists("OUTPUT/form.txt") else "unknown"
if form != "testfast":
    print("planar nozzle: this verify checks the short form; the full form is analysed for the V&V page")
ok = True
print("Planar nozzle, short form: mirror symmetry after 300 iterations (largest relative asymmetry)")
for fam in ("O", "H"):
    for model in ("SA", "SST"):
        f = "OUTPUT/field_%s-L0-%s.tec" % (fam, model)
        if not os.path.exists(f):
            print("  %s %s: FAIL - %s missing (run './MOSE.sh testfast')" % (fam, model, f))
            ok = False
            continue
        a = mirror_asymmetry(f)
        if a is None:
            print("  %s %s: FAIL - the grid is not mirror-symmetric" % (fam, model))
            ok = False
            continue
        worst = max(a.values())
        print("  %-2s %-3s: %s  -> %s" % (fam, model, "  ".join("%s %.1e" % (n, v) for n, v in a.items()),
                                          "PASS" if worst < SYM_TOL else "FAIL"))
        ok &= worst < SYM_TOL
f = "OUTPUT/field_ZZ-L0-SA.tec"
if os.path.exists(f):
    _, V = read_cells(f)
    finite = all(np.all(np.isfinite(v)) for v in V.values())
    positive = bool(np.all(V["rho(1)"] > 0.0) and np.all(V["p"] > 0.0))
    print("  ZZ SA: finite %s, positive density and pressure %s -> %s" % (finite, positive, "PASS" if finite and positive else "FAIL"))
    ok &= finite and positive
else:
    print("  ZZ SA: FAIL - %s missing" % f)
    ok = False
status = "PASS" if ok else "FAIL"
print("  STATUS : %s (tolerance %.0e)" % (status, SYM_TOL))
sys.exit(0 if ok else 1)
