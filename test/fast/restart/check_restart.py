"""Restart round-trip check (I3/F6).

Compares the field from an uninterrupted N-iteration run against the field from
a run that stopped at N/2, wrote a solution, restarted (newrun=false) and ran
the remaining N/2 iterations.  They must agree to within solution-file (Tecplot
ASCII) round-trip precision.

Usage: python check_restart.py <field_full.tec> <field_restart.tec>
Exit 0 = pass, 1 = fail.
"""
import sys
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or install ORION."
    )

TOL = 1.0e-10  # relative tolerance (Tecplot ASCII keeps ~15 significant digits)

full = read_TEC(sys.argv[1])[3][0]
rest = read_TEC(sys.argv[2])[3][0]

# Normalize the difference by a single global field scale (dominated by rho/p),
# not by each variable's own magnitude — otherwise the ~zero transverse
# velocities (v, w in a 1-D problem) divide noise by noise.
names = ["rho", "u", "v", "w", "p"]
global_ref = 0.0
worst_abs = 0.0
worst_var = None
for k, name in enumerate(names):
    a = np.asarray(full[k][:, :, :], dtype=float)
    b = np.asarray(rest[k][:, :, :], dtype=float)
    if not np.all(np.isfinite(b)):
        print(f"Restart: FAIL — non-finite {name} in restarted field")
        sys.exit(1)
    global_ref = max(global_ref, float(np.abs(a).max()))
    d = float(np.abs(a - b).max())
    if d > worst_abs:
        worst_abs, worst_var = d, name

worst_rel = worst_abs / max(global_ref, 1e-30)
if worst_rel > TOL:
    print(f"Restart: FAIL — max difference {worst_abs:.3e} ({worst_rel:.3e} "
          f"relative) in {worst_var} (tol {TOL:.0e})")
    sys.exit(1)

print(f"Restart: ok — full vs split+restart agree to {worst_rel:.3e} relative "
      f"(max abs {worst_abs:.3e} in {worst_var})")
sys.exit(0)
