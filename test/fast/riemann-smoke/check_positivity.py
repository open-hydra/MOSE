"""Positivity / finiteness check for the Riemann smoke-matrix test.

Reads OUTPUT/field.tec and asserts that density and pressure are strictly
positive and finite everywhere.  No reference solution is required — this is a
crash / robustness smoke test, not an accuracy test.

Exit code 0 = pass, 1 = fail.
"""
import sys
from pathlib import Path

import numpy as np

# Make the script work even when conda is not activated: locate ORION in-repo.
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
        "ORION package not found. Run 'conda activate base' or install ORION "
        "in the current Python environment."
    )

solver = sys.argv[1] if len(sys.argv) > 1 else "?"

[x_, y_, z_, var_, vnames] = read_TEC("OUTPUT/field.tec")
rho = var_[0][0][:, 0, 0]
p = var_[0][4][:, 0, 0]

bad = []
if not np.all(np.isfinite(rho)):
    bad.append("rho has non-finite values (NaN/Inf)")
if not np.all(np.isfinite(p)):
    bad.append("p has non-finite values (NaN/Inf)")
if np.nanmin(rho) <= 0.0:
    bad.append(f"min(rho)={np.nanmin(rho):.3e} <= 0")
if np.nanmin(p) <= 0.0:
    bad.append(f"min(p)={np.nanmin(p):.3e} <= 0")

if bad:
    print(f"  [{solver}] FAIL: " + "; ".join(bad))
    sys.exit(1)

print(f"  [{solver}] ok  (min rho={np.nanmin(rho):.4f}, min p={np.nanmin(p):.4f})")
sys.exit(0)
