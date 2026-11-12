"""Mass-conservation check for the closed-box test (I4).

A shock-tube initial state is evolved inside a box closed by slip walls
(BC type 300) on every face, so there is no flux through the boundary and the
total mass must stay constant to round-off, independent of the internal wave
dynamics.

Compares total mass from the initial condition (INPUT/ic.tec) against the final
solution (OUTPUT/field.tec).  The mesh is uniform, so cell volume is constant
and cancels in the relative change — we sum density directly.

Exit code 0 = pass, 1 = fail.
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
        "ORION package not found. Run 'conda activate base' or install ORION "
        "in the current Python environment."
    )

TOL = 1.0e-10  # relative mass-conservation tolerance

# Density is the first cell-centred variable in both files.
[_, _, _, ivar, _] = read_TEC("INPUT/ic.tec")
[_, _, _, fvar, _] = read_TEC("OUTPUT/field.tec")
rho0 = ivar[0][0][:, 0, 0]
rho1 = fvar[0][0][:, 0, 0]

m0 = float(np.sum(rho0))
m1 = float(np.sum(rho1))
rel = abs(m1 - m0) / abs(m0)

if not np.all(np.isfinite(rho1)):
    print("Conservation: FAIL — non-finite density in final field")
    sys.exit(1)

if rel > TOL:
    print(f"Conservation: FAIL — total mass changed by {rel:.3e} "
          f"(initial={m0:.10e}, final={m1:.10e}, tol={TOL:.0e})")
    sys.exit(1)

print(f"Conservation: ok — total mass conserved to {rel:.3e} "
      f"(initial={m0:.6e}, final={m1:.6e})")
sys.exit(0)
