"""
Shu-Osher problem (Mach-3 shock / sine-wave interaction) — convergence check.

There is no closed-form solution, so the coarse run (INPUT/ic.tec, N=200) is
compared against a well-resolved MOSE reference (reference/reference.dat,
N=1600), both with the same scheme.  The metric is the normalised L1 error of
density over the active region; it passes if the coarse solution stays close to
the resolved one (i.e. the scheme converges and did not blow up).

Usage:
    python verify.py            # report error, exit 0/1
    python verify.py --plot     # also show rho(x)

Exit 0 = pass, 1 = fail.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

TOL = 5.0e-2   # normalised L1 density error vs the N=1600 reference

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

parser = argparse.ArgumentParser(description="Shu-Osher convergence check")
parser.add_argument("--plot", action="store_true", help="show rho(x)")
args = parser.parse_args()

# --- coarse MOSE solution ---------------------------------------------------
x_, _, _, var_, _ = read_TEC("OUTPUT/field.tec")
xc = 0.5 * (x_[0][:-1, 0, 0] + x_[0][1:, 0, 0])
rho = var_[0][0][:, 0, 0]

if not np.all(np.isfinite(rho)):
    print("Shu-Osher: FAIL — non-finite density")
    sys.exit(1)

# --- resolved reference -----------------------------------------------------
ref = np.loadtxt("reference/reference.dat")
x_ref, rho_ref = ref[:, 0], ref[:, 1]

# Compare over the region the waves have reached (exclude the inflow plateau and
# the far-field ahead of the shock where both are trivially identical).
lo, hi = -4.5, 3.5
m = (xc >= lo) & (xc <= hi)
rho_ref_i = np.interp(xc[m], x_ref, rho_ref)
l1 = np.sum(np.abs(rho[m] - rho_ref_i)) / np.sum(np.abs(rho_ref_i))

print(f"Shu-Osher (N={len(xc)} vs reference N={len(x_ref)})")
print(f"  normalised L1 density error : {l1*100:.2f} %")
status = "PASS" if l1 < TOL else "FAIL"
print(f"  STATUS : {status}  (tolerance {TOL*100:.1f} %)")

if args.plot:
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 5))
    plt.plot(x_ref, rho_ref, "k-", lw=1.5, label=f"reference N={len(x_ref)}")
    plt.plot(xc, rho, "ro", ms=3, label=f"MOSE N={len(xc)}")
    plt.xlabel("x"); plt.ylabel("rho"); plt.legend(); plt.grid(alpha=0.3)
    plt.tight_layout(); plt.show()

sys.exit(0 if status == "PASS" else 1)
