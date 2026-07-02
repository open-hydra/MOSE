"""
2D isentropic Euler vortex (N1) — order-of-accuracy verification.

The isentropic vortex is a smooth (C-infinity) exact solution of the Euler
equations convected by a uniform mean flow.  After one full pass across the
periodic box it returns to its initial position, so the exact solution equals
the initial condition.  The density error of each grid is measured against that
exact field, and the *observed order of accuracy* is extracted by refinement:

        p = log2( e(N) / e(2N) )     (refinement ratio 2)

A smooth solution + MUSCL reconstruction should give p -> 2 (this is why the
kinked Gresho vortex cannot be used for N1).

The study passes if the error decreases monotonically and the observed order on
the finest pair is at least ORDER_MIN.

A transparent, theme-aware SVG (log-log error vs N with the fitted slope) is
written to OUTPUT/isentropic-vortex.svg for the V&V documentation.

Usage:
    python verify.py            # report, write SVG, exit 0/1
    python verify.py --plot     # also show the figure

Exit 0 = pass, 1 = fail.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

ORDER_MIN = 1.6                 # minimum acceptable observed order (finest pair)
GRIDS = [32, 64, 128]           # refinement sequence (ratio 2)
# vortex parameters — must match build_ic.py
GAMMA, BETA, X0, Y0 = 1.4, 5.0, 5.0, 5.0

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

parser = argparse.ArgumentParser(description="Isentropic-vortex order study")
parser.add_argument("--plot", action="store_true", help="show the figure")
args = parser.parse_args()


def exact_rho(xc, yc):
    r2 = (xc - X0) ** 2 + (yc - Y0) ** 2
    dT = -((GAMMA - 1.0) * BETA * BETA) / (8.0 * GAMMA * np.pi ** 2) * np.exp(1.0 - r2)
    return (1.0 + dT) ** (1.0 / (GAMMA - 1.0))


Ns, errs = [], []
for N in GRIDS:
    path = f"reference/field_{N}.tec"
    if not os.path.exists(path):
        print(f"isentropic-vortex: FAIL — missing {path} (run './MOSE.sh test' first)")
        sys.exit(1)
    x_, y_, _, var_, _ = read_TEC(path)
    xn = x_[0][:, :, 0]
    yn = y_[0][:, :, 0]
    xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
    yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])
    rho = var_[0][0][:, :, 0]
    if not np.all(np.isfinite(rho)):
        print(f"isentropic-vortex: FAIL — non-finite density on N={N}")
        sys.exit(1)
    re = exact_rho(xc, yc)
    e = np.sqrt(np.mean((rho - re) ** 2)) / np.sqrt(np.mean(re ** 2))   # normalised L2
    Ns.append(N)
    errs.append(e)

Ns = np.array(Ns, dtype=float)
errs = np.array(errs)
orders = np.log2(errs[:-1] / errs[1:])          # between successive grids

print("2D isentropic vortex — order of accuracy (normalised L2 density error)")
for N, e in zip(Ns.astype(int), errs):
    print(f"  N={N:4d}^2   error : {e*100:7.4f} %")
for i, p in enumerate(orders):
    print(f"  observed order N={int(Ns[i])}->{int(Ns[i+1])} : {p:.2f}")

monotone = np.all(np.diff(errs) < 0.0)
order_ok = orders[-1] >= ORDER_MIN
status = "PASS" if (monotone and order_ok) else "FAIL"
print(f"  finest-pair order >= {ORDER_MIN} : {order_ok};  monotone : {monotone}")
print(f"  STATUS : {status}")

# --- figure -----------------------------------------------------------------
import matplotlib as mpl
if not args.plot:
    mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

plt.figure(figsize=(7, 6))
plt.loglog(Ns, errs, "bo-", lw=2.0, ms=9, label="MOSE (normalised $L_2$)")
ref2 = errs[0] * (Ns[0] / Ns) ** 2
plt.loglog(Ns, ref2, "k--", lw=1.5, alpha=0.6, label="2nd-order slope")
plt.xlabel("N (cells per side)", fontsize=20)
plt.ylabel(r"normalised $L_2$ density error", fontsize=20)
plt.grid(True, which="both", alpha=0.3)
plt.legend(fontsize=20)
plt.tick_params(labelsize=20)
plt.tight_layout()
os.makedirs("OUTPUT", exist_ok=True)
plt.savefig(os.path.join("OUTPUT", "isentropic-vortex.svg"),
            bbox_inches="tight", transparent=True)

if args.plot:
    plt.show()

sys.exit(0 if status == "PASS" else 1)
