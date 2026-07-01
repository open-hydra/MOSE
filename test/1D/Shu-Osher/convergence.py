"""
Shu-Osher problem (Mach-3 shock / sine-wave interaction) — Richardson
grid-convergence study.

There is no closed-form solution, so a grid-independent ("exact") density field
is estimated by **Richardson extrapolation** and every grid is measured against
that extrapolant.  Four MOSE solutions with refinement ratio r = 2 take part in
the study:

    N = 200, 400, 800   -> OUTPUT/field_x{1,2,4}.tec   (produced by ./MOSE.sh test)
    N = 1600            -> reference/reference.dat      (the finest grid)

All solutions are interpolated onto a common set of points in the active region
-4.5 <= x <= 3.5.  From the three finest grids the observed order of accuracy is

        p = ln( ||f_4h - f_2h||_1 / ||f_2h - f_h||_1 ) / ln(r)

and the Richardson extrapolant of the two finest grids is

        f_ext = f_h + (f_h - f_2h) / (r^p - 1).

The normalised L1 density error of each grid is then measured against f_ext.
The study passes if the observed order is sane (p > 0.5) and the error decreases
monotonically under refinement.

A two-panel, transparent, theme-aware SVG is written to OUTPUT/Shu-Osher.svg
(density profiles + the extrapolant on the left; log-log L1 error vs N with the
observed-order slope on the right) for the V&V documentation.

Usage:
    python convergence.py            # report, write SVG, exit 0/1
    python convergence.py --plot     # also show the figure

Exit 0 = pass, 1 = fail.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

R = 2.0             # grid refinement ratio
LO, HI = -4.5, 3.5  # active region the waves have reached
P_MIN = 0.5         # minimum acceptable observed order of accuracy
# Grids taking part in the study (coarse -> fine); the last is the reference file.
FIELD_GRIDS = [("field_x1.tec", 200), ("field_x2.tec", 400), ("field_x4.tec", 800)]
REF_N = 1600

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

parser = argparse.ArgumentParser(description="Shu-Osher Richardson convergence study")
parser.add_argument("--plot", action="store_true", help="show the figure")
args = parser.parse_args()

# --- load every grid (N, xc, rho), coarse -> fine ---------------------------
grids = []   # (N, xc, rho)
for fname, N in FIELD_GRIDS:
    path = os.path.join("OUTPUT", fname)
    if not os.path.exists(path):
        print(f"Shu-Osher: FAIL — missing {path} (run './MOSE.sh test' first)")
        sys.exit(1)
    x_, _, _, var_, _ = read_TEC(path)
    xc = 0.5 * (x_[0][:-1, 0, 0] + x_[0][1:, 0, 0])
    rho = var_[0][0][:, 0, 0]
    grids.append((N, xc, rho))

ref = np.loadtxt("reference/reference.dat")            # finest grid: N=1600
grids.append((REF_N, ref[:, 0], ref[:, 1]))

for (N, _, rho) in grids:
    if not np.all(np.isfinite(rho)):
        print(f"Shu-Osher: FAIL — non-finite density on N={N}")
        sys.exit(1)

# --- interpolate all grids onto a common fine set of points -----------------
# Use a uniform sampling of the active region at the finest resolution so the
# L1 sums are directly comparable across grids.
Xc = np.linspace(LO, HI, 4 * REF_N)
f = {N: np.interp(Xc, xc, rho) for (N, xc, rho) in grids}
Ns = np.array([g[0] for g in grids], dtype=float)      # [200, 400, 800, 1600]

# --- observed order from the three finest grids (h = 1/1600, 2h, 4h) --------
f_h, f_2h, f_4h = f[1600], f[800], f[400]
d_coarse = np.sum(np.abs(f_4h - f_2h))                 # ||f_4h - f_2h||_1
d_fine = np.sum(np.abs(f_2h - f_h))                    # ||f_2h - f_h||_1
p = np.log(d_coarse / d_fine) / np.log(R)

# --- Richardson extrapolant of the two finest grids -------------------------
f_ext = f_h + (f_h - f_2h) / (R ** p - 1.0)

# --- normalised L1 error of each grid vs the extrapolant --------------------
denom = np.sum(np.abs(f_ext))
errs = np.array([np.sum(np.abs(f[N] - f_ext)) / denom for N in Ns.astype(int)])

# Grid-Convergence Index on the finest solved grid (Roache, Fs = 1.25)
gci_fine = 1.25 * errs[-1] / (R ** p - 1.0)

print(f"Shu-Osher Richardson grid convergence (r={R:.0f})")
for N, e in zip(Ns.astype(int), errs):
    tag = " (reference grid)" if N == REF_N else ""
    print(f"  N={N:4d}   L1 error vs Richardson extrapolant : {e*100:6.3f} %{tag}")
print(f"  observed order p            : {p:.2f}")
print(f"  GCI (finest grid, Fs=1.25)  : {gci_fine*100:.3f} %")

monotone = np.all(np.diff(errs) < 0.0)
order_ok = np.isfinite(p) and p > P_MIN
status = "PASS" if (monotone and order_ok) else "FAIL"
print(f"  monotone error decrease : {monotone};  order > {P_MIN} : {order_ok}")
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

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))

# left: density profiles + Richardson extrapolant
colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
for (N, xc, rho), c in zip(grids, colors):
    axL.plot(xc, rho, "-", color=c, lw=1.0,
             label=f"MOSE (N={N})" + (" — ref" if N == REF_N else ""))
axL.plot(Xc, f_ext, "k--", lw=1.4, label="Richardson extrapolant")
axL.set_xlabel("x"); axL.set_ylabel(r"$\rho$")
axL.set_xlim(-5, 5); axL.grid(True, alpha=0.3); axL.legend()

# right: error vs the extrapolant
axR.loglog(Ns, errs, "ko-", lw=1.2, ms=6, label="L1 error vs extrapolant")
slope = errs[0] * (Ns[0] / Ns) ** p
axR.loglog(Ns, slope, "k--", lw=1.0, alpha=0.6, label=f"slope $p$ = {p:.2f}")
axR.set_xlabel("N (cells)"); axR.set_ylabel(r"normalised $L_1$ density error")
axR.grid(True, which="both", alpha=0.3); axR.legend()
axR.set_title(f"observed order $p$ ≈ {p:.2f}")

plt.tight_layout()
os.makedirs("OUTPUT", exist_ok=True)
plt.savefig(os.path.join("OUTPUT", "Shu-Osher.svg"),
            bbox_inches="tight", transparent=True)

if args.plot:
    plt.show()

sys.exit(0 if status == "PASS" else 1)
