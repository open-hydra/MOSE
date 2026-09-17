"""
2D viscous Method of Manufactured Solutions (N2) — order-of-accuracy verification.

A smooth, periodic field is made an EXACT steady solution of the compressible
Navier-Stokes equations by an analytic source term (injected by the MOSE_mms
driver, src/app/mms.f90).  Because the exact solution is known everywhere, the
discretisation error of each grid is measured directly for every primitive
variable and the *observed order of accuracy* is extracted by refinement:

        p = log2( e(N) / e(2N) )     (refinement ratio 2)

Unlike the inviscid isentropic vortex (N1), this exercises the full viscous
operator — velocity/temperature gradients and mesh metrics — so it catches
gradient/metric bugs the vortex cannot.  With mu = 10 (Re ~ 5) the viscous
fluxes are ~20 % of the convective ones, and the mean Mach number ~0.47 keeps
the compressible scheme out of the low-Mach accuracy-degradation regime.  A
second-order scheme gives p -> 2 in both the L2 and Linf norms of rho, u, v, p.

The study passes if, for every variable, the L2 error decreases monotonically
and the observed order on the finest available pair is at least ORDER_MIN.

A transparent, theme-aware two-panel SVG (L2 and Linf error vs N, all variables,
with an ideal second-order slope) is written to OUTPUT/mms.svg for the V&V page.

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
# refinement sequence (ratio 2). Only the grids whose field_<N>.tec is present
# are used, so the fast CTest subset (32/64) and the full documented study
# (32/64/128) share one verifier.
ALL_GRIDS = [32, 64, 128]

# manufactured-field parameters — must match build_ic.py and src/app/mms.f90
OM = 2.0 * np.pi
RHO0, ARHO = 1.0, 0.1
U0, AU = 40.0, 8.0
V0, AV = 30.0, 8.0
P0, AP = 8.0e3, 8.0e2

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

parser = argparse.ArgumentParser(description="Viscous MMS order study")
parser.add_argument("--plot", action="store_true", help="show the figure")
args = parser.parse_args()

VARS = ["rho", "u", "v", "p"]
LABELS = {"rho": r"$\rho$", "u": r"$u$", "v": r"$v$", "p": r"$p$"}


def exact(xc, yc):
    return {
        "rho": RHO0 + ARHO * np.sin(OM * xc) * np.cos(OM * yc),
        "u":   U0   + AU   * np.sin(OM * xc) * np.sin(OM * yc),
        "v":   V0   + AV   * np.cos(OM * xc) * np.cos(OM * yc),
        "p":   P0   + AP   * np.cos(OM * xc) * np.sin(OM * yc),
    }


GRIDS = [N for N in ALL_GRIDS if os.path.exists(f"reference/field_{N}.tec")]
if len(GRIDS) < 2:
    print("mms: FAIL — need at least two field_<N>.tec (run './MOSE.sh test' or "
          "'./MOSE.sh testfast' first)")
    sys.exit(1)

# normalised L2 and Linf error of each variable on each grid
L2 = {v: [] for v in VARS}
LI = {v: [] for v in VARS}
Ns = []
for N in GRIDS:
    x_, y_, _, var_, _ = read_TEC(f"reference/field_{N}.tec")
    xn, yn = x_[0][:, :, 0], y_[0][:, :, 0]
    xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
    yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])
    num = {"rho": var_[0][0][:, :, 0], "u": var_[0][1][:, :, 0],
           "v": var_[0][2][:, :, 0], "p": var_[0][4][:, :, 0]}
    ex = exact(xc, yc)
    if not np.all([np.all(np.isfinite(num[v])) for v in VARS]):
        print(f"mms: FAIL — non-finite field on N={N}")
        sys.exit(1)
    for v in VARS:
        d = num[v] - ex[v]
        L2[v].append(np.sqrt(np.mean(d ** 2)) / np.sqrt(np.mean(ex[v] ** 2)))
        LI[v].append(np.max(np.abs(d)) / np.max(np.abs(ex[v])))
    Ns.append(N)

Ns = np.array(Ns, dtype=float)
for v in VARS:
    L2[v] = np.array(L2[v]); LI[v] = np.array(LI[v])


def orders(e):
    return np.log2(e[:-1] / e[1:])


print("2D viscous MMS — order of accuracy (normalised errors)")
print(f"  grids: {', '.join(f'{int(N)}^2' for N in Ns)}")
for v in VARS:
    o2, oi = orders(L2[v]), orders(LI[v])
    print(f"  {v:3s}  L2   : " + "  ".join(f"{e:.3e}" for e in L2[v]) +
          "   orders " + " ".join(f"{o:.2f}" for o in o2))
    print(f"       Linf : " + "  ".join(f"{e:.3e}" for e in LI[v]) +
          "   orders " + " ".join(f"{o:.2f}" for o in oi))

# pass criterion: every variable monotone in L2 and finest-pair L2 order >= ORDER_MIN
monotone = all(np.all(np.diff(L2[v]) < 0.0) for v in VARS)
order_ok = all(orders(L2[v])[-1] >= ORDER_MIN for v in VARS)
status = "PASS" if (monotone and order_ok) else "FAIL"
print(f"  all vars: monotone L2 = {monotone};  finest-pair L2 order >= {ORDER_MIN} = {order_ok}")
print(f"  STATUS : {status}")

# --- figure: L2 and Linf error vs N, all variables, ideal 2nd-order slope -----
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

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for ax, norm, data in ((axes[0], "L_2", L2), (axes[1], "L_\\infty", LI)):
    for v in VARS:
        ax.loglog(Ns, data[v], "o-", lw=2.5, ms=9, label=LABELS[v])
    ref = data["rho"][0] * (Ns[0] / Ns) ** 2
    ax.loglog(Ns, ref, "k--", lw=2.0, alpha=0.6)

    ax.set_xlabel("N (cells per side)", fontsize=20)
    ax.set_ylabel(rf"normalised ${norm}$ error", fontsize=20)
    ax.grid(True, which="both", alpha=0.3)
    #ax.set_title(rf"${norm}$ norm", fontsize=22)
    ax.tick_params(labelsize=20)

# Single legend for the whole figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels,
           loc="upper center",
           ncol=4,
           fontsize=20,
           bbox_to_anchor=(0.5, 1.05))

plt.tight_layout(rect=[0, 0, 1, 0.93])  # leave room for the legend

os.makedirs("OUTPUT", exist_ok=True)
plt.savefig(os.path.join("OUTPUT", "mms.svg"),
            bbox_inches="tight",
            transparent=True)

if args.plot:
    plt.show()

sys.exit(0 if status == "PASS" else 1)
