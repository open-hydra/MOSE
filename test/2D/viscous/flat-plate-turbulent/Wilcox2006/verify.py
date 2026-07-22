"""
Turbulent flat-plate skin-friction validation — Wilcox2006.

Reads the wall shear stress from OUTPUT/wall.tec, forms the skin-friction
coefficient Cf(x), and compares it against the reference curve in
reference/cf.dat (NASA turbulence-model verification data).

Usage:
    python verify.py            # report PASS/FAIL, save OUTPUT/cf.svg
    python verify.py --plot     # also show the figure interactively

Exit code 0 = pass (RMS relative error below tolerance), 1 = fail.
"""

# ---- model-specific configuration -------------------------------------------
MODEL        = "Wilcox2006"
REF_MATCH    = ["CFL3D"]     # substrings (all must be in the zone label)
TOLERANCE    = 8.0          # percent RMS relative error
MU           = 1.1858685985e-5 # laminar viscosity [Pa.s] (mil)
X_LO         = 0.02            # lower x bound for the error metric [m]
LABEL_FS     = 15
# -----------------------------------------------------------------------------

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

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

parser = argparse.ArgumentParser(description=f"{MODEL} flat-plate validation")
parser.add_argument("--plot", action="store_true", help="show interactive figure")
args = parser.parse_args()


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.tick_params(labelsize=LABEL_FS - 2)
    ax.legend(loc="best", fontsize=LABEL_FS - 2)
    ax.grid(True, alpha=0.3)


def read_cf_dat(path):
    """Return {zone_label: (x, cf)} from a (multi-zone) Tecplot ASCII file."""
    zones, label, xs, cfs = {}, None, [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.upper().startswith("VARIABLES"):
                continue
            if s.upper().startswith("ZONE") or s.upper().startswith("ZONE,"):
                if label is not None and xs:
                    zones[label] = (np.array(xs), np.array(cfs))
                m = re.search(r'T\s*=\s*"([^"]+)"', s, re.IGNORECASE)
                label, xs, cfs = (m.group(1) if m else s), [], []
            else:
                p = s.split()
                if len(p) >= 2:
                    try:
                        xs.append(float(p[0])); cfs.append(float(p[1]))
                    except ValueError:
                        pass
    if label is not None and xs:
        zones[label] = (np.array(xs), np.array(cfs))
    return zones


def pick_reference(zones, must_contain):
    for label, data in zones.items():
        low = label.lower()
        if all(tok.lower() in low for tok in must_contain):
            return label, data
    raise SystemExit(
        f"verify.py: no reference zone matching {must_contain} in cf.dat "
        f"(have: {list(zones)})"
    )


def read_wall_cf(path, rho_inf, U_inf):
    """Return (x_cell, Cf) on the plate (x > 0). wall.tec var 1 is tauX."""
    xw, _, _, vw, _ = read_TEC(path)
    xn = xw[0][:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    tauX = np.abs(vw[0][1][:, 0, 0])
    m = xc > 0.0
    return xc[m], tauX[m] / (0.5 * rho_inf * U_inf ** 2)


# ---- freestream state from the field --------------------------------------
_, _, _, fv, _ = read_TEC("OUTPUT/field.tec")
U_inf   = fv[0][1][0, -1, 0]
rho_inf = fv[0][0][0, -1, 0]
nu      = MU / rho_inf
print("Freestream:  U_inf=%.4f m/s  rho_inf=%.4f kg/m^3  Re/L=%.1f /m"
      % (U_inf, rho_inf, U_inf / nu))

# ---- MOSE skin friction ----------------------------------------------------
x_m, Cf_m = read_wall_cf("OUTPUT/wall.tec", rho_inf, U_inf)

# ---- reference & error metric ---------------------------------------------
zones = read_cf_dat(Path("reference/cf.dat"))
ref_label, (x_r, cf_r) = pick_reference(zones, REF_MATCH)

x_hi = min(x_m.max(), x_r.max())
mask = (x_m >= X_LO) & (x_m <= x_hi)
cf_ref_i = interp1d(x_r, cf_r, kind="linear")(x_m[mask])
err = np.abs(Cf_m[mask] - cf_ref_i) / cf_ref_i * 100.0
rms, emax = float(np.sqrt(np.mean(err ** 2))), float(err.max())

print(f"\n{MODEL} skin friction vs '{ref_label}'")
print(f"  RMS relative error : {rms:.2f} %")
print(f"  Max relative error : {emax:.2f} %")
status = "PASS" if rms < TOLERANCE else "FAIL"
print(f"  STATUS : {status}  (tolerance {TOLERANCE:.1f} % RMS)")

# ---- figure (optional; matplotlib only needed here) -----------------------
try:
    import matplotlib
    if not args.plot:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_m, Cf_m, "bd-", ms=12, markevery=5, label=f"MOSE", color="#0072B2", linewidth=4.0)
    ax.plot(x_r[x_r > 0], cf_r[x_r > 0], label=ref_label, color="#E69F00", linestyle="--", linewidth=4.0)
    ax.set_xlim(0.0, 1.8); ax.set_ylim(0.002, 0.006)
    style(ax, r"$x$  [m]", r"$C_f$")
    ax.set_xlim(0.0, 1.8); ax.set_ylim(0.002, 0.006)
    fig.tight_layout()
    os.makedirs("OUTPUT", exist_ok=True)
    if args.plot:
        plt.show()
    else:
        fig.savefig("OUTPUT/cf.svg")
        print("  figure : OUTPUT/cf.svg")
except ImportError:
    print("  (matplotlib not available — skipping figure)")

sys.exit(0 if status == "PASS" else 1)
