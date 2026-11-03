"""
Turbulent flat-plate validation — SSG/LRR-omega Reynolds-stress model (SSGLRR).

Reads OUTPUT/wall.tec and OUTPUT/field.tec and compares with the NASA Turbulence
Modeling Resource (TMR) results for SSGLRR-RSM-w2012,
https://tmbwg.github.io/turbmodels/flatplate_ssglrrrsm.html :

  1. skin friction Cf(x) against reference/cf.dat (CFL3D and TAU, finest
     545 x 385 grid);
  2. velocity profile u+(y+) at x = 0.97 against reference/uplus.dat (CFL3D,
     finest grid);
  3. the Reynolds stresses in the log layer of that profile: a Reynolds-stress
     model must return u'u' > w'w' > v'v' and -u'v'/k close to 0.3.

reference/cf-convergence.dat is the TMR grid-convergence table of Cf at
x = 0.97: its 3264-cell row is the grid of this test (69 x 49 nodes).

Usage:
    python verify.py            # report PASS/FAIL, save OUTPUT/cf.svg, OUTPUT/uplus.svg
    python verify.py --plot     # also show the figures interactively

Exit code 0 = pass, 1 = fail.
"""

# ---- configuration ----------------------------------------------------------
MODEL        = "SSGLRR"
REF_MATCH    = ["CFL3D"]     # substrings (all must be in the zone label)
TOLERANCE    = 8.0           # percent RMS relative error on Cf(x)
TOL_PROFILE  = 5.0           # percent RMS relative error on u+(y+)
MU           = 1.1858685985e-5 # laminar viscosity [Pa.s] (mil)
X_LO         = 0.02            # lower x bound for the error metric [m]
X_PROFILE    = 0.97008         # station of the velocity profile [m]
YPLUS_RANGE  = (2.0, 1000.0)   # y+ range of the profile error metric
YPLUS_LOG    = (50.0, 300.0)   # log layer, for the Reynolds-stress check
CF_GRID      = {"CFL3D": 0.0026541927, "TAU": 0.0026305596}   # TMR, Cf(x=0.97), 3264 cells
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


def read_zones(path):
    """Return {zone_label: (col1, col2)} from a (multi-zone) Tecplot ASCII file."""
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
        f"verify.py: no reference zone matching {must_contain} "
        f"(have: {list(zones)})"
    )


# ---- MOSE solution ---------------------------------------------------------
xn, yn, _, fv, names = read_TEC("OUTPUT/field.tec")
var = {n: fv[0][i] for i, n in enumerate(names[3:])}
rho, u = var["rho(1)"][:, :, 0], var["u"][:, :, 0]
U_inf, rho_inf = u[0, -1], rho[0, -1]
print("Freestream:  U_inf=%.4f m/s  rho_inf=%.4f kg/m^3  Re/L=%.1f /m"
      % (U_inf, rho_inf, U_inf * rho_inf / MU))

xw, _, _, vw, wnames = read_TEC("OUTPUT/wall.tec")
wall = {n: vw[0][i][:, 0, 0] for i, n in enumerate(wnames[3:])}
xnod = xw[0][:, 0, 0]
xc = 0.5 * (xnod[:-1] + xnod[1:])
tauw = np.abs(wall["tauX"])
plate = xc > 0.0
x_m, Cf_m = xc[plate], tauw[plate] / (0.5 * rho_inf * U_inf ** 2)

# ---- 1. skin friction ------------------------------------------------------
zones = read_zones(Path("reference/cf.dat"))
ref_label, (x_r, cf_r) = pick_reference(zones, REF_MATCH)
_, (x_t, cf_t) = pick_reference(zones, ["TAU"])

x_hi = min(x_m.max(), x_r.max())
mask = (x_m >= X_LO) & (x_m <= x_hi)
cf_ref_i = interp1d(x_r, cf_r, kind="linear")(x_m[mask])
err = np.abs(Cf_m[mask] - cf_ref_i) / cf_ref_i * 100.0
rms, emax = float(np.sqrt(np.mean(err ** 2))), float(err.max())
cf_097 = float(interp1d(x_m, Cf_m)(X_PROFILE))

print(f"\n{MODEL} skin friction vs '{ref_label}' (finest grid)")
print(f"  RMS relative error : {rms:.2f} %")
print(f"  Max relative error : {emax:.2f} %")
print(f"  Cf at x = 0.97     : {cf_097:.6f}   (TMR on this grid: CFL3D "
      f"{CF_GRID['CFL3D']:.6f}, TAU {CF_GRID['TAU']:.6f}; "
      f"finest grid: {float(interp1d(x_r, cf_r)(X_PROFILE)):.6f})")
ok_cf = rms < TOLERANCE
print(f"  STATUS : {'PASS' if ok_cf else 'FAIL'}  (tolerance {TOLERANCE:.1f} % RMS)")

# ---- 2. velocity profile at x = 0.97 --------------------------------------
i0 = int(np.argmin(np.abs(xc - X_PROFILE)))
ynod = yn[0][i0, :, 0]
yc = 0.5 * (ynod[:-1] + ynod[1:])
nu_w = MU / rho[i0, 0]
u_tau = np.sqrt(tauw[i0] / rho[i0, 0])
yplus, uplus = yc * u_tau / nu_w, u[i0, :] / u_tau

pz = read_zones(Path("reference/uplus.dat"))
_, (logy_r, up_r) = pick_reference(pz, ["0.97008"])
sel = (yplus >= YPLUS_RANGE[0]) & (yplus <= YPLUS_RANGE[1])
up_ref_i = interp1d(logy_r, up_r, kind="linear")(np.log10(yplus[sel]))
errp = np.abs(uplus[sel] - up_ref_i) / up_ref_i * 100.0
rmsp = float(np.sqrt(np.mean(errp ** 2)))

print(f"\n{MODEL} velocity profile u+(y+) at x = {xc[i0]:.4f} vs CFL3D (finest grid)")
print(f"  first cell y+      : {yplus[0]:.2f}")
print(f"  RMS relative error : {rmsp:.2f} %  ({YPLUS_RANGE[0]:g} < y+ < {YPLUS_RANGE[1]:g}, {int(sel.sum())} cells)")
ok_up = rmsp < TOL_PROFILE
print(f"  STATUS : {'PASS' if ok_up else 'FAIL'}  (tolerance {TOL_PROFILE:.1f} % RMS)")

# ---- 3. Reynolds stresses in the log layer --------------------------------
R = {c: var[f"r{c}"][i0, :, 0] / rho[i0, :] for c in ("u'u'", "v'v'", "w'w'", "u'v'")}
k = 0.5 * (R["u'u'"] + R["v'v'"] + R["w'w'"])
log = (yplus >= YPLUS_LOG[0]) & (yplus <= YPLUS_LOG[1])
a = {c: float(np.mean(R[c][log] / k[log])) for c in R}
ok_rs = bool(np.all(R["u'u'"][log] > R["w'w'"][log]) and np.all(R["w'w'"][log] > R["v'v'"][log])
             and 0.2 < -a["u'v'"] < 0.4)

print(f"\n{MODEL} Reynolds stresses in the log layer ({YPLUS_LOG[0]:g} < y+ < {YPLUS_LOG[1]:g}, {int(log.sum())} cells)")
print("  u'u'/k = %.3f   v'v'/k = %.3f   w'w'/k = %.3f   -u'v'/k = %.3f"
      % (a["u'u'"], a["v'v'"], a["w'w'"], -a["u'v'"]))
print(f"  STATUS : {'PASS' if ok_rs else 'FAIL'}  (u'u' > w'w' > v'v' and 0.2 < -u'v'/k < 0.4)")

status = "PASS" if (ok_cf and ok_up and ok_rs) else "FAIL"
print(f"\n  OVERALL : {status}")

# ---- figures (optional; matplotlib only needed here) ----------------------
try:
    import matplotlib
    if not args.plot:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs("OUTPUT", exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_m, Cf_m, "bd-", ms=12, markevery=5, label="MOSE", color="#0072B2", linewidth=4.0)
    ax.plot(x_r[x_r > 0], cf_r[x_r > 0], label=ref_label, color="#E69F00", linestyle="--", linewidth=4.0)
    ax.plot(x_t[x_t > 0], cf_t[x_t > 0], label="TAU", color="#009E73", linestyle=":", linewidth=4.0)
    ax.set_xlim(0.0, 1.8); ax.set_ylim(0.002, 0.006)
    style(ax, r"$x$  [m]", r"$C_f$")
    fig.tight_layout()

    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    ax1.semilogx(yplus, uplus, "bd", ms=10, label="MOSE", color="#0072B2")
    ax1.semilogx(10.0 ** logy_r, up_r, label="CFL3D", color="#E69F00", linestyle="--", linewidth=4.0)
    ax1.set_xlim(0.5, 2.0e4); ax1.set_ylim(0.0, 30.0)
    style(ax1, r"$y^+$", r"$u^+$")
    for c, col, lab in (("u'u'", "#0072B2", r"$\overline{u'u'}^+$"), ("v'v'", "#D55E00", r"$\overline{v'v'}^+$"),
                        ("w'w'", "#009E73", r"$\overline{w'w'}^+$")):
        ax2.semilogx(yplus, R[c] / u_tau ** 2, "o-", ms=6, color=col, label=lab, linewidth=2.5)
    ax2.semilogx(yplus, -R["u'v'"] / u_tau ** 2, "o-", ms=6, color="#CC79A7", label=r"$-\overline{u'v'}^+$", linewidth=2.5)
    ax2.set_xlim(0.5, 2.0e4); ax2.set_ylim(0.0, 4.0)
    style(ax2, r"$y^+$", r"Reynolds stresses / $u_\tau^2$")
    fig2.tight_layout()

    if args.plot:
        plt.show()
    else:
        fig.savefig("OUTPUT/cf.svg")
        fig2.savefig("OUTPUT/uplus.svg")
        print("  figures : OUTPUT/cf.svg  OUTPUT/uplus.svg")
except ImportError:
    print("  (matplotlib not available — skipping figures)")

sys.exit(0 if status == "PASS" else 1)
