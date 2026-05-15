"""
Flat-plate turbulent boundary layer — Spalart-Allmaras validation

Reads wall shear stress from OUTPUT/wall.tec and compares the skin
friction coefficient against CFL3D, FUN3D, and TAU reference data
from reference/cf.dat (NASA turbulence model validation cases).

Usage:
    python verify.py           # saves figure to OUTPUT/verify_Cf.png
    python verify.py --plot    # shows interactive figure
"""

import argparse
import sys
from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

output_dir = "OUTPUT"

# -------------------------------------------------------------
# Locate ORION Python package
# -------------------------------------------------------------
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
        "ORION package not found. Run 'conda activate base' or "
        "install ORION in the current Python environment."
    )

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Turbulent flat-plate validation')
parser.add_argument('--plot', action='store_true', help='Show interactive plots')
args = parser.parse_args()

x_, y_, z_, var_, _ = read_TEC("OUTPUT/field.tec")

gamma  = 1.40
R_gas  = 287.0
mu     = 1.1858685985e-5          # mil from input.ini

U_inf   = var_[0][1][0,-1,0]
p_inf   = var_[0][4][0,-1,0]
rho_inf = var_[0][0][0,-1,0]
nu      = mu / rho_inf

print("Freestream conditions")
print(f"  U_inf   = {U_inf:.4f} m/s")
print(f"  rho_inf = {rho_inf:.4f} kg/m³")
print(f"  mu      = {mu:.2e} Pa·s")
print(f"  nu      = {nu:.4e} m²/s")
print(f"  Re/L    = {U_inf/nu:.1f} m⁻¹")

# -------------------------------------------------------------
# Read reference data from reference/cf.dat
# Format: multi-zone Tecplot ASCII  VARIABLES = "x","cf"
#         zones named "<solver>, grid level <n>"  (n=1 finest)
# -------------------------------------------------------------
def _read_cf_dat(path):
    """Return dict  {zone_label: (x_array, cf_array)}."""
    zones = {}
    current_label = None
    xs, cfs = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.upper().startswith('VARIABLES'):
                continue
            if line.upper().startswith('ZONE'):
                if current_label is not None and xs:
                    zones[current_label] = (np.array(xs), np.array(cfs))
                # parse label from  T="..."
                import re
                m = re.search(r'T\s*=\s*"([^"]+)"', line, re.IGNORECASE)
                current_label = m.group(1) if m else line
                xs, cfs = [], []
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        xs.append(float(parts[0]))
                        cfs.append(float(parts[1]))
                    except ValueError:
                        pass
    if current_label is not None and xs:
        zones[current_label] = (np.array(xs), np.array(cfs))
    return zones

ref_all = _read_cf_dat(Path("reference/cf.dat"))
ref_sst   = _read_cf_dat(Path("reference/cf-sst.dat"))
ref_wilcox = _read_cf_dat(Path("reference/cf-wilcox.dat"))

# Keep only finest grid (level 1) of each solver for the main comparison
ref_solvers = {}
for label, data in ref_all.items():
    if 'level 1' in label.lower():
        solver = label.split(',')[0].strip()
        if solver != 'TAU':
            ref_solvers[solver] = (label, data)

print("\nReference zones loaded:")
for label, (x, cf) in ref_all.items():
    print(f"  {label:35s}  ({len(x)} pts)")
for label, (x, cf) in ref_sst.items():
    print(f"  [SST]   {label:35s}  ({len(x)} pts)")
for label, (x, cf) in ref_wilcox.items():
    print(f"  [Wilcox]{label:35s}  ({len(x)} pts)")

# -------------------------------------------------------------
# Helper: read a wall.tec-style file and return (x_cell, Cf)
# Cell-centered variables (indices 0-6):
#   0:y+  1:tauX  2:tauY  3:tauZ  4:pw  5:Tw  6:qw
# tauX is negative (wall drags fluid in -x), so take abs()
# -------------------------------------------------------------
def _read_wall_cf(path, rho_inf, U_inf):
    """Return (x_cell, Cf) arrays for x > 0 (plate region)."""
    xw_, yw_, zw_, varw_, _ = read_TEC(path)
    xw_n = xw_[0][:, 0, 0]
    xw_c = 0.5 * (xw_n[:-1] + xw_n[1:])
    tauX = np.abs(varw_[0][1][:, 0, 0])
    mask = xw_c > 0.0
    x_c  = xw_c[mask]
    Cf   = tauX[mask] / (0.5 * rho_inf * U_inf**2)
    return x_c, Cf

# Fine SA solution (used for reference comparison)
x_plt, Cf_mose = _read_wall_cf("OUTPUT/wall-sa-fine.tec", rho_inf, U_inf)

# -------------------------------------------------------------
# Error metric vs CFL3D level 1 (interpolated to MOSE x points)
# -------------------------------------------------------------
cfl3d_label, (x_ref, cf_ref) = ref_solvers['CFL3D']
# restrict interpolation to overlapping x range
x_lo = 0.0125
x_hi = min(x_plt.max(), x_ref.max())
mask_err = (x_plt >= x_lo) & (x_plt <= x_hi)

cf_ref_interp = interp1d(x_ref, cf_ref, kind='linear')(x_plt[mask_err])
err_rel = np.abs(Cf_mose[mask_err] - cf_ref_interp) / cf_ref_interp * 100.0
err_rms = np.sqrt(np.mean(err_rel**2))
err_max = err_rel.max()

print(f"\nSkin friction vs {cfl3d_label}")
print(f"  RMS relative error : {err_rms:.2f} %")
print(f"  Max relative error : {err_max:.2f} %")
status = "PASS" if err_rms < 5.0 else "FAIL"
print(f"  STATUS : {status}  (tolerance 5 % RMS)")

# -------------------------------------------------------------
# Plots
# -------------------------------------------------------------
solver_styles = {
    'CFL3D': dict(color='black', ls='-', lw=3.0),
    'FUN3D': dict(color='grey', ls='--', lw=3.0),
}

# --- Cf vs x ---
# Figure 1: Skin friction coefficient
fig1 = plt.figure(figsize=(10, 6))
ax_cf = fig1.add_subplot(111)
ax_cf.plot(x_plt, Cf_mose, 'go', ms=6, lw=2.0, label='MOSE')
for solver, (label, (x_r, cf_r)) in ref_solvers.items():
    style = solver_styles.get(solver, dict(color='gray', ls='--', lw=2.0))
    ax_cf.plot(x_r, cf_r, label=label, **style)
# Sample the 
ax_cf.set_xlim(0.0, 1.8)
ax_cf.set_ylim(0.002, 0.006)
ax_cf.set_xlabel(r'$x$  [m]')
ax_cf.set_ylabel(r'$C_f$')
ax_cf.legend(loc='best')
ax_cf.grid(True, alpha=0.3)

plt.tight_layout()

# -------------------------------------------------------------
# Figure 2: coarse grid comparison (SA-coarse, SST, Wilcox)
# -------------------------------------------------------------
coarse_cases = [
    ("OUTPUT/wall-sa-coarse.tec",  "MOSE SA ",   dict(color='green',  ls='-',  lw=2.0, marker='o', ms=4)),
    ("OUTPUT/wall-sst.tec",        "MOSE SST",   dict(color='blue',   ls='--', lw=2.0, marker='s', ms=4)),
    ("OUTPUT/wall-wilcox.tec",     "MOSE Wilcox",dict(color='orange', ls=':',  lw=2.0, marker='^', ms=4)),
]

fig2 = plt.figure(figsize=(10, 6))
ax_co = fig2.add_subplot(111)
for path, label, style in coarse_cases:
    x_c, cf_c = _read_wall_cf(path, rho_inf, U_inf)
    ax_co.plot(x_c, cf_c, label=label, **style)
ax_co.set_xlim(0.0, 1.8)
ax_co.set_ylim(0.002, 0.006)
ax_co.set_xlabel(r'$x$  [m]')
ax_co.set_ylabel(r'$C_f$')
ax_co.legend(loc='best')
ax_co.grid(True, alpha=0.3)
plt.tight_layout()

# -------------------------------------------------------------
# Figure 3: MOSE SST vs SST reference (FUN3D triangles)
# -------------------------------------------------------------
x_sst, Cf_sst = _read_wall_cf("OUTPUT/wall-sst.tec", rho_inf, U_inf)

solver_styles = {
    'FUN3D': dict(color='grey', ls='--', lw=3.0),
}

fig3 = plt.figure(figsize=(10, 6))
ax_sst = fig3.add_subplot(111)
ax_sst.plot(x_sst, Cf_sst, 'bs', ms=5, lw=2.0, label='MOSE')
for label, (x_r, cf_r) in ref_sst.items():
    solver = label.split(',')[0].strip()
    style = solver_styles.get(solver, dict(color='gray', ls='--', lw=3.0))
    ax_sst.plot(x_r, cf_r, label=label, **style)
ax_sst.set_xlim(0.0, 1.8)
ax_sst.set_ylim(0.002, 0.006)
ax_sst.set_xlabel(r'$x$  [m]')
ax_sst.set_ylabel(r'$C_f$')
ax_sst.legend(loc='best')
ax_sst.grid(True, alpha=0.3)
plt.tight_layout()

# -------------------------------------------------------------
# Figure 4: MOSE Wilcox vs Wilcox reference (CFL3D / FUN3D)
# -------------------------------------------------------------
x_wlx, Cf_wlx = _read_wall_cf("OUTPUT/wall-wilcox.tec", rho_inf, U_inf)

solver_styles = {
    'CFL3D': dict(color='black', ls='-', lw=3.0),
    'FUN3D': dict(color='grey', ls='--', lw=3.0),
}

fig4 = plt.figure(figsize=(10, 6))
ax_wlx = fig4.add_subplot(111)
ax_wlx.plot(x_wlx, Cf_wlx, 'o', color='red', ms=5, lw=2.0, label='MOSE')
for label, (x_r, cf_r) in ref_wilcox.items():
    solver = label.split(',')[0].strip()
    style = solver_styles.get(solver, dict(color='gray', ls='--', lw=3.0))
    ax_wlx.plot(x_r[x_r > 0.0], cf_r[x_r > 0.0], label=label, **style)
ax_wlx.set_xlim(0.0, 1.8)
ax_wlx.set_ylim(0.002, 0.006)
ax_wlx.set_xlabel(r'$x$  [m]')
ax_wlx.set_ylabel(r'$C_f$')
ax_wlx.legend(loc='best')
ax_wlx.grid(True, alpha=0.3)
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_fig1 = Path(os.path.join(output_dir, "TP-cf-sa.svg"))
    out_fig3 = Path(os.path.join(output_dir, "TP-cf-sst.svg"))
    out_fig4 = Path(os.path.join(output_dir, "TP-cf-wilcox.svg"))
    fig1.savefig(out_fig1, dpi=150)
    fig3.savefig(out_fig3, dpi=150)
    fig4.savefig(out_fig4, dpi=150)
    print(f"\nFigures saved to {out_fig1}, {out_fig3}, {out_fig4}")
