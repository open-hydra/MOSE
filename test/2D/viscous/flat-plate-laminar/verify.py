"""
Flat-plate laminar boundary layer — Blasius validation

Reads wall shear stress from OUTPUT/wall.tec and compares the skin
friction coefficient against the Blasius similarity solution computed
on-the-fly by solving the Blasius ODE  f''' + f·f'' = 0.

Also compares the velocity profile u/U_inf vs η at a chosen x station
against the Blasius self-similar profile, and plots the flow field solution.

Generates three figures:
  1. Skin friction coefficient vs Blasius
  2. Velocity profile comparison at x_ref
  3. Flow field velocity magnitude contours

Usage:
    python verify.py           # saves figures to OUTPUT/verify_*.png
    python verify.py --plot    # shows interactive plots
"""

import argparse
import sys
from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
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
parser = argparse.ArgumentParser(description='Flat-plate Blasius validation')
parser.add_argument('--plot', action='store_true', help='Show interactive plots')
args = parser.parse_args()

x_, y_, z_, var_, _ = read_TEC("OUTPUT/field.tec")

gamma  = 1.40
R_gas  = 287.0
mu     = 1e-3
T_inf  = var_[0][5][0,0,0]

a_inf   = np.sqrt(gamma * R_gas * T_inf)
U_inf   = var_[0][1][0,0,0]
p_inf   = var_[0][4][0,0,0]
rho_inf = p_inf / (R_gas * T_inf)
nu      = mu / rho_inf

print("Freestream conditions")
print(f"  U_inf   = {U_inf:.4f} m/s")
print(f"  rho_inf = {rho_inf:.4f} kg/m³")
print(f"  mu      = {mu:.2e} Pa·s")
print(f"  nu      = {nu:.4e} m²/s")
print(f"  Re/L    = {U_inf/nu:.1f} m⁻¹")

# -------------------------------------------------------------
# Blasius reference — solved in-place via shooting + RK45
# ODE: f''' + f·f'' = 0  with  η = y√(U∞/2νx)
# State vector: y = [f'', f', f]
# BCs: f(0)=0, f'(0)=0, f'(η→∞)=1
# → Cf = √2·f''(0)/√Re_x ≈ 0.664/√Re_x
# -------------------------------------------------------------
def _blasius_rhs(eta, y):
    return [-y[0] * y[2], y[0], y[1]]

def _fp_at_inf(fpp0_guess):
    sol = solve_ivp(_blasius_rhs, [0.0, 10.0], [fpp0_guess, 0.0, 0.0],
                    method='RK45', rtol=1e-10, atol=1e-12, max_step=0.01,
                    dense_output=True)
    return sol.y[1, -1]   # f'(10) → should equal 1

fpp0 = brentq(lambda g: _fp_at_inf(g) - 1.0, 0.4, 0.5, xtol=1e-12)

print(f"\nBlasius f''(0) = {fpp0:.6f}  →  Cf = {np.sqrt(2)*fpp0:.4f} / √Re_x")

# -------------------------------------------------------------
# Read wall.tec
# Cell-centered variables (indices 0-6):
#   0:y+  1:tauX  2:tauY  3:tauZ  4:pw  5:Tw  6:qw
# tauX is negative (wall drags fluid in -x), so take abs()
# -------------------------------------------------------------
xw_, yw_, zw_, varw_, _ = read_TEC("OUTPUT/wall.tec")

# Helpers to extract 1D nodal or cell arrays from a 3D block where the
# varying direction can be I (axis 0) or J (axis 1). We pick the first
# axis with size > 1 — this preserves original behaviour for typical
# I-varying data but also supports J-varying files.
def _nodal_1d(arr):
    s = arr.shape
    # consider primary 2D axes (0,1). If both vary, pick the one that
    # shows the largest variation for this coordinate (typical for
    # structured grids where x varies along one axis and y along the
    # other).
    if s[0] > 1 and s[1] > 1:
        col = arr[:, 0, 0]
        row = arr[0, :, 0]
        vcol = float(np.nanstd(col))
        vrow = float(np.nanstd(row))
        if vcol >= vrow:
            return col, 0
        else:
            return row, 1
    # fallback: pick first axis with size > 1
    axes = [i for i in range(len(s)) if s[i] > 1]
    if not axes:
        return arr.flatten(), 0
    ax = axes[0]
    if ax == 0:
        return arr[:, 0, 0], 0
    elif ax == 1:
        return arr[0, :, 0], 1
    else:
        return arr[0, 0, :], 2

def _cell_1d(arr, axis):
    # arr is a cell-centred field (shape ~ (ni-1, nj-1, nk-1)). Return
    # a 1D array along the given axis.
    if axis == 0:
        return arr[:, 0, 0]
    elif axis == 1:
        return arr[0, :, 0]
    else:
        return arr[0, 0, :]

def _cell_centres_from_nodal(nodal):
    if nodal.size > 1:
        return 0.5 * (nodal[:-1] + nodal[1:])
    return nodal.copy()

# Extract wall x nodal positions and detect which axis varies
xw_n, ax_xw = _nodal_1d(xw_[0])
xw_c = _cell_centres_from_nodal(xw_n)

yplus = _cell_1d(varw_[0][0], ax_xw)
tauX  = np.abs(_cell_1d(varw_[0][1], ax_xw))  # skin friction magnitude
Tw    = _cell_1d(varw_[0][5], ax_xw)

# Plate region only (adiabatic patch: x > 0)
mask    = xw_c > 0.0
x_plt   = xw_c[mask]
tau_plt = tauX[mask]
yp_plt  = yplus[mask]
Tw_plt  = Tw[mask]

# -------------------------------------------------------------
# Skin friction coefficient vs Blasius
# -------------------------------------------------------------
Re_x    = U_inf * x_plt / nu
Cf_mose = tau_plt / (0.5 * rho_inf * U_inf**2)
Cf_blas = np.sqrt(2.0) * fpp0 / np.sqrt(Re_x)

err_rel = np.abs(Cf_mose - Cf_blas) / Cf_blas * 100.0
err_rms = np.sqrt(np.mean(err_rel**2))
err_max = err_rel.max()
i_max_err = np.argmax(err_rel)
x_max_err = x_plt[i_max_err]

print(f"\nSkin friction vs Blasius")
print(f"  RMS relative error (all)    : {err_rms:.2f} %")
print(f"  Max relative error : {err_max:.2f} % at x = {x_max_err:.4f} m")

# Exclude leading edge (x < 0.05 m) from tolerance due to finite-difference
# effects near stagnation and coarse initial mesh. Blasius theory assumes
# infinitesimal leading edge; real numerics struggle there.
leading_edge_threshold = 0.05
valid_idx = x_plt >= leading_edge_threshold
leading_edge_idx = x_plt < leading_edge_threshold
downstream_idx = x_plt >= 0.1

if leading_edge_idx.sum() > 0:
    err_rms_le = np.sqrt(np.mean(err_rel[leading_edge_idx]**2))
    print(f"  → Leading edge (x < {leading_edge_threshold} m): RMS = {err_rms_le:.2f} %  [EXCLUDED]")
if valid_idx.sum() > 0:
    err_rms_valid = np.sqrt(np.mean(err_rel[valid_idx]**2))
    print(f"  → Valid region (x ≥ {leading_edge_threshold} m):   RMS = {err_rms_valid:.2f} %  [CHECKED]")
if downstream_idx.sum() > 0:
    err_rms_ds = np.sqrt(np.mean(err_rel[downstream_idx]**2))
    print(f"    ├─ Fully developed (x ≥ 0.1 m): RMS = {err_rms_ds:.2f} %")

status_excl = "PASS" if (valid_idx.sum() > 0 and err_rms_valid < 2.5) else "FAIL"
print(f"  STATUS (x ≥ {leading_edge_threshold} m) : {status_excl}  (tolerance 2.5 % RMS)")
print("\n  NOTE: Leading-edge errors are expected due to numerical stagnation")
print("        effects and coarse initial mesh. Theory assumes infinitesimal")
print("        leading edge; practical CFD has finite resolution.")

# -------------------------------------------------------------
# Velocity profile at x = x_ref vs Blasius
# -------------------------------------------------------------
x_ref = 0.29
# Find which axis carries the x variation in the field data
x_nodes, ax_x = _nodal_1d(x_[0])
x_field = _cell_centres_from_nodal(x_nodes)
i_x = int(np.argmin(np.abs(x_field - x_ref)))
x_act = float(x_field[i_x])
print(f"\nVelocity profile extracted at x = {x_act:.4f} m  (requested {x_ref} m)")

# Slice u profile and y nodal positions at the chosen x-cell index.
# For 2D data the other axis is the wall-normal direction.
def _slice_at_axis(arr, axis, idx):
    # arr expected shape (ni, nj, nk) or similar; return the slice where
    # the given axis is fixed to idx and the remaining primary axis is
    # returned as a 1D array (taking first k-index).
    if axis == 0:
        sl = arr[idx, :, 0]
    elif axis == 1:
        sl = arr[:, idx, 0]
    else:
        sl = arr[:, :, idx].flatten()
    return sl

u_block = var_[0][1]
u_prof = _slice_at_axis(u_block, ax_x, i_x)

y_block = y_[0]
y_nodes = _slice_at_axis(y_block, ax_x, i_x)
y_prof = _cell_centres_from_nodal(y_nodes)

# Detect wall location: wall is where velocity is minimum (no-slip condition).
# Sort so wall is at index 0, freestream at the end.
wall_idx = np.argmin(u_prof)
wall_y = y_prof[wall_idx]
# Create distance array from wall
sort_idx = np.argsort(np.abs(y_prof - wall_y))  # closest to wall first
y_prof   = y_prof[sort_idx]
u_prof   = u_prof[sort_idx]

# wall-normal distance measured from the wall
y_wall  = y_prof[0]
d_prof  = np.abs(y_prof - y_wall)

# Use the requested x_ref when computing the similarity variable η so the
# profile is compared at the nominal x location even if the data is from
# a nearby cell-centre. This applies the user-specified `x_ref` to array
# values used for the Blasius similarity scaling.
eta_num = d_prof * np.sqrt(U_inf / (2.0 * nu * x_ref))
u_norm  = u_prof / U_inf

# Dense Blasius profile via solve_ivp with dense output
_sol_blas = solve_ivp(_blasius_rhs, [0.0, 10.0], [fpp0, 0.0, 0.0],
                      method='RK45', rtol=1e-10, atol=1e-12, max_step=0.01,
                      dense_output=True)
eta_blas = np.linspace(0.0, max(float(eta_num.max()), 6.0), 400)
fp_blas  = _sol_blas.sol(np.clip(eta_blas, 0.0, 10.0))[1]

# -------------------------------------------------------------
# Plots
# -------------------------------------------------------------
plt.close('all')

# Figure 1: Skin friction coefficient
fig1 = plt.figure(figsize=(10, 6))
ax_cf = fig1.add_subplot(111)
ax_cf.plot(x_plt, Cf_blas, 'k-', lw=2.5, label='Blasius')
ax_cf.plot(x_plt, Cf_mose, 'go', ms=4, lw=2, label='MOSE')
#ax_cf.set_xscale('log')
#ax_cf.set_yscale('log')
ax_cf.set_xlim(-0.05, 0.35)
ax_cf.set_ylim(0.0, 0.08)
ax_cf.set_xlabel(r'$x$  [m]')
ax_cf.set_ylabel(r'$C_f$')
ax_cf.grid(True, alpha=0.3)
ax_cf.legend(loc='best')

if args.plot:
    plt.show(block=False)
else:
    fig1_path = Path(os.path.join(output_dir, "LP-Cf.svg"))
    fig1.savefig(fig1_path, dpi=150, bbox_inches='tight', transparent=True)
    print(f"Figure 1 saved to {fig1_path}")

# Figure 2: Velocity profile
fig2 = plt.figure(figsize=(10, 6))
ax_vel = fig2.add_subplot(111)
ax_vel.plot(fp_blas, eta_blas, 'k-', lw=2.5, label='Blasius')
ax_vel.plot(u_norm, eta_num, 'go', ms=5, lw=2, label=f'MOSE')
ax_vel.set_xlabel(r'$u\,/\,U_\infty$', fontsize=12)
ax_vel.set_ylabel(r'$\eta = y\,\sqrt{U_\infty\,/\,2\nu x}$', fontsize=12)
ax_vel.legend(loc='best', fontsize=11)
ax_vel.grid(True, alpha=0.3)
ax_vel.set_xlim(0.0, 1.05)
ax_vel.set_ylim(0.0, 9.0)

if args.plot:
    plt.show(block=False)
else:
    fig2_path = Path(os.path.join(output_dir, "LP-slice.svg"))
    fig2.savefig(fig2_path, dpi=150, bbox_inches='tight', transparent=True)
    print(f"Figure 2 saved to {fig2_path}")

# Figure 3: Flow field solution (velocity contours)
fig3 = plt.figure(figsize=(14, 6))
ax_field = fig3.add_subplot(111)

# Extract cell-centre coordinates and velocity magnitude
x_c = 0.5 * (x_[0][:-1, :-1, 0] + x_[0][1:, :-1, 0])
y_c = 0.5 * (y_[0][:-1, :-1, 0] + y_[0][:-1, 1:, 0])
u = var_[0][1][:, :, 0]  # streamwise velocity (cell-centred)
v = var_[0][2][:, :, 0]  # wall-normal velocity (cell-centred)
u_mag = np.sqrt(u**2 + v**2)/U_inf

# Create contour plot
levels = np.linspace(0.0, 1.01, 25)
cf = ax_field.contourf(x_c, y_c, u_mag, levels=levels, cmap='jet')
cbar = plt.colorbar(cf, ax=ax_field, orientation='horizontal', pad=-0.40, shrink=0.4, label=r'$|\vec{u}|$ [m/s]')
ax_field.set_xlabel('$x$ [m]', fontsize=12)
ax_field.set_ylabel('$y$ [m]', fontsize=12)
ax_field.grid(True, alpha=0.2)
ax_field.set_aspect('equal')

if args.plot:
    plt.show()
else:
    fig3_path = Path(os.path.join(output_dir, "LP-field.svg"))
    fig3.savefig(fig3_path, dpi=150, bbox_inches='tight', transparent=True)
    print(f"Figure 3 saved to {fig3_path}")
