"""
Woodward-Colella forward step — MOSE vs OpenFOAM comparison

Loads:
  - MOSE output : OUTPUT/field.tec  (Tecplot BLOCK, 2 structured zones)
  - OpenFOAM ref: reference/OpenFOAM/forwardStep.0000.{mesh,p,T,U}
                  (EnSight Gold binary, single-layer hexa8 mesh)

Interpolates the OpenFOAM fields onto the MOSE cell centres with a
nearest-neighbour scheme, then computes relative L2 and L-inf errors
for pressure, temperature, and velocity magnitude.

Usage:
    python verify.py           # prints error table, saves figure to OUTPUT/verify.png
    python verify.py --plot    # shows interactive figure instead
"""

import argparse
import struct
import sys
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import NearestNDInterpolator
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
parser = argparse.ArgumentParser(
    description="Woodward-Colella forward step — MOSE vs OpenFOAM"
)
parser.add_argument("--plot", action="store_true", help="Show interactive plots")
args = parser.parse_args()

# =============================================================
# 1.  Load MOSE solution (OUTPUT/field.tec)
# =============================================================
# Variable layout (x,y,z stripped by read_TEC):
#   0=rho  1=u   2=v   3=w   4=p   5=T   6=gamma  7=R
x_, y_, z_, var_, _ = read_TEC("OUTPUT/field.tec")

mose_xc_parts, mose_yc_parts = [], []
mose_p_parts, mose_T_parts, mose_umag_parts = [], [], []

for zone in range(len(x_)):
    xn = np.asarray(x_[zone])      # (I, J, 2) – nodal coordinates
    yn = np.asarray(y_[zone])
    vn = np.asarray(var_[zone])    # (8, I-1, J-1, 1) – cell-centred

    # 2-D cell centres (face-averaged from k=0 nodal layer)
    xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                 xn[:-1,  1:, 0] + xn[1:,  1:, 0])   # (I-1, J-1)
    yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                 yn[:-1,  1:, 0] + yn[1:,  1:, 0])

    mose_xc_parts.append(xc.ravel())
    mose_yc_parts.append(yc.ravel())
    mose_p_parts.append(vn[4, :, :, 0].ravel())
    mose_T_parts.append(vn[5, :, :, 0].ravel())
    mose_umag_parts.append(
        np.sqrt(vn[1, :, :, 0]**2 + vn[2, :, :, 0]**2 + vn[3, :, :, 0]**2).ravel()
    )

mose_xc   = np.concatenate(mose_xc_parts)
mose_yc   = np.concatenate(mose_yc_parts)
mose_p    = np.concatenate(mose_p_parts)
mose_T    = np.concatenate(mose_T_parts)
mose_umag = np.concatenate(mose_umag_parts)

# =============================================================
# 2.  Load OpenFOAM EnSight Gold reference
# =============================================================
REF_DIR = Path("reference/OpenFOAM")


def _es_str(f):
    """Read one 80-byte EnSight Gold string record."""
    return f.read(80).decode("ascii").strip("\x00").strip()


def read_ensight_mesh(path):
    """Parse EnSight Gold binary geometry -> node coords + hexa8 connectivity."""
    with open(path, "rb") as f:
        for _ in range(5):
            _es_str(f)                          # C binary / desc / author / node_id / elem_id
        _es_str(f)                              # "part"
        struct.unpack("i", f.read(4))[0]        # part number
        _es_str(f)                              # part description
        _es_str(f)                              # "coordinates"
        nn = struct.unpack("i", f.read(4))[0]
        xn   = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        yn   = np.frombuffer(f.read(4 * nn), dtype=np.float32)
        f.read(4 * nn)                          # z coords — ignored (2-D problem)
        _es_str(f)                              # "hexa8"
        ne   = struct.unpack("i", f.read(4))[0]
        conn = np.frombuffer(f.read(8 * 4 * ne), dtype=np.int32).reshape(ne, 8)
    return xn, yn, conn, ne


def read_ensight_scalar(path, ne):
    """Parse EnSight Gold per-element scalar field."""
    with open(path, "rb") as f:
        _es_str(f)                              # field description
        _es_str(f)                              # "part"
        struct.unpack("i", f.read(4))[0]        # part number
        _es_str(f)                              # element type
        data = np.frombuffer(f.read(4 * ne), dtype=np.float32)
    return data.astype(float)


def read_ensight_vector(path, ne):
    """Parse EnSight Gold per-element vector field -> (3, ne)."""
    with open(path, "rb") as f:
        _es_str(f)
        _es_str(f)
        struct.unpack("i", f.read(4))[0]
        _es_str(f)
        data = np.frombuffer(f.read(3 * 4 * ne), dtype=np.float32).reshape(3, ne)
    return data.astype(float)


xn_of, yn_of, conn_of, ne_of = read_ensight_mesh(REF_DIR / "forwardStep.0000.mesh")

of_p    = read_ensight_scalar(REF_DIR / "forwardStep.0000.p", ne_of)
of_T    = read_ensight_scalar(REF_DIR / "forwardStep.0000.T", ne_of)
of_U    = read_ensight_vector(REF_DIR / "forwardStep.0000.U", ne_of)

# Element centroids (connectivity is 1-indexed)
of_xc   = xn_of[conn_of - 1].mean(axis=1)
of_yc   = yn_of[conn_of - 1].mean(axis=1)
of_umag = np.sqrt(of_U[0]**2 + of_U[1]**2 + of_U[2]**2)

# =============================================================
# 3.  Interpolate OpenFOAM fields onto MOSE cell centres
# =============================================================
of_pts = np.column_stack([of_xc, of_yc])

interp_p    = NearestNDInterpolator(of_pts, of_p)
interp_T    = NearestNDInterpolator(of_pts, of_T)
interp_umag = NearestNDInterpolator(of_pts, of_umag)

mose_pts = np.column_stack([mose_xc, mose_yc])
of_p_at_mose    = interp_p(mose_pts)
of_T_at_mose    = interp_T(mose_pts)
of_umag_at_mose = interp_umag(mose_pts)

# =============================================================
# 4.  Error metrics
# =============================================================
def l2_rel(a, b):
    """Relative L2 error: ||a-b||_2 / ||b||_2."""
    return np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2))


def linf_rel(a, b):
    """Relative L-inf error: max|a-b| / RMS(b)."""
    return np.max(np.abs(a - b)) / np.sqrt(np.mean(b ** 2))


print("Woodward-Colella forward step — MOSE vs OpenFOAM")
print("=" * 52)
print(f"  MOSE cells   : {len(mose_xc):>6d}")
print(f"  OpenFOAM cells: {ne_of:>5d}")
print()
print(f"  {'Variable':<10}  {'L2 rel [%]':>12}  {'Linf rel [%]':>14}")
print("  " + "-" * 40)
for vname, mv, ov in [
    ("p",   mose_p,    of_p_at_mose),
    ("T",   mose_T,    of_T_at_mose),
    ("|U|", mose_umag, of_umag_at_mose),
]:
    print(f"  {vname:<10}  {100*l2_rel(mv,ov):>12.2f}  {100*linf_rel(mv,ov):>14.2f}")
print()

# =============================================================
# 5.  Plots
# =============================================================
fig, axes = plt.subplots(3, 2, figsize=(14, 9))

# Helper: rebuild structured arrays for one MOSE zone
def _zone_arrays(zone_idx, var_idx):
    xn = np.asarray(x_[zone_idx])
    yn = np.asarray(y_[zone_idx])
    vn = np.asarray(var_[zone_idx])
    return xn[:, :, 0], yn[:, :, 0], vn[var_idx, :, :, 0]


STEP_X = [0.0, 0.6, 0.6]
STEP_Y = [0.2, 0.2, 0.0]

for row, (vname, var_idx, of_vals) in enumerate([
    ("p",   4, of_p),
    ("T",   5, of_T),
    ("|U|", None, of_umag),
]):
    # Collect MOSE cell centres and values for scatter plot
    mose_xc_all, mose_yc_all, mose_vals_all = [], [], []
    for zone in range(len(x_)):
        xn = np.asarray(x_[zone])
        yn = np.asarray(y_[zone])
        vn = np.asarray(var_[zone])
        xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                     xn[:-1,  1:, 0] + xn[1:,  1:, 0])
        yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                     yn[:-1,  1:, 0] + yn[1:,  1:, 0])
        if var_idx is not None:
            vals = vn[var_idx, :, :, 0]
        else:
            vals = np.sqrt(vn[1,:,:,0]**2 + vn[2,:,:,0]**2 + vn[3,:,:,0]**2)
        mose_xc_all.append(xc.ravel())
        mose_yc_all.append(yc.ravel())
        mose_vals_all.append(vals.ravel())
    
    mose_xc_plot = np.concatenate(mose_xc_all)
    mose_yc_plot = np.concatenate(mose_yc_all)
    mose_vals_plot = np.concatenate(mose_vals_all)
    
    # shared colour scale
    vmin = min(mose_vals_plot.min(), of_vals.min())
    vmax = max(mose_vals_plot.max(), of_vals.max())

    # ---- MOSE scatter plot ----
    ax = axes[row, 0]
    sc = ax.scatter(mose_xc_plot, mose_yc_plot, c=mose_vals_plot, s=4, cmap="turbo", vmin=vmin, vmax=vmax,
                    linewidths=0, rasterized=True)
    fig.colorbar(sc, ax=ax, label=vname)
    ax.axvline(0.6, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(f"MOSE — {vname}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)

    # ---- OpenFOAM scatter plot ----
    ax = axes[row, 1]
    sc = ax.scatter(of_xc, of_yc, c=of_vals, s=4, cmap="turbo", vmin=vmin, vmax=vmax,
                    linewidths=0,rasterized=True)
    fig.colorbar(sc, ax=ax, label=vname)
    ax.axvline(0.6, color="gray", lw=0.8, ls=":", alpha=0.7)
    ax.set_title(f"OpenFOAM — {vname}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)

plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"WC-fields.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True)
    print(f"  Figure saved to {out_path}")

# =============================================================
# 6.  Density profile along y = 0.5 horizontal line
# =============================================================
Y_SLICE = 0.5

# --- MOSE: collect cells from both zones nearest to y = Y_SLICE ---
mose_x_slice, mose_rho_slice = [], []
for zone in range(len(x_)):
    xn = np.asarray(x_[zone])
    yn = np.asarray(y_[zone])
    vn = np.asarray(var_[zone])
    xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                 xn[:-1,  1:, 0] + xn[1:,  1:, 0])   # (I-1, J-1)
    yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                 yn[:-1,  1:, 0] + yn[1:,  1:, 0])
    # For each column i, pick the row j closest to Y_SLICE
    j_best = np.argmin(np.abs(yc[0, :] - Y_SLICE))
    # Only include if y is within half a cell height
    dy = np.abs(yc[0, j_best] - Y_SLICE)
    dy_cell = np.abs(yc[0, min(j_best + 1, yc.shape[1] - 1)] - yc[0, j_best])
    if dy <= dy_cell:
        mose_x_slice.append(xc[:, j_best])
        mose_rho_slice.append(vn[0, :, j_best, 0])

mose_x_slice   = np.concatenate(mose_x_slice)
mose_rho_slice = np.concatenate(mose_rho_slice)
sort_idx = np.argsort(mose_x_slice)
mose_x_slice   = mose_x_slice[sort_idx]
mose_rho_slice = mose_rho_slice[sort_idx]

# --- OpenFOAM: cells in a band around Y_SLICE, then sort by x ---
R_gas = 1.0 / 1.4          # same non-dimensional R used by both solvers
of_rho = of_p / (R_gas * of_T)

# use a band of half the typical cell height (~1/128 of domain height)
band = 1.0 / 128.0
mask_of = np.abs(of_yc - Y_SLICE) <= band
# average duplicate x positions (cells at same x, different z layer)
of_x_raw   = of_xc[mask_of]
of_rho_raw = of_rho[mask_of]
sort_of = np.argsort(of_x_raw)
of_x_slice   = of_x_raw[sort_of]
of_rho_slice = of_rho_raw[sort_of]

fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.plot(mose_x_slice, mose_rho_slice, label="MOSE",     lw=2.5)
ax2.plot(of_x_slice,   of_rho_slice,   label="OpenFOAM", lw=2.5, ls="--")
ax2.set_xlabel("x")
ax2.set_ylabel(r"$\rho$")
ax2.set_xlim(0.0, 3.0)
ax2.legend()
ax2.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path2 = os.path.join(output_dir, f"WC-slice.svg")
    fig2.savefig(out_path2, bbox_inches="tight", transparent=True)
    print(f"  Figure saved to {out_path2}")