import argparse
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

output_dir = "./"

# Try to make the script work even when conda is not activated
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
        "ORION package not found. Run 'conda activate base' or install ORION in the current Python environment."
    )

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Validation script')
parser.add_argument('--plot', action='store_true', help='Plot pressure contours and ratio comparison')
args = parser.parse_args()


[x_, y_, z_, var_, vnames] = read_TEC("frozen/reference/field.tec")
frozen = {
    'x': x_[0],
    'y': y_[0],
    'p': var_[0][16],
    'rho': sum(var_[0][0:12]),
    'u': var_[0][13],
    'v': var_[0][14],
    'T': var_[0][17],
    'g': var_[0][18],
    'R': var_[0][19],
}

[x_, y_, z_, var_, vnames] = read_TEC("finite-rate/reference/field.tec")
finite_rate = {
    'x': x_[0],
    'y': y_[0],
    'p': var_[0][22],
    'rho': sum(var_[0][0:18]),
    'u': var_[0][19],
    'v': var_[0][20],
    'T': var_[0][23],
    'g': var_[0][24],
    'R': var_[0][25],
}

[x_, y_, z_, var_, vnames] = read_TEC("equilibrium/reference/field.tec")
equilibrium = {
    'x': x_[0],
    'y': y_[0],
    'p': var_[0][16],
    'rho': sum(var_[0][0:12]),
    'u': var_[0][13],
    'v': var_[0][14],
    'T': var_[0][17],
    'g': var_[0][18],
    'R': var_[0][19],
}


# ----------------------------------------------------------
# Quasi-1D comparison: area-weighted averages vs CEA
# ----------------------------------------------------------
# CEA reference: Merlin 1D, RP-1/LOX O/F=2.35
# Source: CEA.out (Pin=108 bar, O/F=2.35)
# Stations: sub Ae/At=2&1.5, throat, sup Ae/At=2,4,8,16
P0_cea = 108.00    # BAR (stagnation)
T0_cea = 3644.12   # K   (stagnation)
#                     sub2.0    throat   sup2.0   sup4.0    sup8.0   sup16.0
AeAt_cea  = np.array([ 2.0,        1.0,     2.0,     4.0,      8.0,      16.0])
M_cea_fz  = np.array([ 0.312,    1.000,   2.066,   2.650,    3.184,     3.715])
M_cea_eq  = np.array([ 0.314,    1.000,   2.018,    2.549,   3.037,     3.526])
P_cea_fz  = np.array([101.84,   60.606,  12.676,   4.4664,  1.6840,   0.65259])
P_cea_eq  = np.array([102.07,   62.026,  13.783,   5.0021,  1.9243,   0.76119])
T_cea_fz  = np.array([3606.55, 3289.48, 2477.34,  2037.39,  1686.54, 1393.69])
T_cea_eq  = np.array([3621.79, 3428.71, 2862.18,  2472.84,  2113.89,  1797.26])
mw_cea_fz = np.array([22.468,   22.468,  22.468,   22.468,   22.468,   22.468])
mw_cea_eq = np.array([22.495,   22.725,  23.252,   23.424,   23.480,   23.492])
# Cp values from CEA (KJ/(kg·K)) at reference stations
cp_cea_fz = np.array([2.0970,   2.0768,  2.0043,   1.9430,   1.8757,   1.8023])
# For equilibrium case, use "WITH FROZEN REACTIONS" Cp from CEA (frozen chemistry path)
cp_cea_eq = np.array([2.0981,   2.0877,  2.0482,   2.0094,   1.9621,   1.9090])
# Compute gamma = cp/cv = cp/(cp - R) where R = Runi/M
Runi = 8314.51  # J/(kmol·K)
R_cea_fz = Runi / mw_cea_fz  # J/(kg·K)
R_cea_eq = Runi / mw_cea_eq  # J/(kg·K)
g_cea_fz = (cp_cea_fz * 1000) / (cp_cea_fz * 1000 - R_cea_fz)  # convert Cp to J/(kg·K)
g_cea_eq = (cp_cea_eq * 1000) / (cp_cea_eq * 1000 - R_cea_eq)
# (These p/p0 and T/T0 ratios are no longer used since we replaced pressure subplot)
# but kept here for potential future use
pp0_cea_fz   = P_cea_fz / P0_cea
TT0_cea_fz   = T_cea_fz / T0_cea
pp0_cea_eq   = P_cea_eq / P0_cea
TT0_cea_eq   = T_cea_eq / T0_cea

# -------------------------------------------------------------
# Geometry and mass-flowrate-weighted averaging
# -------------------------------------------------------------

xn = frozen['x'][:, :, 0]; yn = frozen['y'][:, :, 0]

# Cell-face radii (axisymmetric: y = radius).
# Nodes: (I, J) = (121, 41); cells: (120, 40).
r_lo = 0.5 * (yn[:120, :40] + yn[1:, :40])   # inner radial face, shape (120, 40)
r_hi = 0.5 * (yn[:120, 1:]  + yn[1:, 1:])    # outer radial face, shape (120, 40)
w_ring = r_hi**2 - r_lo**2                    # proportional to ring area

def mfr_avg_1d(rho, u, q2d):
    """Mass-flow-rate weighted average (ρ·u weight); q2d shape (Ncx, Ncy) -> (Ncx,)."""
    w = rho * u * w_ring
    return np.sum(q2d * w, axis=1) / np.sum(w, axis=1)

# Cell-centre x coordinate
x_c = 0.5 * (xn[:120, 0] + xn[1:, 0])   # shape (120,)

# 2D field slices (k=0)
rho2d_fz = frozen['rho'][:, :, 0]
p2d_fz   = frozen['p'][:, :, 0]
T2d_fz   = frozen['T'][:, :, 0]
g2d_fz   = frozen['g'][:, :, 0]
u2d_fz   = frozen['u'][:, :, 0]
v2d_fz   = frozen['v'][:, :, 0]
a2d_fz   = np.sqrt(g2d_fz * p2d_fz / rho2d_fz)
M2d_fz   = np.sqrt(u2d_fz**2 + v2d_fz**2) / a2d_fz
mw2d_fz  = 8314.51/frozen['R'][:, :, 0]

rho2d_fr = finite_rate['rho'][:, :, 0]
p2d_fr   = finite_rate['p'][:, :, 0]
T2d_fr   = finite_rate['T'][:, :, 0]
g2d_fr   = finite_rate['g'][:, :, 0]
u2d_fr   = finite_rate['u'][:, :, 0]
v2d_fr   = finite_rate['v'][:, :, 0]
a2d_fr   = np.sqrt(g2d_fr * p2d_fr / rho2d_fr)
M2d_fr   = np.sqrt(u2d_fr**2 + v2d_fr**2) / a2d_fr
mw2d_fr  = 8314.51/finite_rate['R'][:, :, 0]

rho2d_eq = equilibrium['rho'][:, :, 0]
p2d_eq   = equilibrium['p'][:, :, 0]
T2d_eq   = equilibrium['T'][:, :, 0]
g2d_eq   = equilibrium['g'][:, :, 0]
u2d_eq   = equilibrium['u'][:, :, 0]
v2d_eq   = equilibrium['v'][:, :, 0]
a2d_eq   = np.sqrt(g2d_eq * p2d_eq / rho2d_eq)
M2d_eq   = np.sqrt(u2d_eq**2 + v2d_eq**2) / a2d_eq
mw2d_eq  = 8314.51/equilibrium['R'][:, :, 0]

# Cell-centre coordinates for Mach contour
x_cent = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
y_cent = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])

# Mach field plot
fig1, ax1 = plt.subplots(figsize=(10, 4))
contour = ax1.contourf(x_cent, y_cent, M2d_fr, levels=50, cmap='turbo')
cax = inset_axes(ax1, width="40%", height="5%", loc='upper center')
cbar = plt.colorbar(contour, cax=cax, orientation='horizontal')
cbar.set_label("Mach number")
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_aspect('equal', adjustable='box')

# 1D averages and axis values
p_avg_fz  = mfr_avg_1d(rho2d_fz, u2d_fz, p2d_fz)
T_avg_fz  = mfr_avg_1d(rho2d_fz, u2d_fz, T2d_fz)
M_avg_fz  = mfr_avg_1d(rho2d_fz, u2d_fz, M2d_fz)
g_avg_fz  = mfr_avg_1d(rho2d_fz, u2d_fz, g2d_fz)
mv_avg_fz = mfr_avg_1d(rho2d_fz, u2d_fz, mw2d_fz)

p_avg_fr  = mfr_avg_1d(rho2d_fr, u2d_fr, p2d_fr)
T_avg_fr  = mfr_avg_1d(rho2d_fr, u2d_fr, T2d_fr)
M_avg_fr  = mfr_avg_1d(rho2d_fr, u2d_fr, M2d_fr)
g_avg_fr  = mfr_avg_1d(rho2d_fr, u2d_fr, g2d_fr)
mv_avg_fr = mfr_avg_1d(rho2d_fr, u2d_fr, mw2d_fr)

p_avg_eq  = mfr_avg_1d(rho2d_eq, u2d_eq, p2d_eq)
T_avg_eq  = mfr_avg_1d(rho2d_eq, u2d_eq, T2d_eq)
M_avg_eq  = mfr_avg_1d(rho2d_eq, u2d_eq, M2d_eq)
g_avg_eq  = mfr_avg_1d(rho2d_eq, u2d_eq, g2d_eq)
mv_avg_eq = mfr_avg_1d(rho2d_eq, u2d_eq, mw2d_eq)

# Inlet stagnation recovered via isentropic relations
T0_sim = T0_cea

# CEA x-positions: locate each Ae/At station in the mesh geometry.
# Subsonic stations live in the converging section, supersonic in the diverging.
R_wall      = r_hi[:, -1]
i_throat    = int(np.argmin(R_wall))
R_t         = R_wall[i_throat]
AeAt_1d     = R_wall**2 / R_t**2            # local area ratio at each axial cell
sub_idx     = np.where(x_c < x_c[i_throat])[0]
sup_idx     = np.where(x_c >= x_c[i_throat])[0]

x_cea_list = []
for k, AeAt in enumerate(AeAt_cea):
    if M_cea_fz[k] < 1.0:                      # subsonic → converging section
        region = sub_idx
    else:                                    # sonic/supersonic → diverging section
        region = sup_idx
    i = region[np.argmin(np.abs(AeAt_1d[region] - AeAt))]
    x_cea_list.append(x_c[i])
x_cea = np.array(x_cea_list)

fig3, axes = plt.subplots(2, 2, figsize=(15, 11), sharex=True)

# Mach
ax = axes[0, 0]
ax.plot(x_c, M_avg_fz,  label='MOSE Frozen', linewidth=2)
ax.plot(x_c, M_avg_fr,  color='green', label='MOSE Finite-rate', linewidth=2)
ax.plot(x_c, M_avg_eq,  label='MOSE Equilibrium', linewidth=2)
ax.scatter(x_cea, M_cea_fz, color='blue', zorder=5, label='CEA Frozen', s=60)
ax.scatter(x_cea, M_cea_eq, color='red', zorder=5, label='CEA Equilibrium', s=60)
ax.set_ylabel('Mach')
ax.set_xlabel('x  [m]')
ax.legend()
ax.grid(True, alpha=0.3)

# gamma
ax = axes[1, 0]
ax.plot(x_c, g_avg_fz, label='MOSE Frozen', linewidth=2)
ax.plot(x_c, g_avg_fr, color='green', label='MOSE Finite-rate', linewidth=2)
ax.plot(x_c, g_avg_eq, color='orange', label='MOSE Equilibrium', linewidth=2)
ax.scatter(x_cea, g_cea_fz, color='blue', zorder=5, label='CEA Frozen', s=60)
ax.scatter(x_cea, g_cea_eq, color='red', zorder=5, label='CEA Equilibrium', s=60)
ax.set_ylabel('Specific heat ratio')
ax.set_xlabel('x  [m]')
ax.legend()
ax.grid(True, alpha=0.3)

# T / T0
ax = axes[0, 1]
ax.plot(x_c, T_avg_fz  / T0_sim, label='MOSE Frozen', linewidth=2)
ax.plot(x_c, T_avg_fr  / T0_sim, color='green', label='MOSE Finite-rate', linewidth=2)
ax.plot(x_c, T_avg_eq  / T0_sim, color='orange', label='MOSE Equilibrium', linewidth=2)
ax.scatter(x_cea, TT0_cea_fz, color='blue', zorder=5, label='CEA Frozen', s=60)
ax.scatter(x_cea, TT0_cea_eq, color='red', zorder=5, label='CEA Equilibrium', s=60)
ax.set_ylabel('Temperature ratio T/T0')
ax.set_xlabel('x  [m]')
ax.legend()
ax.grid(True, alpha=0.3)

# Molecular weight
ax = axes[1, 1]
ax.plot(x_c, mv_avg_fz, label='MOSE Frozen', linewidth=2)
ax.plot(x_c, mv_avg_fr, color='green', label='MOSE Finite-rate', linewidth=2)
ax.plot(x_c, mv_avg_eq, color='orange', label='MOSE Equilibrium', linewidth=2)
ax.scatter(x_cea, mw_cea_fz, color='blue', zorder=5, label='CEA Frozen', s=60)
ax.scatter(x_cea, mw_cea_eq, color='red', zorder=5, label='CEA Equilibrium', s=60)
ax.set_ylabel('Molecular weight [g/mol]')
ax.set_xlabel('x  [m]')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig1.tight_layout()

if args.plot:
    plt.show()
else:
    os.makedirs(output_dir, exist_ok=True)

    out_path = os.path.join(output_dir, "Merlin-field.svg")
    fig1.savefig(out_path, transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")

    out_path = os.path.join(output_dir, "Merlin-XY.svg")
    plt.savefig(out_path, bbox_inches="tight",dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")
