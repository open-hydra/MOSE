import argparse
from numpy import linspace
import numpy as np
from pathlib import Path
import sys
import os
import warnings
warnings.filterwarnings("ignore")

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

# Read the exact solution in reference/exact.txt (x, rho, u, p, sie)
ref = np.loadtxt('reference/exact.txt')
xe   = ref[:, 0]
rhoe = ref[:, 1]
ue   = ref[:, 2]
pe   = ref[:, 3]
siee = ref[:, 4]

# Read the CFD solution
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field.tec')
xc  = 0.5 * (x_[0][:-1] + x_[0][1:])
rho = var_[0][0][:, 0, 0]
u   = var_[0][1][:, 0, 0]
p   = var_[0][4][:, 0, 0]
T   = var_[0][5][:, 0, 0]
gam = var_[0][6][:, 0, 0]
cp  = var_[0][7][:, 0, 0]
cv  = cp / (gam - 1)
sie = cv * T

# Interpolate exact solution onto CFD grid (handles different mesh sizes)
rhoe_i = np.interp(xc[:, 0, 0], xe, rhoe)
ue_i   = np.interp(xc[:, 0, 0], xe, ue)
pe_i   = np.interp(xc[:, 0, 0], xe, pe)
siee_i = np.interp(xc[:, 0, 0], xe, siee)

# Probe positions (away from discontinuities)
ncells = len(rho)
iprobe = [int(ncells * f) for f in [0.20, 0.50, 0.60, 0.80, 0.90]]
threshold = 1.0  # percent

# Check density
for i in iprobe:
    perr = abs((rho[i] - rhoe_i[i]) / rhoe_i[i]) * 100
    if perr > threshold:
        print(f'Density at cell {i}:  CFD={rho[i]:.6f}  Exact={rhoe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

# Check velocity
for i in iprobe:
    perr = abs((u[i] - ue_i[i]) / (ue_i[i] + 1e-10)) * 100
    if perr > threshold:
        print(f'Velocity at cell {i}:  CFD={u[i]:.6f}  Exact={ue_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

# Check pressure
for i in iprobe:
    perr = abs((p[i] - pe_i[i]) / pe_i[i]) * 100
    if perr > threshold:
        print(f'Pressure at cell {i}:  CFD={p[i]:.6f}  Exact={pe_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

# Check SIE
for i in iprobe:
    perr = abs((sie[i] - siee_i[i]) / siee_i[i]) * 100
    if perr > threshold:
        print(f'SIE at cell {i}:  CFD={sie[i]:.6f}  Exact={siee_i[i]:.6f}  Error={perr:.2f}%  threshold={threshold}%')
        sys.exit(1)

print('Correct!')
