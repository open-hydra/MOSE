import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import linspace
import numpy as np
from pathlib import Path
from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver
import sys
import os
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

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

# Define command-line arguments
parser = argparse.ArgumentParser(description='Plot CFD and Exact data')
parser.add_argument('--plot', action='store_true', help='Specify whether to plot the data or not')

# Parse the command-line arguments
args = parser.parse_args()
output_dir = "OUTPUT"
folder_name = os.path.basename(os.getcwd())

# Evaluate exact solution via ExactPack
xmin, xd0, xmax, t = 0.0, 0.5, 1.0, 0.15
rl, ul, pl, gl = 1.0,-2.0, 0.4, 1.4
rr, ur, pr, gr = 1.0, 2.0, 0.4, 1.4
solver_ig = IGEOS_Solver(xmin=xmin, xd0=xd0, xmax=xmax, t=t,
                            rl=rl, ul=ul, pl=pl, gl=gl,
                            rr=rr, ur=ur, pr=pr, gr=gr)
t = 0.15
nx = 100
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln1 = solver_ig(xe, t)
nx = 200
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln2 = solver_ig(xe, t)
nx = 400
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln4 = solver_ig(xe, t)
nx = 800
xe = linspace(0.5*1.0/nx, 1.0-0.5*1.0/nx, nx)
soln8 = solver_ig(xe, t)

# Read the CFD solution
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x1.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x1 = xc; var1 = var_[0]
cv1 = var1[7]/(var1[6]-1)
T1 = var1[5]
ei1 = cv1*var1[5]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x2.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x2 = xc; var2 = var_[0]
cv2 = var2[7]/(var2[6]-1)
T2 = var2[5]
ei2 = cv2*var2[5]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x4.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x4 = xc; var4 = var_[0]
cv4 = var4[7]/(var4[6]-1)
T4 = var4[5]
ei4 = cv4*var4[5]
[x_,y_,z_,var_,vnames] = read_TEC('OUTPUT/field_x8.tec')
xc = 0.5 * (x_[0][:-1] + x_[0][1:])
x8 = xc; var8 = var_[0]
cv8 = var8[7]/(var8[6]-1)
T8 = var8[5]
ei8 = cv8*var8[5]

# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Density
axs[0, 0].plot(xe, soln8['density'], label='Exact', color='black')
axs[0, 0].plot(x1[:, 0, 0], var1[0][:, 0, 0], label='CFD x1', color='red')
axs[0, 0].plot(x8[:, 0, 0], var8[0][:, 0, 0], label='CFD x8', color='green')
axs[0, 0].set_xlabel('x')
axs[0, 0].set_ylabel('Density')
axs[0, 0].grid(True, alpha=0.3)
axs[0, 0].legend()

# Velocity
axs[0, 1].plot(xe, soln8['velocity'], label='Exact', color='black')
axs[0, 1].plot(x1[:, 0, 0], var1[1][:, 0, 0], label='CFD x1', color='red')
axs[0, 1].plot(x8[:, 0, 0], var8[1][:, 0, 0], label='CFD x8', color='green')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('Velocity')
axs[0, 1].grid(True, alpha=0.3)
axs[0, 1].legend()

# Pressure
axs[1, 0].plot(xe, soln8['pressure'], label='Exact', color='black')
axs[1, 0].plot(x1[:, 0, 0], var1[4][:, 0, 0], label='CFD x1', color='red')
axs[1, 0].plot(x8[:, 0, 0], var8[4][:, 0, 0], label='CFD x8', color='green')
axs[1, 0].set_xlabel('x')
axs[1, 0].set_ylabel('Pressure')
axs[1, 0].grid(True, alpha=0.3)
axs[1, 0].legend()

# SIE
axs[1, 1].plot(xe, soln8['specific_internal_energy'], label='Exact', color='black')
axs[1, 1].plot(x1[:, 0, 0], cv1[:,0,0]*T1[:, 0, 0], label='CFD x1', color='red')
axs[1, 1].plot(x8[:, 0, 0], cv8[:,0,0]*T8[:, 0, 0], label='CFD x8', color='green')
axs[1, 1].set_xlabel('x')
axs[1, 1].set_ylabel('SIE')
axs[1, 1].grid(True, alpha=0.3)
axs[1, 1].legend()

# Adjust layout to prevent overlap of subplots
plt.tight_layout()

out_file = os.path.join(output_dir, f"{folder_name}.svg")
plt.savefig(out_file, bbox_inches="tight", transparent=True)

if args.plot:
    # Display the figure with four distinct plots
    plt.show()

iprobe = [8*10,8*20,8*30]
l2_norm = np.zeros((4, 4))
linf_norm = np.zeros((4, 4))

# Density
l2_norm[0,0] = np.linalg.norm(var1[0][:,0,0] - soln1['density'])
l2_norm[0,1] = np.linalg.norm(var2[0][:,0,0] - soln2['density'])
l2_norm[0,2] = np.linalg.norm(var4[0][:,0,0] - soln4['density'])
l2_norm[0,3] = np.linalg.norm(var8[0][:,0,0] - soln8['density'])
linf_norm[0,0] = np.max(np.abs(var1[0][:,0,0] - soln1['density']))
linf_norm[0,1] = np.max(np.abs(var2[0][:,0,0] - soln2['density']))
linf_norm[0,2] = np.max(np.abs(var4[0][:,0,0] - soln4['density']))
linf_norm[0,3] = np.max(np.abs(var8[0][:,0,0] - soln8['density']))
for i in range(len(iprobe)):
    perr = np.abs((var8[0][iprobe[i],0,0] - soln8['density'][iprobe[i]]) / soln8['density'][iprobe[i]]) * 100
    if (perr > 1.0):
        print('Density')
        print(iprobe[i])
        print(var8[0][iprobe[i],0,0], soln8['density'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)

# Velocity
l2_norm[1,0] = np.linalg.norm(var1[1][:,0,0] - soln1['velocity'])
l2_norm[1,1] = np.linalg.norm(var2[1][:,0,0] - soln2['velocity'])
l2_norm[1,2] = np.linalg.norm(var4[1][:,0,0] - soln4['velocity'])
l2_norm[1,3] = np.linalg.norm(var8[1][:,0,0] - soln8['velocity'])
linf_norm[1,0] = np.max(np.abs(var1[0][:,0,0] - soln1['velocity']))
linf_norm[1,1] = np.max(np.abs(var2[0][:,0,0] - soln2['velocity']))
linf_norm[1,2] = np.max(np.abs(var4[0][:,0,0] - soln4['velocity']))
linf_norm[1,3] = np.max(np.abs(var8[0][:,0,0] - soln8['velocity']))
for i in range(len(iprobe)):
    perr = np.abs((var8[1][iprobe[i],0,0] - soln8['velocity'][iprobe[i]]) / (soln8['velocity'][iprobe[i]]+1e-10)) * 100
    if (perr > 1.0):
        print('Velocity')
        print(iprobe[i])
        print(var8[1][iprobe[i],0,0], soln8['velocity'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)

# Pressure
l2_norm[2,0] = np.linalg.norm(var1[4][:,0,0] - soln1['pressure'])
l2_norm[2,1] = np.linalg.norm(var2[4][:,0,0] - soln2['pressure'])
l2_norm[2,2] = np.linalg.norm(var4[4][:,0,0] - soln4['pressure'])
l2_norm[2,3] = np.linalg.norm(var8[4][:,0,0] - soln8['pressure'])
linf_norm[2,0] = np.max(np.abs(var1[0][:,0,0] - soln1['pressure']))
linf_norm[2,1] = np.max(np.abs(var2[0][:,0,0] - soln2['pressure']))
linf_norm[2,2] = np.max(np.abs(var4[0][:,0,0] - soln4['pressure']))
linf_norm[2,3] = np.max(np.abs(var8[0][:,0,0] - soln8['pressure']))
for i in range(len(iprobe)):
    perr = np.abs((var8[4][iprobe[i],0,0] - soln8['pressure'][iprobe[i]]) / soln8['pressure'][iprobe[i]]) * 100
    if (perr > 1.5):
        print('Pressure')
        print(iprobe[i])
        print(var8[4][iprobe[i],0,0], soln8['pressure'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.5)
        sys.exit(1)

# SIE
l2_norm[3,0] = np.linalg.norm(ei1[:,0,0] - soln1['specific_internal_energy'])
l2_norm[3,1] = np.linalg.norm(ei2[:,0,0] - soln2['specific_internal_energy'])
l2_norm[3,2] = np.linalg.norm(ei4[:,0,0] - soln4['specific_internal_energy'])
l2_norm[3,3] = np.linalg.norm(ei8[:,0,0] - soln8['specific_internal_energy'])
linf_norm[3,0] = np.max(np.abs(var1[0][:,0,0] - soln1['specific_internal_energy']))
linf_norm[3,1] = np.max(np.abs(var2[0][:,0,0] - soln2['specific_internal_energy']))
linf_norm[3,2] = np.max(np.abs(var4[0][:,0,0] - soln4['specific_internal_energy']))
linf_norm[3,3] = np.max(np.abs(var8[0][:,0,0] - soln8['specific_internal_energy']))
for i in range(len(iprobe)):
    perr = np.abs((ei8[iprobe[i],0,0] - soln8['specific_internal_energy'][iprobe[i]]) / soln8['specific_internal_energy'][iprobe[i]]) * 100
    if (perr > 1.0):
        print('SIE')
        print(iprobe[i])
        print(ei8[iprobe[i],0,0], soln8['specific_internal_energy'][iprobe[i]])
        print('Error = ',perr, 'threshold = ',1.0)
        sys.exit(1)


# # Set x-axis values
# x_values = [1, 2, 4, 8]

# # Plot each component of l2_norm
# for i in range(1):
#     plt.plot(x_values, linf_norm[i], marker='o', label=f'Component {i+1}')

# # Add labels and legend
# plt.xlabel('X-axis (log scale)')
# plt.ylabel('Y-axis')
# plt.title('Plot of Components of l2_norm')
# plt.legend()

# # Set logarithmic scale on the x-axis
# plt.xscale('log')

# # Set x-ticks and labels
# plt.xticks(x_values, x_values)

# # Show the plot
# plt.show()
