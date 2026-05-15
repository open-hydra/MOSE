import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore")

#Configure mpl for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

# r range
r = np.linspace(0, 3, 600)

# ----- Limiters -----
def minmod(r):
    return np.maximum(0, np.minimum(1, r))

def van_leer(r):
    return (r + np.abs(r)) / (1 + np.abs(r))

def van_albada(r):
    return (r**2 + r) / (1 + r**2)

def mc(r):
    return np.maximum(0, np.minimum(np.minimum(2*r, 0.5*(1+r)), 2))

def superbee(r):
    return np.maximum(0, np.maximum(np.minimum(2*r, 1), np.minimum(r, 2)))

# ----- TVD bounds -----
phi_lower = minmod(r)
phi_upper = superbee(r)

# Plot
plt.figure(figsize=(7, 5))
ax = plt.gca()

# Shade TRUE TVD region
ax.fill_between(r, phi_lower, phi_upper, alpha=0.1)

# Plot limiters
ax.plot(r, minmod(r), label="MinMod")
ax.plot(r, van_leer(r), label="Van Leer")
ax.plot(r, van_albada(r), label="Van Albada")
ax.plot(r, mc(r), label="MC")
ax.plot(r, superbee(r), label="Superbee")

# Labels
ax.set_xlabel("r")
# Latex rendering for phi
ax.set_ylabel("$\phi(r)$")

# Limits
ax.set_xlim(0, 3)
ax.set_ylim(0, 2.5)

# Grid & legend
ax.grid(True, alpha=0.3)
ax.legend()

out_path = "sweby.svg"
plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)