#!/usr/bin/env python3
"""
Analyze the low-Mach Gresho vortex output of several Riemann solvers
(LMRoe, HLLC, AUSM+M).

The Gresho vortex is a STEADY solution of the Euler equations (centrifugal
balance), so the exact answer is the initial vortex preserved unchanged.
Any departure of the final field from the analytic profile is pure numerical
error. Theory (Rieper 2011) predicts:

  - the low-Mach-fixed LMRoe removes the O(1/Ma) momentum dissipation
    (~ -1/2 rho a dU) and should preserve the vortex,
  - plain HLLC keeps that term (we proved HLLC == Roe at leading order),
    so it should over-diffuse: lower peak velocity, smeared u_phi(r), KE decay.

Usage:  python analyze_gresho.py
Outputs: printed metrics + gresho_mach_matrix.png + gresho_profiles_by_mach.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ---- global plot styling -----------------------------------------------------
plt.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 25,
    "axes.titleweight": "bold",
    "axes.labelsize": 25,
    "xtick.labelsize": 23,
    "ytick.labelsize": 23,
    "legend.fontsize": 22,
    "legend.title_fontsize": 22,
    "axes.linewidth": 1.0,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

# ---- analytic Gresho (must match build_ic.py) -------------------------------
XC, YC = 0.5, 0.5
GAMMA = 1.4

def mach_field(d, Mref):
    a = np.sqrt(GAMMA * d["p"] / d["rho"])
    M = np.sqrt(d["u"]**2 + d["v"]**2) / a
    return M / Mref

def uphi_exact(r):
    return np.where(r < 0.2, 5.0*r, np.where(r < 0.4, 2.0 - 5.0*r, 0.0))

# ---- Tecplot BLOCK reader (this specific zone layout) -----------------------
NI, NJ, NK = 101, 101, 2            # nodes
NX, NY, NZ = 100, 100, 1            # cells
NNODE = NI*NJ*NK
NCELL = NX*NY*NZ

def read_field(path):
    vals = []
    with open(path) as f:
        f.readline(); f.readline()           # VARIABLES, ZONE
        for line in f:
            line = line.strip()
            if line:
                vals.append(float(line))
    v = np.array(vals)
    off = 0
    def take(n):
        nonlocal off
        a = v[off:off+n]; off += n; return a
    # nodal x,y,z  (ordering: i fastest, then j, then k)
    xn = take(NNODE).reshape(NK, NJ, NI)[0]   # (j,i) on first k-plane
    yn = take(NNODE).reshape(NK, NJ, NI)[0]
    _z = take(NNODE)
    # cell-centered scalars (ordering: i fastest, then j) on single k-plane
    def cell(): return take(NCELL).reshape(NY, NX)   # (jc,ic)
    rho = cell(); u = cell(); vv = cell(); w = cell()
    p = cell(); T = cell(); g = cell(); R = cell()
    # cell-center coordinates from the 4 surrounding nodes
    xc = 0.25*(xn[:-1,:-1] + xn[:-1,1:] + xn[1:,:-1] + xn[1:,1:])
    yc = 0.25*(yn[:-1,:-1] + yn[:-1,1:] + yn[1:,:-1] + yn[1:,1:])
    return dict(xc=xc, yc=yc, rho=rho, u=u, v=vv, p=p, T=T)

def analyze(name, d):
    xc, yc, u, v, rho, p = d["xc"], d["yc"], d["u"], d["v"], d["rho"], d["p"]
    dx = xc - XC; dy = yc - YC
    r = np.sqrt(dx**2 + dy**2)
    rsafe = np.where(r > 0, r, 1.0)
    # azimuthal velocity:  u_phi = (-u dy + v dx)/r
    uphi = (-u*dy + v*dx)/rsafe
    urad = ( u*dx + v*dy)/rsafe            # should be ~0 for a clean vortex
    ue = uphi_exact(r)

    area = (1.0/NX)*(1.0/NY)
    KE  = np.sum(0.5*rho*(u**2 + v**2))*area
    KEe = np.sum(0.5*1.0*ue**2)*area

    far = r > 0.45
    p0 = np.median(p[far])
    Mach = np.sqrt(1.0/(GAMMA*p0))

    inb = r < 0.4                          # inside the vortex
    l2 = np.sqrt(np.mean((uphi[inb]-ue[inb])**2))

    print(f"\n=== {name} ===")
    print(f"  far-field p0           : {p0:.6g}   -> inferred Mach = {Mach:.4g}")
    print(f"  peak |u_phi|           : {np.max(np.abs(uphi)):.4f}   (exact 1.0000)  "
          f"-> {100*(1-np.max(np.abs(uphi))):.1f}% decay")
    print(f"  kinetic energy         : {KE:.6e}  (exact {KEe:.6e})  "
          f"-> {100*(1-KE/KEe):.1f}% lost")
    print(f"  L2 error u_phi(r<0.4)  : {l2:.4e}")
    print(f"  spurious radial vel    : rms {np.sqrt(np.mean(urad[inb]**2)):.3e}  "
          f"max {np.max(np.abs(urad[inb])):.3e}")
    print(f"  density  min/max       : {rho.min():.6f} / {rho.max():.6f}  "
          f"(exact 1.0)")
    print(f"  pressure fluct dp/p0   : {(p.max()-p.min())/p0:.3e}  (~M^2 expected {Mach**2:.1e})")

    diverged = (np.max(np.abs(uphi)) > 2.0) or (KE/KEe > 5.0) or (rho.min() < 0.9)
    if diverged:
        print("  ** DIVERGED — excluded from profile comparison **")
    return r, uphi, ue, diverged

RUNS = [

    # solver, mach label, Mref, folder, color

    ("HLLC",  "0.1",   1e-1, "HLLC_1",  "tab:red"),
    ("HLLC",  "0.01",  1e-2, "HLLC_2",  "tab:red"),
    ("HLLC",  "0.001", 1e-3, "HLLC_3",  "tab:red"),

    ("LMRoe", "0.1",   1e-1, "LMRoe_1", "tab:blue"),
    ("LMRoe", "0.01",  1e-2, "LMRoe_2", "tab:blue"),
    ("LMRoe", "0.001", 1e-3, "LMRoe_3", "tab:blue"),

    ("AUSM+M", "0.1",   1e-1, "AUSM+M_1", "tab:green"),
    ("AUSM+M", "0.01",  1e-2, "AUSM+M_2", "tab:green"),
    ("AUSM+M", "0.001", 1e-3, "AUSM+M_3", "tab:green"),
]

# LMRoe Mach-cutoff sweep at Ma=0.001 (title, Mref, folder)
LMROE_SWEEP = [
    ("$M_{co}=0.500$", 1e-3, "LMRoe_3_high"),
    ("$M_{co}=0.005$", 1e-3, "LMRoe_3"),
    ("$M_{co}=0.001$", 1e-3, "LMRoe_3_low"),
]


def plot_lmroe_sweep():
    """Mach fields for the LMRoe dissipation variants (single row, like the
    first row of the matrix)."""
    loaded = []
    for title, Mref, folder in LMROE_SWEEP:
        fpath = os.path.join("OUTPUT", "field_"+folder+".tec")
        if not os.path.exists(fpath):
            print(f"(skipping {folder}: missing)")
            continue
        try:
            d = read_field(fpath)
            loaded.append((title, d, mach_field(d, Mref)))
        except Exception as e:
            print(f"(skipping {folder}: {e})")

    if not loaded:
        print("\nNo LMRoe sweep runs found — nothing to plot.")
        return

    fig, axs = plt.subplots(
        1, len(loaded),
        figsize=(max(4.6*len(loaded) + 0.9, 7.8), 5.2),
        squeeze=False,
        constrained_layout=True,
    )
    axs = axs[0]

    pcm = None
    for k, (ax, (title, d, Mnorm)) in enumerate(zip(axs, loaded)):
        pcm = ax.pcolormesh(
            d["xc"], d["yc"], Mnorm,
            shading="auto",
            cmap="turbo",
            vmin=0.0,
            vmax=1.1,
        )
        ax.set_aspect("equal")
        ax.set_title(title, pad=8)
        ax.set_xlabel(r"$x$")
        if k == 0:
            ax.set_ylabel(r"$y$")

    cbar = fig.colorbar(
        pcm,
        ax=axs,
        orientation="horizontal",
        fraction=0.05,
        pad=0.03,
        shrink=0.75,
        aspect=15,
    )
    cbar.set_label(r"$M/M_{\rm ref}$", fontsize=20)
    cbar.ax.tick_params(labelsize=20)

    fig.savefig("gresho_lmroe_sweep.svg", dpi=300, transparent=True)

def main():
    runs = {}
    res = {}
    colors = {}

    solvers = []
    machs = []

    for solver, mach_lbl, Mref, folder, color in RUNS:

        fpath = os.path.join("OUTPUT", "field_"+folder+".tec")
        print(fpath)

        if not os.path.exists(fpath):
            print(f"(skipping {solver} Ma={mach_lbl}: missing)")
            continue

        try:
            d = read_field(fpath)

            runs[(solver, mach_lbl)] = (d, Mref)
            res[(solver, mach_lbl)] = analyze(
                f"{solver}_{mach_lbl}", d
            )

            colors[solver] = color

            if solver not in solvers:
                solvers.append(solver)

            if mach_lbl not in machs:
                machs.append(mach_lbl)

        except Exception as e:
            print(f"(skipping {solver} Ma={mach_lbl}: {e})")

    if not solvers or not machs:
        print("\nNo runs found — nothing to plot.")
        return

    fig, axs = plt.subplots(
        len(solvers),
        len(machs),
        figsize=(max(4.6*len(machs) + 0.9, 7.8), 4.3*len(solvers) + 0.5),
        squeeze=False,
        constrained_layout=True,
    )

    pcm = None

    for i, solver in enumerate(solvers):

        for j, mach_lbl in enumerate(machs):

            ax = axs[i,j]

            key = (solver, mach_lbl)

            if key not in runs:

                ax.text(
                    0.5, 0.5,
                    "MISSING",
                    ha="center",
                    va="center",
                    transform=ax.transAxes
                )
                ax.set_axis_off()
                continue

            d, Mref = runs[key]

            diverged = res[key][3]

            Mnorm = mach_field(d, Mref)

            pcm = ax.pcolormesh(
                d["xc"],
                d["yc"],
                Mnorm,
                shading="auto",
                cmap="turbo",
                vmin=0.0,
                vmax=1.1
            )

            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])

            if i == 0:
                ax.set_title(rf"$M_{{ref}}={mach_lbl}$", pad=8, fontsize=28)

            if j == 0:
                ax.set_ylabel(solver, fontsize=28, fontweight="bold", labelpad=8)

            if diverged:

                ax.text(
                    0.5, 0.95,
                    "DIVERGED",
                    color="white",
                    ha="center",
                    va="top",
                    fontsize=11,
                    fontweight="bold",
                    transform=ax.transAxes,
                    bbox=dict(facecolor="red", alpha=0.85, boxstyle="round,pad=0.3")
                )

    if pcm is not None:
        cbar = fig.colorbar(
            pcm,
            ax=axs,
            fraction=0.025,
            pad=0.02,
            shrink=0.65,
            aspect=20,
        )
        cbar.set_label(r"$M/M_{\rm ref}$", fontsize=20)
        cbar.ax.tick_params(labelsize=20)

    fig.savefig("gresho_mach_matrix.svg", dpi=100, transparent=True)

    # genuine horizontal slice through the vortex centre (y = YC),
    # normalised by the half-domain so r* spans [-1, 1].
    RHALF = 0.5
    ss = np.linspace(-1.0, 1.0, 600)
    vexact = np.sign(ss)*uphi_exact(RHALF*np.abs(ss))

    markers = {s: m for s, m in zip(solvers, ["o", "s", "^", "D", "v"])}

    fig, axs = plt.subplots(
        1,
        len(machs),
        figsize=(max(5.6*len(machs), 7.8), 5.2),
        sharey=True,
        constrained_layout=True,
    )

    if len(machs) == 1:
        axs = [axs]

    for ax, mach_lbl in zip(axs, machs):

        ax.axhline(0.0, color="0.6", lw=0.8)

        ax.plot(
            ss,
            vexact,
            'k-',
            lw=4,
            label='Exact'
        )

        for solver in solvers:

            key = (solver, mach_lbl)

            if key not in runs or res[key][3]:   # missing or diverged
                continue

            d, _ = runs[key]
            yc = d["yc"]
            jrow = np.argmin(np.abs(yc[:, yc.shape[1]//2] - YC))

            x = d["xc"][jrow, :]
            # on the horizontal centre line the tangential velocity is v,
            # signed so it flips across r = 0 (genuine vortex slice).
            vline = d["v"][jrow, :]
            rstar = (x - XC)/RHALF

            ax.plot(
                rstar,
                vline,
                linestyle="none",
                marker=markers[solver],
                markersize=12,
                markerfacecolor="none",
                markeredgewidth=3,
                markevery=4,
                color=colors[solver],
                label=solver
            )

        ax.set_title(rf"$M_{{ref}}={mach_lbl}$", pad=8)
        ax.set_xlim(-1.0,1.0)
        ax.set_ylim(-1.25,1.25)
        ax.set_xlabel(r"$r/R$")
        ax.xaxis.set_major_locator(mticker.MultipleLocator(0.5))
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
        ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.6)

    axs[0].set_ylabel(r"$u_\phi$")

    handles, labels = axs[-1].get_legend_handles_labels()
    leg = fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.14),
        ncol=len(labels),
        frameon=False,
    )

    fig.savefig("gresho_profiles_by_mach.svg", dpi=300, bbox_inches="tight", bbox_extra_artists=(leg,), transparent=True)

    plot_lmroe_sweep()

if __name__ == "__main__":
    main()