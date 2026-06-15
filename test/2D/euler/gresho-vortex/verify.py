#!/usr/bin/env python3
"""
Analyze the low-Mach Gresho vortex output of several Riemann solvers
(LMRoe, HLLC, SLAU2, HLLC-PC).

The Gresho vortex is a STEADY solution of the Euler equations (centrifugal
balance), so the exact answer is the initial vortex preserved unchanged.
Any departure of the final field from the analytic profile is pure numerical
error. Theory (Rieper 2011) predicts:

  - the low-Mach-fixed LMRoe removes the O(1/Ma) momentum dissipation
    (~ -1/2 rho a dU) and should preserve the vortex,
  - plain HLLC keeps that term (we proved HLLC == Roe at leading order),
    so it should over-diffuse: lower peak velocity, smeared u_phi(r), KE decay.

Usage:  python analyze_gresho.py
Outputs: printed metrics + gresho_analysis.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- analytic Gresho (must match build_ic.py) -------------------------------
XC, YC = 0.4, 0.4
GAMMA = 1.4

def uphi_exact(r):
    return np.where(r < 0.2, 5.0*r, np.where(r < 0.4, 2.0 - 5.0*r, 0.0))

# ---- Tecplot BLOCK reader (this specific zone layout) -----------------------
NI, NJ, NK = 81, 81, 2            # nodes
NX, NY, NZ = 80, 80, 1            # cells
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

    area = (0.8/NX)*(0.8/NY)
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

# name -> (output folder, plot color)
RUNS = [
    ("HLLC",    "OUT_HLLC",    "tab:blue"),
    ("LMRoe",   "OUT_LMRoe",   "tab:orange"),
    ("AUSM+M",  "OUT_AUSM+M",  "tab:red"),
    ("HLLC-PC", "OUT_HLLC-PC", "tab:green"),
]

def main():
    import os
    runs, colors = {}, {}
    for name, folder, c in RUNS:
        fpath = os.path.join(folder, "field.tec")
        if not os.path.exists(fpath):
            print(f"  (skipping {name}: {fpath} not found)")
            continue
        runs[name] = read_field(fpath); colors[name] = c

    res = {name: analyze(name, d) for name, d in runs.items()}
    ok = [name for name in runs if not res[name][3]]   # converged runs only

    rr = np.linspace(0, 0.45, 300)
    bins = np.linspace(0, 0.45, 31); bc = 0.5*(bins[:-1]+bins[1:])

    # ---- profile figure: scatter + binned mean/std (converged runs only) ----
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    ax[0].plot(rr, uphi_exact(rr), 'k-', lw=2, label='exact', zorder=5)
    for name in ok:
        r, uphi, ue, _ = res[name]
        ax[0].scatter(r.ravel(), uphi.ravel(), s=2, alpha=0.18, color=colors[name], label=name)
    ax[0].set_xlabel("r"); ax[0].set_ylabel(r"$u_\phi$"); ax[0].set_ylim(-0.2, 1.35)
    ax[0].legend(); ax[0].set_xlim(0,0.45)

    ax[1].plot(rr, uphi_exact(rr), 'k-', lw=2, label='exact')
    for name in ok:
        r, uphi, ue, _ = res[name]
        rf, uf = r.ravel(), uphi.ravel()
        idx = np.digitize(rf, bins)
        mean = np.array([uf[idx==k].mean() if np.any(idx==k) else np.nan for k in range(1,len(bins))])
        std  = np.array([uf[idx==k].std()  if np.any(idx==k) else np.nan for k in range(1,len(bins))])
        ax[1].plot(bc, mean, color=colors[name], label=name)
        ax[1].fill_between(bc, mean-std, mean+std, color=colors[name], alpha=0.15)
    ax[1].set_xlabel("r"); ax[1].set_ylabel(r"$u_\phi$"); ax[1].set_ylim(-0.2, 1.35)
    ax[1].legend(); ax[1].set_xlim(0,0.45)
    fig.tight_layout(); fig.savefig("gresho_profiles.png", dpi=130)
    print("\nSaved gresho_profiles.png")

    # ---- field figure: |u| for every solver (diverged panels auto-scaled) ---
    n = len(runs); ncol = min(n, 4)
    figf, axf = plt.subplots(1, ncol, figsize=(4.2*ncol, 4.2), squeeze=False)
    for axi, name in zip(axf[0], runs):
        d = runs[name]; div = res[name][3]
        vmag = np.sqrt(d["u"]**2 + d["v"]**2)
        vmax = vmag.max() if div else 1.0
        pc = axi.pcolormesh(d["xc"], d["yc"], vmag, cmap="viridis", vmin=0, vmax=vmax)
        axi.set_aspect("equal")
        axi.set_title(f"|u| — {name}" + ("  (DIVERGED)" if div else ""))
        figf.colorbar(pc, ax=axi, fraction=0.046)
    figf.tight_layout(); figf.savefig("gresho_fields.png", dpi=130)
    print("Saved gresho_fields.png")

if __name__ == "__main__":
    main()
