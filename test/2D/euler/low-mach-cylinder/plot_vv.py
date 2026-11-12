"""Figures of the V&V page of the low-Mach cylinder (docs/vv/2D-low-mach-cylinder.md), from the runs of run.sh.

    ./run.sh && python3 -B plot_vv.py          writes docs/vv/images/LMC-*.svg and prints the numbers of the page

Runs read (runs/<tag>/field.tec, runs/<tag>/residual-history.dat): prec-m2 (prec, HLLC-PC, Mach 0.01), prec-m3 (the
same at Mach 0.001) and cons-m2 (cons, HLLC, Mach 0.01).
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import LogNorm

import make_case as mc

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "docs", "vv", "images")
RUNS = {"prec-m2": 0.01, "prec-m3": 0.001, "cons-m2": 0.01}
LABEL = {"prec-m2": "prec + HLLC-PC, M = 0.01", "prec-m3": "prec + HLLC-PC, M = 0.001",
         "cons-m2": "cons + HLLC, M = 0.01"}
COLOR = {"prec-m2": "#1f77b4", "prec-m3": "#ff7f0e", "cons-m2": "#d62728"}
R = mc.R

plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12, "axes.linewidth": 1.0,
    "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none", "svg.fonttype": "none",
})


def load(tag):
    """Nodes (ni, nj, x, y) and cell fields (rho, u, v, p) of runs/<tag>/field.tec, cells as (nj-1, ni-1) arrays."""
    with open(os.path.join("runs", tag, "field.tec")) as f:
        names = f.readline().split("=", 1)[1].replace('"', " ").split()
        zone = f.readline()
        vals = np.array(f.read().split("DATASETAUXDATA")[0].split(), dtype=float)
    ni, nj, nk = (int(zone.split(" %s=" % d)[1].split(",")[0]) for d in "IJK")
    nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * (nk - 1)
    x = vals[:nn].reshape(nk, nj, ni)[0]
    y = vals[nn:2 * nn].reshape(nk, nj, ni)[0]
    cell = lambda n: vals[3 * nn + (names.index(n) - 3) * nc:3 * nn + (names.index(n) - 2) * nc].reshape(nj - 1, ni - 1)
    return dict(x=x, y=y, rho=cell("rho(1)"), u=cell("u"), v=cell("v"), p=cell("p"))


def centres(d):
    x, y = d["x"], d["y"]
    return 0.25 * (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:]), 0.25 * (y[:-1, :-1] + y[:-1, 1:] + y[1:, :-1] + y[1:, 1:])


def q_inf(mach):
    g = mc.GAMMA
    return mc.P_INF * ((1.0 + 0.5 * (g - 1.0) * mach ** 2) ** (g / (g - 1.0)) - 1.0)


def cp(d, mach):
    return (d["p"] - mc.P_INF) / q_inf(mach)


def cp_potential(cx, cy):
    r, th = np.hypot(cx, cy), np.arctan2(cy, cx)
    return 2.0 * (R / r) ** 2 * np.cos(2.0 * th) - (R / r) ** 4


def field_axes(ax, lim=3.0):
    ax.set_rasterization_zorder(0.5)        # the filled contours (zorder 0) as an image inside the SVG
    ax.set_aspect("equal")
    ax.set_xlim(-lim * R, lim * R)
    ax.set_ylim(-lim * R, lim * R)
    ax.add_patch(plt.Circle((0, 0), R, facecolor="0.6", edgecolor="0.3", lw=0.8, zorder=5))
    ax.set_xlabel("x [m]")


def tri(d):
    cx, cy = centres(d)
    # close the O-grid: repeat the first column of cells after the last
    cx = np.hstack([cx, cx[:, :1]])
    cy = np.hstack([cy, cy[:, :1]])
    nj, ni = cx.shape
    tris = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            a, b, c, e = j * ni + i, j * ni + i + 1, (j + 1) * ni + i, (j + 1) * ni + i + 1
            tris += [(a, b, e), (a, e, c)]
    return mtri.Triangulation(cx.ravel(), cy.ravel(), tris)


def closed(a):
    return np.hstack([a, a[:, :1]]).ravel()


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", os.path.relpath(path))


def fig_grid(d):
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))
    for a, lim, step in ((ax[0], 40.0, 2), (ax[1], 2.0, 1)):
        x, y = d["x"], d["y"]
        for j in range(0, x.shape[0], step):
            a.plot(x[j], y[j], color="0.45", lw=0.4)
        for i in range(0, x.shape[1], step):
            a.plot(x[:, i], y[:, i], color="0.45", lw=0.4)
        a.set_aspect("equal")
        a.set_xlim(-lim * R, lim * R)
        a.set_ylim(-lim * R, lim * R)
        a.set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_title("whole domain (every second line)")
    ax[1].set_title("near the cylinder")
    save(fig, "LMC-grid.svg")


def fig_cp_fields(data):
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    for a, tag in zip(ax, ("prec-m2", "prec-m3", "cons-m2")):
        d = data[tag]
        c = cp(d, RUNS[tag])
        lim = (-3.0, 1.0) if tag != "cons-m2" else (-8.0, 8.0)
        cs = a.tricontourf(tri(d), closed(c), levels=np.linspace(*lim, 41), cmap="RdBu_r", extend="both", zorder=0)
        field_axes(a)
        a.set_title(LABEL[tag])
        fig.colorbar(cs, ax=a, shrink=0.8, label="Cp", ticks=np.linspace(*lim, 5))
    ax[0].set_ylabel("y [m]")
    save(fig, "LMC-cp-fields.svg")


def wall_row(d, mach):
    cx, cy = centres(d)
    th = np.degrees(np.arctan2(cy[0], cx[0]))
    o = np.argsort(th)
    return th[o], cp(d, mach)[0][o], cp_potential(cx[0], cy[0])[o]


def fig_cp_wall(data):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    th, _, pot = wall_row(data["prec-m2"], 0.01)
    for a in ax:
        a.plot(th, pot, "-", color="0.2", lw=1.6, label="potential flow (cell centres)")
    for tag, mk in (("prec-m2", "o"), ("prec-m3", "x")):
        t, c, _ = wall_row(data[tag], RUNS[tag])
        ax[0].plot(t, c, mk, ms=5, mfc="none", color=COLOR[tag], label=LABEL[tag])
    t, c, _ = wall_row(data["cons-m2"], 0.01)
    ax[1].plot(t, c, "s", ms=4, mfc="none", color=COLOR["cons-m2"], label=LABEL["cons-m2"])
    t, c, _ = wall_row(data["prec-m2"], 0.01)
    ax[1].plot(t, c, "o", ms=3, color=COLOR["prec-m2"], label=LABEL["prec-m2"])
    for a in ax:
        a.set_xlabel(r"$\theta$ [deg] (0 = downstream)")
        a.set_xlim(-180, 180)
        a.set_xticks(range(-180, 181, 60))
        a.grid(True, alpha=0.3, ls="--")
    for a in ax:
        a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False)
    ax[0].set_ylabel("Cp, first row of cells")
    ax[0].set_ylim(-3.5, 1.5)
    ax[0].set_title("preconditioned: Mach 0.01 and 0.001")
    ax[1].set_title("unpreconditioned against preconditioned")
    save(fig, "LMC-cp-wall.svg")


def fig_mach_independence(data):
    a2, a3 = data["prec-m2"], data["prec-m3"]
    diff = np.abs(cp(a2, 0.01) - cp(a3, 0.001))
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    cs = ax.tricontourf(tri(a2), closed(np.maximum(diff, 1e-7)), levels=np.logspace(-7, -3, 17), norm=LogNorm(),
                        cmap="viridis", extend="both", zorder=0)
    field_axes(ax)
    ax.set_ylabel("y [m]")
    ax.set_title("|Cp(M = 0.01) - Cp(M = 0.001)|")
    fig.colorbar(cs, ax=ax, shrink=0.8, ticks=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3])
    save(fig, "LMC-mach-independence.svg")


def fig_total_pressure(data):
    d = data["prec-m2"]
    p0 = mc.P_INF + q_inf(0.01)
    pt = d["p"] + 0.5 * d["rho"] * (d["u"] ** 2 + d["v"] ** 2)
    loss = (p0 - pt) / q_inf(0.01)
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    cs = ax.tricontourf(tri(d), closed(loss), levels=np.linspace(0, 1, 21), cmap="magma_r", extend="both", zorder=0)
    field_axes(ax, 3.0)
    ax.set_xlim(-2 * R, 6 * R)
    ax.set_ylim(-2 * R, 2 * R)
    ax.set_ylabel("y [m]")
    ax.set_title("total-pressure loss (p0 - pt)/q, prec + HLLC-PC, M = 0.01")
    fig.colorbar(cs, ax=ax, shrink=0.8)
    save(fig, "LMC-total-pressure.svg")


def fig_residuals():
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for tag in RUNS:
        h = np.loadtxt(os.path.join("runs", tag, "residual-history.dat"))
        ax.semilogy(h[::10, 0], h[::10, 2] / h[0, 2], color=COLOR[tag], lw=1.4, ls="--" if tag == "prec-m3" else "-",
                    label=LABEL[tag])
    ax.set_xlabel("iteration")
    ax.set_ylabel("density residual / first iteration")
    ax.grid(True, which="both", alpha=0.3, ls="--")
    ax.legend()
    save(fig, "LMC-residuals.svg")


def numbers(data):
    a2, a3, c2 = data["prec-m2"], data["prec-m3"], data["cons-m2"]
    cx, cy = centres(a2)
    near = np.hypot(cx, cy) < 1.5 * R
    print("  |Cp(0.01) - Cp(0.001)| max over r < 1.5 R: %.2e" % np.abs(cp(a2, 0.01) - cp(a3, 0.001))[near].max())
    front = cx[0] < 0
    for tag, d in (("prec-m2", a2), ("prec-m3", a3), ("cons-m2", c2)):
        e = np.abs(cp(d, RUNS[tag]) - cp_potential(cx, cy))
        h = np.loadtxt(os.path.join("runs", tag, "residual-history.dat"))
        print("  %s: front-half first row max |Cp - potential| %.3f, r < 1.5 R max %.3f, residual fall %.2e"
              % (tag, e[0][front].max(), e[near].max(), h[-1, 2] / h[0, 2]))
    print("  |Cp(cons) - Cp(prec)| max over r < 1.5 R: %.2f" % np.abs(cp(c2, 0.01) - cp(a2, 0.01))[near].max())
    p0 = mc.P_INF + q_inf(0.01)
    pt = a2["p"] + 0.5 * a2["rho"] * (a2["u"] ** 2 + a2["v"] ** 2)
    print("  largest total-pressure loss / q (prec-m2): %.2f" % ((p0 - pt) / q_inf(0.01)).max())


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {tag: load(tag) for tag in RUNS}
    fig_grid(data["prec-m2"])
    fig_cp_fields(data)
    fig_cp_wall(data)
    fig_mach_independence(data)
    fig_total_pressure(data)
    fig_residuals()
    numbers(data)


if __name__ == "__main__":
    main()
