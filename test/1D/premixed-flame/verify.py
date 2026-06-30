"""Verify the MOSE 1D premixed CH4/air laminar flame against a Cantera reference.

Merges the former ``verify.py`` (spatial profiles aligned on the flame front) and
``plot_progress.py`` (journal-style comparison against the progress variable
``c = (T - T_u)/(T_b - T_u)``). Produces three theme-aware, transparent SVG
figures. The two documented figures go to ``docs/vv/images/``:

    flame-thermo.svg   : density (left) + temperature (right) vs c
    flame-species.svg  : major species (left) + minor species (right) vs c

and the spatial-profile figure is kept as a local diagnostic in ``OUTPUT/``:

    flame-profiles.svg : T, u, OH vs x aligned on the T = 1500 K crossing

Field reading uses ORION ``read_TEC``. This case carries 25 species, so the
field.tec variable layout is::

    var[0..24] = rho(1)..rho(25)   (partial densities, INPUT/phase.txt order)
    var[25]=u   var[28]=p   var[29]=T

Density rho = sum_k rho(k);  species mass fraction Y_i = rho(i)/rho. The 1-based
species positions follow INPUT/phase.txt (OH is 7th -> var index 6).
"""
import argparse
import sys
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# Theme-aware transparent SVG styling (matches the rest of the V&V figures)
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    "font.size": 13,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "lines.linewidth": 2.2,
})

T_REF = 1500.0          # K, iso-temperature used to align the flame fronts
N_SPECIES = 25
U_VAR, T_VAR = 25, 29
# 1-based species position in INPUT/phase.txt -> field var index is (pos - 1)
SPECIES_IDX = {"CH4": 17, "O2": 6, "OH": 7, "H2O": 8, "CO": 11, "CO2": 12}
MS, MEW = 7, 1.4        # Cantera marker size / edge width


def _setup_orion_import() -> None:
    root = Path(__file__).resolve().parent
    for parent in [root, *root.parents]:
        candidate = parent / "lib" / "ORION" / "src" / "python"
        if candidate.exists():
            sys.path.insert(0, str(candidate))
            return


_setup_orion_import()
try:
    from ORION import read_TEC  # noqa: E402
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or install ORION "
        "in the current Python environment."
    )


def progress_variable(T: np.ndarray) -> np.ndarray:
    Tu, Tb = T.min(), T.max()
    return (T - Tu) / (Tb - Tu)


def load_mose(field_path: Path) -> dict:
    """Return x-centres plus profiles and progress variable from a field.tec."""
    x_arr, _y, _z, var, _names = read_TEC(str(field_path))
    x_nodes = np.asarray(x_arr[0])[:, 0, 0]
    xc = 0.5 * (x_nodes[:-1] + x_nodes[1:])

    rho = [np.asarray(var[0][i]).reshape(-1) for i in range(N_SPECIES)]
    rho_tot = np.sum(rho, axis=0)
    u = np.asarray(var[0][U_VAR]).reshape(-1)
    T = np.asarray(var[0][T_VAR]).reshape(-1)

    n = min(xc.size, u.size, T.size, rho_tot.size)
    out = {"x": xc[:n], "u": u[:n], "T": T[:n], "rho": rho_tot[:n],
           "c": progress_variable(T[:n])}
    for name, pos in SPECIES_IDX.items():
        out[f"Y_{name}"] = (rho[pos - 1] / rho_tot)[:n]
    return out


def load_cantera(csv_path: Path) -> dict:
    d = np.genfromtxt(csv_path, delimiter=",", names=True)
    out = {"x": d["x"], "u": d["u"], "T": d["T"], "c": progress_variable(d["T"])}
    if "rho" in d.dtype.names:
        out["rho"] = d["rho"]
    for name in SPECIES_IDX:
        col = f"Y_{name}"
        if col in d.dtype.names:
            out[col] = d[col]
    return out


def align_front(x, T, t_ref=T_REF):
    """Shift x so that the T = t_ref crossing sits at x = 0 (T increasing)."""
    order = np.argsort(T)
    x0 = np.interp(t_ref, np.asarray(T)[order], np.asarray(x)[order])
    return np.asarray(x) - x0


def save(fig, output_dir: Path, name: str, show: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{name}.svg"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    print(f"Saved {out}")
    if not show:
        plt.close(fig)


def figure_profiles(m, c, output_dir, show):
    """T, u, OH vs flame-aligned x (former verify.py)."""
    xs_m = align_front(m["x"], m["T"]) * 1e3
    xs_c = align_front(c["x"], c["T"]) * 1e3

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    panels = [
        (axs[0], m["T"], c["T"], "Temperature [K]"),
        (axs[1], m["u"], c["u"], "Velocity [m/s]"),
        (axs[2], m["Y_OH"] * 1e3, c["Y_OH"] * 1e3,
         r"OH mass fraction [$\times 10^{-3}$]"),
    ]
    for ax, ym, yc, ylabel in panels:
        ax.plot(xs_m, ym, "-", color="C3", lw=2.2, label="MOSE")
        ax.plot(xs_c, yc, "--o", color="k", lw=1.6, ms=3, mfc="none",
                markevery=6, label="Cantera")
        ax.set_xlabel("x − x(T=1500 K) [mm]")
        ax.set_ylabel(ylabel)
        ax.set_xlim(-50, 50)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, output_dir, "flame-profiles", show)


def figure_thermo(m, c, mev, output_dir, show):
    """Density (left axis) + temperature (right axis) vs c (former progress_a)."""
    fig, axl = plt.subplots(figsize=(7.5, 6.0))
    axr = axl.twinx()
    cden, cT = "#1f5fbf", "#c0392b"

    axl.plot(m["c"], m["rho"], "-", color=cden)
    axr.plot(m["c"], m["T"], "-", color=cT)
    if "rho" in c:
        axl.plot(c["c"], c["rho"], "s", color=cden, ms=MS, mfc="none",
                 mew=MEW, markevery=mev)
    axr.plot(c["c"], c["T"], "o", color=cT, ms=MS, mfc="none", mew=MEW,
             markevery=mev)

    axl.set_xlabel(r"$c\;(-)$")
    axl.set_ylabel(r"Density $\rho\;\mathrm{[kg\,m^{-3}]}$", color=cden)
    axr.set_ylabel(r"Temperature $T\;\mathrm{[K]}$", color=cT)
    axl.tick_params(axis="y", colors=cden)
    axr.tick_params(axis="y", colors=cT)
    axl.set_xlim(0, 1)

    handles = [
        Line2D([0], [0], color="k", lw=2.2, label="MOSE"),
        Line2D([0], [0], color="k", marker="o", mfc="none", mew=MEW, ms=MS,
               ls="none", label="Cantera"),
        Line2D([0], [0], color=cden, lw=2.2, label=r"$\rho$"),
        Line2D([0], [0], color=cT, lw=2.2, label=r"$T$"),
    ]
    axl.legend(handles=handles, loc="center right", frameon=False)
    fig.tight_layout()
    save(fig, output_dir, "flame-thermo", show)


def figure_species(m, c, mev, output_dir, show):
    """Major species (left) + minor species (right) vs c (former progress_b)."""
    fig, axl = plt.subplots(figsize=(7.5, 6.0))
    axr = axl.twinx()

    major = [("Y_CH4", "#1f77b4", "s", "CH$_4$"),
             ("Y_O2",  "#2ca02c", "o", "O$_2$"),
             ("Y_CO2", "#9467bd", "^", "CO$_2$"),
             ("Y_H2O", "#17becf", "D", "H$_2$O")]
    minor = [("Y_OH",  "#d62728", "v", "OH"),
             ("Y_CO",  "#ff7f0e", "P", "CO")]

    leg = []
    for key, col, mk, lab in major:
        axl.plot(m["c"], m[key], "-", color=col)
        if key in c:
            axl.plot(c["c"], c[key], mk, color=col, ms=MS, mfc="none",
                     mew=MEW, markevery=mev)
        leg.append(Line2D([0], [0], color=col, lw=2.2, marker=mk, mfc="none",
                          mew=MEW, ms=MS, label=lab))
    for key, col, mk, lab in minor:
        axr.plot(m["c"], m[key], "--", color=col)
        if key in c:
            axr.plot(c["c"], c[key], mk, color=col, ms=MS, mfc="none",
                     mew=MEW, markevery=mev)
        leg.append(Line2D([0], [0], color=col, lw=2.2, ls="--", marker=mk,
                          mfc="none", mew=MEW, ms=MS, label=lab))

    axl.set_xlabel(r"$c\;(-)$")
    axl.set_ylabel(r"Major species mass fraction")
    axr.set_ylabel(r"Minor species mass fraction")
    axl.set_xlim(0, 1)
    axl.set_ylim(bottom=0)
    axr.set_ylim(bottom=0)

    style = [Line2D([0], [0], color="k", lw=2.2, label="MOSE"),
             Line2D([0], [0], color="k", marker="o", mfc="none", mew=MEW,
                    ms=MS, ls="none", label="Cantera")]
    axl.legend(handles=leg + style, loc="upper center", frameon=False,
               ncol=2, fontsize=11)
    fig.tight_layout()
    save(fig, output_dir, "flame-species", show)


def summary(m, c) -> None:
    print("\nGlobal comparison (MOSE / Cantera):")
    print(f"  S_L  (unburnt u):  {m['u'][np.argmin(m['T'])]:.4f} / "
          f"{c['u'][np.argmin(c['T'])]:.4f} m/s")
    print(f"  T_burnt:           {m['T'].max():.1f} / {c['T'].max():.1f} K")
    print(f"  u_burnt:           {m['u'].max():.4f} / {c['u'].max():.4f} m/s")
    print(f"  Y_OH peak:         {m['Y_OH'].max():.4e} / {c['Y_OH'].max():.4e}")


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", default="OUTPUT/field.tec")
    parser.add_argument("--cantera", default="reference/flame_phi_1.0.csv")
    parser.add_argument("--output-dir", default="../../../docs/vv/images",
                        help="directory for the SVG figures (relative to the test dir)")
    parser.add_argument("--plot", action="store_true",
                        help="also display the figures interactively")
    args = parser.parse_args()

    m = load_mose(root / args.field)
    c = load_cantera(root / args.cantera)
    mev = max(1, len(c["c"]) // 22)
    output_dir = (root / args.output_dir).resolve()

    # Documented figures go to docs/vv/images; the spatial-profile figure is a
    # local diagnostic (velocity, OH in physical space) kept in OUTPUT/.
    figure_thermo(m, c, mev, output_dir, args.plot)
    figure_species(m, c, mev, output_dir, args.plot)
    figure_profiles(m, c, (root / "OUTPUT").resolve(), args.plot)
    summary(m, c)

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
