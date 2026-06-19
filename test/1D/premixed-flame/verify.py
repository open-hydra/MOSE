"""Compare the MOSE 1D premixed-flame solution against a Cantera reference.

Plots velocity, temperature and OH mass fraction. Profiles are aligned on the
flame front (T = T_REF crossing) since the two solvers anchor the flame at
different absolute x-locations and use different domain lengths.

Field reading recycled from verify.py (ORION read_TEC). NB: this case carries
25 species, so the variable layout in field.tec is:
    var[0]   = rho(1) ... var[24] = rho(25)   (partial densities)
    var[25]  = u,  [28] = p,  [29] = T
OH is the 7th species in INPUT/phase.txt -> rho(7) -> var index 6.
"""
import sys
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")

mpl.rcParams.update({"font.size": 13, "axes.labelsize": 14, "legend.fontsize": 12})

T_REF = 1500.0          # K, iso-temperature used to align the flame fronts
OH_SPECIES_INDEX = 7    # 1-based position of OH in INPUT/phase.txt


def _setup_orion_import() -> None:
    root = Path(__file__).resolve().parent
    for parent in [root, *root.parents]:
        candidate = parent / "lib" / "ORION" / "src" / "python"
        if candidate.exists():
            sys.path.insert(0, str(candidate))
            return


_setup_orion_import()
from ORION import read_TEC  # noqa: E402


def load_mose(field_path: Path):
    """Return x_centers, u, T, Y_OH from a MOSE field.tec file."""
    x_arr, _y, _z, var, _names = read_TEC(str(field_path))
    x_nodes = np.asarray(x_arr[0])[:, 0, 0]
    xc = 0.5 * (x_nodes[:-1] + x_nodes[1:])

    rho = [np.asarray(var[0][i]).reshape(-1) for i in range(25)]
    rho_tot = np.sum(rho, axis=0)
    u = np.asarray(var[0][25]).reshape(-1)
    T = np.asarray(var[0][29]).reshape(-1)
    Y_OH = rho[OH_SPECIES_INDEX - 1] / rho_tot

    n = min(xc.size, u.size, T.size, Y_OH.size)
    return xc[:n], u[:n], T[:n], Y_OH[:n]


def load_cantera(csv_path: Path):
    c = np.genfromtxt(csv_path, delimiter=",", names=True)
    return c["x"], c["u"], c["T"], c["Y_OH"]


def align_front(x, T, T_ref=T_REF):
    """Shift x so that the T = T_ref crossing sits at x = 0 (T increasing)."""
    order = np.argsort(T)
    x0 = np.interp(T_ref, T[order], x[order])
    return x - x0


def main() -> None:
    root = Path(__file__).resolve().parent
    x_m, u_m, T_m, oh_m = load_mose(root / "OUTPUT" / "field.tec")
    x_c, u_c, T_c, oh_c = load_cantera(root / "flame_phi_1.0.csv")

    # align on the flame front, convert to millimetres
    xs_m = align_front(x_m, T_m) * 1e3
    xs_c = align_front(x_c, T_c) * 1e3

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    panels = [
        (axs[0], T_m, T_c, "Temperature [K]"),
        (axs[1], u_m, u_c, "Velocity [m/s]"),
        (axs[2], oh_m * 1e3, oh_c * 1e3, r"OH mass fraction [$\times 10^{-3}$]"),
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

    fig.suptitle(r"1D premixed CH$_4$/air laminar flame ($\phi=1.0$): MOSE vs Cantera",
                 fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out = root / "mose_vs_cantera.png"
    fig.savefig(out, dpi=150)
    print(f"Saved figure to: {out}")

    # quick quantitative summary
    print("\nGlobal comparison:")
    print(f"  S_L (unburnt u):  MOSE {u_m[np.argmin(T_m)]:.4f}  Cantera {u_c[np.argmin(T_c)]:.4f} m/s")
    print(f"  T_burnt:          MOSE {T_m.max():.1f}    Cantera {T_c.max():.1f} K")
    print(f"  u_burnt:          MOSE {u_m.max():.4f}   Cantera {u_c.max():.4f} m/s")
    print(f"  Y_OH peak:        MOSE {oh_m.max():.4e}  Cantera {oh_c.max():.4e}")


if __name__ == "__main__":
    main()
