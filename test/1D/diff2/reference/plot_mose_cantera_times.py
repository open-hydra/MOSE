"""Overlay MOSE (lines) vs the Cantera multicomponent reference (markers) at several times.

Produces three figures -- temperature, CH4 and N2 mass fraction -- each with one MOSE
line and one set of Cantera markers per snapshot (50, 100, 150, 200 ms). The reference
is the low-Mach mixture-averaged solver in reference_multicomponent.py using MOSE's own
(1 - X_k) diffusion rule, so the markers measure pure numerical error.

    python plot_mose_cantera_times.py            # uses ../OUTPUT/field{1..4}.tec
    python plot_mose_cantera_times.py a.tec b.tec ...

SVGs are written to docs/vv/images/ (multicomp_T.svg, multicomp_CH4.svg, multicomp_N2.svg).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import reference_multicomponent as ref
from compare_mose_cantera import load_mose_full, solution_time

# Transparent, theme-aware SVGs for the docs
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

HERE = Path(__file__).resolve().parent
IMG_DIR = HERE.parents[3] / "docs" / "vv" / "images"   # .../MOSE/docs/vv/images

# species index in reference SP order [CH4, H2O, N2, O2]
I_CH4, I_N2 = 0, 2


def main():
    if len(sys.argv) > 1:
        fields = [Path(p) for p in sys.argv[1:]]
    else:
        fields = sorted((HERE.parent / "OUTPUT").glob("field*.tec"))
    if not fields:
        sys.exit("No field*.tec files found in OUTPUT/.")
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    snaps = []
    for fp in fields:
        t = solution_time(fp)
        xm, Tm, Ym = load_mose_full(fp)             # Ym: (n, nsp)
        xr, Tr, Yr = ref.solve(t, N=len(xm), return_full=True)
        snaps.append((t, xm * 1e3, Tm, Ym, Tr, Yr))
    snaps.sort(key=lambda s: s[0])

    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(snaps)))
    every = max(1, len(snaps[0][1]) // 18)          # thin Cantera markers for clarity

    # (ylabel, MOSE getter, reference getter, output file)
    plots = [
        ("Temperature [K]",       lambda Tm, Ym: Tm,           lambda Tr, Yr: Tr,           "multicomp_T.svg"),
        (r"CH$_4$ mass fraction", lambda Tm, Ym: Ym[:, I_CH4], lambda Tr, Yr: Yr[:, I_CH4], "multicomp_CH4.svg"),
        (r"N$_2$ mass fraction",  lambda Tm, Ym: Ym[:, I_N2],  lambda Tr, Yr: Yr[:, I_N2],  "multicomp_N2.svg"),
    ]

    for ylabel, mose_of, ref_of, outname in plots:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        for (t, x, Tm, Ym, Tr, Yr), c in zip(snaps, colors):
            lab = f"{t*1e3:.0f} ms"
            ax.plot(x, mose_of(Tm, Ym), '-', color=c, lw=2.0, label=f"MOSE {lab}")
            ax.plot(x[::every], ref_of(Tr, Yr)[::every], 'o', color=c, ms=5,
                    mfc='none', mew=1.4, label=f"Cantera {lab}")
        ax.set_xlabel("x [mm]"); ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, ncol=2, fontsize=9)
        fig.tight_layout()
        out = IMG_DIR / outname
        plt.savefig(out, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print(f"saved {out}")


if __name__ == "__main__":
    main()
