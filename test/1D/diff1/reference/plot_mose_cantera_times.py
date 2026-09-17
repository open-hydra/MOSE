"""Overlay MOSE (lines) vs Cantera low-Mach reference (markers) at several times.

Produces two figures -- temperature and CH4 mass fraction -- each with one MOSE
line and one set of Cantera markers per snapshot (50, 100, 150, 200 ms).

    python plot_mose_cantera_times.py            # uses OUTPUT/field{1..4}.tec
    python plot_mose_cantera_times.py a.tec b.tec ...
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os

import reference_unity_lewis as ref
from compare_mose_cantera import load_mose, solution_time

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

HERE = Path(__file__).resolve().parent
output_dir = "OUTPUT"

def main():
    if len(sys.argv) > 1:
        fields = [Path(p) for p in sys.argv[1:]]
    else:
        fields = sorted((HERE / "OUTPUT").glob("field*.tec"))
    if not fields:
        sys.exit("No field*.tec files found in OUTPUT/.")

    snaps = []
    for fp in fields:
        t = solution_time(fp)
        xm, Tm, Ym = load_mose(fp)
        xr, Tr, Yr = ref.solve(t, N=len(xm))
        snaps.append((t, xm * 1e3, Tm, Ym, Tr, Yr))
    snaps.sort(key=lambda s: s[0])

    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(snaps)))
    every = max(1, len(snaps[0][1]) // 18)   # thin Cantera markers for clarity

    plots = [
        ("Temperature [K]", 2, 4, "diffusion_T.svg"),
        (r"CH$_4$ mass fraction", 3, 5, "diffusion_CH4.svg"),
    ]

    for ylabel, im, ir, outname in plots:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        for (t, x, Tm, Ym, Tr, Yr), c in zip(snaps, colors):
            mose = Tm if im == 2 else Ym
            cant = Tr if im == 2 else Yr
            lab = f"{t*1e3:.0f} ms"
            ax.plot(x, mose, '-', color=c, lw=2.0, label=f"MOSE {lab}")
            ax.plot(x[::every], cant[::every], 'o', color=c, ms=5,
                    mfc='none', mew=1.4, label=f"Cantera {lab}")
        ax.set_xlabel("x [mm]"); ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, ncol=2, fontsize=9)
        fig.tight_layout()
        out = os.path.join(output_dir, outname)
        plt.savefig(out, bbox_inches="tight", transparent=True)
        print(f"saved {out}")


if __name__ == "__main__":
    main()
