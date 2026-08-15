#!/usr/bin/env python3
"""Cf comparison: MOSE constant-mu vs Sutherland viscosity, one figure per model.

Overlays the two MOSE runs with the experiment / Wind-US references and
prints the separation-bubble length of each MOSE curve. Theme-aware SVG output
in the same style as compare.py (recolours with the docs light/dark theme).

Usage:
    python3 plot_sutherland.py                 # both SA and SST
    python3 plot_sutherland.py --model SST     # one model
    python3 plot_sutherland.py --suth OUTPUT/wall.tec   # custom Sutherland wall file

Inputs (per <MODEL>/):
    reference/MOSE-<MODEL>.tec   constant-mu MOSE wall  (baseline)
    <suth>                       Sutherland  MOSE wall  (default OUTPUT/wall.tec)
    reference/{wind}-<MODEL>.dat, reference/schulein.dat   references
Output:  SWBLI-cf-<model>-sutherland.svg  (in this folder)
"""
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from compare import mose_wall, sep_reatt, style, MOSE_C   # noqa: E402

CONST_C = "0.45"                       # grey: old constant-mu curve
REF_STYLE = {"Wind-US": ("b", "-.")}
REFS = {"SA":  [("Wind-US", "wind-SA.dat")],
        "SST": [("Wind-US", "wind-SST.dat")]}


def _xy(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


def _bubble(x, cf):
    s, r = sep_reatt(np.asarray(x), np.asarray(cf))
    return f"{1e3 * (r - s):.1f} mm" if s else "attached"


def figure(model, suth_rel):
    m = model.lower()
    case = ROOT / model
    xc, cfc, _ = mose_wall(case / "reference" / f"MOSE-{model}.tec")   # constant mu
    xs, cfs, _ = mose_wall(case / suth_rel)                            # Sutherland

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    exp = case / "reference" / "schulein.dat"
    if exp.exists():
        xe, ce = _xy(exp)
        ax.plot(xe, 1e3 * ce, "ko", ms=4, mfc="none", label="experiment (Schülein)")
    for name, fn in REFS[model]:
        f = case / "reference" / fn
        if f.exists():
            xr, cr = _xy(f)
            c, ls = REF_STYLE[name]
            ax.plot(xr, 1e3 * cr, c, ls=ls, lw=1.6, alpha=0.9, label=name)

    ax.plot(xc, 1e3 * cfc, color=CONST_C, ls="--", lw=2.6,
            label=fr"MOSE, constant $\mu$  (bubble {_bubble(xc, cfc)})")
    ax.plot(xs, 1e3 * cfs, MOSE_C, lw=4.0,
            label=f"MOSE, Sutherland  (bubble {_bubble(xs, cfs)})")

    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlim(0.32, 0.41)
    ax.set_ylim(-0.002*1000, 0.0075*1000)
    style(ax, "x [m]", r"$C_f \times 10^{3}$")
    ax.set_title(fr"SWBLI {model}: effect of Sutherland viscosity on $C_f$", fontsize=14)

    out = ROOT / f"SWBLI-cf-{m}-sutherland.svg"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}   "
          f"(const-mu bubble {_bubble(xc, cfc)}, Sutherland {_bubble(xs, cfs)})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=("SA", "SST"), help="default: both")
    ap.add_argument("--suth", default="OUTPUT/wall.tec",
                    help="Sutherland wall file, relative to <MODEL>/ "
                         "(default: OUTPUT/wall.tec)")
    args = ap.parse_args()
    for model in ([args.model] if args.model else ["SA", "SST"]):
        figure(model, args.suth)


if __name__ == "__main__":
    main()
