#!/usr/bin/env python3
"""Figure of the test mr-moving-levels (simple wave with moving level boundaries) for its V&V page (docs/vv/multirate.md), written from
the files <label>.json of verify_moving_levels.py in RUNS_DIR:
  multirate-moving-levels-order.svg  left, the L1 error of the density against the number of cells n for the multirate with the levels
                          recomputed at every macro step (mr), with the map of the first macro step kept (mrf), the single
                          rate at the same macro step (srH) and at dt_min (srdt), with the second-order slope; right, the
                          Linf error of the density of the same runs, with the first- and the second-order slopes; the
                          observed orders of the last pair in the legend.
The test has one form, run in CI; the legend says so. Style of the V&V figures of MOSE (mr-common/mrplot.py). Without
matplotlib the figure is skipped; the test does not depend on it.
Usage: plot_moving_levels.py RUNS_DIR [OUT_DIR]   (OUT_DIR defaults to RUNS_DIR)"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrplot as mp       # noqa: E402

RUNS = [  # label, legend, colour, line, marker
    ('mr', 'multirate, levels recomputed at every macro step', mp.ORANGE, '-', 'o'),
    ('mrf', 'multirate, map of the first macro step kept', mp.GREEN, '-', 's'),
    ('srH', 'single rate, the same macro step H', mp.BLUE, '-', '^'),
    ('srdt', 'single rate, step dt_min', '0.35', '--', 'D'),
]
FIGSIZE2 = (10, 5)
FORM = 'the form run in CI (the test has no other)'


def main():
    runs = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else runs
    plt = mp.setup('mose-multirate-moving-levels')
    if plt is None:
        print('  plot_moving_levels: matplotlib not available, figure skipped')
        return 0
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE2, sharey=False)
    handles = []
    n = None
    for lab, legend, colour, ls, m in RUNS:
        p = os.path.join(runs, lab + '.json')
        if not os.path.isfile(p):
            continue
        d = json.load(open(p))
        n = np.array([r['n'] for r in d])
        e1 = np.array([r['e1'][0] for r in d])
        ei = np.array([r['ei'][0] for r in d])
        o1, oi = d[-1]['o1'][0], d[-1]['oi'][0]
        lw = 3.0 if lab == 'mr' else 2.0
        ln, = axes[0].loglog(n, e1, ls, marker=m, ms=7, lw=lw, color=colour, label='%s (%.2f, %.2f)' % (legend, o1, oi))
        axes[1].loglog(n, ei, ls, marker=m, ms=7, lw=lw, color=colour)
        handles.append(ln)
    if n is None:
        print('  plot_moving_levels: no results, figure skipped')
        return 0
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    e1_mr = json.load(open(os.path.join(runs, 'mr.json')))
    a1 = e1_mr[0]['e1'][0]
    ai = e1_mr[0]['ei'][0]
    s2, = axes[0].loglog(n, a1 * (n[0] / n) ** 2, 'k--', lw=1.5, label='second-order slope')
    axes[1].loglog(n, ai * (n[0] / n) ** 2, 'k--', lw=1.5)
    s1, = axes[1].loglog(n, ai * (n[0] / n), 'k:', lw=1.5, label='first-order slope')
    for ax, ylabel in ((axes[0], 'L¹ error of ρ'), (axes[1], 'L∞ error of ρ')):
        mp.style(ax, 'cells n (Δx and H halved together)', ylabel, legend=False)
        ax.xaxis.set_major_locator(FixedLocator(n))
        ax.xaxis.set_major_formatter(FixedFormatter([str(int(v)) for v in n]))
        ax.xaxis.set_minor_locator(NullLocator())
    fig.tight_layout()
    fig.legend(handles=handles + [s2, s1], loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=1,
               fontsize=mp.LABEL_FS - 2, title_fontsize=mp.LABEL_FS - 2,
               title='simple wave on 256 to 2048 cells, T = 0.5 (observed order of the last pair: L¹, L∞);\n' + FORM)
    mp.save(plt, fig, os.path.join(out, 'multirate-moving-levels-order.svg'))
    print('  plot_moving_levels: %d runs in multirate-moving-levels-order.svg' % len(handles))
    return 0


if __name__ == '__main__':
    sys.exit(main())
