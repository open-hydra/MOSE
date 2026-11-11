#!/usr/bin/env python3
"""Figures of the multirate order test (mr-order) for its V&V page (docs/vv/multirate.md), written from the results of
run.sh in RUNS_DIR (the files <case>-order.json of verify_order.py and the runs themselves):
  multirate-order.svg        the L1 error of the density against the macro step H for every case, relative to the error
                             with the largest macro step H0, and the interface error (Linf on the cells within 4 cells
                             of a face between levels, with the floor of 5 % of the global Linf), with the slope of
                             second order; the mesh and the observed orders of the last pair in the legend;
  multirate-wall-s-levels.svg the level map of mr2d-wall-s (the run with H0 and trace = .true.) on its skewed mesh: the
                             levels 0 | 1 | 2 | 1 | 0 along the wall, the faces between levels crossing the wall.
The figures come from the full form (MR_FORM=full); the legends say so. Style of the V&V figures of MOSE (mr-common/mrplot.py). Without
matplotlib the figures are skipped; the test does not depend on them.
Usage: plot_order.py RUNS_DIR [OUT_DIR]   (OUT_DIR defaults to RUNS_DIR)"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrcheck as mc      # noqa: E402
import mrplot as mp       # noqa: E402

CASES = [  # case, label, group (colour), marker
    ('mr1d-acu', '1D acoustic wave, L = 3', 0, 'o'),
    ('mr1d-ent', '1D entropy wave, L = 3', 0, 's'),
    ('mr2d-vortex-3b', '2D vortex, three blocks, L = 1', 1, '^'),
    ('mr2d-plate', 'laminar flat plate, L = 3', 2, 'D'),
    ('mr2d-hot', 'skewed hot square, L = 1', 2, 'v'),
    ('mr2d-wall-s', 'wall on skewed cells, L = 2', 2, 'P'),
    ('mr2d-chan-sa', 'channel, Spalart-Allmaras, L = 3', 3, 'X'),
    ('mr2d-chan-sst', 'channel, SST, L = 3', 3, '*'),
]
COLOURS = (mp.BLUE, mp.GREEN, mp.ORANGE, '#E69F00')   # compare.py, and the CFL3D colour of the turbulent plate
LEVEL_COLOURS = ('0.75', mp.BLUE, mp.GREEN, mp.ORANGE)
FIGSIZE2 = (10, 5)        # two panels in the width of a single-panel figure (10, 6): the same text size on the page
FORM = 'full form'


def cells(runs, case):
    """number of cells of a case (its run with H0), or None"""
    run = os.path.join(runs, case + '-h0')
    if not os.path.isdir(run):
        return None
    return int(sum(np.prod(b['cell']['p'].shape) for b in mc.read_run(run)))


def order_figure(plt, runs, out):
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE2, sharey=True)
    handles = []
    h = np.array([1.0, 0.5, 0.25, 0.125])
    found = 0
    for case, label, g, m in CASES:
        p = os.path.join(runs, case + '-order.json')
        if not os.path.isfile(p):
            continue
        d = json.load(open(p))['vars'].get('rho(1)')
        if d is None:
            continue
        found += 1
        n = cells(runs, case)
        for ax, key in ((axes[0], 'L1'), (axes[1], 'Lint4')):
            e = np.array(d[key])
            ln, = ax.loglog(h, e / e[0], '-', marker=m, ms=7, lw=2.5, color=COLOURS[g],
                            label='%s, %s (%.2f, %.2f)' % (label, '%d cells' % n if n else 'cells not read',
                                                         d['order_L1'], d['order_int']))
        handles.append(ln)
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    for ax, ylabel in ((axes[0], 'L¹ error of ρ, all cells / at H₀'),
                       (axes[1], 'L∞ error of ρ near the faces / at H₀')):
        ref, = ax.loglog(h, h ** 2, 'k--', lw=1.5, label='second-order slope')
        mp.style(ax, 'macro step H/H₀', ylabel, legend=False)
        ax.xaxis.set_major_locator(FixedLocator(h))
        ax.xaxis.set_major_formatter(FixedFormatter(['1', '1/2', '1/4', '1/8']))
        ax.xaxis.set_minor_locator(NullLocator())
    axes[1].yaxis.set_tick_params(labelleft=False)
    fig.tight_layout()
    fig.legend(handles=handles + [ref], loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=1,
               fontsize=mp.LABEL_FS - 2, title_fontsize=mp.LABEL_FS - 2,
               title='case, cells (observed order of the last pair: L¹, near the faces between levels);\n' + FORM)
    mp.save(plt, fig, os.path.join(out, 'multirate-order.svg'))
    return found


def wall_figure(plt, runs, out):
    run = os.path.join(runs, 'mr2d-wall-s-h0')
    if not os.path.isdir(run):
        return 0
    blocks = mc.read_run(run)
    lev = mc.read_levels(run, [b['cell']['p'].shape for b in blocks])[0][0]          # (nj, ni) of the only block
    P = blocks[0]['points'][0]                                                       # (nj+1, ni+1, 3) nodes
    x, y = 1e3 * P[..., 0], 1e3 * P[..., 1]
    from matplotlib.colors import BoundaryNorm, ListedColormap
    L = int(lev.max())
    cmap = ListedColormap(list(LEVEL_COLOURS[:L + 1]))
    norm = BoundaryNorm(np.arange(-0.5, L + 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(10, 4.2))
    pc = ax.pcolormesh(x, y, lev, cmap=cmap, norm=norm, edgecolors='white', linewidth=0.3, shading='flat')
    ax.plot(x[0, :], y[0, :], 'k-', lw=3.0)                                           # the no-slip wall (301)
    ax.set_aspect('equal')
    mp.style(ax, 'x [mm]', 'y [mm]', legend=False)
    cb = fig.colorbar(pc, ax=ax, ticks=range(L + 1), fraction=0.03, pad=0.03)
    cb.set_label('level', fontsize=mp.LABEL_FS)
    cb.ax.tick_params(labelsize=mp.LABEL_FS - 2)
    nj, ni = lev.shape
    mp.label_box(ax, 'mr2d-wall-s, %d × %d cells, the run with H₀; %s' % (ni, nj, FORM), above=True)
    mp.save(plt, fig, os.path.join(out, 'multirate-wall-s-levels.svg'))
    return 1


def main():
    runs = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else runs
    plt = mp.setup('mose-multirate-order')
    if plt is None:
        print('  plot_order: matplotlib not available, figures skipped')
        return 0
    n = order_figure(plt, runs, out)
    w = wall_figure(plt, runs, out)
    print('  plot_order: %d cases in multirate-order.svg; level map of mr2d-wall-s %s' % (n, 'written' if w else 'skipped'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
