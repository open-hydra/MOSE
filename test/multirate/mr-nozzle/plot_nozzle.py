#!/usr/bin/env python3
"""Figures of the start-up of a nozzle with the multirate integrator (mr-nozzle) for its V&V page (docs/vv/multirate.md),
from the file of verify_nozzle.py --json:
  multirate-nozzle-short-shock.svg   the position of the moving shock (face of largest compression in the divergent) at
                                     the output times where it is a shock (the outputs judged by the test) for ref, sr,
                                     sr2, mra and mra2, with the faces between the bands of cells (between levels);
  multirate-nozzle-short-errors.svg  the density error in the throat band of sr, sr2, mrb, mra and mra2 against ref (the
                                     values above round-off, 1e-12), and
                                     the ratio e(mrb)/e(sr) that the test bounds by 1.05 where mrb has levels.
The form run in CI; every figure names its mesh and form in its legend. Style of the V&V figures of MOSE
(mr-common/mrplot.py). Without matplotlib the figures are skipped; the test does not depend on them.
Usage: plot_nozzle.py JSON [OUT_DIR]   (OUT_DIR defaults to the directory of JSON)"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mr-common'))
import numpy as np        # noqa: E402
import mrplot as mp       # noqa: E402

LINE = {'ref': dict(mp.REFERENCE, zorder=2.6),
        'sr': dict(color=mp.BLUE, ls='-.', lw=2.5, zorder=2.4),
        'sr2': dict(color=mp.BLUE, ls=':', lw=2.5, zorder=2.4),
        'mrb': dict(color='#E69F00', ls='-', lw=2.5, zorder=2.3),
        'mra': dict(color=mp.ORANGE, ls='-', lw=4.0, zorder=2.0),
        'mra2': dict(color=mp.GREEN, ls='-', lw=2.0, zorder=2.2)}
MARK = {'ref': 'o', 'sr': 's', 'sr2': 'v', 'mrb': 'P', 'mra': 'o', 'mra2': 'D'}
LABEL = {'ref': 'single rate, dt₀/4 (reference)', 'sr': 'single rate, dt₀', 'sr2': 'single rate, dt₀/2',
         'mrb': 'multirate, H = dt₀', 'mra': 'multirate, production H', 'mra2': 'multirate, production H/2'}
BAND_FACES = (0.700, 0.912)       # faces between the bands of 1.6, 4 and 10 mm in the divergent (between levels)
FIGSIZE2 = (10, 5)


def shock_figure(plt, J, head, out):
    """the shock at the outputs where it is one (the jump of ref across it at least 1.5, the outputs judged by D4)"""
    fig, ax = plt.subplots(figsize=(10, 6))
    k = J['judged_d4']
    t = 1e3 * np.array(J['t'])[k]
    for xb in BAND_FACES:
        ax.axhline(xb, color='0.6', lw=0.8)
    for r in ('ref', 'sr', 'sr2', 'mra', 'mra2'):
        ax.plot(t, np.array(J['shock_x'][r])[k], marker=MARK[r], ms=7, label=LABEL[r], **LINE[r])
    ax.set_xlim(0.0, 1e3 * J['t'][-1] * 1.03)
    ax.set_ylim(0.5, 1.0)
    mp.style(ax, 't [ms]', 'shock position x [m]', legend_loc='upper left',
             legend_title=head + '\ngrey lines: faces between bands of cells')
    mp.save(plt, fig, os.path.join(out, 'multirate-nozzle-short-shock.svg'))


def errors_figure(plt, J, head, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=FIGSIZE2)
    t = 1e3 * np.array(J['t'])
    for r in ('sr', 'sr2', 'mrb', 'mra', 'mra2'):
        e = np.array(J['e'][r])
        a1.semilogy(t[e > 1e-12], e[e > 1e-12], marker=MARK[r], ms=6, label=LABEL[r], **LINE[r])   # above round-off
    a1.set_xlim(0.0, t[-1] * 1.03)
    mp.style(a1, 't [ms]', 'density error in the throat', legend=False)
    k = [i for i in range(10) if J['mrb_fine'][i] > 0]
    ratio = [J['e']['mrb'][i] / J['e']['sr'][i] for i in k]
    a2.plot(t[k], ratio, marker=MARK['mrb'], ms=7, label='e(mrb) / e(sr)', **LINE['mrb'])
    a2.axhline(J['limit'], color='k', ls=':', lw=2.0, label='acceptance limit %.2f' % J['limit'])
    a2.set_xlim(a1.get_xlim())
    mp.style(a2, 't [ms]', 'e(mrb) / e(sr)', legend_loc='lower right')
    fig.tight_layout()
    h, lab = a1.get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=mp.LABEL_FS - 2,
               title=head, title_fontsize=mp.LABEL_FS - 2)
    mp.save(plt, fig, os.path.join(out, 'multirate-nozzle-short-errors.svg'))
    return max(ratio) if ratio else float('nan')


def main():
    jpath = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(jpath))
    plt = mp.setup('mose-multirate-nozzle')
    if plt is None:
        print('  plot_nozzle: matplotlib not available, figures skipped')
        return 0
    J = json.load(open(jpath))
    head = '%s, %d cells, short form (CI)' % (J['case'], J['cells'])
    shock_figure(plt, J, head, out)
    worst = errors_figure(plt, J, head, out)
    print('  plot_nozzle (%s): largest e(mrb)/e(sr) where mrb has levels %.4f' % (head, worst))
    return 0


if __name__ == '__main__':
    sys.exit(main())
