#!/usr/bin/env python3
"""Figures of the viscous shock tube of Daru and Tenaud with the multirate integrator (mr-daru) for its V&V page
(docs/vv/multirate.md), from the runs of run.sh in RUNS_DIR (ref, sr, mra, mrb) and the file of verify_daru.py --json.
Every form:
  multirate-daru-<form>-errors.svg  the density error e of sr, mra and mrb against the reference ref at the ten output
                                    times, and the ratio e(mrb)/e(sr) with its acceptance limit 1.05;
  multirate-daru-<form>-levels.svg  the share of the cells on the levels 1 to L at every macro step of mra and mrb,
                                    stacked (lines "cells per level" of their logs), and their cell updates against the
                                    single rate at dt0 (the work).
Full form only (500 x 250 cells, the resolution of the published solutions):
  multirate-daru-full-wall.svg      the density along the bottom wall at the final time, the four runs and table 1 of
                                    Zhou, Xu and Liu;
  multirate-daru-full-density.svg   the density field of mra at the final time.
<form> is short (50 x 25 cells, the form run in CI), medium (100 x 50) or full; every figure names its mesh and form.
Density as
rho/(rho_R/1.2), rho_R the density of the low-pressure state (1e5 Pa, 300 K), so that the low-pressure state reads 1.2
as in the non-dimensional problem of the references. Style of the V&V figures of MOSE (mr-common/mrplot.py). Without
matplotlib the figures are skipped; the test does not depend on them.
Usage: plot_daru.py RUNS_DIR CASE_DIR JSON [OUT_DIR]   (OUT_DIR defaults to RUNS_DIR)"""
import json
import os
import re
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'mr-common'))
import numpy as np                # noqa: E402
import mrplot as mp               # noqa: E402
import verify_daru as vd          # noqa: E402

LINE = {'ref': dict(mp.REFERENCE, zorder=2.6),          # the thick line below the thin ones, as compare.py draws them
        'sr': dict(color=mp.BLUE, ls='-.', lw=2.5, zorder=2.4),
        'mra': dict(color=mp.GREEN, ls='-', lw=2.0, zorder=2.2),
        'mrb': dict(color=mp.ORANGE, ls='-', lw=4.0, zorder=2.0)}
MARK = {'sr': 's', 'mra': 'D', 'mrb': 'o'}
LABEL = {'ref': 'single rate, dt₀/4 (reference)', 'sr': 'single rate, dt₀', 'mra': 'multirate, H = 4 dt₀',
         'mrb': 'multirate, H = dt₀'}                   # plain text: no mathtext pieces placed glyph by glyph
LEVEL_COLOURS = ('0.75', mp.BLUE, mp.GREEN, mp.ORANGE)     # level 0 (the rest) is not drawn
RHO_LABEL = 'ρ (non-dimensional)'
FIGSIZE2 = (10, 5)           # two panels in the width of a single-panel figure (10, 6): the same text size on the page
MAX_POINTS = 500             # points of a curve over the macro steps (the full form has up to 129 280 of them)


def solution(run, k):
    """the solution k (10 = final) of a run: one block (the runs of run.sh)"""
    return vd.solution(run, k)


def errors_figure(plt, J, head, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=FIGSIZE2)
    t = 1e3 * np.array(J['t'])
    for r in ('sr', 'mra', 'mrb'):
        a1.semilogy(t, J['e'][r], marker=MARK[r], ms=7, label=LABEL[r], **LINE[r])
    mp.style(a1, 't [ms]', 'density error e', legend_title=head)
    ratio = np.array(J['e']['mrb']) / np.array(J['e']['sr'])
    a2.plot(t, ratio, marker=MARK['mrb'], ms=7, label='e(mrb) / e(sr)', **LINE['mrb'])
    a2.axhline(J['limit'], color='k', ls=':', lw=2.0, label='acceptance limit %.2f' % J['limit'])
    a2.set_ylim(top=J['limit'] + 0.25 * (J['limit'] - min(1.0, float(ratio.min()))))
    mp.style(a2, 't [ms]', 'e(mrb) / e(sr)', legend_loc='center left')
    fig.tight_layout()
    mp.save(plt, fig, os.path.join(out, 'multirate-daru-%s-errors.svg' % J['form']))
    return float(ratio.max())


def level_counts(log):
    return [list(map(int, m.split())) for m in re.findall(r'Multirate: cells per level\s+([0-9 ]+)', log)]


def levels_figure(plt, runs, J, head, out):
    """share of the cells on the finer levels 1 to L at every macro step, stacked (level 0 holds the rest), and the cell
    updates of the run against those of the single rate at dt0 over the same time (the work)"""
    work = {}
    tfin = 1e3 * J['t'][-1]
    top = 0.0
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE2, sharey=True)
    for ax, r, f in zip(axes, ('mra', 'mrb'), (4, 1)):
        c = level_counts(open(os.path.join(runs, r, 'solve.log'), errors='replace').read())
        if len(c) < 2:
            plt.close(fig)
            return None
        L = max(len(x) for x in c)
        a = np.array([x + [0] * (L - len(x)) for x in c], dtype=float)
        n = a.sum(axis=1, keepdims=True)
        tt = tfin * (np.arange(len(c)) + 0.5) / float(len(c))   # [ms], the middle of every macro step
        share = a / n
        work[r] = float((a * (2.0 ** np.arange(L))).sum()) / (f * float(n.sum()))
        s = max(1, -(-len(c) // MAX_POINTS))                     # every s-th macro step: a figure of a few hundred kB
        if L > 1:
            ax.stackplot(tt[::s], share[::s, 1:].T, colors=LEVEL_COLOURS[1:L], labels=['level %d' % l for l in range(1, L)])
            top = max(top, float(share[:, 1:].sum(axis=1).max()))
        ax.set_xlim(0, tfin)
        mp.style(ax, 't [ms]', 'share of the cells on levels ≥ 1' if r == 'mra' else '', legend_loc='upper right',
                 legend_title='%s, work %.2f' % (LABEL[r], work[r]) + ('\n' + head if r == 'mra' else ''))
    axes[0].set_ylim(0.0, 1.6 * top if top > 0 else 1.0)
    fig.tight_layout()
    mp.save(plt, fig, os.path.join(out, 'multirate-daru-%s-levels.svg' % J['form']))
    return work


def wall_figure(plt, runs, J, head, out):
    fig, ax = plt.subplots(figsize=(10, 6))
    zx, zy = zip(*J['zhou'])
    ax.plot(zx, zy, label='Zhou, Xu and Liu (table 1)', zorder=3, **mp.PUBLISHED)
    for r in ('ref', 'sr', 'mra', 'mrb'):
        B = solution(os.path.join(runs, r), 10)
        xc = vd.mc.cell_centres(B['points'])[0][0, :, 0]
        ax.plot(xc, B['cell']['rho(1)'][0, 0, :] / J['rho_scale'], label=LABEL[r], **LINE[r])
    ax.set_xlim(0.3, 1.0)
    mp.style(ax, 'x / L', 'wall density, ' + RHO_LABEL, legend_loc='upper left', legend_title=head + ', final time')
    mp.save(plt, fig, os.path.join(out, 'multirate-daru-full-wall.svg'))


def density_figure(plt, runs, J, head, out):
    B = solution(os.path.join(runs, 'mra'), 10)
    x, y = B['points'][0, :, :, 0], B['points'][0, :, :, 1]
    q = B['cell']['rho(1)'][0] / J['rho_scale']
    fig, ax = plt.subplots(figsize=(10, 4.6))
    pc = ax.pcolormesh(x, y, q, shading='flat', cmap='turbo', rasterized=True)
    cb = fig.colorbar(pc, ax=ax, fraction=0.025, pad=0.03)
    cb.set_label(RHO_LABEL, fontsize=mp.LABEL_FS)
    cb.ax.tick_params(labelsize=mp.LABEL_FS - 2)
    ax.set_xlim(0.3, 1.0)
    ax.set_ylim(0.0, 0.3)
    ax.set_aspect('equal')
    mp.style(ax, 'x / L', 'y / L', legend=False)
    mp.label_box(ax, '%s, final time; %s' % (LABEL['mra'], head), above=True)
    mp.save(plt, fig, os.path.join(out, 'multirate-daru-full-density.svg'))
    return float(q.min()), float(q.max())


def main():
    runs, case, jpath = sys.argv[1], sys.argv[2], sys.argv[3]
    out = sys.argv[4] if len(sys.argv) > 4 else runs
    plt = mp.setup('mose-multirate-daru')
    if plt is None:
        print('  plot_daru: matplotlib not available, figures skipped')
        return 0
    J = json.load(open(jpath))
    full = J['form'] == 'full'
    head = mp.form_text(J['mesh'][0], J['mesh'][1], J['form'])
    worst = errors_figure(plt, J, head, out)
    work = levels_figure(plt, runs, J, head, out)
    msg = 'largest e(mrb)/e(sr) %.4f; cell updates against the single rate at dt0: %s' % (
        worst, ', '.join('%s %.3f' % kv for kv in sorted(work.items())) if work else 'not logged')
    if full:
        wall_figure(plt, runs, J, head, out)
        lo, hi = density_figure(plt, runs, J, head, out)
        msg += '; density of mra at t_fin %.1f to %.1f' % (lo, hi)
    print('  plot_daru (%s, %s): %s' % (J['case'], head, msg))
    return 0


if __name__ == '__main__':
    sys.exit(main())
