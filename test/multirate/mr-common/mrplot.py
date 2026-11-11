"""Style of the figures of the multirate tests for the V&V pages (docs/vv/multirate.md): the style of the V&V figures of
MOSE, as test/2D/viscous/swbli/compare.py draws them:
  - transparent background and text kept as text (svg.fonttype none), so that the pages draw the text in their own
    colour (docs/stylesheets/extra.css);
  - axis labels at 17 pt, tick labels and legend at 15 pt, legend where it hides the least, grid at alpha 0.3;
  - 150 dpi, tight bounding box, transparent;
  - the colours of compare.py: orange for the main MOSE curve (thick), green and blue for the other runs, a MOSE
    reference run in grey dashes, measured or published values as black dots.
Like those figures, no figure has a title: every figure names its mesh and its form in the title of its legend (or in
a box drawn as a legend, for a field), and the page names them in the caption; the short form of a test is the one run
in CI, the medium and the full forms run outside it, and the full form is the one compared with published results.
Text is plain (Unicode, as dt₀): with svg.fonttype none
a mathtext expression is written as pieces placed one by one, which overlap when the page font has other metrics. The
SVG files carry no date and a fixed hash salt, so that the same results give the same files."""

LABEL_FS = 17
ORANGE, GREEN, BLUE = '#eb6834', '#1baf7a', '#2a78d6'          # MOSE_C, OF_C and WIND_C of compare.py
REFERENCE = dict(color='0.35', ls='--', lw=2.0)                # "MOSE reference" of compare.py
PUBLISHED = dict(ls='none', marker='o', color='k', ms=6)       # "Schulein (exp)" of compare.py


def setup(salt):
    """matplotlib with the Agg backend and the settings of the V&V figures of MOSE; None without matplotlib"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    matplotlib.rcParams.update({
        'figure.facecolor': 'none',
        'axes.facecolor': 'none',
        'savefig.facecolor': 'none',
        'svg.fonttype': 'none',
        'svg.hashsalt': salt,
    })
    return plt


def style(ax, xlabel, ylabel, legend=True, legend_loc='best', legend_title=None):
    ax.set_xlabel(xlabel, fontsize=LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.tick_params(labelsize=LABEL_FS - 2)
    if legend:
        ax.legend(loc=legend_loc, fontsize=LABEL_FS - 2, title=legend_title, title_fontsize=LABEL_FS - 2)
    ax.grid(True, which='both', alpha=0.3)


def form_text(ni, nj, form):
    """the mesh and the form (short, medium or full), for the title of a legend"""
    return '%d × %d cells, %s' % (ni, nj, {'short': 'short form (CI)', 'medium': 'medium form', 'full': 'full form'}[form])


def label_box(ax, text, loc='upper left', above=False):
    """a box drawn as a legend (frame, font and place), for a figure without lines to label, as a field; above=True
    puts it just above the axes, so that it hides nothing of the field"""
    from matplotlib.offsetbox import AnchoredText
    if above:
        box = AnchoredText(text, loc='lower left', prop=dict(size=LABEL_FS - 2), frameon=True, borderpad=0.0,
                           bbox_to_anchor=(0.0, 1.02), bbox_transform=ax.transAxes)
    else:
        box = AnchoredText(text, loc=loc, prop=dict(size=LABEL_FS - 2), frameon=True)
    box.patch.set_boxstyle('round,pad=0.3,rounding_size=0.2')
    box.patch.set_facecolor('white')
    box.patch.set_alpha(0.8)
    box.patch.set_edgecolor('0.8')
    ax.add_artist(box)


def save(plt, fig, path):
    fig.savefig(path, dpi=150, bbox_inches='tight', transparent=True, metadata={'Date': None})
    plt.close(fig)
    print('  wrote %s' % path)
