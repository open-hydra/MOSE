#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wall_pressure.py
================
Single-panel WALL PRESSURE figure for the HyShot II combustor V&V page:
MOSE, the Case 2 experiment with its uncertainty band, and the Spalart-Allmaras
RANS (TAU) of Karl's thesis, all along the between-injector wall line.

Two presets, selected with --preset:
  docs   (default)  transparent, theme-aware SVG written into docs/vv/images,
                    MOSE palette and sans-serif type -- matches every other V&V
                    figure in the documentation.
  paper             Times/STIX serif, vector PDF + raster PNG next to this
                    script, fonts embedded as TrueType (Type 42) because
                    publishers reject the matplotlib default (Type 3).

Everything the eye sees is in the STYLE dict below (and the DOCS overrides that
follow it): nothing is left to a matplotlib default.  Edit, re-run, nothing else.

DATA
  wall.tec                    MOSE wall solution (this case, OUTPUT/wall.tec)
  thesis_karl_exp[|_up|_down].txt
                              Case 2 experiment + uncertainty bounds, digitised
                              from Karl's thesis: x in mm from the COMBUSTOR
                              ENTRANCE, p normalised by p0 = 17.73 MPa.
  thesis_karl_rans_sa.txt     TAU Spalart-Allmaras, same thesis, same units.
  paper_fig3a_pressure_EXP.dat
                              the same measurements as digitised from Fig. 3a of
                              Chapuis et al., Proc. Combust. Inst. 34 (2013)
                              2101-2109, in (x[m], p[kPa]).  Used only by
                              check_exp() to audit the unit/origin conversion.

COORDINATE CONVENTION (from MESH/hyshot_mesh.py)
  x = streamwise, y = spanwise (y=0 through the injector axis, y=W half-pitch),
  z = wall-normal (z=0 bottom wall).  Injector axis at x = 0.408 m, y = 0, and
  W = 9.375 mm, so the between-injector symmetry plane is the model centreline.

MESH INDEPENDENCE
  Nothing here is tied to a particular block layout: the wall line is stitched
  from whatever bottom-wall zones wall.tec happens to contain (see extract_line),
  so a remesh -- more blocks, a different notch, the exit ramp that extends the
  domain from x = 0.65 to 0.76 m -- needs no edit.

REQUIREMENTS : numpy, matplotlib
USAGE        : python wall_pressure.py [--preset docs|paper] [--plot]
"""

import argparse
import os
import re

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
WALL_TEC = os.path.join(HERE, "wall.tec")
# The V&V figures are documentation figures: write them straight into the docs.
DOCS_IMAGES = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..",
                                            "docs", "vv", "images"))

# geometry / spanwise sampling line ---------------------------------------
W       = 0.009375          # spanwise half-pitch [m] (between-injector plane)
Y_PRESS = W                 # wall pressure is read on the between-injector line
X_INJ   = 0.408             # injector axis [m]
# row-acceptance tolerances, in LOCAL spanwise cells (see extract_line). They are
# relative because the nearest cell centre to a line ON the domain edge (y=W) is
# always half a cell away, whatever the mesh.
MAX_OFF, MAX_SPREAD = 1.0, 2.0

# Karl's thesis plots x in mm from the combustor entrance and p/p0; the paper
# (and the mesh) put x in m on the vehicle axis.  The entrance therefore sits at
# X0_THESIS, which is fixed by the geometry: the injector axis is 58 mm into the
# combustor, i.e. X_INJ - 0.058.
P0_THESIS = 17.73e6                 # [Pa]  normalising stagnation pressure
X0_THESIS = X_INJ - 0.058           # [m]   combustor entrance -> paper origin


# =============================================================================
# STYLE -- the whole appearance of the figure is here
# =============================================================================
STYLE = dict(
    # -- canvas ---------------------------------------------------------------
    # Elsevier column widths: 90 mm single, 140 mm 1.5-column, 190 mm double.
    # A shock train needs streamwise room, so 1.5-column is the honest minimum.
    width_mm      = 140.0,
    aspect        = 0.58,          # height / width
    png_dpi       = 600,           # raster preview only; the PDF is vector
    pad_inches    = 0.02,
    formats       = ("pdf", "png"),
    outdir        = HERE,
    stem          = "wall_pressure",
    transparent   = False,

    # -- type -----------------------------------------------------------------
    # First font on the list that exists wins; STIX is bundled with matplotlib
    # and is the Times-metric fallback, so this renders the same on any machine.
    font_family   = "serif",
    font_serif    = ["Times New Roman", "Nimbus Roman", "STIXGeneral", "DejaVu Serif"],
    mathtext      = "stix",
    fs_axis       = 12.0,          # axis labels
    fs_tick       = 11.0,
    fs_legend     = 11.0,
    use_tex       = False,         # True needs a working LaTeX install

    # -- axes chrome ----------------------------------------------------------
    spine_lw      = 0.6,
    tick_dir      = "out",
    tick_major    = 3.5,
    tick_minor    = 2.0,
    tick_lw       = 0.6,
    grid          = False,         # journals usually set this off; solid if on
    grid_lw       = 0.4,
    grid_color    = "#e6e6e6",

    # -- data limits and ticks (x in the unit chosen by x_origin, below) -------
    # x_origin: "paper"     -> x [m], vehicle axis, as in Chapuis et al.
    #           "combustor" -> x [mm] from the combustor entrance, as in the thesis
    x_origin      = "paper",
    xlim          = (0.355, 0.760),
    xtick_major   = 0.05,
    xtick_minor   = 0.01,
    ylim          = (0.0, 360.0),
    ytick_major   = 50.0,
    ytick_minor   = 10.0,

    # -- series ---------------------------------------------------------------
    # Dark red / forest green.  NOTE: this is the one pair a red-green colourblind
    # reader cannot separate by hue (they converge to two near-identical dull
    # olives), so the line STYLE carries the distinction -- solid vs dashed --
    # which also survives greyscale printing.
    cfd_color     = "#8B0000",     # dark red
    cfd_lw        = 1.5,
    cfd_ls        = "-",
    cfd_z         = 3,

    rans_color    = "#228B22",     # forest green
    rans_lw       = 1.1,
    rans_ls       = (0, (5.5, 2.2)),
    rans_z        = 2,

    # measured data is ink, not a categorical hue: it is told apart by marker
    # shape and error bars, so it does not compete with the two computed curves
    exp_color     = "#000000",
    exp_ms        = 4.2,
    exp_mew       = 0.8,
    exp_marker    = "s",
    exp_elw       = 0.7,
    exp_capsize   = 1.8,
    exp_capthick  = 0.7,
    exp_z         = 4,

    # -- legend ---------------------------------------------------------------
    # One row ABOVE the axes.  Inside is not an option here: the error bars run
    # to 338 kPa at x = 0.647 m and the post-exit expansion fills the lower
    # right, so every in-axes corner collides with data at some ylim.  Out of
    # the axes it cannot collide at all, and the row costs less height than the
    # headroom an in-axes legend would need.  Set legend_above = False to go
    # back inside, at legend_loc.
    legend_above  = True,
    legend_bbox   = (0.0, 1.005),  # axes fraction, used when legend_above
    legend_loc    = "upper right",  # used when not legend_above
    legend_frame  = False,
    legend_ncol   = 3,             # 3 -> single row above; use 1 when inside
    legend_handle = 2.6,           # handle length [font units]
    legend_labelspacing = 0.35,
    legend_colspacing   = 1.6,

    # -- labels ---------------------------------------------------------------
    lbl_cfd       = "MOSE",
    lbl_rans      = "TAU (Karl)",
    lbl_exp       = "Experiment",
    xlabel_paper     = r"$x$ [m]",
    xlabel_combustor = r"$x$ [mm]",
    ylabel        = r"$p_w$ [kPa]",

    # -- provenance -----------------------------------------------------------
    # The paper file holds one extra measurement (x = 0.7572 m) that the thesis
    # figure does not.  Off by default: a figure should not mix two digitisations
    # of the same experiment under one symbol, and that point sits outside the
    # region the figure is about.  See check_exp().
    include_paper_only_exp = False,
)

# -----------------------------------------------------------------------------
# DOCS preset: what has to change so the figure sits among the other V&V ones.
# Transparent canvas and real <text> elements, so docs/stylesheets/extra.css can
# recolor the labels for the page theme; MOSE orange is the hero colour and the
# reference code takes the categorical blue, as in the SWBLI figures.  Type is
# sans (the docs body face) and larger, because these are read on screen at
# ~700 px wide, not printed at 140 mm.
# -----------------------------------------------------------------------------
DOCS = dict(
    formats       = ("svg",),
    outdir        = DOCS_IMAGES,
    stem          = "HyShot-pw",
    transparent   = True,
    width_mm      = 254.0,          # 10 in, the figsize the other V&V scripts use
    aspect        = 0.6,
    font_family   = "sans-serif",
    mathtext      = "dejavusans",
    fs_axis       = 17.0,
    fs_tick       = 15.0,
    fs_legend     = 15.0,
    spine_lw      = 1.0,
    tick_lw       = 1.0,
    tick_major    = 5.0,
    tick_minor    = 3.0,
    grid          = True,
    grid_lw       = 0.8,
    grid_color    = "#9a9a9a",
    cfd_color     = "#eb6834",      # MOSE orange (hero)
    cfd_lw        = 3.5,
    rans_color    = "#2a78d6",      # reference code blue
    rans_lw       = 2.2,
    rans_ls       = (0, (6.0, 2.5)),
    exp_ms        = 6.0,
    exp_mew       = 1.3,
    exp_elw       = 1.1,
    exp_capsize   = 3.0,
    exp_capthick  = 1.1,
    lbl_rans      = "TAU S-A (Karl)",
    lbl_exp       = "Experiment (Karl)",
    xlabel_paper     = r"$x$  [m]",
    xlabel_combustor = r"$x$  [mm]",
    ylabel        = r"$p_w$  [kPa]",
)


# =============================================================================
# 1) MOSE  --  parse wall.tec and extract the wall-pressure line
# =============================================================================
def parse_tecplot(path):
    """Tecplot reader: DATAPACKING=BLOCK, nodal x,y,z + cell-centred fields."""
    with open(path) as f:
        lines = f.readlines()
    variables = re.findall(r'"([^"]+)"', lines[0])

    zones = []
    for i, ln in enumerate(lines):
        if ln.strip().startswith("ZONE"):
            m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", ln)
            T = re.search(r"T\s*=\s*(\S+?),", ln).group(1)
            v = re.search(r"\[1-(\d+)\]\s*=\s*NODAL", ln)      # e.g. [1-3]=NODAL
            nod = int(v.group(1)) if v else 3
            zones.append((i, T, int(m.group(1)), int(m.group(2)), int(m.group(3)), nod))
    zones.append((len(lines), None, 0, 0, 0, 0))  # sentinel

    def read_floats(seg):
        out = []
        for ln in seg:
            s = ln.strip()
            if s:
                out.extend(float(v) for v in s.split())
        return np.array(out)

    parsed = {}
    for zi in range(len(zones) - 1):
        start, T, I, J, K, nod = zones[zi]
        end = zones[zi + 1][0]
        data = read_floats(lines[start + 1:end])
        nnode = I * J * K
        # a surface zone is FLAT in one index, and that one is not always K: the
        # injector pipe walls come out as I=1 or J=1. Every degenerate direction
        # holds one cell, not zero.
        ncell = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)
        idx, arrs = 0, {}
        for vi, v in enumerate(variables):
            n = nnode if vi < nod else ncell       # x,y,z nodal, fields cell-centred
            arrs[v] = data[idx:idx + n]
            idx += n
        assert idx == len(data), f"{T}: parsed {idx} != {len(data)}"
        parsed[T] = (I, J, K, arrs)
    return parsed


def cell_centers(parsed, T):
    """Return (xc, yc, pw, qw, Tw) as (J-1, I-1) cell-centred arrays."""
    I, J, K, a = parsed[T]
    X = a["x"].reshape(J, I); Y = a["y"].reshape(J, I)
    xc = 0.25 * (X[:-1, :-1] + X[:-1, 1:] + X[1:, :-1] + X[1:, 1:])
    yc = 0.25 * (Y[:-1, :-1] + Y[:-1, 1:] + Y[1:, :-1] + Y[1:, 1:])
    pw = a["pw"].reshape(J - 1, I - 1)
    qw = a["qw"].reshape(J - 1, I - 1)
    Tw = a["Tw"].reshape(J - 1, I - 1)
    return xc, yc, pw, qw, Tw


def bottom_wall_zones(parsed, tol=1e-9):
    """Names of the surface zones lying in the flat bottom wall (z = const, the
    lowest such plane).  Found from the coordinates, so no zone list to maintain:
    it drops the ceiling (which the exit ramp makes non-flat anyway) and the
    injector pipe walls (which span z) on their own."""
    flat = {T: a["z"][0] for T, (I, J, K, a) in parsed.items()
            if K == 1 and np.ptp(a["z"]) < tol}
    if not flat:
        raise RuntimeError(f"{WALL_TEC}: no flat K=1 surface zone found")
    z0 = min(flat.values())
    return [T for T, z in flat.items() if abs(z - z0) < tol]


def extract_line(parsed, y_target, verbose=True):
    """Stitch a streamwise bottom-wall line at fixed y across the blocks.

    Every bottom-wall zone offers the cell row whose mean centre-y is closest to
    y_target, and is kept only if that row really is a sample of the requested
    line: its centre within MAX_OFF local spanwise cells of y_target, and flat
    (y-spread over the row) to within MAX_SPREAD of them.  That test is what
    makes this mesh-agnostic -- it accepts the x-aligned rows of the background
    and notch grids and rejects the O-grid ring blocks around the injector, whose
    "rows" run RADIALLY (y sweeps the whole block along one row) and would
    otherwise smear the profile.  Here it keeps the rows at 0.1-0.5 cells and
    throws out the radial ones at 6-78, so the margin is an order of magnitude.

    Surviving rows are laid down closest-first and each only fills the x-gaps its
    predecessors left, so overlapping chimera grids are never double counted.
    """
    cand = []
    for T in bottom_wall_zones(parsed):
        xc, yc, pw, qw, Tw = cell_centers(parsed, T)
        ym = yc.mean(axis=1)
        jc = int(np.argmin(np.abs(ym - y_target)))
        j2 = jc + 1 if jc + 1 < len(ym) else jc - 1          # local spanwise size
        dy = abs(ym[j2] - ym[jc])
        if dy <= 0.0:
            continue
        off, spread = abs(ym[jc] - y_target) / dy, np.ptp(yc[jc]) / dy
        if off > MAX_OFF or spread > MAX_SPREAD:
            continue
        cand.append((off, T, xc[jc], pw[jc], qw[jc], Tw[jc]))
    cand.sort(key=lambda c: c[0])

    xs, ps, qs, ts, covered = [], [], [], [], []
    for err, T, xrow, p, q, t in cand:
        keep = np.ones(xrow.shape, bool)
        for lo, hi in covered:                  # only fill what is still missing
            keep &= (xrow < lo) | (xrow > hi)
        if not keep.any():
            continue
        covered.append((xrow.min(), xrow.max()))
        xs.append(xrow[keep]); ps.append(p[keep]); qs.append(q[keep]); ts.append(t[keep])
        if verbose:
            print(f"    y={y_target*1e3:6.3f} mm : {T:6s} {keep.sum():5d} pts  "
                  f"x[{xrow[keep].min():.4f},{xrow[keep].max():.4f}]  {err:.2f} cells off")
    if not xs:
        raise RuntimeError(f"no bottom-wall row lies along y={y_target}")
    x = np.concatenate(xs); o = np.argsort(x)
    return (x[o], np.concatenate(ps)[o], np.concatenate(qs)[o], np.concatenate(ts)[o])


# =============================================================================
# 2) REFERENCE DATA  --  thesis .txt (converted on load) + paper .dat
# =============================================================================
def _load(name, sort=True):
    a = np.loadtxt(os.path.join(HERE, name))       # 2 cols: x, value ; '#' header
    if sort:
        a = a[np.argsort(a[:, 0])]
    return a


def _load_thesis(name):
    """Thesis file (x[mm] from combustor entrance, p/p0) -> (x[m], p[kPa])."""
    a = _load(name)
    return np.column_stack([a[:, 0] * 1e-3 + X0_THESIS, a[:, 1] * P0_THESIS / 1e3])


def _load_thesis_exp(base="thesis_karl_exp", tol_x=1.0):
    """Thesis experiment + uncertainty -> (points, yerr) with yerr as (2, N).

    ``base``, ``base_up`` and ``base_down`` were digitised as three separate
    curves, so their x columns are the same measurement stations to within the
    digitisation error only.  The bounds are taken ROW-WISE -- point i of _up
    belongs to point i of the centre file -- which is what the tol_x check
    guards: if the files ever differ in length or drift apart in x by more than
    tol_x mm, the rows no longer line up and the bars would be attached to the
    wrong stations.  x comes from the centre file.
    """
    c, up, dn = (_load_thesis(f"{base}{s}.txt") for s in ("", "_up", "_down"))
    if not len(c) == len(up) == len(dn):
        raise RuntimeError(f"{base}: {len(c)}/{len(up)}/{len(dn)} rows, cannot pair")
    for tag, b in (("_up", up), ("_down", dn)):
        off = np.abs(b[:, 0] - c[:, 0]).max() * 1e3
        if off > tol_x:
            raise RuntimeError(f"{base}{tag}: rows drift {off:.2f} mm > {tol_x} mm in x")
    yerr = np.vstack([c[:, 1] - dn[:, 1], up[:, 1] - c[:, 1]])
    if (yerr < 0).any():
        raise RuntimeError(f"{base}: _up/_down do not bracket the measurement")
    return c, yerr


def load_reference():
    return dict(
        A_exp=_load("paper_fig3a_pressure_EXP.dat"),
        A_karl_rans=_load_thesis("thesis_karl_rans_sa.txt"),
        **dict(zip(("A_karl_exp", "A_karl_exp_err"), _load_thesis_exp())),
    )


def check_exp(ref, tol_x=2.0, tol_p=2.0):
    """Confirm the thesis experiment is the paper experiment in other units.

    Both sets are digitisations of the same Case 2 measurements, so this is a
    units/origin check, not a physics one: it pairs each thesis point with its
    nearest paper point and reports the worst x- and p-disagreement.  Anything
    above a digitisation error (tol_x mm, tol_p %) means X0_THESIS or P0_THESIS
    is wrong, and it warns rather than raises so the figure still gets drawn.

    Returns the paper points with no thesis partner (the thesis figure stops
    earlier).  They are not errors -- see include_paper_only_exp in STYLE.
    """
    karl, paper = ref["A_karl_exp"], ref["A_exp"]
    j = np.array([np.argmin(np.abs(paper[:, 0] - x)) for x in karl[:, 0]])
    dx = (karl[:, 0] - paper[j, 0]) * 1e3                      # [mm]
    dp = 100 * (karl[:, 1] - paper[j, 1]) / paper[j, 1]        # [%]
    extra = sorted(set(range(len(paper))) - set(j.tolist()))
    print(f"  exp check: {len(karl)}/{len(paper)} paper points matched, "
          f"max |dx| = {np.abs(dx).max():.2f} mm, max |dp| = {np.abs(dp).max():.2f} %")
    if extra:
        print("    paper-only points at x = "
              + ", ".join(f"{paper[i, 0]:.4f} m" for i in extra))
    if np.abs(dx).max() > tol_x or np.abs(dp).max() > tol_p:
        print("    WARNING: beyond digitisation scatter -- check X0_THESIS/P0_THESIS")
    return paper[extra]


# =============================================================================
# 3) FIGURE
# =============================================================================
def apply_style(s):
    mpl.rcParams.update({
        "font.family": s["font_family"],
        "font.serif": s["font_serif"],
        "mathtext.fontset": s["mathtext"],
        "text.usetex": s["use_tex"],
        "axes.labelsize": s["fs_axis"],
        "axes.linewidth": s["spine_lw"],
        "xtick.labelsize": s["fs_tick"], "ytick.labelsize": s["fs_tick"],
        "xtick.direction": s["tick_dir"], "ytick.direction": s["tick_dir"],
        "xtick.bottom": True, "ytick.left": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": s["tick_major"], "ytick.major.size": s["tick_major"],
        "xtick.minor.size": s["tick_minor"], "ytick.minor.size": s["tick_minor"],
        "xtick.major.width": s["tick_lw"], "ytick.major.width": s["tick_lw"],
        "xtick.minor.width": s["tick_lw"], "ytick.minor.width": s["tick_lw"],
        "legend.fontsize": s["fs_legend"],
        "legend.frameon": s["legend_frame"],
        "legend.handlelength": s["legend_handle"],
        "legend.labelspacing": s["legend_labelspacing"],
        "lines.solid_capstyle": "round",
        # Type 42 = TrueType.  Publishers reject the matplotlib default (Type 3).
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.bbox": "tight", "savefig.pad_inches": s["pad_inches"],
    })
    if s["transparent"]:
        # Transparent canvas and real <text> elements, so the docs CSS can
        # recolor them for the page theme (see docs/stylesheets/extra.css).
        mpl.rcParams.update({
            "figure.facecolor": "none", "axes.facecolor": "none",
            "savefig.facecolor": "none", "svg.fonttype": "none",
        })


def to_x(x_m, s):
    """Metres on the vehicle axis -> the abscissa the figure is drawn in."""
    return (x_m - X0_THESIS) * 1e3 if s["x_origin"] == "combustor" else x_m


def build(s=STYLE):
    apply_style(s)

    # -- data ------------------------------------------------------------------
    x_cfd, p_cfd, _, _ = extract_line(parse_tecplot(WALL_TEC), Y_PRESS, verbose=False)
    p_cfd = p_cfd / 1e3                                     # Pa -> kPa
    ref = load_reference()
    extra = check_exp(ref)                                  # also prints the audit
    exp, exp_err = ref["A_karl_exp"], ref["A_karl_exp_err"]
    rans = ref["A_karl_rans"]

    w = s["width_mm"] / 25.4
    fig, ax = plt.subplots(figsize=(w, w * s["aspect"]))

    ax.plot(to_x(rans[:, 0], s), rans[:, 1], ls=s["rans_ls"], color=s["rans_color"],
            lw=s["rans_lw"], zorder=s["rans_z"], label=s["lbl_rans"])
    ax.plot(to_x(x_cfd, s), p_cfd, ls=s["cfd_ls"], color=s["cfd_color"],
            lw=s["cfd_lw"], zorder=s["cfd_z"], label=s["lbl_cfd"])
    # yerr is (2, N): the thesis bounds are not symmetric about the measurement
    ax.errorbar(to_x(exp[:, 0], s), exp[:, 1], yerr=exp_err,
                marker=s["exp_marker"], ls="none", mfc="none", mec=s["exp_color"],
                ms=s["exp_ms"], mew=s["exp_mew"], ecolor=s["exp_color"],
                elinewidth=s["exp_elw"], capsize=s["exp_capsize"],
                capthick=s["exp_capthick"], zorder=s["exp_z"], label=s["lbl_exp"])
    if s["include_paper_only_exp"] and len(extra):
        ax.plot(to_x(extra[:, 0], s), extra[:, 1], marker=s["exp_marker"], ls="none",
                mfc="none", mec=s["exp_color"], ms=s["exp_ms"], mew=s["exp_mew"],
                zorder=s["exp_z"])

    xlim = s["xlim"] if s["x_origin"] == "paper" else tuple(to_x(np.array(s["xlim"]), s))
    ax.set_xlim(*xlim); ax.set_ylim(*s["ylim"])
    ax.xaxis.set_major_locator(MultipleLocator(s["xtick_major"]))
    ax.xaxis.set_minor_locator(MultipleLocator(s["xtick_minor"]))
    ax.yaxis.set_major_locator(MultipleLocator(s["ytick_major"]))
    ax.yaxis.set_minor_locator(MultipleLocator(s["ytick_minor"]))
    ax.set_xlabel(s["xlabel_paper"] if s["x_origin"] == "paper" else s["xlabel_combustor"])
    ax.set_ylabel(s["ylabel"])
    if s["grid"]:
        ax.grid(True, ls="-", lw=s["grid_lw"], color=s["grid_color"], alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    # order the legend as the reader meets the curves: data, then the two models
    h, l = ax.get_legend_handles_labels()
    order = [l.index(s[k]) for k in ("lbl_exp", "lbl_rans", "lbl_cfd")]
    leg = dict(ncol=s["legend_ncol"], columnspacing=s["legend_colspacing"])
    if s["legend_above"]:
        leg.update(loc="lower left", bbox_to_anchor=s["legend_bbox"], borderaxespad=0.0)
    else:
        leg.update(loc=s["legend_loc"])
    ax.legend([h[i] for i in order], [l[i] for i in order], **leg)
    return fig


def save(fig, s):
    os.makedirs(s["outdir"], exist_ok=True)
    for ext in s["formats"]:
        out = os.path.join(s["outdir"], f"{s['stem']}.{ext}")
        fig.savefig(out, dpi=s["png_dpi"] if ext == "png" else None,
                    transparent=s["transparent"])
        print(f"  wrote {os.path.relpath(out, os.getcwd())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--preset", choices=("docs", "paper"), default="docs",
                    help="docs: theme-aware SVG into docs/vv/images (default); "
                         "paper: Times PDF + PNG next to this script")
    ap.add_argument("--plot", action="store_true", help="also show the figure")
    args = ap.parse_args()

    s = dict(STYLE)
    if args.preset == "docs":
        s.update(DOCS)
    if not args.plot:
        mpl.use("Agg")

    fig = build(s)
    save(fig, s)
    if args.plot:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
