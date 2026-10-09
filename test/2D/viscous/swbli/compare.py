"""
SWBLI: compare a MOSE run with the reference wall solutions, and regenerate the
V&V documentation figures. MOSE and OpenFOAM (rhoCentralFoam) are run on the same
237,440-cell grid with the same freestream, thermodynamics and Sutherland transport.

Two modes (per model, --model SA|SST):

  python compare.py --model SST
      Your run: compares <MODEL>/OUTPUT/wall.tec with the MOSE and OpenFOAM
      reference solutions for the same setup (for SST, the omega-wall-bc and
      sst-production of <MODEL>/input.ini pick the reference files). Prints
      separation/reattachment and the mean C_f and p_w differences, and writes
      OUTPUT/compare-cf.svg and OUTPUT/compare-pw.svg. Nothing outside OUTPUT/
      is touched.

  python compare.py --model SST --docs
      Regenerates the V&V figures in docs/vv/images from the reference files
      only (no run needed):
        SWBLI-cf-<model>-openfoam.svg     verification: MOSE vs OpenFOAM, skin friction
        SWBLI-pw-<model>-openfoam.svg     verification: MOSE vs OpenFOAM, wall pressure
        SWBLI-cf-<model>-validation.svg   validation: MOSE, OpenFOAM, Wind-US, experiment
        SWBLI-cf-sst-omegaBC-mose.svg     (SST) MOSE, practical vs asymptotic omega BC
        SWBLI-cf-sst-production-mose.svg  (SST) MOSE, incompressible vs compressible production

  --fields  also print the boundary-layer and incident-shock metrics of the V&V
            page; needs OUTPUT/field.tec and an OpenFOAM run (a time directory
            with p, U, rho and constant/polyMesh) copied into <MODEL>/OPENFOAM/.
  --plot    show the figures.

Reference wall solutions (<MODEL>/reference/):
    MOSE-SA.tec, OF-SA.xy                                SA
    MOSE-SST-<omega-wall-bc>-<sst-production>.tec        MOSE SST
        asymptotic-compressible     base case (SST/input.ini as committed)
        asymptotic-incompressible   sst-production = incompressible
        practical-compressible      omega-wall-bc = practical
    OF-SST-<wall>-<production>.xy                        OpenFOAM SST
        asymptotic-compressible     same model as MOSE's base case (system/omegaFaceBC)
        asymptotic-incompressible   same model as MOSE-SST-asymptotic-incompressible
        wallfunction-compressible   OpenFOAM default (kOmegaSST + omegaWallFunction)
        wallfunction-incompressible omegaWallFunction + kOmegaSSTMOSE
    wind-<MODEL>.dat, SU2-<MODEL>.dat, schulein.dat      literature

OpenFOAM walls are raw bottom-wall patch samples (see the README for the
postProcess commands); MOSE walls are OUTPUT/wall.tec files.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Same theme-aware SVG setup as the verify.py scripts: transparent canvas and
# real <text> elements, so the docs CSS can recolor them for the page theme.
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
from ORION import read_TEC   # noqa: E402

# Figures are the V&V documentation figures: write them straight into the docs.
OUT = root.parents[3] / "docs" / "vv" / "images"
LABEL_FS = 17

# Freestream (<MODEL>/input.ini): M=5, p=4000 Pa, T=68.3 K, air
R, GAMMA = 287.0, 1.4
T_INF, P_INF, MACH = 68.3, 4000.0, 5.0
RHO_INF = P_INF / (R * T_INF)
U_INF = MACH * np.sqrt(GAMMA * R * T_INF)
Q_INF = 0.5 * RHO_INF * U_INF ** 2

# Categorical trio (validated CVD-safe in light & dark, dataviz palette slots 1-3):
# MOSE = orange (hero), OpenFOAM = aqua, Wind-US = blue.
MOSE_C, OF_C, WIND_C = "#eb6834", "#1baf7a", "#2a78d6"
# For verification the two codes must solve the *same* SST model: MOSE runs
# sst-production = compressible (OpenFOAM's kOmegaSST production) and OpenFOAM
# applies MOSE's omega wall condition on the wall face (OPENFOAM/system/omegaFaceBC,
# instead of omegaWallFunction, which fixes the wall-adjacent cell).


# ----------------------------------------------------------------------
# OpenFOAM readers (binary case files, no OpenFOAM installation needed)
# ----------------------------------------------------------------------
def _next_list(buf, pos, itemsize, ncomp=1, dtype="<f8"):
    """Read the next '<count>\\n(<binary blob>)' list starting at/after pos."""
    m = re.compile(rb"\n(\d+)\s*\n\(").search(buf, pos)
    if m is None:
        raise ValueError("no binary list found")
    n = int(m.group(1))
    start = m.end()
    nbytes = n * ncomp * itemsize
    a = np.frombuffer(buf, dtype=dtype, count=n * ncomp, offset=start)
    end = start + nbytes
    return (a.reshape(n, ncomp) if ncomp > 1 else a), end


def of_cell_centres(case):
    """Cell centres of an OpenFOAM polyMesh, from points/faces/owner/neighbour.

    For a hex, every vertex is shared by 3 of its faces, so summing the points
    of all the faces of a cell and dividing by the count gives the vertex
    centroid exactly.
    """
    pm = case / "constant" / "polyMesh"

    pts, _ = _next_list((pm / "points").read_bytes(), 0, 8, ncomp=3)

    buf = (pm / "faces").read_bytes()
    off, pos = _next_list(buf, 0, 4, dtype="<i4")          # compact-list offsets
    fpts, _ = _next_list(buf, pos, 4, dtype="<i4")         # face point indices

    owner, _ = _next_list((pm / "owner").read_bytes(), 0, 4, dtype="<i4")
    neigh, _ = _next_list((pm / "neighbour").read_bytes(), 0, 4, dtype="<i4")

    face_sum = np.add.reduceat(pts[fpts], off[:-1], axis=0)
    face_cnt = np.diff(off).astype(float)

    ncells = int(max(owner.max(), neigh.max())) + 1
    csum = np.zeros((ncells, 3))
    ccnt = np.zeros(ncells)
    np.add.at(csum, owner, face_sum)
    np.add.at(ccnt, owner, face_cnt)
    np.add.at(csum, neigh, face_sum[: neigh.size])
    np.add.at(ccnt, neigh, face_cnt[: neigh.size])
    return csum / ccnt[:, None]


def of_field(path):
    """internalField of a binary volScalarField."""
    buf = path.read_bytes()
    m = re.compile(rb"internalField\s+nonuniform\s+List<scalar>").search(buf)
    if m is None:
        m = re.compile(rb"internalField\s+uniform\s+([-\d.eE+]+)").search(buf)
        raise ValueError(f"{path.name}: uniform field ({float(m.group(1))})")
    a, _ = _next_list(buf, m.end(), 8)
    return a


def of_vfield(path):
    """internalField of a binary volVectorField, shape (ncells, 3)."""
    buf = path.read_bytes()
    m = re.compile(rb"internalField\s+nonuniform\s+List<vector>").search(buf)
    if m is None:
        raise ValueError(f"{path.name}: not a nonuniform vector field")
    a, _ = _next_list(buf, m.end(), 8, ncomp=3)
    return a


def of_latest_time(case):
    times = [d for d in case.iterdir()
             if d.is_dir() and re.fullmatch(r"\d+", d.name) and d.name != "0"]
    if not times:
        raise FileNotFoundError(f"no time directory in {case} — is the run copied over?")
    return max(times, key=lambda d: int(d.name))


def of_wall(path):
    """(x, Cf, p) on the bottom wall from a raw patchSurface sample."""
    d = np.loadtxt(path)
    o = np.argsort(d[:, 0])
    # wallShearStress is the traction on the fluid, hence negative under an
    # attached BL: flip the sign to match the MOSE tauX convention.
    return d[o, 0], -d[o, 3] / Q_INF, d[o, 6]


def ref_xy(path):
    """(x, Cf) from a literature reference file, or (None, None) if absent."""
    if not path.exists():
        return None, None
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


# ----------------------------------------------------------------------
# MOSE readers
# ----------------------------------------------------------------------
def mose_wall(path):
    """(x, Cf, p) from wall.tec. Vars are ['y+','tauX','tauY','tauZ','pw',...]."""
    xw, _, _, vw, _ = read_TEC(str(path))
    xn = xw[0][:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    return xc, vw[0][1][:, 0, 0] / Q_INF, vw[0][4][:, 0, 0]


def mose_fields(path):
    """Per block: cell-centre x, y and the cell-centred rho, u, p from field.tec.

    Coordinates are nodal, variables are cell-centred: average the four nodes
    of each cell to get the cell centres.
    """
    xb, yb, _, vb, names = read_TEC(str(path))
    irho = [k - 3 for k, n in enumerate(names) if n.startswith("rho(")]
    iu, ip = (names.index(n) - 3 for n in ("u", "p"))
    out = []
    for x, y, v in zip(xb, yb, vb):
        xn, yn = x[:, :, 0], y[:, :, 0]
        out.append({
            "x": 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:]),
            "y": 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:]),
            "rho": sum(v[k][:, :, 0] for k in irho),
            "u": v[iu][:, :, 0],
            "p": v[ip][:, :, 0],
        })
    return out


# ----------------------------------------------------------------------
# Field metrics (--fields): incoming boundary layer and incident shock
# ----------------------------------------------------------------------
def bl_metrics(y, u, rho):
    """delta99, compressible theta and incompressible H of one wall-normal profile.

    The edge is the velocity maximum of the profile (the leading-edge shock keeps
    it below U_inf); integrals run from the wall to delta99.
    """
    ke = int(np.argmax(u))
    ue, re_ = u[ke], rho[ke]
    j = int(np.argmax(u >= 0.99 * ue))
    d99 = np.interp(0.99 * ue, u[: j + 1], y[: j + 1])
    m = y <= d99
    f = u[m] / ue
    theta = np.trapz(rho[m] / re_ * f * (1 - f), y[m])
    h_i = np.trapz(1 - f, y[m]) / np.trapz(f * (1 - f), y[m])
    return d99, theta, h_i


def shock_metrics(x, p, lo=0.20, hi=0.33):
    """Incident-shock position (mid-pressure), strength p2/p1 and 10-90% width in cells.

    Levels are located by walking out from the steepest point to the first
    crossing, so a post-shock over/undershoot (central-upwind flux) is not
    counted as shock thickness.
    """
    m = (x > lo) & (x < hi)
    x, p = x[m], p[m]
    k = int(np.argmax(np.gradient(p, x)))
    p1 = np.median(p[(x > x[k] - 0.02) & (x < x[k] - 0.008)])
    p2 = np.median(p[(x > x[k] + 0.008) & (x < x[k] + 0.02)])

    def crossing(f):
        lvl = p1 + f * (p2 - p1)
        i = k
        if p[k] >= lvl:
            while p[i - 1] >= lvl:
                i -= 1
            i -= 1
        else:
            while p[i + 1] < lvl:
                i += 1
        return x[i] + (lvl - p[i]) * (x[i + 1] - x[i]) / (p[i + 1] - p[i])

    dx = np.mean(np.diff(x[k - 3: k + 4]))
    return crossing(0.5), p2 / p1, (crossing(0.9) - crossing(0.1)) / dx


def along_y(blk, f, yv):
    """(x, f) along the horizontal line y = yv through a structured block."""
    xs, fs = [], []
    for i in range(blk["x"].shape[0]):
        yy = blk["y"][i, :]
        if not yy[0] <= yv <= yy[-1]:
            continue
        j = min(max(int(np.searchsorted(yy, yv)), 1), yy.size - 1)
        w = (yv - yy[j - 1]) / (yy[j] - yy[j - 1])
        xs.append((1 - w) * blk["x"][i, j - 1] + w * blk["x"][i, j])
        fs.append((1 - w) * f[i, j - 1] + w * f[i, j])
    return np.array(xs), np.array(fs)


def field_report(case):
    """Print the BL and shock metrics of MOSE and OpenFOAM on the shared grid.

    Needs OUTPUT/field.tec and an OpenFOAM time directory with p, U and rho
    (copied over from the run, with constant/polyMesh).
    """
    from scipy.spatial import cKDTree

    blocks = mose_fields(case / "OUTPUT" / "field.tec")
    ofc = case / "OPENFOAM"
    t = of_latest_time(ofc)
    tree = cKDTree(of_cell_centres(ofc)[:, :2])
    of_p, of_rho, of_u = of_field(t / "p"), of_field(t / "rho"), of_vfield(t / "U")[:, 0]

    def mapped(blk):
        """OpenFOAM cell of every MOSE cell: same grid, matched by centre."""
        d, idx = tree.query(np.c_[blk["x"].ravel(), blk["y"].ravel()])
        if d.max() > 1e-6:
            raise RuntimeError(f"grids differ: cell centres up to {d.max():.2e} m apart")
        return idx.reshape(blk["x"].shape)

    print(f"\nField metrics (OpenFOAM time {t.name}):")
    b = blocks[0]
    o = mapped(b)
    i = int(np.argmin(np.abs(b["x"][:, 0] - 0.25)))
    mm = bl_metrics(b["y"][i, :], b["u"][i, :], b["rho"][i, :])
    mo = bl_metrics(b["y"][i, :], of_u[o[i, :]], of_rho[o[i, :]])
    for lbl, a, c, s in (("delta99 [mm]", mm[0], mo[0], 1e3), ("theta [mm]", mm[1], mo[1], 1e3),
                         ("H_i", mm[2], mo[2], 1.0)):
        print(f"  x = 0.25  {lbl:13s} MOSE {s * a:8.4f}  OpenFOAM {s * c:8.4f}  ({100 * (c / a - 1):+5.1f}%)")
    for yv in (0.020, 0.030):
        b = next(bb for bb in blocks if bb["y"].min() < yv < bb["y"].max())
        o = mapped(b)
        sm = shock_metrics(*along_y(b, b["p"], yv))
        so = shock_metrics(*along_y(b, of_p[o], yv))
        print(f"  y = {yv:.3f} incident shock: x MOSE {sm[0]:.4f}  OpenFOAM {so[0]:.4f}"
              f" | p2/p1 {sm[1]:.3f} {so[1]:.3f} | 10-90% width {sm[2]:.1f} {so[2]:.1f} cells")


# ----------------------------------------------------------------------
def sep_reatt(x, cf, lo=0.30, hi=0.42):
    m = (x > lo) & (x < hi)
    xs, cs = np.asarray(x)[m], np.asarray(cf)[m]
    n = np.where(cs < 0)[0]
    return (xs[n[0]], xs[n[-1]]) if n.size else (None, None)


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.tick_params(labelsize=LABEL_FS - 2)
    ax.legend(loc="best", fontsize=LABEL_FS - 2)
    ax.grid(True, alpha=0.3)


def save(fig, name):
    """Write a V&V figure into docs/vv/images."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=150, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote docs/vv/images/{name}")


def save_local(fig, path):
    """Write a user-comparison figure (never into the docs)."""
    fig.savefig(path, dpi=150, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {path.relative_to(root)}")


# ----------------------------------------------------------------------
# Reference selection and metrics
# ----------------------------------------------------------------------
def sst_setup(case):
    """(omega-wall-bc, sst-production) of <case>/input.ini, with MOSE's defaults."""
    import configparser
    ini = configparser.ConfigParser(strict=False, inline_comment_prefixes=("!", "#", ";"))
    ini.read(case / "input.ini")
    turb = ini["MOSE-Turbulence"] if ini.has_section("MOSE-Turbulence") else {}
    return (turb.get("omega-wall-bc", "practical").strip().lower(),
            turb.get("sst-production", "incompressible").strip().lower())


def references(case, M, tag=None):
    """Paths of the MOSE and OpenFOAM reference walls (None where there is none).

    SA has one of each; SST has one per (omega-wall-bc, sst-production) tag.
    """
    ref = case / "reference"
    if M == "SA":
        mose, of = ref / "MOSE-SA.tec", ref / "OF-SA.xy"
    else:
        mose, of = ref / f"MOSE-SST-{tag}.tec", ref / f"OF-SST-{tag}.xy"
    return (mose if mose.exists() else None), (of if of.exists() else None)


def report(curves, base):
    """Print separation/reattachment of each curve and its differences from `base`."""
    print("\nSeparation / reattachment (Cf = 0):")
    for name, (x, cf, _) in curves.items():
        s, r = sep_reatt(x, cf)
        L = f"{1e3 * (r - s):.2f} mm" if s else "attached"
        print(f"  {name:22s} x_sep = {s:.5f}  x_reatt = {r:.5f}  L = {L}" if s
              else f"  {name:22s} attached")
    xb, cfb, pb = curves[base]
    for name, (x, cf, p) in curves.items():
        if name == base:
            continue
        print(f"\nMean |difference|, {name} vs {base}:")
        for lo, hi, lbl in ((0.10, 0.30, "upstream   "),
                            (0.35, 0.40, "post-reatt."),
                            (0.45, 0.52, "downstream ")):
            xs = np.linspace(lo, hi, 300)
            cb, pbi = np.interp(xs, xb, cfb), np.interp(xs, xb, pb)
            ec = 100 * np.mean(np.abs(np.interp(xs, x, cf) - cb) / np.abs(cb))
            ep = 100 * np.mean(np.abs(np.interp(xs, x, p) - pbi) / pbi)
            print(f"  {lbl} ({lo:.2f}-{hi:.2f}):  Cf {ec:5.1f}%   p_w {ep:5.2f}%")


def cf_axes(ax, xlim=(0.32, 0.41)):
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.002, 0.0075)
    style(ax, r"$x$  [m]", r"$C_f$")


# ----------------------------------------------------------------------
# Modes
# ----------------------------------------------------------------------
def compare_run(case, M):
    """Your run (OUTPUT/wall.tec) against the references for the same setup."""
    tag = None
    if M == "SST":
        bc, prod = sst_setup(case)
        tag = f"{bc}-{prod}"
        print(f"Setup (input.ini): omega-wall-bc = {bc}, sst-production = {prod}")
    mose_ref, of_ref = references(case, M, tag)
    print(f"MOSE reference:     {mose_ref.relative_to(root) if mose_ref else 'none for this setup'}")
    print(f"OpenFOAM reference: {of_ref.relative_to(root) if of_ref else 'none for this setup'}")

    wall = case / "OUTPUT" / "wall.tec"
    if not wall.exists():
        print(f"\nNo {wall.relative_to(root)}: run MOSE first (./MOSE.sh solve).")
        print("Use --docs to regenerate the V&V figures from the references alone.")
        return

    curves = {"MOSE (OUTPUT)": mose_wall(wall)}
    if mose_ref:
        curves["MOSE reference"] = mose_wall(mose_ref)
    if of_ref:
        curves["OpenFOAM reference"] = of_wall(of_ref)
    report(curves, "MOSE (OUTPUT)")

    style_of = {"MOSE (OUTPUT)": dict(color=MOSE_C, lw=4.0),
                "MOSE reference": dict(color="0.35", ls="--", lw=2.0),
                "OpenFOAM reference": dict(color=OF_C, lw=2.0)}
    fig, ax = plt.subplots(figsize=(10, 6))
    xe, cfe = ref_xy(case / "reference" / "schulein.dat")
    if xe is not None:
        ax.plot(xe, cfe, "ko", ms=6, label="Schulein (exp)")
    for name, (x, cf, _) in curves.items():
        ax.plot(x, cf, label=name, **style_of[name])
    cf_axes(ax, (0.32, 0.45))
    save_local(fig, case / "OUTPUT" / "compare-cf.svg")

    fig, ax = plt.subplots(figsize=(10, 6))
    for name, (x, _, p) in curves.items():
        ax.plot(x, p / P_INF, label=name, **style_of[name])
    ax.set_xlim(0.05, 0.52)
    style(ax, r"$x$  [m]", r"$p_w / p_\infty$")
    save_local(fig, case / "OUTPUT" / "compare-pw.svg")


def docs_figures(case, M):
    """The V&V page figures, from the reference files only."""
    m = M.lower()
    ref = case / "reference"
    mose_ref, of_ref = references(case, M, None if M == "SA" else "asymptotic-compressible")
    xm, cfm, pwm = mose_wall(mose_ref)
    xo, cfo, pwo = of_wall(of_ref)
    print(f"MOSE:     {mose_ref.relative_to(root)}\nOpenFOAM: {of_ref.relative_to(root)}")

    curves = {"MOSE": (xm, cfm, pwm), "OpenFOAM": (xo, cfo, pwo)}
    report(curves, "MOSE")

    # --- verification: MOSE vs OpenFOAM, same grid and model -----------
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xm, cfm, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xo, cfo, OF_C, lw=2.5, label="OpenFOAM")
    cf_axes(ax)
    save(fig, f"SWBLI-cf-{m}-openfoam.svg")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xm, pwm / P_INF, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xo, pwo / P_INF, OF_C, lw=2.5, label="OpenFOAM")
    ax.set_xlim(0.05, 0.52)
    style(ax, r"$x$  [m]", r"$p_w / p_\infty$")
    save(fig, f"SWBLI-pw-{m}-openfoam.svg")

    # --- validation: MOSE and OpenFOAM vs the literature ----------------
    fig, ax = plt.subplots(figsize=(10, 6))
    xe, cfe = ref_xy(ref / "schulein.dat")
    xw, cfw = ref_xy(ref / f"wind-{M}.dat")
    if xe is not None:
        ax.plot(xe, cfe, "ko", ms=6, label="Schulein (exp)")
    if xw is not None:
        ax.plot(xw, cfw, WIND_C, ls="-.", lw=2.5, label="Wind-US")
    ax.plot(xm, cfm, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xo, cfo, OF_C, lw=2.0, label="OpenFOAM")
    cf_axes(ax)
    save(fig, f"SWBLI-cf-{m}-validation.svg")

    if M != "SST":
        return

    # --- MOSE sensitivity studies: same grid, one input line changed ------
    # Same colour (MOSE) for both curves; the line style carries the option.
    for other, name, labels in (
            ("practical-compressible", "SWBLI-cf-sst-omegaBC-mose.svg",
             (r"MOSE, practical $\omega$ BC  ($800\,\nu/y^2$)",
              r"MOSE, asymptotic $\omega$ BC  ($80\,\nu/y^2$)")),
            ("asymptotic-incompressible", "SWBLI-cf-sst-production-mose.svg",
             (r"MOSE, incompressible production  ($\mu_t S^2$)",
              r"MOSE, compressible production  ($\tau_{ij}\,\partial u_i/\partial x_j$)"))):
        xs_, cfs_, _ = mose_wall(ref / f"MOSE-SST-{other}.tec")
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(xs_, cfs_, MOSE_C, ls="--", lw=2.5, label=labels[0])
        ax.plot(xm, cfm, MOSE_C, lw=3.5, label=labels[1])
        cf_axes(ax, (0.32, 0.45))
        save(fig, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", choices=("SA", "SST"), default="SA")
    ap.add_argument("--docs", action="store_true",
                    help="regenerate the V&V figures in docs/vv/images from the references")
    ap.add_argument("--fields", action="store_true",
                    help="also print BL and shock metrics (needs OUTPUT/field.tec "
                         "and the OpenFOAM run copied into OPENFOAM/)")
    ap.add_argument("--plot", action="store_true", help="show the figures")
    args = ap.parse_args()

    M = args.model
    case = root / M
    print(f"Model: {M}")
    print(f"Freestream: U={U_INF:.2f} m/s  rho={RHO_INF:.4f} kg/m3  q={Q_INF:.0f} Pa")

    if args.docs:
        docs_figures(case, M)
    else:
        compare_run(case, M)

    if args.fields:
        field_report(case)

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
