#!/usr/bin/env python3
"""HyShot II injector mesh -- three independent, overlapping structured grids
coupled by CHIMERA overlap (no conformal interfaces, no gmsh).

The domain is split into three regions that overlap at their seams; ATLAS BCB
computes the donor/overlap (type-102) connectivity from the overlaps:

    x = 360 ...... 408 (injector) ...... 650 (X_RAMP) ...... 760   (streamwise)
                                                                _.-'
  +---------------------------------------------------+---_.-''      ceiling
  |  UPSTREAM      |            CHAMBER  (background)  |_.-'   12 deg
  |  coarse box    |   homogeneous  dx = dy ,  dz graded (BL)     |  y = 0..w
  |  dx_up         |          +----------------+                  |  z = 0..H(x)
  +----------------+----------| INJECTION      |------------------+
        overlap               |  O-grid + pipe |   overset patch
                              +----------------+
Past X_RAMP = 650 the duct is no longer constant-section: the ceiling turns up
by ramp_deg and the channel height becomes H(x) (bottom wall stays flat).  Only
the downstream chamber block gets there; X_OUT <= 650 gives the plain tube.

Coordinates:  x = streamwise, y = spanwise (y=0 symmetry through the injector
axis, y = w half pitch), z = wall-normal (z=0 bottom wall, z=H ceiling).

The injector is the upper half (y>=0) of a circle of radius R centred at
(x_inj, 0) on the wall; the pipe is the same half-disk extruded below the wall.

FULL CONTROL is exposed as the knobs in the "PARAMETERS" block below:
  dx, dy         chamber streamwise / spanwise spacing (homogeneous)
  dx_up          upstream coarse spacing
  dz1, rbl, nbl  wall boundary layer (first height, growth ratio, #cells)
  dz_core        wall-normal spacing of the uniform core (between the BLs)
  refine_inj     injection grid fineness vs. the chamber (d_inj = dx/refine_inj)
  rim_refine     extra refinement AT the injector rim (1 = none, 4 = 4x finer)
  core_size      size of the O-grid core square = size of its cells (bore centre)
  X_RAMP, ramp_deg, dx_exit, exit_stretch   exit expansion (ramp past X_RAMP)
  x_up_end, x_ch_start, ob2, z_inj, Lpipe   region extents / overlaps

Usage:
    python3 hyshot_mesh.py                 -> mesh.p3d  (PLOT3D, for MOSE)
    python3 hyshot_mesh.py mesh.tec        -> Tecplot   (for ATLAS BCB chimera)
    python3 hyshot_mesh.py mesh.p3d mesh.tec   -> both
"""

import sys
from math import cos, sin, tan, pi, sqrt, log, ceil, radians

# ============================ PARAMETERS ==================================== #
fs     = 0.001
# ---- domain --------------------------------------------------------------- #
X_IN   = 360.0     # inlet plane
X_OUT  = 800.0     # outlet plane (clamped to X_MAX, see the expansion below)
X_INJ  = 408.0     # injector axis
W      = 9.375     # spanwise half-pitch (y = 0..W)
H      = 9.8       # channel height (z = 0..H) of the CONSTANT-SECTION duct
R      = 1.0       # injector radius (d/2)
LPIPE  = 2.0       # injector pipe length below the wall

# ---- homogeneous chamber (background) ------------------------------------- #
dx     = 0.25      # streamwise spacing        (== dy for a truly square grid)
dy     = 0.25      # spanwise spacing
dz_core = 0.25     # wall-normal spacing of the uniform core

# ---- exit expansion (constant-section duct -> nozzle) --------------------- #
# The duct is constant-section only up to X_RAMP.  Past it the CEILING turns up
# by ramp_deg while the bottom wall stays flat:
#       H(x) = H + max(0, x - X_RAMP) * tan(ramp_deg)
# so the channel is H tall at X_RAMP and H + (X_MAX-X_RAMP)*tan(ramp_deg) at the
# very end.  X_OUT <= X_RAMP -> nothing changes, it is the plain tube as before.
# ONLY the downstream chamber block reaches past X_RAMP (everything else ends
# near the injector), so the shared chamber node lattice stays conformal.
#   wall-normal: the wall BL keeps dz1 and the ceiling BL rides up rigidly with
#                the ceiling; the uniform core between them is STRETCHED.  The
#                cell count is fixed by the shared lattice, so the core spacing
#                grows with the duct -- the run prints the exit value.
#   streamwise:  a node is forced onto X_RAMP so the kink stays sharp, and the
#                expansion GROWS OUT OF the chamber spacing instead of jumping to
#                its own: cells start at dx right after X_RAMP, grow by
#                exit_stretch per cell, stop growing once they reach dx_exit and
#                stay uniform to the outlet (the graded() law again).  dx is
#                therefore CONTINUOUS across the kink and dx_exit is only the CAP
#                on how coarse the nozzle may ever get -- not a step change.
#                exit_stretch <= 1 -> old behaviour, uniform dx_exit throughout.
#                The run prints the cell count and the dx range in the nozzle.
X_RAMP   = 650.0   # start of the expansion
X_MAX    = 760.0   # end of the modelled hardware; X_OUT is clamped to it
ramp_deg = 12.0    # ceiling inclination past X_RAMP, degrees (0 = flat)
dx_exit  = 2.0     # COARSEST streamwise cell in the expansion (the cap)
exit_stretch = 1.2 # growth per cell from dx at X_RAMP up to dx_exit

# ---- upstream block (block 1) --------------------------------------------- #
# x is geometrically STRETCHED: the finest cell (up_dxmin) sits at the injection
# side (x = CX-ob2, to match the finer chamber/injection) and grows toward the
# inlet.  y is uniform (dy_up); z keeps the shared wall BL (dz1/rbl/nbl) with its
# own core spacing (dz_up).
up_dxmin   = 0.25   # finest streamwise cell, at the injection side
up_stretch = 1.02   # streamwise growth ratio toward the inlet (1.0 = uniform)
dy_up      = 0.40   # spanwise spacing
dz_up      = 0.40   # wall-normal core spacing (wall BL retained)

# ---- wall boundary layer (shared by every region) ------------------------- #
dz1    = 0.001      # first cell height at the wall
rbl    = 1.20      # geometric growth ratio
nbl    = 28        # number of BL cells
CEIL_BL = True     # mirror the BL at the ceiling (z = H) as well

# ---- injection overset ---------------------------------------------------- #
ob2   = W / 2 * (sqrt(2) / 2) * 1.4   # half-size of the injection footprint box
z_inj = 4.0        # top of the injection patch (must be < H, overlaps chamber)

# O-GRID CORE SQUARE (block S3, the middle of the injector bore).  Its cell
# COUNT is fixed at 2N x N by the azimuthal resolution N of the arc, so its cell
# size is exactly  hl/N : ENLARGE THE SQUARE AND ITS CELLS GROW WITH IT, one for
# one.  This is the knob for the small cells in the middle of the bore (they are
# what limits the explicit time step there).  Footprint on the wall: 2*hl x hl,
# i.e. area 2*hl^2 -- doubling core_size quadruples the area.
#   core_size  half-width of the square as a fraction of the injector radius R.
#              What limits it is the PETALS: they take up what is left, and the
#              square's corners approach the circle faster than its edges, so the
#              petal thins much faster at 45 deg than on the axis --
#                   0.25 -> 0.75 .. 0.65      0.5 -> 0.50 .. 0.29
#                   0.6  -> 0.40 .. 0.15      0.7 -> 0.30 .. 0.01  (pinched off)
#              Cell ANGLES stay ~48 deg throughout, so it is the thinning, not
#              skew, that sets the cap core_max (0.707 = corners on the circle).
#              At 0.6 the core cells are within ~25 % of the arc cells at the rim
#              -- about as uniform as an O-grid core gets.  The run prints the
#              core cell size, its ratio to the arc cell, and the petal thickness
#              range, so you can dial it by eye.
core_size = 0.50   # 0.25 = the original square
core_max  = 0.60   # cap; > 0.707 is geometrically impossible
hl = min(core_size, core_max) * R     # half-size of the O-grid core square

# In-plane (i,j) discretization of the O-grid cross-section.
#   refine_inj  is the master fineness; every count derives from it UNLESS you
#               override that direction with an explicit cell count below.
#   n_arc  circumferential resolution AROUND the O (cells per 45-deg quadrant;
#          the 90-deg top arc gets 2*n_arc, the core gets n_arc x 2*n_arc)
#   n_rad  radial cells from the circle out to the outer box (ring thickness)
#   n_pet  radial cells from the circle in to the core square (petal thickness)
refine_inj = 4.0
n_arc = None       # None -> round((R*pi/4)/(dx/refine_inj))
n_rad = None       # None -> auto (follows rim_refine, see below)
n_pet = None       # None -> auto (follows rim_refine, see below)

# RADIAL REFINEMENT AT THE INJECTOR RIM (the circle r = R, where the pipe/petal
# blocks meet the outer ring blocks -- the shear layer of the jet).
#   rim_refine   plain refinement FACTOR, like refine_inj: the radial cell that
#                touches the rim is  d_inj / rim_refine.
#                   1 = uniform (no clustering)     2 = twice as fine at the rim
#                   4 = four times finer            ... no units, no thresholds.
#   rim_growth   how fast cells grow away from the rim (1.2 = +20 % per cell).
#                They stop growing at d_inj and stay uniform after that, so
#                refining the rim ADDS cells instead of coarsening the far side
#                (same monotone law as graded() uses in z).  n_rad / n_pet are
#                derived from this unless you set them explicitly.
#   rim_refine_in / rim_refine_out   refine only one side: inside the circle
#                (petals/pipe bore) / outside it (rings).  None = use rim_refine.
# The outer side is limited by the top-ring (S7) trapezoid: too strong a value
# folds it.  That is now handled automatically -- the outer factor is backed off
# toward 1 until the blocks are valid, and the value actually used is printed.
rim_refine     = 8.0
rim_growth     = 1.2
rim_refine_in  = None
rim_refine_out = None

# Wall-normal (k) discretization of the injection column, INDEPENDENT of the
# surrounding chamber.  Unlike the chamber's fixed nbl boundary layer, the
# injection column is MONOTONE (see graded()): cells grow from dz1 at the wall,
# capped at dz_core_inj, then uniform -- so refining the core (smaller
# dz_core_inj) automatically resolves the BL better, with no fine-coarse-fine
# bump.  nbl does not apply here; the near-wall count follows from dz1 -> dz_core_inj.
dz_core_inj = 0.25 * dz_core   # smaller -> finer wall-normal AND better-resolved BL

# ---- region seams -------------------------------------------------------- #
# The chamber is carved into THREE blocks around the injection "notch" so the
# injection grid is NOT embedded in a background volume (no pure volume
# overlap): chimera only couples the differing discretizations across the
# injection's outer surfaces.  The injection footprint is the box
#   x in [CX-ob2, CX+ob2],  y in [0, ob2],  z in [0, z_inj].
# The three chamber blocks share one node lattice so they connect conformally
# to each other; they abut the injection box on its outer faces:
#   cham DOWNSTREAM  x in [CX+ob2, X_OUT], full y, full z
#   cham ABOVE       x in [CX-ob2, CX+ob2], y in [0, ob2], z in [z_inj, H]
#   cham LATERAL     x in [CX-ob2, CX+ob2], y in [ob2, W], full z
front_overlap = 0.0   # upstream extends this far past the injection front (x=CX-ob2)

# ---- multigrid compatibility --------------------------------------------- #
# If the solver coarsens the grid (MG_LEVEL-1) times (2:1 per level), every
# block's CELL count in each direction must be divisible by 2^(MG_LEVEL-1).
# With MG_LEVEL > 1 the cell counts derived from dx/dy/dz.../n_* (and each
# internal segment: BL+core columns, chamber splits at xR/ob2/z_inj, the pipe
# part of the merged hole blocks) are rounded UP to satisfy this -- spacings
# shift slightly, split planes still survive coarsening.  MG_LEVEL = 1 = off.
MG_LEVEL = 3
# =========================================================================== #

SR2 = sqrt(2) / 2
OB  = R * SR2                       # 45-degree point of the circle
CX  = X_INJ
TAN_RAMP = tan(radians(ramp_deg))   # ceiling slope past X_RAMP


def z_top(x):
    """Local ceiling height: H through the constant-section duct, then rising at
    ramp_deg past X_RAMP (the bottom wall stays flat at z = 0)."""
    return H + max(0.0, x - X_RAMP) * TAN_RAMP


# --------------------------- 1-D node distributions ------------------------ #
def mg(n):
    """Round a cell count UP to a multiple of the multigrid factor 2^(MG_LEVEL-1)
    (so the block can be coarsened MG_LEVEL-1 times). MG_LEVEL=1 -> unchanged."""
    n = max(1, n)
    f = 2 ** (MG_LEVEL - 1)
    return n if f <= 1 else ((n + f - 1) // f) * f


def linspace(a, b, n):
    """n+1 evenly spaced nodes from a to b (n cells)."""
    return [a + (b - a) * i / n for i in range(n + 1)]


def uniform(a, b, dtarget):
    """Node coords from a to b at ~dtarget spacing; cell count MG-rounded."""
    return linspace(a, b, mg(round(abs(b - a) / dtarget)))


def bl_heights():
    """Cumulative BL node offsets [0, dz1, ...] -> nbl+1 nodes, top = hbl."""
    z, acc = [0.0], 0.0
    for i in range(nbl):
        acc += dz1 * rbl ** i
        z.append(acc)
    return z            # z[-1] == hbl


def z_column(z0, z1, wall=True, ceil=True, dcore=None):
    """Wall-normal node coords from z0 to z1: geometric BL at z0 (if wall),
    uniform core (spacing dcore, default dz_core), mirrored BL at z1 (if ceil)."""
    bl = bl_heights()
    hbl = bl[-1]
    dc = dcore if dcore else dz_core
    lo = z0 + hbl if wall else z0
    hi = z1 - hbl if ceil else z1
    base = (nbl if wall else 0) + (nbl if ceil else 0)   # fixed BL cells
    ncore = mg(base + max(1, round((hi - lo) / dc))) - base   # so total is MG-divisible
    nodes = [z0 + h for h in bl] if wall else [z0]
    nodes += linspace(lo, hi, ncore)[1:]     # drop duplicate first node
    if ceil:
        for h in reversed(bl[:-1]):          # mirror BL, skip duplicate hbl top
            nodes.append(z1 - h)
    return nodes


def graded(z0, z1, d_wall, ratio, d_core, fine='low'):
    """Monotone nodes z0..z1 with NO coarse overshoot: cells start at d_wall on
    the fine end (fine='low' -> z0, 'high' -> z1), grow geometrically by `ratio`,
    and are CAPPED at d_core (uniform thereafter). Unlike the fixed nbl boundary
    layer, a refined core (small d_core) automatically resolves the BL instead of
    leaving a coarse bump between the fine core and the fine wall.  Used for the
    wall-normal columns and, in x, to grow the nozzle out of the chamber spacing."""
    L = z1 - z0
    sizes, s, acc = [], d_wall, 0.0
    while s < d_core and acc + s < L:        # geometric ramp, only while growing
        sizes.append(s)
        acc += s
        s *= ratio
    rem = L - acc                            # single uniform batch fills the rest
    if rem > 1e-9:
        ncore = max(1, mg(len(sizes) + max(1, round(rem / d_core))) - len(sizes))
        sizes += [rem / ncore] * ncore       # ncore chosen so total cells are MG-divisible
    off = [0.0]
    for cs in sizes:
        off.append(off[-1] + cs)
    off[-1] = L
    if fine == 'low':
        return [z0 + p for p in off]
    return [z1 - p for p in reversed(off)]   # fine end at z1


def stretch(a, b, d_fine, strength, fine='high'):
    """Nodes a..b geometrically STRETCHED: the cell at the `fine` end has size
    d_fine and each next cell grows toward the other end.  The cell COUNT auto-
    fits the length and is MG-rounded; d_fine is held exact and the growth ratio
    is solved to fit (so it lands slightly below `strength`).  strength<=1 or an
    unusable d_fine -> uniform at d_fine."""
    L = b - a
    if strength <= 1.0 + 1e-9 or d_fine <= 0.0 or d_fine >= L:
        return uniform(a, b, d_fine if d_fine > 0.0 else L)
    n = mg(max(1, ceil(log(1.0 + L * (strength - 1.0) / d_fine) / log(strength))))
    lo, hi = 1.0 + 1e-12, 16.0               # solve growth r so n cells fill L exactly
    for _ in range(200):
        r = 0.5 * (lo + hi)
        try:
            s = d_fine * (r ** n - 1.0) / (r - 1.0)
        except OverflowError:
            s = float('inf')
        if s > L:
            hi = r
        else:
            lo = r
    r = 0.5 * (lo + hi)
    off = [0.0]
    for i in range(n):
        off.append(off[-1] + d_fine * r ** i)
    off[-1] = L                              # absorb the residual into the coarse end
    return [a + p for p in off] if fine == 'low' else [b - p for p in reversed(off)]


def extrude(face_x, face_y, zs):
    """Extrude a 2-D face (lists indexed [j][i]) along z -> (X,Y,Z) i-fastest."""
    nj = len(face_x)
    ni = len(face_x[0])
    X, Y, Z = [], [], []
    for z in zs:
        for j in range(nj):
            for i in range(ni):
                X.append(face_x[j][i])
                Y.append(face_y[j][i])
                Z.append(z)
    return (ni, nj, len(zs)), X, Y, Z


# ------------------------------- box region -------------------------------- #
def block_from_nodes(xs, ys, zs):
    """Rectangular block from explicit x, y, z node lists."""
    fx = [[x for x in xs] for _ in ys]
    fy = [[y for _ in xs] for y in ys]
    return extrude(fx, fy, zs)


def block_ramped(xs, ys, zs, top):
    """Same block, but with the ceiling following top(x) instead of being flat.
    The column zs (built for the constant-section height zs[-1]) is warped per x:
    the wall BL is untouched, the ceiling BL is translated RIGIDLY so it keeps
    dz1 at the moving wall, and only the uniform core between them stretches.
    top(x) == zs[-1] everywhere -> identical to block_from_nodes()."""
    hbl = bl_heights()[-1]
    lo = zs[0] + hbl                                # top of the wall BL (fixed)
    hi = zs[-1] - hbl if CEIL_BL else zs[-1]        # foot of the ceiling BL (rides up)
    rise = [max(0.0, top(x) - zs[-1]) for x in xs]  # how far the ceiling has risen
    X, Y, Z = [], [], []
    for z in zs:
        # fraction of the rise this node follows: 0 in the wall BL, 1 from the
        # foot of the ceiling BL upward, linear across the core
        s = min(1.0, max(0.0, (z - lo) / (hi - lo))) if hi > lo else 1.0
        for y in ys:
            for i, x in enumerate(xs):
                X.append(x)
                Y.append(y)
                Z.append(z + s * rise[i])
    return (len(xs), len(ys), len(zs)), X, Y, Z


# --------------------------- O-grid (Coons/TFI) ---------------------------- #
def ramp_count(length, d_fine, ratio, d_coarse):
    """Cells needed to span `length` under the ramp law of ramp_fracs -- i.e.
    how many cells the requested rim refinement costs."""
    n, s, acc = 0, d_fine, 0.0
    while s < d_coarse and acc + s < length:
        acc += s
        n += 1
        s *= ratio
    if length - acc > 1e-9:
        n += max(1, round((length - acc) / d_coarse))
    return max(2, n)


def ramp_fracs(n, length, d_fine, ratio, d_coarse):
    """n+1 monotone fractions 0..1 for a direction of extent `length`: the first
    cell is `d_fine`, cells grow by `ratio` until they reach `d_coarse` and are
    uniform after that (graded()'s law, as FRACTIONS so the same profile can be
    reused on the shorter edges of the same block family).  d_fine <= 0, or too
    coarse to cluster (>= length/n), -> uniform."""
    if d_fine <= 0.0 or d_fine >= length / n:
        return [i / n for i in range(n + 1)]
    sizes, s, acc = [], d_fine, 0.0
    while len(sizes) < n - 1 and s < d_coarse and acc + s < length:
        sizes.append(s)
        acc += s
        s *= ratio
    # the leftover cells split the rest uniformly; give back ramp cells until
    # that uniform size is no smaller than the last ramp cell (stays monotone)
    while sizes and (length - acc) / (n - len(sizes)) < sizes[-1]:
        acc -= sizes.pop()
    sizes += [(length - acc) / (n - len(sizes))] * (n - len(sizes))
    tot = sum(sizes)
    f = [0.0]
    for s in sizes:
        f.append(f[-1] + s / tot)
    f[-1] = 1.0
    return f


def line(p, q, fr, cluster='start'):
    """Points from p to q at the monotone fractions `fr` (0..1); an int means
    that many uniform cells.  cluster='end' mirrors the profile so its fine end
    lands at q instead of p."""
    if isinstance(fr, int):
        fr = [i / fr for i in range(fr + 1)]
    elif cluster == 'end':
        fr = [1.0 - f for f in reversed(fr)]
    return [(p[0] + (q[0] - p[0]) * f, p[1] + (q[1] - p[1]) * f) for f in fr]


def arc(th0, th1, n):
    """n+1 points on the injector circle, angle th0 -> th1 (radians)."""
    return [(CX + R * cos(th), R * sin(th))
            for th in (th0 + (th1 - th0) * i / n for i in range(n + 1))]


def coons(bottom, top, left, right):
    """Transfinite (Coons) interpolation from 4 matching edge point-lists.
    bottom/top have ni+1 pts (v=0/1); left/right have nj+1 pts (u=0/1)."""
    ni = len(bottom) - 1
    nj = len(left) - 1
    c00, c10 = bottom[0], bottom[-1]
    c01, c11 = top[0], top[-1]
    fx = [[0.0] * (ni + 1) for _ in range(nj + 1)]
    fy = [[0.0] * (ni + 1) for _ in range(nj + 1)]
    for j in range(nj + 1):
        v = j / nj
        for i in range(ni + 1):
            u = i / ni
            for c, B, T, L, Rt, C00, C10, C01, C11 in (
                (0, bottom[i][0], top[i][0], left[j][0], right[j][0],
                 c00[0], c10[0], c01[0], c11[0]),
                (1, bottom[i][1], top[i][1], left[j][1], right[j][1],
                 c00[1], c10[1], c01[1], c11[1])):
                val = ((1 - v) * B + v * T + (1 - u) * L + u * Rt
                       - ((1 - u) * (1 - v) * C00 + u * (1 - v) * C10
                          + (1 - u) * v * C01 + u * v * C11))
                (fx if c == 0 else fy)[j][i] = val
    return fx, fy


def n_arc_cells():
    """Cells per 45-deg quadrant of the injector circle -- the O-grid's azimuthal
    resolution.  It also fixes the core square at 2N x N cells, which is why the
    core's cell size is hl/N and scales with core_size alone."""
    return mg(n_arc if n_arc else max(2, round((R * pi / 4) / (dx / refine_inj))))


def injection_faces(r_in, r_out):
    """Return the seven O-grid cross-section faces (each a (fx,fy) pair), with
    the radial direction refined by r_in / r_out at the injector rim (1 = off).
    The four half-disk faces are also the ones used for the pipe."""
    d_inj = dx / refine_inj
    # cell counts (shared edges are forced equal by construction)
    N    = n_arc_cells()                         # 45-deg arc -> ny_core
    nxc  = 2 * N                                 # 90-deg top arc -> nx_core (uniform arc)
    # radial: rim cell = d_inj/rim_refine, growing back to d_inj -> the count
    # follows the refinement (unless n_pet / n_rad pin it)
    d_in  = d_inj / max(1.0, r_in)               # radial cell AT the rim, inside
    d_out = d_inj / max(1.0, r_out)              # radial cell AT the rim, outside
    npet = mg(n_pet if n_pet else ramp_count(R - hl, d_in, rim_growth, d_inj))
    nrad = mg(n_rad if n_rad else ramp_count(ob2 - R, d_out, rim_growth, d_inj))
    # one profile per family, reused on every edge of it (keeps shared edges equal)
    g_in  = ramp_fracs(npet, R - hl, d_in if r_in > 1.0 else 0.0, rim_growth, d_inj)
    g_out = ramp_fracs(nrad, ob2 - R, d_out if r_out > 1.0 else 0.0, rim_growth, d_inj)

    # key points ------------------------------------------------------------ #
    A_l0, A_r0 = (CX - R, 0.0), (CX + R, 0.0)          # arc ends on the wall
    A_l45, A_r45 = (CX - OB, OB), (CX + OB, OB)        # 45-deg arc points
    c3, c5 = (CX - hl, 0.0), (CX + hl, 0.0)            # core bottom corners
    c7, c8 = (CX - hl, hl), (CX + hl, hl)              # core top corners
    o_l0, o_r0 = (CX - ob2, 0.0), (CX + ob2, 0.0)      # outer box wall corners
    o_l, o_r = (CX - ob2, ob2), (CX + ob2, ob2)        # outer box top corners

    faces = {}

    # S3 core square  (i:+x  j:+y)
    faces['S3'] = coons(line(c3, c5, nxc), line(c7, c8, nxc),
                        line(c3, c7, N),  line(c5, c8, N))

    # S2 hole-left  (3001,3,7,3004) : petal radial, circle at p -> fine at p
    faces['S2'] = coons(line(A_l0, c3, g_in, 'start'), line(A_l45, c7, g_in, 'start'),
                        arc(pi, 0.75 * pi, N), line(c3, c7, N))
    # S5 hole-right (5,3002,3003,8) : petal radial, circle at q -> fine at q
    faces['S5'] = coons(line(c5, A_r0, g_in, 'end'), line(c8, A_r45, g_in, 'end'),
                        line(c5, c8, N), arc(0.0, 0.25 * pi, N))
    # S4 hole-top   (7,8,3003,3004) : left/right radial, circle at q -> fine at q
    faces['S4'] = coons(line(c7, c8, nxc), arc(0.75 * pi, 0.25 * pi, nxc),
                        line(c7, A_l45, g_in, 'end'), line(c8, A_r45, g_in, 'end'))

    # S1 ring-left  (12,3001,3004,6) : radial bottom/top, circle at q -> fine at q
    faces['S1'] = coons(line(o_l0, A_l0, g_out, 'end'), line(o_l, A_l45, g_out, 'end'),
                        line(o_l0, o_l, N), arc(pi, 0.75 * pi, N))
    # S6 ring-right (3002,11,4,3003) : radial, circle at p -> fine at p
    faces['S6'] = coons(line(A_r0, o_r0, g_out, 'start'), line(A_r45, o_r, g_out, 'start'),
                        arc(0.0, 0.25 * pi, N), line(o_r0, o_r, N))
    # S7 ring-top   (3004,3003,4,6) : left/right radial, circle at p -> fine at p
    faces['S7'] = coons(arc(0.75 * pi, 0.25 * pi, nxc), line(o_l, o_r, nxc),
                        line(A_l45, o_l, g_out, 'start'), line(A_r45, o_r, g_out, 'start'))

    return faces


def folded(face):
    """True if the Coons blend turned the face inside out anywhere (a cell whose
    orientation differs from the corner one)."""
    fx, fy = face
    ref = 0.0
    for j in range(len(fx) - 1):
        for i in range(len(fx[0]) - 1):
            ax, ay = fx[j][i + 1] - fx[j][i], fy[j][i + 1] - fy[j][i]
            bx, by = fx[j + 1][i] - fx[j][i], fy[j + 1][i] - fy[j][i]
            det = ax * by - ay * bx
            ref = ref or det
            if det * ref <= 0.0:
                return True
    return False


def build_faces():
    """injection_faces() with the rim refinement auto-limited per side.  The two
    sides are independent (r_in shapes S2/S4/S5, r_out shapes S1/S6/S7), and the
    outer one is capped by the top-ring trapezoid: past some clustering the Coons
    blend folds it.  Rather than let that reach the mesh, an over-ambitious
    factor is bisected down to the strongest value that still gives valid cells,
    and what was actually used is printed."""
    ask_in  = rim_refine_in  if rim_refine_in  is not None else rim_refine
    ask_out = rim_refine_out if rim_refine_out is not None else rim_refine

    def largest(r, ok):
        """Strongest factor <= r (>= 1) that ok() accepts."""
        if r <= 1.0 or ok(r):
            return r
        lo, hi = 1.0, r                        # lo taken as valid, hi known folded
        while hi - lo > 0.05:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if ok(mid) else (lo, mid)
        return lo

    r_in  = largest(ask_in,  lambda r: not any(folded(injection_faces(r, 1.0)[t])
                                               for t in ('S2', 'S4', 'S5')))
    r_out = largest(ask_out, lambda r: not any(folded(injection_faces(1.0, r)[t])
                                               for t in ('S1', 'S6', 'S7')))
    d_inj = dx / refine_inj
    N, arc_cell = n_arc_cells(), (R * pi / 4) / n_arc_cells()
    print(f"O-grid core (S3): half-width {hl:.4g} = {hl / R:g} R, {2 * N} x {N} cells "
          f"of {hl / N:.4g} = {hl / (R * pi / 4):.2f}x the arc cell {arc_cell:.4g}")
    print(f"  petals left with {R - hl:.3f} (on the axis) .. {(OB - hl) * sqrt(2):.3f} "
          f"(at 45 deg) of thickness")
    if core_size > core_max:
        print(f"  core_size {core_size:g} capped at core_max {core_max:g} "
              f"(0.707 would pinch the 45-deg petal to nothing)")
    print(f"rim refinement: in={r_in:g} out={r_out:g}  (radial cell at the rim "
          f"{d_inj / max(1.0, r_in):.4g} / {d_inj / max(1.0, r_out):.4g})")
    if (r_in, r_out) != (ask_in, ask_out):
        print(f"  relaxed from in={ask_in:g} out={ask_out:g} -- stronger folds a block")
    return injection_faces(r_in, r_out)


# ------------------------------- writers ----------------------------------- #
def write_p3d(path, blocks):
    with open(path, 'w') as f:
        f.write(f"{len(blocks)}\n")
        for _, dims, *_ in blocks:
            f.write(f"{dims[0]} {dims[1]} {dims[2]}\n")
        for _, _, X, Y, Z in blocks:
            for arr in (X, Y, Z):
                f.write("\n".join(f"{v*fs:.15g}" for v in arr) + "\n")


def write_tec(path, blocks):
    with open(path, 'w') as f:
        f.write(' TITLE     = "HyShot chimera mesh"\n')
        f.write(' VARIABLES = "X", "Y", "Z"\n')
        for k, (name, dims, X, Y, Z) in enumerate(blocks, 1):
            f.write(f' ZONE T="BLOCCO {k}"\n')
            f.write(f' I= {dims[0]}, J= {dims[1]}, K= {dims[2]}, ZONETYPE=Ordered\n')
            f.write(' DATAPACKING=BLOCK\n')
            for arr in (X, Y, Z):
                f.write("\n".join(f"{v*fs:.15g}" for v in arr) + "\n")


# ------------------------------- validation -------------------------------- #
def inverted_cells(dims, X, Y, Z):
    """Count hexahedra with non-positive Jacobian (left-handed / folded)."""
    ni, nj, nk = dims
    def P(i, j, k):
        q = (k * nj + j) * ni + i
        return (X[q], Y[q], Z[q])
    bad = 0
    for k in range(nk - 1):
        for j in range(nj - 1):
            for i in range(ni - 1):
                o = P(i, j, k)
                a = [p - q for p, q in zip(P(i + 1, j, k), o)]
                b = [p - q for p, q in zip(P(i, j + 1, k), o)]
                c = [p - q for p, q in zip(P(i, j, k + 1), o)]
                det = (a[0] * (b[1] * c[2] - b[2] * c[1])
                       - a[1] * (b[0] * c[2] - b[2] * c[0])
                       + a[2] * (b[0] * c[1] - b[1] * c[0]))
                if det <= 0.0:
                    bad += 1
    return bad


# --------------------------------- main ------------------------------------ #
def main():
    outs = sys.argv[1:] or ["mesh.p3d"]

    xL, xR = CX - ob2, CX + ob2                            # injection x-extent

    # ---- one node lattice shared by the 3 chamber blocks ----------------- #
    # split nodes fall exactly on xR (x), ob2 (y) and z_inj (z) so the three
    # chamber blocks connect conformally to each other.
    x_out = min(X_OUT, X_MAX)                              # never past the hardware
    if x_out < X_OUT:
        print(f"note: X_OUT {X_OUT:g} clamped to the modelled length {X_MAX:g}")
    xseg = uniform(xL, xR, dx)
    xdn  = uniform(xR, min(x_out, X_RAMP), dx)             # constant-section part
    if x_out > X_RAMP:                                     # + expansion, node at X_RAMP
        # grow out of the chamber spacing (dx at the kink -> capped at dx_exit)
        xdn = xdn[:-1] + (graded(X_RAMP, x_out, dx, exit_stretch, dx_exit, fine='low')
                          if exit_stretch > 1.0 else uniform(X_RAMP, x_out, dx_exit))
    xs   = xseg[:-1] + xdn                                 # node at xR = xs[iXR]
    iXR  = len(xseg) - 1
    yseg = uniform(0.0, ob2, dy)
    ys   = yseg[:-1] + uniform(ob2, W, dy)                 # node at ob2 = ys[iOB]
    iOB  = len(yseg) - 1
    zlo  = z_column(0.0, z_inj, wall=True, ceil=False)
    zs   = zlo[:-1] + z_column(z_inj, H, wall=False, ceil=CEIL_BL)   # node at z_inj
    iZI  = len(zlo) - 1

    z_ovr  = graded(0.0, z_inj, dz1, rbl, dz_core_inj, fine='low')    # injection: BL blends into refined core
    z_pipe = graded(-LPIPE, 0.0, dz1, rbl, dz_core_inj, fine='high')  # pipe: fine at the exit (z=0)

    blocks = []

    # 1) upstream block: x stretched (fine up_dxmin at the injection side x=xL),
    #    own spanwise (dy_up) and wall-normal core (dz_up); wall BL retained.
    x_up = stretch(X_IN, xL + front_overlap, up_dxmin, up_stretch, fine='high')
    z_up = z_column(0.0, H, wall=True, ceil=CEIL_BL, dcore=dz_up)
    blocks.append(("upstream       ",
                   *block_from_nodes(x_up, uniform(0.0, W, dy_up), z_up)))

    # 2) chamber, split into three blocks around the injection notch.  Only the
    #    downstream one reaches the exit expansion, so it is the only block with
    #    a non-flat ceiling (z_top is the identity while x <= X_RAMP).
    blocks.append(("cham downstream",
                   *block_ramped(xs[iXR:], ys, zs, z_top)))              # x[xR,X_OUT]
    blocks.append(("cham above     ",
                   *block_from_nodes(xs[:iXR + 1], ys[:iOB + 1], zs[iZI:])))
    blocks.append(("cham lateral   ",
                   *block_from_nodes(xs[:iXR + 1], ys[iOB:], zs)))

    # 3) injection overset : 7 O-grid blocks.  The four half-disk blocks (S2-S5,
    #    the injector bore) run through the pipe AND above the wall as ONE block
    #    each -- z_hole spans [-LPIPE, z_inj], fine at the exit plane z=0 from
    #    both sides.  The ring blocks (S1,S6,S7) sit on solid wall below z=0, so
    #    they exist above the wall only.
    faces = build_faces()
    HOLE = ('S2', 'S3', 'S4', 'S5')
    z_hole = z_pipe[:-1] + z_ovr           # merge pipe + above-wall (shared node at z=0)
    for tag in ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7'):
        fx, fy = faces[tag]
        zcol = z_hole if tag in HOLE else z_ovr
        label = f"inj {tag} (hole)" if tag in HOLE else f"inj {tag}       "
        blocks.append((label, *extrude(fx, fy, zcol)))

    npts = 0
    total_bad = 0
    fmg = 2 ** (MG_LEVEL - 1)
    mg_bad = 0
    if MG_LEVEL > 1:
        print(f"multigrid: level {MG_LEVEL} -> cell counts divisible by {fmg}")
    if x_out > X_RAMP and TAN_RAMP:
        hbl = bl_heights()[-1]
        core = H - hbl - (hbl if CEIL_BL else 0.0)         # stretchable core height
        f_str = (core + z_top(x_out) - H) / core           # core stretch at the outlet
        dzc = max(b - a for a, b in zip(zs, zs[1:]))       # actual core spacing
        print(f"exit expansion: ceiling {ramp_deg:g} deg from x={X_RAMP:g}, "
              f"height {H:g} -> {z_top(x_out):.3f} at x={x_out:g}"
              f"  (core dz {dzc:.4g} -> {dzc * f_str:.4g}, wall/ceiling BL kept)")
        nz = [b - a for a, b in zip(xs, xs[1:]) if a > X_RAMP - 1e-9]
        print(f"  streamwise: {len(nz)} cells, dx {nz[0]:.4g} -> {nz[-1]:.4g} "
              f"(chamber {dx:g}, growth {exit_stretch:g}/cell, cap {dx_exit:g})")
    print(f"{len(blocks)} blocks:")
    for name, dims, X, Y, Z in blocks:
        ni, nj, nk = dims
        bad = inverted_cells(dims, X, Y, Z)
        total_bad += bad
        cells = (ni - 1, nj - 1, nk - 1)          # coarsen-able dimension is #cells
        div = fmg > 1 and any(c % fmg for c in cells)
        mg_bad += div
        flag = f"  <-- {bad} INVERTED CELLS" if bad else ""
        flag += "  <-- NOT MG-DIVISIBLE" if div else ""
        print(f"  {name}  {ni:4d} x{nj:4d} x{nk:3d} = {ni*nj*nk:>9d}"
              f"  (cells {cells[0]}x{cells[1]}x{cells[2]}){flag}")
        npts += ni * nj * nk
    print(f"total grid points: {npts}")
    if total_bad:
        print(f"*** WARNING: {total_bad} inverted cells -- reduce rim_refine "
              f"or coarsen the injection ***")
    if mg_bad:
        print(f"*** WARNING: {mg_bad} block(s) not divisible by {fmg} ***")

    for out in outs:
        (write_tec if out.endswith(".tec") else write_p3d)(out, blocks)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
