"""Generate the Shu-Osher initial condition (INPUT/ic.tec) and boundary file
(INPUT/bc.txt) for a given number of cells.

Shu-Osher problem (Mach-3 shock into a sine-perturbed density field):
    x < -4 :  rho = 3.857143,  u = 2.629369,  p = 10.33333   (post-shock)
    x >= -4:  rho = 1 + 0.2 sin(5x),  u = 0,  p = 1
    domain x in [-5, 5],  evolve to t = 1.8

The mesh is a thin 1-D slab (I = N+1 nodes, J = K = 2) in the same Tecplot
BLOCK layout MOSE reads for the 1-D test cases (nodal coordinates, cell-centred
state, variable order x y z rho u v w p). Values are written one per line in
list-directed-readable form.

Usage:  python gen_ic.py [N] [tag]
    N     number of cells               (default 200)
    tag   if given, write ic_x<tag>.tec / bc_x<tag>.txt instead of the
          default ic.tec / bc.txt       (used to stage the convergence grids)
"""
import sys
from math import sin

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
TAG = sys.argv[2] if len(sys.argv) > 2 else ""
IC_FILE = f"INPUT/ic_x{TAG}.tec" if TAG else "INPUT/ic.tec"
BC_FILE = f"INPUT/bc_x{TAG}.txt" if TAG else "INPUT/bc.txt"
XMIN, XMAX = -5.0, 5.0

I, J, K = N + 1, 2, 2
dx = (XMAX - XMIN) / N
xn = [XMIN + i * dx for i in range(I)]
yn = [0.0, dx]
zn = [-0.5 * dx, 0.5 * dx]


def state(x):
    if x < -4.0:
        return 3.857143, 2.629369, 10.33333
    return 1.0 + 0.2 * sin(5.0 * x), 0.0, 1.0


def w(fh, v):
    fh.write(f" {v: .15E}\n")


# ---- ic.tec ---------------------------------------------------------------
with open(IC_FILE, "w") as fh:
    fh.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
    fh.write(f' ZONE  T = B1-IG, I={I}, J={J}, K={K}, DATAPACKING=BLOCK, '
             'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n')
    # nodal coordinates: i fastest, then j, then k
    for coord in (lambda i, j, k: xn[i],
                  lambda i, j, k: yn[j],
                  lambda i, j, k: zn[k]):
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    w(fh, coord(i, j, k))
    # cell-centred state: rho, u, v, w, p
    xc = [0.5 * (xn[i] + xn[i + 1]) for i in range(I - 1)]
    cells = [state(x) for x in xc]
    for comp in (lambda s: s[0],          # rho
                 lambda s: s[1],          # u
                 lambda s: 0.0,           # v
                 lambda s: 0.0,           # w
                 lambda s: s[2]):         # p
        for k in range(K - 1):
            for j in range(J - 1):
                for s in cells:
                    w(fh, comp(s))

# ---- bc.txt : transmissive (extrapolation) ends, mirrors Sod layout -------
# faces: 1=imin 2=imax 3=jmin 4=jmax 5=kmin 6=kmax
with open(BC_FILE, "w") as fh:
    def line(i, j, k, face, typ):
        fh.write(f"{1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{typ:8d}\n")
    line(1, 1, 1, 1, 400)                     # imin (cell i=1)
    line(N, 1, 1, 2, 400)                     # imax (cell i=N)
    for i in range(1, N + 1):                 # jmin
        line(i, 1, 1, 3, 400)
    for i in range(1, N + 1):                 # jmax
        line(i, 1, 1, 4, 400)
    for i in range(1, N + 1):                 # kmin (degenerate)
        line(i, 1, 1, 5, 0)
    for i in range(1, N + 1):                 # kmax (degenerate)
        line(i, 1, 1, 6, 0)

print(f"wrote {IC_FILE} and {BC_FILE} for N={N} cells (x in [{XMIN},{XMAX}])")
