# 3D Viscous Manufactured Solutions on Skewed Grids

The 3D counterpart of the [skewed-grid manufactured solutions](2D-mms-skewed.md):
the same driver and checks with the fields `periodic3d` and `wall3d` of the
[manufactured-solution driver](2D-mms.md#methodology), on grids whose cells are
skewed in the three directions or change size abruptly, made of two blocks that
meet with permuted and reversed indices.

## Grids and cases

The generator of the 2D case (`test/2D/viscous/mms/skewed/make_case.py`) builds,
on the unit cube with $N$ cells per side:

- **zigzag**: the nodes of every other $i$ line moved by $\pm0.3\,h$ in $y$, also
  in $z$ alternating with $j$ and in $x$ alternating with $k$, weighted so that
  the boundary stays plane: every cell is skewed in the three directions, and its
  neighbours the opposite way;
- **iface**: a cell-size jump of 1.5 at $y = 1/2$ and grid lines kinked there
  ($x$ moved by $0.3\,y$ below and $0.3\,(1-y)$ above, $z$ by half of it).

The two blocks of each grid meet with permuted and reversed indices (on zigzag
$i$ and $k$ exchanged; on iface a $j$ face meets a $k$ face).

| Case | Field (`MOSE_MMS_FIELD`) | Boundaries | Grids |
|---|---|---|---|
| zigzag, periodic | `periodic3d` | periodic on every side | short form $N = 8, 12$; full form $N = 16, 24, 36$ |
| iface, periodic | `periodic3d` | periodic on every side | full form $N = 20, 30, 45$ |
| iface, walls | `wall3d` | periodic in $x$ and $z$; no-slip walls at $y = 0$ (`302`, $T = T_w(x, z)$) and $y = 1$ (`301`, the exact heat flux) | full form $N = 20, 30, 45$ |

The short form, run by CTest, leaves out the iface grids: on $N = 10, 15$ the case
with walls passes the checks below, but its convergence costs twenty to
twenty-five times the zigzag pair in CPU time.

## Setup

As in the [2D case](2D-mms-skewed.md#setup): calorically perfect air with
$\mu = 10$ and $\mathrm{Pr} = 0.72$ (the same gas files, written by ATLAS GPB), laminar Navier–Stokes, RK3 with local time
stepping until the residual has fallen by $5\cdot10^{-4}$, CFL 0.2 and VNN 0.1,
MUSCL with van Leer, HLLC, `diffusive-metric = stencil` (default).

## Checks

`verify.py` applies to every case run the four checks of the
[2D case](2D-mms-skewed.md#checks), with $w$ among the variables of check 1.

## Results

- **Short form** (`MMSViscous3D`, validation tier of CTest, zigzag on
  $N = 8, 12$): passes checks 1, 2 and 4 (check 3 needs walls, which the short
  form does not run).
- **Full form** (`./MOSE.sh test`): the refinement study of the cases;
  `verify.py` then gives the errors and the observed orders of every case and
  the figures `OUTPUT/mms3d-grids.svg` (the face z = 0 of the coarsest grid of
  each family) and `OUTPUT/mms3d-<case>.svg` (the $L_2$ and $L_\infty$ errors
  against $N$).

## Running the case

```bash
cd test/3D/viscous/mms
./MOSE.sh testfast        # or: ./MOSE.sh test (full form)
python verify.py
```

## References

See the [viscous manufactured solution](2D-mms.md#references).
