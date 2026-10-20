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
- **Full form** (`./MOSE.sh test`; the CTest test `MMSViscous3DFull`, with the
  `verification` label, outside the CI of the pull requests): the refinement
  study of the three cases. With `MOSE_mms` of commit `17367503009f` and the
  inputs of this folder, `verify.py` passes every check of every case, gives the
  errors and the observed orders below and draws the figures
  `OUTPUT/mms3d-grids.svg` (the face $z = 0$ of the coarsest grid of each family)
  and `OUTPUT/mms3d-<case>.svg` (the $L_2$ and $L_\infty$ errors against $N$).

The face $z = 0$ of the coarsest grid of each family, the two blocks in two colours:

<figure>
  {% include "vv/images/mms3d-grids.svg" %}
</figure>

**Normalised $L_2$ error (%) of the full form**

| Case | Variable | $N = 16$ | $N = 24$ | $N = 36$ | order $16\!\to\!24$ | order $24\!\to\!36$ |
|---|---|---|---|---|---|---|
| zigzag, periodic | $\rho$ | 0.570 | 0.370 | 0.248 | 1.07 | 0.99 |
| | $u$ | 0.434 | 0.231 | 0.131 | 1.55 | 1.39 |
| | $v$ | 0.456 | 0.241 | 0.130 | 1.57 | 1.52 |
| | $w$ | 0.892 | 0.482 | 0.273 | 1.52 | 1.40 |
| | $p$ | 0.584 | 0.379 | 0.252 | 1.07 | 1.01 |

| Case | Variable | $N = 20$ | $N = 30$ | $N = 45$ | order $20\!\to\!30$ | order $30\!\to\!45$ |
|---|---|---|---|---|---|---|
| iface, periodic | $\rho$ | 0.229 | 0.111 | 0.0548 | 1.78 | 1.74 |
| | $u$ | 0.105 | 0.0390 | 0.0145 | 2.45 | 2.44 |
| | $v$ | 0.182 | 0.0832 | 0.0375 | 1.94 | 1.96 |
| | $w$ | 0.218 | 0.0802 | 0.0296 | 2.47 | 2.46 |
| | $p$ | 0.209 | 0.102 | 0.0511 | 1.76 | 1.71 |
| iface, walls | $\rho$ | 0.192 | 0.0927 | 0.0462 | 1.80 | 1.72 |
| | $u$ | 0.270 | 0.119 | 0.0523 | 2.02 | 2.03 |
| | $v$ | 1.41 | 0.555 | 0.231 | 2.30 | 2.16 |
| | $w$ | 0.349 | 0.151 | 0.0654 | 2.06 | 2.07 |
| | $p$ | 0.191 | 0.0929 | 0.0467 | 1.78 | 1.70 |

With walls, the heat flux integrated over the wall $y = 0$ differs from the
exact one by 0.214 %, 0.0923 % and 0.0366 % (orders 2.07 and 2.28). On the
finest pair the $L_2$ orders are 0.99 to 1.52 on zigzag and 1.70 to 2.46 on
iface.

<figure>
  {% include "vv/images/mms3d-zigzag.svg" %}
</figure>

<figure>
  {% include "vv/images/mms3d-iface.svg" %}
</figure>

<figure>
  {% include "vv/images/mms3d-iface-walls.svg" %}
</figure>

## Running the case

```bash
cd test/3D/viscous/mms
./MOSE.sh testfast        # or: ./MOSE.sh test (full form)
python verify.py
```

## References

See the [viscous manufactured solution](2D-mms.md#references).
