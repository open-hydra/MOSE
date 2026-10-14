# Viscous Manufactured Solutions on Skewed Grids

The [viscous manufactured solution](2D-mms.md) measures the order of accuracy of
the viscous Navier–Stokes discretisation on a Cartesian mesh. This case runs the
same driver on two-block grids whose cells are skewed or change size abruptly,
where the face gradient of the viscous fluxes is exact for a linear field only if
its metric takes the same differences of the cell centres as of the field: the
`stencil` value of `diffusive-metric`, the default
([Numerics](../theory/numerics.md)).

!!! note
    The convective fluxes take the MUSCL reconstruction along the grid lines,
    which is not exact for linear fields where every cell is skewed the other way
    of its neighbours: on such grids the observed order of the whole scheme drops
    to about one, whatever the viscous metric.

## Grids and cases

Two grid families on the unit square, $N$ cells per side, one cell across in
$z$ with symmetry (`300`) on the $k$ faces; `make_case.py` in the test folder
builds the mesh, the initial field and the boundary records of every case:

- **zigzag**: the nodes of every other $i$ line moved by $\pm0.3\,h$ in $y$,
  weighted by $4s(1-s)$ so that the boundary stays straight: every cell is
  skewed, and its neighbours the opposite way;
- **iface**: a cell-size jump of 1.5 at $y = 1/2$ ($3N/5$ cells below, $2N/5$
  above) and grid lines kinked there ($x$ moved by $0.3\,y$ below and
  $0.3\,(1-y)$ above).

Each grid is made of two blocks that meet with permuted indices: on zigzag the
second half is turned by 180 degrees about $k$, on iface the upper part is turned
so that its $j$ face meets an $i$ face.

| Case | Field (`MOSE_MMS_FIELD`) | Boundaries |
|---|---|---|
| zigzag, periodic | `periodic2d` | periodic on every side |
| iface, periodic | `periodic2d` | periodic on every side |
| iface, walls | `wall2d` | periodic in $x$; no-slip walls at $y = 0$ (`302`, $T = T_w(x)$) and $y = 1$ (`301`, the exact heat flux) |

The fields and their sources are those of the
[manufactured-solution driver](2D-mms.md#methodology); the initial state is the
exact field at the cell centres.

## Setup

| Parameter | Value |
|---|---|
| Gas | calorically perfect air (`test/common/Air`), $\mu = 10$ (flat transport table), $\mathrm{Pr} = 0.72$, one species |
| Equations | Navier–Stokes (laminar) |
| Time scheme | RK3, local time stepping until the residual has fallen by $5\cdot10^{-4}$ from the first iteration |
| CFL / VNN | 0.2 / 0.1 |
| Space reconstruction | MUSCL, van Leer |
| Riemann solver | HLLC |
| Face metric | `diffusive-metric = stencil` (default) |
| Grids | short form (CTest): $N = 10, 20$ cells per side; full form: $N = 20, 40, 80$ |

## Checks

`verify.py` passes if, for every case:

1. the normalised $L_2$ error of $\rho$, $u$, $v$ and $p$ decreases on the finest
   pair of grids;
2. the observed $L_2$ order on the finest pair is at least 70 % of the order the
   family is expected to show: 0.7 on zigzag (order about one, see the note
   above) and 1.4 on iface (order two);
3. with walls, the heat flux integrated over the wall $y = 0$ converges to the
   exact one with an order of at least that of check 2;
4. on every grid, the totals $\sum \mathbf{U} V$ of the conserved variables whose
   boundary flux is zero (all of them with `periodic2d`, the mass with
   `wall2d`) differ from those of the initial state by less than the smallest
   $L_2$ error of the grid. The discrete steady problem of these equations fixes
   their solution only up to a mean, which the local time steps move during the
   iteration; this check makes sure that the mean so chosen does not dominate the
   measured error.

## Results

- **Short form** (`MMSViscousSkewed`, validation tier of CTest): passes its
  checks on every case: the four on the case with walls, checks 1, 2 and 4 on
  the two periodic ones.
- **Full form** (`./MOSE.sh test`): the refinement study of the three cases;
  `verify.py` then gives the errors and the observed orders of every case and
  the figures `OUTPUT/mms-skewed-grids.svg` (the coarsest grid of each family)
  and `OUTPUT/mms-skewed-<case>.svg` (the $L_2$ and $L_\infty$ errors against
  $N$).

## Running the case

```bash
cd test/2D/viscous/mms/skewed
./MOSE.sh testfast        # or: ./MOSE.sh test (full form)
python verify.py
```

The case runs `bin/MOSE_mms`, the variant of the solver that adds the analytic
source at every Runge–Kutta stage.

## References

See the [viscous manufactured solution](2D-mms.md#references).
