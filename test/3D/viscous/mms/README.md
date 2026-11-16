# 3D Viscous Manufactured Solutions on Skewed Grids

The results are discussed in the documentation, [3D Viscous Manufactured Solutions on Skewed Grids](../../../../docs/vv/3D-mms.md); this file holds the setup of the case and how to run it.

## Setup

As in the [2D case](../../../../docs/vv/2D-mms-skewed.md#setup): calorically perfect air with
$\mu = 10$ and $\mathrm{Pr} = 0.72$ (the same gas files, those of ATLAS GPB at ATLAS commit `ebeaf4f`), laminar Navier–Stokes, RK3 with local time
stepping until the residual has fallen by $5\cdot10^{-4}$, CFL 0.2 and VNN 0.1,
MUSCL with van Leer, HLLC, `diffusive-metric = stencil` (default).

## Running the case

```bash
cd test/3D/viscous/mms
./MOSE.sh testfast        # or: ./MOSE.sh test (full form)
python verify.py
```
