# Viscous Manufactured Solutions on Skewed Grids

The results are discussed in the documentation, [Viscous Manufactured Solutions on Skewed Grids](../../../../docs/vv/2D-mms-skewed.md); this file holds the setup of the case and how to run it.

## Setup

| Parameter | Value |
|---|---|
| Gas | calorically perfect air, $\gamma = 1.4$, $R = 287$ J/(kg·K), $\mu = 10$, $\mathrm{Pr} = 0.72$, one species |
| Equations | Navier–Stokes (laminar) |
| Time scheme | RK3, local time stepping until the residual has fallen by $5\cdot10^{-4}$ from the first iteration |
| CFL / VNN | 0.2 / 0.1 |
| Space reconstruction | MUSCL, van Leer |
| Riemann solver | HLLC |
| Face metric | `diffusive-metric = stencil` (default) |
| Grids | short form (CTest): $N = 10, 20$ cells per side; full form: $N = 20, 40, 80$ |

The gas is that of `test/common/Air`, the files that ATLAS GPB writes from `test/common/Air/input.ini` at ATLAS commit `ebeaf4f`, with the flat transport
table `test/2D/viscous/mms/INPUT/transport.dat` (the values of ATLAS GPB with `mil = 10` and `Pr = 0.72`).

## Running the case

```bash
cd test/2D/viscous/mms/skewed
./MOSE.sh testfast        # or: ./MOSE.sh test (full form)
python verify.py
```

The case runs `MOSE_mms`, the variant of the solver that adds the analytic
source at every Runge–Kutta stage. It is built with the test programs in
`build/test/` (CTest passes its path; by hand, `MOSE.sh` takes `build/test/MOSE_mms`
or the file named by `MOSE_MMS_BIN`).
