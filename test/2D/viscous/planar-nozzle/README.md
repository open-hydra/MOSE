# Planar Nozzle on Grid Families

The results are discussed in the documentation, [Planar Nozzle on Grid Families](../../../../docs/vv/2D-planar-nozzle-grids.md); this file holds the setup of the case and how to run it.

## Setup

| Parameter | Value |
|---|---|
| Gas | air: perfect gas with $\gamma = 1.4$ and $R = 287$ J/(kg·K) (`test/common/Air-viscous`, the files that ATLAS GPB writes at ATLAS commit `ebeaf4f` for an ideal gas `Air` with `gamma = 1.4`, `R = 287.0`, `mil = 1.1858685985e-5` and `Pr = 0.72`), viscosity of Sutherland's law $\mu = 1.458\times10^{-6}\,T^{3/2}/(T + 110.4)$ Pa·s (`test/common/Air-sutherland/transport.dat`), $\mathrm{Pr} = 0.72$, $\mathrm{Pr}_t = 0.9$ |
| Turbulence | Spalart–Allmaras; SST (`omega-wall-bc = practical`) |
| Time scheme | RK3, local time stepping (steady) |
| CFL / VNN | 0.8 / 0.3 |
| Space reconstruction | MUSCL, van Leer |
| Riemann solver | HLLC |

## Running the case

```bash
cd test/2D/viscous/planar-nozzle
./MOSE.sh testfast        # or: ./MOSE.sh -p <threads> test (full form, long)
python verify.py
python plot_grids.py      # the figure of the grid families (no run needed)
```
