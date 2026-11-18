# 2D Flat Plate – Laminar Boundary Layer with the FAS Multigrid Cycle

The results are discussed in the documentation, [Multigrid Cycles (FAS)](../../../../docs/vv/fas-cycles.md); this file holds the setup of the case and how to run it.

The case is the laminar flat plate of [`../flat-plate-laminar`](../flat-plate-laminar/README.md), with the same mesh, gas, initial field and scheme, solved with the FAS multigrid cycle at its default settings. Its `input.ini` is that of the plate with one more section:

```ini
[MOSE-Multigrid]
levels      = 3
level3-iter = 2000
level2-iter = 2000
level1-iter = 100000
cycle       = w
```

`INPUT/bc2.txt` and `INPUT/bc3.txt` are the boundary records of levels 2 and 3, written by ATLAS BCB from the `[BCB-Block1]` section of `input.ini`.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK2 |
| CFL | 0.5 |
| VNN | 0.5 |
| Space reconstruction | MUSCL |
| Flux limiter | Van Leer |
| Riemann solver | HLLC |
| Grid levels | 3, 2000 iterations on each coarser level |
| Multigrid cycle | `cycle = w`, the other keys at their defaults |

## Running the case

```bash
cd test/2D/viscous/flat-plate-laminar-fas
./MOSE.sh solve       # the run of input.ini, output in OUTPUT/
./MOSE.sh test        # CTest FlatPlateLaminar.FAS
python3 compare.py    # after ./MOSE.sh test: skin friction of the cycle against the grid sequencing
```

`test` runs the plate to convergence twice, in `RUN/fas` with the cycle and in `RUN/seq` with grid sequencing (`cycle = none`): at most 100 000 iterations of the finest level, with a residual threshold of $10^{-10}$ of the first residual of each level. The `verify.py` of `../flat-plate-laminar` must find the skin friction of both runs within 2.5 % RMS of the Blasius solution for $x \ge 0.05$ m, and `compare.py` the two skin frictions within 0.5 % RMS of each other: the two methods reach the same solution. The cycle must also reach a density residual of $10^{-7}$ (`OUTPUT/residual-history.dat`) in at most half the iterations of the finest level of grid sequencing.

The runs of `test` write the solution at the end only. `RUN/` and `bin/` are not tracked.
