# Verification & Validation

This section documents the Verification & Validation (V&V) test suite for MOSE.  The purpose is to demonstrate that the solver produces correct results across a range of physical regimes, by comparison against analytical or well-established reference solutions.

## Test suite

| Test | Dim | Mach | Physics | Verification | Solver(s) | Ref |
|---|---|---|---|---|---|---|
| [Sod Shock Tube](1D.md#sod-shock-tube) | 1D | subsonic–supersonic | Shock, contact, rarefaction | Analytical | SLAU | Sod (1978) |
| [Einfeldt Double Rarefaction](1D.md#einfeldt-double-rarefaction) | 1D | subsonic | Near-vacuum, two rarefactions | Analytical | HLLE | Einfeldt et al. (1991) |
| [Noh Implosion](1D.md#noh-implosion-problem) | 1D | supersonic | Strong shock, density jump | Analytical | HLLE | Noh (1987) |
| [Toro Test 3](1D.md#toro-test-case-3) | 1D | subsonic–supersonic | Compound wave (shock + contact + shock) | Analytical | HLLE++ | Toro (1999) |
| [Woodward-Colella Step](2D-woodward-colella.md) | 2D | 3.0 | Oblique shock, expansion fan, interactions | OpenFOAM reference | HLLE | Woodward & Colella (1984) |
| [Oblique Shock](2D-oblique-shock.md) | 2D | 4.0 | Oblique shock, post-shock instabilities | Analytical | HLLC, HLLC+, HLLE, SLAU | oblique shock theory |
| [Hypersonic Cylinder](2D-hypersonic-cylinder.md) | 2D | 8.1 | Bow shock, carbuncle | Analytical | HLLC, HLLC+, HLLE, SLAU | normal shock theory |
| [Rocket Nozzle](2D-nozzle.md) | 2D | transonic–supersonic | Multi-species frozen expansion | CEA 1D reference | HLLC | CEA/NASA |
| [Laminar Flat Plate](2D-flat-plate-laminar.md) | 2D | 0.2 | Laminar boundary layer, viscous effects | Blasius similarity | HLLC | Blasius (1908) |
| [Turbulent Flat Plate](2D-flat-plate-turbulent.md) | 2D | 0.2 | Turbulent boundary layer, turbulence models | NASA solver comparison (CFL3D, FUN3D) | HLLC | NASA |
| [Shock Wave-Boundary Layer Interaction](2D-swbli.md) | 2D | 5.0 | Shock-boundary-layer interaction, separation/reattachment | Schulein/SU2/Wind-US comparison | HLLC (SA) | Schulein + code-to-code |

## Running the tests

Each test case lives under `test/` and contains a `verify.py` script.  The script reads the solver output, compares it against the exact solution, and writes output figures.

```bash
cd test/1D/<TestName>
./MOSE.sh solve
python verify.py          # check errors and export figure
python verify.py --plot   # as above, also display figure interactively
```
