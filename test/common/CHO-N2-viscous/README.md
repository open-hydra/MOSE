# CHO-N2-viscous

The eight species of the surface reactions of a carbon wall (BC 505 and 506: O2, O, OH, H, H2, H2O, CO, CO2) and N2 as a
diluent: `phase.txt`, `thermo.dat` (1..3000 K, 1 K step), `transport.dat` and `composition.txt` (elements C H N O of each
species), for the columns of these walls in `test/2D/viscous/ablating-wall`, which link `phase.txt`, `thermo.dat` and
`transport.dat` (`composition.txt` completes the product of GPB; no test reads it).

Provenance: ATLAS GPB, `input.ini` of this folder, ATLAS `1b3c3d5` (its ideal-gas sources as at `34c128e`), Python 3.13.2
with Cantera 3.1.0 (conda environment ct-env); thermodynamics NASA-9 (`database/thermo/nasa9.yaml` of ATLAS): absolute
enthalpy on the CEA datum, the elements at 0 at 298.15 K, J/kg; transport of CEA. Molar masses of `phase.txt` = the sums of
Cantera's element weights (C 12.011, H 1.008, O 15.999, N 14.007). A run of GPB on `input.ini` in a copy of this folder
writes the same four files, byte for byte, into `fromATLAStoSolver/`.
