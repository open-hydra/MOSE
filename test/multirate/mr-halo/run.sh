#!/usr/bin/env bash
#===============================================================================
#  Halo from the true stencils: safety 2, the halo cells and the ghost records of every
#  level written by MOSE (trace = .true.) are identical to an independent construction in Python with the same stencil
#  offsets. Cases: 1D acoustic wave (L = 3), deep hierarchy (L = 8), 2D box (L = 3); deep hierarchies with bands
#  of one cell (L = 8 and L = 12), where the halo of a level l also holds cells of levels m <= l - 2, required for every
#  l >= 2; acoustic wave in three blocks, deep hierarchy with 13 levels (L = 12) in three blocks
#  and 2D box in 2 x 2 blocks, whose halos cross the block connections. Laminar Navier-Stokes, with the
#  diagonals of the diffusive faces and the tangential neighbours of the partner cells of the connections (Pg(:,3:6)):
#  flat plate in 2 x 2 blocks (L = 3) and skewed hot square in one block and in 2 x 2 blocks (L = 1), whose square of
#  level 1 has sharp corners where some halo cells are reached only through the diagonals (required); its 3D twin, the
#  hot cube (L = 1, diagonals in k too). Spalart-Allmaras and SST, with the stencil of the
#  turbulence sources (+-1 along the axes) and the ghost cells of the wall records: the plate with Spalart-Allmaras in
#  2 x 2 blocks and the wall channel with SST (L = 3). The 2-D axisymmetric wedge with levels (mr2d-wedge-band, laminar,
#  L = 2: faces between levels on cells with the records 200 of the wedge). In the viscous cases, at every stage of the first macro step, the
#  viscous variables that the viscous fluxes of the stage read (blk%Visc) are those of the primitive variables of their
#  cells at that point (lines MR-VISC): no halo cell, finer cell or ghost cell, diagonals included, read stale.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu mr1d-deep mr2d-box mr1d-deep1 mr1d-deep12-1 mr1d-acu-3b mr1d-deep12-3b mr2d-box-4b \
               mr2d-plate-4b mr2d-hot mr2d-hot-4b mr3d-hot mr2d-plate-sa-4b mr2d-chan-sst mr2d-wedge-band \
  || { echo "mr-halo: FAIL - cannot build the cases"; exit 1; }
status=0
for cl in "mr1d-acu 3" "mr1d-deep 8" "mr2d-box 3" "mr1d-deep1 8 --require-deep" "mr1d-deep12-1 12 --require-deep" \
          "mr1d-acu-3b 3" "mr1d-deep12-3b 12" "mr2d-box-4b 3" "mr2d-plate-4b 3" "mr2d-hot 1 --require-diag" \
          "mr2d-hot-4b 1 --require-diag" "mr3d-hot 1 --require-diag" "mr2d-plate-sa-4b 3" "mr2d-chan-sst 3" \
          "mr2d-wedge-band 2"; do
  set -- $cl; c=$1; L=$2; shift 2
  mr_run "runs/$c" "cases/$c" "MOSE-Parameters|iter-threshold|1" "MOSE-Multirate|enabled|.true." \
         "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|trace|.true." || { status=1; continue; }
  echo "  $c:"
  python3 -B verify_halo.py "runs/$c" "cases/$c" "$@" || status=1
done
if [[ $status -eq 0 ]]; then echo "mr-halo: ok"; else echo "mr-halo: FAIL"; fi
exit $status
