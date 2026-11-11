#!/usr/bin/env bash
#===============================================================================
#  Mirror symmetry (metamorphic test): 2D box with the banded mesh of mr2d-box (refined
#  cross, levels 0..3 with safety 2) and the pressure pulse at the centre of symmetry of the mesh, in one block and in
#  2 x 2 blocks; multirate (safety 2, max-level 3) and single rate, cfl 0.2, 100 macro steps: the fields are symmetric
#  in x, in y and about the diagonal (u <-> v) to 1e-10 of the excursion of each variable, and the multirate within 10
#  times the asymmetry of the single rate on the same case (round-off floor: the single rate alone reaches 3-4e-12).
#  The steps of the first macro step stay below the slab thickness / a (wall mode across the slab), checked from
#  the logs.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr2d-box-sym mr2d-box-sym-4b || { echo "mr-sym: FAIL - cannot build the cases"; exit 1; }
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0")
status=0
for c in mr2d-box-sym mr2d-box-sym-4b; do
  mr_run "runs/$c-mr" "cases/$c" "${MR[@]}" || status=1
  mr_run "runs/$c-sr" "cases/$c" || status=1
done
if [[ $status -eq 0 ]]; then
  python3 -B verify_sym.py runs/mr2d-box-sym-mr runs/mr2d-box-sym-sr runs/mr2d-box-sym-4b-mr runs/mr2d-box-sym-4b-sr || status=1
  python3 -B ../mr-common/verify_slab.py cases/mr2d-box-sym runs/mr2d-box-sym-mr runs/mr2d-box-sym-4b-mr --h-max || status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-sym: ok"; else echo "mr-sym: FAIL"; fi
exit $status
