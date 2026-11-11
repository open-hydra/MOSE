#!/usr/bin/env bash
#===============================================================================
#  Multirate, 1 against 3 and 4 OpenMP threads: Sod tube with the levels
#  recomputed at every macro step, 2D box with a pressure pulse, 1D acoustic wave in three blocks and 2D box in 2 x 2
#  blocks (max-level 3); laminar Navier-Stokes: flat plate in 2 x 2 blocks (max-level 3) and skewed hot square in
#  2 x 2 blocks (max-level 1); its 3D twin, the hot cube (max-level 1); Spalart-Allmaras and SST:
#  the turbulent flat plate with Spalart-Allmaras in 2 x 2 blocks and the wall channel with SST (max-level 3);
#  safety 2, 100 macro steps: fields, wall output and residual history (written at every step, res-diter = 1) identical
#  byte by byte. Three threads:
#  with schedule(static) the split of the lists changes with the number of threads.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-sod mr2d-box mr1d-acu-3b mr2d-box-4b mr2d-plate-4b mr2d-hot-4b mr3d-hot mr2d-plate-sa-4b mr2d-chan-sst \
  || { echo "mr-threads: FAIL - cannot build the cases"; exit 1; }
status=0
for spec in "mr1d-sod 3" "mr2d-box 3" "mr1d-acu-3b 3" "mr2d-box-4b 3" "mr2d-plate-4b 3" "mr2d-hot-4b 1" "mr3d-hot 1" \
            "mr2d-plate-sa-4b 3" "mr2d-chan-sst 3"; do
  set -- $spec; c=$1; L=$2
  MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0" "MOSE-IO|res-diter|1")
  lv=(); [[ $c == mr1d-sod ]] && lv=("MOSE-Multirate|levels|recompute")
  for t in 1 3 4; do
    OMP_NUM_THREADS=$t mr_run "runs/$c-t$t" "cases/$c" "${MR[@]}" "${lv[@]}" || { status=1; continue 2; }
  done
  for t in 3 4; do
    if mr_same "runs/$c-t1" "runs/$c-t$t"; then echo "  $c: 1 and $t threads identical"; else echo "  $c: 1 and $t threads DIFFERENT"; status=1; fi
  done
done
if [[ $status -eq 0 ]]; then echo "mr-threads: ok"; else echo "mr-threads: FAIL"; fi
exit $status
