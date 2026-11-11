#!/usr/bin/env bash
#===============================================================================
#  Sod tube with levels: 1D tube closed by slip walls, banded mesh, levels
#  recomputed at every macro step (safety 2, max-level 3), up to t >= 1e-3 s: the run completes with positive density
#  and pressure; mass and energy drift <= 1e-12; on x >= 0.15 m (outside the reflection of the rarefaction head):
#  L1 error against the exact Riemann solution <= 1.05 x the single rate's, shock within one local cell of the exact
#  position, no more local extrema of the density than the single rate.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-sod || { echo "mr-sod: FAIL - cannot build the case"; exit 1; }
T=("MOSE-Parameters|time-threshold|1.0e-3" "MOSE-Parameters|iter-threshold|1000000")
status=0
mr_run runs/sr cases/mr1d-sod "${T[@]}" || status=1
mr_run runs/mr cases/mr1d-sod "${T[@]}" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" \
       "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" || status=1
if [[ $status -eq 0 ]]; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-sod runs/mr --quantities mass,energy --max 1e-12 || status=1
  python3 -B verify_sod.py runs/mr runs/sr --xmin 0.15 || status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-sod: ok"; else echo "mr-sod: FAIL"; fi
exit $status
