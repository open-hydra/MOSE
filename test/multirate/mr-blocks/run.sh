#!/usr/bin/env bash
#===============================================================================
#  Multirate on several blocks, safety 2.
#  - equal levels across the connections: 1D acoustic wave in three blocks (one connection inside a band, one on the
#    edge of two bands, where the rule raises a cell), deep hierarchy in three blocks (L = 12) and 2D box in 2 x 2 blocks (L = 3): first macro step with
#    trace = .true.; every connection record has the level of its partner cell, level jumps <= 1 across every face;
#    against the same mesh in one block the map only rises (raised cells listed).
#  - one block against three blocks: acoustic wave in one block and in three blocks with both connections inside bands
#    (same level map), max-level 3, 100 macro steps: largest relative difference of rho, u, p between the two runs
#    <= 1e-11 and <= 10 x the same difference of the single-rate step.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu mr1d-acu-3b mr1d-acu-3bi mr1d-deep12 mr1d-deep12-3b mr2d-box mr2d-box-4b \
  || { echo "mr-blocks: FAIL - cannot build the cases"; exit 1; }
status=0
echo "  Equal levels across the connections:"
for spec in "mr1d-acu-3b mr1d-acu 3" "mr1d-deep12-3b mr1d-deep12 12" "mr2d-box-4b mr2d-box 3"; do
  set -- $spec; c=$1; c1=$2; L=$3
  T=("MOSE-Parameters|iter-threshold|1" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" \
     "MOSE-Multirate|safety|2.0" "MOSE-Multirate|trace|.true.")
  mr_run "runs/$c-map" "cases/$c" "${T[@]}" || { status=1; continue; }
  mr_run "runs/$c1-map" "cases/$c1" "${T[@]}" || { status=1; continue; }
  echo "  $c:"
  python3 -B verify_levels.py "runs/$c-map" "cases/$c" --single "runs/$c1-map" || status=1
done
echo "  One block against three blocks, same level map:"
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0")
st=0
for c in mr1d-acu mr1d-acu-3bi; do
  mr_run "runs/$c-mr" "cases/$c" "${MR[@]}" || st=1
  mr_run "runs/$c-sr" "cases/$c" || st=1
done
if [[ $st -eq 0 ]]; then
  python3 -B verify_one_block.py runs/mr1d-acu-mr runs/mr1d-acu-3bi-mr runs/mr1d-acu-sr runs/mr1d-acu-3bi-sr || status=1
else status=1; fi
if [[ $status -eq 0 ]]; then echo "mr-blocks: ok"; else echo "mr-blocks: FAIL"; fi
exit $status
