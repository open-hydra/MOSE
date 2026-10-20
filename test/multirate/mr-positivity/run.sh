#!/usr/bin/env bash
#===============================================================================
#  Multirate with explicit turbulence sources, positivity: the positivity cases of the test
#  TimeIntegration/TurbulenceSources (near-wall turbulence far above its equilibrium: the
#  destruction acts faster than the transport) with the multirate, Spalart-Allmaras and SST, max-level 3, the levels
#  recomputed at every macro step (the stiffness of the current state sets them), safety 2, and max-level 0 with safety 1
#  (the stiffest cells, which set the step, run at h = dtcap: the worst case),
#  200 macro steps. Checks (verify_positivity.py): the run ends; no stage makes nu~, k or omega
#  negative (count printed by the solver); the local step was limited by the turbulence sources. Positive control:
#  fixed-macro-step 48 with the map of that step (max-level 3, safety 2) takes the finest cells to 12 times their limit:
#  the run must end and print the count, otherwise a missing count does not show that no stage went negative. With 32 to
#  64 both models end and print it (Spalart-Allmaras prints nothing with 24 and diverges with 96, SST prints nothing with 16
#  and diverges with 128, with gfortran and ifx alike).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr2d-chanpos-sa mr2d-chanpos-sst || { echo "mr-positivity: FAIL - cannot build the cases"; exit 1; }
status=0
for c in mr2d-chanpos-sa mr2d-chanpos-sst; do
  for s in 2 1; do
    L=3; [[ $s == 1 ]] && L=0
    mr_run "runs/$c-s$s" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" \
      "MOSE-Multirate|safety|$s.0" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" \
      "MOSE-Parameters|iter-threshold|200" "MOSE-IO|shell-diter|200" || { status=1; continue; }
    python3 -B verify_positivity.py "cases/$c" "runs/$c-s$s" || status=1
  done
done
for c in mr2d-chanpos-sa mr2d-chanpos-sst; do
  r=runs/$c-control
  if mr_run "$r" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" \
       "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" "MOSE-Multirate|fixed-macro-step|48.0" \
       "MOSE-Parameters|iter-threshold|200" "MOSE-IO|shell-diter|200"; then
    n=$(command grep -a -m1 'Turbulence variables below zero after a stage' "$r/solve.log" | sed -E 's/.*: *//')
    if [[ -n $n ]]; then echo "  $c positive control: negative nu~, k or omega after a stage: $n"
    else echo "  $c positive control: no count printed, the check cannot see a negative value"; status=1; fi
  else status=1; fi
done
if [[ $status -eq 0 ]]; then echo "mr-positivity: ok"; else echo "mr-positivity: FAIL"; fi
exit $status
