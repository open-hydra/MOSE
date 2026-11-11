#!/usr/bin/env bash
#===============================================================================
#  Planar shock through the faces between levels: mesh of mr2d-box (refined cross,
#  levels 0..3 with safety 2), Riemann problem at x = 0.2 m (1.2 kg/m3, 1e5 Pa | 0.15 kg/m3, 1e4 Pa), slip walls,
#  cfl 0.2, up to t >= 1e-3 s: the shock crosses the refined cross and the faces between levels in x and in y. Two
#  multirate runs, levels fixed from the initial state and recomputed at every macro step, and the single rate at dt_min.
#  Each multirate run: exit 0, rho and p > 0; mass and energy drift <= 1e-12; L1(rho - exact) on x >= 0.25 m <= 1.05 x
#  the single rate's; on every row no more local extrema of rho than the single rate. The steps of the first macro
#  step stay below the slab thickness / a (wall mode across the slab), checked from the logs.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr2d-shock || { echo "mr-shock: FAIL - cannot build the case"; exit 1; }
T=("MOSE-Parameters|time-threshold|1.0e-3" "MOSE-Parameters|iter-threshold|1000000")
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0")
status=0
mr_run runs/sr cases/mr2d-shock "${T[@]}" || status=1
mr_run runs/mr-fixed cases/mr2d-shock "${T[@]}" "${MR[@]}" "MOSE-Multirate|levels|fixed" || status=1
mr_run runs/mr-recompute cases/mr2d-shock "${T[@]}" "${MR[@]}" "MOSE-Multirate|levels|recompute" \
       "MOSE-Multirate|recompute-every|1" || status=1
if [[ $status -eq 0 ]]; then
  for r in runs/mr-fixed runs/mr-recompute; do
    python3 -B ../mr-common/verify_conservation.py cases/mr2d-shock "$r" --quantities mass,energy --max 1e-12 || status=1
  done
  python3 -B ../mr-common/verify_conservation.py cases/mr2d-shock runs/sr --quantities mass,energy --report
  python3 -B verify_shock.py runs/sr runs/mr-fixed runs/mr-recompute || status=1
  python3 -B ../mr-common/verify_slab.py cases/mr2d-shock runs/mr-fixed runs/mr-recompute --h-max || status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-shock: ok"; else echo "mr-shock: FAIL"; fi
exit $status
