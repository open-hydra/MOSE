#!/usr/bin/env bash
#===============================================================================
#  Free-stream preservation with levels: uniform flow u = 100, v = 50 m/s on
#  a curvilinear, stretched and skewed 2D mesh periodic in x and y, in one block and in 2 x 2 blocks (connections 101,
#  periodic 201), safety 2, max-level 2, levels fixed from the mesh, cfl 0.1, res-threshold 1e-300, 100 macro steps:
#  every conserved variable stays uniform to machine zero (<= 1e-13). The slab is one cell (1 mm) thick between slip
#  walls: the steps of the first macro step stay <= dz/a, below the onset of the wall mode across the slab, checked
#  from the log.
#  Informative (no threshold): the same at cfl 0.8, beyond the TVD limit of MUSCL with Heun on cells close to
#  isotropic (sum of the directional Courant numbers up to 0.8 > 0.5).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr2d-free mr2d-free-4b || { echo "mr-free: FAIL - cannot build the cases"; exit 1; }
MR=("MOSE-Parameters|res-threshold|1e-300" "MOSE-Parameters|iter-threshold|100" "MOSE-Multirate|enabled|.true."
    "MOSE-Multirate|max-level|2" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|fixed")
status=0
for c in mr2d-free mr2d-free-4b; do
  mr_run "runs/$c-mr" "cases/$c" "${MR[@]}" "MOSE-Numerics|cfl|0.1" || { status=1; continue; }
  mr_run "runs/$c-sr" "cases/$c" "MOSE-Parameters|res-threshold|1e-300" "MOSE-Parameters|iter-threshold|300" "MOSE-Numerics|cfl|0.1" \
    || { status=1; continue; }
  python3 -B verify_free.py "cases/$c" "runs/$c-mr" --max 1e-13 || status=1
  python3 -B verify_free.py "cases/$c" "runs/$c-sr" --report
  python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-mr" --h-max || status=1
done
echo "  informative, cfl 0.8 (expected to lose the uniform flow):"
if mr_run runs/mr2d-free-mr-cfl08 cases/mr2d-free "${MR[@]}" "MOSE-Numerics|cfl|0.8"; then
  python3 -B verify_free.py cases/mr2d-free runs/mr2d-free-mr-cfl08 --report
else
  echo "  mr2d-free-mr-cfl08: did not complete (see runs/mr2d-free-mr-cfl08/solve.log)"
fi
if [[ $status -eq 0 ]]; then echo "mr-free: ok"; else echo "mr-free: FAIL"; fi
exit $status
