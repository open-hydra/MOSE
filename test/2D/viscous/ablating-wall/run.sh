#!/usr/bin/env bash
#===============================================================================
#  Ablative wall boundary conditions -- BC 503 (melting) and 504 (pyrolysis)
#
#  Runs the blown-column case once per wall model and checks each result against
#  the analytic solution (see build_case.py and verify.py). Every configuration
#  is generated from scratch, so nothing but the scripts needs to be stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

CONFIGS=${CONFIGS:-"melting htpb hdpe pp melting-top hdpe-top melting-inert srm"}
status=0

for cfg in $CONFIGS; do
  python3 -B build_case.py "$cfg" > /dev/null || { echo "AblatingWall: FAIL - could not build '$cfg'"; exit 1; }

  rm -f OUTPUT/wall.tec OUTPUT/field.tec
  if ! ./MOSE.sh solve > OUTPUT/solve-"$cfg".log 2>&1; then
    echo "AblatingWall [$cfg]: FAIL - solver exited non-zero"
    tail -20 OUTPUT/solve-"$cfg".log
    status=1
    continue
  fi
  # The burning-grain BC is not a viscous wall, so it writes no wall file; every
  # other configuration must produce both.
  need_wall=1
  [[ "$cfg" == srm ]] && need_wall=0
  if [[ ! -f OUTPUT/field.tec ]] || [[ $need_wall -eq 1 && ! -f OUTPUT/wall.tec ]]; then
    echo "AblatingWall [$cfg]: FAIL - solver produced no wall/field output"
    tail -20 OUTPUT/solve-"$cfg".log
    status=1
    continue
  fi

  python3 -B verify.py || status=1
done

echo
if [[ $status -ne 0 ]]; then
  echo "AblatingWall: FAILED"
else
  echo "AblatingWall: PASSED"
fi
exit $status
