#!/usr/bin/env bash
#===============================================================================
#  Reynolds-stress model -- decay of isotropic turbulence in a uniform stream
#
#  Runs the case once per SSG/LRR variant and checks each result against the
#  analytic solution (see build_case.py and verify.py). Every configuration is
#  generated from scratch, so nothing but the scripts needs to be stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MODELS=${MODELS:-"SSGLRR SSGLRR-SD SSGLRR-w2019 SSGLRR-w2019-SD"}
status=0

for model in $MODELS; do
  python3 -B build_case.py "$model" > /dev/null || { echo "RSMDecay: FAIL - could not build '$model'"; exit 1; }

  rm -f OUTPUT/field.tec
  if ! ./MOSE.sh solve > OUTPUT/solve-"$model".log 2>&1; then
    echo "RSMDecay [$model]: FAIL - solver exited non-zero"
    tail -20 OUTPUT/solve-"$model".log
    status=1
    continue
  fi
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "RSMDecay [$model]: FAIL - solver produced no field output"
    tail -20 OUTPUT/solve-"$model".log
    status=1
    continue
  fi

  python3 -B verify.py || status=1
done

echo
if [[ $status -ne 0 ]]; then
  echo "RSMDecay: FAILED"
else
  echo "RSMDecay: PASSED"
fi
exit $status
