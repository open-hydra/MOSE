#!/usr/bin/env bash
#===============================================================================
#  Boundary conditions (fast tier)
#
#  Boundary conditions: inlets, the records of bc.txt, walls, and walls in a
#  rotating frame.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    Inlets (inlets): inlet boundary conditions
#    BCInput (bc-input): property lines of bc.txt
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(inlets bc-input)
if [[ $# -gt 0 ]]; then
  PARTS=("$@")
fi

status=0
for part in "${PARTS[@]}"; do
  echo "=== $part"
  if (cd "$HERE/$part" && ./run.sh); then
    echo "=== $part: PASS"
  else
    echo "=== $part: FAIL"
    status=1
  fi
done
exit $status
