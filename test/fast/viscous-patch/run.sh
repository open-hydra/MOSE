#!/usr/bin/env bash
#===============================================================================
#  Patch test of the viscous fluxes (fast tier)
#
#  Linear fields have zero viscous force and heat conduction: every cell must
#  give zero to round-off, inside the blocks and next to walls; next to inlets
#  and outlets, which have no viscous flux, minus the flux through them.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    ViscousLinearPatch (viscous-linear-patch): patch test of the viscous
#      fluxes
#    ViscousWallPatch (viscous-wall-patch): the same next to walls of
#      prescribed temperature
#    ViscousOpenPatch (viscous-open-patch): the same next to inlets and
#      outlets (401-408, 420)
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(viscous-linear-patch viscous-wall-patch viscous-open-patch)
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
