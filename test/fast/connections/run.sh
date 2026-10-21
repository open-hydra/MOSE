#!/usr/bin/env bash
#===============================================================================
#  Faces between blocks (fast tier)
#
#  The faces between blocks: connections (101, 201), chimera faces (102) and
#  rotational periodicity (200), with their inviscid and viscous fluxes.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    ChimeraFace (chimera-face): viscous stencil of chimera faces at the
#      edges of a block face
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(chimera-face)
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
