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
#    ConnectionFlux (connection-flux): one flux per face between blocks
#    ConnectionOrder (connection-order): faces between blocks, whatever the
#      order of bc.txt
#    SectorRadial (sector-radial): rotational periodicity of a 3-D sector
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(chimera-face connection-flux connection-order sector-radial)
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
