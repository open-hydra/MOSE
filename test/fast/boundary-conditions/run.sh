#!/usr/bin/env bash
#===============================================================================
#  Boundary conditions (fast tier)
#
#  Boundary conditions: inlets, the choked nozzle, the records of bc.txt,
#  walls, walls in a rotating frame, and extrapolation faces and the gas
#  entering through them.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    Inlets (inlets): inlet boundary conditions
#    Nozzle (nozzle): the face of a choked nozzle (420) below and above
#      p_shock
#    BCInput (bc-input): property lines of bc.txt
#    RotatingFrame (rotating-frame): stationary walls of a rotating frame
#    WallHeatSign (wall-heat-sign): sign of the wall heat flux and wall
#      stress
#    Extrapolation400 (extrapolation-400): extrapolation faces (400) and the
#      lines of the gas entering through them
#    SupersonicInlet (supersonic-inlet): the supersonic inlet (405) with the
#      flow direction normal
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(inlets nozzle bc-input rotating-frame wall-heat-sign extrapolation-400 supersonic-inlet)
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
