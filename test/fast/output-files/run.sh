#!/usr/bin/env bash
#===============================================================================
#  Output files (fast tier)
#
#  The output files: wall, field, diagnostic and probe outputs, their names,
#  formats, times and values.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    WallOutput (wall-output): wall output files
#    VtkNames (vtk-names): names and times of the VTK outputs
#    VtkRawFormat (vtk-raw-format): vtk raw solutions as raw appended data
#    WallOutputResidual (wall-output-residual): wall output and residual
#    ProbePosition (probe-position): probes given by position
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(wall-output vtk-names vtk-raw-format wall-output-residual probe-position)
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
