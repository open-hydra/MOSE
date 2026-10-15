#!/usr/bin/env bash
#===============================================================================
#  Time integration (fast tier)
#
#  The explicit step and the time of a run: the stages and their boundary data,
#  time tables, the stop rules, the local time step and source terms.
#
#  One part per folder, each with its own case and checks, described at
#  the head of its run.sh; each part runs in its folder, as CTest runs a
#  test:
#    TimeControl (time-control): boundary data in time, time tables and stop
#      rules
#    LocalTimeStep (local-time-step): stability of the local time step
#    TurbulenceSources (turbulence-sources): turbulence sources in
#      time-accurate runs
#  ./run.sh runs every part; ./run.sh <part> ... runs those parts only.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PARTS=(time-control local-time-step turbulence-sources)
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
