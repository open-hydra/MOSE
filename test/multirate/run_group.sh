#!/usr/bin/env bash
#===============================================================================
#  Runs a group of the checks of the multirate integrator, one after the other, each in its own directory
#  (mr-<check>/run.sh). The validation and needs-mpi entries of the multirate integrator in test/CMakeLists.txt
#  (MultirateIdentity, MultirateConservation, MultirateOrder, MultirateOrderInletsMoving, MultirateShocks,
#  MultirateInput, MultirateHybrid) are groups of these checks; see test/TESTING-PLAN.md and docs/vv/multirate.md.
#  Every check of the group runs, also after a failure; the group fails if one of them fails, and its last line names
#  the checks that failed. Each check prints its own lines and its time in seconds.
#  Usage: ./run_group.sh GROUP CHECK...      e.g.  ./run_group.sh input reject fixedstep trace
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
[[ $# -ge 2 ]] || { echo "usage: $0 GROUP CHECK..."; exit 2; }
group=$1; shift
failed=()
for c in "$@"; do
  if [[ ! -f mr-$c/run.sh ]]; then
    echo "$group: no check mr-$c"; failed+=("$c"); continue
  fi
  t0=$(date +%s)
  bash "mr-$c/run.sh" || failed+=("$c")
  echo "  ($group: mr-$c took $(( $(date +%s) - t0 )) s)"
done
if [[ ${#failed[@]} -eq 0 ]]; then
  echo "$group: ok ($*)"
  exit 0
fi
echo "$group: FAIL (${failed[*]})"
exit 1
