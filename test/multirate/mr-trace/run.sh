#!/usr/bin/env bash
#===============================================================================
#  Stage sequence of the multirate integrator: 1D acoustic wave on the
#  banded mesh with levels 0..3 (safety 2, max-level 3), trace of the first macro step: level, stage, abscissa and
#  theta of every coarser level of each of the 30 stages equal to the reference sequence of the recursion
#  (bench_sequence_L3.txt).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu || { echo "mr-trace: FAIL - cannot build the case"; exit 1; }
status=0
mr_run runs/mr cases/mr1d-acu "MOSE-Parameters|iter-threshold|1" "MOSE-Multirate|enabled|.true." \
       "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|trace|.true." || status=1
[[ $status -eq 0 ]] && { python3 -B verify_trace.py runs/mr/solve.log bench_sequence_L3.txt || status=1; }
if [[ $status -eq 0 ]]; then echo "mr-trace: ok"; else echo "mr-trace: FAIL"; fi
exit $status
