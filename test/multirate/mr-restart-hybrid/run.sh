#!/usr/bin/env bash
#===============================================================================
#  Multirate restart with MPI processes and OpenMP threads: mr2d-box-4b (Euler,
#  2 x 2 blocks) and mr2d-plate-4b (laminar plate in 2 x 2 blocks, a vertex of four blocks in the boundary layer), multirate
#  (max-level 3, safety 2, levels recomputed at every macro step), solutions every 30 macro steps, residual at every
#  macro step, on 2 processes x 2 threads: 100 macro steps in one go against 50 + restart (newrun = false, vtk raw) + 50,
#  OUTPUT trees identical byte by byte (test/checks/restart/check_restart_vtk.py), on the threads asked.
#  The launcher is the one CMake found (MOSE_MPIEXEC, MOSE_MPIEXEC_NUMPROC_FLAG; mpirun -n by default).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
CHECK=../../checks/restart/check_restart_vtk.py
rm -rf cases runs
MPIEXEC=${MOSE_MPIEXEC:-mpirun}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}
command -v "$MPIEXEC" > /dev/null || { echo "mr-restart-hybrid: FAIL - MPI launcher $MPIEXEC not found"; exit 1; }
mr_cases cases mr2d-box-4b mr2d-plate-4b || { echo "mr-restart-hybrid: FAIL - cannot build the cases"; exit 1; }
status=0
KEYS=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0"
      "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" "MOSE-IO|res-diter|1"
      "MOSE-IO|sol-overwrite|false" "MOSE-IO|sol-diter|30")

# part DIR -> one run on 2 processes x 2 threads, appending to solve.log
part() {
  ( cd "$1" && OMP_NUM_THREADS=2 "$MPIEXEC" "$NPFLAG" 2 "$MR_MASTER" >> solve.log 2>&1 ) \
    || { echo "  run $1: exit non-zero"; tail -3 "$1/solve.log"; return 1; }
  mr_threads_ok "$1" 2
}

for c in mr2d-box-4b mr2d-plate-4b; do
  d=runs/$c
  mr_prepare "$d/full" "cases/$c" "${KEYS[@]}" "MOSE-Parameters|iter-threshold|100" && part "$d/full" || { status=1; continue; }
  mr_prepare "$d/split" "cases/$c" "${KEYS[@]}" "MOSE-Parameters|iter-threshold|50" && part "$d/split" || { status=1; continue; }
  python3 -B "$MR_COMMON/set_ini.py" "$d/split/input.ini" "MOSE-Parameters|newrun|false" && part "$d/split" \
    || { status=1; continue; }
  if [[ -s $d/full/OUTPUT/residual-history.dat && -s $d/split/OUTPUT/residual-history.dat ]] \
     && python3 -B "$CHECK" "$d/full/OUTPUT" "$d/split/OUTPUT" > "$d/compare.log" 2>&1; then
    echo "  $c 2x2: 100 = 50 + restart + 50, identical ($(cd "$d/full/OUTPUT" && find . -type f ! -name '*.log' | wc -l) files)"
  else
    echo "  $c 2x2: DIFFERENT or empty"; sed 's/^/    /' "$d/compare.log" 2>/dev/null | head -6; status=1
  fi
done
if [[ $status -eq 0 ]]; then echo "mr-restart-hybrid: ok"; else echo "mr-restart-hybrid: FAIL"; fi
exit $status
