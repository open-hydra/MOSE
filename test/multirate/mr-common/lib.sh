#!/usr/bin/env bash
#===============================================================================
#  Common helpers of the multirate tests (test/multirate/mr-*): build the small cases with gen_mr_cases.py, run MOSE in a
#  separate directory per run, compare outputs byte by byte. Sourced by each run.sh (cwd = the test directory).
#  MOSE_BIN overrides the solver (default: bin/MOSE of this tree), MOSE_MMS_BIN the MMS driver.
#===============================================================================
MR_COMMON="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MR_MASTER="${MOSE_BIN:-$MR_COMMON/../../../bin/MOSE}"
MR_MMS="${MOSE_MMS_BIN:-$MR_COMMON/../../../bin/MOSE_mms}"
MR_AIR="$MR_COMMON/../../common/Air"
ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# mr_cases DIR CASE...  -> DIR/<case>/{input.ini,INPUT/*} (only if missing)
mr_cases() {
  local dir=$1; shift
  mkdir -p "$dir"
  local c
  for c in "$@"; do
    [[ -d "$dir/$c" ]] || python3 -B "$MR_COMMON/gen_mr_cases.py" --air "$MR_AIR" --out "$dir" "$c" > /dev/null || return 1
  done
}

# mr_prepare RUNDIR CASEDIR [ 'Section|key|value' ... ]  -> fresh run directory with the case linked and the keys set
# (the files of the case directory besides input.ini, e.g. the time tables p0.dat and q2d.dat of the inlets 402 and
# 410, are read by MOSE from the run directory: they are linked there too)
mr_prepare() {
  local run=$1 case=$2; shift 2
  rm -rf "$run"; mkdir -p "$run/INPUT" "$run/OUTPUT"
  local f
  for f in "$case"/INPUT/*; do ln -s "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" "$run/INPUT/$(basename "$f")"; done
  for f in "$case"/*; do
    [[ -f $f && $(basename "$f") != input.ini ]] && ln -s "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" "$run/$(basename "$f")"
  done
  cp "$case/input.ini" "$run/input.ini"
  python3 -B "$MR_COMMON/set_ini.py" "$run/input.ini" "$@"
}

# mr_threads_ok RUNDIR NT -> 0 if the log of the run says that MOSE ran on NT OpenMP threads (src/app/main.f90)
mr_threads_ok() {
  command grep -a -qE "Number of threads --> *$2\$" "$1/solve.log" || { echo "  run $1: not on $2 threads"; return 1; }
}

# mr_run RUNDIR CASEDIR [ 'Section|key|value' ... ]  -> runs MOSE in RUNDIR; 0 if it ends with exit 0, a field and the
# number of threads asked (OMP_NUM_THREADS)
mr_run() {
  local run=$1
  mr_prepare "$@" || return 1
  ( cd "$run" && "$MR_MASTER" > solve.log 2>&1 )
  local rc=$?
  [[ $rc -eq 0 && -f "$run/OUTPUT/field.vtm" ]] || { echo "  run $run: exit $rc, field.vtm $( [[ -f $run/OUTPUT/field.vtm ]] && echo present || echo missing)"; tail -3 "$run/solve.log"; return 1; }
  mr_threads_ok "$run" "$OMP_NUM_THREADS" || return 1
  return 0
}

# mr_same RUN_A RUN_B -> 0 if the two runs wrote the same output files (logs aside), none of the compared files is empty,
# and every .vts file of OUTPUT/vtk, field.vtm and the residual history are identical byte by byte (the other files of
# OUTPUT are compared by name). An empty or a missing file is never
# «identical»: a run that stopped early leaves an empty residual-history.dat (the runs compared must write the residual,
# res-diter = 1).
mr_same() {
  local a=$1 b=$2 f g n=0 la lb
  la=$(cd "$a/OUTPUT" 2>/dev/null && find . -type f ! -name '*.log' | sort)
  lb=$(cd "$b/OUTPUT" 2>/dev/null && find . -type f ! -name '*.log' | sort)
  [[ -n $la && $la == "$lb" ]] || { echo "  differ: the output files of $a and $b are not the same"; return 1; }
  for f in "$a"/OUTPUT/vtk/*.vts "$a"/OUTPUT/field.vtm "$a"/OUTPUT/residual-history.dat; do
    [[ -f $f ]] || continue
    g=$b/OUTPUT/${f#$a/OUTPUT/}
    [[ -s $f && -s $g ]] || { echo "  empty: ${f#$a/} or ${g#$b/}"; return 1; }
    cmp -s "$f" "$g" || { echo "  differ: ${f#$a/}"; return 1; }
    n=$((n + 1))
  done
  [[ $n -gt 0 ]] || return 1
  return 0
}
