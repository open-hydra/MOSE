#!/usr/bin/env bash
#===============================================================================
#  Implicit residual smoothing with time-accurate runs (fast tier)
#
#  Implicit residual smoothing (IRS) is a convergence accelerator for steady
#  (pseudo-time) runs: it replaces the residual by a smoothed one, so the
#  computed transient is no longer that of the space discretization. With
#  time-accurate = true it must be refused at set-up:
#    1. time-accurate + irs-beta > 0  -> "[ERROR] Implicit residual
#                                        smoothing ...", no solution written;
#    2. time-accurate + irs = true    -> the same;
#    3. steady run + irs-beta > 0     -> runs (IRS stays available there);
#    4. time-accurate without IRS     -> runs.
#  Short coarse Sod shock tube; inputs shared with test/1D/Sod79.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79
COMMON=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"    INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"   INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec" INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt" INPUT/bc.txt

fail=0

# check <status> <description>: status 0 = pass
check() {
  if [[ $1 -eq 0 ]]; then
    echo "  [ok  ] $2"
  else
    echo "  [FAIL] $2"
    fail=1
  fi
}

# run_case <tag> <time-accurate> <extra Numerics lines...>
#   solver output in OUTPUT/<tag>.log; final field (if any) kept as OUTPUT/<tag>.tec
run_case() {
  local tag="$1" accurate="$2"
  shift 2
  {
    printf '[MOSE-Parameters]\n'
    if [[ "$accurate" == true ]]; then
      printf 'time-threshold = 0.04\n\n'
    else
      printf 'iter-threshold = 20\n\n'
    fi
    printf '[MOSE-Physics]\nequations = euler\n\n'
    printf '[MOSE-Numerics]\ncfl = 0.5\ntime-accurate = %s\ntime-scheme = RK3\n' "$accurate"
    printf 'integration-variables = cons\nspace-reconstruction = MUSCL\n'
    printf 'flux-limiter = vanleer\nriemann-solver = HLLC\n'
    local line
    for line in "$@"; do
      printf '%s\n' "$line"
    done
  } > input.ini
  rm -f OUTPUT/field.tec "OUTPUT/$tag.tec"
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  if [[ -f OUTPUT/field.tec ]]; then
    mv OUTPUT/field.tec "OUTPUT/$tag.tec"
  fi
}

# refused <tag>: the set-up error on residual smoothing and no solution written
refused() {
  [[ "$(< "OUTPUT/$1.log")" == *"[ERROR] Implicit residual smoothing"* && ! -f "OUTPUT/$1.tec" ]]
}
# ran <tag>: a solution was written and no [ERROR] printed
ran() {
  [[ -f "OUTPUT/$1.tec" && "$(< "OUTPUT/$1.log")" != *"[ERROR]"* ]]
}

run_case accurate-beta true  "irs-beta = 0.5"
run_case accurate-irs  true  "irs = true"
run_case steady-beta   false "irs-beta = 0.5"
run_case accurate      true

echo "-- time-accurate runs with implicit residual smoothing --"
refused accurate-beta
check $? "time-accurate + irs-beta = 0.5: refused at set-up, no solution written"
refused accurate-irs
check $? "time-accurate + irs = true: refused at set-up, no solution written"

echo "-- controls --"
ran steady-beta
check $? "steady run + irs-beta = 0.5: runs"
ran accurate
check $? "time-accurate run without IRS: runs"

if [[ $fail -ne 0 ]]; then
  echo "IRSTimeAccurate: FAILED"
  exit 1
fi
echo "IRSTimeAccurate: all checks passed"
exit 0
