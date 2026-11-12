#!/usr/bin/env bash
#===============================================================================
#  Scope of the preconditioning
#
#  integration-variables = prec (Weiss-Smith preconditioning) is for steady
#  (pseudo-time) runs, and the preconditioned flux HLLC-PC needs it. At set-up:
#    1. time-accurate + prec           -> "[ERROR] integration-variables=prec
#                                         ...", no solution written;
#    2. steady + HLLC-PC + cons        -> "[ERROR] Preconditioned HLL solver
#                                         ...", no solution written;
#    3. steady + HLLC-PC + prec        -> passes the set-up (its summary names
#                                         the preconditioned variables). The
#                                         Sod tube is no case for it (its
#                                         shocks stop the run): test/2D/euler/
#                                         low-mach-cylinder runs it.
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

# run_case <tag> <time-accurate> <integration variables> <riemann solver>
#   solver output in OUTPUT/<tag>.log; final field (if any) kept as OUTPUT/<tag>.tec
run_case() {
  local tag="$1" accurate="$2" vars="$3" rs="$4"
  {
    printf '[MOSE-Parameters]\n'
    if [[ "$accurate" == true ]]; then
      printf 'time-threshold = 0.04\n\n'
    else
      printf 'iter-threshold = 20\n\n'
    fi
    printf '[MOSE-Physics]\nequations = euler\n\n'
    printf '[MOSE-Numerics]\ncfl = 0.5\ntime-accurate = %s\ntime-scheme = RK3\n' "$accurate"
    printf 'integration-variables = %s\nspace-reconstruction = MUSCL\n' "$vars"
    printf 'flux-limiter = vanleer\nriemann-solver = %s\n' "$rs"
  } > input.ini
  rm -f OUTPUT/field.tec "OUTPUT/$tag.tec"
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  if [[ -f OUTPUT/field.tec ]]; then
    mv OUTPUT/field.tec "OUTPUT/$tag.tec"
  fi
}

# refused <tag> <message>: the set-up error and no solution written
refused() {
  [[ "$(< "OUTPUT/$1.log")" == *"$2"* && ! -f "OUTPUT/$1.tec" ]]
}
# past_setup <tag>: the set-up summary names the preconditioned variables, and neither refusal is in the log
past_setup() {
  local log
  log="$(< "OUTPUT/$1.log")"
  [[ "$log" =~ Integration\ variables\ +Preconditioned && "$log" != *"[ERROR] integration-variables=prec"* &&
     "$log" != *"[ERROR] Preconditioned HLL solver"* ]]
}

run_case accurate-prec true  prec HLLC
run_case pc-cons       false cons HLLC-PC
run_case steady-pc     false prec HLLC-PC

echo "-- the scope of integration-variables = prec and of HLLC-PC --"
refused accurate-prec "[ERROR] integration-variables=prec is currently supported only for steady"
check $? "time-accurate + prec: refused at set-up, no solution written"
refused pc-cons "[ERROR] Preconditioned HLL solver selected. integration-variables must be set to \"prec\""
check $? "HLLC-PC with cons: refused at set-up, no solution written"

echo "-- control --"
past_setup steady-pc
check $? "steady run, prec + HLLC-PC: passes the set-up"

if [[ $fail -ne 0 ]]; then
  echo "PrecScope: FAILED"
  exit 1
fi
echo "PrecScope: all checks passed"
exit 0
