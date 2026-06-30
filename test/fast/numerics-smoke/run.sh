#!/usr/bin/env bash
#===============================================================================
#  Numerics smoke-matrix (fast tier, N4 + scheme/recon sweeps)
#
#  Sweeps the numerical-method options one axis at a time on a short coarse Sod
#  shock tube and asserts density/pressure stay positive and finite.  Guarantees
#  that every flux limiter, time scheme, reconstruction and integration-variable
#  path at least executes without producing NaN/p<0.  Reuses the coarse Sod
#  inputs (test/1D/Sod79) and the positivity checker from ../riemann-smoke.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air
CHECK=../riemann-smoke/check_positivity.py

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

# run_case <label> <recon> <limiter> <time-scheme> <int-vars>
run_case() {
  local label="$1" recon="$2" limiter="$3" tscheme="$4" ivars="$5"
  cat > input.ini <<EOF
[MOSE-Parameters]
time-threshold = 0.04

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = $tscheme
integration-variables = $ivars
space-reconstruction = $recon
flux-limiter = $limiter
riemann-solver = HLLC
EOF
  rm -f OUTPUT/field.tec
  if ! "$LOCAL" > "OUTPUT/solve_${label//[^A-Za-z0-9]/_}.log" 2>&1; then
    echo "  [$label] FAIL: solver exited non-zero"; fail=1; return
  fi
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [$label] FAIL: no OUTPUT/field.tec produced"; fail=1; return
  fi
  python3 -B "$CHECK" "$label" || fail=1
}

echo "-- flux limiters (MUSCL, RK3, cons) --"
for lim in vanleer vanalbada minmod superbee mc; do
  run_case "limiter:$lim" MUSCL "$lim" RK3 cons
done

echo "-- time schemes (MUSCL, vanleer, cons) --"
for ts in euler RK2 RK3; do
  run_case "time:$ts" MUSCL vanleer "$ts" cons
done

echo "-- reconstruction (HLLC, RK3, cons) --"
run_case "recon:first-order" first-order vanleer RK3 cons
run_case "recon:MUSCL"       MUSCL       vanleer RK3 cons

echo "-- integration variables (MUSCL, RK3) --"
run_case "ivars:cons" MUSCL vanleer RK3 cons
run_case "ivars:prim" MUSCL vanleer RK3 prim

if [[ $fail -ne 0 ]]; then
  echo "NumericsSmoke: FAILED"
  exit 1
fi
echo "NumericsSmoke: all configurations passed"
exit 0
