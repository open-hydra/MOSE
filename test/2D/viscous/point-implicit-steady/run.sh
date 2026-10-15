#!/usr/bin/env bash
#===============================================================================
#  Steady state with the point-implicit factor (validation tier)
#
#  The steady state must not depend on the time step. A short hypersonic
#  plate segment (make_case.py) with Spalart-Allmaras and with SST, RK3, local
#  time stepping, point-implicit factor on:
#    1. 10000 iterations at cfl 0.8, vnn 0.4, until the wall shear changes by
#       less than 1e-11 over the last 2000;
#    2. 2000 more iterations from that state at cfl 0.4, vnn 0.2;
#    3. 2000 more iterations from that state with implicit residual smoothing
#       (irs-beta = 0.5), which smooths the turbulence residual too.
#  The wall shear and heat flux of runs 2 and 3 must agree with those of run 1
#  to 1e-9.
#  Gases: test/common/Air-viscous, viscosity of test/common/Air-sutherland.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../common/Air-viscous
SUTH=../../../common/Air-sutherland

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SUTH/transport.dat" INPUT/transport.dat

fail=0

# solve <model> <cfl> <vnn> <iterations> <ramp> <output every> <tag> [numerics line]: steady run; the solver output goes to
# OUTPUT/<tag>.log
solve() {
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $4
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = $1
schmidt = 0.7

[MOSE-Turbulence]
point-implicit = .true.
omega-wall-bc = practical

[MOSE-Numerics]
cfl = $2
vnn = $3
cfl-rise-threshold = $5
time-accurate = false
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
${8:-}

[MOSE-IO]
sol-diter = $6
sol-overwrite = false
shell-diter = $6
res-diter = $6
EOF
  rm -f OUTPUT/field*.tec OUTPUT/wall*.tec OUTPUT/diagnostic*.tec
  "$LOCAL" > "OUTPUT/$7.log" 2>&1
  if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$7.log"; then
    echo "  [FAIL] the run blew up (see OUTPUT/$7.log)"
    fail=1
    return 1
  fi
  local done_steps
  done_steps=$(grep -a -oE 'Iter = +[0-9]+' "OUTPUT/$7.log" | tail -1 | grep -a -oE '[0-9]+$')
  if [[ "$done_steps" != "$4" ]]; then
    echo "  [FAIL] the run stopped after ${done_steps:-0} of $4 iterations (see OUTPUT/$7.log)"
    fail=1
    return 1
  fi
}

for model in SA SST; do
  echo "-- $model: steady state at cfl 0.8, vnn 0.4, then continued at cfl 0.4, vnn 0.2 --"
  python3 -B make_case.py case "$model" || exit 1
  solve "$model" 0.8 0.4 10000 500 2000 "$model-A" || continue
  cp OUTPUT/wall.tec "OUTPUT/$model-A-wall.tec"
  if ! python3 -B make_case.py converged OUTPUT/wall4.tec OUTPUT/wall.tec; then
    echo "  [FAIL] the first run did not converge: the comparison would not be meaningful"
    fail=1
    continue
  fi
  rm -f INPUT/ic.tec
  cp OUTPUT/field.tec INPUT/ic.tec
  solve "$model" 0.4 0.2 2000 0 2000 "$model-B" || continue
  cp OUTPUT/wall.tec "OUTPUT/$model-B-wall.tec"
  if python3 -B make_case.py compare "OUTPUT/$model-A-wall.tec" "OUTPUT/$model-B-wall.tec"; then
    echo "  [ok  ] the steady state does not depend on the time step"
  else
    echo "  [FAIL] the steady state depends on the time step"
    fail=1
  fi
  solve "$model" 0.8 0.4 2000 0 2000 "$model-C" "irs-beta = 0.5" || continue
  cp OUTPUT/wall.tec "OUTPUT/$model-C-wall.tec"
  if python3 -B make_case.py compare "OUTPUT/$model-A-wall.tec" "OUTPUT/$model-C-wall.tec"; then
    echo "  [ok  ] the steady state does not change with implicit residual smoothing"
  else
    echo "  [FAIL] the steady state changes with implicit residual smoothing"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "PointImplicitSteady: FAILED"
  exit 1
fi
echo "PointImplicitSteady: all checks passed"
exit 0
