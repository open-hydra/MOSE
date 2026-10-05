#!/usr/bin/env bash
#===============================================================================
#  Spalart-Allmaras with variable density (fast tier)
#
#  A constructed periodic field at rest with variable temperature and nu~ and
#  no wall (make_case.py): one explicit Euler step with local time stepping
#  gives, in OUTPUT/diagnostic.tec, the right-hand side of the equation for
#  rho nu~ that the solver evaluates. It must converge to the exact right-hand
#  side of the standard conservation form (NASA Turbulence Modeling Resource,
#  Spalart-Allmaras page) with second order, on 32 x 32 and 64 x 64 cells.
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../../common/Air-viscous

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
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = SA
schmidt = 0.7

[MOSE-Turbulence]
point-implicit = false

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF

errors=()
for n in 32 64; do
  python3 -B make_case.py case "$n" || exit 1
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/n$n.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/n$n.log"; then
    echo "  [FAIL] the run on $n x $n cells did not complete (see OUTPUT/n$n.log)"
    echo "MMSResiduals: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-n$n.tec"
  errors+=("$(python3 -B make_case.py error "OUTPUT/diagnostic-n$n.tec")")
done

echo "-- right-hand side of rho nu~ against the standard conservation form --"
if python3 -B make_case.py check "${errors[0]}" "${errors[1]}"; then
  echo "  [ok  ] second-order agreement"
  echo "MMSResiduals: all checks passed"
  exit 0
fi
echo "  [FAIL] the solver does not evaluate the standard form"
echo "MMSResiduals: FAILED"
exit 1
