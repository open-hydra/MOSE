#!/usr/bin/env bash
#===============================================================================
#  Wall output and residual (fast tier)
#
#  The wall output must not change the residual. A short laminar hypersonic
#  plate segment (make_case.py), RK3, local time stepping, is converged in
#  4000 iterations at cfl 0.8, vnn 0.1, with an output every 400; the wall
#  shear must change by less than 1e-11 over the last 400. The wall output is
#  written just before the diagnostic one: in OUTPUT/diagnostic.tec the
#  residual of every wall cell must then be round-off small too, below 1e-6
#  of the wall shear force (rho u) and of the wall heat flow (rho e) of its
#  face.
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

python3 -B make_case.py case || exit 1
cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 4000
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.8
vnn = 0.1
cfl-rise-threshold = 500
time-accurate = false
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-diter = 400
sol-overwrite = false
shell-diter = 400
res-diter = 400
EOF

echo "-- laminar plate segment: steady state, then the residual of the wall cells --"
rm -f OUTPUT/field*.tec OUTPUT/wall*.tec OUTPUT/diagnostic*.tec
"$LOCAL" > OUTPUT/run.log 2>&1
if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' OUTPUT/run.log; then
  echo "  [FAIL] the run blew up or wrote no diagnostic output (see OUTPUT/run.log)"
  echo "WallOutputResidual: FAILED"
  exit 1
fi
done_steps=$(grep -a -oE 'Iter = +[0-9]+' OUTPUT/run.log | tail -1 | grep -a -oE '[0-9]+$')
if [[ "$done_steps" != "4000" ]]; then
  echo "  [FAIL] the run stopped after ${done_steps:-0} of 4000 iterations (see OUTPUT/run.log)"
  echo "WallOutputResidual: FAILED"
  exit 1
fi

fail=0
if ! python3 -B make_case.py converged OUTPUT/wall9.tec OUTPUT/wall.tec; then
  echo "  [FAIL] the run did not converge: the check would not be meaningful"
  fail=1
elif python3 -B make_case.py residual OUTPUT/diagnostic.tec OUTPUT/wall.tec OUTPUT/field.tec; then
  echo "  [ok  ] the wall output leaves the residual of the wall cells alone"
else
  echo "  [FAIL] the residual of the wall cells holds the wall fluxes a second time"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "WallOutputResidual: FAILED"
  exit 1
fi
echo "WallOutputResidual: all checks passed"
exit 0
