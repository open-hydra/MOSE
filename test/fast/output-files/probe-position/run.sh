#!/usr/bin/env bash
#===============================================================================
#  Probes given by position (fast tier)
#
#  A probe given by position records the cell whose centre is the nearest.
#  A planar 2-D duct with a different pressure in every cell (make_case.py)
#  runs two steps with three probes, each at the centre of a cell in the
#  mid-plane of the one-cell layer, and a solution at each step: the probe
#  rows of the last step must hold the values of their cells in the solution
#  of that step (field2), never those of a ghost cell.
#  Gas: test/common/Air.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"  INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat

cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 2
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
sol-diter = 1
sol-overwrite = false
INI

python3 -B make_case.py case || exit 1
rm -f OUTPUT/field*.tec OUTPUT/c*.txt
echo "-- three probes at cell centres of the mid-plane, two steps --"
"$LOCAL" > OUTPUT/run.log 2>&1
if [[ $? -ne 0 || ! -f OUTPUT/field2.tec ]]; then
  echo "  [FAIL] the run did not complete, or wrote no solution at its last step (see OUTPUT/run.log)"
  echo "ProbePosition: FAILED"
  exit 1
fi
if python3 -B make_case.py check OUTPUT/field2.tec; then
  echo "  [ok  ] each probe holds the values of its cell"
  echo "ProbePosition: all checks passed"
  exit 0
fi
echo "  [FAIL] the probes do not hold the values of their cells"
echo "ProbePosition: FAILED"
exit 1
