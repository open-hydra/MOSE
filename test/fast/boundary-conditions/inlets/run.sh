#!/usr/bin/env bash
#===============================================================================
#  Inlet boundary conditions (fast tier)
#
#  1. Far-field inlet with the flow direction "normal". The flow direction of
#     an inlet record is an angle or the token normal, the face normal
#     (docs/user/boundary-conditions.md, common fields). A planar duct with an
#     ambient boundary (407) where the gas enters at 50 m/s runs 20 steps twice
#     (farfield_normal.py): with alpha = beta = normal, and with the angles the
#     token stands for on that face, alpha = beta = 0. The two runs must give
#     the same field to the bit, with no transverse velocity.
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

fail=0

echo "-- 1. far-field inlet with the flow direction normal --"
cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 20
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
INI

ffn=0
for dir in normal angle; do
  python3 -B farfield_normal.py case "$dir" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$dir.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] the run with the direction $dir did not complete (see OUTPUT/$dir.log)"
    ffn=1
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/field-$dir.tec"
  echo "   inlet direction $dir: largest |v| / u = $(python3 -B farfield_normal.py check "OUTPUT/field-$dir.tec")"
done

if [[ $ffn -eq 0 ]]; then
  if cmp -s OUTPUT/field-normal.tec OUTPUT/field-angle.tec; then
    echo "  [ok  ] normal and the equivalent angles give the same field to the bit"
  else
    echo "  [FAIL] normal and the equivalent angles give different fields"
    ffn=1
  fi
fi
[[ $ffn -eq 0 ]] || fail=1

if [[ $fail -ne 0 ]]; then
  echo "Inlets: FAILED"
  exit 1
fi
echo "Inlets: all checks passed"
exit 0
