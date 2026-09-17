#!/usr/bin/env bash
#===============================================================================
#  Closed-box mass conservation (fast tier, I4)
#
#  Evolves a shock-tube initial state inside a box closed by slip walls on every
#  face (BC type 300) and checks that total mass is conserved to round-off.
#  The closed-box BC file is derived from the coarse Sod BC (test/1D/Sod79) by
#  turning every extrapolation face (400) into a slip wall (300); the initial
#  condition is reused as-is.  No new mesh/IC data is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air

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
# Closed box: every extrapolation face (400) -> slip wall (300).
sed 's/400$/300/' "$SOD/INPUT/bc_x1.txt" > INPUT/bc.txt

cat > input.ini <<'EOF'
[MOSE-Parameters]
time-threshold = 0.10

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
EOF

rm -f OUTPUT/field.tec
if ! "$LOCAL" > OUTPUT/solve.log 2>&1; then
  echo "Conservation: FAIL — solver exited non-zero"
  tail -5 OUTPUT/solve.log
  exit 1
fi
if [[ ! -f OUTPUT/field.tec ]]; then
  echo "Conservation: FAIL — no OUTPUT/field.tec produced"
  exit 1
fi

python3 -B check_conservation.py
