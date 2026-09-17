#!/usr/bin/env bash
#===============================================================================
#  Restart round-trip (fast tier, I3/F6)
#
#  Verifies that stopping a run, writing a solution and restarting (newrun=false,
#  which reads OUTPUT/field.tec) reproduces the uninterrupted run:
#     A : N iterations in one go
#     B : N/2 iterations -> restart -> N/2 more
#  A and B must match to solution-file round-trip precision.  Reuses the coarse
#  Sod inputs (test/1D/Sod79); no new mesh/IC data is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air
N=50          # total iterations
HALF=25       # iterations before the restart

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt

write_ini() {  # <iter-threshold> <newrun>
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = $2

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.9
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
}

# --- A: uninterrupted N-iteration run --------------------------------------
write_ini "$N" true
rm -f OUTPUT/field*.tec
"$LOCAL" > OUTPUT/solve_full.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from full run"; exit 1; }
cp OUTPUT/field.tec field_full.tec

# --- B: N/2 then restart for the remaining N/2 -----------------------------
write_ini "$HALF" true
rm -f OUTPUT/field*.tec
"$LOCAL" > OUTPUT/solve_half.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from first half"; exit 1; }

write_ini "$HALF" false
"$LOCAL" > OUTPUT/solve_restart.log 2>&1
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from restart"; exit 1; }
cp OUTPUT/field.tec field_restart.tec

python3 -B check_restart.py field_full.tec field_restart.tec
