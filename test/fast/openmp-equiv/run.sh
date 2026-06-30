#!/usr/bin/env bash
#===============================================================================
#  Serial vs OpenMP equivalence (fast tier, I5)
#
#  Runs the same coarse Sod shock tube with OMP_NUM_THREADS=1 and =4 and asserts
#  the two solution files are bit-for-bit identical.  Catches threading bugs
#  (races, non-deterministic reductions, unguarded shared state).  Reuses the
#  coarse Sod inputs (test/1D/Sod79); no new mesh/IC data is stored here.
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

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x4.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x4.txt"  INPUT/bc.txt

cat > input.ini <<'EOF'
[MOSE-Parameters]
time-threshold = 0.10

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

run_threads() {
  local n="$1" dest="$2"
  rm -f OUTPUT/field.tec
  OMP_NUM_THREADS="$n" "$LOCAL" > "OUTPUT/solve_omp$n.log" 2>&1
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "OpenMP-equiv: FAIL — no field.tec for OMP=$n"; exit 1
  fi
  cp OUTPUT/field.tec "$dest"
}

run_threads 1 OUTPUT/field_omp1.tec
run_threads 4 OUTPUT/field_omp4.tec

if cmp -s OUTPUT/field_omp1.tec OUTPUT/field_omp4.tec; then
  echo "OpenMP-equiv: ok — 1 vs 4 threads bit-identical"
  exit 0
else
  echo "OpenMP-equiv: FAIL — 1 vs 4 threads differ"
  exit 1
fi
