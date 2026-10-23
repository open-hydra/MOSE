#!/usr/bin/env bash
#===============================================================================
#  Serial vs OpenMP equivalence (fast tier, I5)
#
#  Runs the same coarse Sod shock tube with OMP_NUM_THREADS=1 and =4 and asserts
#  the two solution files are bit-for-bit identical.  Catches threading bugs
#  (races, non-deterministic reductions, unguarded shared state).  Reuses the
#  coarse Sod inputs (test/1D/Sod79); no new mesh/IC data is stored here.
#
#  Chimera: a block of 40 x 40 cells overset on a background block of 48 x 48
#  (make_case.py), whose four side faces are chimera boundaries meeting at its
#  edges, air of test/common/Air-viscous in a smooth state, Navier-Stokes, 10
#  RK3 steps, with OMP_NUM_THREADS=1 and =4: the two solution files must be
#  bit-for-bit identical. Where two chimera faces of a block meet, both write
#  the ghost cell of their common edge and read it back for the viscous
#  stencil of their face; in one parallel loop those would be unordered writes
#  of one location by different threads.
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
rm -f INPUT/*
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

status=0
if cmp -s OUTPUT/field_omp1.tec OUTPUT/field_omp4.tec; then
  echo "OpenMP-equiv: ok — 1 vs 4 threads bit-identical"
else
  echo "OpenMP-equiv: FAIL — 1 vs 4 threads differ"
  status=1
fi

# --- chimera: four chimera faces of a block meeting at its edges ------------
rm -f INPUT/*
for f in phase.txt thermo.dat transport.dat; do
  ln -sf "../../../common/Air-viscous/$f" "INPUT/$f"
done
python3 -B make_case.py chimera || exit 1
cat > input.ini <<'EOF'
[MOSE-Parameters]
iter-threshold = 10
newrun = true

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
EOF
run_threads 1 OUTPUT/chimera_omp1.tec
run_threads 4 OUTPUT/chimera_omp4.tec
if cmp -s OUTPUT/chimera_omp1.tec OUTPUT/chimera_omp4.tec; then
  echo "OpenMP-equiv: ok — chimera, 1 vs 4 threads bit-identical"
else
  echo "OpenMP-equiv: FAIL — chimera, 1 vs 4 threads differ"
  status=1
fi
exit $status
