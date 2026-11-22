#!/usr/bin/env bash
#===============================================================================
#  The manifold (BC 501) across MPI ranks (needs-mpi tier)
#
#  The two-block manifold cases of test/checks/boundary-conditions/manifold
#  (straight and turned back, from the exact state and from a perturbed one),
#  with each block on its own rank: the two sides of the manifold sit on two
#  ranks and the sums of each side come from one rank each, so two ranks must
#  give one rank bit for bit (one reduction per stage adds the zeros of the
#  other rank). Gas: test/common/Air.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../common/Air
CASE=../../checks/boundary-conditions/manifold/make_case.py
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}
ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=1
mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat
ln -sf "../$AIR/phase.txt"  INPUT/phase.txt
ln -sf "../$AIR/thermo.dat" INPUT/thermo.dat

write_ini() {
  cat > input.ini <<EOT
[MOSE-Parameters]
iter-threshold = $1
newrun = true
res-threshold = $2

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = $3
ini-diter = 1000000000
EOT
}

# run <ranks> <tag>: solve the case in INPUT/ on <ranks> ranks, keep OUTPUT as runs/<tag>
run() {
  rm -rf OUTPUT "runs/$2"
  mkdir -p OUTPUT/vtk
  if ! "$MPIEXEC" "$NPFLAG" "$1" "$LOCAL" > "runs/$2.log" 2>&1; then
    echo "  [FAIL] the run $2 exited non-zero (see runs/$2.log)"; return 1
  fi
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|\[ERROR\]|ERROR STOP|STOP: |ABORT: |MPI error code|Fortran runtime error|Program received signal|forrtl: ' "runs/$2.log"; then
    echo "  [FAIL] the run $2 did not complete (see runs/$2.log)"; return 1
  fi
  mv OUTPUT "runs/$2"
  [[ $2 == *reverse* ]] && touch "runs/$2/reverse"
  return 0
}

fail=0
for start in exact perturbed; do
  iters=500; res=1e-300; shell=1; extra=""
  [[ $start == perturbed ]] && { iters=60000; res=1e-13; shell=500; extra=perturbed; }
  write_ini $iters $res $shell
  for layout in manifold reverse; do
    echo "-- $layout, $start start --"
    python3 -B "$CASE" case "$layout" $extra || exit 1
    run 1 "$start-$layout-np1" || { fail=1; continue; }
    run 2 "$start-$layout-np2" || { fail=1; continue; }
    echo "  2 ranks against 1 rank:"
    python3 -B "$CASE" same "runs/$start-$layout-np1" "runs/$start-$layout-np2" || fail=1
    echo "  the line of the manifold on 2 ranks:"
    python3 -B "$CASE" carried "runs/$start-$layout-np2.log" "runs/$start-$layout-np2" 1e-12 || fail=1
  done
done

if [[ $fail -ne 0 ]]; then
  echo "MPIManifold: FAILED"
  exit 1
fi
echo "MPIManifold: all checks passed"
exit 0
