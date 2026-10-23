#!/usr/bin/env bash
#===============================================================================
#  One flux on the faces between blocks, across MPI ranks (needs-mpi tier)
#
#  The two-block periodic cases of test/fast/connections/connection-flux (shock tube; low
#  Mach number wave for HLLC-PC), with each block on its own rank: the block
#  cut and the periodic connection are faces between ranks. With one face
#  flux, identical on the two sides:
#    - two ranks give one rank bit for bit, and the one-block run;
#    - mass, momentum and energy are conserved to round-off (time-accurate runs).
#  Same Riemann solvers as the fast test. Gas: test/common/Air.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../common/Air
CASE=../../fast/connections/connection-flux/make_case.py
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
  local accurate=true vars=cons extra=""
  if [[ $1 == HLLC-PC ]]; then accurate=false; vars=prec; fi
  if [[ $1 == HLLC+Chen ]]; then extra="riemann-options-Mco = 0.2"; fi
  cat > input.ini <<EOT
[MOSE-Parameters]
iter-threshold = 40
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = $accurate
time-scheme = RK3
integration-variables = $vars
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = $1
$extra

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
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
  # MOSE's own error lines only ([ERROR] lines, the stop of a diverging run, the stops of every rank, the runtime errors
  # of the compilers): the log also holds the output of mpiexec, whose lines (ORTE_ERROR_LOG, ...) are not MOSE's.
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|\[ERROR\]|ERROR STOP|STOP: |ABORT: |MPI error code|Fortran runtime error|Program received signal|forrtl: ' "runs/$2.log"; then
    echo "  [FAIL] the run $2 did not complete (see runs/$2.log)"; return 1
  fi
  mv OUTPUT "runs/$2"
}

fail=0
for rs in HLLC HLLC+Tramel HLLC+Chen HLLC-PC; do
  echo "-- $rs --"
  write_ini "$rs"
  kind=sod; [[ $rs == HLLC-PC ]] && kind=wave
  python3 -B "$CASE" case one "$kind" || exit 1
  run 1 "$rs-one" || { fail=1; continue; }
  python3 -B "$CASE" case two "$kind" || exit 1
  run 1 "$rs-two-np1" || { fail=1; continue; }
  run 2 "$rs-two-np2" || { fail=1; continue; }
  echo "  two blocks, 2 ranks against 1 rank:"
  python3 -B "$CASE" same "runs/$rs-two-np1" "runs/$rs-two-np2" || fail=1
  echo "  two blocks on 2 ranks against one block:"
  python3 -B "$CASE" same "runs/$rs-one" "runs/$rs-two-np2" || fail=1
  if [[ $rs != HLLC-PC ]]; then
    echo "  conservation, two blocks on 2 ranks:"
    python3 -B "$CASE" conserved "runs/$rs-two-np2" || fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "MPIConnectionFlux: FAILED"
  exit 1
fi
echo "MPIConnectionFlux: all checks passed"
exit 0
