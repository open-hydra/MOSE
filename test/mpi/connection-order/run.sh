#!/usr/bin/env bash
#===============================================================================
#  One flux per face between blocks across MPI ranks, whatever the order of
#  bc.txt (needs-mpi tier)
#
#  The 2-D case of test/fast/connections/connection-order: two blocks cut at y = 0.5 with
#  their index directions permuted across the cut (face 4 of one block with
#  face 2 of the other, two faces of the same parity), Navier-Stokes with
#  HLLC+Tramel and a shock leaving the cut. With each block on its own rank
#  the cut is a face between ranks, with its inviscid and viscous shares sent
#  across. Two ranks with the records of bc.txt sorted and in reverse order
#  must give the solution of one rank with the sorted file, bit for bit.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../common/Air-viscous
CASE=../../fast/connections/connection-order/make_case.py
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=1

mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

cat > input.ini <<'EOT'
[MOSE-Parameters]
iter-threshold = 20
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC+Tramel

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
EOT

# run <ranks> <order>: solve with the bc.txt records in that order, keep OUTPUT as runs/<order>-np<ranks>
run() {
  local label="$2-np$1"
  python3 -B "$CASE" case "$2" > /dev/null || exit 1
  rm -rf OUTPUT "runs/$label"
  mkdir -p OUTPUT/vtk
  "$MPIEXEC" "$NPFLAG" "$1" "$LOCAL" > "runs/$label.log" 2>&1
  # MOSE's own error lines only ([ERROR] lines, the stop of a diverging run, the stops of every rank, the runtime errors
  # of the compilers): the log also holds the output of mpiexec, whose lines (ORTE_ERROR_LOG, ...) are not MOSE's.
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|\[ERROR\]|ERROR STOP|STOP: |ABORT: |MPI error code|Fortran runtime error|Program received signal|forrtl: ' "runs/$label.log"; then
    echo "  [FAIL] the run $label did not complete (see runs/$label.log)"
    echo "MPIConnectionOrder: FAILED"
    exit 1
  fi
  mv OUTPUT "runs/$label"
}

# same <label-a> <label-b>: the two block files byte-identical
same() {
  local f n=0 bad=0
  for f in runs/"$1"/vtk/*.vts; do
    n=$((n + 1))
    cmp -s "$f" "runs/$2/${f#runs/$1/}" || { echo "  [FAIL] ${f#runs/$1/} differs ($1 vs $2)"; bad=1; }
  done
  [[ $n -eq 2 ]] || { echo "  [FAIL] $n block files, 2 expected"; bad=1; }
  [[ $bad -eq 0 ]] && echo "  [ok  ] $1 = $2 bit for bit"
  return $bad
}

run 1 sorted
run 2 sorted
run 2 reversed
echo "-- two blocks with permuted indices, one block per rank --"
fail=0
same sorted-np1 sorted-np2 || fail=1
same sorted-np1 reversed-np2 || fail=1
if [[ $fail -ne 0 ]]; then
  echo "MPIConnectionOrder: FAILED"
  exit 1
fi
echo "MPIConnectionOrder: all checks passed"
exit 0
