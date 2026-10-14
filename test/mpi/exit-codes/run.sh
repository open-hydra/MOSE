#!/usr/bin/env bash
#===============================================================================
#  Exit codes of the error stops under MPI (needs-mpi tier)
#
#  An error that rank 0 alone finds at set-up must stop every rank: the
#  100-cell Sod tube of test/1D/Sod79 with cfl = -1, on 2 ranks, must end
#  within 60 s with exit code 1 to 128 and the message of the rule; and so
#  must the run without input.ini, which every rank finds it cannot read
#  (the ranks stop together, with MPI_FINALIZE). A rank
#  that ends alone, with STOP and without MPI_FINALIZE, leaves the others to
#  the launcher: the mpiexec of Intel MPI 2021.12 kills them and ends with
#  exit code 255. The inputs come from ../../fast/input-checks/exit-codes/make_case.py
#  (the Sod ones); nothing else is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79/INPUT
COMMON=../../common/Air
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/*
ln -sf "../$COMMON/phase.txt" INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat
python3 -B ../../fast/input-checks/exit-codes/make_case.py sod "$SOD" || exit 1
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 50
newrun = true

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = -1
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOI

fail=0
# run2 <tag> <message>: 2 ranks must end within 60 s with exit code 1 to 128 and the message
run2() {
  local tag="$1" message="$2"
  timeout -k 10 60 "$MPIEXEC" "$NPFLAG" 2 "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  if [[ $rc -ge 1 && $rc -le 128 && $rc -ne 124 ]] && grep -a -q -- "$message" "OUTPUT/$tag.log"; then
    echo "  [ok  ] $tag, 2 ranks, exit code $rc, \"$message\""
  else
    echo "  [FAIL] $tag, 2 ranks, exit code $rc (required 1 to 128 within 60 s; 124 is the time limit), message \"$message\" $(grep -a -q -- "$message" "OUTPUT/$tag.log" && echo found || echo missing) (see OUTPUT/$tag.log)"
    fail=1
  fi
}
echo "-- an error that rank 0 finds at set-up stops every rank (2 ranks) --"
run2 mpi-input "cfl must be > 0"
echo "-- an input.ini that cannot be read stops every rank (2 ranks) --"
rm -f input.ini
run2 mpi-ini-missing "input.ini cannot be opened"
if [[ $fail -eq 0 ]]; then
  echo "MPIExitCodes: all checks passed"
  exit 0
fi
echo "MPIExitCodes: FAILED"
exit 1
