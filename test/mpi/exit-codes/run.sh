#!/usr/bin/env bash
#===============================================================================
#  Exit codes of the error stops under MPI (needs-mpi tier)
#
#  An error that rank 0 alone finds at set-up must stop every rank: the
#  100-cell Sod tube of test/1D/Sod79 with cfl = -1, on 2 ranks, must end
#  within 60 s with exit code 1, the message of the rule and the STOP line
#  of both ranks: rank 0 broadcasts the error and the ranks stop together,
#  with MPI_FINALIZE, once rank 0 has written its message (with MPI_ABORT
#  from rank 0 alone, Intel MPI 2021.12 ended once in ten runs with exit
#  code 141, SIGPIPE, and lines lost). The run without input.ini, which rank 0
#  cannot read and whose outcome rank 0 passes to every rank, must end within
#  60 s with exit code 1 to 128 and its message (the ranks stop together, with
#  MPI_FINALIZE). A rank that ends alone, with STOP and without MPI_FINALIZE,
#  leaves the others to the launcher: the mpiexec of Intel MPI 2021.12 kills
#  them and ends with exit code 255. The inputs come from ../../fast/input-checks/exit-codes/make_case.py
#  (the Sod ones); nothing else is stored here.
#
#  A value of the initial condition that is not finite must stop every rank,
#  also when one rank alone holds it: the Sod tube with a NaN velocity in one
#  cell (make_case.py nan), its one block on rank 1 (INPUT/decomposition.map,
#  ../write_map.py), so that rank 1 alone holds the state of its cells and
#  rank 0 only its nodes, and the tube with an infinite coordinate of one
#  node (make_case.py inf), which both ranks hold, must end within 60 s with
#  exit code 1, the place of the value, written once, and the STOP line of
#  both ranks.
#
#  A value out of its rule written during the run must stop every rank when
#  input.ini is read again: two ranks launched as two programs (MPMD), each
#  in its own folder with its own input.ini (mpmd/runs/A for rank 0,
#  mpmd/runs/B for rank 1), on the 2x2-block viscous box of RuntimeInput
#  (../make_case.py, 48 x 48 cells per block), input.ini read again every 20
#  iterations. Once the solution of iteration 10 is written, the input.ini
#  of rank 0 gets cfl = -1: rank 0 reads input.ini and passes its text to
#  every rank, so at iteration 20 every rank finds the value out of its
#  rule (the input.ini of rank 1, unchanged, is not read: MPIRuntimeInput
#  checks that). The run must end within 60 s with exit code 1, at
#  iteration 20, with the rule in the log. A rank that stopped alone would
#  wait for the others in the barrier of its stop, and the run would hang
#  until the time limit (exit code 124). The change must come before the
#  run writes its solution of iteration 20; if it does not (a host too slow
#  to poll), the test says so and fails. The MOSE processes of this run that
#  the launcher leaves behind (this copy of the solver, working in
#  mpmd/runs/A or mpmd/runs/B) are named, stopped, and fail the test.
#
#  An initial condition that one rank alone cannot read must stop every
#  rank: the Sod tube in two folders (MPMD, mpmd/ic/A for rank 0, mpmd/ic/B
#  for rank 1), the initial condition of rank 1 cut in the middle of its
#  values, must end within 60 s with exit code 1, the line "Initial
#  conditions FAIL" of the set-up summary, the error of the read written
#  once with the rank that found it, and the STOP line of both ranks: the
#  ranks decide together on the error of the read. The processes left
#  behind are handled as above.
#
#  A rank that stops every rank with MPI_ABORT (an error that it finds
#  alone, here a probe of an unknown variable) writes out what it holds for
#  its output and its error output first: the Sod tube in air with
#  chemistry = finite-rate, whose set-up writes a warning of FLINT on the
#  error output, and a probe with the variables "p foo", run as one process
#  without the launcher with its output on a file (as CTest runs the fast
#  tests in an MPI build), must end within 60 s with exit code 1 to 128 and
#  keep in its log the message, the line "Boundary conditions OK" of the
#  set-up summary, written before it, and the warning of FLINT. Built with
#  gfortran and run with Open MPI 4.1.4, the process ended without writing
#  them out and the log held only the lines of Open MPI; with Intel MPI
#  2021.12 nothing was lost.
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
# run2 <tag> <message> [<stop>]: 2 ranks must end within 60 s with exit code 1 to 128 and the message; with <stop>, with
# exit code 1 and the line "STOP: <stop>" of both ranks (the ranks stop together)
run2() {
  local tag="$1" message="$2" stop="${3:-}"
  timeout -k 10 60 "$MPIEXEC" "$NPFLAG" 2 "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  local ok=1 nstop=0 need="1 to 128 within 60 s; 124 is the time limit"
  [[ $rc -ge 1 && $rc -le 128 && $rc -ne 124 ]] || ok=0
  grep -a -q -- "$message" "OUTPUT/$tag.log" || ok=0
  if [[ -n $stop ]]; then
    need="1 within 60 s, and the STOP line of both ranks"
    nstop=$(grep -a -c -F -- "] STOP: $stop" "OUTPUT/$tag.log")
    [[ $rc -eq 1 && $nstop -eq 2 ]] || ok=0
  fi
  if [[ $ok -eq 1 ]]; then
    echo "  [ok  ] $tag, 2 ranks, exit code $rc, \"$message\"${stop:+, the STOP line of both ranks}"
  else
    echo "  [FAIL] $tag, 2 ranks, exit code $rc (required $need), message \"$message\" $(grep -a -q -- "$message" "OUTPUT/$tag.log" && echo found || echo missing)${stop:+, STOP lines $nstop of 2} (see OUTPUT/$tag.log)"
    fail=1
  fi
}
echo "-- an error that rank 0 finds at set-up stops every rank (2 ranks) --"
run2 mpi-input "cfl must be > 0" "input file error (see above)"
echo "-- an input.ini that cannot be read stops every rank (2 ranks) --"
rm -f input.ini
run2 mpi-ini-missing "input.ini cannot be opened"

echo "-- a value of the initial condition that is not finite stops every rank, also when one rank alone holds it (2 ranks) --"
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 50
newrun = true

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
EOI
python3 -B ../write_map.py --layout 1x1 --cells 100 1 --owners list:1 > /dev/null || exit 1
for tag in mpi-nan-ic mpi-inf-node; do
  if [[ $tag == mpi-nan-ic ]]; then
    python3 -B ../../fast/input-checks/exit-codes/make_case.py nan "$SOD" || exit 1
    place="variable 2 (u) of block 1, cell i, j, k = 50, 1, 1"
  else
    python3 -B ../../fast/input-checks/exit-codes/make_case.py inf "$SOD" || exit 1
    place="coordinate y of block 1, node i, j, k = 20, 1, 1"
  fi
  run2 $tag "holds a value that is not finite (NaN or infinity): $place" "initial condition error (see above)"
  n=$(grep -a -c -F -- "holds a value that is not finite" "OUTPUT/$tag.log")
  if [[ $n -eq 1 ]]; then
    echo "  [ok  ] $tag, the place of the value written once"
  else
    echo "  [FAIL] $tag, the place of the value written $n times (required once, by the first rank that found it)"
    fail=1
  fi
done
rm -f INPUT/decomposition.map

echo "-- a value out of its rule written in the input.ini of rank 0 stops every rank when it is read again (2 ranks, two folders) --"
SOLVER="$(cd bin && pwd -P)/MOSE"
M=mpmd
write_box_ini() {  # <folder> <cfl>: the input.ini of the cfl case of RuntimeInput, written whole, then moved in place
  cat > "$1/input.ini.new" <<EOF
[MOSE-Parameters]
iter-threshold = 40
res-threshold  = 1e-300
newrun         = true

[MOSE-Physics]
equations = navier-stokes
schmidt   = 1.0

[MOSE-Numerics]
cfl                   = $2
time-accurate         = false
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format    = vtk raw
sol-diter     = 10
sol-overwrite = false
ini-diter     = 20
shell-diter   = 1000000000
EOF
  mv "$1/input.ini.new" "$1/input.ini"
}
left_behind() {  # the processes of this copy of the solver whose working folder is one of the two of the run
  local p
  for p in /proc/[0-9]*; do
    [[ $(readlink "$p/exe" 2>/dev/null) == "$SOLVER" ]] || continue
    case $(readlink "$p/cwd" 2>/dev/null) in
      "$RUNS/A" | "$RUNS/B") echo "${p#/proc/}" ;;
    esac
  done
}
rm -rf "$M"
mkdir -p "$M/INPUT" "$M/runs"
for f in phase.txt thermo.dat transport.dat; do
  ln -s "../../../../common/Air-viscous/$f" "$M/INPUT/$f"
done
if ! python3 -B ../make_case.py --layout 2x2 --cells 48 48 --bc box --zfaces 300 --outdir "$M" > /dev/null; then
  echo "  [FAIL] mpi-read-again, case generation"
  fail=1
else
  RUNS="$(cd "$M/runs" && pwd -P)"
  for d in A B; do  # two levels under mpmd/, where the links of their INPUT lead
    mkdir -p "$M/runs/$d/INPUT" "$M/runs/$d/OUTPUT"
    for f in phase.txt thermo.dat transport.dat ic.tec bc.txt; do
      ln -s "../../../INPUT/$f" "$M/runs/$d/INPUT/$f"
    done
    write_box_ini "$M/runs/$d" 0.5
  done
  A=$M/runs/A
  B=$M/runs/B
  timeout -k 10 60 "$MPIEXEC" "$NPFLAG" 1 -wdir "$RUNS/A" "$SOLVER" : "$NPFLAG" 1 -wdir "$RUNS/B" "$SOLVER" \
    > OUTPUT/mpi-read-again.log 2>&1 &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $A/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $A/OUTPUT/field1.vtm ]] && first=yes
  [[ -f $A/OUTPUT/field2.vtm ]] && late=yes
  write_box_ini "$A" -1
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$A/OUTPUT/residual-history.dat" 2>/dev/null)
  left=$(left_behind | paste -sd ' ' -)
  if [[ -n $left ]]; then
    kill -TERM $left 2>/dev/null
    sleep 2
    kill -KILL $left 2>/dev/null
  fi
  rule=missing
  grep -a -q -- "cfl must be > 0" OUTPUT/mpi-read-again.log && rule=found
  if [[ $first != yes ]]; then
    echo "  [FAIL] mpi-read-again, the run ended before its solution of iteration 10, exit code $rc (see OUTPUT/mpi-read-again.log)"
    fail=1
  elif [[ $late == yes ]]; then
    echo "  [FAIL] mpi-read-again, input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"
    fail=1
  elif [[ $rc -eq 1 && ${last:-} == 20 && $rule == found ]] && grep -a -q -- "read again" OUTPUT/mpi-read-again.log; then
    echo "  [ok  ] mpi-read-again, 2 ranks, cfl = -1 in the input.ini of rank 0, read again: exit code 1 at iteration 20, \"cfl must be > 0\""
  else
    echo "  [FAIL] mpi-read-again, 2 ranks, cfl = -1 in the input.ini of rank 0, read again: exit code $rc (required 1 within 60 s; 124 is" \
         "the time limit), last row of the residual history ${last:-none} (required 20), rule \"cfl must be > 0\" $rule" \
         "(see OUTPUT/mpi-read-again.log)"
    fail=1
  fi
  if [[ -n $left ]]; then
    echo "  [FAIL] mpi-read-again, MOSE processes of this run left behind by the launcher, now stopped: $left"
    fail=1
  fi
fi

echo "-- an initial condition that one rank alone cannot read stops every rank (2 ranks, two folders) --"
python3 -B ../../fast/input-checks/exit-codes/make_case.py sod "$SOD" || exit 1
mkdir -p "$M/ic"
for d in A B; do
  rm -rf "$M/ic/$d"
  mkdir -p "$M/ic/$d/INPUT" "$M/ic/$d/OUTPUT"
  cp INPUT/phase.txt INPUT/thermo.dat INPUT/bc.txt INPUT/ic.tec "$M/ic/$d/INPUT/"
  cat > "$M/ic/$d/input.ini" <<EOI
[MOSE-Parameters]
iter-threshold = 50
newrun = true

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
EOI
done
# the initial condition of rank 1 cut in the middle of its values: the reader meets the end of the file
head -c $(( $(wc -c < INPUT/ic.tec) / 2 )) INPUT/ic.tec > "$M/ic/B/INPUT/ic.tec"
RUNS="$(cd "$M/ic" && pwd -P)"
timeout -k 10 60 "$MPIEXEC" "$NPFLAG" 1 -wdir "$RUNS/A" "$SOLVER" : "$NPFLAG" 1 -wdir "$RUNS/B" "$SOLVER" \
  > OUTPUT/mpi-ic-read.log 2>&1
rc=$?
left=$(left_behind | paste -sd ' ' -)
if [[ -n $left ]]; then
  kill -TERM $left 2>/dev/null
  sleep 2
  kill -KILL $left 2>/dev/null
fi
nmsg=$(grep -a -c -F -- "[ERROR] reading input file INPUT/ic.tec (rank 1)" OUTPUT/mpi-ic-read.log)
nstop=$(grep -a -c -F -- "] STOP: initial condition error (see above)" OUTPUT/mpi-ic-read.log)
summary=missing
grep -a -q -E 'Initial conditions +FAIL' OUTPUT/mpi-ic-read.log && summary=found
if [[ $rc -eq 1 && $nmsg -eq 1 && $nstop -eq 2 && $summary == found ]]; then
  echo "  [ok  ] mpi-ic-read, 2 ranks, the initial condition of rank 1 cut: exit code 1, \"Initial conditions FAIL\"," \
       "the error of the read of rank 1 written once, the STOP line of both ranks"
else
  echo "  [FAIL] mpi-ic-read, 2 ranks, the initial condition of rank 1 cut: exit code $rc (required 1 within 60 s; 124 is" \
       "the time limit), \"Initial conditions FAIL\" $summary, the error of the read of rank 1 written $nmsg times (required" \
       "once), STOP lines $nstop of 2 (see OUTPUT/mpi-ic-read.log)"
  fail=1
fi
if [[ -n $left ]]; then
  echo "  [FAIL] mpi-ic-read, MOSE processes of this run left behind by the launcher, now stopped: $left"
  fail=1
fi

echo "-- a rank that stops every rank with MPI_ABORT writes out its output first (one process without the launcher) --"
python3 -B ../../fast/input-checks/exit-codes/make_case.py sod "$SOD" || exit 1
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2
newrun = true

[MOSE-Physics]
equations = euler
chemistry = finite-rate

[MOSE-Numerics]
cfl = 0.9
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-Chemistry]
ode-solver = sdirk4b

[MOSE-Probes]
probe1 = pA

[pA]
variables = p foo
diter = 1
position = 0.5 0.0 0.0
EOI
timeout -k 10 60 "$LOCAL" > OUTPUT/mpi-abort-log.log 2>&1
rc=$?
message=missing; summary=missing; warning=missing
grep -a -q -F -- 'unknown probe variable "foo"' OUTPUT/mpi-abort-log.log && message=found
grep -a -q -E 'Boundary conditions +OK' OUTPUT/mpi-abort-log.log && summary=found
grep -a -q -F -- 'FLINT Assign_Mechanism: explicit procedure for' OUTPUT/mpi-abort-log.log && warning=found
if [[ $rc -ge 1 && $rc -le 128 && $rc -ne 124 && $message == found && $summary == found && $warning == found ]]; then
  echo "  [ok  ] mpi-abort-log, one process without the launcher, stopped with MPI_ABORT: exit code $rc; in its log the" \
       "message, the line of the set-up summary written before it and the warning of FLINT on the error output"
else
  echo "  [FAIL] mpi-abort-log, one process without the launcher, stopped with MPI_ABORT: exit code $rc (required 1 to" \
       "128 within 60 s; 124 is the time limit), the message (unknown probe variable) $message, the line Boundary" \
       "conditions OK of the set-up summary $summary, the warning of FLINT on the error output $warning (see" \
       "OUTPUT/mpi-abort-log.log)"
  fail=1
fi
if [[ $fail -eq 0 ]]; then
  echo "MPIExitCodes: all checks passed"
  exit 0
fi
echo "MPIExitCodes: FAILED"
exit 1
