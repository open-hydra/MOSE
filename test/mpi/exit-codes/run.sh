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
#  code 141, SIGPIPE, and lines lost). The run without input.ini, which every
#  rank finds it cannot read, must end within 60 s with exit code 1 to 128
#  and its message (the ranks stop together, with MPI_FINALIZE). A rank
#  that ends alone, with STOP and without MPI_FINALIZE, leaves the others to
#  the launcher: the mpiexec of Intel MPI 2021.12 kills them and ends with
#  exit code 255. The inputs come from ../../fast/input-checks/exit-codes/make_case.py
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
#  A value out of its rule that one rank alone reads again must stop every
#  rank as well: two ranks launched as two programs (MPMD), each in its own
#  folder with its own input.ini (mpmd/runs/A for rank 0, mpmd/runs/B for
#  rank 1), on the 2x2-block viscous box of RuntimeInput (../make_case.py,
#  48 x 48 cells per block), input.ini read again every 20 iterations. Once
#  the solution of iteration 10 is written, the input.ini of rank 1 alone
#  gets cfl = -1: at iteration 20 rank 1 finds the value out of its rule
#  and rank 0 does not. The run must end within 60 s with exit code 1, at
#  iteration 20, with the rule in the log. A rank that stopped alone would
#  wait for the others in the barrier of its stop, and the run would hang
#  until the time limit (exit code 124). The change must come before the
#  run writes its solution of iteration 20; if it does not (a host too slow
#  to poll), the test says so and fails. The MOSE processes of this run that
#  the launcher leaves behind (this copy of the solver, working in
#  mpmd/runs/A or mpmd/runs/B) are named, stopped, and fail the test.
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

echo "-- a value out of its rule that one rank alone reads again stops every rank (2 ranks, two folders) --"
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
  echo "  [FAIL] mpi-one-rank, case generation"
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
    > OUTPUT/mpi-one-rank.log 2>&1 &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $A/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $A/OUTPUT/field1.vtm ]] && first=yes
  [[ -f $A/OUTPUT/field2.vtm ]] && late=yes
  write_box_ini "$B" -1
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
  grep -a -q -- "cfl must be > 0" OUTPUT/mpi-one-rank.log && rule=found
  if [[ $first != yes ]]; then
    echo "  [FAIL] mpi-one-rank, the run ended before its solution of iteration 10, exit code $rc (see OUTPUT/mpi-one-rank.log)"
    fail=1
  elif [[ $late == yes ]]; then
    echo "  [FAIL] mpi-one-rank, input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"
    fail=1
  elif [[ $rc -eq 1 && ${last:-} == 20 && $rule == found ]] && grep -a -q -- "read again" OUTPUT/mpi-one-rank.log; then
    echo "  [ok  ] mpi-one-rank, 2 ranks, cfl = -1 read again by rank 1 alone: exit code 1 at iteration 20, \"cfl must be > 0\""
  else
    echo "  [FAIL] mpi-one-rank, 2 ranks, cfl = -1 read again by rank 1 alone: exit code $rc (required 1 within 60 s; 124 is" \
         "the time limit), last row of the residual history ${last:-none} (required 20), rule \"cfl must be > 0\" $rule" \
         "(see OUTPUT/mpi-one-rank.log)"
    fail=1
  fi
  if [[ -n $left ]]; then
    echo "  [FAIL] mpi-one-rank, MOSE processes of this run left behind by the launcher, now stopped: $left"
    fail=1
  fi
fi
if [[ $fail -eq 0 ]]; then
  echo "MPIExitCodes: all checks passed"
  exit 0
fi
echo "MPIExitCodes: FAILED"
exit 1
