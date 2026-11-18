#!/usr/bin/env bash
#===============================================================================
#  Multigrid cycle: the checks of the set-up and of input.ini read again, and
#  the grid level in the shell line
#
#  A run whose cycle (cycle = v, w) cannot work stops with the reason, at the
#  set-up and when input.ini is read again during the run; with
#  ini-error = warn, the default, a reading that breaks a check is not taken.
#  1. set-up, on the Sod tube of test/1D/Sod79 (100 cells, steady):
#     a. cycle = w with levels = 1: no coarser level for the cycle;
#     b. cycle = w on two levels with pre-iter = 0 and post-iter = 0: the
#        running level would not iterate;
#     c. cycle = w and cycle-turbulence = true with SSG/LRR (the TMR flat
#        plate of test/2D/viscous/flat-plate-turbulent-tmr, 35x25, on two
#        levels): the cycle solves the RANS equations of one- and
#        two-equation models only;
#     each run must exit non-zero, at the set-up, with the reason.
#  2. input.ini read again, on the Spalart-Allmaras box of the MPI tests
#     (../../../mpi/bit-identity/make_case.py SA --levels 2, 2 x 2 blocks of
#     16 x 16 x 4 cells) on two levels, no iteration on the coarse one, 240 on
#     the finest, input.ini read again at its iteration 200 and changed once
#     the residual history has its row of iteration 10:
#     a. cycle = w written in a run started with cycle = none, ini-error =
#        stop: the run must stop at the reading, the cycle not being
#        allocated;
#     b. the same with ini-error left out (warn): the reading is not taken,
#        the run ends normally with the residual history of the same run
#        without the change, byte for byte;
#     c. pre-iter = 0 and post-iter = 0 written in a run with cycle = w,
#        ini-error = stop: the run must stop at the reading.
#     If the change comes after the reading (a host too slow to poll), the
#     test says so and fails.
#  3. the shell line, on the Sod tube with cycle = w on two levels, a line
#     every 5 iterations: with level2-iter = 0 the coarse level serves the
#     cycle only, and every line is that of a single grid ("MOSE | Iter");
#     with level2-iter = 10 the lines name the level ("MOSE Grid Level 2",
#     then "MOSE Grid Level 1").
#  Reuses the inputs of the cases named; nothing else is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79
TMR=../../../2D/viscous/flat-plate-turbulent-tmr
AIR=../../../common/Air
AIRV=../../../common/Air-viscous
SOLVER="$HERE/bin/MOSE"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

status=0
mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

# Sod tube: its inputs, and the boundary records of the coarse grid (50 cells)
mkdir -p INPUT/sod
ln -sf "../../$AIR/phase.txt"     INPUT/sod/phase.txt
ln -sf "../../$AIR/thermo.dat"    INPUT/sod/thermo.dat
ln -sf "../../$SOD/INPUT/ic_x1.tec" INPUT/sod/ic.tec
ln -sf "../../$SOD/INPUT/bc_x1.txt" INPUT/sod/bc.txt
awk '{ if ($5 == 2) $2 = 50; if ($5 >= 3 && $2 > 50) next; print }' "$SOD/INPUT/bc_x1.txt" > INPUT/sod/bc2.txt
# Spalart-Allmaras box on two levels
python3 -B ../../../mpi/bit-identity/make_case.py SA --levels 2 --outdir INPUT/sa > /dev/null \
  || { echo "MultigridChecks: FAIL — case generation (SA box)"; exit 1; }
for f in phase.txt thermo.dat transport.dat; do ln -sf "../../../$AIRV/$f" "INPUT/sa/INPUT/$f"; done
# TMR flat plate, SSG/LRR, 35x25
python3 -B "$TMR/make_case.py" --model SSGLRR --grid 35x25 --out INPUT/ssg > /dev/null \
  || { echo "MultigridChecks: FAIL — case generation (TMR plate)"; exit 1; }
for f in phase.txt thermo.dat transport.dat; do ln -sf "../../../$AIRV/$f" "INPUT/ssg/INPUT/$f"; done

# sod <name> <lines of [MOSE-Multigrid]>: input.ini of the Sod tube in runs/<name>, 30 iterations of the finest level
sod() {
  local d=runs/$1; shift
  rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt bc2.txt; do ln -s "../../../INPUT/sod/$f" "$d/INPUT/$f"; done
  printf '%s\n' '[MOSE-Parameters]' 'res-threshold = 1e-300' 'newrun = true' 'iter-threshold = 30' '' \
    '[MOSE-Physics]' 'equations = euler' '' \
    '[MOSE-Numerics]' 'cfl = 0.5' 'time-accurate = false' 'time-scheme = RK3' 'integration-variables = cons' \
    'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC' '' \
    '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 1000' 'shell-diter = 5' 'ini-diter = 1000000000' '' \
    '[MOSE-Multigrid]' "$@" > "$d/input.ini"
}

# refused <directory> <case> <reason>: the run exits non-zero at the set-up (no row of residual history) with the reason
refused() {
  local d=$1 tag=$2 why=$3 rc rows
  ( cd "$d" && exec "$SOLVER" > solve.log 2>&1 )
  rc=$?
  rows=$(awk 'NF' "$d/OUTPUT/residual-history.dat" 2>/dev/null | wc -l)
  if [[ $rc -ne 0 && $rows -eq 0 ]] && command grep -a -q -F "$why" "$d/solve.log"; then
    echo "  ok — $tag: stopped at the set-up: $why"
  else
    echo "  FAIL — $tag: exit code $rc, rows of residual history $rows (required: non-zero, 0), reason \"$why\"" \
         "(see $d/solve.log)"; status=1
  fi
}

echo "-- 1. checks of the set-up --"
sod setup-levels 'levels = 1' 'level1-iter = 30' 'cycle = w'
refused runs/setup-levels "a. cycle = w with levels = 1" \
  'The multigrid cycle (cycle = v, w) needs the coarser grid levels'
sod setup-smooth 'levels = 2' 'level2-iter = 10' 'level1-iter = 30' 'cycle = w' 'pre-iter = 0' 'post-iter = 0'
refused runs/setup-smooth "b. cycle = w with pre-iter = 0 and post-iter = 0" \
  'The multigrid cycle needs pre-iter + post-iter of at least 1'
d=runs/setup-ssglrr
rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
for f in phase.txt thermo.dat transport.dat ic.tec bc.txt; do ln -s "../../../INPUT/ssg/INPUT/$f" "$d/INPUT/$f"; done
sed -e 's/^iter-threshold = .*/iter-threshold = 10/' -e 's/^sol-diter .*/sol-diter   = 1000/' "$TMR/input-SSGLRR.ini" > "$d/input.ini"
printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' 'level1-iter = 10' 'cycle = w' \
  'cycle-turbulence = true' >> "$d/input.ini"
refused "$d" "c. cycle = w and cycle-turbulence = true with SSG/LRR" \
  'cycle-turbulence = true needs a RANS model of one or two equations'

# box <directory> <ini-error, or empty> <cycle>: input.ini of the SA box, written whole, then moved in place
box() {
  local d=$1 err=$2 cyc=$3 smooth=${4:-}
  {
    printf '%s\n' '[MOSE-Parameters]' 'res-threshold = 1e-300' 'newrun = true' '' \
      '[MOSE-Physics]' 'equations = navier-stokes' 'turbulence = SA' '' \
      '[MOSE-Numerics]' 'cfl = 0.5' 'vnn = 0.3' 'time-accurate = false' 'time-scheme = RK3' \
      'integration-variables = cons' 'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC' '' \
      '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 1000' 'shell-diter = 1000000000' 'ini-diter = 200'
    if [[ -n $err ]]; then echo "ini-error = $err"; fi
    printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 0' 'level1-iter = 240' "cycle = $cyc"
    if [[ -n $smooth ]]; then printf '%s\n' 'pre-iter = 0' 'post-iter = 0'; fi
  } > "$d/input.ini.new"
  mv "$d/input.ini.new" "$d/input.ini"
}

prepare_box() {  # <directory>
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat transport.dat ic.tec bc.txt bc2.txt; do ln -s "../../../INPUT/sa/INPUT/$f" "$1/INPUT/$f"; done
}

# reread <directory> <ini-error, or empty> <cycle at the start> <cycle written> [smooth]: the run of the box, input.ini
# written again once the residual history has its row of iteration 10; sets first, late, rc and last
reread() {
  local d=$1 err=$2 c0=$3 c1=$4 smooth=${5:-} pid n
  prepare_box "$d"; box "$d" "$err" "$c0"
  ( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 10 ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 10 ]] && first=yes
  box "$d" "$err" "$c1" "$smooth"
  # the change came too late if the run had already passed the reading of level-1 iteration 200
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 200 ]] && late=yes
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
}

# stopped <directory> <case> <reason>: the run stopped at the reading, before its last iteration, with the reason
stopped() {
  local d=$1 tag=$2 why=$3
  if [[ $first != yes || $late == yes ]]; then
    echo "  FAIL — $tag: the change did not come between iterations 10 and 200 (host too slow to poll)"; status=1
  elif [[ $rc -ne 0 && ${last:-0} -lt 240 ]] && command grep -a -q -F "$why" "$d/solve.log" \
       && command grep -a -q -E '(ABORT|STOP): ' "$d/solve.log"; then
    echo "  ok — $tag: the run stopped at the reading: $why"
  else
    echo "  FAIL — $tag: exit code $rc, last iteration ${last:-none}, reason \"$why\" (see $d/solve.log)"; status=1
  fi
}

echo "-- 2. checks of input.ini read again --"
d=runs/box-ref
prepare_box "$d"; box "$d" "" none
( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) || { echo "  FAIL — the run without a change failed (see $d/solve.log)"; status=1; }
reread runs/reread-stop stop none w
stopped runs/reread-stop "a. cycle = w in a run started without it, ini-error = stop" \
  'The run started without the multigrid cycle, which cannot be switched on during the run'
reread runs/reread-warn "" none w
nw=$(command grep -a -c -F 'is not taken' runs/reread-warn/solve.log)
nr=$(command grep -a -c -F 'The run started without the multigrid cycle' runs/reread-warn/solve.log)
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — b. ini-error = warn: the change did not come between iterations 10 and 200 (host too slow to poll)"; status=1
elif [[ $rc -eq 0 && $nw -eq 1 && $nr -eq 1 ]] && cmp -s runs/box-ref/OUTPUT/residual-history.dat runs/reread-warn/OUTPUT/residual-history.dat; then
  echo "  ok — b. cycle = w in a run started without it, ini-error = warn: the reading is not taken, the run goes on as without the change"
else
  echo "  FAIL — b. ini-error = warn: exit code $rc, warnings $nw, reasons $nr (required 0, 1, 1), residual history" \
       "of the run without the change (see runs/reread-warn/solve.log)"; status=1
fi
reread runs/reread-smooth stop w w smooth
stopped runs/reread-smooth "c. pre-iter = 0 and post-iter = 0 in a run with cycle = w, ini-error = stop" \
  'The multigrid cycle needs pre-iter + post-iter of at least 1'

echo "-- 3. the grid level in the shell line --"
sod shell-cycle 'levels = 2' 'level2-iter = 0' 'level1-iter = 30' 'cycle = w'
( cd runs/shell-cycle && exec "$SOLVER" > solve.log 2>&1 ) || { echo "  FAIL — shell-cycle failed"; status=1; }
nl=$(command grep -a -c 'Grid Level' runs/shell-cycle/solve.log)
np=$(command grep -a -c '^MOSE | Iter' runs/shell-cycle/solve.log)
if [[ $nl -eq 0 && $np -gt 0 ]]; then
  echo "  ok — coarse level for the cycle only: $np lines of a single grid, none with the grid level"
else
  echo "  FAIL — coarse level for the cycle only: $nl lines with the grid level, $np of a single grid (required 0, > 0)"; status=1
fi
sod shell-sequence 'levels = 2' 'level2-iter = 10' 'level1-iter = 30' 'cycle = w'
( cd runs/shell-sequence && exec "$SOLVER" > solve.log 2>&1 ) || { echo "  FAIL — shell-sequence failed"; status=1; }
n2=$(command grep -a -c '^MOSE Grid Level 2 | Iter' runs/shell-sequence/solve.log)
n1=$(command grep -a -c '^MOSE Grid Level 1 | Iter' runs/shell-sequence/solve.log)
if [[ $n2 -gt 0 && $n1 -gt 0 ]]; then
  echo "  ok — grid sequence: $n2 lines of level 2, $n1 of level 1"
else
  echo "  FAIL — grid sequence: $n2 lines of level 2, $n1 of level 1 (required > 0 each)"; status=1
fi

if [[ $status -eq 0 ]]; then
  echo "MultigridChecks: ok"
else
  echo "MultigridChecks: FAIL"
fi
exit $status
