#!/usr/bin/env bash
#===============================================================================
#  Restart counters (fast tier, I3/F6)
#
#  A restart must carry on the counters that decide what the run does and
#  writes next. The case below solves the Sod tube of test/1D/Sod79 (100
#  cells) in two parts, with a restart in between (newrun = false):
#
#   transient  a time-accurate restart from a steady vtk raw solution starts
#              the transient at t = 0: after 5 steps the solution time of the
#              last output is the sum of their time steps, in (0, 0.1).
#
#  Reuses the Sod inputs (test/1D/Sod79); nothing else is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79
COMMON=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

fail() { echo "RestartCounters: FAIL — $*"; exit 1; }

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt

# prepare <run directory>
prepare() {
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt; do
    ln -s "../../../../INPUT/$f" "$1/INPUT/$f"
  done
}

# solve <run directory> <label>: the solver must exit 0 and report no error
solve() {
  ( cd "$1" && "../../../$LOCAL" > "solve_$2.log" 2>&1 ) || { tail -5 "$1/solve_$2.log"; fail "solver exited non-zero ($1, $2)"; }
  if command grep -q -E 'ERROR|NaN or p<0' "$1/solve_$2.log"; then
    command grep -E 'ERROR|NaN or p<0' "$1/solve_$2.log" | head -3
    fail "the solver stopped with an error ($1, $2)"
  fi
}

# write_ini <run directory> <case> <iter-threshold> <newrun>
write_ini() {
  local dir=$1 case=$2 iters=$3 newrun=$4
  {
    printf '%s\n' '[MOSE-Parameters]' "iter-threshold = $iters" 'res-threshold = 1e-300' "newrun = $newrun" ''
    printf '%s\n' '[MOSE-Physics]' 'equations = euler' ''
    printf '%s\n' '[MOSE-Numerics]' 'time-scheme = RK3' 'integration-variables = cons' \
                  'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC'
    case $case in
      steady)    printf '%s\n' 'cfl = 0.5' 'cfl-rise-threshold = 30' 'time-accurate = false' ;;
      dtime)     printf '%s\n' 'cfl = 0.9' 'time-accurate = true' ;;
      multigrid|multigrid-tecplot) printf '%s\n' 'cfl = 0.5' 'time-accurate = false' ;;
      transient) if [[ $newrun == true ]]; then printf '%s\n' 'cfl = 0.5' 'time-accurate = false'
                 else printf '%s\n' 'cfl = 0.9' 'time-accurate = true'; fi ;;
    esac
    local fmt='vtk raw'
    [[ $case == multigrid-tecplot ]] && fmt='tecplot ascii'
    printf '\n%s\n' '[MOSE-IO]' "sol-format = $fmt" 'sol-overwrite = false' 'shell-diter = 1000000000' \
                     'ini-diter = 1000000000'
    case $case in
      steady)    printf '%s\n' 'sol-diter = 7' '' '[MOSE-Probes]' 'probe1 = p1' '' '[p1]' \
                               'variables = rho(1) u p' 'diter = 3' 'index-position = 1 50 1 1' ;;
      dtime)     printf '%s\n' 'sol-dtime = 0.0083' ;;
      multigrid) printf '%s\n' 'sol-diter = 4' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' \
                               "level1-iter = $iters" ;;
      multigrid-tecplot) printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' "level1-iter = $iters" ;;
    esac
  } > "$dir/input.ini"
}

mkdir -p runs
status=0
# transient: 10 steady iterations, then a time-accurate restart for 5 steps
prepare runs/transient/split
write_ini runs/transient/split transient 10 true;  solve runs/transient/split steady
write_ini runs/transient/split transient 5 false;  solve runs/transient/split transient
python3 -B - runs/transient/split/OUTPUT/field.vtm <<'EOF' || status=1
import sys
from vtk_time import time_of
t = time_of(sys.argv[1])
ok = 0.0 < t < 0.1
print('  transient: solution time after 5 steps from the steady solution = %.17g (%s)' % (t, 'ok' if ok else 'should be in (0, 0.1)'))
sys.exit(0 if ok else 1)
EOF

if [[ $status -ne 0 ]]; then
  fail "a restart does not carry on the counters"
fi
echo "RestartCounters: ok — a transient from a steady vtk raw solution starts at t = 0"
exit 0
