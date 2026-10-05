#!/usr/bin/env bash
#===============================================================================
#  Restart counters (fast tier, I3/F6)
#
#  An exact restart from vtk raw must also carry on the counters that decide
#  what the run does and writes next. Each case below solves the Sod tube of
#  test/1D/Sod79 (100 cells) twice, in one go ("full") and in two parts with a
#  restart in between ("split", newrun = false), and the two OUTPUT/ trees
#  must be identical byte for byte (../check_restart_vtk.py): numbered
#  solutions, diagnostic files, residual-history.dat, probe files.
#
#   steady     local time stepping, CFL ramp over 30 iterations, solutions
#              every 7 iterations (sol-overwrite = false) with the diagnostic
#              output, a probe every 3 iterations; 20 + 20 against 40. Carries
#              the iteration count, the global iteration count (residual
#              history, probe rows, steady solution time) and the output
#              number.
#   dtime      time-accurate, a solution every 0.0083 s of physical time,
#              25 + 25 against 50 iterations with the first part ending
#              between two outputs. Carries the time of the last output.
#   multigrid  two grid levels, 10 iterations on the coarse one, then 30 on
#              the fine one, the restart after 15 of them; solutions every 4
#              iterations of the level. Carries the level.
#  and
#   transient  a time-accurate restart from a steady solution starts the
#              transient at t = 0: after 5 steps the solution time of the
#              last output is the sum of their time steps, in (0, 0.1).
#  and
#   probes     time-accurate, a solution every 10 iterations and three probes:
#              p1 every 0.01 s of physical time, p2 every 4 iterations, p3
#              every 7 iterations and every 0.03 s. 25 iterations, then the
#              final solution of the first part (field.vtm) is removed, as if
#              that run had stopped before writing it, and the restart reads
#              field2 (iteration 20) and runs 30 iterations, against 50 in one
#              go: residual-history.dat and the probe files of the first part
#              also hold rows 21 to 25, which the restart must not keep twice.
#              In the run in one go, p1 writes one row per multiple of its
#              dtime, from the first, and p2 only rows of iterations.
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
CHECK=../check_restart_vtk.py

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
# Boundary conditions of the coarse grid (50 cells) for the multigrid case:
# the same faces, with the cell indices of the coarse grid
awk '{ if ($5 == 2) $2 = 50; if ($5 >= 3 && $2 > 50) next; print }' "$SOD/INPUT/bc_x1.txt" > INPUT/bc2.txt

# prepare <run directory>
prepare() {
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt bc2.txt; do
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
      dtime|probes) printf '%s\n' 'cfl = 0.9' 'time-accurate = true' ;;
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
      probes)    printf '%s\n' 'sol-diter = 10' '' '[MOSE-Probes]' 'probe1 = p1' 'probe2 = p2' 'probe3 = p3' '' \
                               '[p1]' 'variables = rho(1) u p' 'dtime = 0.01' 'index-position = 1 50 1 1' '' \
                               '[p2]' 'variables = rho(1) u p' 'diter = 4' 'index-position = 1 30 1 1' '' \
                               '[p3]' 'variables = rho(1) u p' 'diter = 7' 'dtime = 0.03' 'index-position = 1 70 1 1' ;;
      multigrid) printf '%s\n' 'sol-diter = 4' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' \
                               "level1-iter = $iters" ;;
      multigrid-tecplot) printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' "level1-iter = $iters" ;;
    esac
  } > "$dir/input.ini"
}

# run_case <case> <total> <first part>: full run against split run + restart
run_case() {
  local case=$1 total=$2 first=$3
  prepare "runs/$case/full"
  write_ini "runs/$case/full" "$case" "$total" true;          solve "runs/$case/full" full
  prepare "runs/$case/split"
  write_ini "runs/$case/split" "$case" "$first" true;         solve "runs/$case/split" first
  write_ini "runs/$case/split" "$case" "$((total - first))" false; solve "runs/$case/split" restart
  python3 -B "$CHECK" "runs/$case/full/OUTPUT" "runs/$case/split/OUTPUT" > "runs/$case/compare.log" 2>&1
  local rc=$?
  sed "s/^Restart (vtk raw)/  $case/" "runs/$case/compare.log"
  return $rc
}

mkdir -p runs
status=0
run_case steady 40 20 || status=1

# dtime: the premise is that the first part ends strictly between two outputs of
# the full run, with outputs before and after it (checked on the full run)
run_case dtime 50 25 || status=1
python3 -B - runs/dtime/full/OUTPUT 25 <<'EOF' || status=1
import sys
from pathlib import Path
from vtk_time import time_of
out, split = Path(sys.argv[1]), int(sys.argv[2])
times = sorted(time_of(p) for p in out.glob('field[0-9]*.vtm'))
rows = [l.split() for l in (out / 'residual-history.dat').read_text().splitlines() if l.strip()]
t_split = float(rows[split - 1][1])
ok = any(t < t_split for t in times) and any(t > t_split for t in times) and all(t != t_split for t in times)
print('  dtime: %d time-based outputs in the full run, the first part ends at t = %.6g %s'
      % (len(times), t_split, 'between two of them' if ok else '(premise not met)'))
sys.exit(0 if ok else 1)
EOF

run_case multigrid 30 15 || status=1   # with levels, level1-iter counts the iterations on the fine level

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

# probes: 25 iterations, then the restart from the solution of iteration 20 (the final solution of the first
# part removed, as after a run that stopped before writing it); residual-history.dat and the probe files
# must be those of the run in one go
prepare runs/probes/full
write_ini runs/probes/full probes 50 true;    solve runs/probes/full full
prepare runs/probes/split
write_ini runs/probes/split probes 25 true;   solve runs/probes/split first
rm -f runs/probes/split/OUTPUT/field.vtm
write_ini runs/probes/split probes 30 false;  solve runs/probes/split restart
python3 -B "$CHECK" runs/probes/full/OUTPUT runs/probes/split/OUTPUT > runs/probes/compare.log 2>&1 || status=1
sed "s/^Restart (vtk raw)/  probes/" runs/probes/compare.log
# the time-based rows of the run in one go: p1 one row per multiple of dtime = 0.01 from the first (none at the first
# step, which ends at t < 0.01), p2 (diter only) no time row at all. A probe row is the iteration or the time and its
# three variables, written list-directed: some compilers break it over two lines
python3 -B - runs/probes/full/OUTPUT <<'EOF' || status=1
import sys
from pathlib import Path
out = Path(sys.argv[1])
def rows(name, n):
    v = (out / name).read_text().split()
    return [v[k:k + n] for k in range(0, len(v), n)]
is_time = lambda tok: any(c in tok for c in '.eEdD')
p1, p2 = rows('p1.txt', 4), rows('p2.txt', 4)
t_end = float([l.split() for l in (out / 'residual-history.dat').read_text().splitlines() if l.strip()][-1][1])
t1 = [float(r[0]) for r in p1 if is_time(r[0])]
ok1 = len(t1) == len(p1) == int(t_end / 0.01) and t1[0] >= 0.01
ok2 = len(p2) == 12 and not any(is_time(r[0]) for r in p2)
print('  probes: p1 %d time rows from t = %.6g (%s), p2 %d rows, %d of them by time (%s)'
      % (len(t1), t1[0] if t1 else -1, 'ok' if ok1 else 'one per multiple of 0.01 from the first expected',
         len(p2), sum(is_time(r[0]) for r in p2), 'ok' if ok2 else '12 iteration rows expected'))
sys.exit(0 if ok1 and ok2 else 1)
EOF

if [[ $status -ne 0 ]]; then
  fail "a restart does not carry on the counters"
fi
echo "RestartCounters: ok — steady, dtime, multigrid and probes restarts reproduce the uninterrupted runs bit for bit, and a transient from a steady solution starts at t = 0"
exit 0
