#!/usr/bin/env bash
#===============================================================================
#  Restart counters (fast tier, I3/F6)
#
#  An exact restart must also carry on the counters that decide what the run
#  does and writes next, in every format: the solution stores them, as vtk
#  field data or as Tecplot data set auxiliary data. Each case below solves the
#  Sod tube of test/1D/Sod79 (100 cells) twice, in one go ("full") and in two
#  parts with a restart in between ("split", newrun = false), and the two
#  OUTPUT/ trees must be identical byte for byte (../check_restart_vtk.py):
#  numbered solutions, diagnostic files, residual-history.dat, probe files.
#  The cases run with vtk raw; steady, dtime and multigrid also with vtk ascii,
#  Tecplot ASCII and, when MOSE is built with TecIO (MOSE_TECIO = 1, which the
#  CTest definitions set), Tecplot binary (.szplt): the text formats write 17
#  significant digits, which read back as the 64-bit values written.
#
#   steady     local time stepping, CFL ramp over 30 iterations, solutions
#              every 7 iterations (sol-overwrite = false) with the diagnostic
#              output, a probe every 3 iterations; 20 + 20 against 40. Carries
#              the iteration count, the global iteration count (residual
#              history, probe rows, steady solution time), the output
#              number and the residual of the first iteration (the field data
#              MOSE_RES_REFERENCE of the solutions; res-threshold is off here).
#   dtime      time-accurate, a solution every 0.0083 s of physical time,
#              25 + 25 against 50 iterations with the first part ending
#              between two outputs. Carries the time of the last output.
#   multigrid  two grid levels, 10 iterations on the coarse one, then 30 on
#              the fine one, the restart after 15 of them; solutions every 4
#              iterations of the level. Carries the level.
#  and, with Tecplot ASCII,
#   multigrid-tecplot  the same two levels: the global iteration count of the
#              restart goes on from the magnitude of the solution time read on
#              the finest level, minus the count in a steady solution (the
#              coarse levels hold none); the solution of the coarse level
#              (field-level2), written when the run moves to the fine one, has
#              the solution time of the prolongated solution written at the
#              same iteration, minus the count, and no STRANDID (static zones)
#  and
#   transient  a time-accurate restart from a steady solution starts the
#              transient at t = 0: after 5 steps the solution time of the
#              last output is the sum of their time steps, in (0, 0.1); from
#              vtk raw and from Tecplot ASCII (transient-tecplot), where the
#              steady solution is the one with a negative SOLUTIONTIME.
#   time-only  a steady restart from a vtk ascii solution of 4 iterations that
#              holds TIME only (CYCLE and the MOSE_* field data removed) goes
#              on from the iteration count |TIME|: the restart writes the row
#              of iteration 5 of residual-history.dat.
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
#   probes-cut the first part of probes, then the files cut by hand after the
#              rows of iteration 20, without a newline at the end:
#              residual-history.dat ends with its complete row of iteration
#              20, p2.txt with the first characters of its next row ("2" of
#              24). The restart drops the broken row, ends each file with a
#              newline and goes on as the run in one go of probes.
#   torque     time-accurate, a solution every 10 iterations, with the
#              rotating frame on (rigid-body, omega = 1 rad/s about x):
#              OUTPUT/torque.dat gets a header line and a row at every
#              iteration, as residual-history.dat. 25 iterations, field.vtm
#              removed, then the restart from field2 (iteration 20) for 30
#              iterations, against 50 in one go: torque.dat must keep its
#              header and rows 1 to 20, drop rows 21 to 25 of the first part
#              and go on as the run in one go. The Sod tube has no wall
#              faces, so the torque and power columns are 0: the rows differ
#              by their iteration and time.
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
    printf '%s\n' '[MOSE-Physics]' 'equations = euler'
    [[ $case == torque ]] && printf '%s\n' 'rotational-frame = rigid-body'
    printf '\n'
    printf '%s\n' '[MOSE-Numerics]' 'time-scheme = RK3' 'integration-variables = cons' \
                  'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC'
    case $case in
      steady)    printf '%s\n' 'cfl = 0.5' 'cfl-rise-threshold = 30' 'time-accurate = false' ;;
      dtime|probes|torque) printf '%s\n' 'cfl = 0.9' 'time-accurate = true' ;;
      multigrid|multigrid-tecplot|time-only) printf '%s\n' 'cfl = 0.5' 'time-accurate = false' ;;
      transient|transient-tecplot|transient-szplt) if [[ $newrun == true ]]; then printf '%s\n' 'cfl = 0.5' 'time-accurate = false'
                 else printf '%s\n' 'cfl = 0.9' 'time-accurate = true'; fi ;;
    esac
    local fmt=${FMT:-vtk raw}
    [[ $case == multigrid-tecplot || $case == transient-tecplot ]] && fmt='tecplot ascii'
    [[ $case == transient-szplt ]] && fmt='tecplot binary'
    [[ $case == time-only ]] && fmt='vtk ascii'
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
      torque)    printf '%s\n' 'sol-diter = 10' '' '[MOSE-Rotating-Frame]' 'omega = 1.0' 'axis = 1.0 0.0 0.0' \
                               'origin = 0.0 0.0 0.0' ;;
      multigrid) printf '%s\n' 'sol-diter = 4' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' \
                               "level1-iter = $iters" ;;
      multigrid-tecplot) printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' "level1-iter = $iters" ;;
    esac
  } > "$dir/input.ini"
}

# run_case <case> <total> <first part> [<sol-format>]: full run against split run + restart, with vtk raw or with
# the format given (in runs/<case>-<format>/)
run_case() {
  local case=$1 total=$2 first=$3 rc
  FMT=${4:-vtk raw}
  local d=runs/$case label=$case
  [[ $FMT != 'vtk raw' ]] && d=runs/$case-${FMT// /-} && label="$case ($FMT)"
  prepare "$d/full"
  write_ini "$d/full" "$case" "$total" true;          solve "$d/full" full
  prepare "$d/split"
  write_ini "$d/split" "$case" "$first" true;         solve "$d/split" first
  write_ini "$d/split" "$case" "$((total - first))" false; solve "$d/split" restart
  python3 -B "$CHECK" "$d/full/OUTPUT" "$d/split/OUTPUT" > "$d/compare.log" 2>&1
  rc=$?
  sed "s/^Restart (vtk raw)/  $label/" "$d/compare.log"
  FMT='vtk raw'
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

# The same three cases in the other formats
formats=('vtk ascii' 'tecplot ascii')
if [[ ${MOSE_TECIO:-0} == 1 ]]; then
  formats+=('tecplot binary')
else
  echo "  (Tecplot binary not run: MOSE_TECIO is not 1, MOSE built without TecIO)"
fi
for fmt in "${formats[@]}"; do
  run_case steady 40 20 "$fmt" || status=1
  run_case dtime 50 25 "$fmt" || status=1
  run_case multigrid 30 15 "$fmt" || status=1
done

# multigrid-tecplot: the same two levels with Tecplot ASCII and the final solution only. The solution of
# the coarse level, written when the run moves to the fine one, has the time of its iteration, and the
# iteration column of residual-history.dat keeps increasing across the restart (10 + 15 iterations, then
# the restart)
prepare runs/multigrid-tecplot/split
write_ini runs/multigrid-tecplot/split multigrid-tecplot 15 true;   solve runs/multigrid-tecplot/split first
python3 -B - runs/multigrid-tecplot/split/OUTPUT <<'EOF' || status=1
import re
import sys
from pathlib import Path
out = Path(sys.argv[1])


def stamp(name):
    """SOLUTIONTIME and STRANDID of the first zone of a Tecplot ASCII file (None when absent)."""
    zone = next(l for l in (out / name).read_text().splitlines() if l.lstrip().upper().startswith('ZONE'))
    t = re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)', zone)
    s = re.search(r'STRANDID\s*=\s*(-?\d+)', zone)
    return (float(t.group(1).replace('D', 'E')) if t else None, int(s.group(1)) if s else None)


coarse, fine = stamp('field-level2.tec'), stamp('field-prolongated.tec')
ok = coarse == fine and fine[0] is not None and fine[0] < 0 and fine[1] is None
print('  multigrid-tecplot: field-level2 SOLUTIONTIME, STRANDID %s, field-prolongated %s (%s)'
      % (coarse, fine, 'ok' if ok else 'the coarse level does not have the time of its iteration'))
sys.exit(0 if ok else 1)
EOF
write_ini runs/multigrid-tecplot/split multigrid-tecplot 15 false;  solve runs/multigrid-tecplot/split restart
python3 -B - runs/multigrid-tecplot/split/OUTPUT/residual-history.dat <<'EOF' || status=1
import sys
its = [int(l.split()[0]) for l in open(sys.argv[1]) if l.strip()]
ok = len(its) > 25 and all(b > a for a, b in zip(its, its[1:]))
print('  multigrid-tecplot: residual history rows %d to %d around the restart: %s (%s)'
      % (25, 26, its[24:26], 'ok' if ok else 'the iteration count starts again'))
sys.exit(0 if ok else 1)
EOF

# transient: 10 steady iterations, then a time-accurate restart for 5 steps, from vtk raw, from Tecplot ASCII and,
# with TecIO, from Tecplot binary (whose time is read from the last row of residual-history.dat)
cases=(transient transient-tecplot)
[[ ${MOSE_TECIO:-0} == 1 ]] && cases+=(transient-szplt)
for case in "${cases[@]}"; do
  prepare runs/$case/split
  write_ini runs/$case/split $case 10 true;  solve runs/$case/split steady
  write_ini runs/$case/split $case 5 false;  solve runs/$case/split transient
  sol=runs/$case/split/OUTPUT/field.vtm
  [[ $case == transient-tecplot ]] && sol=runs/$case/split/OUTPUT/field.tec
  [[ $case == transient-szplt ]] && sol=runs/$case/split/OUTPUT/residual-history.dat
  python3 -B - "$sol" "$case" <<'EOF' || status=1
import sys
from vtk_time import time_of
if sys.argv[1].endswith('.dat'):
    t = float([l.split() for l in open(sys.argv[1]) if l.strip()][-1][1])
else:
    t = time_of(sys.argv[1])
ok = 0.0 < t < 0.1
print('  %s: solution time after 5 steps from the steady solution = %.17g (%s)'
      % (sys.argv[2], t, 'ok' if ok else 'should be in (0, 0.1)'))
sys.exit(0 if ok else 1)
EOF
done

# time-only: 4 steady iterations in vtk ascii, the field data of the solution cut down to TIME (minus the iteration count),
# then a steady restart of one iteration, which must go on from |TIME| = 4: its row of residual-history.dat is iteration 5
prepare runs/time-only/split
write_ini runs/time-only/split time-only 4 true;  solve runs/time-only/split first
python3 -B - runs/time-only/split/OUTPUT/field.vtm <<'EOF' || status=1
import re
import sys
from pathlib import Path
vtm = Path(sys.argv[1])
for entry in re.findall(r'file="([^"]+)"', vtm.read_text()):
    f = vtm.parent / entry
    f.write_text(re.sub(r'\s*<DataArray[^>]*Name="(CYCLE|MOSE_[A-Z_]+)"[^>]*>.*?</DataArray>', '', f.read_text(), flags=re.S))
EOF
write_ini runs/time-only/split time-only 1 false;  solve runs/time-only/split restart
python3 -B - runs/time-only/split/OUTPUT/residual-history.dat <<'EOF' || status=1
import sys
its = [int(l.split()[0]) for l in open(sys.argv[1]) if l.strip()]
ok = its[-1] == 5
print('  time-only: the restart from a solution with TIME only writes the row of iteration %d (%s)'
      % (its[-1], 'ok' if ok else 'expected 5, |TIME| = 4 plus one'))
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
rm -rf runs/probes-cut; mkdir -p runs/probes-cut; cp -R runs/probes/split runs/probes-cut/split   # first part of probes-cut
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

# probes-cut: the first part of probes (25 iterations, field.vtm removed, copied above before the restart of probes),
# then residual-history.dat and p2.txt cut by hand after the rows of iteration 20, as a user may cut them, with no
# newline at the end: the history ends with its complete row of iteration 20, p2.txt (rows at iterations 4, 8, ...,
# 24) with the first characters of its row of iteration 24, up to the "2". The restart must give the files of the
# run in one go of probes.
python3 -B - runs/probes-cut/split/OUTPUT <<'EOF' || status=1
import re
import sys
from pathlib import Path
out = Path(sys.argv[1])
h = out / 'residual-history.dat'
rows = [l for l in h.read_text().splitlines() if l.strip()]
keep = [l for l in rows if int(l.split()[0]) <= 20]
h.write_text('\n'.join(keep))
p = out / 'p2.txt'
text = p.read_text()
values = list(re.finditer(r'\S+', text))
ok = len(keep) == 20 and len(values) == 24 and values[20].group() == '24'
if ok:
    p.write_text(text[:values[20].start() + 1])
print('  probes-cut: residual-history.dat cut after iteration 20 without a newline, p2.txt after "%s" (%s)'
      % (p.read_text()[-12:].replace('\n', '\\n'), 'ok' if ok else 'unexpected first part'))
sys.exit(0 if ok else 1)
EOF
write_ini runs/probes-cut/split probes 30 false;  solve runs/probes-cut/split restart
python3 -B "$CHECK" runs/probes/full/OUTPUT runs/probes-cut/split/OUTPUT > runs/probes-cut/compare.log 2>&1 || status=1
sed "s/^Restart (vtk raw)/  probes-cut/" runs/probes-cut/compare.log

# torque: the rotating frame writes OUTPUT/torque.dat, a header line and a row per residual row. 25 iterations, then
# the restart from the solution of iteration 20 (field.vtm removed, as in probes): torque.dat, compared with the other
# outputs, must be that of the run in one go, which must hold the header and the rows of iterations 1 to 50
prepare runs/torque/full
write_ini runs/torque/full torque 50 true;    solve runs/torque/full full
prepare runs/torque/split
write_ini runs/torque/split torque 25 true;   solve runs/torque/split first
rm -f runs/torque/split/OUTPUT/field.vtm
write_ini runs/torque/split torque 30 false;  solve runs/torque/split restart
python3 -B "$CHECK" runs/torque/full/OUTPUT runs/torque/split/OUTPUT > runs/torque/compare.log 2>&1 || status=1
sed "s/^Restart (vtk raw)/  torque/" runs/torque/compare.log
python3 -B - runs/torque/full/OUTPUT/torque.dat <<'EOF' || status=1
import sys
from pathlib import Path
p = Path(sys.argv[1])
lines = [l for l in p.read_text().splitlines() if l.strip()] if p.is_file() else []
heads = [k for k, l in enumerate(lines) if l.lstrip().startswith('#')]
its = [int(l.split()[0]) for l in lines if not l.lstrip().startswith('#')]
ok = heads == [0] and its == list(range(1, 51))
print('  torque: torque.dat of the run in one go, %d header line(s) and %d rows, iterations %s (%s)'
      % (len(heads), len(its), '%d to %d' % (its[0], its[-1]) if its else '-',
         'ok' if ok else 'one header line first and the rows of iterations 1 to 50 expected'))
sys.exit(0 if ok else 1)
EOF

if [[ $status -ne 0 ]]; then
  fail "a restart does not carry on the counters"
fi
echo "RestartCounters: ok — steady, dtime, multigrid, probes and rotating-frame (torque.dat) restarts reproduce the uninterrupted runs bit for bit, also from history files cut by hand without a final newline, a transient from a steady solution starts at t = 0 and a steady restart from a solution with TIME only goes on from |TIME|"
exit 0
