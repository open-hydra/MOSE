#!/usr/bin/env bash
#===============================================================================
#  Multirate: exact restart and exact output instants.
#  Every run writes vtk raw and its residual at every macro step (res-diter = 1). A "split" run makes the first part
#  with newrun = true, then restarts in the same directory (newrun = false, from OUTPUT/field.vtm) and makes the rest;
#  its OUTPUT tree must be identical byte by byte to that of the run in one go (test/checks/restart/check_restart_vtk.py:
#  numbered solutions, final solution, residual history), and neither may have an empty residual history.
#   R    mr1d-sod, mr2d-box-4b and the Spalart-Allmaras plate mr2d-plate-sa, multirate (max-level 3, safety 2) and the
#        one-level path (max-level 0, safety 1), levels recomputed at every macro step, solutions every 30 macro
#        steps: 100 = 50 + restart + 50.
#   R-tec mr1d-sod and mr2d-box-4b as in R, with sol-format = tecplot ascii: the restart reads the Tecplot solution
#        OUTPUT/field.tec (17 significant digits, the counters as data set auxiliary data); the OUTPUT trees (numbered
#        and final .tec solutions, residual history) identical byte by byte, and no run writes a vtk file.
#   R2   mr2d-shock, the same; mr1d-sod with recompute-every 4: 100 = 48 + restart + 52 (48 is a multiple of 4: the
#        restart rebuilds the map where the run in one go does); 50 + restart + 50 is reported, without a threshold.
#        recompute-every reaches the map: mr1d-sod with recompute-every 4 and nonlocal-cfl = off (no redo raises a map),
#        every macro step printed: the cells per level change, and only at the macro steps 4 m + 1, where the map is rebuilt.
#   U    mr1d-sod and mr2d-shock (multirate), solutions every sol-dtime = dt_out and time-threshold = 10 dt_out +
#        dt_out/3: the TIME of the k-th numbered solution is k*dt_out bit for bit (the product in double precision)
#        and that of the final solution is time-threshold bit for bit; a restart between two solutions gives the
#        OUTPUT tree of the run in one go.
#   C    in the runs of U: one row of the residual history per macro step, the cut ones included (the iteration of
#        the last row = the number of rows).
#   CONS mr2d-shock, recompute-every 5, 400 macro steps with restarts after 100, 200 and 300: mass drift within 1e-12
#        (energy reported).
#   D    mr1d-blast with max-level 10, whose depth follows the octaves of the local steps and changes as the front moves
#        (at least two depths printed in the run in one go, and cells at level 0 in each of its macro steps):
#        100 = 50 + restart + 50.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
CHECK=../../checks/restart/check_restart_vtk.py
rm -rf cases runs
mr_cases cases mr1d-sod mr2d-box-4b mr2d-plate-sa mr2d-shock mr1d-blast || { echo "mr-restart: FAIL - cannot build the cases"; exit 1; }
status=0
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute")
SR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" "MOSE-Multirate|levels|recompute")
OUT=("MOSE-IO|res-diter|1" "MOSE-IO|sol-overwrite|false")

# part DIR -> one run of MOSE in DIR (input.ini set), appending to solve.log; 0 if exit 0 and the threads asked
part() {
  ( cd "$1" && "$MR_MASTER" >> solve.log 2>&1 ) || { echo "  run $1: exit non-zero"; tail -3 "$1/solve.log"; return 1; }
  mr_threads_ok "$1" "$OMP_NUM_THREADS"
}

# same LABEL FULL SPLIT [FINAL] -> the two OUTPUT trees identical, residual histories not empty; FINAL is the final
# solution of the run in one go (default field.vtm)
same() {
  local lab=$1 a=$2 b=$3 fin=${4:-field.vtm}
  [[ -s $a/OUTPUT/residual-history.dat && -s $b/OUTPUT/residual-history.dat && -f $a/OUTPUT/$fin ]] \
    || { echo "  $lab: missing final solution or empty residual history"; return 1; }
  if python3 -B "$CHECK" "$a/OUTPUT" "$b/OUTPUT" > "$b/../compare.log" 2>&1; then
    echo "  $lab: identical ($(cd "$a/OUTPUT" && find . -type f ! -name '*.log' | wc -l) files)"
  else
    echo "  $lab: DIFFERENT"; sed 's/^/    /' "$b/../compare.log" | head -6; return 1
  fi
}

# split LABEL CASE N FIRST KEYS... -> N macro steps in one go against FIRST + restart + (N - FIRST); FINAL names the
# final solution (field.vtm; field.tec in the Tecplot ASCII runs of R-tec)
FINAL=field.vtm
split() {
  local lab=$1 case=$2 n=$3 first=$4; shift 4
  local d=runs/$lab
  mr_prepare "$d/full" "cases/$case" "$@" "${OUT[@]}" "MOSE-Parameters|iter-threshold|$n" || return 1
  part "$d/full" || return 1
  mr_prepare "$d/split" "cases/$case" "$@" "${OUT[@]}" "MOSE-Parameters|iter-threshold|$first" || return 1
  part "$d/split" || return 1
  python3 -B "$MR_COMMON/set_ini.py" "$d/split/input.ini" "MOSE-Parameters|newrun|false" \
    "MOSE-Parameters|iter-threshold|$((n - first))" || return 1
  part "$d/split" || return 1
  same "$lab ($n = $first + restart + $((n - first)))" "$d/full" "$d/split" "$FINAL"
}

echo "  R: restart in the middle of the run"
for c in mr1d-sod mr2d-box-4b mr2d-plate-sa; do
  split "R-$c-mr" "$c" 100 50 "${MR[@]}" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-diter|30" || status=1
  split "R-$c-sr" "$c" 100 50 "${SR[@]}" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-diter|30" || status=1
done

echo "  R-tec: restart from a Tecplot ASCII solution"
FINAL=field.tec
for c in mr1d-sod mr2d-box-4b; do
  for p in mr sr; do
    if [[ $p == mr ]]; then keys=("${MR[@]}"); else keys=("${SR[@]}"); fi
    d=runs/R-tec-$c-$p
    split "R-tec-$c-$p" "$c" 100 50 "${keys[@]}" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-diter|30" \
      "MOSE-IO|sol-format|tecplot ascii" || { status=1; continue; }
    # the second part started from the Tecplot solution of the first one, and no run wrote a vtk file
    nread=$(command grep -a -cE 'Initial conditions file +OUTPUT/field\.tec *$' "$d/split/solve.log")
    nvtk=$(find "$d/full/OUTPUT" "$d/split/OUTPUT" \( -name '*.vtm' -o -name '*.vts' \) | wc -l)
    if [[ $nread -eq 1 && $nvtk -eq 0 ]]; then
      echo "    the restart read OUTPUT/field.tec; no vtk file"
    else
      echo "  R-tec-$c-$p: restarts from OUTPUT/field.tec $nread (1 expected), vtk files $nvtk (0 expected)"; status=1
    fi
  done
done
FINAL=field.vtm

echo "  R2: map recomputed every k macro steps"
split R2-mr2d-shock mr2d-shock 100 50 "${MR[@]}" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-diter|30" || status=1
split R2-sod-k4-48 mr1d-sod 100 48 "${MR[@]}" "MOSE-Multirate|recompute-every|4" "MOSE-IO|sol-diter|30" || status=1
if split R2-sod-k4-50 mr1d-sod 100 50 "${MR[@]}" "MOSE-Multirate|recompute-every|4" "MOSE-IO|sol-diter|30" > runs/r2-50.txt; then
  echo "  R2-sod-k4-50 (restart not on a multiple of 4, informative): identical"
else
  echo "  R2-sod-k4-50 (restart not on a multiple of 4, informative): different, as declared"
fi
mr_run runs/R2-sod-k4-map cases/mr1d-sod "${MR[@]}" "MOSE-Multirate|recompute-every|4" "MOSE-Multirate|nonlocal-cfl|off" "MOSE-IO|shell-diter|1" \
       "MOSE-Parameters|iter-threshold|100" || status=1
if awk '/Multirate: cells per level/ {n++; if (n > 1 && $0 != p) {c++; if ((n - 1) % 4) b++}; p = $0} END {exit !(c > 0 && !b)}' runs/R2-sod-k4-map/solve.log; then
  echo "  R2-sod-k4-map: the cells per level change, and only at the macro steps 4 m + 1"
else
  echo "  R2-sod-k4-map: FAIL - the cells per level do not follow recompute-every 4"; status=1
fi

echo "  U and C: output instants k*dt_out and time-threshold, residual rows"
for spec in "mr1d-sod 2.5e-5" "mr2d-shock 2.0e-5"; do
  set -- $spec; c=$1; dto=$2
  tt=$(python3 -c "print(repr(10 * $dto + $dto / 3))")
  keys=("${MR[@]}" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-dtime|$dto" "MOSE-Parameters|time-threshold|$tt")
  d=runs/U-$c
  mr_prepare "$d/full" "cases/$c" "${keys[@]}" "${OUT[@]}" "MOSE-Parameters|iter-threshold|1000000" || { status=1; continue; }
  part "$d/full" || { status=1; continue; }
  # the restart point: the first macro step from the middle of the run on that does not end on an output instant
  first=$(python3 - "$d/full/OUTPUT/residual-history.dat" "$dto" <<'EOF'
import sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip()]
dto = float(sys.argv[2])
for i in range(len(rows) // 2, len(rows)):
    t = float(rows[i - 1][1])
    if all(abs(t - k * dto) > 1e-9 * t for k in range(1, 12)):
        print(i)
        break
EOF
)
  [[ -n $first ]] || { echo "  U-$c: no restart point between two solutions"; status=1; continue; }
  mr_prepare "$d/split" "cases/$c" "${keys[@]}" "${OUT[@]}" "MOSE-Parameters|iter-threshold|$first" || { status=1; continue; }
  part "$d/split" || { status=1; continue; }
  python3 -B "$MR_COMMON/set_ini.py" "$d/split/input.ini" "MOSE-Parameters|newrun|false" \
    "MOSE-Parameters|iter-threshold|1000000" || { status=1; continue; }
  part "$d/split" || { status=1; continue; }
  same "U-$c (restart after $first macro steps)" "$d/full" "$d/split" || status=1
  python3 -B - "$d/full/OUTPUT" "$dto" "$tt" "$first" <<'EOF' || status=1
import sys
from pathlib import Path
sys.path.insert(0, '../../checks/restart/restart-counters')
from vtk_time import time_of
out, dto, tt, first = Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
num = sorted(out.glob('field[0-9]*.vtm'), key=lambda p: int(p.stem[5:]))
ok = len(num) == 10
for p in num:
    k = int(p.stem[5:])
    if time_of(p) != float(k) * dto:
        ok = False
        print('    %s: TIME %.17g, k*dt_out %.17g' % (p.name, time_of(p), float(k) * dto))
tf = time_of(out / 'field.vtm')
ok = ok and tf == tt
rows = [l.split() for l in (out / 'residual-history.dat').read_text().splitlines() if l.strip()]
nrow, last = len(rows), int(rows[-1][0])
t_first = float(rows[first - 1][1])
between = all(abs(t_first - float(k) * dto) > 0 for k in range(1, 11)) and dto < t_first < 10 * dto
print('  U: %d numbered solutions at k*dt_out bit for bit: %s; final TIME = time-threshold bit for bit: %s; the restart '
      'falls between two solutions (t = %.6g): %s' % (len(num), 'yes' if len(num) == 10 and ok else 'NO',
      'yes' if tf == tt else 'NO (%.17g)' % tf, t_first, 'yes' if between else 'NO'))
print('  C: residual rows %d, iteration of the last row %d: %s' % (nrow, last, 'equal' if nrow == last else 'DIFFERENT'))
sys.exit(0 if ok and between and nrow == last else 1)
EOF
done

echo "  CONS: 400 macro steps with three restarts"
d=runs/CONS
mr_prepare "$d" cases/mr2d-shock "${MR[@]}" "MOSE-Multirate|recompute-every|5" "${OUT[@]}" "MOSE-Parameters|iter-threshold|100" \
  && part "$d" || status=1
for r in 1 2 3; do
  python3 -B "$MR_COMMON/set_ini.py" "$d/input.ini" "MOSE-Parameters|newrun|false" && part "$d" || status=1
done
rows=$(command grep -a -c . "$d/OUTPUT/residual-history.dat" 2>/dev/null)
echo "  CONS: residual rows $rows (400 expected)"
[[ $rows -eq 400 ]] || status=1
python3 -B ../mr-common/verify_conservation.py cases/mr2d-shock "$d" --quantities mass --max 1e-12 || status=1
python3 -B ../mr-common/verify_conservation.py cases/mr2d-shock "$d" --quantities energy --report

echo "  D: a depth that changes from one macro step to the next"
split D-mr1d-blast mr1d-blast 100 50 "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|10" "MOSE-Multirate|safety|2.0" \
  "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" "MOSE-IO|sol-diter|30" "MOSE-IO|shell-diter|1" || status=1
depths=$(command grep -a 'Multirate: H =' runs/D-mr1d-blast/full/solve.log | sed -E 's/.*L = +([0-9]+),.*/\1/' | sort -n | uniq | tr '\n' ' ')
if [[ $(wc -w <<< "$depths") -ge 2 ]]; then echo "  D-mr1d-blast: depths printed: $depths"
else echo "  D-mr1d-blast: one depth only ($depths): the depth does not change in the run"; status=1; fi
n0=$(command grep -a -c 'Multirate: cells per level' runs/D-mr1d-blast/full/solve.log)
z0=$(command grep -a 'Multirate: cells per level' runs/D-mr1d-blast/full/solve.log | awk '$5 + 0 == 0' | wc -l)
if [[ $n0 -gt 0 && $z0 -eq 0 ]]; then echo "  D-mr1d-blast: level 0 occupied in the $n0 macro steps of the run in one go"
else echo "  D-mr1d-blast: $z0 of $n0 macro steps of the run in one go with level 0 empty"; status=1; fi

if [[ $status -eq 0 ]]; then echo "mr-restart: ok"; else echo "mr-restart: FAIL"; fi
exit $status
