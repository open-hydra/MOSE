#!/usr/bin/env bash
#===============================================================================
#  Restart round-trip (fast tier, I3/F6)
#
#  Verifies that stopping a run, writing a solution and restarting (newrun=false)
#  reproduces the uninterrupted run:
#     A : N iterations in one go
#     B : N/2 iterations -> restart -> N/2 more
#  1. Tecplot ASCII (OUTPUT/field.tec): A and B must match to the round-trip
#     precision of the file (17 significant digits).
#  2. vtk raw: A and B must be identical bit for bit, counters included: the
#     runs write numbered solutions every 10 iterations (sol-overwrite = false)
#     and ramp the CFL number over the first 40 iterations, so the restart must
#     also carry on the iteration count, the output numbering and the
#     iteration column of residual-history.dat.
#  3. A restart from a solution that holds no time must stop at setup with a
#     message, in a time-accurate and in a steady run: the Tecplot solution of
#     part 1 without SOLUTIONTIME, and the vtk raw solution of part 2 without
#     the TIME field data of its block file.
#  4. Initial condition from a solution written by MOSE, in vtk raw, vtk ascii
#     and Tecplot ASCII: a new run (newrun = true) from the solution copied into
#     INPUT/ (field.vtm and its vtk/ folder as ic.vtm and vtk/, field.tec as
#     ic.tec) must have, after N/2 more steps, the cell values of the restart
#     from the same solution, bit for bit (check_initial.py).
#  Reuses the coarse Sod inputs (test/1D/Sod79); no new mesh/IC data is stored
#  here.
#     C : probes at a restart. Two probes on cells of the tube, by index (pA on
#         cell 50, pB on cell 45: the nodes of this 1-D mesh lie on the planes
#         j = 0 and k = 0), one row per iteration: 20 iterations in one go,
#         then 10 with the file of pA removed before a restart for 10 more.
#         The restart must end normally, create pA.txt again with the rows of
#         the restarted part and append those rows to pB.txt, and the rows of
#         both files must be those of the run in one go (check_probes.py; the
#         restarted run goes on with the iteration count of its solution). A
#         restart did not open the probes from a missing file on, and the run
#         ended at its first probe output with a run-time error.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air
N=50          # total iterations
HALF=25       # iterations before the restart

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt
SOLVER="$HERE/bin/MOSE"

# solve <directory> <log file> <label>: the solver must exit 0 and report no
# error (it exits 0 also when it stops on NaN or p<0, so the log is checked too)
solve() {
  ( cd "$1" && "$SOLVER" > "$2" 2>&1 ) || { tail -5 "$1/$2"; echo "Restart: FAIL — solver exited non-zero ($3)"; exit 1; }
  if command grep -q -E 'ERROR|NaN or p<0' "$1/$2"; then
    command grep -E 'ERROR|NaN or p<0' "$1/$2" | head -3
    echo "Restart: FAIL — the solver stopped with an error ($3)"; exit 1
  fi
}

write_ini() {  # <iter-threshold> <newrun>
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = $2

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
EOF
}

# --- A: uninterrupted N-iteration run --------------------------------------
write_ini "$N" true
rm -f OUTPUT/field*.tec
solve . OUTPUT/solve_full.log "Tecplot ASCII, full run"
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from full run"; exit 1; }
cp OUTPUT/field.tec field_full.tec

# --- B: N/2 then restart for the remaining N/2 -----------------------------
write_ini "$HALF" true
rm -f OUTPUT/field*.tec
solve . OUTPUT/solve_half.log "Tecplot ASCII, first half"
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from first half"; exit 1; }

write_ini "$HALF" false
solve . OUTPUT/solve_restart.log "Tecplot ASCII, restart"
[[ -f OUTPUT/field.tec ]] || { echo "Restart: FAIL — no field.tec from restart"; exit 1; }
cp OUTPUT/field.tec field_restart.tec

python3 -B check_restart.py field_full.tec field_restart.tec || exit 1

# --- 2. vtk raw, bit for bit -----------------------------------------------
write_ini_vtk() {  # <iter-threshold> <newrun>
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = $2

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.9
cfl-rise-threshold = 40
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = vtk raw
sol-diter = 10
sol-overwrite = false
EOF
}

prepare() {  # <run directory>
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt; do
    ln -s "../../../INPUT/$f" "$1/INPUT/$f"
  done
}

prepare runs/full
( cd runs/full && write_ini_vtk "$N" true ); solve runs/full solve_full.log "vtk raw, full run"
[[ -f runs/full/OUTPUT/field.vtm ]] || { echo "Restart (vtk raw): FAIL — no field.vtm from full run"; exit 1; }

prepare runs/restart
( cd runs/restart && write_ini_vtk "$HALF" true ); solve runs/restart solve_half.log "vtk raw, first half"
[[ -f runs/restart/OUTPUT/field.vtm ]] || { echo "Restart (vtk raw): FAIL — no field.vtm from first half"; exit 1; }
( cd runs/restart && write_ini_vtk "$HALF" false ); solve runs/restart solve_restart.log "vtk raw, restart"

python3 -B check_restart_vtk.py runs/full/OUTPUT runs/restart/OUTPUT || exit 1

# --- 3. A restart from a solution that holds no time -----------------------
status=0
for fmt in tecplot vtk; do
  for kind in time-accurate steady; do
    d=runs/no-time-$fmt-$kind
    prepare "$d"
    if [[ $fmt == tecplot ]]; then
      sed -e 's/, *SOLUTIONTIME *= *[^,]*//' -e 's/, *STRANDID *= *[0-9]*//' field_full.tec > "$d/OUTPUT/field.tec"
      ( cd "$d" && write_ini 1 false )
      file=OUTPUT/field.tec
    else
      # the vtk raw solution of part 2, with the TIME field data taken out of the XML header of its block
      # file (the other arrays keep their offsets into the appended data)
      python3 -B - runs/full/OUTPUT "$d/OUTPUT" <<'EOF' || { echo "  vtk: FAIL — the solution without TIME was not written"; status=1; continue; }
import re
import sys
from pathlib import Path
src, dst = Path(sys.argv[1]), Path(sys.argv[2])
vtm = (src / 'field.vtm').read_text()
(dst / 'field.vtm').write_text(vtm)
for entry in re.findall(r'file="([^"]+)"', vtm):
    head, sep, tail = (src / entry).read_bytes().partition(b'<AppendedData')
    head, n = re.subn(rb'\s*<DataArray[^>]*Name="TIME"[^>]*/>', b'', head)
    if n != 1:
        sys.exit(f'{entry}: {n} TIME field data arrays instead of 1')
    (dst / entry).parent.mkdir(parents=True, exist_ok=True)
    (dst / entry).write_bytes(head + sep + tail)
EOF
      ( cd "$d" && write_ini_vtk 1 false )
      file=OUTPUT/field.vtm
    fi
    if [[ $kind == steady ]]; then sed -i 's/^time-accurate = true$/time-accurate = false/' "$d/input.ini"; fi
    ( cd "$d" && "$SOLVER" > solve.log 2>&1 ); rc=$?
    if [[ $rc -ne 0 ]] && command grep -q "ERROR\] Restart from $file: the solution holds no time" "$d/solve.log" \
       && command grep -q 'STOP: restart from a solution without a time' "$d/solve.log" \
       && ! command grep -q -E 'NaN|Program received signal' "$d/solve.log"; then
      echo "  $fmt, $kind: ok — the restart from a solution without a time stops at setup"
    else
      echo "  $fmt, $kind: FAIL — the restart from a solution without a time did not stop with the message, exit code $rc (see $d/solve.log)"
      status=1
    fi
  done
done
command grep -h 'ERROR\]' runs/no-time-*/solve.log 2>/dev/null | sort -u | sed 's/^/   /'
[[ $status -eq 0 ]] || { echo "Restart: FAIL — a restart from a solution without a time does not stop"; exit 1; }

# --- 4. Initial condition from a solution written by MOSE ------------------
# The solution of N/2 steps, written in FORMAT, starts a new run (copied into INPUT/) and a restart; both go on for N/2
# steps at a constant CFL number, so they must have the same cell values. In vtk the .vtm lists its block files relative
# to its own folder (vtk/fieldBlock1.vts), and the copy keeps that layout.
write_ini_ic() {  # <newrun> <ic-format> <sol-format>
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $HALF
newrun = $1

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

[MOSE-IO]
ic-format = $2
sol-format = $3
EOF
}

# solved <directory> <log file>: like solve, but returns non-zero instead of stopping the test
solved() {
  ( cd "$1" && "$SOLVER" > "$2" 2>&1 ) || { tail -3 "$1/$2"; return 1; }
  ! command grep -q -E 'ERROR|NaN or p<0' "$1/$2"
}

ic_status=0
for fmt in "vtk raw" "vtk ascii" "tecplot ascii"; do
  tag=${fmt// /-}
  src=runs/ic-$tag-src
  new=runs/ic-$tag-new
  prepare "$src"
  ( cd "$src" && write_ini_ic true "tecplot ascii" "$fmt" ); solve "$src" solve_first.log "$fmt, first half"
  rm -rf "$new"; mkdir -p "$new/INPUT" "$new/OUTPUT"
  for f in phase.txt thermo.dat bc.txt; do
    ln -s "../../../INPUT/$f" "$new/INPUT/$f"
  done
  if [[ $fmt == vtk* ]]; then
    cp "$src/OUTPUT/field.vtm" "$new/INPUT/ic.vtm"
    cp -r "$src/OUTPUT/vtk" "$new/INPUT/vtk"
  else
    cp "$src/OUTPUT/field.tec" "$new/INPUT/ic.tec"
  fi
  ( cd "$new" && write_ini_ic true "$fmt" "$fmt" )
  if ! solved "$new" solve_new.log; then
    echo "  $fmt: FAIL — the new run from the solution did not complete (see $new/solve_new.log)"
    ic_status=1
    continue
  fi
  ( cd "$src" && write_ini_ic false "tecplot ascii" "$fmt" ); solve "$src" solve_restart.log "$fmt, restart"
  python3 -B check_initial.py "${fmt%% *}" "$src/OUTPUT" "$new/OUTPUT" "$fmt" || ic_status=1
done
if [[ $ic_status -ne 0 ]]; then
  echo "Restart: FAIL — a new run from a solution written by MOSE does not read it"
  exit 1
fi
echo "Restart: ok — Tecplot ASCII to its precision, vtk raw bit for bit, a solution without a time refused, and new runs from solutions written by MOSE"
status=$?

# --- C: probes at a restart, the file of one of them removed ----------------
add_probes() {
  cat >> input.ini <<EOF

[MOSE-Probes]
probe1 = pA
probe2 = pB

[pA]
variables = p u
diter = 1
index-position = 1 50 1 1

[pB]
variables = p
diter = 1
index-position = 1 45 1 1
EOF
}
write_ini 20 true
add_probes
rm -f OUTPUT/pA.txt OUTPUT/pB.txt
"$LOCAL" > OUTPUT/probes_full.log 2>&1
mv -f OUTPUT/pA.txt OUTPUT/pA_full.txt 2>/dev/null
mv -f OUTPUT/pB.txt OUTPUT/pB_full.txt 2>/dev/null
write_ini 10 true
add_probes
rm -f OUTPUT/pA.txt OUTPUT/pB.txt
"$LOCAL" > OUTPUT/probes_first.log 2>&1
rows() { if [[ -f $1 ]]; then wc -l < "$1"; else echo 0; fi; }
nA=$(rows OUTPUT/pA.txt)
nB=$(rows OUTPUT/pB.txt)
rm -f OUTPUT/pA.txt
write_ini 10 false
add_probes
"$LOCAL" > OUTPUT/probes_restart.log 2>&1
rc=$?
mA=$(rows OUTPUT/pA.txt)
mB=$(rows OUTPUT/pB.txt)
python3 -B check_probes.py 0 OUTPUT/pA_full.txt OUTPUT/pA.txt 0 OUTPUT/pB_full.txt OUTPUT/pB.txt "$nB"
same=$?
if [[ $rc -eq 0 && $nA -gt 0 && $mA -eq $nA && $mB -eq $((2 * nB)) && $same -eq 0 ]]; then
  echo "  [ok  ] probes at a restart: pA.txt created again ($mA rows), pB.txt appended ($nB + $nB rows), the rows of the run in one go"
else
  echo "  [FAIL] probes at a restart: exit status $rc; pA.txt $mA rows (before the restart $nA), pB.txt $mB rows" \
       "(before $nB), rows of the run in one go: $([[ $same -eq 0 ]] && echo yes || echo no); see OUTPUT/probes_restart.log"
  status=1
fi
# Parts of this test in their own folders, each with its case and checks
# described at the head of its run.sh, run in their folders as CTest runs
# this script:
for part in restart-counters; do
  echo "=== $part"
  if (cd "$HERE/$part" && ./run.sh); then
    echo "=== $part: PASS"
  else
    echo "=== $part: FAIL"
    status=1
  fi
done
exit $status
