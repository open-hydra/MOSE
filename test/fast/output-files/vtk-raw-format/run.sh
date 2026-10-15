#!/usr/bin/env bash
#===============================================================================
#  vtk raw solutions are raw (fast tier)
#
#  With sol-format = vtk raw every .vts file must hold its data as raw binary
#  appended to the file (format="appended"), as sol-format asks and as the wall
#  and diagnostic outputs already do. Two runs of the Sod tube of
#  test/1D/Sod79, 10 iterations each, a numbered solution every 3:
#     steady     local time stepping, with the diagnostic output after each
#                solution: the solutions before and after the first diagnostic
#     transient  time-accurate, no diagnostic output
#  check_vtk_raw_format.py reads the XML header of every block file of every
#  .vtm in OUTPUT/ and asks for format="appended" in every DataArray.
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

fail() { echo "VtkRawFormat: FAIL — $*"; exit 1; }

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt

status=0
for case in steady transient; do
  dir=runs/$case
  rm -rf "$dir"; mkdir -p "$dir/INPUT" "$dir/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt; do ln -s "../../../INPUT/$f" "$dir/INPUT/$f"; done
  if [[ $case == steady ]]; then accurate=false; cfl=0.5; else accurate=true; cfl=0.9; fi
  printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 10' 'res-threshold = 1e-300' 'newrun = true' '' \
                '[MOSE-Physics]' 'equations = euler' '' \
                '[MOSE-Numerics]' "cfl = $cfl" "time-accurate = $accurate" 'time-scheme = RK3' \
                'integration-variables = cons' 'space-reconstruction = MUSCL' 'flux-limiter = vanleer' \
                'riemann-solver = HLLC' '' \
                '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 3' 'sol-overwrite = false' \
                'shell-diter = 1000000000' > "$dir/input.ini"
  ( cd "$dir" && "../../$LOCAL" > solve.log 2>&1 ) || { tail -5 "$dir/solve.log"; fail "solver exited non-zero ($case)"; }
  if command grep -q -E 'ERROR|NaN or p<0' "$dir/solve.log"; then
    command grep -E 'ERROR|NaN or p<0' "$dir/solve.log" | head -3; fail "the solver stopped with an error ($case)"
  fi
  python3 -B check_vtk_raw_format.py "$dir/OUTPUT" "$case" || status=1
done

[[ $status -eq 0 ]] || fail "a vtk raw solution is not raw"
echo "VtkRawFormat: ok — every vtk raw output of the steady and transient runs is raw appended data"
exit 0
