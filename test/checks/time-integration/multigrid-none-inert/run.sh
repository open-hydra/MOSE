#!/usr/bin/env bash
#===============================================================================
#  Multigrid cycle off (fast tier)
#
#  cycle = none must be the solver without the FAS cycle, whatever the other
#  keys of the cycle say. The Sod tube of test/1D/Sod79 (100 cells), steady,
#  on two grid levels: 10 iterations of the coarse one, then 30 of the fine
#  one, solutions every 4 iterations in vtk raw.
#    plain   no key of the cycle;
#    inert   cycle = none, pre-iter = 3, post-iter = 2, coarse-iter = 7,
#            damp-restriction = 0.5, damp-prolongation = 0.5:
#            every output file identical to that of plain, byte for byte;
#    cycle   the same keys with cycle = w: residual-history.dat must differ
#            from that of plain, so that the keys reach the cycle when it is on.
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

fail() { echo "MultigridNoneInert: FAIL — $*"; exit 1; }

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt
# Boundary conditions of the coarse grid (50 cells): the same faces, with the cell indices of the coarse grid
awk '{ if ($5 == 2) $2 = 50; if ($5 >= 3 && $2 > 50) next; print }' "$SOD/INPUT/bc_x1.txt" > INPUT/bc2.txt

# solve <name> <lines of [MOSE-Multigrid] after the levels>: the run in runs/<name>
solve() {
  local d=runs/$1; shift
  rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
  for f in phase.txt thermo.dat ic.tec bc.txt bc2.txt; do ln -s "../../../INPUT/$f" "$d/INPUT/$f"; done
  printf '%s\n' '[MOSE-Parameters]' 'res-threshold = 1e-300' 'newrun = true' '' \
    '[MOSE-Physics]' 'equations = euler' '' \
    '[MOSE-Numerics]' 'cfl = 0.5' 'time-accurate = false' 'time-scheme = RK3' 'integration-variables = cons' \
    'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC' '' \
    '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 4' 'sol-overwrite = false' 'shell-diter = 1000000000' \
    'ini-diter = 1000000000' '' \
    '[MOSE-Multigrid]' 'levels = 2' 'level2-iter = 10' 'level1-iter = 30' "$@" > "$d/input.ini"
  ( cd "$d" && "../../$LOCAL" > solve.log 2>&1 ) || { tail -5 "$d/solve.log"; fail "solver exited non-zero ($d)"; }
}

KEYS=('pre-iter = 3' 'post-iter = 2' 'coarse-iter = 7' 'damp-restriction = 0.5' 'damp-prolongation = 0.5')
solve plain
solve inert 'cycle = none' "${KEYS[@]}"
solve cycle 'cycle = w' "${KEYS[@]}"

n=0
while IFS= read -r f; do
  n=$((n + 1))
  cmp -s "runs/plain/OUTPUT/$f" "runs/inert/OUTPUT/$f" || fail "cycle = none: ${f#./} differs from the run without the keys"
done < <(cd runs/plain/OUTPUT && find . -type f | sort)
diff <(cd runs/plain/OUTPUT && find . -type f | sort) <(cd runs/inert/OUTPUT && find . -type f | sort) > /dev/null \
  || fail "cycle = none: other output files than the run without the keys"
[[ $n -gt 0 ]] || fail "no output file compared"
grep -q 'w-cycle' runs/cycle/solve.log || fail "cycle = w: the set-up summary does not list the cycle"
cmp -s runs/plain/OUTPUT/residual-history.dat runs/cycle/OUTPUT/residual-history.dat \
  && fail "cycle = w: the residual history is that of the run without the cycle"
echo "MultigridNoneInert: ok — cycle = none with the other keys of the cycle set writes the $n output files of the run without them, byte for byte; cycle = w changes the residual history"
exit 0
