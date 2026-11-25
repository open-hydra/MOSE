#!/usr/bin/env bash
#===============================================================================
#  The FAS cycle, and Spalart-Allmaras in it, from a given iteration
#
#  cycle-start sets the iteration of the running level from which the cycle
#  acts; cycle-turbulence = .true. solves the RANS equations on the coarse
#  levels of the cycle too, from cycle-turbulence-start iterations after
#  cycle-start. The case of sa-fas of MPIBitIdentity
#  (../../../mpi/bit-identity/make_case.py SA --levels 3: a 3-D box of 2 x 2
#  blocks of 16 x 16 x 4 cells, walls, inlets, outlets, Spalart-Allmaras) on
#  three grid levels, 5 iterations on each coarse level, then 30 on the finest
#  with cycle = w and the other keys of the cycle at their defaults (a cycle
#  every second iteration), one thread, the solutions of each level at its
#  end in vtk raw. Start of the cycle:
#    none    cycle = none;
#    cnever  cycle-start = 1000000, beyond the run: every output file
#            identical to that of none, byte for byte;
#    clate   cycle-start = 12: no cycle on level 2 (5 iterations) and in the
#            iterations 1 to 11 of the finest level, a cycle from 12 on; the
#            residual history differs from that of none;
#    csplit  cycle = none for 10 iterations of the finest level, then a
#            restart from its solution with cycle = w for 20 more: its first
#            cycle is that of 12, and every output file must be identical to
#            that of clate, byte for byte.
#  Start of the RANS equations in the cycle (cycle-start left out):
#    frozen  cycle-turbulence = .false.;
#    never   cycle-turbulence = .true., cycle-turbulence-start = 1000000,
#            beyond the run: every output file identical to that of frozen,
#            byte for byte;
#    first   cycle-turbulence = .true., no start key: the residual history
#            differs from that of frozen, so that the RANS equations reach
#            the coarse levels;
#    zero    cycle-turbulence = .true., cycle-turbulence-start = 0 (the
#            default): every output file identical to that of first;
#    late    cycle-turbulence = .true., cycle-turbulence-start = 12: frozen
#            in the cycles of the iterations 2 to 10 of the finest level (and
#            in those of level 2, 5 iterations), solved from the cycle of 12
#            on, the start itself (a start acting one cycle late, from 14,
#            would differ from split); the residual history differs from those
#            of frozen and first;
#    split   cycle-turbulence = .false. for 10 iterations of the finest level,
#            then a restart from its solution with cycle-turbulence = .true.
#            and no start key for 20 more (level1-iter counts the iterations
#            since the launch): its first cycle with the RANS equations is that
#            of 12, as in late, and every output file must be identical to
#            that of late, byte for byte;
#    trel    cycle-start = 12, cycle-turbulence = .true.,
#            cycle-turbulence-start = 4: the RANS equations in the cycle from
#            the iteration 16 (12 + 4);
#    tsplit  cycle-start = 12, cycle-turbulence = .false. for 14 iterations of
#            the finest level, then a restart with cycle-turbulence = .true.
#            and no start key for 16 more: every output file identical to that
#            of trel, byte for byte (a start counted from the first iteration
#            of the level, 4, would make trel differ from it).
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
CASE=../../../mpi/bit-identity/make_case.py
COMMON=../../../common/Air-viscous

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=1

fail() { echo "MultigridStart: FAIL — $*"; exit 1; }

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$COMMON/transport.dat" INPUT/transport.dat
python3 -B "$CASE" SA --outdir . --levels 3 > /dev/null || fail "case generation"

# write_ini <dir> <newrun> <level1-iter> <lines of [MOSE-Multigrid] after the levels>
write_ini() {
  local d=$1 newrun=$2 iters=$3; shift 3
  printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 1000000000' 'res-threshold = 1e-300' "newrun = $newrun" '' \
    '[MOSE-Physics]' 'equations = navier-stokes' 'turbulence = SA' '' \
    '[MOSE-Numerics]' 'cfl = 0.5' 'time-accurate = false' 'time-scheme = RK3' \
    'integration-variables = cons' 'space-reconstruction = MUSCL' 'flux-limiter = vanleer' \
    'riemann-solver = HLLC' 'irs-beta = 0.5' '' \
    '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 1000000000' 'shell-diter = 1000000000' 'ini-diter = 1000000000' '' \
    '[MOSE-Multigrid]' 'levels = 3' "level1-iter = $iters" 'level2-iter = 5' 'level3-iter = 5' \
    "$@" > "$d/input.ini"
}

# solve <dir> <log>: the run in <dir>
solve() {
  ( cd "$1" && "../../$LOCAL" > "$2" 2>&1 ) || { tail -5 "$1/$2"; fail "solver exited non-zero ($1)"; }
  ls "$1"/OUTPUT/vtk/fieldBlock*.vts > /dev/null 2>&1 || fail "no field output ($1)"
}

# fresh <name>: an empty run folder runs/<name> with the inputs of the case
fresh() {
  local d=runs/$1 f
  rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT/vtk"
  for f in phase.txt thermo.dat transport.dat ic.tec bc.txt bc2.txt bc3.txt; do ln -s "../../../INPUT/$f" "$d/INPUT/$f"; done
}

# same <a> <b> <what>: every output file of runs/<a> identical to that of runs/<b>, and the same list of files
same() {
  local n=0 f
  diff <(cd "runs/$1/OUTPUT" && find . -type f | sort) <(cd "runs/$2/OUTPUT" && find . -type f | sort) > /dev/null \
    || fail "$3: other output files than $2"
  while IFS= read -r f; do
    n=$((n + 1))
    cmp -s "runs/$1/OUTPUT/$f" "runs/$2/OUTPUT/$f" || fail "$3: ${f#./} differs from that of $2"
  done < <(cd "runs/$1/OUTPUT" && find . -type f | sort)
  [[ $n -gt 0 ]] || fail "$3: no output file compared"
  echo "  $1 = $2: $n files identical"
}

# differ <a> <b> <what>: the residual history of runs/<a> not that of runs/<b>
differ() {
  cmp -s "runs/$1/OUTPUT/residual-history.dat" "runs/$2/OUTPUT/residual-history.dat" \
    && fail "$3: the residual history is that of $2"
  echo "  $1 != $2: the residual histories differ"
}

W='cycle = w'
T='cycle-turbulence = .true.'
for r in none cnever clate frozen never first zero late trel; do fresh "$r"; done
write_ini runs/none   true 30 'cycle = none';                                  solve runs/none   solve.log
write_ini runs/cnever true 30 "$W" 'cycle-start = 1000000';                    solve runs/cnever solve.log
write_ini runs/clate  true 30 "$W" 'cycle-start = 12';                         solve runs/clate  solve.log
fresh csplit
write_ini runs/csplit true 10 'cycle = none';                                  solve runs/csplit solve_none.log
write_ini runs/csplit false 20 "$W";                                           solve runs/csplit solve_restart.log
write_ini runs/frozen true 30 "$W" 'cycle-turbulence = .false.';               solve runs/frozen solve.log
write_ini runs/never  true 30 "$W" "$T" 'cycle-turbulence-start = 1000000';    solve runs/never  solve.log
write_ini runs/first  true 30 "$W" "$T";                                       solve runs/first  solve.log
write_ini runs/zero   true 30 "$W" "$T" 'cycle-turbulence-start = 0';          solve runs/zero   solve.log
write_ini runs/late   true 30 "$W" "$T" 'cycle-turbulence-start = 12';         solve runs/late   solve.log
fresh split
write_ini runs/split  true 10 "$W" 'cycle-turbulence = .false.';               solve runs/split  solve_frozen.log
write_ini runs/split false 20 "$W" "$T";                                       solve runs/split  solve_restart.log
write_ini runs/trel   true 30 "$W" 'cycle-start = 12' "$T" 'cycle-turbulence-start = 4'
solve runs/trel solve.log
fresh tsplit
write_ini runs/tsplit true 14 "$W" 'cycle-start = 12' 'cycle-turbulence = .false.'
solve runs/tsplit solve_frozen.log
write_ini runs/tsplit false 16 "$W" 'cycle-start = 12' "$T"
solve runs/tsplit solve_restart.log

same   cnever none  'cycle-start beyond the run'
differ clate  none  'cycle-start = 12'
same   csplit clate 'cycle = none to 10, then restarted with cycle = w'
same   never frozen 'start beyond the run'
differ first frozen 'cycle-turbulence = .true. from the first iteration'
same   zero  first  'cycle-turbulence-start = 0'
differ late  frozen 'cycle-turbulence-start = 12'
differ late  first  'cycle-turbulence-start = 12'
same   split late   'frozen to 10, then restarted with cycle-turbulence = .true.'
same   tsplit trel  'cycle-start = 12, frozen to 14, then restarted with cycle-turbulence = .true.'
echo "MultigridStart: ok — cycle-start beyond the run is cycle = none, and 12 the run without the cycle to 10 restarted with it, byte for byte; cycle-turbulence-start beyond the run is cycle-turbulence = .false., 0 the first iteration, 12 the run frozen to 10 restarted with the RANS equations in the cycle, byte for byte, and it counts from cycle-start"
exit 0
