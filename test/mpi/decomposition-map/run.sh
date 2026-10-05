#!/usr/bin/env bash
#===============================================================================
#  Block-to-rank assignment from INPUT/decomposition.map under MPI
#  (needs-mpi tier)
#
#  The optional file INPUT/decomposition.map, in the format the ATLAS mesh
#  decomposer (MDB) writes, gives the MPI rank of every block; without it
#  MOSE assigns the blocks itself, largest first, each to the least loaded
#  rank. The assignment must not change the solution: every run below gives
#  the field of one rank without a map, bit for bit.
#   1. A strip of 4 x 1 blocks: 2 ranks without a map; 2 ranks with the map
#      of MOSE's own assignment, written as MDB writes it; 2 ranks with the
#      strip cut in halves (0 0 1 1); 2 ranks with three blocks on rank 0
#      (the log reports the balance, 66.7 % of ideal); 2 ranks with every
#      block on rank 1 (rank 0 idle, with a warning); 4 ranks in reverse
#      order (3 2 1 0).
#   2. A square of 2 x 2 blocks around one vertex: 2 ranks with the diagonal
#      blocks together (0 1 1 0); 4 ranks in reverse order.
#   3. 110 blocks (11 x 10 of 4 x 4 cells): 2 ranks cyclic (0 1 0 1 ...) and
#      4 ranks with MOSE's own assignment. Output in Tecplot ASCII, so this
#      comparison holds to its 15 significant digits.
#  A map for other ranks stops the run on every rank with the same message:
#  a block on rank 2 of a 2-rank run; every block on rank 0 of a 2-rank run.
#  The residual history is not compared: its sums over the ranks follow the
#  assignment. The cases come from ../make_case.py, the maps from
#  ../write_map.py; nothing else is stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../common/Air-viscous
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}
MAP=INPUT/decomposition.map

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$COMMON/transport.dat" INPUT/transport.dat

# write_ini <sol-format>
write_ini() {
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 20
res-threshold  = 1e-300
newrun         = true

[MOSE-Physics]
equations = navier-stokes
schmidt   = 1.0

[MOSE-Numerics]
cfl                   = 0.3
vnn                   = 0.3
time-accurate         = true
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format  = $1
sol-diter   = 1000000000
shell-diter = 1000000000
ini-diter   = 1000000000
EOF
}

fail() { echo "MPIDecompositionMap: FAIL — $*"; exit 1; }

# run <ranks> <label>: solve and keep OUTPUT as runs/<label>; status in RC
run() {
  local np="$1" label="$2"
  rm -rf OUTPUT "runs/$label"
  mkdir -p OUTPUT/vtk runs
  OMP_NUM_THREADS=1 "$MPIEXEC" "$NPFLAG" "$np" "$LOCAL" > "solve_$label.log" 2>&1
  RC=$?
  mv OUTPUT "runs/$label"
}

# same <label-a> <label-b> <glob...>: every listed output file byte-identical
same() {
  local a="$1" b="$2" f x n=0; shift 2
  for f in "$@"; do
    for x in runs/"$a"/$f; do
      [[ -f "$x" ]] || continue
      n=$((n + 1))
      cmp -s "$x" "runs/$b/${x#runs/$a/}" || { echo "  differs: ${x#runs/$a/} ($a vs $b)"; return 1; }
    done
  done
  [[ $n -gt 0 ]] || { echo "  nothing compared ($a vs $b)"; return 1; }
  echo "  identical: $a vs $b ($n files)"
}

# with_map <ranks> <label> <owners> <reference> <glob> [<text expected in the log>]
with_map() {
  local np="$1" label="$2" owners="$3" ref="$4" glob="$5" text="${6:-}"
  python3 -B ../write_map.py $CASE --owners "$owners" > /dev/null || fail "map generation ($owners)"
  run "$np" "$label"
  rm -f "$MAP"
  if [[ $RC -ne 0 ]]; then
    echo "  the run $label failed (solve_$label.log)"; status=1; return
  fi
  same "$ref" "$label" "$glob" || status=1
  if ! grep -aq "MPI partition read from $MAP" "solve_$label.log"; then
    echo "  $label: the setup log does not name $MAP"; status=1
  fi
  if [[ -n "$text" ]] && ! grep -aqF "$text" "solve_$label.log"; then
    echo "  $label: the setup log does not say \"$text\""; status=1
  fi
}

# stops <ranks> <label> <owners> <text>: the run must stop, every rank printing <text>
stops() {
  local np="$1" label="$2" owners="$3" text="$4" k
  python3 -B ../write_map.py $CASE --owners "$owners" > /dev/null || fail "map generation ($owners)"
  run "$np" "$label"
  rm -f "$MAP"
  k=$(grep -acF "$text" "solve_$label.log")
  if [[ $RC -ne 0 && $k -ge $np ]]; then
    echo "  stops on every rank: $owners on $np ranks -> \"$text\""
  else
    echo "  $label: exit $RC, \"$text\" printed $k times for $np ranks (solve_$label.log)"; status=1
  fi
}

status=0
rm -f "$MAP"

# --- 1. strip of 4 x 1 blocks ------------------------------------------------
CASE="--layout 4x1 --cells 16 16"
write_ini "vtk raw"
python3 -B ../make_case.py $CASE --bc channel --zfaces 300 > /dev/null || fail "case generation (strip)"
run 1 strip-ref
[[ $RC -eq 0 ]] || fail "the reference run failed (solve_strip-ref.log)"
run 2 strip-nomap
[[ $RC -eq 0 ]] || fail "the run without a map failed (solve_strip-nomap.log)"
same strip-ref strip-nomap 'vtk/*.vts' || status=1
with_map 2 strip-lpt    lpt:2          strip-ref 'vtk/*.vts'
with_map 2 strip-halves list:0,0,1,1   strip-ref 'vtk/*.vts'
with_map 2 strip-3to1   list:0,0,0,1   strip-ref 'vtk/*.vts' "balance  66.7% of ideal"
with_map 2 strip-idle0  list:1,1,1,1   strip-ref 'vtk/*.vts' "WARNING: the decomposition map gives no block to 1 of 2 ranks"
with_map 4 strip-rev4   reverse:4      strip-ref 'vtk/*.vts'
stops 2 strip-rank2 list:0,1,2,0 "block 3 is on rank 2, outside the ranks of the run (0 to 1)"
stops 2 strip-for1  list:0,0,0,0 "the largest rank is 0: the map is for 1 rank, the run has 2"

# --- 2. square of 2 x 2 blocks around one vertex -------------------------------
CASE="--layout 2x2 --cells 16 16"
python3 -B ../make_case.py $CASE --bc box --zfaces 300 > /dev/null || fail "case generation (square)"
run 1 square-ref
[[ $RC -eq 0 ]] || fail "the reference run failed (solve_square-ref.log)"
with_map 2 square-diag list:0,1,1,0 square-ref 'vtk/*.vts'
with_map 4 square-rev4 reverse:4    square-ref 'vtk/*.vts'

# --- 3. 110 blocks -------------------------------------------------------------
CASE="--layout 11x10 --cells 4 4"
write_ini "tecplot ascii"
python3 -B ../make_case.py $CASE --bc channel --zfaces 300 > /dev/null || fail "case generation (110 blocks)"
run 1 many-ref
[[ $RC -eq 0 ]] || fail "the reference run failed (solve_many-ref.log)"
with_map 2 many-cyclic cyclic:2 many-ref 'field*.tec'
with_map 4 many-lpt    lpt:4    many-ref 'field*.tec'

if [[ $status -ne 0 ]]; then
  echo "MPIDecompositionMap: FAIL"
  exit 1
fi
echo "MPIDecompositionMap: ok — every assignment gives the one-rank field; maps for other ranks stop every rank"
exit 0
