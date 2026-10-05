#!/usr/bin/env bash
#===============================================================================
#  Block-to-rank assignment from INPUT/decomposition.map (fast tier, 1 rank)
#
#  The optional file INPUT/decomposition.map, in the format the ATLAS mesh
#  decomposer (MDB) writes, gives the MPI rank of every block. On one rank:
#   1. a map with every block on rank 0 gives the bits of the run without a
#      map, and the setup log names the file the assignment comes from;
#   2. a map that does not fit the mesh or the run stops the solver at setup
#      with a message that names the fault: a block line missing or one too
#      many, blocks out of order, sizes or cell counts other than the mesh's,
#      a rank the run does not have, a map without the rank column.
#  Case: 4 x 1 blocks of 8 x 8 cells (../../../mpi/make_case.py); maps written
#  by ../../../mpi/write_map.py. The MPI tier runs the assignments on 2 and 4
#  ranks (test/mpi/decomposition-map).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../../common/Air-viscous
CASE="--layout 4x1 --cells 8 8"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=1

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$COMMON/transport.dat" INPUT/transport.dat

cat > input.ini <<'EOF'
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
sol-format  = vtk raw
sol-diter   = 1000000000
shell-diter = 1000000000
ini-diter   = 1000000000
EOF

fail() { echo "DecompositionMap: FAIL — $*"; exit 1; }
MAP=INPUT/decomposition.map

# run <label>: solve and keep OUTPUT as runs/<label>; status in RC
run() {
  local label="$1"
  rm -rf OUTPUT "runs/$label"
  mkdir -p OUTPUT/vtk runs
  "$LOCAL" > "solve_$label.log" 2>&1
  RC=$?
  mv OUTPUT "runs/$label"
}

# same <label-a> <label-b>: field outputs and residual history byte-identical
same() {
  local a="$1" b="$2" x n=0
  for x in runs/"$a"/vtk/*.vts runs/"$a"/residual-history.dat; do
    [[ -f "$x" ]] || continue
    n=$((n + 1))
    cmp -s "$x" "runs/$b/${x#runs/$a/}" || { echo "  differs: ${x#runs/$a/} ($a vs $b)"; return 1; }
  done
  [[ $n -gt 1 ]] || { echo "  nothing compared ($a vs $b)"; return 1; }
  echo "  identical: $a vs $b ($n files)"
}

status=0
python3 -B ../../../mpi/make_case.py $CASE --bc channel --zfaces 300 > /dev/null || fail "case generation"

# 1. no map, then every block on rank 0
rm -f "$MAP"
run nomap
[[ $RC -eq 0 ]] || fail "the run without a map failed (solve_nomap.log)"
python3 -B ../../../mpi/write_map.py $CASE --owners lpt:1 > /dev/null || fail "map generation"
run map
if [[ $RC -ne 0 ]]; then
  echo "  the run with a map on rank 0 failed (solve_map.log)"; status=1
else
  same nomap map || status=1
  if grep -aq "MPI partition read from $MAP" solve_map.log; then
    echo "  the setup log names $MAP"
  else
    echo "  the setup log does not name $MAP"; status=1
  fi
fi

# 2. maps that do not fit: the solver must stop at setup with the message
#    <label>|<write_map.py arguments>|<text expected in the log>
while IFS='|' read -r label args expected; do
  python3 -B ../../../mpi/write_map.py $CASE $args > /dev/null || fail "map generation ($args)"
  run "$label"
  if [[ $RC -ne 0 ]] && grep -aqF "$expected" "solve_$label.log"; then
    echo "  stops: $args -> \"$expected\""
  else
    echo "  does not stop with \"$expected\": $args (exit $RC, solve_$label.log)"; status=1
  fi
done <<'EOF'
missing-row|--owners lpt:1 --bad missing-row|3 block lines for the 4 blocks of the mesh
extra-row|--owners lpt:1 --bad extra-row|more block lines than the 4 blocks of the mesh
order|--owners lpt:1 --bad order|line 5 lists block 2 where block 1 is expected
dims|--owners lpt:1 --bad dims|block 2 is 9 x 8 x 1 cells in the map, 8 x 8 x 1 in the mesh
cells|--owners lpt:1 --bad cells|block 2 has 65 cells in the map, 64 in the mesh
rank|--owners list:0,1,0,0|block 2 is on rank 1, outside the ranks of the run (0 to 0)
no-rank|--owners lpt:1 --bad no-rank|line 5 is not 13 integers
EOF
rm -f "$MAP"

if [[ $status -ne 0 ]]; then
  echo "DecompositionMap: FAIL"
  exit 1
fi
echo "DecompositionMap: ok — a map on rank 0 gives the bits of the run without it; maps that do not fit stop the run"
exit 0
