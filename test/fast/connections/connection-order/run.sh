#!/usr/bin/env bash
#===============================================================================
#  One flux per face between blocks, whatever the order of bc.txt (fast tier)
#
#  Two 2-D blocks cut at y = 0.5 with their index directions permuted across
#  the cut, so that the connection pairs face 4 of one block with face 2 of
#  the other: two faces of the same parity (make_case.py). Of the two records
#  of such a face, the one that computes its flux must not depend on the order
#  of the records in bc.txt: with the shock sensor each side reads the switch
#  of its own cell, so another choice gives other bits. The same case with
#  the records of bc.txt sorted and in reverse order, Navier-Stokes with
#  HLLC+Tramel and a shock leaving the cut, 20 RK3 steps: the two solutions
#  must be the same bit for bit. The ghost cells of a connection face are
#  the first two cells of the block across it: with block 1 one cell thick
#  across the cut (thin), or two cells thick with two multigrid levels, one
#  at level 2 (thinmg), the run must stop at set-up with a non-zero exit code
#  and a message that names the record, the block across and its cells.
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../common/Air-viscous

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=1

mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

cat > input.ini <<'EOT'
[MOSE-Parameters]
iter-threshold = 20
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC+Tramel

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
EOT

# run <order>: solve the case with the bc.txt records in that order, keep OUTPUT as runs/<order>
run() {
  python3 -B make_case.py case "$1" > /dev/null || exit 1
  rm -rf OUTPUT "runs/$1"
  mkdir -p OUTPUT/vtk
  "$LOCAL" > "runs/$1.log" 2>&1
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|ERROR' "runs/$1.log"; then
    echo "  [FAIL] the run $1 did not complete (see runs/$1.log)"
    echo "ConnectionOrder: FAILED"
    exit 1
  fi
  mv OUTPUT "runs/$1"
}

run sorted
run reversed
echo "-- two blocks with permuted indices, records of bc.txt sorted and reversed --"
n=0; fail=0
for f in runs/sorted/vtk/*.vts; do
  n=$((n + 1))
  if ! cmp -s "$f" "runs/reversed/${f#runs/sorted/}"; then
    echo "  [FAIL] ${f#runs/sorted/} differs"; fail=1
  fi
done
[[ $n -eq 2 ]] || { echo "  [FAIL] $n block files, 2 expected"; fail=1; }
[[ $fail -eq 0 ]] && echo "  [ok  ] the same solution bit for bit ($n block files)"

# thin <case> <level>: block 1 has one cell across the cut (at the multigrid level), fewer than the two ghost cells of
# a connection read from it; the run must stop at set-up with a non-zero exit code and the message that names it
thin() {
  python3 -B make_case.py case "$1" > /dev/null || exit 1
  rm -rf OUTPUT
  mkdir -p OUTPUT/vtk
  "$LOCAL" > "runs/$1.log" 2>&1
  local rc=$?
  local want="block 1 across it has 1 cell across its face 4$2, fewer than the 2 a connection needs"
  if [[ $rc -ne 0 ]] && grep -a -qF "$want" "runs/$1.log"; then
    echo "  [ok  ] $1: stopped at set-up, exit code $rc, \"$want\""
  else
    echo "  [FAIL] $1: exit code $rc, the run was not stopped with the message (see runs/$1.log)"
    fail=1
  fi
}
echo "-- block 1 one cell thick across the cut, on the mesh and at multigrid level 2 --"
thin thin ''
printf '%s\n' '' '[MOSE-Multigrid]' 'levels = 2' 'level1-iter = 20' 'level2-iter = 10' >> input.ini
thin thinmg ' at multigrid level 2'

if [[ $fail -ne 0 ]]; then
  echo "ConnectionOrder: FAILED"
  exit 1
fi
echo "ConnectionOrder: all checks passed"
exit 0
