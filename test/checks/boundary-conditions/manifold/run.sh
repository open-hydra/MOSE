#!/usr/bin/env bash
#===============================================================================
#  The manifold (BC 501): faces joined by a lumped plenum
#
#  A duct cut in two blocks (make_case.py): joined by a block connection
#  (101) or by a manifold (501), straight or with the second block turned
#  back (the two face-2 ends joined, the turn-around of a counter-flow element
#  in small). The uniform state of the inlet (T, g) and of the back pressure is
#  the exact steady solution of every layout, so:
#    - started from it, every layout keeps it (the manifold carries g, T0 and
#      p of a uniform flow exactly): the two-block runs equal the one-block run
#      to round-off, and the log line of the manifold shows the same mass flow
#      and total enthalpy on its two sides;
#    - started from a perturbed state (the pressure of block 1 raised by 2 %),
#      the manifold runs converge to the solution of the one-block run;
#    - the ratios of the link act: mdot-ratio = 0.5 halves the mass flux of
#      block 2, dh0 lowers its total temperature by dh0/cp, and the log line
#      shows the mass flow and total enthalpy in the relations of the record.
#  Euler, HLLC, MUSCL van Leer, steady (local time steps), gas: test/common/Air.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../../common/Air
ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$AIR/phase.txt"  INPUT/phase.txt
ln -sf "../$AIR/thermo.dat" INPUT/thermo.dat

# write_ini <iterations> <residual> <shell interval>: steady, local time steps, to the residual or the iterations (the
# exact state has a residual of zero: the run stops at once), the manifold line every <shell interval> iterations
write_ini() {
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = true
res-threshold = $2

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = $3
ini-diter = 1000000000
EOF
}

# run <tag>: solve the case in INPUT/, keep OUTPUT as runs/<tag> (with a marker for a reversed second block)
run() {
  rm -rf OUTPUT "runs/$1"
  mkdir -p OUTPUT/vtk
  "$LOCAL" > "runs/$1.log" 2>&1
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|ERROR' "runs/$1.log"; then
    echo "  [FAIL] the run $1 did not complete (see runs/$1.log)"
    return 1
  fi
  mv OUTPUT "runs/$1"
  [[ $1 == reverse* ]] && touch "runs/$1/reverse"
  return 0
}

fail=0

echo "-- the exact state kept: 500 iterations from it --"
write_ini 500 1e-300 1
for layout in one two manifold reverse; do
  python3 -B make_case.py case "$layout" || exit 1
  run "exact-$layout" || fail=1
done
for layout in two manifold reverse; do
  echo "  $layout against one block (round-off):"
  python3 -B make_case.py same runs/exact-one "runs/exact-$layout" 1e-12 || fail=1
done
for layout in manifold reverse; do
  echo "  the line of the manifold, $layout:"
  python3 -B make_case.py carried "runs/exact-$layout.log" "runs/exact-$layout" 1e-12 || fail=1
done

echo "-- the steady state found from a perturbed one: to a residual of 1e-13 --"
write_ini 60000 1e-13 500
for layout in one manifold reverse; do
  python3 -B make_case.py case "$layout" perturbed || exit 1
  run "perturbed-$layout" || fail=1
done
for layout in manifold reverse; do
  echo "  $layout against one block (1e-9):"
  python3 -B make_case.py same runs/perturbed-one "runs/perturbed-$layout" 1e-9 || fail=1
done

echo "-- the ratios of the link: mdot-ratio 0.5, K 0.5, p0-ratio 0.9, dh0 20 cp --"
write_ini 60000 1e-13 500
dh0=$(python3 -c "print(20.0 * 1004.5)")
python3 -B make_case.py case manifold mr=0.5 K=0.5 pr=0.9 "dh0=$dh0" || exit 1
run ratios || fail=1
echo "  the mass flux of block 2 is half that of block 1:"
python3 -B make_case.py split runs/ratios 0.5 || fail=1
echo "  the total temperature of block 2 is 20 K below that of block 1:"
python3 -B make_case.py drop runs/ratios "$dh0" || fail=1
echo "  the line of the manifold:"
python3 -B make_case.py carried runs/ratios.log runs/ratios 1e-8 mr=0.5 K=0.5 pr=0.9 "dh0=$dh0" || fail=1

if [[ $fail -ne 0 ]]; then
  echo "Manifold: FAILED"
  exit 1
fi
echo "Manifold: all checks passed"
exit 0
