#!/usr/bin/env bash
#===============================================================================
#  Choked nozzle (420) below and above p_shock: design state, mass flux g (fast)
#
#  The pressure p of the boundary cell selects the state of the face of a
#  nozzle 420 (docs/user/boundary-conditions.md): below p_shock, the pressure
#  behind a normal shock standing at the exit of the design flow, the
#  supersonic design exit state of its record, the isentropic state of the
#  reservoir with the mass flux g, on the tables of the gas; from p_shock up
#  to p_sub the state of mass flux g and total temperature T0 that the
#  invariant of the cell reaches. The records are those that the
#  boundary-condition builder of ATLAS writes for real nozzles, from their
#  exit-to-throat area ratio (design_state.py). A 1-D tube with the nozzle on
#  face 1, its gas moving in; one forward Euler step with first-order
#  reconstruction: the fluxes of mass, momentum and energy through the nozzle
#  face, from the change of its cell, must be those of the state the face
#  takes, within 1e-9 of each. With the gas of the tube in the design state
#  at 0.98 of its pressure, the face takes the design state, g, p_e + g u_e
#  and g h0:
#    mach25, mach3, mach4  air, T0 = 300 K, p0 = 2 bar, nozzles of exit Mach
#              number 2.5, 3 and 4;
#    rocket45  the gas of the inlet of the frozen rocket nozzle of
#              test/2D/euler/rocket-nozzle (13 species, T0 = 3644.12 K,
#              p0 = 108 bar), a nozzle of area ratio 40, exit Mach number 4.46;
#    sonic     the air of mach25 at the throat (area ratio 1) with g 0.5e-9
#              above the choked mass flux: the sonic state, whose momentum
#              flux alone is checked (design_state.py).
#  The air of mach25 in a nozzle of area ratio 4 (exit Mach number 2.94),
#  the gas of the tube on the line of the states of mass flux g and total
#  enthalpy h0, at a pressure between p_sup and p_sub:
#    overexp-sup  2 p_sup, supersonic: the design state;
#    overexp-sub  midway between p_son, the pressure of the sonic state of the
#              line, and p_shock, subsonic: the design state;
#    shock-below, shock-above  1e-6 below and above p_shock: the design state,
#              then the state of the gas of the tube, which holds the p_shock
#              of the solver within 1e-6 of the pressure behind a normal shock
#              at the Mach number of the design state;
#    shock-in  1.02 p_shock, a normal shock in the divergent part: the state of
#              the gas of the tube.
#  Gas: test/common/Air, and the gas of test/2D/euler/rocket-nozzle/frozen.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../../common/Air
ROCKET=../../../2D/euler/rocket-nozzle/frozen/INPUT

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC
INI

fail=0
for name in mach25 mach3 mach4 rocket45 sonic overexp-sup overexp-sub shock-below shock-above shock-in; do
  echo "   $name"
  gas=$AIR
  [[ $name == rocket45 ]] && gas=$ROCKET
  ln -sf "../$gas/phase.txt"  INPUT/phase.txt
  ln -sf "../$gas/thermo.dat" INPUT/thermo.dat
  python3 -B design_state.py case "$name" || exit 1
  rm -f OUTPUT/field.tec "OUTPUT/$name.tec"
  timeout 60 "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] the run did not complete (see OUTPUT/$name.log)"
    fail=1
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/$name.tec"
  if python3 -B design_state.py check "$name" "OUTPUT/$name.tec"; then
    echo "  [ok  ] $name: the fluxes through the nozzle face are those of the state the face takes"
  else
    echo "  [FAIL] $name: the fluxes through the nozzle face are not those of the state the face takes"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "Nozzle: FAILED"
  exit 1
fi
echo "Nozzle: all checks passed"
exit 0
