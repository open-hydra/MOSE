#!/usr/bin/env bash
#===============================================================================
#  Choked nozzle (420) below p_sup: the supersonic design exit state (fast)
#
#  Below p_sup the face of a nozzle 420 takes the supersonic design exit state
#  of its record, the isentropic state of the reservoir with the mass flux g,
#  on the tables of the gas (docs/user/boundary-conditions.md). The records are
#  those that the boundary-condition builder of ATLAS writes for real nozzles,
#  from their exit-to-throat area ratio (design_state.py). A 1-D tube with the
#  nozzle on face 1, its gas the design state at 0.98 of its pressure, moving
#  in; one forward Euler step with first-order reconstruction: the fluxes of
#  mass, momentum and energy through the nozzle face, from the change of its
#  cell, must be those of the design state, g, p_e + g u_e and g h0, within
#  1e-9 of each:
#    mach25, mach3, mach4  air, T0 = 300 K, p0 = 2 bar, nozzles of exit Mach
#              number 2.5, 3 and 4;
#    rocket45  the gas of the inlet of the frozen rocket nozzle of
#              test/2D/euler/rocket-nozzle (13 species, T0 = 3644.12 K,
#              p0 = 108 bar), a nozzle of area ratio 40, exit Mach number 4.46;
#    sonic     the air of mach25 at the throat (area ratio 1) with g 0.5e-9
#              above the choked mass flux: the sonic state, whose momentum
#              flux alone is checked (design_state.py).
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
for name in mach25 mach3 mach4 rocket45 sonic; do
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
    echo "  [ok  ] $name: the fluxes through the nozzle face are those of the design state"
  else
    echo "  [FAIL] $name: the fluxes through the nozzle face are not those of the design state"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "Nozzle: FAILED"
  exit 1
fi
echo "Nozzle: all checks passed"
exit 0
