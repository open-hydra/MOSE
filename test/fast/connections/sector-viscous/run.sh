#!/usr/bin/env bash
#===============================================================================
#  Viscous flux of the rotational periodicity of a 3-D sector (fast)
#
#  A 3-D sector of an annulus closed in k by type 200: its two k faces are
#  interior faces of the whole annulus, and their viscous flux must be that of
#  an interior face (make_case.py). Viscous air at 10 Pa (the pressure only
#  sets the round-off of the viscous force, make_case.py), one explicit Euler
#  step; the viscous force of a cell is r/dt of a Navier-Stokes step minus
#  r/dt of an Euler step from the same state (OUTPUT/diagnostic.tec).
#  - ring: the sector from -45 to 45 degrees against the sector from -90 to 90
#    degrees that holds it, where the faces at +-45 degrees are interior, with
#    diffusive-metric = inverse-mean; a velocity periodic for a turn by 90
#    degrees with normal stress and shear on the meridian planes (radial
#    field, potential vortex, axial harmonic field). In every cell, momentum
#    and energy, the viscous forces must agree within 1e-9 of their scale.
#  - mirror: the radial field alone, k faces 200 against 300 (the field is its
#    own mirror image across the planes); away from the j faces the residuals
#    must agree within 0.02 of the normal stress (2/3) mu U A/V.
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
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

# write_ini <equations> [diffusive-metric]: one explicit Euler step with local time stepping, first order
write_ini() {
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = $1
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC
EOF
  if [[ -n ${2:-} ]]; then echo "diffusive-metric = $2" >> input.ini; fi
}

# step <equations> <tag> [diffusive-metric]: run it, keep OUTPUT/diagnostic.tec as OUTPUT/diagnostic-<tag>.tec
step() {
  write_ini "$1" "${3:-}"
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$2.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$2.log"; then
    echo "  [FAIL] the $1 step of $2 did not complete (see OUTPUT/$2.log)"
    echo "SectorViscous: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$2.tec"
}

status=0

echo "-- ring: the sector of 90 degrees against the same cells inside the sector of 180 degrees --"
for m in 45 90; do
  python3 -B make_case.py case "ring$m" || exit 1
  step navier-stokes "ring$m-ns" inverse-mean
  step euler "ring$m-eu" inverse-mean
done
if python3 -B make_case.py ring OUTPUT/diagnostic-ring45-ns.tec OUTPUT/diagnostic-ring45-eu.tec \
                                OUTPUT/diagnostic-ring90-ns.tec OUTPUT/diagnostic-ring90-eu.tec; then
  echo "  [ok  ] the periodic faces carry the viscous flux of the interior faces of the annulus"
else
  echo "  [FAIL] the viscous force next to the periodic faces is not that of the whole annulus"
  status=1
fi

echo "-- mirror: radial flow, k faces 200 against symmetry planes 300 --"
for t in 200 300; do
  python3 -B make_case.py case "mirror$t" || exit 1
  step navier-stokes "mirror$t"
done
if python3 -B make_case.py mirror OUTPUT/diagnostic-mirror200.tec OUTPUT/diagnostic-mirror300.tec; then
  echo "  [ok  ] the periodic faces carry the normal viscous stress of the symmetry planes"
else
  echo "  [FAIL] the periodic faces and the symmetry planes give different viscous forces"
  status=1
fi

if [[ $status -eq 0 ]]; then
  echo "SectorViscous: all checks passed"
  exit 0
fi
echo "SectorViscous: FAILED"
exit 1
