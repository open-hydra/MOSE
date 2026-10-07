#!/usr/bin/env bash
#===============================================================================
#  Sector radial: rotational periodicity of a 3-D sector (fast)
#
#  A 3-D sector of an annulus closed in k by type 200: the cells next to one k
#  face go on in those next to the other, turned by the sector angle. The
#  sector spans -45 to 45 degrees and the flow is radial, with uniform density
#  and pressure: its own mirror image across z = 0 (make_case.py). After 20
#  steps of RK3 the columns k and NK + 1 - k must hold the same rho, u, v and
#  p, and opposite w, within 1e-10 of their scales.
#  Gas: test/common/Air.
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

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$AIR/phase.txt"  INPUT/phase.txt
ln -sf "../$AIR/thermo.dat" INPUT/thermo.dat

cat > input.ini <<EOT
[MOSE-Parameters]
iter-threshold = 20
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOT

echo "-- radial flow in a 3-D sector closed by the rotational periodicity --"
python3 -B make_case.py case || exit 1
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/sector.log 2>&1
if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed' OUTPUT/sector.log; then
  echo "  [FAIL] the run did not complete (see OUTPUT/sector.log)"
  echo "SectorRadial: FAILED"
  exit 1
fi
read -r dev where < <(python3 -B make_case.py check OUTPUT/field.tec)
if python3 -c "import sys; sys.exit(0 if float('$dev') <= 1e-10 else 1)"; then
  echo "  [ok  ] largest departure from the mirror symmetry $dev ($where)"
  echo "SectorRadial: all checks passed"
  exit 0
fi
echo "  [FAIL] largest departure from the mirror symmetry $dev ($where)"
echo "SectorRadial: FAILED"
exit 1
