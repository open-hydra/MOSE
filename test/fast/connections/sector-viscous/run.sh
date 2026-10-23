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
#    diffusive-metric = stencil; a velocity periodic for
#    a turn by 90 degrees with normal stress and shear on the meridian planes
#    (radial field, potential vortex, axial harmonic field). In every cell,
#    momentum and energy, the viscous forces must agree within 1e-9 of their
#    scale. (With mean the faces of type 200 take the inverse of the mean of the
#    inverse metrics of the cell and its image, the interior faces the mean of
#    the metrics: the two sectors then differ by 2.2e-2 of the scale, so mean
#    is not checked here.)
#  - mirror: the radial field alone, k faces 200 against 300 (the field is its
#    own mirror image across the planes); away from the j faces the residuals
#    must agree within 0.02 of the normal stress (2/3) mu U A/V.
#  - shift: a sector with no angle, a block periodic in z closed in k by type
#    200 (its last k node plane is the first one translated) against the same
#    block closed by periodic pairs (201); a shear flow u(z), Navier-Stokes with
#    Spalart-Allmaras. Both put the cell at the other end behind each k face,
#    at the same place (that cell translated by the period), so in every cell
#    the step of rho, momentum, energy and rho nu~ must agree within 1e-9 of
#    its largest value.
#  - stretch: the mirror sector with stretched radii (as the ring, and growing
#    by 1.2), k faces 200, the default metric (stencil), the radial field, a
#    linear field that the rotation allows: away from the j faces the viscous
#    force must vanish below 1e-11 of (2/3) mu U A/V, next to the periodic
#    faces too.
#  - wall: a translated sector with tilted k faces next to a wall 301, cells
#    growing from the wall and bent lines, so that the ghost cells of the wall
#    hold the velocity and the temperature at two positions: linear velocity
#    zero on the wall (wallmom) and linear temperature with its heat flux on
#    the wall (wallene); a Navier-Stokes step against one with mu = k = 1e-30.
#    In every cell with no face 400, the wall cells included, the viscous force
#    and the heat conduction must vanish below 1e-11 of their scales.
#  - split: type 200 joins the two k faces of one block (make_case.py), one
#    Euler step. A sector of 45 degrees cut along j, each block closed by
#    type 200, gives the step of the whole sector within 1e-12; cut along k,
#    with type 200 on face 4, with a record outside its block or with no
#    record on the other k face, and a 2-D wedge with type 200 on face 4 or
#    on one k face only, stop at set-up, before the mesh type and any step,
#    with exit code 1 and the message that names the record. A ring of four
#    blocks joined by 101, the junction at 2 pi included (the remedy of the
#    message), gives in block 1 the step of the sector closed by type 200:
#    Euler within 1e-12, viscous force within 1e-9.
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

# write_ini <equations> [diffusive-metric] [turbulence]: one explicit Euler step with local time stepping, first
# order; a turbulence model with explicit sources
write_ini() {
  local turb=""
  if [[ -n ${3:-} ]]; then
    turb="turbulence = $3

[MOSE-Turbulence]
point-implicit = false"
  fi
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = $1
schmidt = 0.7
$turb

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

# step <equations> <tag> [diffusive-metric] [turbulence]: run it, keep OUTPUT/diagnostic.tec as OUTPUT/diagnostic-<tag>.tec
step() {
  write_ini "$1" "${3:-}" "${4:-}"
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
for metric in stencil; do
  sfx="-stencil"
  for m in 45 90; do
    python3 -B make_case.py case "ring$m" || exit 1
    step navier-stokes "ring$m-ns$sfx" $metric
    step euler "ring$m-eu$sfx" $metric
  done
  if python3 -B make_case.py ring OUTPUT/diagnostic-ring45-ns$sfx.tec OUTPUT/diagnostic-ring45-eu$sfx.tec \
                                  OUTPUT/diagnostic-ring90-ns$sfx.tec OUTPUT/diagnostic-ring90-eu$sfx.tec; then
    echo "  [ok  ] $metric: the periodic faces carry the viscous flux of the interior faces of the annulus"
  else
    echo "  [FAIL] $metric: the viscous force next to the periodic faces is not that of the whole annulus"
    status=1
  fi
done

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

echo "-- shift: a sector with no angle, k faces 200 against periodic pairs 201 (Spalart-Allmaras) --"
for t in 200 201; do
  python3 -B make_case.py case "shift$t" || exit 1
  step navier-stokes "shift$t" "" SA
done
if python3 -B make_case.py shift OUTPUT/diagnostic-shift200.tec OUTPUT/diagnostic-shift201.tec; then
  echo "  [ok  ] the k faces of a translated sector give the step of the periodic pairs"
else
  echo "  [FAIL] the k faces of a translated sector and the periodic pairs give different steps"
  status=1
fi

echo "-- stretch: the radial field on stretched radii, k faces 200, the default metric --"
for law in pot geo; do
  python3 -B make_case.py case "stretch$law" || exit 1
  step navier-stokes "stretch$law-ns"
  step euler "stretch$law-eu"
  if python3 -B make_case.py linear "stretch$law" OUTPUT/diagnostic-stretch$law-ns.tec OUTPUT/diagnostic-stretch$law-eu.tec; then
    echo "  [ok  ] stretch$law: no viscous force of the linear field, next to the periodic faces too"
  else
    echo "  [FAIL] stretch$law: a viscous force of the linear field"
    status=1
  fi
done

# wall_step <run> <tag>: a Navier-Stokes step with INPUT/bc-<tag>.txt; for the reference (ref) mu = k = 1e-30
wall_step() {
  write_ini navier-stokes
  cp "INPUT/bc-$2.txt" INPUT/bc.txt
  rm -f INPUT/transport.dat
  if [[ $2 == ref ]]; then
    ln -sf transport-ref.dat INPUT/transport.dat
  else
    ln -sf "../$VISC/transport.dat" INPUT/transport.dat
  fi
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$1-$2.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$1-$2.log"; then
    echo "  [FAIL] the $2 step of $1 did not complete (see OUTPUT/$1-$2.log)"
    echo "SectorViscous: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$1-$2.tec"
}

echo "-- wall: a translated sector next to a wall, one metric per group of quantities at the periodic faces --"
for run in wallmom wallene; do
  python3 -B make_case.py case $run || exit 1
  wall_step $run ns
  wall_step $run ref
  if python3 -B make_case.py wall $run OUTPUT/diagnostic-$run-ns.tec OUTPUT/diagnostic-$run-ref.tec; then
    echo "  [ok  ] $run: the linear field is exact in every cell, next to the wall and the periodic faces too"
  else
    echo "  [FAIL] $run: the linear field is not exact next to the wall or the periodic faces"
    status=1
  fi
done
rm -f INPUT/transport.dat
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

echo "-- split: type 200 joins the two k faces of one block; other layouts stop at set-up --"
for c in whole cutj; do
  python3 -B make_case.py case "split-$c" || exit 1
  step euler "split-$c"
done
if python3 -B make_case.py split OUTPUT/diagnostic-split-whole.tec OUTPUT/diagnostic-split-cutj.tec; then
  echo "  [ok  ] a sector cut along j, each block closed by type 200, gives the step of the whole sector"
else
  echo "  [FAIL] the sector cut along j and the whole sector give different steps"
  status=1
fi
# stop_case <case> <message>: the run must stop at set-up, before the type of the mesh and any step: exit code 1, the
# line "Boundary conditions" with FAIL, the message, the STOP line of the set-up, no line "Mesh type", no diagnostic output
stop_case() {
  python3 -B make_case.py case "$1" || exit 1
  write_ini euler
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$1.log" 2>&1
  local rc=$?
  if [[ $rc -eq 1 ]] && grep -a -qE '^ +Boundary conditions +FAIL' "OUTPUT/$1.log" && grep -a -qF -- "$2" "OUTPUT/$1.log" \
     && grep -a -qF 'STOP: boundary condition error (see above)' "OUTPUT/$1.log" \
     && ! grep -a -q 'Mesh type' "OUTPUT/$1.log" && [[ ! -f OUTPUT/diagnostic.tec ]]; then
    echo "  [ok  ] $1 stops at set-up: $2"
  else
    echo "  [FAIL] $1 (exit code $rc) does not stop at set-up with: $2 (see OUTPUT/$1.log)"
    status=1
  fi
}
rule="type 200 joins the k faces of one block:"
ring="$rule cut sectors along i or j; join rings by 101 at each junction"
k_only="type 200 is for the k faces (5 and 6) of a wedge or a sector only"
stop_case split-cutk "[ERROR] bc.txt, record 49: type 200 on face 5 of block 1 at i, j = 1 1, but face 6 there is of type 101: $ring"
stop_case split-face4 "[ERROR] bc.txt, record 81: type 200 on face 4 of block 1 at i, j, k = 1 4 1: $k_only"
stop_case split-outside "[ERROR] bc.txt, record 97: type 200 on face 5 of block 1 at i, j, k = 3 1 1: the cell is outside the block, of 2 x 4 x 8 cells"
stop_case split-missing "[ERROR] bc.txt, record 105: type 200 on face 6 of block 1 at i, j = 1 1, but face 5 there has no record: $rule give face 5 there type 200 too"
stop_case split-wedge4 "[ERROR] bc.txt, record 13: type 200 on face 4 of block 1 at i, j, k = 1 4 1: $k_only"
stop_case split-wedge6 "[ERROR] bc.txt, record 17: type 200 on face 5 of block 1 at i, j = 1 1, but face 6 there is of type 300: $rule give face 6 there type 200 too"
python3 -B make_case.py case ring101 || exit 1
step navier-stokes ring101-ns-stencil stencil
step euler ring101-eu-stencil stencil
if python3 -B make_case.py ring101 OUTPUT/diagnostic-ring45-ns-stencil.tec OUTPUT/diagnostic-ring45-eu-stencil.tec \
                                   OUTPUT/diagnostic-ring101-ns-stencil.tec OUTPUT/diagnostic-ring101-eu-stencil.tec; then
  echo "  [ok  ] a ring of four blocks joined by 101, the junction at 2 pi included, gives the step of the sector closed by 200"
else
  echo "  [FAIL] the ring of four blocks joined by 101 and the sector closed by type 200 give different steps"
  status=1
fi

if [[ $status -eq 0 ]]; then
  echo "SectorViscous: all checks passed"
  exit 0
fi
echo "SectorViscous: FAILED"
exit 1
