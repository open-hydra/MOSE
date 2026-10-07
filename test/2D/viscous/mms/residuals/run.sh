#!/usr/bin/env bash
#===============================================================================
#  Right-hand sides on constructed fields: Spalart-Allmaras and the viscous
#  stress of the axisymmetric wedge (fast tier)
#
#  1. Variable density: a constructed periodic field at rest with variable
#     temperature and nu~ and no wall (make_case.py): one explicit Euler step
#     with local time stepping gives, in OUTPUT/diagnostic.tec, the right-hand
#     side of the equation for rho nu~ that the solver evaluates. It must
#     converge to the exact right-hand side of the standard conservation form
#     (NASA Turbulence Modeling Resource, Spalart-Allmaras page) with second
#     order, on 32 x 32 and 64 x 64 cells.
#  2. Rotation tensor of SA-RC and SST-RC: W_ij = (du_i/dx_j - du_j/dx_i)/2
#     (same page) is a tensor, so the right-hand side of rho nu~, and of rho k
#     and rho omega with SST-RC, must not depend on how the flow is turned in
#     space. A cube of 8^3 cells, one block (-RC needs one), with a 3-D
#     velocity whose vorticity has three components and whose strain varies
#     in space, so that the term r~ of the rotation-curvature correction is
#     not zero, and the same field turned by the rotation that takes x to y,
#     y to z and z to x; one explicit step each, with SA-RC and with SST-RC
#     (k = 1 m2/s2, omega = 3000 1/s). The change of rho nu~ (rho k and rho
#     omega) must be the same in corresponding cells within 1e-9 of its
#     largest value. A rotation tensor whose W_23 has the sign opposite to
#     that of W_12 and W_13 (+omega_x/2 instead of -omega_x/2) makes the two
#     changes differ by their own size.
#  3. SA-RC in a rotating frame and where the vorticity vanishes: f_r1
#     multiplies the production of nu~, and its r~ takes, in a rotating frame,
#     the rotation of the frame next to the material derivative of the strain
#     (same page, SA-RC). A planar block of 8 x 8 cubic cells, one explicit
#     step: a uniform shear in a frame turning about z, where the change of
#     rho nu~ must be dt rho c_b1 f_r1 S~ nu~ in every cell within 1e-6 of it,
#     with r~ from the rotation of the frame alone (f_r1 = 4.65; 1.5 without
#     the rotation term); and a uniform flow in a fixed frame, where the
#     vorticity is zero and r* = S/omega divides by zero: the run must
#     complete and rho nu~ must not change.
#     And air at rest next to a wall (an isothermal wall on one face, rho nu~
#     = 0.5 mu, so that nu~ is produced next to the wall), with SA and with
#     SA-RC: strain and vorticity vanish, f_r1 takes the value 1, and SA-RC
#     must change rho nu~ as SA does, within 1e-12 of the largest change.
#     And the same air flowing at 10 m/s along the wall, on cells whose widths
#     along the flow grow by 1.15 from one to the next: off the wall row the
#     velocity differences are zero, so that strain, vorticity and r~ come from
#     the 1e-40 of the velocity gradient alone; f_r1 must take the value 1
#     there too (r~ = 0), and SA-RC must change rho nu~ in those rows as SA
#     does, within 1e-12 of the largest change.
#  4. Viscous stress on the 2-D axisymmetric wedge (k faces of type 200,
#     1 degree, axis at j = 0; wedge_viscous.py): Homann's axisymmetric
#     stagnation-point flow, an exact solution of the Navier-Stokes
#     equations, divergence free, and a polynomial field with dilatation
#     and swirl, u_x = b x^2, u_r = c r^3, u_theta = s r^3 (1 + x/L), with a
#     closed-form viscous force. One explicit Euler step with local time
#     stepping with the Navier-Stokes equations and one with the Euler
#     equations, from the same state: the difference of the right-hand
#     sides in OUTPUT/diagnostic.tec is the viscous force the solver
#     evaluates. On 32 x 32 and 64 x 64 cells it must converge to the exact
#     one with second order: homann, the radial (on the cells with r >= L/2;
#     next to the axis it must match the hoop stress 2 mu u_r/r^2) and axial
#     components, the azimuthal one staying zero; swirl, the radial and
#     azimuthal components on the cells with r >= L/2, the axial one exact.
#     Axis at r = 0 on the wedge: the faces of the wedge would get the hoop
#     stress of the first row wrong, so a viscous run must stop at setup and
#     say why, with diffusive-metric = stencil and inverse-mean. Polar line of
#     a 3-D mesh (360 degrees around the axis, axis at r = 0): with symmetry on
#     the axis the viscous force of a linear velocity field must vanish to
#     round-off, next to the axis too; with a wall on the axis the run must
#     stop at setup with the message that names the face and the remedy
#     (symmetry on the axis), without NaN before it.
#  5. Spalart-Allmaras next to a wall (wall_gradient.py): a column of fluid
#     at rest over an isothermal wall with nu~ = a y, chi below 1e-3: the
#     linear profile is an equilibrium of the Spalart-Allmaras equation
#     (cw1 = cb1/kappa^2 + (1 + cb2)/sigma). One explicit Euler step with
#     local time stepping; the change of rho nu~ over the step is the
#     right-hand side, which must vanish in every cell next to the wall
#     too: there the source terms must see grad nu~ = a, from the wall
#     value nu~ = 0 (ghost cell -nu~), not a/2 from a ghost cell that
#     copies the cell.
#     The same with cells growing by 1.2 from the wall, and on skewed cells
#     (2-D and 3-D) with the density linear along the wall, where the ghost
#     cell of the wall copies the density and the terms in grad rho must stay
#     zero. SA-RC: the uniform column runs; the same column as two blocks must
#     stop at setup and say why (its result would depend on the division in
#     blocks).
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../../common/Air-viscous

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

cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = SA
schmidt = 0.7

[MOSE-Turbulence]
point-implicit = false

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF

fail=0
errors=()
for n in 32 64; do
  python3 -B make_case.py case "$n" || exit 1
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/n$n.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/n$n.log"; then
    echo "  [FAIL] the run on $n x $n cells did not complete (see OUTPUT/n$n.log)"
    echo "MMSResiduals: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-n$n.tec"
  errors+=("$(python3 -B make_case.py error "OUTPUT/diagnostic-n$n.tec")")
done

echo "-- 1. right-hand side of rho nu~ against the standard conservation form --"
if python3 -B make_case.py check "${errors[0]}" "${errors[1]}"; then
  echo "  [ok  ] second-order agreement"
else
  echo "  [FAIL] the solver does not evaluate the standard form"
  fail=1
fi

for model in SA-RC SST-RC; do
  echo "-- 2. $model: a 3-D flow and the same flow turned, one step --"
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = $model

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
EOF
  for name in a b; do
    python3 -B make_case.py rotation "$name" "$model" || exit 1
    rm -f OUTPUT/field.tec "OUTPUT/field-$model-$name.tec"
    "$LOCAL" > "OUTPUT/rotation-$model-$name.log" 2>&1
    if [[ -f OUTPUT/field.tec ]]; then
      cp OUTPUT/field.tec "OUTPUT/field-$model-$name.tec"
    else
      echo "  (rotation $model $name: no OUTPUT/field.tec, see OUTPUT/rotation-$model-$name.log)"
    fi
  done
  if python3 -B make_case.py rotation-check "OUTPUT/field-$model-a.tec" "OUTPUT/field-$model-b.tec"; then
    echo "  [ok  ] $model: the change of the turbulence variables does not depend on how the flow is turned"
  else
    echo "  [FAIL] $model: the change of the turbulence variables depends on how the flow is turned"
    fail=1
  fi
done

echo "-- 3. SA-RC in a rotating frame and where the vorticity vanishes, one step --"
# write_rc_ini <rotational-frame>
write_rc_ini() {
  {
    printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 1' 'newrun = true' 'res-threshold = 1e-300' ''
    printf '%s\n' '[MOSE-Physics]' 'equations = navier-stokes' 'turbulence = SA-RC' "rotational-frame = $1" ''
    printf '%s\n' '[MOSE-Numerics]' 'cfl = 0.1' 'time-accurate = true' 'time-scheme = euler' \
                  'integration-variables = cons' 'space-reconstruction = first-order' 'riemann-solver = HLLC' ''
    printf '%s\n' '[MOSE-IO]' 'sol-format = tecplot ascii'
    if [[ $1 == rigid-body ]]; then
      printf '\n%s\n' '[MOSE-Rotating-Frame]' 'omega = 200.0' 'axis = 0.0 0.0 1.0' 'origin = 0.0 0.0 0.0'
    fi
  } > input.ini
}
for kind in shear uniform; do
  if [[ $kind == shear ]]; then write_rc_ini rigid-body; else write_rc_ini none; fi
  python3 -B make_case.py rc "$kind" || exit 1
  rm -f OUTPUT/field.tec "OUTPUT/field-$kind.tec"
  "$LOCAL" > "OUTPUT/rc-$kind.log" 2>&1
  if [[ -f OUTPUT/field.tec ]]; then
    cp OUTPUT/field.tec "OUTPUT/field-$kind.tec"
  else
    echo "  (rc $kind: no OUTPUT/field.tec, see OUTPUT/rc-$kind.log)"
  fi
  if python3 -B make_case.py rc-check "$kind" "OUTPUT/field-$kind.tec"; then
    echo "  [ok  ] $kind: the change of rho nu~ is the one of f_r1"
  else
    echo "  [FAIL] $kind: the change of rho nu~ is not the one of f_r1"
    fail=1
  fi
done

# write_rest_ini <turbulence model>: the cases next to a wall, fixed frame, one step
write_rest_ini() {
  {
    printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 1' 'newrun = true' 'res-threshold = 1e-300' ''
    printf '%s\n' '[MOSE-Physics]' 'equations = navier-stokes' "turbulence = $1" 'rotational-frame = none' ''
    printf '%s\n' '[MOSE-Numerics]' 'cfl = 0.1' 'time-accurate = true' 'time-scheme = euler' \
                  'integration-variables = cons' 'space-reconstruction = first-order' 'riemann-solver = HLLC' ''
    printf '%s\n' '[MOSE-IO]' 'sol-format = tecplot ascii'
  } > input.ini
}
for kind in rest stream; do
  python3 -B make_case.py "$kind" || exit 1
  for model in SA SA-RC; do
    write_rest_ini "$model"
    rm -f OUTPUT/field.tec "OUTPUT/field-$kind-$model.tec"
    "$LOCAL" > "OUTPUT/$kind-$model.log" 2>&1
    if [[ -f OUTPUT/field.tec ]]; then
      cp OUTPUT/field.tec "OUTPUT/field-$kind-$model.tec"
    else
      echo "  ($kind $model: no OUTPUT/field.tec, see OUTPUT/$kind-$model.log)"
    fi
  done
  where="next to the wall"
  [[ $kind == stream ]] && where="off the wall row of a uniform flow on stretched cells"
  if python3 -B make_case.py "$kind-check" "OUTPUT/field-$kind-SA.tec" "OUTPUT/field-$kind-SA-RC.tec"; then
    echo "  [ok  ] $kind: SA-RC changes rho nu~ $where as SA does"
  else
    echo "  [FAIL] $kind: SA-RC does not change rho nu~ $where as SA does"
    fail=1
  fi
done

# wedge_ini <equations> [diffusive-metric]: one explicit Euler step with local time stepping
wedge_ini() {
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
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
  if [[ -n ${2:-} ]]; then echo "diffusive-metric = $2" >> input.ini; fi
}

# wedge_step <equations> <tag>: run it, keep OUTPUT/diagnostic.tec as OUTPUT/diagnostic-<tag>.tec
wedge_step() {
  wedge_ini "$1"
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$2.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$2.log"; then
    echo "  [FAIL] the $1 step did not complete (see OUTPUT/$2.log)"
    echo "MMSResiduals: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$2.tec"
}

echo "-- 4. viscous stress on the 2-D axisymmetric wedge, 32 x 32 and 64 x 64 cells --"
for field in homann swirl; do
  errors=()
  for n in 32 64; do
    python3 -B wedge_viscous.py case "$field" "$n" || exit 1
    wedge_step navier-stokes "$field-ns-n$n"
    wedge_step euler "$field-eu-n$n"
    errors+=("$(python3 -B wedge_viscous.py error "$field" "OUTPUT/diagnostic-$field-ns-n$n.tec" \
                                                    "OUTPUT/diagnostic-$field-eu-n$n.tec")")
  done

  if [[ $field == homann ]]; then
    echo "-- viscous force on Homann's stagnation-point flow, 2-D axisymmetric wedge --"
  else
    echo "-- viscous force on a polynomial field with dilatation and swirl, 2-D axisymmetric wedge --"
  fi
  if python3 -B wedge_viscous.py check "$field" "${errors[0]}" "${errors[1]}"; then
    echo "  [ok  ] second-order agreement"
  else
    echo "  [FAIL] the viscous force on the wedge is not the Navier-Stokes one"
    fail=1
  fi
done

echo "-- the axis at r = 0 on the wedge (homann, N = 32): a declared stop of the viscous runs --"
python3 -B wedge_viscous.py case homann 32 axis0 || exit 1
for metric in stencil inverse-mean; do
  wedge_ini navier-stokes "$metric"
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/axis0-$metric.log" 2>&1
  if python3 -B wedge_viscous.py axisstop "OUTPUT/axis0-$metric.log"; then
    echo "  [ok  ] $metric: the run stops at setup and says why"
  else
    echo "  [FAIL] $metric: no stop with the message on the axis at r = 0 (see OUTPUT/axis0-$metric.log)"
    fail=1
  fi
done

echo "-- polar line of a 3-D mesh, symmetry on the axis: a linear velocity field --"
python3 -B wedge_viscous.py polar 300 || exit 1
wedge_step navier-stokes polar-300-ns
wedge_step euler polar-300-eu
if python3 -B wedge_viscous.py linear OUTPUT/diagnostic-polar-300-ns.tec OUTPUT/diagnostic-polar-300-eu.tec; then
  echo "  [ok  ] no viscous force, next to the axis too"
else
  echo "  [FAIL] the viscous force of a linear field does not vanish on the polar mesh"
  fail=1
fi

echo "-- polar line of a 3-D mesh, wall on the axis: a declared limit of the stencil metric --"
python3 -B wedge_viscous.py polar 302 || exit 1
wedge_ini navier-stokes
rm -f OUTPUT/diagnostic.tec
"$LOCAL" > OUTPUT/polar-302.log 2>&1
if python3 -B wedge_viscous.py stopped OUTPUT/polar-302.log; then
  echo "  [ok  ] the run stops at setup, names the face and the remedy"
else
  echo "  [FAIL] the run does not stop with the message on the wall faces of the axis (see OUTPUT/polar-302.log)"
  fail=1
fi

cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = SA
schmidt = 0.7

[MOSE-Turbulence]
point-implicit = false

[MOSE-Numerics]
cfl = 0.5
vnn = 0.1
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF

# solve <case> [uni]: one explicit step, outputs kept as OUTPUT/wall-{field,diagnostic}-<case>[-uni].tec
solve() {
  local tag="$1${2:+-$2}"
  python3 -B wall_gradient.py case "$@" || exit 1
  rm -f OUTPUT/field.tec OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/wall-run-$tag.log" 2>&1
  if [[ ! -f OUTPUT/field.tec || ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/wall-run-$tag.log"; then
    echo "  [FAIL] the run $tag did not complete (see OUTPUT/wall-run-$tag.log)"
    echo "MMSResiduals: FAILED"
    exit 1
  fi
  mv OUTPUT/field.tec "OUTPUT/wall-field-$tag.tec"
  mv OUTPUT/diagnostic.tec "OUTPUT/wall-diagnostic-$tag.tec"
}

echo "-- 5. Spalart-Allmaras next to a wall: right-hand side of rho nu~ for a linear profile --"
solve uniform
if python3 -B wall_gradient.py check OUTPUT/wall-field-uniform.tec OUTPUT/wall-diagnostic-uniform.tec; then
  echo "  [ok  ] the source terms see the wall value nu~ = 0 in the wall cells"
else
  echo "  [FAIL] the source terms of the wall cells see a wrong gradient of nu~"; fail=1
fi

echo "-- the same on cells growing by 1.2 from the wall --"
solve stretched
if python3 -B wall_gradient.py check OUTPUT/wall-field-stretched.tec OUTPUT/wall-diagnostic-stretched.tec stretched; then
  echo "  [ok  ] the cell gradients take the differences of the cell centres"
else
  echo "  [FAIL] the cell gradients of the stretched cells are wrong"; fail=1
fi

echo "-- skewed cells, rho linear along the wall: the terms in grad rho (2-D and 3-D) --"
for case in skewed skewed3d; do
  solve "$case"
  solve "$case" uni
  if python3 -B wall_gradient.py checkdiff "$case" "OUTPUT/wall-field-$case.tec" "OUTPUT/wall-diagnostic-$case.tec" \
       "OUTPUT/wall-field-$case-uni.tec" "OUTPUT/wall-diagnostic-$case-uni.tec"; then
    echo "  [ok  ] $case: the copied density of the wall ghost cells sits at the cell itself"
  else
    echo "  [FAIL] $case: grad rho of the wall cells takes a wrong normal component"; fail=1
  fi
done

echo "-- SA-RC: one block runs, two blocks stop at setup --"
sed -i 's/^turbulence = SA$/turbulence = SA-RC/' input.ini
python3 -B wall_gradient.py case uniform || exit 1
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/wall-run-sarc-1.log 2>&1
if [[ -f OUTPUT/field.tec ]] && ! grep -a -qE 'NaN|Integration failed|ABORT' OUTPUT/wall-run-sarc-1.log; then
  echo "  [ok  ] one block: the run completes"
else
  echo "  [FAIL] SA-RC on one block did not complete (see OUTPUT/wall-run-sarc-1.log)"; fail=1
fi
python3 -B wall_gradient.py case uniform2 || exit 1
"$LOCAL" > OUTPUT/wall-run-sarc-2.log 2>&1
if python3 -B wall_gradient.py sarcstop OUTPUT/wall-run-sarc-2.log; then
  echo "  [ok  ] two blocks: the run stops at setup and says why"
else
  echo "  [FAIL] SA-RC on two blocks does not stop with the message (see OUTPUT/wall-run-sarc-2.log)"; fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "MMSResiduals: FAILED"
  exit 1
fi
echo "MMSResiduals: all checks passed"
exit 0
