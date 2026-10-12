#!/usr/bin/env bash
#===============================================================================
#  Exit codes of the error stops (fast tier)
#
#  A run that MOSE stops on an error must end with a non-zero exit code, so
#  that scripts and CI see the failure, and with its message, not with a
#  crash signal; it writes no solution. Reuses the inputs of the 100-cell Sod
#  tube (test/1D/Sod79) and the gases of test/common and test/1D/diff1; no
#  new mesh or IC data is stored here.
#
#  Unknown probe variable. A probe that asks for a variable MOSE does not
#  know must stop the run at set-up with exit status 1 and a message naming
#  the variable, not crash: the probe pointer stayed undefined and the run
#  ended on a segmentation fault (exit status 139). The Sod tube of
#  test/1D/Sod79 in air, one probe with variables "p foo".
#
#  Default Schmidt number. schmidt defaults to 0.0, i.e. mixture-averaged
#  multicomponent diffusion, which needs the binary diffusion table
#  INPUT/diffusion.dat when the gas has more than one species (one species
#  does not diffuse: AblatingWall runs a single-species column without
#  schmidt). The Sod tube as a viscous mixture of four species (those of
#  test/1D/diff1, with the mass fractions of air), without schmidt and
#  without a table, must stop at set-up with the diffusion-table error and
#  write no solution.
#
#  Diagnostic output of a single species without thermo or transport. The
#  diagnostic output of a steady run reuses the buffer of the solution
#  output, which with a single species and sol-variables naming neither
#  thermo nor transport is one variable short: the Sod tube in air, steady,
#  with sol-variables = none, must stop at its first output with exit status
#  1 and a message that says to name thermo or transport.
#
#  Diagnostic output with a soot model. A steady run with the soot model
#  LL91 must end normally, and its diagnostic output must name every
#  variable it writes, the residuals of rho, rho u, rho v, rho w, rho E, of
#  the soot number density and of the turbulence variable, then dt and beta,
#  with the local time step under dt (the one of the same run without soot,
#  within 1e-6). Six species with the names the soot model looks for, in a
#  uniform state on the mesh of the Sod tube (make_case.py soot), with type 0
#  on its faces 3 and 4, with and without Spalart-Allmaras, two iterations.
#  With Spalart-Allmaras the output filled the turbulence column from arrays
#  of different sizes (a run-time error with bounds checking); without it,
#  it named one variable fewer than it wrote and put the density of C(gr)
#  under dt.
#
#  turbulence = SAcomp. Its compressibility correction takes c from
#  c^2 f(c) = 25 Omega nu_t / a, a number against a speed, so its result
#  depends on the units: on the state of the soot runs with
#  Spalart-Allmaras and no soot, which runs to the end with SA, the run must
#  stop at set-up with exit code 1 to 128 and the message, and write no
#  diagnostic output.
#
#  Error stops. Each case below breaks the Sod tube in one way: a run that
#  diverges (CFL 20; stop inside the OpenMP loop of the state update), an
#  initial condition with a NaN velocity in one cell (the same stop: the
#  check of the state update must catch a value that is not finite in any
#  build; -ffast-math, an option of the RELEASE build, lets gfortran fold
#  isnan to false), an initial condition with fewer variables than the model
#  needs (set-up), a rotating frame with an axis of zero length (set-up).
#  The exit code must be between 1 and 128 (128 is what ifx gives to ERROR
#  STOP; above 128 the shell reports a signal) and the log must hold the
#  message.
#
#  And a valid run must end with exit code 0: with finite-rate chemistry and
#  exclude-blocks = 1 the same tube must run, its set-up summary listing
#  block 1 under "Inactive blocks" (the setup deallocated that list before it
#  was ever allocated: exit code 153 with ifx, "allocatable array or pointer
#  is not allocated"). The air has no reaction (no mechanism file): the
#  chemistry leaves the state as it is.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79/INPUT
MIX=../../../1D/diff1/INPUT
COMMON=../../../common/Air
VISC=../../../common/Air-viscous
ERR='Sc<=0 (multicomponent diffusion) requires a binary diffusion table'

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

fail=0

# check <status> <description>: status 0 = pass
check() {
  if [[ $1 -eq 0 ]]; then
    echo "  [ok  ] $2"
  else
    echo "  [FAIL] $2"
    fail=1
  fi
}

echo "-- probe with an unknown variable --"
rm -f INPUT/*
ln -sf "../../../../common/Air/phase.txt" INPUT/phase.txt
ln -sf "../../../../common/Air/thermo.dat" INPUT/thermo.dat
ln -sf "../$SOD/ic_x1.tec" INPUT/ic.tec
ln -sf "../$SOD/bc_x1.txt" INPUT/bc.txt
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2
newrun = true

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.9
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-Probes]
probe1 = pA

[pA]
variables = p foo
diter = 1
position = 0.5 0.0 0.0
EOI
rm -f OUTPUT/field.tec OUTPUT/pA.txt
"$LOCAL" > OUTPUT/probe-variable.log 2>&1
rc=$?
[[ $rc -ne 0 && $rc -lt 128 ]]
check $? "the solver stops with an error, not a crash (exit status $rc)"
grep -a -qF 'unknown probe variable "foo"' OUTPUT/probe-variable.log
check $? "it names the unknown variable"
[[ ! -f OUTPUT/field.tec ]]
check $? "no solution is written"

echo "-- default schmidt: four species, no diffusion table --"
rm -f INPUT/*
for f in phase.txt thermo.dat transport.dat; do
  ln -sf "../$MIX/$f" "INPUT/$f"
done
ln -sf "../$SOD/bc_x1.txt" INPUT/bc.txt
python3 -B make_case.py mixture "$SOD/ic_x1.tec" INPUT/ic.tec || exit 1
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2
newrun = true

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOI
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/schmidt-mixture.log 2>&1
rc=$?
[[ $rc -ne 0 ]]
check $? "the solver stops with an error (exit status $rc)"
grep -a -qF "$ERR" OUTPUT/schmidt-mixture.log
check $? "it reports the missing diffusion table"
[[ ! -f OUTPUT/field.tec ]]
check $? "no solution is written"

echo "-- single species, steady, sol-variables = none: the diagnostic output --"
rm -f INPUT/*
ln -sf "../../../../common/Air/phase.txt" INPUT/phase.txt
ln -sf "../../../../common/Air/thermo.dat" INPUT/thermo.dat
ln -sf "../$SOD/ic_x1.tec" INPUT/ic.tec
ln -sf "../$SOD/bc_x1.txt" INPUT/bc.txt
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2
newrun = true

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
sol-variables = none
EOI
"$LOCAL" > OUTPUT/diagnostic-single.log 2>&1
rc=$?
[[ $rc -ne 0 && $rc -lt 128 ]]
check $? "the solver stops with an error, not a crash (exit status $rc)"
grep -a -qF 'needs sol-variables to include thermo (the default) or transport' OUTPUT/diagnostic-single.log
check $? "it says which sol-variables to set"

echo "-- soot model LL91, steady, two iterations: the diagnostic output --"
PREMIX=../../../1D/premixed-flame/INPUT
# soot_run <kind: sa, soot, plain, sa-plain or sacomp>: the run of make_case.py soot <kind> (sacomp: the state of sa-plain,
# with turbulence = SAcomp), its exit status in OUTPUT/soot-<kind>.rc
soot_run() {
  rm -f INPUT/*
  local state=$1
  [[ $1 == sacomp ]] && state=sa-plain
  python3 -B make_case.py soot "$PREMIX" "$SOD/ic_x1.tec" INPUT "$state" || exit 1
  # faces 3 and 4 of the tube, one cell thick in j: a copy (type 0), not the extrapolation 400, which takes
  # three cells and with one cell reads ghost cells that no record fills
  awk '$5 == 3 || $5 == 4 { $6 = 0 } { print }' "$SOD/bc_x1.txt" > INPUT/bc.txt
  local turb=none soot=LL91
  [[ $1 == sa || $1 == sa-plain ]] && turb=SA
  [[ $1 == plain || $1 == sa-plain ]] && soot=none
  [[ $1 == sacomp ]] && turb=SAcomp && soot=none
  {
    printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 2' 'newrun = true' ''
    printf '%s\n' '[MOSE-Physics]' 'equations = navier-stokes' "turbulence = $turb" 'schmidt = 0.7' "soot-generation = $soot" ''
    printf '%s\n' '[MOSE-Numerics]' 'cfl = 0.5' 'time-accurate = false' 'time-scheme = euler' \
                  'integration-variables = cons' 'space-reconstruction = first-order' 'riemann-solver = HLLC' ''
    printf '%s\n' '[MOSE-IO]' 'sol-format = tecplot ascii'
  } > input.ini
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/soot-$1.log" 2>&1
  echo $? > "OUTPUT/soot-$1.rc"
  [[ -f OUTPUT/diagnostic.tec ]] && cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$1.tec"
}
rm -f OUTPUT/diagnostic-*.tec
# the references without soot: plain for the soot run, sa-plain (Spalart-Allmaras) for the sa run
soot_run plain
soot_run sa-plain
for kind in sa soot; do
  soot_run $kind
  rc=$(cat "OUTPUT/soot-$kind.rc")
  [[ $rc -eq 0 ]]
  check $? "$kind: the run ends normally (exit status $rc)"
  names="rho rhou rhov rhow rhoe rho*Np"
  ref=OUTPUT/diagnostic-plain.tec
  [[ $kind == sa ]] && names="$names rhoNuTilde" && ref=OUTPUT/diagnostic-sa-plain.tec
  python3 -B make_case.py diagnostic-check "OUTPUT/diagnostic-$kind.tec" "$ref" $names dt beta
  check $? "$kind: the diagnostic output names every variable it writes, dt under dt"
done

echo "-- turbulence = SAcomp stops at set-up --"
soot_run sacomp
rc=$(cat OUTPUT/soot-sacomp.rc)
SACOMP='turbulence = SAcomp: its compressibility correction solves'
if [[ $rc -ge 1 && $rc -le 128 && ! -f OUTPUT/diagnostic-sacomp.tec ]] && grep -a -q -- "$SACOMP" OUTPUT/soot-sacomp.log; then
  echo "  [ok  ] sacomp: exit code $rc, \"$SACOMP\""
else
  echo "  [FAIL] sacomp: exit code $rc (required 1 to 128, no diagnostic output), message \"$SACOMP\" $(grep -a -q -- "$SACOMP" OUTPUT/soot-sacomp.log && echo found || echo missing) (see OUTPUT/soot-sacomp.log)"
  fail=1
fi

# write_case <tag> <physics lines> [other sections]: the Sod tube and its input.ini
write_case() {
  local tag="$1" physics="$2" extra="${3:-}"
  rm -f INPUT/*
  if [[ "$physics" == *navier-stokes* ]]; then
    ln -sf "../$VISC/phase.txt" INPUT/phase.txt; ln -sf "../$VISC/thermo.dat" INPUT/thermo.dat
    ln -sf "../$VISC/transport.dat" INPUT/transport.dat
  else
    ln -sf "../$COMMON/phase.txt" INPUT/phase.txt; ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat
  fi
  if [[ "$tag" == nan-ic ]]; then
    python3 -B make_case.py nan "$SOD" || exit 1
  else
    python3 -B make_case.py sod "$SOD" || exit 1
  fi
  local cfl=0.9
  [[ "$tag" == diverge ]] && cfl=20
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 50
newrun = true

[MOSE-Physics]
$physics

[MOSE-Numerics]
cfl = $cfl
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
$extra
EOI
}

# run_case <tag> <message expected in the log> <physics lines> [other sections]
run_case() {
  local tag="$1" message="$2"
  write_case "$tag" "$3" "${4:-}"
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  if [[ $rc -ge 1 && $rc -le 128 ]] && grep -a -q -- "$message" "OUTPUT/$tag.log"; then
    echo "  [ok  ] $tag: exit code $rc, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1 to 128), message \"$message\" $(grep -a -q -- "$message" "OUTPUT/$tag.log" && echo found || echo missing) (see OUTPUT/$tag.log)"
    fail=1
  fi
}

echo "-- exit codes of the error stops --"
run_case diverge  "Integration failed"              "equations = euler"
run_case nan-ic   "Integration failed"              "equations = euler"
run_case ic-vars  "Number of variables in IOfield"  "equations = navier-stokes
turbulence = SA
schmidt = 0.7"
run_case rot-axis "Rotation axis has zero length"   "equations = euler
rotational-frame = rigid-body" "
[MOSE-Rotating-Frame]
omega = 1.0
axis = 0.0 0.0 0.0"

echo "-- a valid run ends with exit code 0 --"
write_case chem-excl "equations = euler
chemistry = finite-rate" "
[MOSE-Chemistry]
exclude-blocks = 1
ode-solver = sdirk4b"
"$LOCAL" > OUTPUT/chem-excl.log 2>&1
rc=$?
if [[ $rc -eq 0 ]] && grep -a -q -E 'Inactive blocks +1 *$' OUTPUT/chem-excl.log; then
  echo "  [ok  ] chem-excl: exit code 0, block 1 under Inactive blocks"
else
  echo "  [FAIL] chem-excl: exit code $rc (required 0), block 1 $(grep -a -q -E 'Inactive blocks +1 *$' OUTPUT/chem-excl.log && echo listed || echo not listed) under Inactive blocks (see OUTPUT/chem-excl.log)"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "ExitCodes: FAILED"
  exit 1
fi
echo "ExitCodes: all checks passed"
exit 0
