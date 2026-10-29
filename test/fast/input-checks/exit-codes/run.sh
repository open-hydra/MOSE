#!/usr/bin/env bash
#===============================================================================
#  Exit codes of the error stops (validation tier)
#
#  A run that MOSE stops on an error must end with a non-zero exit code, so
#  that scripts and CI see the failure, and with its message, not with a
#  crash signal; it writes no solution. Reuses the inputs of the 100-cell Sod
#  tube (test/1D/Sod79) and the gases of test/common and test/1D/diff1; no
#  new mesh or IC data is stored here.
#
#  Unknown probe variable. A probe that asks for a variable MOSE does not
#  know must stop the run at set-up with exit status 1 and a message naming
#  the variable, not crash: without the check the probe pointer stays undefined
#  and the run ends on a segmentation fault (exit status 139). The Sod tube of
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
#  An output that fills the turbulence column from arrays of different sizes
#  stops with a run-time error under bounds checking (with Spalart-Allmaras),
#  and one that names a variable fewer than it writes puts the density of
#  C(gr) under dt (without Spalart-Allmaras).
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
#  initial condition with fewer variables than the model needs (set-up), a
#  rotating frame with an axis of zero length (set-up),
#  chemistry = equilibrium without the elements of the species,
#  INPUT/composition.txt (set-up; without the stop the run would go on with an
#  equilibrium of no element), and a value outside the rule of the registry
#  for three forms of rule: preconditioning-eps-min = 1 (the
#  interval in (0, 1)), ode-relative-tol-species = 0 (a list of reals, > 0)
#  and index-position = 1 -1 1 1 of a probe (a list of integers, >= 0); without
#  the check the three would run with no message (set-up). The exit code must be
#  between 1 and 128 (128 is what ifx gives to ERROR STOP; above 128 the
#  shell reports a signal) and the log must hold the message.
#
#  A key of input.ini whose value does not read as the type of its parameter
#  stops the run at set-up, in every build, with exit code 1, the file, the
#  section, the key and the value as written, the STOP line of the set-up,
#  and no solution: the Sod tube with sol-dtime = NaN, Infinity and 0,5 (a
#  decimal comma), sol-diter = seven and 2.5, sol-overwrite = yes, and
#  levels = 2.0 in [MOSE-Multigrid], which the reader takes outside the
#  registry. Without the check the run would go on with the default of the
#  key, with exit code 0 and "Input file OK"; an integer written 2.5 or 2.0
#  is refused by gfortran and cut to 2 by ifx, so the text of an integer is
#  checked too.
#
#  A rotating frame whose axis or origin in input.ini holds a value that is
#  not finite (axis = NaN 0.0 0.0, origin = 0.0 Infinity 0.0) stops the run
#  at set-up, in every build, with exit code 1 and a message that names the
#  key: read as numbers, the values would go into every term of the frame.
#
#  A value of the initial condition that is not finite stops the run at
#  set-up, in every build: the tube with a NaN velocity in one cell, and the
#  tube with an infinite coordinate of one node, must end with exit code 1,
#  the place of the value (the variable or the coordinate, the block, the
#  cell or the node) and the STOP line of the set-up, and write no solution.
#  The check reads the bits of each value, so that it holds in a DEBUG
#  build, which traps the invalid operations (a test by comparison would end
#  the NaN case there with SIGFPE, exit code 136, at the end of the set-up),
#  and with -ffast-math, an option of the RELEASE build, which lets gfortran
#  fold isnan to false.
#
#  A time table of a boundary condition with a time or a value that is not
#  finite stops the run at set-up, in every build, with the file, the line
#  and the column: the duct of BCInput, an inlet 402 whose total pressure
#  comes from the table p0.dat (../../boundary-conditions/bc-input/make_case.py case p0.dat), with
#  the pressure of line 3 written as NaN and as Infinity, and the time of
#  line 4 as NaN, must end with exit code 1 and the message, and write no
#  solution. The check reads the bits: every comparison with a NaN is false,
#  so the check of decreasing times lets a NaN time through, and a NaN
#  pressure would close the face of the inlet in silence.
#
#  A value of a property line of bc.txt that is not finite stops the run at
#  set-up, in every build, with the line, the record, its type and the name of
#  the value: the duct of BCInput with its inlet made an inlet 401 whose total
#  pressure is NaN and Infinity, or whose flow direction alpha is NaN, with its
#  inlet 402 of an infinite total temperature, and, with the Navier-Stokes
#  equations, with a wall 302 on a face cell whose roughness, a value with a
#  default, is NaN, must end with exit code 1 and the message, and write no
#  solution. The check of an empty value lets a NaN through, as every
#  comparison with a NaN is false.
#
#  A value of the Q2D file of an assigned state (410) that is not finite stops
#  the run at set-up, in every build, with the file, the zone, the line of its
#  header, and the variable and the cell of the value: the duct of BCInput
#  with its inlet made a face of type 410 that reads q2d.dat, two zones of
#  the face at t = 0 and t = 1 s, with the pressure of the second zone NaN,
#  or its SOLUTIONTIME Infinity, must end with exit code 1 and the message,
#  and write no solution.
#  So does a Q2D file that cannot be read, with the file and, for a zone, the
#  zone and the line of its header: the same duct with q2d.dat missing, cut
#  in the values of the second zone, with a word among them, with I=two in
#  the header of the second zone, and with one value more in the first zone
#  (the line after its data is not the header of the next zone), and so does
#  a file whose zones go back in time (the first zone at t = 2 s, the second
#  at t = 1 s). Without these checks the values not read would keep those of
#  the zone before, or the zones after it would be read shifted; and the times
#  of the zones would be taken in the order of the file, where the
#  interpolation wants them increasing. The run would go on.
#  So does a zone header that gives no I= or no J=, the only form of I and J
#  that the reader takes (upper case, with no blank before =): the same duct
#  with I = 3, or j=2, in the header of the second zone. Without the check I
#  or J would keep its default, 2, without a word, and the run would go on.
#  So does a file that lacks a part of the state of the cell across the face,
#  which the ghost cells of a face of type 410 take whole (the density rho(N)
#  of at least one species, u, v, w and p), with the file and the variables
#  it lacks: the same duct with p, with u and w, or with rho(1), the density
#  of its only species, taken out of the file. Without the check the ghost
#  cells would take the value of the cell for a missing p, u or w, and the run
#  would end with exit code 0; without any rho(N) they would have zero
#  density, and the run would end at its first step with a run-time error
#  (gfortran with bounds checking).
#  So does a file that leaves a cell of a face of type 410 without a value at
#  one of its times, or whose zones give nodes that are not those of the
#  faces, with the file, the block and the record, and for the nodes the zone
#  and the node: the duct of BCInput in two blocks with faces of type 410 on
#  its six kinds of faces, read from three files in the form that ATLAS
#  writes (make_case.py duct410), the file of its faces 5 and 6 with the
#  zones of the two blocks under each other's title, with no zone of block 2,
#  with no zone of block 2 at the second time, with the zones of block 2 cut
#  to their first 25 cells, with the zone of block 1 at the second time
#  holding that of block 2, and with its nodes moved by 2e-2 of the shortest
#  edge of these faces. Without these checks such cells would keep the zero
#  state, and the run would end at its first step with a run-time error that
#  names no file, or take that state at the second time; and the values of a
#  zone of another block, of another time or of another mesh would go to the
#  face without a word. The files as written must run, with exit code
#  0, and so must the same files with the nodes of the faces 5 and 6 moved by
#  0.5e-2 of that edge, or with the coordinates written with 15 significant
#  digits, with the solution of the files as written, bit for bit; and so
#  must a restart from that solution, whose nodes the check reads again.
#
#  An input.ini that cannot be read stops the run at set-up with its reason:
#  input.ini missing, empty, and of comment lines only (no section, an error
#  of the INI reader). Without this check the run would stop on the first
#  required key, naming that key and not the file, or, with the empty file,
#  end inside the reader with a run-time error.
#
#  A restart from a solution whose time or counter is not finite stops the
#  run at set-up, in every build, with exit code 1, the file, the name of the
#  value and the STOP line of the set-up, and leaves the solution as it was:
#  the Sod tube written as OUTPUT/field.tec (Tecplot ASCII), its SOLUTIONTIME
#  made NaN, or its MOSE_LAST_OUTPUT_TIME made Infinity, then read by a
#  restart (newrun = false).
#
#  And a valid run must end with exit code 0: with finite-rate chemistry and
#  exclude-blocks = 1 the tube holding the gas of test/1D/Fer14 with its
#  reaction mechanism (its INPUT/chemistry-info.txt and
#  INPUT/chemistry-Arrhenius.dat; make_case.py fer14) must run, its set-up
#  summary listing block 1 under "Inactive blocks" (a set-up that deallocates
#  that list before it is allocated ends with exit code 153 with ifx,
#  "allocatable array or pointer is not allocated"); the excluded block takes
#  no chemistry.
#  exclude-blocks with a word that is not a block number (x 1, 1 x) or with a
#  block beyond the mesh (1000) must stop the run at set-up with exit code 1
#  and a message that names the key, the value and the word or the block,
#  and write no solution: the same tube with its mechanism; the words after
#  one that does not read are not read (a read beyond the list would be a
#  run-time error with bounds checks), and a block beyond the mesh names no
#  block.
#  chemistry = finite-rate needs the reaction mechanism: the tube holding the
#  gas of test/1D/Fer14 without its INPUT/chemistry-info.txt must stop at
#  set-up with exit code 1, the message that names the file and the STOP line
#  of the set-up, before any iteration, without the warning that FLINT writes
#  for a mechanism it does not know, and write no solution. Without the stop
#  the run would go on with no rate table: FLINT answers every call of the
#  source terms with its bail-out value, -1 in every component, which the
#  integrators of OSLO other than cvode take as a rate of -1 for every
#  density and for the temperature: the run goes on with its composition and
#  its temperature changed, in silence.
#  A failed integration of the finite-rate chemistry stops the run: the same
#  tube with the mechanism of test/1D/Fer14 and a relative tolerance of the
#  species of 1e-30, which the default integrator, H-radau5, refuses at every
#  call (RADAU5 takes a relative tolerance above 10 times its rounding unit,
#  1.1e-19, and returns -1, input not consistent, without integrating), must
#  end at its first step with exit code 1, one line that gives the
#  integrator, its return code and what it means, the block, the cell (every
#  cell fails: the line gives the failure of the smallest block, k, j, i, cell
#  1, 1, 1) and the state of the cell, the STOP line, and no solution; with 4
#  threads the line is the same as with one, and on two multigrid levels
#  (steady, with the INPUT/bc2.txt of the steam case below) the line names the
#  coarse level, where the run starts. Without the stop no cell would take its
#  chemistry and the run would end with exit code 0.
#  exclude-blocks holds for chemistry = equilibrium too, on every multigrid
#  level: the tube filled with steam far from equilibrium (0.2 % of hydrogen
#  and 1.6 % of oxygen, at 1500 and 1200 K; the 13 C-H-O species of
#  test/2D/euler/rocket-nozzle/equilibrium; make_case.py steam), its one
#  block excluded, must give the bits of the run without chemistry, on one
#  grid (time-accurate) and on two multigrid levels (steady); with the block
#  not excluded the equilibrium changes the state.
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
CHO=../../../2D/euler/rocket-nozzle/equilibrium/INPUT
FER=../../../1D/Fer14/INPUT
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
  # faces 3 and 4 of the tube, one cell thick in j: the null type 0, no flux through them
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
  case "$tag" in
    nan-ic)   python3 -B make_case.py nan "$SOD" || exit 1 ;;
    inf-node) python3 -B make_case.py inf "$SOD" || exit 1 ;;
    *)        python3 -B make_case.py sod "$SOD" || exit 1 ;;
  esac
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

# write_fer14 <mechanism> <tag> <physics lines> [other sections]: the tube of write_case holding the gas of test/1D/Fer14
# (make_case.py fer14); with mechanism = yes, also its reaction mechanism, INPUT/chemistry-info.txt and
# INPUT/chemistry-Arrhenius.dat
write_fer14() {
  local mechanism="$1"
  shift
  write_case "$@"
  rm -f INPUT/phase.txt INPUT/thermo.dat
  ln -sf "../$FER/phase.txt" INPUT/phase.txt; ln -sf "../$FER/thermo.dat" INPUT/thermo.dat
  if [[ $mechanism == yes ]]; then
    ln -sf "../$FER/chemistry-info.txt" INPUT/chemistry-info.txt
    ln -sf "../$FER/chemistry-Arrhenius.dat" INPUT/chemistry-Arrhenius.dat
  fi
  python3 -B make_case.py fer14 "$SOD" INPUT/phase.txt || exit 1
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

# run_value <tag> <place>: the tube whose initial condition holds a value that is not finite at <place> (make_case.py
# nan, inf) must stop at set-up with exit code 1, the message that names the place, the STOP line of the set-up and no
# solution
run_value() {
  local tag="$1" place="$2"
  local message="holds a value that is not finite (NaN or infinity): $place" found=missing stop=missing solution=none
  write_case "$tag" "equations = euler"
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  grep -a -q -F -- "STOP: initial condition error (see above)" "OUTPUT/$tag.log" && stop=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $stop == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, STOP line $stop, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}

echo "-- exit codes of the error stops --"
run_case diverge  "Integration failed"              "equations = euler"
run_value nan-ic   "variable 2 (u) of block 1, cell i, j, k = 50, 1, 1"
run_value inf-node "coordinate y of block 1, node i, j, k = 20, 1, 1"
run_case ic-vars  "Number of variables in IOfield"  "equations = navier-stokes
turbulence = SA
schmidt = 0.7"
run_case rot-axis "Rotation axis has zero length"   "equations = euler
rotational-frame = rigid-body" "
[MOSE-Rotating-Frame]
omega = 1.0
axis = 0.0 0.0 0.0"
run_case rot-axis-nan "Rotating frame: a value of axis in input.ini is not finite (NaN or infinity)" "equations = euler
rotational-frame = rigid-body" "
[MOSE-Rotating-Frame]
omega = 1.0
axis = NaN 0.0 0.0"
run_case rot-origin-inf "Rotating frame: a value of origin in input.ini is not finite (NaN or infinity)" "equations = euler
rotational-frame = rigid-body" "
[MOSE-Rotating-Frame]
omega = 1.0
origin = 0.0 Infinity 0.0"
run_case chem-elements "chemistry = equilibrium needs the elements of the species" "equations = euler
chemistry = equilibrium"
run_case rule-interval "preconditioning-eps-min must be in (0, 1)" "equations = euler" "preconditioning-eps-min = 1"
run_case rule-reals "ode-relative-tol-species must be > 0 (value 1 of the list)" "equations = euler" "
[MOSE-Chemistry]
ode-relative-tol-species = 0"
run_case rule-integers "index-position must be >= 0 (value 2 of the list)" "equations = euler" "
[MOSE-Probes]
probe1 = middle

[middle]
variables = p
index-position = 1 -1 1 1"

echo "-- a key of input.ini whose value does not read as its type stops the run at set-up --"
# run_key <tag> <lines> <message>: the Sod tube of write_case with the lines added to its input.ini must stop at set-up with
# exit code 1, the message, the STOP line of the set-up, and no solution
run_key() {
  local tag="$1" lines="$2" message="$3" found=missing stop=missing solution=none
  write_case "$tag" "equations = euler" "$lines"
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  grep -a -q -F -- "STOP: input file error (see above)" "OUTPUT/$tag.log" && stop=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $stop == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, STOP line $stop, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
run_key key-nan $'[MOSE-IO]\nsol-dtime = NaN' \
  "[ERROR] input.ini, [MOSE-IO] sol-dtime = NaN: the value is not a finite real number"
run_key key-inf $'[MOSE-IO]\nsol-dtime = Infinity' \
  "[ERROR] input.ini, [MOSE-IO] sol-dtime = Infinity: the value is not a finite real number"
run_key key-comma $'[MOSE-IO]\nsol-dtime = 0,5' \
  "[ERROR] input.ini, [MOSE-IO] sol-dtime = 0,5: the value is not a finite real number"
run_key key-word $'[MOSE-IO]\nsol-diter = seven' \
  "[ERROR] input.ini, [MOSE-IO] sol-diter = seven: the value is not an integer"
run_key key-integer $'[MOSE-IO]\nsol-diter = 2.5' \
  "[ERROR] input.ini, [MOSE-IO] sol-diter = 2.5: the value is not an integer"
run_key key-logical $'[MOSE-IO]\nsol-overwrite = yes' \
  "[ERROR] input.ini, [MOSE-IO] sol-overwrite = yes: the value is not a logical (true or false)"
run_key key-levels $'[MOSE-Multigrid]\nlevels = 2.0' \
  "[ERROR] input.ini, [MOSE-Multigrid] levels = 2.0: the value is not an integer"

echo "-- a time table with a value that is not finite stops the run at set-up --"
# run_table <tag> <line> <column> <token>: the duct of the inlet 402 of BCInput, the field <column> of line <line> of its
# table p0.dat written as <token>, must stop at set-up with exit code 1, the file, the line and the column, and no solution
run_table() {
  local tag="$1" line="$2" col="$3" tok="$4"
  local message="[ERROR] Time table p0.dat, line $line, column $col: a value that is not finite (NaN or infinity)"
  local found=missing solution=none
  rm -f INPUT/*
  ln -sf "../$COMMON/phase.txt" INPUT/phase.txt; ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat
  python3 -B ../../boundary-conditions/bc-input/make_case.py case p0.dat || exit 1
  awk -v n="$line" -v c="$col" -v t="$tok" 'NR == n { $c = t } { print }' p0.dat > p0.new && mv p0.new p0.dat
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 40
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
EOI
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
run_table table-nan  3 2 NaN
run_table table-inf  3 2 Infinity
run_table table-time 4 1 NaN

echo "-- a value of a property line of bc.txt that is not finite stops the run at set-up --"
# run_bc_case_files [equations]: the inputs of the duct of the inlet 402 of BCInput (../../boundary-conditions/bc-input/make_case.py case p0.dat), in
# air (Euler) or in the viscous air of test/common (navier-stokes), and its input.ini
run_bc_case_files() {
  local equations="${1:-euler}" gas=$COMMON
  rm -f INPUT/*
  [[ $equations == navier-stokes ]] && gas=$VISC
  ln -sf "../$gas/phase.txt" INPUT/phase.txt; ln -sf "../$gas/thermo.dat" INPUT/thermo.dat
  [[ $equations == navier-stokes ]] && ln -sf "../$gas/transport.dat" INPUT/transport.dat
  python3 -B ../../boundary-conditions/bc-input/make_case.py case p0.dat || exit 1
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 40
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = $equations

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOI
}
# run_bc <tag> <equations> <awk program> <message>: the duct, its bc.txt changed by the awk program, must stop at set-up with
# exit code 1 and the message, and write no solution
run_bc() {
  local tag="$1" equations="$2" edit="$3" message="$4" found=missing solution=none
  run_bc_case_files "$equations"
  awk "$edit" INPUT/bc.txt > INPUT/bc.new && mv INPUT/bc.new INPUT/bc.txt
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
# the inlet (record 1, lines 1 and 2) made an inlet 401: T0, p0, alpha, beta, rel_fac, mass fraction
run_bc bc-p0-nan euler 'NR == 1 { $6 = 401 } NR == 2 { $0 = "300.0, NaN, normal, normal, 1.0, 1.0" } { print }' \
  "[ERROR] bc.txt line 2 (record 1, type 401): the value of p0 is not finite (NaN or infinity)"
run_bc bc-p0-inf euler 'NR == 1 { $6 = 401 } NR == 2 { $0 = "300.0, Infinity, normal, normal, 1.0, 1.0" } { print }' \
  "[ERROR] bc.txt line 2 (record 1, type 401): the value of p0 is not finite (NaN or infinity)"
run_bc bc-alpha euler 'NR == 1 { $6 = 401 } NR == 2 { $0 = "300.0, 1.2E+05, NaN, normal, 1.0, 1.0" } { print }' \
  "[ERROR] bc.txt line 2 (record 1, type 401): the value of alpha is not finite (NaN or infinity)"
# the inlet 402: T0, p0 file, alpha, beta, rel_fac, mass fraction
run_bc bc-402-t0 euler 'NR == 2 { $0 = "Infinity, p0.dat, normal, normal, 1.0, 1.0" } { print }' \
  "[ERROR] bc.txt line 2 (record 1, type 402): the value of T0 is not finite (NaN or infinity)"
# record 3 (line 5, face 3 of cell 1) made a wall 302 with its property line: T_wall, k_rough, eps_wall
run_bc bc-wall-k navier-stokes 'NR == 5 { $6 = 302; print; print "300.0, NaN, 0.0"; next } { print }' \
  "[ERROR] bc.txt line 6 (record 3, type 302): the value of k_rough is not finite (NaN or infinity)"

echo "-- a value of the Q2D file of a face 410 that is not finite stops the run at set-up --"
# run_q2d <tag> <time of zone 2> <pressure of zone 2> <message> [edit]: the duct of BCInput, its inlet (record 1) made a face of
# type 410 that reads q2d.dat, two zones of that face (block 1, face 1: one cell) at t = 0 and at the time given, the state of
# the second one with the pressure given, the file then changed by the awk program edit (missing: no file), must stop at
# set-up with exit code 1 and the message, and write no solution
run_q2d() {
  local tag="$1" time2="$2" p2="$3" message="$4" edit="${5:-}" found=missing solution=none
  run_bc_case_files
  awk 'NR == 1 { $6 = 410 } NR == 2 { $0 = "q2d.dat" } { print }' INPUT/bc.txt > INPUT/bc.new && mv INPUT/bc.new INPUT/bc.txt
  {
    echo 'VARIABLES = "x" "y" "z" "rho(1)" "u" "v" "w" "p"'
    for zone in 1 2; do
      local t=0.0 p=1.0E+05
      [[ $zone == 2 ]] && t=$time2 && p=$p2
      echo "ZONE T=\"B1_T$zone\", I=2, J=2, SOLUTIONTIME=$t, DATAPACKING=BLOCK, VARLOCATION=([4-8]=CELLCENTERED)"
      printf '%s\n' '0.0 0.0 0.0 0.0' '0.0 0.01 0.0 0.01' '-5.0E-04 -5.0E-04 5.0E-04 5.0E-04' '1.16' '50.0' '0.0' '0.0' "$p"
    done
  } > q2d.dat
  if [[ $edit == missing ]]; then
    rm -f q2d.dat
  elif [[ -n $edit ]]; then
    awk "$edit" q2d.dat > q2d.new && mv q2d.new q2d.dat
  fi
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
run_q2d q2d-nan 1.0 NaN \
  "[ERROR] Q2D file q2d.dat, zone 2 (line 11), variable p, cell 1, 1: a value that is not finite (NaN or infinity)"
run_q2d q2d-time Infinity 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, line 11 (zone 2): SOLUTIONTIME is not finite (NaN or infinity)"

echo "-- a Q2D file of a face 410 that lacks a part of the state across the face stops the run at set-up --"
# the file of run_q2d with finite values (zone 2 at 1 s) without p, without u and w, without rho(1) (the density of the only
# species): the name out of the VARIABLES line, the values out of both zones (zone 1: lines 6, 7, 9, 10 for rho(1), u, w, p;
# zone 2: lines 15, 16, 18, 19), the CELLCENTERED range of both headers shortened
run_q2d q2d-no-p 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat lacks p: the file of a face of type 410 gives the whole state of the cell across the face, the density rho(N) of at least one species N, u, v, w and p" \
  'NR == 1 { sub(/ "p"/, "") } NR == 2 || NR == 11 { sub(/4-8/, "4-7") } NR == 10 || NR == 19 { next } { print }'
run_q2d q2d-no-uw 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat lacks u, w: the file of a face of type 410 gives the whole state of the cell across the face, the density rho(N) of at least one species N, u, v, w and p" \
  'NR == 1 { sub(/ "u"/, ""); sub(/ "w"/, "") } NR == 2 || NR == 11 { sub(/4-8/, "4-6") } NR == 7 || NR == 9 || NR == 16 || NR == 18 { next } { print }'
run_q2d q2d-no-rho 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat lacks rho(N): the file of a face of type 410 gives the whole state of the cell across the face, the density rho(N) of at least one species N, u, v, w and p" \
  'NR == 1 { sub(/ "rho.1."/, "") } NR == 2 || NR == 11 { sub(/4-8/, "4-7") } NR == 6 || NR == 15 { next } { print }'

echo "-- a Q2D file of a face 410 that cannot be read stops the run at set-up --"
# the file of run_q2d with finite values (zone 2 at 1 s), 19 lines (zone 1: lines 2 to 10, its header and the values of x, y,
# z, rho(1), u, v, w, p; zone 2: lines 11 to 19), missing, cut before its last line (the value of p of zone 2), with a word for
# the value of u of zone 2 (line 16), with I=two in the header of zone 2, with I = 3 or j=2 there (not the form I=, J= that the
# reader takes), with the line of p of zone 1 written twice, with the time of zone 1 made 2 s (after that of zone 2)
run_q2d q2d-missing 1.0 1.0E+05 "[ERROR] Q2D file q2d.dat cannot be opened" missing
run_q2d q2d-short 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, zone 2 (line 11): the values of p: the file ends before all the values that I and J give" \
  'NR < 19'
run_q2d q2d-word 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, zone 2 (line 11): the values of u: a value that is not a number" \
  'NR == 16 { $0 = "fifty" } { print }'
run_q2d q2d-header 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, line 11 (zone 2): the value of I cannot be read: two" \
  'NR == 11 { sub(/I=2/, "I=two") } { print }'
run_q2d q2d-header-blank 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, line 11 (zone 2): no I= in the header; give the nodes as I=<n>, J=<n>, in upper case with no blank around =" \
  'NR == 11 { sub(/I=2/, "I = 3") } { print }'
run_q2d q2d-header-case 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, line 11 (zone 2): no J= in the header; give the nodes as I=<n>, J=<n>, in upper case with no blank around =" \
  'NR == 11 { sub(/J=2/, "j=2") } { print }'
run_q2d q2d-long 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, zone 1 (line 2): the line after its data is not the header of zone 2 (line 12)" \
  'NR == 10 { print } { print }'
# the first zone at t = 2 s, the second at t = 1 s: the zones go back in time
run_q2d q2d-order 1.0 1.0E+05 \
  "[ERROR] Q2D file q2d.dat, line 11 (zone 2): SOLUTIONTIME 1.00000E+00 is less than 2.00000E+00, the time of a zone above it: the zones must come in increasing time" \
  'NR == 2 { sub(/SOLUTIONTIME=0.0/, "SOLUTIONTIME=2.0") } { print }'

echo "-- a file of faces of type 410 that leaves a cell without a value, or whose nodes are not those of the faces, stops the run at set-up --"
# run_q2d_mesh <tag> <form> <message>: the duct of BCInput in two blocks with faces of type 410 on its six kinds of faces, read
# from three files in the form that ATLAS writes (make_case.py duct410 <form>), 5 steps. With a message, it must stop at set-up
# with exit code 1 and the message, and write no solution; with "valid", it must end with exit code 0, its solution kept as
# OUTPUT/q2d-mesh-valid.tec; with "same", it must end with exit code 0 and give that solution, bit for bit
run_q2d_mesh() {
  local tag="$1" form="$2" message="$3" found=missing solution=none
  rm -f INPUT/* q2d-ends.dat q2d-jfaces.dat q2d-side.dat
  ln -sf "../$COMMON/phase.txt" INPUT/phase.txt; ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat
  python3 -B make_case.py duct410 "$form" || exit 1
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 5
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
EOI
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  if [[ $message == valid || $message == same ]]; then
    [[ $message == valid && $rc -eq 0 ]] && cp OUTPUT/field.tec OUTPUT/q2d-mesh-valid.tec
    if [[ $rc -eq 0 ]] && cmp -s OUTPUT/field.tec OUTPUT/q2d-mesh-valid.tec; then
      if [[ $message == valid ]]; then
        echo "  [ok  ] $tag: exit code 0"
      else
        echo "  [ok  ] $tag: exit code 0, the solution of the files as written, bit for bit"
      fi
    else
      echo "  [FAIL] $tag: exit code $rc (required 0), solution $(cmp -s OUTPUT/field.tec OUTPUT/q2d-mesh-valid.tec && echo same || echo different) (see OUTPUT/$tag.log)"
      fail=1
    fi
    return
  fi
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  [[ -f OUTPUT/field.tec ]] && solution=written
  if [[ $rc -eq 1 && $found == found && $solution == none ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
rm -f OUTPUT/q2d-mesh-valid.tec
run_q2d_mesh q2d-mesh-valid valid valid
# a restart from that solution: the check reads its nodes again
sed -i 's/^newrun = true$/newrun = false/' input.ini
"$LOCAL" > OUTPUT/q2d-mesh-restart.log 2>&1
rc=$?
if [[ $rc -eq 0 ]] && grep -a -q 'MOSE | Iter' OUTPUT/q2d-mesh-restart.log; then
  echo "  [ok  ] q2d-mesh-restart: exit code 0, a restart from the solution of the files as written"
else
  echo "  [FAIL] q2d-mesh-restart: exit code $rc (required 0, with its steps) (see OUTPUT/q2d-mesh-restart.log)"
  fail=1
fi
run_q2d_mesh q2d-mesh-moved-5e-5 moved-5e-5 same
run_q2d_mesh q2d-mesh-digits-15 digits-15 same
run_q2d_mesh q2d-mesh-swapped swapped \
  "[ERROR] Q2D file q2d-side.dat, zone 1 (line 2): node I=1, J=1 of the zone, a corner of the face of record 205 (block 1, cell 1, 1, 1, face 5), is at"
run_q2d_mesh q2d-mesh-moved-2e-4 moved-2e-4 \
  "[ERROR] Q2D file q2d-side.dat, zone 1 (line 2): node I=1, J=1 of the zone, a corner of the face of record 205 (block 1, cell 1, 1, 1, face 5), is at"
run_q2d_mesh q2d-mesh-no-block-2 no-block-2 \
  "[ERROR] Q2D file q2d-side.dat has no zone of block 2, whose 50 cells of faces of type 410 read it (the first: record 559, cell 1, 1, 2, face 6)"
run_q2d_mesh q2d-mesh-no-time-2 no-time-2 \
  "[ERROR] Q2D file q2d-side.dat: block 2 has no zone at SOLUTIONTIME 1.00000E-03, a time of the file, and its cells of faces of type 410 need a value at every time of the file (the first: record 559, cell 1, 1, 2, face 6)"
run_q2d_mesh q2d-mesh-short-zone short-zone \
  "[ERROR] Q2D file q2d-side.dat: at SOLUTIONTIME 0.00000E+00 no zone of block 2 holds cell 26, 1 of the face of record 584 (cell 26, 1, 2, face 6); the first zone of the block at that time, zone 2 (line 559), has I=26, J=2, that is 25 by 1 cells, along i and j of the face"
run_q2d_mesh q2d-mesh-time-2-other-block time-2-other-block \
  "[ERROR] Q2D file q2d-side.dat, zone 3 (line 1116): node I=1, J=1 of the zone, a corner of the face of record 205 (block 1, cell 1, 1, 1, face 5), is  5.000E-01 from the same node in zone 1 (line 2)"
rm -f q2d-ends.dat q2d-jfaces.dat q2d-side.dat OUTPUT/q2d-mesh-valid.tec

echo "-- an input.ini that cannot be read stops the run at set-up --"
# run_ini <tag> <reason expected in the log>: the Sod tube with the input.ini left by the caller
run_ini() {
  local tag="$1" message="$2"
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  if [[ $rc -ge 1 && $rc -le 128 && ! -f OUTPUT/field.tec ]] && grep -a -q -- "$message" "OUTPUT/$tag.log"; then
    echo "  [ok  ] $tag: exit code $rc, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1 to 128, no solution), message \"$message\" $(grep -a -q -- "$message" "OUTPUT/$tag.log" && echo found || echo missing) (see OUTPUT/$tag.log)"
    fail=1
  fi
}
write_case ini-missing "equations = euler"; rm -f input.ini
run_ini ini-missing "input.ini cannot be opened"
write_case ini-empty "equations = euler"; : > input.ini
run_ini ini-empty "input.ini is empty"
write_case ini-comments "equations = euler"
printf '%s\n' '; the keys of the run were here' '# iter-threshold = 50' > input.ini
run_ini ini-comments "input.ini cannot be read by the INI reader"

echo "-- a restart from a solution whose time or counter is not finite stops the run at set-up --"
# run_restart <tag> <sed program> <message>: the Sod tube of write_case, its solution OUTPUT/field.tec (Tecplot ASCII) changed by
# the sed program, then a restart from it (newrun = false), must stop at set-up with exit code 1, the message and the STOP line
# of the set-up, and leave the solution as it was
run_restart() {
  local tag="$1" edit="$2" message="$3" found=missing stop=missing solution=kept
  write_case "$tag" "equations = euler"
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag-first.log" 2>&1
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] $tag: the run that writes the solution to restart from failed (see OUTPUT/$tag-first.log)"
    fail=1
    return
  fi
  sed "$edit" OUTPUT/field.tec > "OUTPUT/$tag-restart.tec" && cp "OUTPUT/$tag-restart.tec" OUTPUT/field.tec
  sed -i 's/^newrun = true$/newrun = false/' input.ini
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  local rc=$?
  grep -a -q -F -- "$message" "OUTPUT/$tag.log" && found=found
  grep -a -q -F -- "STOP: initial condition error (see above)" "OUTPUT/$tag.log" && stop=found
  cmp -s OUTPUT/field.tec "OUTPUT/$tag-restart.tec" || solution=written
  rm -f OUTPUT/field.tec
  if [[ $rc -eq 1 && $found == found && $stop == found && $solution == kept ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $found, STOP line $stop, solution $solution (see OUTPUT/$tag.log)"
    fail=1
  fi
}
run_restart restart-time 's/SOLUTIONTIME=[^,]*/SOLUTIONTIME=NaN/' \
  "[ERROR] OUTPUT/field.tec: the solution time is not finite (NaN or infinity)"
run_restart restart-counter 's/^DATASETAUXDATA MOSE_LAST_OUTPUT_TIME = .*$/DATASETAUXDATA MOSE_LAST_OUTPUT_TIME = "Infinity"/' \
  "[ERROR] OUTPUT/field.tec: MOSE_LAST_OUTPUT_TIME is not finite (NaN or infinity)"

echo "-- a valid run ends with exit code 0 --"
write_fer14 yes chem-excl "equations = euler
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

echo "-- exclude-blocks with a word that is not a block number, or a block beyond the mesh, stops at set-up --"
for value in "x 1" "1 x" "1000"; do
  tag="chem-excl-$(printf '%s' "$value" | tr ' ' '_')"
  write_fer14 yes "$tag" "equations = euler
chemistry = finite-rate" "
[MOSE-Chemistry]
exclude-blocks = $value
ode-solver = sdirk4b"
  if [[ "$value" == 1000 ]]; then
    message="[ERROR] input.ini, [MOSE-Chemistry] exclude-blocks = 1000: block 1000 is not a block of the mesh (1 to 1)"
  else
    message="[ERROR] input.ini, [MOSE-Chemistry] exclude-blocks = $value: \"x\" is not a block number"
  fi
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  rc=$?
  if [[ $rc -eq 1 ]] && grep -a -q -F -- "$message" "OUTPUT/$tag.log" && [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [ok  ] $tag: exit code 1, \"$message\""
  else
    echo "  [FAIL] $tag: exit code $rc (required 1), message \"$message\" $(grep -a -q -F -- "$message" "OUTPUT/$tag.log" && echo found || echo missing), solution $([[ -f OUTPUT/field.tec ]] && echo written || echo none) (see OUTPUT/$tag.log)"
    fail=1
  fi
done

echo "-- chemistry = finite-rate without its reaction mechanism stops at set-up --"
write_fer14 no chem-no-mechanism "equations = euler
chemistry = finite-rate"
message="[ERROR] chemistry = finite-rate needs the reaction mechanism, INPUT/chemistry-info.txt: the file cannot be opened"
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/chem-no-mechanism.log 2>&1
rc=$?
found=missing; stop=missing; iter=none; warning=none; solution=none
grep -a -q -F -- "$message" OUTPUT/chem-no-mechanism.log && found=found
grep -a -q -F -- "STOP: setup error (see above)" OUTPUT/chem-no-mechanism.log && stop=found
grep -a -q -F -- "MOSE | Iter" OUTPUT/chem-no-mechanism.log && iter=written
grep -a -q -F -- "Explicit procedure for" OUTPUT/chem-no-mechanism.log && warning=written
[[ -f OUTPUT/field.tec ]] && solution=written
if [[ $rc -eq 1 && $found == found && $stop == found && $iter == none && $warning == none && $solution == none ]]; then
  echo "  [ok  ] chem-no-mechanism: exit code 1, \"$message\""
else
  echo "  [FAIL] chem-no-mechanism: exit code $rc (required 1), message \"$message\" $found, STOP line $stop, iteration lines" \
       "$iter (required none), warning of FLINT $warning (required none), solution $solution (see OUTPUT/chem-no-mechanism.log)"
  fail=1
fi

echo "-- a failed integration of the finite-rate chemistry stops the run (1 and 4 threads) --"
# chem_fail <threads>: the tube of write_fer14 with its mechanism and ode-relative-tol-species = 1e-30, which the default
# integrator, H-radau5, refuses at every call, run for at most 3 iterations with <threads> threads; writes
# OUTPUT/chem-integrator-<threads>.{log,rc} and OUTPUT/chem-integrator-<threads>.solution (written or none)
chem_fail() {
  local nt="$1" solution=none
  write_fer14 yes "chem-integrator-$nt" "equations = euler
chemistry = finite-rate" "
[MOSE-Chemistry]
ode-relative-tol-species = 1e-30"
  sed -i 's/^iter-threshold = 50$/iter-threshold = 3/' input.ini
  rm -f OUTPUT/field.tec
  OMP_NUM_THREADS=$nt "$LOCAL" > "OUTPUT/chem-integrator-$nt.log" 2>&1
  echo $? > "OUTPUT/chem-integrator-$nt.rc"
  [[ -f OUTPUT/field.tec ]] && solution=written
  echo $solution > "OUTPUT/chem-integrator-$nt.solution"
  rm -f OUTPUT/field.tec
}
prefix="[ERROR] Chemistry: the integrator H-radau5 failed (return code -1: input is not consistent) in block 1, cell i, j, k = 1, 1, 1, over a chemistry step of"
skip=0
for nt in 1 4; do
  [[ $nt -eq 4 && $skip -eq 1 ]] && break
  chem_fail $nt
  log="OUTPUT/chem-integrator-$nt.log"
  rc=$(cat "OUTPUT/chem-integrator-$nt.rc")
  n=$(grep -a -c -F -- "$prefix" "$log")
  stop=missing; threads=used
  grep -a -q -F -- "STOP: chemistry integration error (see above)" "$log" && stop=found
  [[ $nt -eq 4 ]] && ! grep -a -q "Number of threads --> *4\$" "$log" && threads="not used (see the line Number of threads)"
  if [[ $rc -eq 1 && $n -eq 1 && $stop == found && $(cat "OUTPUT/chem-integrator-$nt.solution") == none && $threads == used ]]; then
    echo "  [ok  ] chem-integrator, $nt thread(s): exit code 1, the line of the failure once, \"$prefix ...\""
  else
    echo "  [FAIL] chem-integrator, $nt thread(s): exit code $rc (required 1), the line \"$prefix ...\" written $n times (required" \
         "once), STOP line $stop, solution $(cat "OUTPUT/chem-integrator-$nt.solution") (required none), $nt threads $threads (see $log)"
    fail=1
  fi
  # a MOSE built without OpenMP writes "Serial execution" in place of the number of threads: it has one thread only
  if [[ $nt -eq 1 ]] && grep -a -q "^ Serial execution\$" "$log"; then
    skip=1
  fi
done
if [[ $skip -eq 1 ]]; then
  echo "  [skip] chem-integrator, 4 threads: MOSE built without OpenMP (OUTPUT/chem-integrator-1.log: Serial execution)"
elif [[ "$(grep -a -F -- "$prefix" OUTPUT/chem-integrator-1.log)" == "$(grep -a -F -- "$prefix" OUTPUT/chem-integrator-4.log)" ]] && \
     [[ -n "$(grep -a -F -- "$prefix" OUTPUT/chem-integrator-1.log)" ]]; then
  echo "  [ok  ] chem-integrator: the line of the failure with 4 threads is the one with 1 thread"
else
  echo "  [FAIL] chem-integrator: the line of the failure with 4 threads is not the one with 1 thread (see OUTPUT/chem-integrator-1.log" \
       "and OUTPUT/chem-integrator-4.log)"
  fail=1
fi
# the same failure on the coarse level of two multigrid levels (steady, INPUT/bc2.txt of 50 cells as for the steam below), where
# the run starts: the line names the grid level
write_fer14 yes chem-integrator-mg "equations = euler
chemistry = finite-rate" "
[MOSE-Chemistry]
ode-relative-tol-species = 1e-30

[MOSE-Multigrid]
levels = 2
level2-iter = 2
level1-iter = 2"
awk '{ if ($5 == 2) $2 = 50; if ($5 >= 3 && $2 > 50) next; print }' INPUT/bc.txt > INPUT/bc2.txt
sed -i -e 's/^cfl = 0.9$/cfl = 0.5/' -e 's/^time-accurate = true$/time-accurate = false/' input.ini
rm -f OUTPUT/field.tec OUTPUT/field-level2.tec
"$LOCAL" > OUTPUT/chem-integrator-mg.log 2>&1
rc=$?
prefix="[ERROR] Chemistry: the integrator H-radau5 failed (return code -1: input is not consistent) in block 1 of grid level 2, cell i, j, k = 1, 1, 1, over a chemistry step of"
n=$(grep -a -c -F -- "$prefix" OUTPUT/chem-integrator-mg.log)
stop=missing; solution=none
grep -a -q -F -- "STOP: chemistry integration error (see above)" OUTPUT/chem-integrator-mg.log && stop=found
[[ -f OUTPUT/field.tec || -f OUTPUT/field-level2.tec ]] && solution=written
if [[ $rc -eq 1 && $n -eq 1 && $stop == found && $solution == none ]]; then
  echo "  [ok  ] chem-integrator-mg, two multigrid levels: exit code 1, the line of the failure once, \"$prefix ...\""
else
  echo "  [FAIL] chem-integrator-mg, two multigrid levels: exit code $rc (required 1), the line \"$prefix ...\" written $n times" \
       "(required once), STOP line $stop, solution $solution (required none) (see OUTPUT/chem-integrator-mg.log)"
  fail=1
fi
rm -f OUTPUT/field.tec OUTPUT/field-level2.tec INPUT/bc2.txt

echo "-- exclude-blocks with chemistry = equilibrium: the excluded block keeps the bits of the run without chemistry --"
# write_steam <chemistry line> <numerics lines> [other sections]: the tube of write_case holding steam far
# from equilibrium (make_case.py steam), with INPUT/bc2.txt for a coarse multigrid level of 50 cells; the files of the
# finite-rate mechanism that the cases above link in INPUT/ are taken out
write_steam() {
  local chemistry="$1" numerics="$2" extra="${3:-}"
  rm -f INPUT/composition.txt INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat INPUT/bc2.txt
  rm -f INPUT/chemistry-info.txt INPUT/chemistry-Arrhenius.dat
  ln -sf "../$CHO/phase.txt" INPUT/phase.txt; ln -sf "../$CHO/thermo.dat" INPUT/thermo.dat
  ln -sf "../$CHO/composition.txt" INPUT/composition.txt
  python3 -B make_case.py steam "$SOD" INPUT/phase.txt || exit 1
  awk '{ if ($5 == 2) $2 = 50; if ($5 >= 3 && $2 > 50) next; print }' INPUT/bc.txt > INPUT/bc2.txt
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 10
newrun = true

[MOSE-Physics]
equations = euler
$chemistry

[MOSE-Numerics]
$numerics
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
$extra
EOI
}
for grid in one multigrid; do
  if [[ $grid == one ]]; then
    label='one grid'
    numerics=$'cfl = 0.9\ntime-accurate = true\ntime-scheme = RK3'; levels=''
  else
    label='two multigrid levels'
    numerics=$'cfl = 0.5\ntime-accurate = false\ntime-scheme = RK3'
    levels=$'\n[MOSE-Multigrid]\nlevels = 2\nlevel2-iter = 5\nlevel1-iter = 10'
  fi
  for chem in none excl all; do
    case $chem in
      none) write_steam "" "$numerics" "$levels" ;;
      excl) write_steam "chemistry = equilibrium" "$numerics" "$levels"$'\n\n[MOSE-Chemistry]\nexclude-blocks = 1' ;;
      all)  write_steam "chemistry = equilibrium" "$numerics" "$levels" ;;
    esac
    rm -f OUTPUT/field.tec "OUTPUT/steam-$grid-$chem.tec"
    "$LOCAL" > "OUTPUT/steam-$grid-$chem.log" 2>&1
    rc=$?
    if [[ $rc -eq 0 && -f OUTPUT/field.tec ]]; then
      mv OUTPUT/field.tec "OUTPUT/steam-$grid-$chem.tec"
    else
      echo "  [FAIL] steam, $label, chemistry $chem: exit code $rc, no field (see OUTPUT/steam-$grid-$chem.log)"
      fail=1
    fi
  done
  if [[ -f OUTPUT/steam-$grid-excl.tec && -f OUTPUT/steam-$grid-none.tec ]] && \
     cmp -s "OUTPUT/steam-$grid-excl.tec" "OUTPUT/steam-$grid-none.tec"; then
    echo "  [ok  ] steam, $label: the block excluded from the equilibrium has the bits of the run without chemistry"
  else
    echo "  [FAIL] steam, $label: with the block excluded the field is not that of the run without chemistry"
    fail=1
  fi
  if [[ -f OUTPUT/steam-$grid-all.tec && -f OUTPUT/steam-$grid-none.tec ]] && \
     ! cmp -s "OUTPUT/steam-$grid-all.tec" "OUTPUT/steam-$grid-none.tec"; then
    echo "  [ok  ] steam, $label: with the block not excluded the equilibrium changes the state"
  else
    echo "  [FAIL] steam, $label: with the block not excluded the field is that of the run without chemistry"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "ExitCodes: FAILED"
  exit 1
fi
echo "ExitCodes: all checks passed"
exit 0
