#!/usr/bin/env bash
#===============================================================================
#  Boundary input: the property lines of bc.txt read right or refused (fast)
#
#  A 1-D duct driven by an inlet 402, whose total pressure comes from a table
#  in time named on the property line of the inlet (make_case.py). The same
#  table under four names must give the same run, bit for bit: p0.dat, a name
#  with a directory (INPUT/tables/p0.dat), a name longer than 32 characters
#  and a quoted name with a blank. A table that does not exist, or whose
#  times decrease, must stop the run at set-up, with a non-zero exit code and
#  the name of the file. A table that starts after the start of the run must
#  give the run of the same table with its first value repeated at t = 0.
#  A property line with an empty required value must stop the run at set-up,
#  naming the line of bc.txt, the record, its type and the value; so must a
#  line short of a value, or with a word for a number, and a flow direction
#  that is neither a number nor normal. An SRM grain boundary (502) with the
#  Spalart-Allmaras model runs (40 steps), and so with SST, its ghost cells
#  hold the wall value nu~ = 0 and it is a wall for the wall distance: one
#  explicit step at rest from a uniform rho nu~ changes it in the cell next to
#  the grain by the flux through the grain and the cb2 term of the gradient to
#  the wall value, plus the source terms with the distance to the grain, and
#  in the other cells by those source terms alone (make_case.py srmcheck).
#  Gas: test/common/Air (test/common/Air-viscous for the grain).
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

cat > input.ini <<EOF
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
EOF

fail=0
run() {   # run <tag> <table name>
  python3 -B make_case.py case "$2" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$1.log" 2>&1
  local rc=$?
  [[ -f OUTPUT/field.tec ]] && cp OUTPUT/field.tec "OUTPUT/field-$1.tec"
  return $rc
}

echo "-- the same table of total pressure under four names --"
run ref p0.dat
if [[ ! -f OUTPUT/field-ref.tec ]] || grep -a -qE 'NaN|Integration failed' OUTPUT/ref.log; then
  echo "  [FAIL] the reference run (p0.dat) did not complete (see OUTPUT/ref.log)"
  echo "BCInput: FAILED"
  exit 1
fi
for tag_name in "dir:INPUT/tables/p0.dat" "long:INPUT/a_directory_with_a_long_name/inlet_total_pressure.dat" \
                "blank:INPUT/p0 table.dat"; do
  tag=${tag_name%%:*}
  name=${tag_name#*:}
  run "$tag" "$name"
  if [[ -f "OUTPUT/field-$tag.tec" ]] && python3 -B make_case.py same OUTPUT/field-ref.tec "OUTPUT/field-$tag.tec"; then
    echo "  [ok  ] $name: the run of p0.dat, bit for bit"
  else
    echo "  [FAIL] $name: not the run of p0.dat (see OUTPUT/$tag.log)"
    fail=1
  fi
done

echo "-- a table that does not exist --"
run missing missing.dat
rc=$?
if [[ $rc -ne 0 ]] && grep -a -q 'missing.dat' OUTPUT/missing.log; then
  echo "  [ok  ] stopped at set-up, exit code $rc, file named"
else
  echo "  [FAIL] exit code $rc, the missing table was not reported (see OUTPUT/missing.log)"
  fail=1
fi

echo "-- a table that starts after the start of the run --"
run late late.dat
run early early.dat
if [[ -f OUTPUT/field-late.tec && -f OUTPUT/field-early.tec ]] && \
   python3 -B make_case.py same OUTPUT/field-early.tec OUTPUT/field-late.tec; then
  echo "  [ok  ] the first value before the first time: the run of early.dat, bit for bit"
else
  echo "  [FAIL] late.dat is not the run of early.dat (see OUTPUT/late.log)"
  fail=1
fi

echo "-- a table whose times decrease --"
run decreasing decreasing.dat
rc=$?
if [[ $rc -ne 0 ]] && grep -a -q 'decreasing.dat' OUTPUT/decreasing.log; then
  echo "  [ok  ] stopped at set-up, exit code $rc, file named"
else
  echo "  [FAIL] exit code $rc, the decreasing times were not reported (see OUTPUT/decreasing.log)"
  fail=1
fi

echo "-- property lines that cannot be read as the values of their type --"
bad() {   # bad <kind> <pattern>: the run must stop with a non-zero code and a message matching the pattern
  python3 -B make_case.py bad "$1" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/bad-$1.log" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]] && grep -a -q -E "$2" "OUTPUT/bad-$1.log"; then
    echo "  [ok  ] $1: stopped at set-up, exit code $rc"
  else
    echo "  [FAIL] $1: exit code $rc, the line was not refused as it should (see OUTPUT/bad-$1.log)"
    fail=1
  fi
}
bad empty 'bc.txt line 4 \(record 2, type 406\): the value of rel_fac is empty'
bad short 'bc.txt line 4 \(record 2, type 406\): the line is short or holds a value that is not a number'
bad word 'bc.txt line 4 \(record 2, type 406\): the line is short or holds a value that is not a number'
bad angle "bc.txt line 2 \\(record 1, type 402\\): alpha is neither a number nor normal: 'nromal'"

echo "-- an SRM grain boundary (502) with a RANS model --"
python3 -B make_case.py srm || exit 1
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../../../../common/Air-viscous/phase.txt"     INPUT/phase.txt
ln -sf "../../../../common/Air-viscous/thermo.dat"    INPUT/thermo.dat
ln -sf "../../../../common/Air-viscous/transport.dat" INPUT/transport.dat
sed -e 's/^equations = euler$/equations = navier-stokes\nturbulence = SA\nschmidt = 0.7/' input.ini > input-srm.ini && mv input-srm.ini input.ini
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/srm-rans.log 2>&1
rc=$?
if [[ $rc -eq 0 && -f OUTPUT/field.tec ]] && ! grep -a -qE 'NaN|Integration failed' OUTPUT/srm-rans.log; then
  echo "  [ok  ] srm-rans: 40 steps, exit code 0"
else
  echo "  [FAIL] srm-rans: exit code $rc, the run with the grain and Spalart-Allmaras did not complete (see OUTPUT/srm-rans.log)"
  fail=1
fi
python3 -B make_case.py srm sst || exit 1
sed -e 's/^turbulence = SA$/turbulence = SST/' input.ini > input-sst.ini
mv input.ini input-sa.ini && mv input-sst.ini input.ini
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/srm-sst.log 2>&1
rc=$?
if [[ $rc -eq 0 && -f OUTPUT/field.tec ]] && ! grep -a -qE 'NaN|Integration failed' OUTPUT/srm-sst.log && \
   grep -a -q 'Turbulence variables below zero after a stage (set to the floor): 0$' OUTPUT/srm-sst.log; then
  echo "  [ok  ] srm-sst: 40 steps with SST, exit code 0, k and omega never below zero"
else
  echo "  [FAIL] srm-sst: exit code $rc, the run with the grain and SST did not complete as it should (see OUTPUT/srm-sst.log)"
  fail=1
fi
mv input-sa.ini input.ini
python3 -B make_case.py srm rest || exit 1
sed -e 's/^iter-threshold = 40$/iter-threshold = 1/' -e 's/^time-accurate = true$/time-accurate = false/' \
    -e 's/^time-scheme = RK3$/time-scheme = euler/' -e 's/^schmidt = 0.7$/schmidt = 0.7\n\n[MOSE-Turbulence]\npoint-implicit = false/' \
    input.ini > input-rest.ini && mv input-rest.ini input.ini
rm -f OUTPUT/field.tec OUTPUT/diagnostic.tec
"$LOCAL" > OUTPUT/srm-rest.log 2>&1
if [[ -f OUTPUT/diagnostic.tec ]] && python3 -B make_case.py srmcheck OUTPUT/diagnostic.tec; then
  echo "  [ok  ] srm-rest: the ghost cells of the grain hold the wall value nu~ = 0, and the grain is a wall for d"
else
  echo "  [FAIL] srm-rest: the wall value nu~ = 0 or the distance to the grain missing (see OUTPUT/srm-rest.log)"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "BCInput: FAILED"
  exit 1
fi
echo "BCInput: all checks passed"
exit 0
