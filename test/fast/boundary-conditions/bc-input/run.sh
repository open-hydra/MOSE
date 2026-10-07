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
#  line short of a value, or with a word for a number.
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

if [[ $fail -ne 0 ]]; then
  echo "BCInput: FAILED"
  exit 1
fi
echo "BCInput: all checks passed"
exit 0
