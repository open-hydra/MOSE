#!/usr/bin/env bash
#===============================================================================
#  input.ini read again during an MPI run: every rank takes the values of the
#  file that rank 0 reads (needs-mpi tier)
#
#  Every ini-diter iterations MOSE reads input.ini again. Rank 0 reads it and
#  passes its text to every rank, which parses the same text: the ranks take
#  the same values. Some values decide collective calls that each rank makes
#  on its own: sol-diter (the solution output, a gather to rank 0),
#  res-diter (the reduction of the residuals) and ini-diter (the reading
#  itself). With different values the ranks would call their collectives in
#  different orders, which MPI forbids (MPI-4.1, section 6.14), and the run
#  would hang: so it did when every rank read input.ini on its own and a file
#  rewritten during the run gave the ranks different contents.
#
#  Two ranks launched as two programs (MPMD), each in its own folder with its
#  own input.ini (runs/<case>/A for rank 0, runs/<case>/B for rank 1), on the
#  2x2-block viscous box of RuntimeInput (../make_case.py, 48 x 48 cells per
#  block), 40 iterations, solutions every 10 iterations, input.ini read again
#  every 20 iterations. Once the solution of iteration 10 is written, the
#  input.ini of rank 0 alone gets
#    sol-rank0  sol-diter = 7;
#  every rank must take it: the run must end within 60 s with exit code 0 at
#  iteration 40, and the numbered solutions in the OUTPUT folder of rank 0
#  must be those of iterations 10, 20, 21, 28 and 35 (the CYCLE of their
#  block files): the reading of iteration 20 comes after the solution of
#  iteration 20, and a rank that kept sol-diter = 10 would leave the run
#  hanging at iteration 21. Then, in the same way, the input.ini of rank 1
#  alone gets
#    sol  sol-diter = 7,
#    res  res-diter = 7,
#    ini  ini-diter = 30,
#    cfl  cfl = -1 (out of its rule);
#  rank 1 must not take it: each run must end within 60 s with exit code 0 at
#  iteration 40, and the files in the OUTPUT folder of rank 0 must be those of
#  the same run with no change, byte for byte (ref). The change must be in
#  place before the run writes its solution of iteration 20, which comes
#  before the reading of iteration 20; if it is not (a host too slow to poll),
#  the test says so and fails. The cases run one after the other and the test
#  stops at the first that fails, so that it ends within its TIMEOUT of 120 s.
#  The MOSE processes of a run that the launcher leaves behind (this copy of
#  the solver, working in runs/<case>/A or runs/<case>/B) are named, stopped,
#  and fail the test.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
SOLVER="$(cd bin && pwd -P)/MOSE"

# write_box_ini <folder> [<key> <value>]: the input.ini of the cfl case of RuntimeInput, with one key changed, written
# whole, then moved in place
write_box_ini() {
  local cfl=0.5 sol=10 ini=20 res=1
  case ${2:-} in
    cfl) cfl=$3 ;;
    sol-diter) sol=$3 ;;
    ini-diter) ini=$3 ;;
    res-diter) res=$3 ;;
  esac
  cat > "$1/input.ini.new" <<EOF
[MOSE-Parameters]
iter-threshold = 40
res-threshold  = 1e-300
newrun         = true

[MOSE-Physics]
equations = navier-stokes
schmidt   = 1.0

[MOSE-Numerics]
cfl                   = $cfl
time-accurate         = false
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format    = vtk raw
sol-diter     = $sol
sol-overwrite = false
res-diter     = $res
ini-diter     = $ini
shell-diter   = 1000000000
EOF
  mv "$1/input.ini.new" "$1/input.ini"
}

rm -rf runs INPUT
mkdir -p INPUT runs
for f in phase.txt thermo.dat transport.dat; do
  ln -s "../../../common/Air-viscous/$f" "INPUT/$f"
done
python3 -B ../make_case.py --layout 2x2 --cells 48 48 --bc box --zfaces 300 > /dev/null || exit 1
RUNS="$(cd runs && pwd -P)"

left_behind() {  # <case>: the processes of this copy of the solver whose working folder is one of the two of the case
  local p
  for p in /proc/[0-9]*; do
    [[ $(readlink "$p/exe" 2>/dev/null) == "$SOLVER" ]] || continue
    case $(readlink "$p/cwd" 2>/dev/null) in
      "$RUNS/$1/A" | "$RUNS/$1/B") echo "${p#/proc/}" ;;
    esac
  done
}

# cycles <OUTPUT folder>: the CYCLE of each numbered solution, from the first block file of each (vtk raw: the field
# data CYCLE, an Int32 in the appended data, after its header of 4 bytes, UInt32, or 8, UInt64)
cycles() {
  python3 -B - "$1" <<'PY'
import re, struct, sys
out, n = [], 1
while True:
    try:
        b = open('%s/vtk/field%dBlock1.vts' % (sys.argv[1], n), 'rb').read()
    except OSError:
        break
    hs = 8 if b'header_type="UInt64"' in b[:400] else 4
    off = int(re.search(rb'Name="CYCLE" format="appended" offset="(\d+)"', b).group(1))
    a = b.index(b'_', b.index(b'<AppendedData encoding="raw">')) + 1
    out.append(str(struct.unpack('<i', b[a + off + hs:a + off + hs + 4])[0]))
    n += 1
print(' '.join(out))
PY
}

fail=0
# run_case <case> [<folder> <key> <value>]: the run of two folders, the change written in the input.ini of one folder
# alone (A of rank 0, B of rank 1) once the solution of iteration 10 is written
run_case() {
  local c=$1 A=$RUNS/$1/A B=$RUNS/$1/B d f n pid rc first late last left cyc
  for d in "$A" "$B"; do
    mkdir -p "$d/INPUT" "$d/OUTPUT"
    for f in phase.txt thermo.dat transport.dat ic.tec bc.txt; do
      ln -s "../../../../INPUT/$f" "$d/INPUT/$f"   # runs/<case>/A/INPUT: four levels below this folder
    done
    write_box_ini "$d"
  done
  timeout -k 10 60 "$MPIEXEC" "$NPFLAG" 1 -wdir "$A" "$SOLVER" : "$NPFLAG" 1 -wdir "$B" "$SOLVER" \
    > "runs/$c.log" 2>&1 &
  pid=$!
  if [[ $# -gt 1 ]]; then
    for ((n = 0; n < 6000; n++)); do
      [[ -f $A/OUTPUT/field1.vtm ]] && break
      kill -0 "$pid" 2>/dev/null || break
      sleep 0.05
    done
    first=no; late=no
    [[ -f $A/OUTPUT/field1.vtm ]] && first=yes
    write_box_ini "$RUNS/$c/$2" "$3" "$4"
    # MOSE reads input.ini again after it writes the solution of the same iteration (MOSE_postprocess): without the
    # solution of iteration 20 once the change is in place, the reading of iteration 20 comes after the change
    [[ -f $A/OUTPUT/field2.vtm ]] && late=yes
  fi
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$A/OUTPUT/residual-history.dat" 2>/dev/null)
  left=$(left_behind "$c" | paste -sd ' ' -)
  if [[ -n $left ]]; then
    kill -TERM $left 2>/dev/null
    sleep 2
    kill -KILL $left 2>/dev/null
    echo "  [FAIL] $c, MOSE processes of this run left behind by the launcher, now stopped: $left"
    fail=1
  fi
  if [[ $# -gt 1 && $first != yes ]]; then
    echo "  [FAIL] $c, the run ended before its solution of iteration 10, exit code $rc (see runs/$c.log)"
    fail=1
    return
  fi
  if [[ $# -gt 1 && $late == yes ]]; then
    echo "  [FAIL] $c, input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"
    fail=1
    return
  fi
  if [[ $rc -ne 0 || ${last:-} != 40 ]]; then
    echo "  [FAIL] $c, exit code $rc (required 0 within 60 s; 124 is the time limit), last row of the residual history" \
         "${last:-none} (required 40) (see runs/$c.log)"
    fail=1
    return
  fi
  if [[ $c == ref ]]; then
    echo "  [ok  ] ref, no change: exit code 0 at iteration 40"
  elif [[ $2 == A ]]; then
    cyc=$(cycles "$A/OUTPUT")
    if [[ $cyc != "10 20 21 28 35" || ! -f $A/OUTPUT/field.vtm ]]; then
      echo "  [FAIL] $c, $3 = $4 in the input.ini of rank 0 alone: numbered solutions of iterations ${cyc:-none}" \
           "(required 10 20 21 28 35), final solution field.vtm: $([[ -f $A/OUTPUT/field.vtm ]] && echo yes || echo no)"
      fail=1
      return
    fi
    echo "  [ok  ] $c, $3 = $4 in the input.ini of rank 0 alone: exit code 0 at iteration 40, numbered solutions of" \
         "iterations $cyc: every rank took it"
  else
    if ! diff <(cd "$RUNS/ref/A/OUTPUT" && find . -type f | sort) <(cd "$A/OUTPUT" && find . -type f | sort) > /dev/null; then
      echo "  [FAIL] $c, the files in the OUTPUT folder of rank 0 differ in their list from those of ref"
      fail=1
      return
    fi
    while IFS= read -r f; do
      if ! cmp -s "$RUNS/ref/A/OUTPUT/$f" "$A/OUTPUT/$f"; then
        echo "  [FAIL] $c, ${f#./} differs from that of ref"
        fail=1
        return
      fi
    done < <(cd "$RUNS/ref/A/OUTPUT" && find . -type f | sort)
    echo "  [ok  ] $c, $3 = $4 in the input.ini of rank 1 alone: exit code 0 at iteration 40, the files of ref byte for byte"
  fi
}

echo "-- input.ini read again: every rank takes the values of the file of rank 0 (2 ranks, two folders) --"
run_case ref
# One case after the other, up to the first that fails: a run that hangs lasts up to 70 s (its time limit)
for args in "sol-rank0 A sol-diter 7" "sol B sol-diter 7" "res B res-diter 7" "ini B ini-diter 30" "cfl B cfl -1"; do
  [[ $fail -eq 0 ]] || break
  run_case $args
done

if [[ $fail -eq 0 ]]; then
  echo "MPIRuntimeInput: all checks passed"
  exit 0
fi
echo "MPIRuntimeInput: FAILED"
exit 1
