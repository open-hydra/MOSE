#!/usr/bin/env bash
#===============================================================================
#  input.ini read again during a run (validation tier)
#
#  Every ini-diter iterations MOSE reads input.ini again and takes its values.
#  They must keep to the rules of the registry, as at setup: a value out of
#  its rule stops the run with the message of the rule.
#  A steady viscous run on a 2x2-block box (the case of VtkNames,
#  ../../../mpi/make_case.py), 48 x 48 cells per block, ini-diter = 20, 40
#  iterations, numbered vtk raw solutions every 10 iterations:
#  1. control: a run that starts with cfl = -1 stops at setup with the
#     message "cfl must be > 0";
#  2. a run that starts with cfl = 0.5 gets cfl = -1 in input.ini once its
#     solution of iteration 10 is written: it must stop when it reads
#     input.ini again, at iteration 20, with the same message, instead of
#     going on with a negative CFL number.
#  The change must come before the run writes its solution of iteration 20;
#  if it does not (a host too slow to poll), the test says so and fails.
#  3. the keys of input.ini keep their values when it is read again
#     unchanged: on the same box, a run whose input.ini holds the required
#     keys only (every other key takes its default), and a run whose
#     input.ini also writes integration-variables = prim, riemann-solver =
#     HLLE and flux-limiter = minmod, each read again every 5 iterations,
#     must reach iteration 10, end normally and write the residual history
#     of the same run without readings, byte for byte. A setup that turned
#     some values into names (integration-variables cons into
#     Conservative, riemann-solver HLLC into HLLC Batten, flux-limiter
#     vanleer into Van Leer) would make the check after each reading find
#     them out of the rules, and the run with the required keys only would
#     stop at its first reading.
#  4. a comment line that holds "=" in [MOSE-Numerics], before cfl, is
#     skipped: with cfl = 0.4, the run with the line writes the residual
#     history of the same run without it, byte for byte, and the setup names
#     no key as not read (FiNeR before 2.1.0 ended the section at that line:
#     cfl and the keys after it kept their defaults).
#  5. the rotating frame read again, on the box with 24 x 24 cells per block
#     (runs/frame-case): the box in a rigid-body frame turning
#     about the x axis, with its wall on face 1 of block 1 (the plane x = 0,
#     a surface of revolution about the axis) at rest in the laboratory
#     (stationary-face-1 = 1 1), solutions every 7 iterations, the change
#     written once the solution of iteration 7 is, before the reading of
#     iteration 20 (else the test says so and fails):
#     a. omega from 1e4 to 2e4 rad/s: from iteration 21 the run must go on as
#        one that restarts from its solution of iteration 20 with omega = 2e4
#        from its set-up for 10 iterations: residual-history.dat, torque.dat,
#        the solutions of iterations 21, 28 and 30 and the wall output of
#        iteration 30 byte for byte
#        (the Coriolis and centrifugal sources and the velocity of the
#        stationary wall all at 2e4 rad/s), and the log must give the omega
#        read again;
#     b. rotational-frame = none, and c. stationary-face-1 = 1 3: the run
#        must stop when it reads input.ini again, at iteration 20, with a
#        message that the model, or the stationary faces, cannot change
#        during the run;
#     d. axis = NaN 0.0 0.0: the run must stop when it reads input.ini again,
#        at iteration 20, with the message that a value of axis is not finite
#        (read as a number, the NaN went into every term of the frame);
#     e. two stationary faces at the set-up, stationary-face-1 = 1 1 and
#        stationary-face-2 = 2 2 (the plane x = x_max of block 2), and omega
#        from 1e4 to 2e4: read again with both keys (e0), with
#        stationary-face-1 alone (e1) and with no stationary-face key (e2).
#        The keys left out keep their faces, as every key left out keeps
#        its value: e1 and e2 must reach iteration 30 and write every output
#        file of e0, byte for byte: the keys left out must not stop the run.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../../common/Air-viscous
SOLVER="$HERE/bin/MOSE"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$COMMON/transport.dat" INPUT/transport.dat
python3 -B ../../../mpi/make_case.py --layout 2x2 --cells 48 48 --bc box --zfaces 300 > /dev/null \
  || { echo "RuntimeInput: FAIL — case generation"; exit 1; }
python3 -B ../../../mpi/make_case.py --layout 2x2 --cells 24 24 --bc box --zfaces 300 --outdir runs/frame-case > /dev/null \
  || { echo "RuntimeInput: FAIL — case generation (part 5)"; exit 1; }

write_ini() {  # <directory> <cfl>: input.ini written whole, then moved in place
  cat > "$1/input.ini.new" <<EOF
[MOSE-Parameters]
iter-threshold = 40
res-threshold  = 1e-300
newrun         = true

[MOSE-Physics]
equations = navier-stokes
schmidt   = 1.0

[MOSE-Numerics]
cfl                   = $2
time-accurate         = false
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format    = vtk raw
sol-diter     = 10
sol-overwrite = false
ini-diter     = 20
shell-diter   = 1000000000
EOF
  mv "$1/input.ini.new" "$1/input.ini"
}

prepare() {  # <directory>
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat transport.dat ic.tec bc.txt; do
    ln -s "../../../INPUT/$f" "$1/INPUT/$f"
  done
}

prepare_frame() {  # <directory>: the box of part 5, 24 x 24 cells per block
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat transport.dat; do
    ln -s "../../../INPUT/$f" "$1/INPUT/$f"
  done
  for f in ic.tec bc.txt; do
    ln -s "../../frame-case/INPUT/$f" "$1/INPUT/$f"
  done
}

write_short() {  # <directory> <ini-diter, or empty> [key = value ...]: the required keys and the keys given, 10 iterations
  local dir=$1 diter=$2
  shift 2
  {
    printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 10' ''
    printf '%s\n' '[MOSE-Physics]' 'equations = navier-stokes' ''
    printf '%s\n' '[MOSE-Numerics]' 'time-scheme = RK3' 'cfl = 0.5' 'time-accurate = false' 'space-reconstruction = MUSCL' "$@"
    if [[ -n $diter ]]; then printf '\n%s\n' '[MOSE-IO]' "ini-diter = $diter"; fi
  } > "$dir/input.ini"
}

stopped_on_cfl() {  # <log>: the run stopped with the message of the rule of cfl
  command grep -q 'ERROR\] cfl must be > 0' "$1" && command grep -q -E '(ABORT|STOP): ' "$1" \
    && ! command grep -q -E 'NaN|Program received signal' "$1"
}

status=0
echo "-- 1. cfl = -1 from the start: the run stops at setup --"
d=runs/setup
prepare "$d"; write_ini "$d" -1
( cd "$d" && "$SOLVER" > solve.log 2>&1 )
if stopped_on_cfl "$d/solve.log"; then
  echo "  ok — $(command grep -m1 'ERROR\]' "$d/solve.log" | sed 's/^ *//')"
else
  echo "  FAIL — the run did not stop at setup with the message (see $d/solve.log)"; status=1
fi

echo "-- 2. cfl = -1 written during the run: the run stops when it reads input.ini again --"
d=runs/during
prepare "$d"; write_ini "$d" 0.5
( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
pid=$!
for ((n = 0; n < 6000; n++)); do
  [[ -f $d/OUTPUT/field1.vtm ]] && break
  kill -0 "$pid" 2>/dev/null || break
  sleep 0.05
done
first=no; late=no
[[ -f $d/OUTPUT/field1.vtm ]] && first=yes
[[ -f $d/OUTPUT/field2.vtm ]] && late=yes
write_ini "$d" -1
wait "$pid"
last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
echo "  solution of iteration 10 written before the change: $first; solution of iteration 20 written before it: $late;" \
     "last row of residual-history.dat: ${last:-none}"
if [[ $first != yes ]]; then
  echo "  FAIL — the run ended before its solution of iteration 10 (see $d/solve.log)"; status=1
elif [[ $late == yes ]]; then
  echo "  FAIL — input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"; status=1
elif stopped_on_cfl "$d/solve.log" && [[ $last == 20 ]] && command grep -q 'read again' "$d/solve.log"; then
  echo "  ok — the run stopped at iteration 20: $(command grep -m1 'ERROR\]' "$d/solve.log" | sed 's/^ *//')"
else
  echo "  FAIL — the run did not stop with the message at iteration 20 (see $d/solve.log)"; status=1
fi

echo "-- 3. input.ini read again unchanged every 5 iterations: the keys left out and the keys written keep their values --"
for c in minimal written; do
  keys=()
  [[ $c == written ]] && keys=('integration-variables = prim' 'riemann-solver = HLLE' 'flux-limiter = minmod')
  d=runs/$c-once
  prepare "$d"; write_short "$d" "" ${keys[@]+"${keys[@]}"}
  ( cd "$d" && "$SOLVER" > solve.log 2>&1 )
  d=runs/$c
  prepare "$d"; write_short "$d" 5 ${keys[@]+"${keys[@]}"}
  ( cd "$d" && "$SOLVER" > solve.log 2>&1 )
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
  same=no
  cmp -s "runs/$c-once/OUTPUT/residual-history.dat" "$d/OUTPUT/residual-history.dat" && same=yes
  if [[ $rc -eq 0 && $last == 10 && $same == yes ]] && ! command grep -q -E 'read again|ERROR|ABORT' "$d/solve.log"; then
    echo "  ok — $c keys: the run reached iteration 10 (exit code 0), with the residual history of the run without readings"
  else
    echo "  FAIL — $c keys: the run stopped at iteration ${last:-none} with exit code $rc, residual history as without" \
         "readings: $same: $(command grep -m1 -E 'ERROR\]|must be' "$d/solve.log" | sed 's/^ *//') (see $d/solve.log)"; status=1
  fi
done

echo "-- 4. a comment line with \"=\" before cfl is skipped: cfl is read, no key is named as not read --"
for v in plain comment; do
  d=runs/comment-$v
  prepare "$d"; write_ini "$d" 0.4
  sed -i 's/^iter-threshold = 40$/iter-threshold = 3/' "$d/input.ini"
  if [[ $v == comment ]]; then
    sed -i 's/^\[MOSE-Numerics\]$/[MOSE-Numerics]\n; a slower start: cfl = 0.1/' "$d/input.ini"
  fi
  ( cd "$d" && "$SOLVER" > solve.log 2>&1 )
done
if [[ -s runs/comment-plain/OUTPUT/residual-history.dat ]] \
   && cmp -s runs/comment-plain/OUTPUT/residual-history.dat runs/comment-comment/OUTPUT/residual-history.dat \
   && ! command grep -q 'is not read' runs/comment-comment/solve.log; then
  echo "  ok — with the comment line the run writes the residual history of the run without it, no key named as not read"
else
  echo "  FAIL — with the comment line: residual history not that of the run without it, or a key named as not read:" \
       "$(command grep -m1 'is not read' runs/comment-comment/solve.log | cut -c1-80) (see runs/comment-*/solve.log)"; status=1
fi

echo "-- 5. the rotating frame read again --"
write_frame() {  # <directory> <iter-threshold> <newrun> <model> <omega> <stationary faces> [axis]: input.ini written whole;
                 # the stationary faces separated by commas, written as stationary-face-1, -2, ... (no key if empty)
  local keys="" rest=$6 n=0
  while [[ -n $rest ]]; do
    n=$((n + 1)); keys+="stationary-face-$n = ${rest%%,*}"$'\n'
    [[ $rest == *,* ]] && rest=${rest#*,} || rest=""
  done
  keys=${keys%$'\n'}
  cat > "$1/input.ini.new" <<FRAME
[MOSE-Parameters]
iter-threshold = $2
res-threshold  = 1e-300
newrun         = $3

[MOSE-Physics]
equations        = navier-stokes
schmidt          = 1.0
rotational-frame = $4

[MOSE-Numerics]
cfl                   = 0.5
time-accurate         = false
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format    = vtk raw
sol-diter     = 7
sol-overwrite = false
ini-diter     = 20
shell-diter   = 1000000000

[MOSE-Rotating-Frame]
omega             = $5
axis              = ${7:-1.0 0.0 0.0}
origin            = 0.0 0.0 0.0
${keys}
FRAME
  mv "$1/input.ini.new" "$1/input.ini"
}
during() {  # <directory> <model> <omega> <faces> [axis] [faces at the start]: start at omega = 1e4 with the faces at the
            # start (1 1 if not given), write the change after iteration 7
  local d=$1 pid n
  prepare_frame "$d"; write_frame "$d" 30 true rigid-body 1.0e4 "${6:-1 1}"
  ( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $d/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $d/OUTPUT/field1.vtm ]] && first=yes
  write_frame "$d" 30 true "$2" "$3" "$4" "${5:-1.0 0.0 0.0}"
  # the change came too late if the run had already written its row of iteration 20, after which it reads input.ini
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 20 ]] && late=yes
  wait "$pid"
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
}
during runs/frame-omega rigid-body 2.0e4 '1 1'
d=runs/frame-restart
prepare_frame "$d"; write_frame "$d" 20 true rigid-body 1.0e4 '1 1'
( cd "$d" && "$SOLVER" > solve1.log 2>&1 )
write_frame "$d" 10 false rigid-body 2.0e4 '1 1'
( cd "$d" && "$SOLVER" > solve2.log 2>&1 )
same=yes
for f in residual-history.dat torque.dat field3.vtm field4.vtm field.vtm vtk/field3Block1.vts vtk/field4Block2.vts \
         vtk/fieldBlock1.vts vtk/fieldBlock2.vts vtk/fieldBlock3.vts vtk/fieldBlock4.vts vtk/wallB1F1.vts; do
  if ! cmp -s "runs/frame-omega/OUTPUT/$f" "$d/OUTPUT/$f"; then
    same=no
    echo "  differs: $f"
  fi
done
if [[ $first != yes ]]; then
  echo "  FAIL — a. the run ended before its solution of iteration 7 (see runs/frame-omega/solve.log)"; status=1
elif [[ $late == yes ]]; then
  echo "  FAIL — a. input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"; status=1
elif [[ $last == 30 && $same == yes ]] && \
     command grep -q 'Rotating frame, input.ini read again: omega =  2.0000000E+04' runs/frame-omega/solve.log; then
  echo "  ok — a. omega read again at iteration 20: the run goes on as the restart with omega = 2e4, byte for byte"
else
  echo "  FAIL — a. omega read again: last iteration ${last:-none}, the same files as the restart: $same" \
       "(see runs/frame-omega/solve.log)"; status=1
fi
for c in "model:none:1 1:the model cannot change" "face:rigid-body:1 3:the stationary faces cannot change"; do
  IFS=: read -r tag model face msg <<< "$c"
  during "runs/frame-$tag" "$model" 1.0e4 "$face"
  if [[ $first != yes || $late == yes ]]; then
    echo "  FAIL — $tag: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
  elif [[ $last == 20 ]] && command grep -q "$msg during the run" "runs/frame-$tag/solve.log"; then
    echo "  ok — $tag: the run stopped at iteration 20: $(command grep -m1 'ERROR\] Rotating frame' "runs/frame-$tag/solve.log" | cut -c1-90)…"
  else
    echo "  FAIL — $tag: the run did not stop at iteration 20 with the message (last iteration ${last:-none}," \
         "see runs/frame-$tag/solve.log)"; status=1
  fi
done
during runs/frame-axis rigid-body 1.0e4 '1 1' 'NaN 0.0 0.0'
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — axis: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
elif [[ $last == 20 ]] && \
     command grep -q 'a value of axis in input.ini read again is not finite (NaN or infinity)' runs/frame-axis/solve.log; then
  echo "  ok — axis: the run stopped at iteration 20: $(command grep -m1 'ERROR\] Rotating frame' runs/frame-axis/solve.log | cut -c1-90)…"
else
  echo "  FAIL — axis: the run did not stop at iteration 20 with the message (last iteration ${last:-none}," \
       "see runs/frame-axis/solve.log)"; status=1
fi
# e. two stationary faces at the set-up (the walls x = 0 of block 1 and x = x_max of block 2), omega to 2e4 with both keys
#    (e0), with stationary-face-1 alone (e1) and with no stationary-face key (e2): the keys left out keep their faces
for c in "e0:1 1,2 2" "e1:1 1" "e2:"; do
  IFS=: read -r tag faces <<< "$c"
  during "runs/frame-$tag" rigid-body 2.0e4 "$faces" '1.0 0.0 0.0' '1 1,2 2'
  if [[ $first != yes || $late == yes ]]; then
    echo "  FAIL — $tag: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
  elif [[ $last != 30 ]] || ! command grep -q 'Rotating frame, input.ini read again: omega =  2.0000000E+04' \
                                 "runs/frame-$tag/solve.log"; then
    echo "  FAIL — $tag: the run stopped at iteration ${last:-none}, or the log does not give omega = 2e4 read again" \
         "(see runs/frame-$tag/solve.log)"; status=1
  elif [[ $tag != e0 ]]; then
    list0=$(cd runs/frame-e0/OUTPUT && find . -type f | sort)
    differ=$(cd runs/frame-e0/OUTPUT && find . -type f ! -exec cmp -s {} "../../frame-$tag/OUTPUT/{}" \; -print)
    if [[ $list0 == "$(cd "runs/frame-$tag/OUTPUT" && find . -type f | sort)" && -z $differ ]]; then
      echo "  ok — $tag: stationary-face keys left out (${faces:-none} read again): the run reaches iteration 30 with" \
           "every output file of e0, byte for byte"
    else
      echo "  FAIL — $tag: the output files are not those of e0:" $differ; status=1
    fi
  fi
done

if [[ $status -eq 0 ]]; then
  echo "RuntimeInput: all checks passed"
  exit 0
fi
echo "RuntimeInput: FAILED"
exit 1
