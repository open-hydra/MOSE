#!/usr/bin/env bash
#===============================================================================
#  input.ini read again during a run (validation tier)
#
#  Every ini-diter iterations MOSE reads input.ini again and takes its values.
#  They must keep to the rules of the registry, as at setup. An error of the
#  reading (a value out of its rule or that does not read, a file that cannot
#  be read, another value of a key fixed for the run, of levels or more
#  probes, a change of the model or of the stationary faces of the rotating
#  frame, a value of its axis that is not finite or does not read, residual
#  smoothing switched on in a run started without it, values that a check of
#  the set-up between keys refuses) does what ini-error says: with stop the
#  run stops with the message, as at setup;
#  with warn, the default, the reading is not taken: a warning gives the
#  error, and the run goes on with every value it had before the reading.
#  A steady viscous run on a 2x2-block box (the case of VtkNames,
#  ../../../mpi/make_case.py), 48 x 48 cells per block, ini-diter = 20, 40
#  iterations, numbered vtk raw solutions every 10 iterations:
#  1. control: a run that starts with cfl = -1 stops at setup with the
#     message "cfl must be > 0";
#  2. a run that starts with cfl = 0.5 and ini-error = stop gets cfl = -1 in
#     input.ini once its solution of iteration 10 is written: it must stop
#     when it reads input.ini again, at iteration 20, with the same message,
#     instead of going on with a negative CFL number.
#     The box with 24 x 24 cells per block of part 5 (runs/frame-case),
#     without the rotating frame, 30 iterations, solutions every 7, the
#     change written once the solution of iteration 7 is: with ini-error left
#     out (warn), a. cfl = -1, b. input.ini emptied, c. cfl = 0.4 with
#     sol-diter = seven (a value that does not read), e. cfl = -1 with
#     ini-error = stop, f. input.ini of comment lines only (the INI reader
#     finds no section): the run must reach iteration 30 with exit code 0,
#     write once that the reading of iteration 20 is not taken, with its
#     error, and no [ERROR], and write the residual history of the same run
#     without the change, byte for byte (with c, cfl = 0.4 taken at iteration
#     20 would change it; with e, the value of ini-error in force when the
#     reading starts decides). With ini-error = stop from the start, b.
#     input.ini emptied and e. cfl = -1 with ini-error = warn must stop the
#     run at iteration 20 with the error.
#  The change must come before the run reads input.ini again at iteration
#  20; if it does not (a host too slow to poll), the test says so and fails.
#  3. the keys of input.ini keep their values when it is read again
#     unchanged: on the same box, a run whose input.ini holds the required
#     keys only (every other key takes its default), and a run whose
#     input.ini also writes integration-variables = prim, riemann-solver =
#     HLLE and flux-limiter = minmod, each read again every 5 iterations,
#     must reach iteration 10, end normally and write the residual history
#     of the same run without readings, byte for byte. A setup that turns
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
#     b. rotational-frame = none, and c. stationary-face-1 = 1 3: with
#        ini-error = stop the run must stop when it reads input.ini again,
#        at iteration 20, with a message that the model, or the stationary
#        faces, cannot change during the run;
#     d. axis = NaN 0.0 0.0: with ini-error = stop the run must stop when it
#        reads input.ini again, at iteration 20, with the message that a
#        value of axis is not finite (read as a number, the NaN would go into
#        every term of the frame);
#     b, c, d and g. axis = 1.0 abc 0.0 (a value that does not read) with
#        ini-error left out (warn): the reading is not taken: the run must
#        reach iteration 30 with exit code 0, write once that the reading of
#        iteration 20 is not taken, with the error, and no [ERROR], and write
#        every output file of the frame run without the change, byte for
#        byte;
#     e. two stationary faces at the set-up, stationary-face-1 = 1 1 and
#        stationary-face-2 = 2 2 (the plane x = x_max of block 2), and omega
#        from 1e4 to 2e4: read again with both keys (e0), with
#        stationary-face-1 alone (e1) and with no stationary-face key (e2).
#        The keys left out keep their faces, as every key left out keeps
#        its value: e1 and e2 must reach iteration 30 and write every output
#        file of e0, byte for byte: the keys left out must not stop the run.
#  6. the class of a key (the column During the run of the registry: taken,
#     fixed for the run, start only), on the box of part 2 with 24 x 24
#     cells per block, 30 iterations, solutions every 7, the change written
#     once the solution of iteration 7 is, before the reading of iteration
#     20 (else the test says so and fails). A reading that writes for a
#     fixed key another value than the start read is an error of the
#     reading:
#     a. time-scheme = RK2 in a run started with RK3, ini-error = stop: the
#        run must stop at the reading of iteration 20 with the message that
#        time-scheme cannot change during the run;
#     b. time-scheme = RK2 and cfl = 0.4 in the same reading, ini-error left
#        out (warn): the reading is not taken, so cfl = 0.4 is not taken
#        either: iteration 30, exit code 0, one warning with the error, no
#        [ERROR], the residual history of the run without the change, byte
#        for byte;
#     c. levels = 2 added in [MOSE-Multigrid], and d. a probe added, with
#        warn: the grid levels and the probes are those of the start, so
#        the reading is not taken, as in b;
#     e. newrun, a key of the start only: newrun = false written in a run
#        in the rotating frame of part 5 with res-diter = 25, whose
#        OUTPUT/torque.dat holds a row of an earlier run before the start:
#        the reading is taken, with no warning, and the first row of the
#        torque, at iteration 25, starts the file again with its header,
#        as newrun = true of the start says (the row of the earlier run is
#        gone);
#     f. ode-relative-tol-temperature, a fixed key, 1e-5 at the start and
#        1.0e-05 in the reading: the same value, so the reading is taken,
#        with no warning, and the run writes the residual history of the
#        run without the change, byte for byte;
#     g. time-accurate = true written during the steady run, ini-error
#        left out (warn): the reading is not taken (a run that took it
#        would advance its time from -1 s), as in b.
#  7. implicit residual smoothing (IRS), on with irs = true or irs-beta > 0
#     at the start and at every reading, on the box of part 6:
#     a. irs-beta = 0.5 written during a run started without IRS, with
#        ini-error = stop: its arrays exist only when it is on at the
#        start, so the run must stop at the reading of iteration 20 with
#        the message that it cannot be switched on during the run;
#     b. irs = false and irs-beta = 0.5 written together in the input.ini
#        of the start, on the box of part 3, read again unchanged every 5
#        iterations: IRS stays on, and the run writes the residual history
#        of the same run without readings, byte for byte;
#     c. irs-beta = 0 written alone during a run started with
#        irs-beta = 0.5 (irs left out): IRS is off from the reading of
#        iteration 20, so from iteration 21 the run must go on as one that
#        restarts from its solution of iteration 20 with irs-beta = 0:
#        residual-history.dat and the solutions of iterations 21, 28 and
#        30, byte for byte.
#  8. the checks of the set-up between keys, made again on the values read,
#     with ini-error = stop: the run must stop at the reading of iteration
#     20 with the key, its value and the reason of the set-up:
#     a. riemann-solver = AUSM+M, a low-Mach solver, with
#        riemann-options-Mco = 0.1 at the start (the box of part 6), and
#        riemann-options-Mco = 0 written during the run;
#     b. the gas of four species of test/1D/diff1, without
#        INPUT/diffusion.dat, in its time-accurate tube of 200 cells, with
#        schmidt = 1.0 at the start and schmidt = 0, multicomponent
#        diffusion, which needs that table, written once the run has
#        written the row of iteration 7 of its residual history; an
#        iteration takes about a millisecond (gfortran TESTING build, one
#        thread), so input.ini is read again at iteration 1000, of 1010,
#        where the run must stop.
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

write_ini() {  # <directory> <cfl> [ini-error] [sol-diter] [iter-threshold]: input.ini written whole, then moved in place
  {
  cat <<EOF
[MOSE-Parameters]
iter-threshold = ${5:-40}
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
sol-diter     = ${4:-10}
sol-overwrite = false
ini-diter     = 20
shell-diter   = 1000000000
EOF
  if [[ -n ${3:-} ]]; then echo "ini-error     = $3"; fi
  } > "$1/input.ini.new"
  mv "$1/input.ini.new" "$1/input.ini"
}

prepare() {  # <directory>
  rm -rf "$1"; mkdir -p "$1/INPUT" "$1/OUTPUT"
  for f in phase.txt thermo.dat transport.dat ic.tec bc.txt; do
    ln -s "../../../INPUT/$f" "$1/INPUT/$f"
  done
}

prepare_frame() {  # <directory>: the box of 24 x 24 cells per block (runs/frame-case), of part 2 and part 5
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

echo "-- 2. cfl = -1 written during the run, ini-error = stop: the run stops when it reads input.ini again --"
d=runs/during
prepare "$d"; write_ini "$d" 0.5 stop
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
write_ini "$d" -1 stop
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

# change <directory> <ini-error at the start, or empty> <change>: the box of 24 x 24 cells per block with cfl = 0.5, 30
# iterations, solutions every 7, input.ini read again at iteration 20, the change written whole once the solution of iteration
# 7 is written; sets first, late, rc and last
change() {
  local d=$1 err=$2 how=$3 pid n
  prepare_frame "$d"; write_ini "$d" 0.5 "$err" 7 30
  ( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $d/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $d/OUTPUT/field1.vtm ]] && first=yes
  case $how in
    cfl-1)        write_ini "$d" -1 "$err" 7 30 ;;
    cfl-1-stop)   write_ini "$d" -1 stop 7 30 ;;
    cfl-1-warn)   write_ini "$d" -1 warn 7 30 ;;
    cfl0.4-seven) write_ini "$d" 0.4 "$err" seven 30 ;;
    empty)        : > "$d/input.ini.new"; mv "$d/input.ini.new" "$d/input.ini" ;;
    comments)     printf '%s\n' '; input.ini of comment lines only' '; cfl = 0.4' > "$d/input.ini.new"
                  mv "$d/input.ini.new" "$d/input.ini" ;;
  esac
  # the change came too late if the run had already written its row of iteration 20, after which it reads input.ini
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 20 ]] && late=yes
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
}

# not_taken <directory> <case> <error>: iteration 30, exit code 0, the warning of the reading of iteration 20 once with the
# error, no [ERROR] in the log, and the residual history of the run without the change, byte for byte
not_taken() {
  local d=$1 tag=$2 err=$3 nw ne same=no
  nw=$(command grep -a -c -F "[WARNING] input.ini read again at iteration 20 is not taken: $err" "$d/solve.log")
  ne=$(command grep -a -c -F '[ERROR]' "$d/solve.log")
  cmp -s runs/box-ref/OUTPUT/residual-history.dat "$d/OUTPUT/residual-history.dat" && same=yes
  if [[ $first != yes ]]; then
    echo "  FAIL — $tag: the run ended before its solution of iteration 7 (see $d/solve.log)"; status=1
  elif [[ $late == yes ]]; then
    echo "  FAIL — $tag: input.ini changed after the run had passed iteration 20: the test cannot tell (host too slow to poll)"; status=1
  elif [[ $rc -eq 0 && $last == 30 && $nw -eq 1 && $ne -eq 0 && $same == yes ]]; then
    echo "  ok — $tag: the reading of iteration 20 is not taken ($err), the run goes on with the values before it"
  else
    echo "  FAIL — $tag: exit code $rc, last iteration ${last:-none} (required 30), warnings $nw (required 1), [ERROR] lines" \
         "$ne, residual history of the run without the change: $same (see $d/solve.log)"; status=1
  fi
}

# stopped <directory> <case> <error>: the run stopped at the reading of iteration 20 with the error
stopped() {
  local d=$1 tag=$2 err=$3
  if [[ $first != yes || $late == yes ]]; then
    echo "  FAIL — $tag: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
  elif [[ $rc -ne 0 && $last == 20 ]] && command grep -a -q -F "$err" "$d/solve.log" \
       && command grep -a -q -E '(ABORT|STOP): ' "$d/solve.log"; then
    echo "  ok — $tag: the run stopped at iteration 20: $err"
  else
    echo "  FAIL — $tag: exit code $rc, last iteration ${last:-none} (required 20), error \"$err\" (see $d/solve.log)"; status=1
  fi
}

echo "-- 2. errors of input.ini read again, ini-error left out (warn): the reading is not taken --"
d=runs/box-ref
prepare_frame "$d"; write_ini "$d" 0.5 "" 7 30
( cd "$d" && "$SOLVER" > solve.log 2>&1 ) || { echo "  FAIL — the run without a change failed (see $d/solve.log)"; status=1; }
change runs/warn-cfl "" cfl-1
not_taken runs/warn-cfl "a. cfl = -1" "cfl must be > 0"
change runs/warn-empty "" empty
not_taken runs/warn-empty "b. input.ini emptied" "input.ini is empty"
change runs/warn-whole "" cfl0.4-seven
not_taken runs/warn-whole "c. cfl = 0.4 with sol-diter = seven" \
  "input.ini, [MOSE-IO] sol-diter = seven: the value is not an integer"
change runs/warn-key "" cfl-1-stop
not_taken runs/warn-key "e. cfl = -1 with ini-error = stop" "cfl must be > 0"
change runs/warn-comments "" comments
not_taken runs/warn-comments "f. comment lines only" "input.ini cannot be read by the INI reader"
echo "-- 2. errors of input.ini read again, ini-error = stop from the start: the run stops --"
change runs/stop-empty stop empty
stopped runs/stop-empty "b. input.ini emptied" "input.ini is empty"
change runs/stop-key stop cfl-1-warn
stopped runs/stop-key "e. cfl = -1 with ini-error = warn" "cfl must be > 0"

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
                 # the stationary faces separated by commas, written as stationary-face-1, -2, ... (no key if empty);
                 # ini-error = $INIERR when INIERR is not empty
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
${INIERR:+ini-error     = $INIERR}

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
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
}
INIERR=
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
INIERR=stop
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
INIERR=
# b, c, d and g with ini-error left out (warn): the reading of iteration 20 is not taken; the frame run without a change
d=runs/frame-ref
prepare_frame "$d"; write_frame "$d" 30 true rigid-body 1.0e4 '1 1'
( cd "$d" && "$SOLVER" > solve.log 2>&1 ) || { echo "  FAIL — the frame run without a change failed (see $d/solve.log)"; status=1; }
list0=$(cd runs/frame-ref/OUTPUT && find . -type f | sort)
for c in "model:none:1 1:1.0 0.0 0.0:the model cannot change during the run" \
         "face:rigid-body:1 3:1.0 0.0 0.0:the stationary faces cannot change during the run" \
         "axis:rigid-body:1 1:NaN 0.0 0.0:a value of axis in input.ini read again is not finite (NaN or infinity)" \
         "axisword:rigid-body:1 1:1.0 abc 0.0:axis = 1.0 abc 0.0 in input.ini read again does not read as three numbers"; do
  IFS=: read -r tag model face axis msg <<< "$c"
  during "runs/frame-warn-$tag" "$model" 1.0e4 "$face" "$axis"
  d=runs/frame-warn-$tag
  nw=$(command grep -a -c -F "is not taken: Rotating frame: " "$d/solve.log")
  ne=$(command grep -a -c -F '[ERROR]' "$d/solve.log")
  differ=$(cd runs/frame-ref/OUTPUT && find . -type f ! -exec cmp -s {} "../../frame-warn-$tag/OUTPUT/{}" \; -print)
  if [[ $first != yes || $late == yes ]]; then
    echo "  FAIL — $tag, warn: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
  elif [[ $rc -eq 0 && $last == 30 && $nw -eq 1 && $ne -eq 0 && -z $differ ]] \
       && [[ $list0 == "$(cd "$d/OUTPUT" && find . -type f | sort)" ]] && command grep -a -q -F "$msg" "$d/solve.log"; then
    echo "  ok — $tag, warn: the reading of iteration 20 is not taken ($msg), every output file of the run without a change"
  else
    echo "  FAIL — $tag, warn: exit code $rc, last iteration ${last:-none} (required 30), warnings $nw (required 1), [ERROR]" \
         "lines $ne, files other than those of the run without a change:" $differ "(see $d/solve.log)"; status=1
  fi
done
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

echo "-- 6. the class of a key: what a reading takes --"
# rewrite <directory> <sed script>: input.ini written again whole with the sed script applied, then moved in place
rewrite() {
  sed -e "$2" "$1/input.ini" > "$1/input.ini.new"
  mv "$1/input.ini.new" "$1/input.ini"
}
# extend <directory> <line>...: input.ini written again whole with the lines added at its end, then moved in place
extend() {
  local d=$1
  shift
  { cat "$d/input.ini"; printf '%s\n' "$@"; } > "$d/input.ini.new"
  mv "$d/input.ini.new" "$d/input.ini"
}
to_rk2()      { rewrite "$1" 's/^time-scheme           = RK3$/time-scheme           = RK2/'; }
to_rk2_cfl()  { rewrite "$1" 's/^time-scheme           = RK3$/time-scheme           = RK2/; s/^cfl                   = 0.5$/cfl                   = 0.4/'; }
add_levels()  { extend "$1" '' '[MOSE-Multigrid]' 'levels = 2'; }
add_probe()   { extend "$1" '' '[MOSE-Probes]' 'probe1 = p1' '' '[p1]' 'variables = p' 'index-position = 1 1 1 1'; }
tol_start()   { extend "$1" '' '[MOSE-Chemistry]' 'ode-relative-tol-temperature = 1e-5'; }
tol_written() { rewrite "$1" 's/^ode-relative-tol-temperature = 1e-5$/ode-relative-tol-temperature = 1.0e-05/'; }
# reread <directory> <ini-error at the start, or empty> <change> [<start>]: the box of part 2 with cfl = 0.5, 30 iterations,
# solutions every 7, input.ini read again at iteration 20, written by write_ini and then, if given, by the function <start>; once
# the solution of iteration 7 is written, the function <change> writes input.ini again; sets first, late, rc and last
reread() {
  local d=$1 err=$2 how=$3 start=${4:-} pid n
  prepare_frame "$d"; write_ini "$d" 0.5 "$err" 7 30
  [[ -n $start ]] && "$start" "$d"
  ( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $d/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $d/OUTPUT/field1.vtm ]] && first=yes
  "$how" "$d"
  # the change came too late if the run had already written its row of iteration 20, after which it reads input.ini
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 20 ]] && late=yes
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
}
fixed='input.ini, [MOSE-Numerics] time-scheme = RK2: the run has time-scheme = RK3 from its set-up, and time-scheme cannot change during the run'
reread runs/fixed-stop stop to_rk2
stopped runs/fixed-stop "a. time-scheme = RK2, ini-error = stop" "$fixed"
reread runs/fixed-warn "" to_rk2_cfl
not_taken runs/fixed-warn "b. time-scheme = RK2 and cfl = 0.4" "$fixed"
reread runs/fixed-levels "" add_levels
not_taken runs/fixed-levels "c. levels = 2 added" \
  'input.ini, [MOSE-Multigrid] levels = 2: the run started without levels in input.ini (its default: 1), and levels cannot change during the run'
reread runs/fixed-probe "" add_probe
not_taken runs/fixed-probe "d. a probe added" \
  'input.ini, [MOSE-Probes] probe keys (probe1, probe2, ... up to the first missing): 1 read again, 0 at the set-up: the probes cannot change during the run'
if [[ -e runs/fixed-probe/OUTPUT/p1.txt ]]; then
  echo "  FAIL — d. a probe added: the run wrote OUTPUT/p1.txt"; status=1
fi
# e. newrun = false written in a rotating-frame run with res-diter = 25 and a row of an earlier run in OUTPUT/torque.dat
d=runs/start-only
prepare_frame "$d"; write_frame "$d" 30 true rigid-body 1.0e4 '1 1'
rewrite "$d" 's/^sol-diter     = 7$/sol-diter     = 7\nres-diter     = 25/'
echo 'a row of an earlier run' > "$d/OUTPUT/torque.dat"
( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
pid=$!
for ((n = 0; n < 6000; n++)); do
  [[ -f $d/OUTPUT/field1.vtm ]] && break
  kill -0 "$pid" 2>/dev/null || break
  sleep 0.05
done
first=no; late=no
[[ -f $d/OUTPUT/field1.vtm ]] && first=yes
write_frame "$d" 30 false rigid-body 1.0e4 '1 1'
rewrite "$d" 's/^sol-diter     = 7$/sol-diter     = 7\nres-diter     = 25/'
# res-diter = 25 leaves no row of iteration 20 in the residual history: the change came in time if the solution of iteration
# 14 was not there when it was written, so the run had not passed iteration 14
[[ -f $d/OUTPUT/field2.vtm ]] && late=yes
wait "$pid"
rc=$?
head1=$(head -n 1 "$d/OUTPUT/torque.dat" 2>/dev/null)
rows=$(awk '$1 !~ /^#/' "$d/OUTPUT/torque.dat" 2>/dev/null | awk '{print $1}' | paste -sd ' ' -)
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — e. newrun: the change did not come between iterations 7 and 14 (host too slow to poll)"; status=1
elif [[ $rc -eq 0 && -f $d/OUTPUT/field.vtm && ${head1:0:1} == '#' && $rows == 25 ]] \
     && ! command grep -a -q -F 'is not taken' "$d/solve.log"; then
  echo "  ok — e. newrun = false read again: no warning, and torque.dat starts again at iteration 25, as newrun = true says"
else
  echo "  FAIL — e. newrun = false read again: exit code $rc, first line of torque.dat \"$head1\", rows \"$rows\" (required" \
       "the header and the row of iteration 25), warning: $(command grep -a -c -F 'is not taken' "$d/solve.log") (see $d/solve.log)"
  status=1
fi
# f. the same value of a fixed key, written in another way: taken
reread runs/fixed-same "" tol_written tol_start
nw=$(command grep -a -c -F 'is not taken' runs/fixed-same/solve.log)
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — f. same value: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
elif [[ $rc -eq 0 && $last == 30 && $nw -eq 0 ]] \
     && cmp -s runs/box-ref/OUTPUT/residual-history.dat runs/fixed-same/OUTPUT/residual-history.dat; then
  echo "  ok — f. ode-relative-tol-temperature = 1.0e-05 read again after 1e-5: taken, the residual history of the run without the change"
else
  echo "  FAIL — f. same value: exit code $rc, last iteration ${last:-none} (required 30), warnings $nw (required 0)," \
       "residual history of the run without the change: $(cmp -s runs/box-ref/OUTPUT/residual-history.dat \
       runs/fixed-same/OUTPUT/residual-history.dat && echo yes || echo no) (see runs/fixed-same/solve.log)"; status=1
fi
# g. time-accurate, a fixed logical key, switched on during the steady run: not taken
to_accurate() { rewrite "$1" 's/^time-accurate         = false$/time-accurate         = true/'; }
reread runs/fixed-accurate "" to_accurate
not_taken runs/fixed-accurate "g. time-accurate = true" \
  'input.ini, [MOSE-Numerics] time-accurate = true: the run has time-accurate = false from its set-up, and time-accurate cannot change during the run'

echo "-- 7. implicit residual smoothing: on with irs = true or irs-beta > 0, at the start and at every reading --"
irs_start() { rewrite "$1" 's/^riemann-solver        = HLLC$/riemann-solver        = HLLC\nirs-beta              = 0.5/'; }
irs_on()    { irs_start "$1"; }
irs_off()   { rewrite "$1" 's/^irs-beta              = 0.5$/irs-beta              = 0/'; }
reread runs/irs-on stop irs_on
stopped runs/irs-on "a. irs-beta = 0.5 in a run started without IRS, ini-error = stop" \
  'input.ini, [MOSE-Numerics] irs-beta = 0.5: the run started without implicit residual smoothing, which cannot be switched on during the run'
# b. irs = false and irs-beta = 0.5 together at the start, read again unchanged every 5 iterations
d=runs/irs-both-once
prepare "$d"; write_short "$d" "" 'irs = false' 'irs-beta = 0.5'
( cd "$d" && "$SOLVER" > solve.log 2>&1 )
d=runs/irs-both
prepare "$d"; write_short "$d" 5 'irs = false' 'irs-beta = 0.5'
( cd "$d" && "$SOLVER" > solve.log 2>&1 )
rc=$?
last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
if [[ $rc -eq 0 && $last == 10 ]] && cmp -s runs/irs-both-once/OUTPUT/residual-history.dat "$d/OUTPUT/residual-history.dat" \
   && ! command grep -a -q -F 'is not taken' "$d/solve.log"; then
  echo "  ok — b. irs = false with irs-beta = 0.5, read again unchanged: IRS stays on, the residual history of the run without readings"
else
  echo "  FAIL — b. irs = false with irs-beta = 0.5, read again unchanged: exit code $rc, last iteration ${last:-none}," \
       "residual history of the run without readings: $(cmp -s runs/irs-both-once/OUTPUT/residual-history.dat \
       "$d/OUTPUT/residual-history.dat" && echo yes || echo no) (see $d/solve.log)"; status=1
fi
# c. irs-beta = 0 written alone in a run started with irs-beta = 0.5: from the reading of iteration 20 the run of the restart
#    from iteration 20 without IRS
reread runs/irs-off "" irs_off irs_start
d=runs/irs-restart
prepare_frame "$d"; write_ini "$d" 0.5 "" 7 20; irs_start "$d"
( cd "$d" && "$SOLVER" > solve1.log 2>&1 )
write_ini "$d" 0.5 "" 7 10; irs_start "$d"; irs_off "$d"
rewrite "$d" 's/^newrun         = true$/newrun         = false/'
( cd "$d" && "$SOLVER" > solve2.log 2>&1 )
same=yes
for f in residual-history.dat field3.vtm field4.vtm field.vtm vtk/field3Block1.vts vtk/field4Block2.vts \
         vtk/fieldBlock1.vts vtk/fieldBlock2.vts vtk/fieldBlock3.vts vtk/fieldBlock4.vts; do
  if ! cmp -s "runs/irs-off/OUTPUT/$f" "$d/OUTPUT/$f"; then
    same=no
    echo "  differs: $f"
  fi
done
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — c. irs-beta = 0: the change did not come between iterations 7 and 20 (host too slow to poll)"; status=1
elif [[ $rc -eq 0 && $last == 30 && $same == yes ]] && ! command grep -a -q -F 'is not taken' runs/irs-off/solve.log; then
  echo "  ok — c. irs-beta = 0 read again at iteration 20: IRS off, the run goes on as the restart without IRS, byte for byte"
else
  echo "  FAIL — c. irs-beta = 0 read again: exit code $rc, last iteration ${last:-none}, the files of the restart without" \
       "IRS: $same (see runs/irs-off/solve.log)"; status=1
fi

echo "-- 8. the checks of the set-up between keys, made again on the values read --"
mco_start() { rewrite "$1" 's/^riemann-solver        = HLLC$/riemann-solver        = AUSM+M\nriemann-options-Mco   = 0.1/'; }
mco_zero()  { rewrite "$1" 's/^riemann-options-Mco   = 0.1$/riemann-options-Mco   = 0/'; }
reread runs/check-mco stop mco_zero mco_start
stopped runs/check-mco "a. riemann-options-Mco = 0 with AUSM+M, ini-error = stop" \
  'input.ini, [MOSE-Numerics] riemann-options-Mco = 0: the run has riemann-options-Mco = 0.1, and the set-up refuses the value read: Low-Mach Riemann solver selected. Mco (riemann-options-Mco) must be defined in input.'
# b. schmidt = 0 written during a run of a gas of four species without INPUT/diffusion.dat
write_tube() {  # <directory> <schmidt>: input.ini of the tube of test/1D/diff1 (time-accurate), written whole, then moved in place
  cat > "$1/input.ini.new" <<TUBE
[MOSE-Parameters]
iter-threshold = 1010
newrun         = true

[MOSE-Physics]
equations = navier-stokes
schmidt   = $2
prandtl   = 1.0

[MOSE-Numerics]
cfl                   = 0.95
vnn                   = 0.95
time-accurate         = true
time-scheme           = RK2
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
sol-format    = vtk raw
ini-diter     = 1000
ini-error     = stop
shell-diter   = 1000000000
TUBE
  mv "$1/input.ini.new" "$1/input.ini"
}
d=runs/check-schmidt
rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
for f in phase.txt thermo.dat transport.dat composition.txt ic.tec bc.txt; do
  ln -s "../../../../../../1D/diff1/INPUT/$f" "$d/INPUT/$f"
done
write_tube "$d" 1.0
( cd "$d" && exec "$SOLVER" > solve.log 2>&1 ) &
pid=$!
for ((n = 0; n < 6000; n++)); do
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 7 ]] && break
  kill -0 "$pid" 2>/dev/null || break
  sleep 0.05
done
first=no; late=no
[[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 7 ]] && first=yes
write_tube "$d" 0
[[ $(awk 'NF {n = $1} END {print n + 0}' "$d/OUTPUT/residual-history.dat" 2>/dev/null) -ge 1000 ]] && late=yes
wait "$pid"
rc=$?
last=$(awk 'NF {n = $1} END {print n}' "$d/OUTPUT/residual-history.dat" 2>/dev/null)
err='input.ini, [MOSE-Physics] schmidt = 0: the run has schmidt = 1.0, and the set-up refuses the value read: Sc<=0 (multicomponent diffusion) requires a binary diffusion table (INPUT/diffusion.dat)'
if [[ $first != yes || $late == yes ]]; then
  echo "  FAIL — b. schmidt = 0: the change did not come between iterations 7 and 1000 (host too slow to poll)"; status=1
elif [[ $rc -ne 0 && $last == 1000 ]] && command grep -a -q -F "$err" "$d/solve.log" \
     && command grep -a -q -E '(ABORT|STOP): ' "$d/solve.log"; then
  echo "  ok — b. schmidt = 0 with four species and no diffusion table, ini-error = stop: the run stopped at iteration 1000: $err"
else
  echo "  FAIL — b. schmidt = 0 with four species and no diffusion table: exit code $rc, last iteration ${last:-none}" \
       "(required 1000), error \"$err\" (see $d/solve.log)"; status=1
fi

if [[ $status -eq 0 ]]; then
  echo "RuntimeInput: all checks passed"
  exit 0
fi
echo "RuntimeInput: FAILED"
exit 1
