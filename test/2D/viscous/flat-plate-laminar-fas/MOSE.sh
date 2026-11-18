#!/bin/bash -
#===============================================================================
#
#          FILE: MOSE.sh
#
#         USAGE: run "./MOSE.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run MOSE on the laminar flat plate
#                (../flat-plate-laminar) with the FAS multigrid cycle at its
#                defaults: three grid levels, 2000 iterations on each coarse
#                level, cycle = w (input.ini). INPUT/bc2.txt and bc3.txt are
#                the boundary records of levels 2 and 3, written by ATLAS BCB
#                from the [BCB-Block1] section of input.ini.
#
#   test (CTest FlatPlateLaminar.FAS): the cycle and the grid sequencing to
#     convergence, RUN/fas and RUN/seq (input.ini, the latter with
#     cycle = none): 100000 iterations of the finest level at most, a residual
#     threshold of 1e-10 of the first residual of each level. verify.py of
#     ../flat-plate-laminar must find the skin friction of both within its
#     tolerance of Blasius, 2.5 % RMS for x >= 0.05 m, and compare.py the two
#     skin frictions within 0.5 % RMS of each other: the same solution. The
#     cycle must reach a density residual of 1e-7 (residual-history.dat) in at
#     most half the iterations of the finest level of the grid sequencing.
#   The runs of test write the solution at the end only.
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     compile accordingly with the preset file"
  echo " solve                 -->     run MOSE"
  echo " kill                  -->     kill the process"
  echo " test                  -->     the cycle and the grid sequencing to convergence, against Blasius and each other"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../../
MASTER=$MASTERDIR/bin/MOSE
LOCAL=./bin/MOSE

# Default Options
BG=0
NTHREADS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    -b | --background)
        BG=1
        shift
        ;;
    -p | --parallel)
        NTHREADS=$2
        shift 2
        ;;
    -h | --help)
        print_usage
        shift
        ;;
    * )
      break
      ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

function fail {
  echo "FlatPlateLaminar FAS: FAIL — $*"
  exit 1
}

# run_case <name> <cycle>: input.ini with that cycle, to convergence (residual threshold 1e-10, 100000 iterations of the
# finest level at most), the solution written at the end only, in RUN/<name>
function run_case {
  local d=RUN/$1 f
  rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
  for f in phase.txt thermo.dat transport.dat ic.tec bc.txt bc2.txt bc3.txt; do ln -s "../../../INPUT/$f" "$d/INPUT/$f"; done
  tr -d '\r' < input.ini | sed -e "s/^iter-threshold = .*/iter-threshold = 1000000000\nres-threshold = 1e-10/" \
      -e "s/^sol-diter = .*/sol-diter = 1000000000/" -e "s/^shell-diter = .*/shell-diter = 1000/" \
      -e "s/^level1-iter = .*/level1-iter = 100000/" -e "s/^cycle       = .*/cycle       = $2/" > "$d/input.ini"
  ( cd "$d" && "../../$LOCAL" > solve.log 2>&1 ) || { tail -5 "$d/solve.log"; fail "solver exited non-zero ($d)"; }
}

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi

if [[ $1 == solve || $1 == test ]]; then
  # INPUT/phase.txt, thermo.dat and transport.dat are the gas of input.ini ([GPB-Phase1]:
  # mil = 1e-3, Pr = 1), as written by ATLAS GPB: not the shared air of test/common
  mkdir -p bin
  [[ $1 == solve ]] && mkdir -p OUTPUT
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
fi

if [[ $1 == solve ]]; then
  # Run the solver
  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

# rms <name>: the RMS error of the skin friction of RUN/<name> against Blasius for x >= 0.05 m (verify.py)
function rms {
  local out
  out=$(cd RUN/$1 && python3 -B ../../../flat-plate-laminar/verify.py 2>&1) || { echo "$out" | tail -5 >&2; fail "verify.py ($1)"; }
  echo "$out" | sed -n 's/.*Valid region (x ≥ 0.05 m): *RMS = \([0-9.]*\) %.*/\1/p'
}

# reached <name>: the first iteration of the finest level whose density residual is below 1e-7, or 0
function reached {
  local coarse
  coarse=$(tr -d '\r' < input.ini | awk -F= '/^level[23]-iter/ {n += $2} END {print n}')
  awk -v c="$coarse" 'NF && $1 > c && $3 < 1e-7 {print $1 - c; found = 1; exit} END {if (!found) print 0}' \
    "RUN/$1/OUTPUT/residual-history.dat"
}

if [[ $1 == test ]]; then
  run_case fas w
  run_case seq none
  efas=$(rms fas); eseq=$(rms seq)
  [[ -n $efas && -n $eseq ]] || fail "no RMS error from verify.py"
  awk -v e="$efas" 'BEGIN {exit !(e <= 2.5)}' || fail "the cycle is $efas % RMS from Blasius, above 2.5 %"
  awk -v e="$eseq" 'BEGIN {exit !(e <= 2.5)}' || fail "the grid sequencing is $eseq % RMS from Blasius, above 2.5 %"
  out=$(python3 -B compare.py 2>&1) || { echo "$out" | tail -5; fail "compare.py"; }
  echo "$out"
  dcf=$(echo "$out" | sed -n 's/.*RMS = \([0-9.]*\) %.*/\1/p')
  awk -v e="$dcf" 'BEGIN {exit !(e != "" && e <= 0.5)}' || fail "the skin frictions of the cycle and of the grid sequencing differ by $dcf % RMS, above 0.5 %"
  nfas=$(reached fas); nseq=$(reached seq)
  [[ $nfas -gt 0 && $nseq -gt 0 ]] || fail "a run does not reach a density residual of 1e-7 (cycle: $nfas, grid sequencing: $nseq)"
  (( 2 * nfas <= nseq )) || fail "the cycle reaches 1e-7 at iteration $nfas, the grid sequencing at $nseq: not in half the iterations"
  echo "FlatPlateLaminar FAS: ok — converged skin friction within 2.5 % RMS of Blasius with both methods, $dcf % RMS from each other; density residual 1e-7 at iteration $nfas of the finest level with the cycle, $nseq with the grid sequencing"
fi

if [[ $1 == kill ]]; then
read PID < .ID && kill $PID
fi
