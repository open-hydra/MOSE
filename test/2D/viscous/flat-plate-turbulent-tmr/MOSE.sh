#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#   DESCRIPTION: turbulent flat plate of the NASA Turbulence Modeling Resource
#                (TMR) on three grids of the TMR family (35x25, 69x49 and
#                137x97 nodes), with the Spalart-Allmaras, SST and Wilcox 2006
#                models
#
#   Every run has its own folder, RUN/<model>-<grid>: make_case.py writes its
#   INPUT/ic.tec and INPUT/bc.txt, input-<model>.ini its input.ini with the
#   iterations of the form, and test/common/Air-viscous its gas.
#
#   test (full form): the nine runs from the uniform field, 300 000 iterations
#     each, a numbered solution every 20 000 (long: about 24.5 h at one thread
#     in a TESTING build; use -p);
#   testfast (short form, for CTest): for each model, a restart of 5000
#     iterations from the converged 35x25 solution of the full form
#     (reference/field-35x25-<model>.tec), in RUN/<model>-35x25-short;
#   solve <model> <grid>: one run of the full form.
#   verify.py reads RUN/form.txt, which this script writes, and checks the form.
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " test                  -->     full form: SA, SST, Wilcox2006 on 35x25, 69x49, 137x97, from the uniform field"
  echo " testfast              -->     short form: SA, SST, Wilcox2006 restarted from reference/field-35x25-<model>.tec"
  echo " solve <model> <grid>  -->     one run of the full form, e.g. solve SA 35x25"
  echo
  echo "Solver options"
  echo " -p | --parallel <n>   -->     launch the solver with <n> threads"
  echo " -i | --iterations <n> -->     iterations of each run (test and solve: 300000, a numbered"
  echo "                               solution every n/15; testfast: 5000)"
  echo
  exit 1
}

MASTERDIR=../../../../
MASTER=$MASTERDIR/bin/MOSE
LOCAL=./bin/MOSE
MODELS="SA SST Wilcox2006"
GRIDS="35x25 69x49 137x97"
ITER_FULL=300000
ITER_SHORT=5000

NTHREADS=1
ITER=""
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in
    -p | --parallel)   NTHREADS=$2; shift 2 ;;
    -i | --iterations) ITER=$2; shift 2 ;;
    -h | --help)       print_usage ;;
    * ) break ;;
  esac
done
[[ $# == 0 ]] && print_usage
[[ -z $ITER || $ITER =~ ^[1-9][0-9]*$ ]] || { echo "MOSE.sh: --iterations must be a positive integer"; exit 1; }

# one_run <model> <grid> <folder> <iterations> <sol-diter> <newrun>
function one_run {
  local model=$1 grid=$2 dir=$3 iter=$4 diter=$5 newrun=$6 f rc
  rm -rf "$dir"
  mkdir -p "$dir/INPUT" "$dir/OUTPUT"
  python3 make_case.py --model "$model" --grid "$grid" --out "$dir" > /dev/null || return 1
  for f in phase.txt thermo.dat transport.dat; do
    ln -sfn ../../../../../../common/Air-viscous/$f "$dir/INPUT/$f"
  done
  sed -e "s/^iter-threshold.*/iter-threshold = $iter/" -e "s/^newrun.*/newrun = $newrun/" \
      -e "s/^sol-diter\( *\)=.*/sol-diter\1= $diter/" "input-$model.ini" > "$dir/input.ini"
  if [[ $newrun == false ]]; then
    cp "reference/field-35x25-$model.tec" "$dir/OUTPUT/field.tec" || return 1
  fi
  ( cd "$dir" && ../../$LOCAL > log.txt 2> err.txt )
  rc=$?
  echo "  $dir: exit code $rc, $(command grep -a -E '^ *Iterations' "$dir/log.txt" | tail -1 | tr -s ' ')"
  if [[ $rc != 0 ]]; then tail -20 "$dir/log.txt" "$dir/err.txt"; return 1; fi
}

if [[ $1 == solve || $1 == test || $1 == testfast ]]; then
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  mkdir -p bin RUN
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi
  status=0
  if [[ $1 == testfast ]]; then
    ITER=${ITER:-$ITER_SHORT}
    echo "testfast $ITER" > RUN/form.txt
    for MODEL in $MODELS; do
      one_run $MODEL 35x25 RUN/$MODEL-35x25-short $ITER 1000000000 false || status=1
    done
  else
    ITER=${ITER:-$ITER_FULL}
    DITER=$(( ITER / 15 > 0 ? ITER / 15 : 1 ))
    if [[ $1 == solve ]]; then
      [[ " $MODELS " == *" $2 "* && " $GRIDS " == *" $3 "* ]] || print_usage
      RUNS="$2:$3"
      echo "solve $ITER" > RUN/form.txt
    else
      RUNS=""
      for MODEL in $MODELS; do for GRID in $GRIDS; do RUNS="$RUNS $MODEL:$GRID"; done; done
      echo "test $ITER" > RUN/form.txt
    fi
    for RUN in $RUNS; do
      IFS=: read -r MODEL GRID <<< "$RUN"
      one_run $MODEL $GRID RUN/$MODEL-$GRID $ITER $DITER true || status=1
    done
  fi
  exit $status
fi
print_usage
