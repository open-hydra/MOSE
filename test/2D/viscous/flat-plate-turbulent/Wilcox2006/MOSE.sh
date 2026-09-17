#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#         USAGE: ./MOSE.sh [task] [options]
#   DESCRIPTION: compile / run MOSE for this turbulent flat-plate case
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     recompile MOSE"
  echo " solve                 -->     run MOSE"
  echo " kill                  -->     kill the background process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../../../../
MASTER=$MASTERDIR/bin/MOSE
LOCAL=./bin/MOSE

BG=0
NTHREADS=1
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in
    -b | --background) BG=1; shift ;;
    -p | --parallel)   NTHREADS=$2; shift 2 ;;
    -h | --help)       print_usage ;;
    * ) break ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi

if [[ $1 == solve ]]; then
  cd INPUT
  ln -sfn ../../../../../common/Air-viscous/phase.txt     phase.txt
  ln -sfn ../../../../../common/Air-viscous/thermo.dat    thermo.dat
  ln -sfn ../../../../../common/Air-viscous/transport.dat transport.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi
  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
