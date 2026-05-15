#!/bin/bash -
#===============================================================================
#
#          FILE: figarolame.sh
#
#         USAGE: run "./figarolame.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run MOSE
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     compile accordingly with the presets file"
  echo " solve                 -->     run MOSE"
  echo " kill                  -->     kill the process" 
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../../../
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
        #echo " -> Background"
        BG=1
        shift
        ;;
    -p | --parallel)
        #echo " -> OpenMP($2)"
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

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi
   
if [[ $1 == solve ]]; then
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  # Run the solver
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
