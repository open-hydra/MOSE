#!/bin/bash -
#===============================================================================
#
#          FILE: test.sh
#
#         USAGE: ./test.sh [options]
#
#   DESCRIPTION: bash srcipt to run the validation tests
#===============================================================================
function print_usage {
  echo "Bash srcipt to run the validation tests"
  echo "Options:"
  echo " -p | --parallel         => launch test simulations in parallel"
  echo
  echo "Usage:"
  echo "   ./test.sh <test>   => run the specific <test> case"
  echo "   ./test.sh 1D       => run the 1D test cases"
  echo "   ./test.sh 2D       => run the 2D test cases"
  echo "   ./test.sh clean    => clean the test directories"
  exit 1
}

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

function check {
  if [[ $esito == "0" ]]; then
    echo -e "${GREEN}success${NC} - $TEST" 
    echo success - $TEST >> $FILE; 
  else 
    echo -e "${RED}fail${NC}    - $TEST" 
    echo fail    - $TEST >> $FILE; 
  fi
}

function clean {
  for TEST in $(ls -d */); do
    cd $TEST
    echo $TEST
    rm -rf bin log* OUTPUT errors_file *.out .ID
    find . -maxdepth 2 -type l -delete
    cd ../
  done
}

# Default
PARAL=0

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case "$1" in
    --parallel | -p )
      shift
      PARAL=1
      ;;

    * )
      break
      ;;
  esac
done
# Check input
[[ $# == 0 ]] && print_usage

MAXNTHREADS=24
MINNTHREADS=1
FILE=$(pwd)/testlog.txt

rm -f $FILE
NTEST=$(ls -l | grep ^d | wc -l)

if [[ $SHELL == *"zsh"* ]]; then
  THREADS=1
  IDLE=1
elif [[ $SHELL == *"bash"* ]]; then
  THREADS=$(nproc --all)
  WORKLOAD=$(awk '{u=$2+$4; t=$2+$4+$5; if (NR==1){u1=u; t1=t;} else print ($2+$4-u1) / (t-t1); }' \
  <(grep 'cpu ' /proc/stat) <(sleep 1;grep 'cpu ' /proc/stat))
  IDLE=$(awk "BEGIN {print $THREADS*(1-$WORKLOAD)}")
  IDLE=${IDLE%.*}
fi

echo 'MOSE TEST SESSION' >> $FILE
date >> $FILE
echo >> $FILE
if [[ $@ == all ]]; then
  echo 'All cases tested' >> $FILE
else
  echo 'Test case:'$@ >> $FILE
fi
echo >> $FILE
echo 'Total number of threads : '$THREADS >> $FILE
echo 'Number of idle threads  : '$IDLE >> $FILE
echo >> $FILE

if [ $PARAL == 0 ] ; then
  echo 'Consecutive test run' >> $FILE
  N=$(( $MAXNTHREADS < $IDLE ? $MAXNTHREADS : $IDLE ))
else
  echo 'Parallel test run' >> $FILE
  N=$(awk "BEGIN {print $IDLE/$NTEST}")
  N=${N%.*}
  N=$(( $MAXNTHREADS < $N ? $MAXNTHREADS : $N ))
  echo ' OpenMP with '$N' threads' >> $FILE
fi
echo >> $FILE

if (( N<MINNTHREADS )); then 
  echo 'Not enough free threads' $N' < '$MINNTHREADS
  echo 'Not enough free threads' $N' < '$MINNTHREADS >> $FILE
  exit
fi

for i in $@; do

if [[ $i == 1D ]]; then

  cd '1D'
  for TEST in $(ls -d */); do 
    cd $TEST
    ./MOSE.sh -b test
    sleep 1
    python3 -B verify.py >> $FILE
    esito=$(echo $?)
    check
    cd ../
  done
  cd ..
  exit

elif [[ $i == clean ]]; then
  cd 1D && clean && cd ../
  cd 2D/euler && clean && cd ../../
  cd 2D/viscous && clean && cd ../../
  rm $FILE
  cd ..
  exit

else

  echo $@ >> $FILE; cd $@
  ./figarolame.sh compile
  ./figarolame.sh -b -p $N solve > /dev/null
  
fi

done
