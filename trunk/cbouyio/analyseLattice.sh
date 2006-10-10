#!/bin/bash

# Shell script to wrap all the lattice simulation procedure.

# Usage.
if [ "$1" == "--help" ] || [ "$1" = "-h" ] ;
then
  echo "Usage: $0 transsys_program lattice_size timesteps"
  echo "This is the minimum requirementsa for the script to run.
For more info about the script check the source code."
  exit 1
fi

if (( $# < 3 )) ;
then
  echo "$0 : The script works with at least 3 parameters check help."
  exit 1
fi

tpName=$1
name=${tpName%\.tra}

latticeSize=$2

timesteps=$3

uBorders=$4

exprValue=$5

if [ $6 ];
then
  identifier=$6
else 
  identifier="factor_activator"
fi

if [ $7 ];
then
  expression=$7
else
  expression="diffusibility"
fi


# Start the process.
./tunningTranssys -i $identifier -v $expression:$exprValue | ./latticeSimulator -n $latticeSize -t $timesteps -u $uBorders $tpName $name\_ftable.dat


# Invoke R (applies only to the patternFormation.tra transsys program)
echo "source('translattice.r')" | R CMD BATCH
echo "lframe <- readTransLattice('$name\_ftable.dat')" | R CMD BATCH
echo "postscript('$name.ps'); plotConcentrationSeries(lframe, 'factor_activator', getConcentrationRange(lframe, 'factor_activator'), oneSecondDelay); dev.off();" | R CMD BATCH

rm $name\_ftable.dat

