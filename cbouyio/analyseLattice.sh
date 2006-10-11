#!/bin/bash

# Shell script to wrap all the lattice simulation procedure.

# Usage.
#if [ "$1" == "--help" ] || [ "$1" = "-h" ] ;
#then
#  echo "Usage: $0 transsys_program lattice_size timesteps"
#  echo "This is the minimum requirement for the script to run.
#For more info about the script check the source code."
#  exit 1
#fi

#if (( $# < 3 )) ;
#then
#  echo "$0 : The script works with at least 3 parameters check help."
#  exit 1
#fi


# Start the process.
./latticeSimulator -n 10x10 -t 100 -u 0:1 patternFormation.tra patternFormation_ftable.dat

# Invoke R (applies only to the patternFormation.tra transsys program)
echo "source("\""translattice.r"\"")" | R CMD BATCH
echo "lframe <- readTransLattice("\""patternFormation_ftable.dat"\"")" | R CMD BATCH
echo "postscript("\""patternFormation.ps"\"", width=10, height=10); plotConcentrationSeries(lframe, "\""factor_activator"\"", getConcentrationRange(lframe, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH


# Run the control.
./zeroTranssysDiffusibility patternFormation.tra patternFormation_zeroDiff.tra
./latticeSimulator -n 10x10 -t 100 -u 0:1 patternFormation_zeroDiff.tra patternFormation_zeroDiff_ftable.dat

# Invoke R.
echo "lframe_zero <- readTransLattice("\""patternFormation_zeroDiff_ftable.dat"\"")" | R CMD BATCH
echo "postscript("\""patternFormation_zeroDiff.ps"\"", width=10, height=10); plotConcentrationSeries(lframe_zero, "\""factor_activator"\"", getConcentrationRange(lframe_zero, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH

