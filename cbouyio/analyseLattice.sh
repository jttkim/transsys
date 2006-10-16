#!/bin/bash

# Shell script to wrap all the lattice simulation procedure.

# Usage.
if test "$1" == "--help"  || test "$1" == "-h" ;
then
  echo "Usage: $0 transsys_program"
  echo "This is the minimum requirement for the script to run.
For more info about the script check the source code."
  exit 1
fi

if test $# -lt 1 ;
then
  echo "$0 : a transsys program file is required to run the script."
  exit 1
fi

tpName=$1;

if ! test -e "$tpName" ;
then
  echo "Input file $tpName does not exist"
  exit 
fi

bname=${tpName%\.tra}

#if test $bname -eq 1 ;
#then
#  echo "Give a transsys program file .tra"
#  exit 1
#fi

#echo $bname

#bname=`basename $1`

#echo $bname

#exit 1


# Start the process.
if ! ./latticeSimulator -n 10x10 -t 10 -u 0:1 $tpName ${bname}_ftable.dat ;
then
  exit $?
fi

# Invoke R (applies only to the patternFormation.tra transsys program)
echo "source("\""translattice.r"\"")" | R CMD BATCH --no-restore
echo "lframe <- readTransLattice("\""${bname}_ftable.dat"\"")" | R CMD BATCH

#echo "postscript("\""patternFormation.ps"\"", width=10, height=10); plotConcentrationSeries(lframe, "\""factor_activator"\"", getConcentrationRange(lframe, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH


# Generate the control.
./zeroTranssysDiffusibility $tpName ${bname}_zeroDiff.tra

if ! ./zeroTranssysDiffusibility $tpName ${bname}_zeroDiff.tra ;
then
  exit $?
fi

# Run the control
if ! ./latticeSimulator -n 10x10 -t 10 -u 0:1 ${bname}_zeroDiff.tra ${bname}_zeroDiff_ftable.dat ;
then
  exit $?
fi

# Invoke R.
echo "lframe_zero <- readTransLattice("\""${bname}_zeroDiff_ftable.dat"\"")" | R CMD BATCH

#echo "postscript("\""patternFormation_zeroDiff.ps"\"", width=10, height=10); plotConcentrationSeries(lframe_zero, "\""factor_activator"\"", getConcentrationRange(lframe_zero, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH

# Call an R procedure which contains all the necessary elements.
R --no-save

