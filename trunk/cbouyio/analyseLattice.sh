#!/bin/bash

# Shell script to wrap all the lattice simulation procedure.

# Variables
TIMESTEPS=500
LATTICESIZE=10x10
UNI_RANGE=0:1

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

# Get the transsys program name.
tpName=$1;

if ! test -e "$tpName" ;
then
  echo "Input file $tpName does not exist."
  exit 
fi

# Extract the basename.
bname=`basename $tpName .tra`

if test "$bname" == "$tpname" ;
then
  echo "  Error in basename proccesing.
  Specify a transsys program (.tra) filename."
  exit
fi

# Start the process.
if ! ./latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE $tpName ${bname}_ftable.dat ;
then
  exit $?
fi

# Produce the .R source file.
echo "source("\""translattice.r"\"")" | cat > ${bname}_RData.R
echo "lframe <- readTransLattice("\""${bname}_ftable.dat"\"")" | cat >> ${bname}_RData.R

# Invoke R
#echo "source("\""translattice.r"\"")" | R CMD BATCH --no-restore
#echo "lframe <- readTransLattice("\""${bname}_ftable.dat"\"")" | R CMD BATCH
#echo "postscript("\""patternFormation.ps"\"", width=10, height=10); plotConcentrationSeries(lframe, "\""factor_activator"\"", getConcentrationRange(lframe, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH


# Generate the control.
if ! ./zeroTranssysDiffusibility $tpName ${bname}_control.tra ;
then
  exit $?
fi

# Run the control
if ! ./latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE ${bname}_control.tra ${bname}_control_ftable.dat ;
then
  exit $?
fi

# Produce the .R source file.
echo "lframe_control <- readTransLattice("\""${bname}_control_ftable.dat"\"")" | cat >> ${bname}_RData.R

# Invoke R.
#echo "lframe_zero <- readTransLattice("\""${bname}_control_ftable.dat"\"")" | R CMD BATCH
#echo "postscript("\""patternFormation_control.ps"\"", width=10, height=10); plotConcentrationSeries(lframe_zero, "\""factor_activator"\"", getConcentrationRange(lframe_zero, "\""factor_activator"\""), oneSecondDelay); dev.off();" | R CMD BATCH

# Remove some files.
rm -rf ${bname}_control.tra 

