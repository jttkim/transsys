#!/bin/sh

# Shell script to wrap all the lattice simulation procedure.


# Check parameters.
if test $# != 3 ;
then
  echo "Usage:"
  echo "  `basename $0` [Latice_Size] [Timesteps] [Transsys_Program]"
  exit $?
fi

# Variable assignment.
LATTICESIZE=$1
TIMESTEPS=$2
UNI_RANGE=0:1
TPNAME=$3

# Check for the existance of the transsys program file.
if ! test -e "$TPNAME" ;
then
  echo "Input file $TPNAME does not exist."
  exit $?
fi

# Extract the basename.
BASENAME=`basename ${TPNAME} .tra`

if test "$TPNAME" == "$BASENAME" ;
then
  echo "  Error in basename proccesing.
  Specify a transsys program (.tra) filename."
  exit $?
fi

# Start the process.
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE $TPNAME ${BASENAME}_ftable.dat ;
then
  exit $?
fi

# Generate the .R source file.
echo "source("\""~/devel/transsys/trunk/cbouyio/translattice.r"\"")" | cat > ${BASENAME}_Rsource.r
echo "lframe <- readTransLattice("\""${BASENAME}_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r


# Generate the control.
if ! zeroTranssysDiffusibility $TPNAME ${BASENAME}_control.tra ;
then
  exit $?
fi

# Run the control
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE ${BASENAME}_control.tra ${BASENAME}_control_ftable.dat ;
then
  exit $?
fi

# Append to the .R source file.
echo "lframeControl <- readTransLattice("\""${BASENAME}_control_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r

# Remove some files.
rm -rf ${BASENAME}_control.tra
