#!/bin/sh


# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit


# Shell script to wrap all the lattice simulation procedure,
# conduct all the control experiments and produce the R source file.

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
UNI_RANGE=0.0:0.1
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



# Run the basic experiment.
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE $TPNAME ${BASENAME}_ftable.dat ;
then
  exit $?
fi

# Generate the .R source file.
echo "source("\""~/devel/transsys/trunk/cbouyio/transsysLattice/translattice.r"\"")" | cat > ${BASENAME}_Rsource.r
echo "lframe <- readTransLattice("\""${BASENAME}_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r



# Generate the zero control tp.
if ! alterTranssysDiffusibility -d 0.0 $TPNAME ${BASENAME}_zeroControl.tra ;
then
  exit $?
fi

# Run the zero controlexperiment.
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE ${BASENAME}_zeroControl.tra ${BASENAME}_zeroControl_ftable.dat ;
then
  exit $?
fi

# Append to the .R source file.
echo "lframeZeroControl <- readTransLattice("\""${BASENAME}_zeroControl_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r


#
## Generate the max control tp.
#if ! alterTranssysDiffusibility -d 1.0 $TPNAME ${BASENAME}_maxControl.tra ;
#then
#  exit $?
#fi
#
## Run the max control experiment.
#if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE ${BASENAME}_maxControl.tra ${BASENAME}_maxControl_ftable.dat ;
#then
#  exit $?
#fi
#
## Append to the .R source file.
#echo "lframeMaxControl <- readTransLattice("\""${BASENAME}_maxControl_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r
#
#
#
## Run the homogenized control experiment.
#if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS $TPNAME ${BASENAME}_homogenControl_ftable.dat ;
#then
#  exit $?
#fi
#
## Append to the .R source file.
#echo "lframeHomogenControl <- readTransLattice("\""${BASENAME}_homogenControl_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.r
#


# Remove generated transsys programs.
rm -rf ${BASENAME}_zeroControl.tra # ${BASENAME}_maxControl.tra

