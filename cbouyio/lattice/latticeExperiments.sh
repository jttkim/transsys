#!/bin/sh

# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit


# Shell script to wrap all the lattice simulation procedure,
# conduct all the control experiments and produce the R source file.

# Check parameters.
if test $# -lt 3 || test $# -gt 4 ;
then
  echo "Usage:"
  echo "  `basename $0` [Latice_Size] [Timesteps] [Transsys_Program] [Sampling_Interval (Optional)]"
  exit $?
fi

# Variable assignment.
LATTICESIZE=$1
TIMESTEPS=$2
SAMPLINGINTERVALS=2
RNDSEED=1
UNI_RANGE=0:0.1
TPNAME=$3

# Assign the sampling interval if it has been specified.
if test $4 ;
then
  SAMPLINGINTERVALS=$4
fi

# Check for the existance of the transsys program file.
if ! test -e "$TPNAME" ;
then
  echo "Input file $TPNAME does not exist."
  exit $?
fi

# Extract the basename.
BASENAME=`basename ${TPNAME} .tra`

if test "$TPNAME" = "$BASENAME" ;
then
  echo "  Error in basename proccesing.
  Specify a transsys program (.tra) filename."
  exit $?
fi


# Run the basic experiment.
echo "Runing the lattice experiment..."
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE -d $SAMPLINGINTERVALS -r $RNDSEED $TPNAME ${BASENAME}_ftable.dat ;
then
  exit $?
  echo "Error on the run of the lattice experiment..."
fi


# Generate the .R source file.
echo "source("\""~/devel/transsys/trunk/cbouyio/transsysLattice/translattice.r"\"")" | cat > ${BASENAME}_Rsource.R
echo "lframe <- readTransLattice("\""${BASENAME}_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.R
echo "Lattice experiment finish succefully!"



# Run the zero diffusibility control experiment.
echo "Runing the individual collection of cells control experiment..."
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE -d $SAMPLINGINTERVALS -r $RNDSEED -o ${BASENAME}.tra ${BASENAME}_zeroControl.dat ;
then
  exit $?
  echo "Error on the run of the individual collection of cells control experiment..."
fi

# Append to the .R source file.
echo "lframeCtrl <- readTransLattice("\""${BASENAME}_zeroControl.dat"\"")" | cat >> ${BASENAME}_Rsource.R
echo "Individual collection of cells experiment finish succefully!"



# Run the well strired experiment.
echo "Runing the well stirred reactor control experiment..."
if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE -d $SAMPLINGINTERVALS -r $RNDSEED -w $TPNAME ${BASENAME}_wellStirred.dat ;
then
  exit $?
  echo "Error on the run of the well stirred reactor control experiment..."
fi

# Append the well stirred results to the .R source file.
echo "lframeWSR <- readTransLattice("\""${BASENAME}_wellStirred.dat"\"")" | cat >> ${BASENAME}_Rsource.R
echo "Well stirred reactor experiment finish succefully!"


# SOME OTHER CONTROLS TO BE TESTED LATER
## Generate the maximum diffusibility control tp.
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
#echo "lframeMaxControl <- readTransLattice("\""${BASENAME}_maxControl_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.R
#
#
## Run the homogenized control experiment (without initial randomisation).
#if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS $TPNAME ${BASENAME}_homogenControl_ftable.dat ;
#then
#  exit $?
#fi
#
## Append to the .R source file.
#echo "lframeHomogenControl <- readTransLattice("\""${BASENAME}_homogenControl_ftable.dat"\"")" | cat >> ${BASENAME}_Rsource.R
#

## Remove the generated transsys programs.
#rm -rf ${BASENAME}_zeroControl.tra    # ${BASENAME}_maxControl.tra

