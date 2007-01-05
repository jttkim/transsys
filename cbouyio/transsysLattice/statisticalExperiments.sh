#!/bin/bash

# statisticalExperiments.sh
# Script to produce data for the whole statistical significant experiments.

# Run the simulator with an number of different random seeds and then calculates some statistics.

# Check parameters.
if test $# != 4 ;
then
  echo "Usage:";
  echo "  `basename $0` [Latice Size] [Timesteps] [Random Seed Times] [Transsys Program]";
  exit $?;
fi

# Variable assignment.
LATTICESIZE=$1;
TIMESTEPS=$2;
UNI_RANGE=0.0:0.1;
RND_SEED_TIMES=$3;
TPNAME=$4;

# Check for the existance of the transsys program file.
if ! test -e "$TPNAME" ;
then
  echo "Input file $TPNAME does not exist.";
  exit $?
fi

# Extract the basename.
BASENAME=`basename ${TPNAME} .tra`;

if test "$TPNAME" = "$BASENAME" ;
then
  echo "  Error in basename proccesing.
  Specify a transsys program (.tra) filename.";
  exit $?;
fi

# Generate the .R source file.
# Source the tarnslattice.r package.
echo "source("\""~/devel/transsys/trunk/cbouyio/transsysLattice/translattice.r"\"");" | cat > ${BASENAME}_Rsource_Stats.R;
# Initialise the frame name vector.
echo "lframeNames <- list();" | cat >> ${BASENAME}_Rsource_Stats.R; 


# Run the basic experiment for $RND_SEED_TIMES times.
for (( rnd=0 ; rnd<=RND_SEED_TIMES ; rnd++))
do
  if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE -i $TIMESTEPS -r $rnd $TPNAME ${BASENAME}_rs${rnd}_ftable.dat;
    then
      exit $?;
  fi

  # Generate the .R source file.
  echo "lframe${rnd} <- readTransLattice("\""${BASENAME}_rs${rnd}_ftable.dat"\"");" | cat >> ${BASENAME}_Rsource_Stats.R;
  echo "lframeNames <- append(lframeNames, lframe${rnd});" | cat >> ${BASENAME}_Rsource_Stats.R
done
