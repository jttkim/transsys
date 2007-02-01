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

## Initialise the frame name vector.
#echo "lframeNames <- vector();" | cat >> ${BASENAME}_Rsource_Stats.R;

# Initialise the frame names variable.
FRAMENAMES=""

# Run the basic experiment for $RND_SEED_TIMES times.
for (( RNDSEED=1 ; RNDSEED<=RND_SEED_TIMES ; RNDSEED++))
do
  if ! latticeSimulator -n $LATTICESIZE -t $TIMESTEPS -u $UNI_RANGE -i $TIMESTEPS -r $RNDSEED $TPNAME ${BASENAME}_rs${RNDSEED}_ftable.dat;
    then
      exit $?;
  fi

  # Write to the .R source file.
  echo "lframe${RNDSEED} <- readTransLattice("\""${BASENAME}_rs${RNDSEED}_ftable.dat"\"");" | cat >> ${BASENAME}_Rsource_Stats.R;
  FRAMENAMES=$FRAMENAMES"lframe${RNDSEED}, "

#  echo "lframeNames <- append(lframeNames, lframe${RNDSEED});" | cat >> ${BASENAME}_Rsource_Stats.R
done

# Generate the list of lframes.
echo "lframeNames <- list($FRAMENAMES);" | cat >> ${BASENAME}_Rsource_Stats.R;

