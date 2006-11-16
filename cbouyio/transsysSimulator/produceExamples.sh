#!/bin/bash

+# $Rev::               $:  Revision of last commit
+# $Author::            $:  Author of last commit
+# $Date$:  Date of last commit


# Changes on 07-11-2005 (r_01)
# Uses the transsys_distances_r_02 script (r_01)
# Automates the procedure of transsys examples
# takes only the .rtp file as an argument
# Minor fixes at 4-11-2005
# CBouyio, UEA, 1-11-2005


if [ $# -eq 0 ]
then
  echo "$0 : The script works only with an .rtp filename"
  exit 1
fi

filename=$1
name=${filename%\.rtp}
./tp_generator_r01 $1 $name.tra

if [ $2 ]
then 
  t=$2
else
  t=200
fi
  
# Invoke the transsys_distances_r_03 program (with n=200 
# as default value)
./transsys_distances_r_03 -n $t $name.tra $name.out
mkdir $name-files

# By something like that i could invoke R also
# echo "postscript('$name.eps'); plot(1:10); dev.off()" | R CMD BATCH

cp $1 $name.tra $name.out $name-files
rm $1 $name.tra $name.out

