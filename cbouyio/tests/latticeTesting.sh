#!/bin/bash

# Script to conduct several tests on the latticeSimulator python program.

# Test with zero transsys program (all the values set to default).
./latticeSimulator exampleZero.tra zeroTMP.out
z='grep --count 0.0 zeroTMP.out'

if z=50
then
  echo 'Zero Test Succesfully passed!'
else
  echo 'Zero Test Failed!!!'
fi
rm zeroTMP.out 
