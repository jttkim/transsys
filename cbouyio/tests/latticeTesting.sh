#!/bin/bash

# Script to conduct several tests on the latticeSimulator python program.

# Test with zero transsys program (all the values set to default).
../latticeSimulator exampleZero.tra zeroTMP.out
if diff zeroTMP.out.expect zeroTMP.out ; then
  true
else
  echo 'zero test: files zeroTMP.out.expect zeroTMP.out differ'
  exit 1
fi
echo 'zero test: successful'
rm zeroTMP.out 

# next test goes here...

if false ; then
  true
else
  echo 'nonexistent test: problem'
  exit 1
fi
echo 'nonexistent test: successful'

echo
echo 'all tests run successfully'

