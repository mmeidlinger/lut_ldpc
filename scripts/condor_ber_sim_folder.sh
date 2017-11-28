#!/bin/bash

FOLDER=${1%/}
NUMJOBS=$2

FILES=$FOLDER/*.ini

echo $FILES

for f in $FILES
do
    echo "Submitting parameter file $f..."
    ./condor_ber_sim $f $NUMJOBS
done
