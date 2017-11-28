#!/bin/sh
# Converts .ens files in a way that the Progressive Edge Growth (PEG) program can read 
# them, runs PEG to generate a code and converts the PEG output to an .alist file

# Binaries
binaryPEG=peg/MainPEG
binaryDAT2ALIST=bin/dat2alist
binaryENS2DEG=bin/ens2deg

# Parameters
numM=$1
numN=$2
codeName=$3
ensFileName=$4
sglConcent=${5:-1}
tgtGirth=${6:-100000}



degFileName="$(basename $ensFileName).deg"

# Uncomment this for debugging outputs
#set -x

# Generate degree file
echo "Generating temporary .deg file from .ens file ..."
$binaryENS2DEG $ensFileName $degFileName

# Run PEG
echo "Running Progressive Edge Growth ..."
$binaryPEG -numM $numM -numN $numN -codeName  ${codeName}.dat  -degFileName  $degFileName -sglConcent $sglConcent -tgtGirth $tgtGirth 

# Create .alist file
echo "Creating .alist file ..."
$binaryDAT2ALIST ${codeName}.dat  $codeName

echo "Deleting temporary files... "
# Delete PEG output
rm -rf ${codeName}.dat
# Delete deg file
rm -rf $degFileName
# Delete log file
rm -rf leftHandGirth.log



