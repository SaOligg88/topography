
## Setup
To compile the exact geodesic calculation script:
  
    g++ -Wall -I../geo exactGeodesicMatrix.cpp -o exactGeodesicMatrix

## Prepare surface mesh
The following script loads in the surface, removes the medial wall, and output the file in the proper format:
  
%    prepSurf.m
    DoDistIndividual2fsa5(sub, hemi);

## Run distance calculation on condor

./x.runAllcondor.sh

#####################################
#!/bin/bash

dir=`pwd`
cd /scr/ilz2/LEMON_LSD/freesurfer/
List=`ls -d ?????`
cd $dir

for sub in $List
do 
    ./x.prepCondor.sh ${sub} ${hemi}
    condor_submit condor_${sub}_${hemi}
done
########################################

## Reassemble the matrix
  
RunAll('lh');

dataAll = ReassembleDistMat(sub, hemi);

% check if file
% rm 
