
## Setup
To compile the exact geodesic calculation script:
  g++ -Wall -I../geo exactGeodesicMatrix.cpp -o exactGeodesicMatrix

## Prepare surface mesh
The following script loads in the surface, removes the medial wall, and output the file in the proper format:
  prepSurf.m

## Run distance calculation on condor
  ./x.prepCondor.sh

## Reassemble the matrix
  ReassembleDistMat.m
