#!/bin/bash


for i in $(seq 1 29930); 
do
	echo "arguments = /scr/litauen1/Dropbox/misc/yeoTopo/lme/geodesic/exactgeodesic/surf.patch.asc ${i}" >> condor
	echo "output = /scr/litauen2/projects/distance/condor/${i}.out" >> condor
	echo "error = /scr/litauen2/projects/distance/condor/${i}.error" >> condor
	echo "queue" >> condor
	echo "" >> condor
done
