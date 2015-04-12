#!/bin/bash

echo "executable = /scr/litauen1/Dropbox/misc/yeoTopo/lme/geodesic/exactgeodesic/geodesicMatrix" > condor
echo "arguments = /scr/litauen1/Dropbox/misc/yeoTopo/lme/geodesic/exactgeodesic/surf.patch.asc 0" >> condor
echo "universe = vanilla" >> condor
echo "output = /scr/litauen2/projects/distance/condor/0.out" >> condor
echo "error = /scr/litauen2/projects/distance/condor/0.error" >> condor
echo "log = /scr/litauen2//projects/distance/condor/0.log" >> condor
echo "request_memory = 2000" >> condor
echo "request_cpus = 1" >> condor
echo "getenv = True" >> condor
echo "notification = Error" >> condor
echo "queue" >> condor

for i in $(seq 1 29930); 
do
	echo "arguments = /scr/litauen1/Dropbox/misc/yeoTopo/lme/geodesic/exactgeodesic/surf.patch.asc ${i}" >> condor
	echo "output = /scr/litauen2/projects/distance/condor/${i}.out" >> condor
	echo "error = /scr/litauen2/projects/distance/condor/${i}.error" >> condor
	echo "queue" >> condor
	echo "" >> condor
done
