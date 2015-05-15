#!/bin/bash

sub=$1
hemi=$2

echo "executable = /scr/liberia1/workspace/topography/utils/dist_matrix/exactGeodesicMatrix" > condor_${sub}_${hemi}
echo "arguments = /scr/liberia1/data/lsd/surface/mesh/${sub}_${hemi}.asc /scr/liberia1/data/lsd/surface/temp/ ${sub}_${hemi} /scr/liberia1/data/lsd/surface/mesh/${hemi}_${sub}_indices.txt 0" >> condor_${sub}_${hemi}
echo "universe = vanilla" >> condor_${sub}_${hemi}
echo "output = /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}/${sub}_${hemi}_0.out" >> condor_${sub}_${hemi}
echo "error = /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}/${sub}_${hemi}_0.error" >> condor_${sub}_${hemi}
echo "log = /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}/${sub}_${hemi}_0.log" >> condor_${sub}_${hemi}
echo "request_memory = 2000" >> condor_${sub}_${hemi}
echo "request_cpus = 1" >> condor_${sub}_${hemi}
echo "getenv = True" >> condor_${sub}_${hemi}
echo "notification = Error" >> condor_${sub}_${hemi}
echo "queue" >> condor_${sub}_${hemi}
echo "" >> condor_${sub}_${hemi}

# fix for right and left...
# r: 9354
# l: 9361
if [[ $hemi == "lh" ]]
then
  length=9353;
fi
if [[ $hemi == "rh" ]]
then
  length=9360;
fi

for i in $(seq 1 ${length}); 
do
	echo "arguments = /scr/liberia1/data/lsd/surface/mesh/${sub}_${hemi}.asc /scr/liberia1/data/lsd/surface/temp/ ${sub}_${hemi} /scr/liberia1/data/lsd/surface/mesh/${hemi}_${sub}_indices.txt ${i}" >> condor_${sub}_${hemi}
	echo "output = /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}/${sub}_${hemi}_${i}.out" >> condor_${sub}_${hemi}
	echo "error = /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}/${sub}_${hemi}_${i}.error" >> condor_${sub}_${hemi}
	echo "queue" >> condor_${sub}_${hemi}
	echo "" >> condor_${sub}_${hemi}
done
