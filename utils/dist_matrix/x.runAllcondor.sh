#!/bin/bash

dir=`pwd`
cd /scr/ilz2/LEMON_LSD/freesurfer/
List=`ls -d ?????`
cd $dir

for sub in $List
do 
  for hemi in lh rh
  do
    mkdir /scr/liberia1/data/lsd/surface/distance/condor/condor_${sub}_${hemi}
    ./x.prepCondor.sh ${sub} ${hemi}
    # condor_submit condor_${sub}_${hemi}
  done
done
