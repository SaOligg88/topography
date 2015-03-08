#!/bin/bash

list="0_10 10_20 20_30 30_40 40_50 50_60 60_70 70_80 80_90"

#00_05 05_10 10_15 15_20 20_25 25_30 30_35 35_40 40_45 45_50 50_55 55_60 60_65 65_70 70_75 75_80 80_85 85_90 90_95"

for i in $list
do
	#applytransform -o distDMN_MNI_1mm.xform distDMN_MNI_1mm_${i}.nii distDMN_MNI_1mm_${i}.nii
	makeroi -min 0.01 distDMN_MNI_1mm_${i}_trans.nii.gz
done