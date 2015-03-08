#!/bin/bash

# /Applications/fsl/bin/flirt -in /Volumes/hd2/Project_Eickhoff/distDMN/Clusters/distDMN_MNI_2mm.nii.gz -applyxfm -init /Applications/fsl/etc/flirtsch/ident.mat -out /Volumes/hd2/Project_Eickhoff/distDMN/Clusters/distDMN_MNI_1mm.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref /Applications/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz


fslmaths distDMN_MNI_1mm.nii.gz -thr 0 -uthr 10 -bin distDMN_MNI_1mm_0_10
fslmaths distDMN_MNI_1mm.nii.gz -thr 10 -uthr 20 -bin distDMN_MNI_1mm_10_20
fslmaths distDMN_MNI_1mm.nii.gz -thr 20 -uthr 30 -bin distDMN_MNI_1mm_20_30
fslmaths distDMN_MNI_1mm.nii.gz -thr 30 -uthr 40 -bin distDMN_MNI_1mm_30_40
fslmaths distDMN_MNI_1mm.nii.gz -thr 40 -uthr 50 -bin distDMN_MNI_1mm_40_50
fslmaths distDMN_MNI_1mm.nii.gz -thr 50 -uthr 60 -bin distDMN_MNI_1mm_50_60
fslmaths distDMN_MNI_1mm.nii.gz -thr 60 -uthr 70 -bin distDMN_MNI_1mm_60_70
fslmaths distDMN_MNI_1mm.nii.gz -thr 70 -uthr 80 -bin distDMN_MNI_1mm_70_80
fslmaths distDMN_MNI_1mm.nii.gz -thr 80 -uthr 90 -bin distDMN_MNI_1mm_80_90


fslmaths distDMN_MNI_1mm.nii.gz -thr 00 -uthr 05 -bin distDMN_MNI_1mm_00_10
fslmaths distDMN_MNI_1mm.nii.gz -thr 05 -uthr 10 -bin distDMN_MNI_1mm_05_10
fslmaths distDMN_MNI_1mm.nii.gz -thr 10 -uthr 15 -bin distDMN_MNI_1mm_10_15
fslmaths distDMN_MNI_1mm.nii.gz -thr 15 -uthr 20 -bin distDMN_MNI_1mm_15_20
fslmaths distDMN_MNI_1mm.nii.gz -thr 20 -uthr 25 -bin distDMN_MNI_1mm_20_25
fslmaths distDMN_MNI_1mm.nii.gz -thr 25 -uthr 30 -bin distDMN_MNI_1mm_25_30
fslmaths distDMN_MNI_1mm.nii.gz -thr 30 -uthr 35 -bin distDMN_MNI_1mm_30_35
fslmaths distDMN_MNI_1mm.nii.gz -thr 35 -uthr 40 -bin distDMN_MNI_1mm_35_40
fslmaths distDMN_MNI_1mm.nii.gz -thr 40 -uthr 45 -bin distDMN_MNI_1mm_40_45
fslmaths distDMN_MNI_1mm.nii.gz -thr 45 -uthr 50 -bin distDMN_MNI_1mm_45_50
fslmaths distDMN_MNI_1mm.nii.gz -thr 50 -uthr 55 -bin distDMN_MNI_1mm_50_55
fslmaths distDMN_MNI_1mm.nii.gz -thr 55 -uthr 60 -bin distDMN_MNI_1mm_55_60
fslmaths distDMN_MNI_1mm.nii.gz -thr 60 -uthr 65 -bin distDMN_MNI_1mm_60_65
fslmaths distDMN_MNI_1mm.nii.gz -thr 65 -uthr 70 -bin distDMN_MNI_1mm_65_70
fslmaths distDMN_MNI_1mm.nii.gz -thr 70 -uthr 75 -bin distDMN_MNI_1mm_70_75
fslmaths distDMN_MNI_1mm.nii.gz -thr 75 -uthr 80 -bin distDMN_MNI_1mm_75_80
fslmaths distDMN_MNI_1mm.nii.gz -thr 80 -uthr 85 -bin distDMN_MNI_1mm_80_85
fslmaths distDMN_MNI_1mm.nii.gz -thr 85 -uthr 90 -bin distDMN_MNI_1mm_85_90
fslmaths distDMN_MNI_1mm.nii.gz -thr 90 -uthr 95 -bin distDMN_MNI_1mm_90_95
