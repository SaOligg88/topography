#!/bin/bash

for sub in subList; do
	dir="/a/documents/connectome/_all/${sub}/XXXXXXXXXXXXX"
	for hemi in RL LR; do
	  for sess in 1 2; do
		  input="XXXXXXXXXXX"
			wb_command -cifti-reduce ${dir}/${input} MEAN ${dir}/mean.dtseries.nii
			wb_command -cifti-reduce ${dir}/${input} STDEV ${dir}/stdev.dtseries.nii
			wb_command -cifti-math '(x - mean) / stdev' -fixnan 0 \
			  -var x ${dir}/${input} -var mean ${dir}/mean.dtseries.nii \
				-select 1 1 -repeat -var stdev ${dir}/stdev.dtseries.nii -select 1 1 -repeat​
			# where is the output filename for the above? 
			rm -f ${dir}/mean.dtseries.nii
			rm -f ${dir}/stdev.dtseries.nii
		done
	done
	output="XXXXXXXX"
	wb_command -cifti-merge ${dir}/${output} #INPUTS...
	rm -f #INPUTS...
	echo "Finished ${sub}"
done
