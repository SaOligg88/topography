# Save matrix in Matlab
	
	SaveMat('/path/to/filename.mat', data);
	
Where *data* is the variable for the nxn connectivity matrix.

# Cluster data using python

	python clustering_embedding_macaque.py -i /path/to/filename.mat data -o /path/to/cluster -e 25 -c 2 20
	
