# Save correlation matrix matrix in Matlab
	
	save('/path/to/filename.mat', '-v7.3', 'data');
	
Where `data` is the variable for the `nxn connectivity matrix`.

# Cluster data using diffusion embedding in python

## Setup ##

Install python packages:

	pip install sys os h5py scipy numpy sklearn --user

Use the `--user` flag to install for local user
	
## Run ##

Set the following variable in the script before running:
    
* n_components_embedding = 25	# number of diffusion embedding components
* comp_min = 2 			# min number of clusters
* comp_max = 20 + 1 		# max number of clusters
* varname = 'data' 		# variable name saved into .mat file 
* filename = '/path/to/filename' 	# Filename *without* suffix '.mat'


	python clustering_embedding_macaque.py

The output will be a .mat file with the variable `results` called: 

	/path/to/filename_results.mat

which can be read into matlab using:

	raw_results = load('/path/to/filename_results.mat');
	results = raw_results.results;

