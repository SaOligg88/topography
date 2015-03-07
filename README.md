# Topography scripts

## Paths through cluster results: ##

### DoFindPaths.m ###

`[clus, countC, randpath] = DoPaths(n);`

where `n` is the number of permutations,
and a config file is loaded in the script.


### embedding.py ###

To create cluster results based on individual connectivity matrices:

`embedding.py -s <subject> -f <output filebasename>`

where input is a connectivity matrix saved as a .mat file.

### find_similar_paths.py ###

To find most common paths through a cluster map.


## Distance maps ##

### Calculating distance from DMN ###
	

