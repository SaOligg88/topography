# Step 1: Clustering #

Modify the variables at the top, including adding the subject list (comma separated), then run:

	python individual_distance_cluster.py

# Step 2: Find DMN cluster and calculate distance to sensory regions #

In Matlab, run: 

	dist2mask = DoFindPaths_LSD(subjects);

Be sure to clone the entire repository, as the necessary functions are located in `../../utils`

