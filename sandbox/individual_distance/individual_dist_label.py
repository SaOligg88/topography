#!/usr/bin/python

import os, numpy as np, scipy as sp, nibabel.freesurfer as fs
from sklearn.utils.arpack import eigsh  

# Set defaults:
base_dir = '/scr/liberia1/LEMON_LSD/LSD_rest_surf'
output_base_dir = '/scr/liberia1'
subjects = [26410]
    
for subject in subjects:
    for hemi in ['lh', 'rh']:

		# read in cortical mask
		cort = # nodes

		dataCorr = # load conn mat and mask out only cortex
		fullsize = # length of full cortex

        embedding = DoFiedler(dataCorr[cort,cort]) # see below for details
        del dataCorr
        # reinsert zeros:
        fiedler = np.zeros(fullsize)
        fiedler[cort] = embedding[1] # check if this reads the first eigenvector correctly
        
        # read in distance matrix        
        distmat = # read in
        # get labels from freesurfer for anatomy
        fs_labels = # read in
        label_parietal = fs_labels == # XXX # grab parietal mask
        for i in [label1, label2, etc]:
            label_dist = np.mean(distmat(fs_labels == i))
            # mask fiedler by parietal to get peak of DMN in parietal
            masked_fiedler = fiedler * label_parietal
            if masked_fiedler > mean(fiedler):
                anat_dist = label_dist(max(masked_fiedler)) # does that compute elementwise product?
            else:
                anat_dist = label_dist(min(masked_fiedler)) # does that compute elementwise product?
        
        # save out anat_dist for subject / hemi / anat label
        # also create images for quality control: fiedler, masked_fiedler


def DoFiedler(conn):
    # prep for embedding
    K = (conn + 1) / 2.  
    v = np.sqrt(np.sum(K, axis=1)) 
    A = K/(v[:, None] * v[None, :])  
    del K
    A = np.squeeze(A * [A > 0])

    # diffusion embedding
    n_components_embedding = 2
    lambdas, vectors = eigsh(A, k=n_components_embedding+1)  
    del A
    lambdas = lambdas[::-1]  
    vectors = vectors[:, ::-1]  
    psi = vectors/vectors[:, 0][:, None]  
    lambdas = lambdas[1:] / (1 - lambdas[1:])  
    embedding = psi[:, 1:(n_components_embedding + 1 + 1)] * lambdas[:n_components_embedding+1][None, :]  

    return embedding
        
